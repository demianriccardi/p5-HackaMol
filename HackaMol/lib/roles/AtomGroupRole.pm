package AtomGroupRole;
#ABSTRACT: Role for a group of atoms   
use Moose::Role;
use Carp;
use MooseX::Storage;
use Math::Vector::Real;
with Storage( 'io' => 'StorableFile' );

my $angste_debye = 4.80320;

has 'atoms' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Atom]',
    default => sub { [] },
    handles => {
        push_atoms            => 'push',
        get_atoms             => 'get',
        set_atoms             => 'set',
        delete_atoms          => 'delete',
        all_atoms             => 'elements',
        count_atoms           => 'count',
        clear_atoms           => 'clear',
    },
    lazy     => 1,
);

sub dipole {
    my $self    = shift;
    return(V(0)) unless ($self->count_atoms);
    my @atoms   = $self->all_atoms;
    my @vectors = grep {defined} map { $_->get_coords(  $_->t ) } @atoms;
    my @charges = grep {defined} map { $_->get_charges( $_->t ) } @atoms;
    my $dipole = V( 0, 0, 0 );
    if ( $#vectors != $#charges ){
      carp "build_dipole> mismatch number of coords and charges. all defined?";
      return $dipole;
    }
    $dipole += $vectors[$_] * $charges[$_] foreach 0 .. $#charges;
    return ($dipole);
}

sub COM {
    my $self      = shift;
    return(V(0)) unless ($self->count_atoms);
    my @atoms     = $self->all_atoms;
    my @m_vectors = map { $_->mass * $_->get_coords( $_->t ) } @atoms;
    my $com       = V( 0, 0, 0 );
    $com += $_ foreach @m_vectors;
    return ($com/$self->total_mass);
}

sub COZ {
    my $self      = shift;
    return(V(0)) unless ($self->count_atoms);
    my @atoms     = $self->all_atoms;
    my @z_vectors = map { $_->Z * $_->get_coords( $_->t ) } @atoms;
    my $coz       = V( 0, 0, 0 );
    $coz += $_ foreach @z_vectors;
    return ($coz/$self->total_Z);
}

sub gt {
#set group time
  my $self = shift;
  my $t    = shift;
  $self->do_forall('t',$t);
}

sub do_forall{
  my $self   = shift;
  my $method = shift;
  do{carp "doing nothing for all"; return} unless(@_);
  my @atoms = $self->all_atoms;
  $_->$method(@_) foreach @atoms;
}

sub total_charge {
    my $self    = shift;
    return(0) unless ($self->count_atoms);
    my @atoms   = $self->all_atoms;
    my @charges = map { $_->get_charges( $_->t ) } @atoms;
    my $sum     = 0;
    $sum += $_ foreach @charges;
    return ($sum);
}

sub total_mass {
    my $self   = shift;
    return(0) unless ($self->count_atoms);
    my @masses = map { $_->mass } $self->all_atoms;
    my $sum    = 0;
    $sum += $_ foreach @masses;
    return ($sum);
}

sub total_Z {
    my $self = shift;
    return(0) unless ($self->count_atoms);
    my @Zs   = map { $_->Z } $self->all_atoms;
    my $sum  = 0;
    $sum    += $_ foreach @Zs;
    return ($sum);
}

sub dipole_moment {
    my $self = shift;
    return ( abs( $self->dipole )*$angste_debye );
}

sub bin_atoms {
    my $self  = shift;
    my $bin_hr = {};
    my $z_hr   = {};
    return ($bin_hr,$z_hr) unless $self->count_atoms;
    foreach my $atom ($self->all_atoms){
      $bin_hr->{$atom->symbol}++;
      $z_hr->{$atom->symbol}=$atom->Z;
    }
    return ($bin_hr,$z_hr);
}

sub count_unique_atoms {
    my $self = shift;
    my ($bin_hr,$z_hr) = $self->bin_atoms;
    return (scalar(keys %{$bin_hr}));
}

sub canonical_name {
    # return something like C4H10 sort in order of descending Z
    my $self = shift;
    my ($bin_hr,$z_hr) = $self->bin_atoms;
    my @names = map { 
                     my $name = $_ . $bin_hr->{$_}; 
                     $name =~ s/(\w+)1$/$1/; $name; # substitue 1 away? 
                    }
                    sort {
                          $z_hr->{$b} <=> $z_hr->{$a}    # sort by Z!  see above...
                         } keys %{$bin_hr};
                return join( '', @names );
}

sub translate {
  my $self = shift;
  my $tvec = shift or croak "pass MVR translation vector";
  my $tf   = shift;

  my @atoms = $self->all_atoms; 
  $tf = $atoms[0]->t unless(defined($tf));

  $_->set_coords($tf, $_->xyz+$tvec) foreach ($self->all_atoms);
}

sub rotate {
  #rotate about origin. having origin allow rotation of subgroup
  #without having to translate everything.  the origin could be taken
  #as that of the atoms, but this may be tricky and needs some thinking
  #wrt to API
  my $self = shift;
  my $rvec = shift or croak "pass MVR rotation vector";
  my $ang  = shift or croak "pass rotation angle";
  my $orig = shift or croak "pass MVR origin";
  my $tf   = shift;

  my @atoms = $self->all_atoms;
  my $t = $atoms[0]->t;
  $tf = $t unless(defined($tf));
  $rvec = $rvec->versor; #unit vector

  my @cor  = map{$_->get_coords($t)-$orig} @atoms; 
  my @rcor = $rvec->rotate_3d( deg2rad($ang), @cor );

  $atoms[$_]->set_coords($tf, $rcor[$_]+$orig) foreach 0 .. $#rcor;
}


no Moose::Role;

1;

__END__

=head1 SYNOPSIS

my $atom1 = HackaMol::Atom->new(
    name    => 'O1',
    coords  => [ V( 2.05274, 0.01959, -0.07701 ) ],
    Z       => 8,
);

my $atom2 = HackaMol::Atom->new(
    name    => 'H1',
    coords  => [ V( 1.08388, 0.02164, -0.12303 ) ],
    Z       => 1,
);

my $atom3 = HackaMol::Atom->new(
    name    => 'H2',
    coords  => [ V( 2.33092, 0.06098, -1.00332 ) ],
    Z       => 1,
);

$atom1->push_charges(-0.834);
$_->push_charges(0.417) foreach ($atom1, $atom2);

# instance of class that consumes the AtomGroupRole 

my $group = Class_with_AtomGroupRole->new(atoms=> [$atom1,$atom2,$atom3]);

print $group->count_atoms . "\n"; #3

print $group->total_charge . "\n"; # 0

print $group->total_mass . "\n";  

my @atoms = $group->all_atoms;

print $group->dipole_moment . "\n";

$group->do_forall('push_charges',0);

$group->do_forall('push_coords',$group->COM);

$group->gt(1); # same as $group->do_forall('t',1);

print $group->dipole_moment . "\n";

print $group->canonical_name . "\n";

print $group->unique_atoms . "\n";

$group->translate(V(10,0,0));

$group->rotate( V(1,0,0),
                     180,
                V(0,0,0));

=head1 DESCRIPTION

The HackaMol AtomGroupRole class provides core methods and attributes for 
consuming classes that use groups of atoms. The original implementation of 
this role relied heavily on attributes, builders, and clearers.  Such an approach
naturally gives fast lookup tables, but the ability to change atoms and coordinates
made the role to difficult.  Such an approach may be pursued again (without changing
the API) in the future after the API has matured.  The AtomGroupRole calculates all
values for atoms using their own t attributes.

=array_method push_atoms, get_atoms, set_atoms, all_atoms, count_atoms, clear_atoms

ARRAY traits for the atoms attribute, respectively: push, get, set, elements, count, clear

=array_method push_atoms

push atom on to atoms array

$group->push_atoms($atom1, $atom2, @otheratoms);

=array_method all_atoms

returns array of all elements in atoms array

print $_->symbol, "\n" foreach $group->all_atoms; 

=array_method get_atoms

return element by index from atoms array

print $group->get_atoms(1); # returns $atom2 from above

=array_method set_atoms

set atoms array by index

$group->set_atoms(1, $atom1);

=array_method count_atoms

return number of atoms in group 
  
print $group->count_atoms; 

=array_method clear_atoms

clears atoms array

=method do_for_all

pass method and arguments down to atoms in group

$group->do_for_all('t',1); #sets t to 1 for all atoms

=method gt
  
integer argument. wraps do_for_all for setting time within group

$group->gt(1);

=method dipole

no arguments. return dipole calculated from charges and coordinates as Math::Vector::Real object  

=method COM

no arguments. return center of mass calculated from masses and coordinates as Math::Vector::Real object  
 
=method COM

no arguments. return center of nuclear charge calculated from Zs and coordinates as Math::Vector::Real object  

=method total_charge

no arguments. return sum of atom charges.

=method total_mass

no arguments. return sum of atom masses.

=method total_Z

no arguments. return sum of Zs.

=method dipole_moment

no arguments. returns the norm of the dipole in debye (assuming charges in electrons, AKMA)

=method bin_atoms

no arguments. returns two hash references. The histogram of atom symbols, and a map from symbol-> Z for the same
keys.  The second hash reference was added, to be able to sort by Z in the absence of Atom objects.

=method count_unique_atoms

no arguments. returns the number of keys in each hash returned by bin_atoms

=method canonical_name

no arguments. returns a string summary of the atoms in the group.  Take the bin_atoms hashes, sorts
by Z and generates something like OH2 for water or O2H2 for peroxide.

=attr atoms

isa ArrayRef[Atom] that is lazy with public ARRAY traits described in ARRAY_METHODS

=method translate

requires Math::Vector::Real vector argument. Optional argument: integer tf.  

Translates all atoms in group by the MVR vector.  Pass tf to the translate method to store new 
coordinates in tf rather than atom->t.

=method rotate

requires Math::Vector::Real vector, an angle (in degrees), and a MVR vector origin as arguments. 
Optional argument: integer tf.  

Rotates all atoms in the group around the MVR vector. Pass tf to the translate method to store new 
coordinates in tf rather than atom->t.

