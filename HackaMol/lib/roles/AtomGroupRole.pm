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

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

use HackaMol::Angle;

my $atom1 = HackaMol::Atom->new(
    name    => 'O1',
    coords  => [ V( 2.05274, 0.01959, -0.07701 ) ],
    Z       => 8,
);

my $atom2 = HackaMol::Atom->new(
    name    => 'H1',
    coords  => [ V( 1.08388, 0.02164, -0.12303 ) ],
    Z       => 8,
);

my $atom3 = HackaMol::Atom->new(
    name    => 'H2',
    coords  => [ V( 2.33092, 0.06098, -1.00332 ) ],
    Z       => 8,
);

$atom1->push_charges(-0.834);
$_->push_charges(0.417) foreach ($atom1, $atom2);

# instance of class that consumes the AtomGroupRole 

my $group = Class_with_AtomGroupRole->new(atoms=> [$atom1,$atom2,$atom3]);


print $group->count_atoms . "\n";

my @atoms = $group->all_atoms;

$group->do_forall('push_charges',0);

$group->do_forall('push_coords',$group->COM);

print $group->dipole_moment . "\n";

$group->gt(1); # same as $group->do_forall('t',1);

print $group->dipole_moment . "\n";




=head1 DESCRIPTION

PhysVecMVR provides the core attributes and methods shared between Atom 
and Molecule classes. Consuming this role gives Classes a place to store 
coordinates, forces, and charges, perhaps, over the course of a simulation 
or for a collection of configurations for which all other object metadata (name, 
mass, etc) remains fixed. As such, the 't' attribute, described below, is 
important to understand. The PhysVecMVR uses Math::Vector::Real (referred to as MVR), 
which has pure Perl and XS implementations.  MVR::XS is fast with many useful/powerful
overloaded methods. PhysVecMVR leaves many attributes rw so that they may be set and 
reset on the fly. This seems most intuitive from the 
perspective of carrying out computational work on molecules.  Thus, HackaMol bravely 
ignores Moose recommendations to use mostly 'ro' attributes and to generate 
objects as needed.  HackaMol may be coerced to be more rigid in future releases.    

Comparing the PhysVec within Atom and Molecule may be helpful. For both, the PhysVecRol 
generates a little metadata (mass, name, etc.) and an array of coordinates, forces, and 
charges.  For an atom, the array of coordinates gives an atom (with fixed metadata) the ability 
to store multiple [x,y,z] positions (as a function of time, symmetry, distribution, etc.). What 
is the array of coordinates for Molecule? Usually, the coordinates for a molecule will likely 
remain empty (because the atoms that Molecule contains have the more useful coordinates), but we
can imagine using the coordinates array to track the center of mass of the molecule if needed. 

In the following:  Methods with mean_foo msd_foo intra_dfoo out front, carries out some analysis
within $self. Methods with inter_ out front carries out some analysis between
two objects that "does" PhysVecMVR (tests for this should be added) at the
$self->t and $obj->t; a warning is carped if the ts are different 

=array_method push_$_, all_$_, get_$_, set_$_, count_$_, clear_$_ foreach qw(charges coords forces)

  ARRAY traits, respectively: push, get, set, all, elements, clear
  Descriptions for charges and coords follows.  forces analogous to coords.
  
=array_method push_charges

  push value on to charges array

  $obj->push_charges($_) foreach (0.15, 0.17, 0.14, 0.13);

=array_method all_charges

  returns array of all elements in charges array

    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 0.14 0.13

=array_method get_charges

  return element by index from charges array

    print $obj->get_charges(1); # prints 0.17

=array_method set_charges

  set value of element by index from charges array

    $obj->set_charges(2, 1.00);
    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 1.00 0.13

=array_method count_charges

  return number of elements in charges array
  
    print $obj->count_charges; # prints 4 

=array_method clear_charges
  
  clears charges array
    
    $obj->clear_charges;
    print $_ . " " foreach $obj->all_charges; # does nothing 
    print $obj->count_charges; # prints 0

=array_method push_coords

  push value on to coords array

  $obj->push_coords($_) foreach ([0,0,0],[1,1,1],[-1.0,2.0,-4.0], [3,3,3]);

=array_method all_coords

  returns array of all elements in coords array

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

    my @new_coords = map {[$_->[0]+1,$_->[1]+1,$_->[2]+1]} $obj->all_coords;

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach @new_coords; 

=array_method get_coords

  return element by index from coords array

    printf ("%10.3f %10.3f %10.3f \n", @{$obj->get_coords(1)}); # 

=array_method set_coords

  set value of element by index from coords array

    $obj->set_coords(2, [100,100,100]);
    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

=array_method count_coords

  return number of elements in coords array
  
    print $obj->count_coords; # prints 4 

=array_method clear_coords
  
  clears coords array
    
    $obj->clear_coords;
    print $_ . " " foreach $obj->all_coords; # does nothing 
    print $obj->count_coords # prints 0

=attr name

isa Str that is rw. useful for labeling, bookkeeping...

=attr t 

isa Int or ScalarRef that is rw with default of 0

t is intended to describe the current "setting" of the object. Objects that 
consume the PhysVecRol have  arrays of coordinates, forces, and charges to 
allow storage with the passing of time (hence t) or the generation alternative 
configurations.  For example, a crystal lattice can be stored as a single 
object that consumes PhysVec (with associated metadata) along with the 
array of 3d-coordinates resulting from lattice vector translations. 

Experimental: Setting t to a ScalarRef allows all objects to share the same t.  
Although, to use this, a $self->t accessor that dereferences the value would 
seem to be required.  There is nothing, currently, in the core to do so. Not 
sure yet if it is a good or bad idea to do so.

    my $t  = 0;
    my $rt = \$t;  
    $_->t($rt) for (@objects);
    $t = 1; # change t for all objects.

=attr mass 

isa Num that is rw and lazy with a default of 0

=attr xyzfree

isa ArrayRef that is rw and lazy with a default value [1,1,1]. Using this array
allows the object to be fixed for calculations that support it.  For example, to
fix the X and Y coordinates:

  $obj->xyzfree([0,0,1]);

fixing any coordinate will set trigger the is_fixed(1) flag.  

=attr is_fixed

isa Bool that is rw and lazy with a default of 0 (false).

=attr charges

isa ArrayRef[Num] that is lazy with public ARRAY traits described in ARRAY_METHODS

Gives atoms and molecules t-dependent arrays of charges. e.g. store and analyze
atomic charges from a quantum mechanical molecule in several intramolecular 
configurations or a fixed configuration in varied external potentials. 

=attr coords forces

isa ArrayRef[Math::Vector::Real] that is lazy with public ARRAY traits described in ARRAY_METHODS

Gives atoms and molecules t-dependent arrays of coordinates and forces,
for the purpose of analysis.  

=method distance

Takes one argument ($obj2) and calculates the distance using Math::Vector::Real

  $obj1->distance($obj2);

=method angle

Takes two arguments ($obj2,$obj3) and calculates the angle (degrees) between 
the vectors with $obj1 as orgin using Math::Vector::Real.  

  $obj1->angle($obj2,$obj3);

=method intra_dcharges

Calculates the change in charge from initial t ($ti) to final t ($tf). I.e.:

$self->get_charges($tf) - $self->get_charges($ti)

=method mean_charges

No arguments.  Calculates the mean of all stored charges.

=method msd_charges

No arguments.  Calculates the mean square deviation of all stored charges.

=method intra_dcoords intra_dforces 

returns the difference (Math::Vector::Real object) from the initial t ($ti) to
the final t ($tf).

$obj1->intra_dcoords($ti,$tf);

$obj1->intra_dforces($ti,$tf);

=method mean_coords mean_forces 

No arguments. Calculates the mean (Math::Vector::Real object) vector.

=method msd_coords msd_forces 

No arguments. Calculates the mean (Math::Vector::Real object) dot product of the
vector difference of the xyz coords from the mean vector.

=method charge

called with no arguments.  returns $self->get_charges($self->t);

=method xyz

called with no arguments.  returns $self->get_coords($self->t);

=method force

called with no arguments.  returns $self->get_forces($self->t);

=method copy_ref_from_t1_through_t2

called as $obj->copy_ref_from_t1_through_t2($qcf,$t,$tf); where qcf is
charges|coords|forces.  Fills the qcf from t+1 to tf with the value at t.
Added this while trying to compare a dipole as a function of changing charges
at a single set of coordinates.

