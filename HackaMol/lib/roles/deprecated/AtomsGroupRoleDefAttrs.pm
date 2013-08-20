package AtomsGroupRole;
use Moose::Role;
use Carp;
use MooseX::Storage;
use Math::Vector::Real;
with Storage( 'io' => 'StorableFile' );

requires '_clear_group_attrs';

my $angste_debye = 4.80320;

has 'atoms' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Atom]',
    default => sub { [] },
    handles => {
        push_atoms            => 'push',
        quick_push_atoms      => 'push',  # for quick pushing with no binning each time
        get_atoms             => 'get',
        set_atoms             => 'set',
        quick_set_atoms       => 'set',
        delete_atoms          => 'delete',
        quick_delete_atoms    => 'delete',
        all_atoms             => 'elements',
    #    atoms                 => 'elements',
        count_atoms           => 'count',
        clear_atoms           => 'clear',
    },
    lazy     => 1,
);

#anytime the group changes, we need to reset the defaults!
after $_ => sub {
  my $self = shift;
  $self->_clear_group_attrs;
  $self->bin_atoms;
} foreach (qw(push_atoms set_atoms delete_atoms clear_atoms));

has $_ => (
    is      => 'rw',
    isa     => 'Math::Vector::Real',
    builder => "_build_$_",
    clearer => "clear_$_",
    lazy => 1,
) foreach (qw(dipole COM COZ));

sub _build_dipole {
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

sub _build_COM {
    my $self      = shift;
    return(V(0)) unless ($self->count_atoms);
    my @atoms     = $self->all_atoms;
    my @m_vectors = map { $_->mass * $_->get_coords( $_->t ) } @atoms;
    my $com       = V( 0, 0, 0 );
    $com += $_ foreach @m_vectors;
    return ($com/$self->total_mass);
}

sub _build_COZ {
    my $self      = shift;
    return(V(0)) unless ($self->count_atoms);
    my @atoms     = $self->all_atoms;
    my @z_vectors = map { $_->Z * $_->get_coords( $_->t ) } @atoms;
    my $coz       = V( 0, 0, 0 );
    $coz += $_ foreach @z_vectors;
    return ($coz/$self->total_Z);
}

sub do_forall{
  my $self   = shift;
  my $method = shift;
  do{carp "doing nothing for all"; return} unless(@_);
  my @atoms = $self->all_atoms;
  $_->$method(@_) foreach @atoms;
}

has $_ => (
    is      => 'rw',
    isa     => 'Num',
    builder => "_build_$_",
    clearer => "clear_$_",
    lazy => 1,
) foreach (qw(dipole_moment total_charge total_mass total_Z));

sub _build_total_charge {
    my $self    = shift;
    return(0) unless ($self->count_atoms);
    my @atoms   = $self->all_atoms;
    my @charges = map { $_->get_charges( $_->t ) } @atoms;
    my $sum     = 0;
    $sum += $_ foreach @charges;
    return ($sum);
}

sub _build_total_mass {
    my $self   = shift;
    return(0) unless ($self->count_atoms);
    my @masses = map { $_->mass } $self->all_atoms;
    my $sum    = 0;
    $sum += $_ foreach @masses;
    return ($sum);
}

sub _build_total_Z {
    my $self = shift;
    return(0) unless ($self->count_atoms);
    my @Zs   = map { $_->Z } $self->all_atoms;
    my $sum  = 0;
    $sum    += $_ foreach @Zs;
    return ($sum);
}

sub _build_dipole_moment {
    my $self = shift;
    return ( abs( $self->dipole )*$angste_debye );
}


has 'atoms_bin' => (
    traits  => ['Hash'],
    isa     => 'HashRef',
    default => sub{{}},
    clearer => 'clear_atoms_bin',
    handles => {
        set_atoms_bin      => 'set',
        get_atoms_bin      => 'get',
        has_empty_bin      => 'is_empty',
        count_unique_atoms => 'count',
        all_unique_atoms   => 'keys',
        atom_counts        => 'kv',
        exists_in_bin      => 'exists',
    },
    lazy => 1,
);

sub bin_atoms {
    my $self  = shift;

    return ({}) unless $self->count_atoms;
    $self->clear_atoms_bin;
    foreach my $atom ($self->all_atoms){
      my $symb = $atom->symbol;
      if ($self->exists_in_bin($symb)){
        my $count_z = $self->get_atoms_bin($symb);
        $count_z->[0]++;
        $self->set_atoms_bin($symb => $count_z);
      }
      else {
        $self->set_atoms_bin($symb => [1, $atom->Z]);
      }
    }

}

sub canonical_name {
    # return something like C4H10 sort in order of descending Z
    my $self = shift;
    my @names = map { 
                     my $name = $_->[0] . $_->[1][0] ; 
                     $name =~ s/(\w+)1$/$1/; $name; 
                    }
                    sort {
                          $b->[1][1] <=> $a->[1][1]    # sort by Z!  see above...
                         } $self->atom_counts;
                return join( '', @names );
}

no Moose::Role;

1;
