package AtomsGroupRole;
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
