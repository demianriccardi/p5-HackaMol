package AtomsGroup;
use Moose::Role;
use Carp;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' );

has 'atoms' => (
    traits  => ['Array'],
    isa     => 'ArrayRef[Atom]',
    default => sub { [] },
    handles => {
        add_atoms    => 'push',
        get_atoms    => 'get',
        set_atoms    => 'set',
        delete_atoms => 'delete',
        all_atoms    => 'elements',
        count_atoms  => 'count',
        clear_atoms  => 'clear',
    },
);

has $_ => (
    is      => 'rw',
    isa     => 'Math::Vector::Real',
    builder => "_build_$_",
    clearer => "_clear_$_",
) foreach (qw(dipole COM COZ));

sub _build_dipole {
    my $self    = shift;
    my @atoms   = $self->all_atoms;
    my @vectors = map { $_->get_coords( $_->t ) } @atoms;
    my @charges = map { $_->get_charges( $_->t ) } @atoms;
    croak "mismatch number of coords and charges" if ( $#vectors != $#charges );
    my $dipole = V( 0, 0, 0 );
    $dipole += $vectors[$_] * $charges[$_] foreach 0 .. $#charges;
    return ($dipole);
}

sub _build_COM {
    my $self      = shift;
    my @atoms     = $self->all_atoms;
    my @m_vectors = map { $_->mass * $_->get_coords( $_->t ) } @atoms;
    my $com       = V( 0, 0, 0 );
    $com += $_ foreach @m_vectors;
    return ($com);
}

sub _build_COZ {
    my $self      = shift;
    my @atoms     = $self->all_atoms;
    my @z_vectors = map { $_->Z * $_->get_coords( $_->t ) } @atoms;
    my $coz       = V( 0, 0, 0 );
    $coz += $_ foreach @z_vectors;
    return ($coz);
}

has $_ => (
    is      => 'rw',
    isa     => 'Num',
    builder => "_build_$_",
    clearer => "_clear_$_",
) foreach (qw(dipole_moment total_charge));

sub _build_total_charge {
    my $self    = shift;
    my @atoms   = $self->all_atoms;
    my @charges = map { $_->get_charges( $_->t ) } @atoms;
    my $sum     = 0;
    $sum += $_ foreach @charges;
    return ($sum);
}

sub _build_dipole_moment {
    my $self = shift;
    return ( abs( $self->dipole ) );
}

has 'atoms_bin' => (
    traits  => ['Hash'],
    isa     => 'HashRef[Str]',
    builder => '_build_atoms_bin',
    clearer => 'clear_bin',
    handles => {
        set_atoms_bin      => 'set',
        get_atoms_bin      => 'get',
        has_empty_bin      => 'is_empty',
        count_unique_atoms => 'count',
        all_unique_atoms   => 'keys',
        atom_counts        => 'kv',
    },
);

sub _build_atoms_bin {
    my $self  = shift;
    my @atoms = $self->all_atoms;
    foreach my $atom (@atoms) {
        my $sym     = $atom->symbol;
        my $count_Z = $self->get_atoms_bin($sym);
        $self->set_atoms_bin( $sym => [ $count_Z->[0]++, $atom->Z ] );
    }
}

sub canonical_name {

    # return something like C4H10 sort in order of descending Z
    my $self = shift;
    my @names = map { $_->[0] . $_->[1] }
      sort {
        $b->[1][1] <=> $a->[1][1]    # sort by Z!  see above...
      } $self->atom_counts;
    return join( '', @names );
}

no Moose::Role;

1;
