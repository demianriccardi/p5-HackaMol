package Atom;
#ABSTRACT: HackaMol Atom Class
use Moose;
use namespace::autoclean;
use lib 'lib/roles', 'lib/HackaMol/lib';
use Carp;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' ), 
     'PhysVecMVRRole', 'BondsAnglesDihedralsRole', 'PdbRole', 'QmRole';
use PeriodicTable
  qw(@ELEMENTS %ELEMENTS %ATOMIC_MASSES @COVALENT_RADII @VDW_RADII %ATOM_MULTIPLICITY);

my @delta_attrs = qw(Z symbol mass vdw_radius covalent_radius);

has 'is_dirty' => (
# when attributes change, the Atom gets dirty. change_symbol, change_Z
# generally, classes that have Atom should decide whether to clean Atom
  is      => 'rw',
  isa     => 'Bool',
  lazy    => 1,
  default => 0, # anytime called, the atom becomes dirty forever!  
);

has 'symbol' => (
    is        => 'rw',
    isa       => 'Str',
    predicate => 'has_symbol',
    clearer   => 'clear_symbol',
    lazy      => 1,
    builder   => '_build_symbol',
);

sub _build_symbol {
    my $self = shift;
    return ( _Z_to_symbol( $self->Z ) );
}

has 'Z' => (
    is        => 'rw',
    isa       => 'Int',
    predicate => 'has_Z',
    clearer   => 'clear_Z',
    lazy      => 1,
    builder   => '_build_Z',
);

sub _build_Z {
    my $self = shift;
    return  ( _symbol_to_Z( $self->symbol ) );
}


has $_ => (
    is        => 'rw',
    isa       => 'Num',
    predicate => "has_$_",
    clearer   => "clear_$_",
    lazy      => 1,
    builder   => "_build_$_",
) foreach (qw(covalent_radius vdw_radius));

sub _build_covalent_radius {
    my $self = shift;
    return ( _Z_to_covalent_radius( $self->Z ) );
}

sub _build_vdw_radius {
    my $self = shift;
    return( _Z_to_vdw_radius( $self->Z ) );
}

sub change_Z {
    my $self = shift;
    my $Z    = shift or croak "pass argument Z to change_Z method";
    $self->_clean_atom;
    $self->Z($Z);
}

sub change_symbol{
    my $self   = shift;
    my $symbol = shift or croak "pass argument symbol to change_Z method";   
    $self->_clean_atom;
    $self->symbol(_fix_symbol($symbol));
}

sub _clean_atom {
    my $self = shift;
    foreach my $clearthis (map {"clear_$_"} @delta_attrs) {
      $self->$clearthis;
    }
    carp "cleaning atom attributes for in place change. setting atom->is_dirty";
    $self->is_dirty(1); 
}

sub BUILD {
    my $self = shift;

    if ($self->has_symbol){
      $self->symbol( _fix_symbol( $self->symbol ) ) ;
      return;
    }  

    return if ($self->has_Z); 

    unless ( $self->has_Z or $self->has_symbol ) {
        croak "Either Z or Symbol must be set when calling Atom->new()";
    }
    return;
}

sub _build_mass {
  my $self = shift;
  return (_symbol_to_mass($self->symbol));
};

sub _symbol_to_Z {
    use PeriodicTable qw(%ELEMENTS);
    my $symbol = shift;
    $symbol = ucfirst( lc($symbol) );
    return $ELEMENTS{$symbol};
}

sub _Z_to_symbol {
    use PeriodicTable qw(@ELEMENTS);
    my $Z = shift;
    return $ELEMENTS[$Z];
}

sub _symbol_to_mass {
    use PeriodicTable qw(%ATOMIC_MASSES);
    my $symbol = shift;
    return $ATOMIC_MASSES{$symbol};
}

sub _fix_symbol {
    return ucfirst( lc(shift) );
}

sub _Z_to_covalent_radius {
    my $Z = shift;

    # index 1 for single bond length..
    return $COVALENT_RADII[$Z][1] / 100;
}

sub _Z_to_vdw_radius {
    my $Z = shift;
    return $VDW_RADII[$Z][1] / 100;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

my $atom1 = Atom->new(
    name    => 'C1',
    charges => [-1.0],
    coords  => [ [ 3.12618, -0.06060, 0.05453 ] ],
    forces  => [ [ 0 , 0, 0 ] ],
    Z       => 6
);

my $atom2 = Atom->new(
    name    => 'Hg1',
    charges => [2.0],
    coords  => [ [ 1.04508, -0.06088, 0.05456  ] ],
    forces  => [ [ 0 , 0, 0 ] ],
    symbol  => 'Hg'
);





