package Atom;
#ABSTRACT: HackaMol Atom Class
use Moose;
use namespace::autoclean;
use lib 'lib/roles', 'lib/HackaMol/lib'; 
use Carp;
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecRole';
use PeriodicTable qw(@ELEMENTS %ELEMENTS %ATOMIC_MASSES @COVALENT_RADII @VDW_RADII %ATOM_MULTIPLICITY);

has 'symbol'  => (
                  is        => 'rw',
                  isa       => 'Str',
                  predicate => 'has_symbol',
                  lazy      => 1,
                  builder   => '_build_symbol',
                 );

sub _build_symbol {
  # if we are building symbol, Z must exist.  BUILD croaks without one of them
  my $self = shift;
  $self->symbol( _Z_to_symbol($self->Z) );
}

has 'Z'       => (
                  is        => 'rw',
                  isa       => 'Int',
                  predicate => 'has_Z',
                  lazy      => 1,
                  builder   => '_build_Z',
                 );

sub _build_Z {
  # if we are building Z, symbol must exist.  BUILD croaks without one of them
  my $self = shift;
  $self->Z( _symbol_to_Z($self->symbol) );
}

has $_        => (
                  is        => 'rw',
                  isa       => 'Num',
                  predicate => "has_$_",
                  lazy      => 1,
                  builder   => "_build_$_",
                 ) foreach (qw(covalent_radius vdw_radius));

sub _build_covalent_radius {
  # if we are building symbol, Z must exist.  BUILD croaks without one of them
  my $self = shift;
  $self->covalent_radius( _Z_to_covalent_radius($self->Z) );
}

sub _build_vdw_radius {
  # if we are building symbol, Z must exist.  BUILD croaks without one of them
  my $self = shift;
  $self->vdw_radius( _Z_to_vdw_radius($self->Z) );
}

has 'bonds' => (
           traits    => [ 'Array' ],
           is        => 'rw',
           isa       => 'ArrayRef[Bond]',
           default   => sub { [] },
           lazy      => 1,
           weak_ref  => 1,
           handles   =>
           {
             "has_bonds" => 'count'   ,
             "add_bonds" => 'push'    ,
             "get_bonds" => 'get'     ,
             "set_bonds" => 'set'     ,
             "all_bonds" => 'elements',
           "count_bonds" => 'count'   ,
           "break_bonds" => 'delete'  ,
           "clear_bonds" => 'clear'   ,
           },
);

sub BUILD {
  my $self = shift;

  unless ($self->has_Z or $self->has_symbol){
    croak "Either Z or Symbol must be set when calling Atom->new()"; 
  }

  #$self->push_charges($self->charge) if $self->has_charge;
  #$self->push_coords($self->coord)   if $self->has_coord;
  #$self->push_forces($self->force)   if $self->has_force;
  $self->symbol(  _fix_symbol($self->symbol) ) if ($self->has_symbol);
}

sub _symbol_to_Z{
  use PeriodicTable qw(%ELEMENTS);
  my $symbol = shift;
  $symbol = ucfirst( lc($symbol) );
  return $ELEMENTS{$symbol};
}

sub _Z_to_symbol{
  use PeriodicTable qw(@ELEMENTS);
  my $Z = shift;
  return $ELEMENTS[$Z];
}

sub _symbol_to_mass{
  use PeriodicTable qw(%ATOMIC_MASSES);
  my $symbol = shift;
  $symbol = ucfirst( lc($symbol) );
  return $ATOMIC_MASSES{$symbol};
}

sub _fix_symbol{
  return ucfirst( lc( shift ) );
}

sub _Z_to_covalent_radius{
  my $Z = shift;
  # index 1 for single bond length..
  return $COVALENT_RADII[$Z][1]/100;
}

sub _Z_to_vdw_radius{
  my $Z = shift;
  return $VDW_RADII[$Z][1]/100;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Atom

=cut
