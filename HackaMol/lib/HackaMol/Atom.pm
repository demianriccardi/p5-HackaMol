package Atom;
#ABSTRACT: HackaMol Atom Class
use Moose;
use namespace::autoclean;
use lib 'lib/roles', 'lib/HackaMol/lib';
use Carp;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' ), 
     'PhysVecMVRRole', 'PdbRole', 'QmRole';
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

has 'bond_count' => (
      traits  => ['Counter'],
      is      => 'ro',
      isa     => 'Num',
      default => 0,
      handles => {
          inc_bond_count   => 'inc',
          dec_bond_count   => 'dec',
          reset_bond_count => 'reset',
      },
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

use HackaMol::Atom;
use Math::Vector::Real;


my $atom1 = HackaMol::Atom->new(
    name    => 'Zinc',
    coords  => [ V( 2.05274, 0.01959, -0.07701 ) ],
    Z       => 30,
);

print $atom->symbol ; #prints "Zn"

print "clean " unless $atom->is_dirty; #prints clean

$atom->change_symbol("Hg");

print $atom->Z ; #prints 80

print "dirty " if $atom->is_dirty; #prints dirty

=head1 DESCRIPTION

Central to HackaMol, the Atom class provides methods and attributes for a given atom. 
The Atom class consumes PhysVecMVRRole, PdbRole, and QmRole. See the documentation of those 
roles for details.  The Atom class builds on the capabilities of those roles to add
attributes (such as symbol, Z, covalent radius) and methods (such as change_symbol) 
specific to atoms. Creating an instance of an Atom object requires either the atomic 
number (Z) or the atomic symbol. The other attributes are lazily built when needed.  
An Atom object is flexible. The atom type can be changed in place 
(e.g. convert a zinc atom to a mercury atom, see SYNOPSIS), but changing the type of atom will set
the is_dirty flag so that other objects using the atom have the ability to know whether
atom-type dependent attributes need to be updated (e.g. forcefield parameters).  Atom data
is generated from the PeriodicTable module that borrows data from PerlMol.  The PeriodicTable 
module may be extended in the future.

=attr is_dirty

isa Bool that is lazy and rw.  Default is 0.  $self->is_dirty(1) called during the change_symbol and 
change_Z methods.

=attr symbol

isa Str that is lazy and rw. _build_symbol builds the default. 

Generating an atom instance with symbol, will run ucfirst(lc ($symbol)) to make sure
the format is correct.  Thus, creating an atom object is slightly slower with symbol than with Z. If Z is 
used to generate the instance of the Atom class (my $atom = Atom->new(Z=>1)), the _build_symbol method generates the symbol from Z only when the symbol attribute is read (symbol attribute is lazy).

=attr Z

isa Int that is lazy and rw. _build_Z builds the default

Z is the Atomic number.

=attr covalent_radius

isa Num that is lazy and rw. _build_covalent_radius builds the default.

the covalent radii are taken from those tabulated in:

=for :list
* P. Pyykkö, M. Atsumi (2009). 
"Molecular Single-Bond Covalent Radii for Elements 1-118". Chemistry: A European Journal 15: 186.
* P. Pyykkö, M. Atsumi (2009). 
"Molecular Double-Bond Covalent Radii for Elements Li–E112". Chemistry: A European Journal 15 (46): 12770.
* P. Pyykkö, S. Riedel, M. Patzschke (2005). 
"Triple-Bond Covalent Radii". Chemistry: A European Journal 11 (12): 3511. 

=attr vdw_radius

isa Num that is lazy and rw. _build_vdw_radius builds the default. 

Atomic vdw radii information will be revisited and revised. Included as reminder for now. See the source of PeriodicTable.pm for more information.

=method change_Z

no arguments.  Changes the atom type using Z.  change_Z calls _clean_atom which clears all attributes and
sets calls is_dirty(1).

=method change_symbol

no arguments.  Changes the atom type using symbol. Similar to change_Z, change_symbol calls _clean_atom which clears all attributes and
sets calls is_dirty(1).

=head1 SEE ALSO

=for :list
* L<PhysVecMVRRole>
* L<PdbRole>
* L<QmRole>
* L<PerlMol>

