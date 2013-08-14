package QmRole;
#ABSRACT: simple role that provides attributes needed for setting up quantum chemistry calculations
use Moose::Role;
use namespace::autoclean;

has 'basis' => (
                is        => 'rw',
                isa       => 'Str',
                predicate => 'has_basis',
                clearer   => 'clear_basis',
                lazy      => 1,
                default   => '6-31+G*', #to provide example; EMSL Str can be many lines...
               );
has 'ecp'   => (
                is        => 'rw',
                isa       => 'Str',
                clearer   => 'clear_ecp',
                predicate => 'has_ecp',
                );

has 'multiplicity' , is => 'rw', isa => 'Int', lazy=>1, default => 1;
has 'basis_geom'   , is => 'rw', isa => 'Str';

has 'dummy' => (
                is          => 'rw',
                isa         => 'Str',
                predicate   => 'is_dummy',
                clearer     => 'clear_dummy',
               );

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

# instance of class that consumes the QmRole.

$obj->multiplicity(1);

$obj->basis_geom('spherical');

$obj->basis_geom('cartesian');

$obj->dummy('bq');

print "dummy! "     if $obj->is_dummy;

$obj->clear_dummy;

print "not dummy! " unless $obj->is_dummy;

$obj->basis('SDD');

$obj->ecp('SDD');

=head1 DESCRIPTION

QmRole provides attributes that will be useful for setting up interfaces to quantum chemistry
packages.  All attributes are 'rw' and lazy, so they will not contaminate the namespace unless 
called upon. The functionality of the QmRole is expected to be extended in the future. I have
used this role by generating 'basis_atoms' from the binned unique atoms in a molecule, which are
then used to append basis_sets and effective core potentials pulled from the EMSL basis set exchange
as a single Str. https://bse.pnl.gov/bse/portal

A dream is to interface with EMSL library directly.

=attr multiplicity

isa Int that is lazy and rw

=attr basis_geom

isa Str that is lazy and rw

=attr dummy
 
isa Str that is lazy and rw

=attr basis

isa Str that is lazy and rw

=attr ecp

isa Str that is lazy and rw

