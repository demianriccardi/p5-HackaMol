package HackaMol::QmMolRole;
#ABSTRACT: provides attributes needed for quantum chemistry calculations
# this will need updating as needs arise
use Moose::Role;

with 'HackaMol::QmAtomRole';

has 'multiplicity' , is => 'rw', isa => 'Int', lazy=>1, default => 1;

my @tscl = qw(
              total_energy electronic_energy nuclear_energy 
              qm_dipole_moment ionization_energy gradient_norm
              heat_of_formation 
              U H G S
              S_t
              S_r
              S_v
              total_energy_mp2
              total_energy_ccsdt
              nonelectrostatic_energy
              );

has "$_" => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Num]',
    default => sub { [] },
    handles => {
        "push_$_"  => 'push',
        "get_$_"   => 'get',
        "all_$_"   => 'elements',
        "clear_$_" => 'clear',
        "count_$_" => 'count',
    },
    lazy   => 1,
) foreach @tscl;

has "$_" => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[Math::Vector::Real]',
    default => sub { [] },
    handles => {
        "push_$_"  => 'push',
        "get_$_"   => 'get',
        "all_$_"   => 'elements',
        "clear_$_" => 'clear',
        "count_$_" => 'count',
    },
    lazy   => 1,
) for qw(qm_dipole frequencies eigvec alpha beta);


no Moose::Role;

1;

__END__

=head1 SYNOPSIS

# instance of class that consumes the QmMolRole.

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

QmMolRole provides attributes that will be useful for setting up interfaces to quantum chemistry
packages.  This role consumes the QmAtomRole so there is some overlap for basis_geom, basis, 
and ecp.  For interfaces, the Molecule should take precedence over the atom; i.e. if a Molecule 
has a basis of 6-31G*, that should be used for all atoms regardless of the basis
set that they may have. All attributes are 'rw' and lazy, so they will not contaminate the 
namespace unless called upon. QmMolRole has 'basis_atoms' that should be generated from the 
binned unique atoms in a molecule.  basis_atoms are just instances of the HackaMol::Atom class 
with all the basis sets and effective core potentials loaded, either as simple
strings supported by the package or the full descriptions pulled from the EMSL basis set 
exchange as a single Str. https://bse.pnl.gov/bse/portal

A dream is to interface with EMSL library directly. Attributes below are
described without much detail; they will contain information mapped from
calculations and are not exhaustive.  This role will probably evolve as 
interfaces are added.

=attr multiplicity

isa Int that is lazy and rw

=attr total_energy electronic_energy nuclear_energy 

each isa ArrayRef[Num] that is lazy with public ARRAY traits: push_$_ get_$_
all_$_ clear_$_ count_$_

=attr qm_dipole_moment ionization_energy gradient_norm heat_of_formation

each isa ArrayRef[Num] that is lazy with public ARRAY traits: push_$_ get_$_
all_$_ clear_$_ count_$_
 
=attr U H G S S_t S_r S_v

each isa ArrayRef[Num] that is lazy with public ARRAY traits: push_$_ get_$_
all_$_ clear_$_ count_$_
              
=attr total_energy_mp2 total_energy_ccsdt nonelectrostatic_energy
 
each isa ArrayRef[Num] that is lazy with public ARRAY traits: push_$_ get_$_
all_$_ clear_$_ count_$_

=attr qm_dipole frequencies eigvec alpha beta

each isa ArrayRef[Math::Vector::Real] that is lazy with public ARRAY traits: push_$_ get_$_
all_$_ clear_$_ count_$_


=head1 SEE ALSO

=for :list
* L<HackaMol::Molecule>
* L<EMSL | https://bse.pnl.gov/bse/portal>
