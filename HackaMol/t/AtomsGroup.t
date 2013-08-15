use Test::Most;
use Test::Warnings;
use Test::Moose;
use MooseX::ClassCompositor;    #use this for testing roles
use lib 'lib/roles','lib/HackaMol';
use Atom;
use AtomsGroup;                # v0.001;#To test for version availability

my @attributes = qw(
atoms dipole COM COZ dipole_moment total_charge atoms_bin
);
my @methods = qw(
_build_dipole _build_COM _build_COZ _build_total_charge _build_dipole_moment 
_build_atoms_bin clear_atoms_bin clear_dipole clear_dipole_moment clear_COM 
clear_COZ clear_total_charge clear_dipole_moment
set_atoms_bin get_atoms_bin has_empty_bin count_unique_atoms all_unique_atoms 
atom_counts canonical_name Rg all_atoms push_atoms get_atoms delete_atoms count_atoms
clear_atoms
);

my $class = MooseX::ClassCompositor->new( { 
                                            class_basename => 'Test', 
                                          } )->class_for('AtomsGroup');

map has_attribute_ok( $class, $_ ), @attributes;
map can_ok( $class, $_ ), @methods;
my $obj;
lives_ok {
    $obj = $class->new();
}
'Test creation of an obj';

my $atom1 = Atom->new(
    name    => 'H',
    charges => [−0.80,−0.82,-0.834],
    coords  => [ V(2.05274,        0.01959,       -0.07701) ],
    Z       => 1
);

my $atom2 = Atom->new(
    name    => 'H',
    charges => [0.4,0.41,0.417],
    coords  => [ V( 1.08388,        0.02164,       -0.12303 ) ],
    Z       => 1
);
my $atom3 = Atom->new(
    name    => 'H',
    charges => [0.4,0.41,0.417],
    coords  => [ V( 2.33092,        0.06098,       -1.00332 ) ],
    Z       => 1
);

$obj->push_atoms($_) foreach ($atom1, $atom2, $atom3);


done_testing();
