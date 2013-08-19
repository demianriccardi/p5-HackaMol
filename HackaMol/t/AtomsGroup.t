use Test::Most;
use Test::Warnings;
use Test::Moose;
use lib 'lib/HackaMol';
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;
use AtomsGroup;
use Atom;

my @attributes = qw(
atoms dipole COM COZ dipole_moment total_charge atoms_bin
);
my @methods = qw(
_build_dipole _build_COM _build_COZ _build_total_charge _build_dipole_moment 
clear_atoms_bin clear_dipole clear_dipole_moment clear_COM 
clear_COZ clear_total_charge clear_dipole_moment
set_atoms_bin get_atoms_bin has_empty_bin count_unique_atoms all_unique_atoms 
atom_counts canonical_name atoms all_atoms push_atoms get_atoms delete_atoms count_atoms clear_atoms Rg
);
my @roles = qw(AtomsGroupRole);

map has_attribute_ok( 'AtomsGroup', $_ ), @attributes;
map can_ok (          'AtomsGroup', $_ ), @methods;
map does_ok(          'AtomsGroup', $_ ), @roles;

my $radius = 16;
my $natoms = int(0.0334*($radius**3)*4*pi/3);

my @atoms = map {Atom->new(Z => 8, charges=> [0], coords => [$_]) } 
            map {$_*$radius} 
            map {Math::Vector::Real->random_in_sphere(3)} 1 .. $natoms;

my $group = AtomsGroup->new(gname => 'biggroup', atoms=> [@atoms]);

is($group->count_atoms, $natoms, "atom count: $natoms");
is($group->count_unique_atoms, 1, 'unique atoms in sphere is 1');
is($group->canonical_name, "O$natoms", "sphere atoms named O$natoms");
cmp_ok(1-abs($group->COM), '>',0, 'center of mass within 1 angstrom of 0,0,0');
my $exp_Rg = sqrt($radius*$radius*3/5);
cmp_ok(abs($exp_Rg-$group->Rg), '<',0.75, 'numerical Rg within 0.75 Angs of theoretical');

done_testing();



