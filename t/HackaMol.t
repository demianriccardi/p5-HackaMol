use strict;
use warnings;
use Test::More;
use Test::Warn;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;
use HackaMol;



my $merc = HackaMol::Atom->new(
    name    => "Mercury",
    Z       => 80,
    charges => [0],
    coords  => [ V( 0, 0, 0 ) ]
);
my $mercg = HackaMol::AtomGroup->new( name => "MercG", atoms => [$merc] );
my $mercm =
  HackaMol::Molecule->new( name => "MercM", atoms => [ $mercg->all_atoms ] );

ok( defined $merc,                      'HackaMol::Atom->new' );
ok( $merc->isa('HackaMol::Atom'),       'isa HackaMol::Atom' );
ok( $merc->symbol eq 'Hg',              'symbol Hg' );
ok( $merc->name eq 'Mercury',           'name is as set' );
ok( defined $mercg,                     'HackaMol::AtomGroup->new' );
ok( $mercg->isa('HackaMol::AtomGroup'), 'isa HackaMol::AtomGroup' );
ok( defined $mercm,                     'HackaMol::Molecule->new' );
ok( $mercm->isa('HackaMol::Molecule'),  'isa HackaMol::Molecule' );

#taken from AtomGroup Test
my $radius = 16;
my $natoms = int( 0.0334 * ( $radius**3 ) * 4 * pi / 3 );

my @atoms =
  map { HackaMol::Atom->new( Z => 8, charges => [0], coords => [$_] ) }
  map { Math::Vector::Real->random_in_sphere( 3, $radius ) } 1 .. $natoms;

my $group = HackaMol::AtomGroup->new( name => 'biggroup', atoms => [@atoms] );
my $mol = HackaMol::Molecule->new(
    name   => 'bg_mol',
    atoms  => [ $group->all_atoms ],
    groups => [$group]
);

is( $group->count_atoms,        $natoms, "group atom count: $natoms" );
is( $mol->count_atoms,          $natoms, "mol atom count: $natoms" );
is( $group->count_unique_atoms, 1,       'group unique atoms in sphere is 1' );
is( $mol->count_unique_atoms,   1,       "mol unique atoms in sphere is 1" );
is( $group->bin_atoms_name, "O$natoms", "group sphere atoms named O$natoms" );
is( $mol->bin_atoms_name,   "O$natoms", "mol sphere atoms named O$natoms" );
cmp_ok( 2 - abs( $group->COM ),
    '>', 0, 'group center of mass within 2 angstrom of 0,0,0' );
cmp_ok( abs( $mol->COM - $group->COM ),
    '<', 1E-10, 'mol com same as group com' );
cmp_ok( abs( $group->COZ - $group->COM ), '<',  1E-6, 'COM ~ COZ' );
cmp_ok( $group->total_charge,             '==', 0,    'group total charges 0' );
cmp_ok( $mol->total_charge,               '==', 0,    'mol total charges 0' );
cmp_ok( $group->dipole_moment, '==', 0,
    'group dipole moment is zero, no charges' );
cmp_ok( $mol->dipole_moment, '==', 0, 'mol dipole moment is zero, no charges' );
my $exp_Rg = sqrt( $radius * $radius * 3 / 5 );
cmp_ok( abs( $exp_Rg - $group->Rg ),
    '<', 0.75, 'group numerical Rg within 0.75 Angs of theoretical' );
cmp_ok( abs( $mol->Rg - $group->Rg ), '<', 1E-10, 'group and Mol Rg same' );


#$mol->push_groups_by_atom_attr('resid');
#is( $mol->count_groups, 22, "group_by_atom_resid yields 22 groups" );
#$mol->clear_groups;
#is( $mol->count_groups, 0, "clear->groups yields 0 groups" );
#$mol->push_groups_by_atom_attr('symbol');
#is( $mol->count_groups, 4, "group_by_atom_symbol yields 4 (ONCH) groups" );
#$mol->clear_groups;
#$mol->push_groups_by_atom_attr('name');
#is( $mol->count_groups, 60, "group_by_atom_name yields 60 groups" );

done_testing();

