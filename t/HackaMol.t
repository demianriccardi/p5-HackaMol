use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Fatal qw(dies_ok);
use Test::Warn;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;
use HackaMol;

my @attributes = qw( 
                    name 
                   );
my @methods = qw(
  build_bonds build_angles build_dihedrals 
  group_by_atom_attr read_file_atoms read_pdb_atoms 
  read_xyz_atoms
);

my @roles = qw(HackaMol::MolReadRole HackaMol::NameRole);

map has_attribute_ok( 'HackaMol', $_ ), @attributes;
map can_ok( 'HackaMol', $_ ), @methods;
map does_ok( 'HackaMol', $_ ), @roles;


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

#test hackamol class

my $hack = HackaMol->new( name => "hackitup" );
is ($hack->name, "hackitup", "HackaMol name attr");

my @atoms1 = $hack->read_file_atoms("t/lib/1L2Y.pdb");
my $mol1 = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms1] );
is ($mol1->count_atoms, 304 , "read atoms in from pdb");

unlink("t/lib/1L2Y.xyz");
my $fh = $mol1->print_xyz( "t/lib/1L2Y.xyz" );
$fh->close;
#foreach my $t ( 1 .. $atoms1[0]->count_coords - 1 ) {
#       $mol1->t($t);
#       $mol1->print_xyz($fh);
#}

my @atoms2 = $hack->read_file_atoms("t/lib/1L2Y.xyz");
my $mol2 = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms2] );

is ($mol2->count_atoms, 304 , "read atoms in from xyz");

my @Z1 = map {$_->Z} $mol1->all_atoms;
my @Z2 = map {$_->Z} $mol2->all_atoms;

is_deeply(\@Z1,\@Z2, "xyz and pdb give same atoms");

dies_ok {$hack->read_file_atoms("bah.mol")} "Croak on unsupported file type";

dies_ok {$hack->read_file_atoms("t/lib/bad1.xyz")} "xyz Croak change symbol";
dies_ok {$hack->read_file_atoms("t/lib/bad2.xyz")} "xyz Croak change Z";
dies_ok {$hack->read_file_atoms("t/lib/bad3.xyz")} "xyz Croak change number of atoms";
dies_ok {$hack->read_file_atoms("t/lib/bad1.pdb")}  "pdb Croak change atom name";
dies_ok {$hack->read_file_atoms("t/lib/bad2.pdb")}  "pdb Croak change element";

my @wats1 = $hack->read_file_atoms("t/lib/byZ.xyz");  
my @wats2 = $hack->read_file_atoms("t/lib/byZSym.xyz");  

my @Zw1 = map {$_->Z} @wats1;
my @Zw2 = map {$_->Z} @wats2;

is_deeply(\@Z1,\@Z2, "different Z/Symbol formatted xyz give same");

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

