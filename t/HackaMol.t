#!/usr/bin/env perl

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
cmp_ok( abs( $group->COZ - $group->COM ), '<',  1E-7, 'COM ~ COZ' );
cmp_ok( $group->total_charge,             '==', 0,    'group total charges 0' );
cmp_ok( $mol->total_charge,               '==', 0,    'mol total charges 0' );
cmp_ok( abs($group->dipole_moment), '<', 1E-7,
    'group dipole moment is zero, no charges' );
cmp_ok( abs($mol->dipole_moment), '<', 1E-7, 'mol dipole moment is zero, no charges' );
my $exp_Rg = sqrt( $radius * $radius * 3 / 5 );
cmp_ok( abs( $exp_Rg - $group->Rg ),
    '<', 0.75, 'group numerical Rg within 0.75 Angs of theoretical' );
cmp_ok( abs( $mol->Rg - $group->Rg ), '<', 1E-7, 'group and Mol Rg same' );

#test hackamol class

my $hack = HackaMol->new( name => "hackitup" );
is( $hack->name, "hackitup", "HackaMol name attr" );

my @atoms1 = $hack->read_file_atoms("t/lib/1L2Y.pdb");
my $mol1 = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms1] );
is( $mol1->count_atoms, 304, "read atoms in from pdb" );

unlink("t/lib/1L2Y.xyz");
my $fh = $mol1->print_xyz("t/lib/1L2Y.xyz");
$fh->close;

#foreach my $t ( 1 .. $atoms1[0]->count_coords - 1 ) {
#       $mol1->t($t);
#       $mol1->print_xyz($fh);
#}

my @atoms2 = $hack->read_file_atoms("t/lib/1L2Y.xyz");
my $mol2 = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms2] );

is( $mol2->count_atoms, 304, "read atoms in from xyz" );

my @Z1 = map { $_->Z } $mol1->all_atoms;
my @Z2 = map { $_->Z } $mol2->all_atoms;

is_deeply( \@Z1, \@Z2, "xyz and pdb give same atoms" );

dies_ok { $hack->read_file_atoms("bah.mol") } "Croak on unsupported file type";

dies_ok { $hack->read_file_atoms("t/lib/bad1.xyz") } "xyz Croak change symbol";
dies_ok { $hack->read_file_atoms("t/lib/bad2.xyz") } "xyz Croak change Z";
dies_ok { $hack->read_file_atoms("t/lib/bad3.xyz") }
"xyz Croak change number of atoms";
dies_ok { $hack->read_file_atoms("t/lib/bad1.pdb") }
"pdb Croak change atom name";
dies_ok { $hack->read_file_atoms("t/lib/bad2.pdb") } "pdb Croak change element";

my @wats1 = $hack->read_file_atoms("t/lib/byZ.xyz");
my @wats2 = $hack->read_file_atoms("t/lib/byZSym.xyz");

my @Zw1 = map { $_->Z } @wats1;
my @Zw2 = map { $_->Z } @wats2;

is_deeply( \@Z1, \@Z2, "different Z/Symbol formatted xyz give same" );

#test group generation
my @gresids = $hack->group_by_atom_attr( 'resid',  $mol1->all_atoms );
my @gsymbls = $hack->group_by_atom_attr( 'symbol', $mol1->all_atoms );
my @gnames  = $hack->group_by_atom_attr( 'name',   $mol1->all_atoms );
$mol1->push_groups(@gresids);
is( $mol1->count_groups, 20, "group_by_atom_resid" );
$mol1->clear_groups;
is( $mol1->count_groups, 0, "clear->groups" );
$mol1->push_groups(@gsymbls);
is( $mol1->count_groups, 4, "group_by_atom_symbol" );
$mol1->clear_groups;
$mol1->push_groups(@gnames);
is( $mol1->count_groups, 74, "group_by_atom_name" );

my @bb = grep { $_->name eq 'N' or $_->name eq 'CA' or $_->name eq 'C' }
  $mol1->all_atoms;

my @bonds     = $hack->build_bonds(@bb);
my @angles    = $hack->build_angles(@bb);
my @dihedrals = $hack->build_dihedrals(@bb);

is( scalar(@bonds),       scalar(@bb) - 1,    "number of bonds generated" );
is( scalar(@angles),      scalar(@bb) - 2,    "number of angles generated" );
is( scalar(@dihedrals),   scalar(@bb) - 3,    "number of dihedrals generated" );
is( $dihedrals[0]->name,  'N1_CA1_C1_N2',     "dihedral name1" );
is( $angles[0]->name,     'N1_CA1_C1',        "angle name1" );
is( $bonds[0]->name,      'N1_CA1',           "bond name1" );
is( $bonds[1]->name,      'CA1_C1',           "bond name2" );
is( $bonds[2]->name,      'C1_N2',            "bond name3" );
is( $angles[1]->name,     'CA1_C1_N2',        "second angle name" );
is( $dihedrals[-1]->name, 'C19_N20_CA20_C20', "last dihedral name" );
is( $bonds[-1]->name,     'CA20_C20',         "last bond name1" );
is( $angles[-1]->name,    'N20_CA20_C20',     "last angle name1" );

dies_ok { $hack->build_dihedrals(@bb[0 .. 2]) } "build_dihedrals croak";
dies_ok { $hack->build_bonds($bb[0]) }          "build_bonds croak";
dies_ok { $hack->build_angles(@bb[0,1]) }    "build_angles croak";


done_testing();

