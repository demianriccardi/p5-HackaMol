use Test::Most;
use Test::Warnings;
use Test::Moose;
use lib 'lib/HackaMol','t/lib';
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;
use AtomGroup;
use Bond;
use Molecule;
use PDBintoAtoms qw(readinto_atoms);

my @attributes = qw(
atoms mass name 
);
my @methods = qw(
push_groups set_groups get_groups all_groups clear_groups
delete_groups count_groups 
Rg 
);
my @roles = qw(PhysVecMVRRole BondsAnglesDihedralsRole AtomGroupRole);

map has_attribute_ok( 'Molecule', $_ ), @attributes;
map can_ok (          'Molecule', $_ ), @methods;
map does_ok(          'Molecule', $_ ), @roles;

my @atoms = readinto_atoms("t/lib/2LL5.pdb");
my $max_t = $atoms[0]->count_coords -1;

my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

is($mol->count_atoms, 283, 'number of atoms: 283');

#dubious test, Rgs COMs generated from HackaMol, double check
my @Rgt = qw(7.19251198554957 7.15700534761354 7.14685267470545 7.00573838170868 7.21240816173855 7.18269096506165 7.090359584332 7.15239515247956 7.16336404927373 7.21275609325176 7.19505186699334 7.10928755300199 7.11941832497372 7.21511627204864 7.3088703725112 7.19876994990288 7.25732589008288 7.23532446379484 7.1890879967172 7.13373953831758 7.12469776280913 6.95147600134586 7.09343947318801 7.25038296705691 7.16489832835357 7.0962667805229 7.05103784479553 7.04501594915213 7.16973657632147 7.12776512790071 7.19457930266239 7.1262398705025 7.21766853448075);

my @COMs = (
 [29.132,   0.299,   0.412],
 [29.135,   0.259,   0.422],
 [29.236,   0.088,   0.452],
 [29.137,   0.298,   0.481],
 [29.122,   0.284,   0.504],
 [29.085,   0.210,   0.419],
 [29.081,   0.316,   0.455],
 [29.116,   0.176,   0.385],
 [29.111,   0.315,   0.454],
 [29.175,   0.168,   0.443],
 [29.117,   0.280,   0.393],
 [29.243,   0.256,   0.432],
 [29.082,   0.300,   0.403],
 [29.117,   0.168,   0.447],
 [29.169,   0.269,   0.477],
 [29.256,   0.250,   0.503],
 [29.181,   0.216,   0.445],
 [29.119,   0.180,   0.473],
 [29.126,   0.186,   0.459],
 [29.137,   0.287,   0.434],
 [29.171,   0.228,   0.517],
 [29.078,   0.128,   0.481],
 [29.222,   0.125,   0.476],
 [29.186,   0.108,   0.372],
 [29.084,   0.345,   0.415],
 [29.077,   0.162,   0.464],
 [29.240,   0.254,   0.445],
 [29.149,   0.204,   0.452],
 [29.079,   0.203,   0.457],
 [29.044,   0.193,   0.350],
 [29.074,   0.256,   0.482],
 [29.134,   0.125,   0.523],
 [29.101,   0.268,   0.427],
);

foreach my $t (0 .. $max_t) {
  $mol->t($t);
  my $ocom = sprintf ("%8.3f %8.3f %8.3f\n",@{$mol->COM});
  my $ecom = sprintf ("%8.3f %8.3f %8.3f\n",@{$COMs[$t]});
  is($ocom,$ecom, "COM at $t");
  cmp_ok(abs($mol->Rg-$Rgt[$t]) ,'<', 1E-7, "Rg at $t");
}

$mol->push_groups_by_atom_attr('resid'); 

is($mol->count_groups, 22, "group_by_atom_resid yields 22 groups");
$mol->clear_groups;
is($mol->count_groups, 0, "clear->groups yields 0 groups");

$mol->push_groups_by_atom_attr('symbol'); 

is($mol->count_groups, 4, "group_by_atom_symbol yields 4 (ONCH) groups");

$mol->clear_groups;
$mol->push_groups_by_atom_attr('name'); 

is($mol->count_groups, 60, "group_by_atom_name yields 60 groups");

my @bonds = map{Bond->new(atoms=>[$atoms[0],$_])} @atoms;  

$mol->push_bonds(@bonds);

my @ncord = grep{$_->bond_length < 5.0 and $_->bond_length > 0.0} $mol->all_bonds;
foreach my $bond (@ncord){
  is($bond->get_atoms(0)->t, $mol->t, "atoms in bond time and mol time are same");
}

is(scalar(@ncord), 23, "23 atoms within 5 angstroms of first atom");

use Math::Trig;

foreach my $t (1 .. 360) {
  #$mol->t($t);

  print $mol->count_groups . "\n\n";
  
  my $v = V(1,0,0);
  my @gs = map{$_->COM} $mol->all_groups;
  my @s = $v->rotate_3d($t/(2*pi), @gs);  
  printf("Hg %8.3f %8.3f %8.3f\n", @{$_}) foreach @s;

#  foreach my $g ($mol->all_groups){
#    printf("Hg %8.3f %8.3f %8.3f\n", @{$g->COM});
#    printf("Zn %8.3f %8.3f %8.3f\n", @{$g->COZ});
#  }
}

done_testing();



