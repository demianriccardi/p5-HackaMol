use Modern::Perl;
use Test::Most;
use Test::Warnings;
use Test::Moose;
use Test::More;
use Math::Vector::Real;
use lib 'lib/HackaMol';
use Time::HiRes qw(time);
use Atom;
use Dihedral;


my @attributes = qw(
atoms dihe_eq dihe_multi dihe_dphase dihe_fc 
);
my @methods = qw(
torsion_energy improper_dihe_energy dihe
);

my @roles = qw(AtomsGroupRole);

map has_attribute_ok( 'Dihedral', $_ ), @attributes;
map can_ok( 'Dihedral', $_ ), @methods;
map does_ok( 'Dihedral', $_ ), @roles;

my $atom0 = Atom->new(
    name    => 'C1',
    charges => [0],
    coords  => [
                V( 1.0, 1.0, 0.0 ),
               ],
    Z => 6,
);
my $atom1 = Atom->new(
    name    => 'S1',
    charges => [0],
    coords  => [
                V( 1.0, 0.0, 0.0 ),
               ],
    Z => 16,
);
my $atom2 = Atom->new(
    name    => 'S2',
    charges => [0],
    coords  => [ 
                V( -1.0, 0.0, 0.0 ), 
               ],
    Z       => 16
);
my $atom3 = Atom->new(
    name    => 'C2',
    charges => [0],
    coords  => [
                V( -1.0, -1.0, 0.0 ),
               ],
    Z => 6,
);

my $dihe= Dihedral->new(atoms => [$atom0,$atom1,$atom2,$atom3]);

cmp_ok($dihe->dihe,'==', 180.0, "180 dihedral");
is_deeply($dihe->COM, V(0,0,0), "COM at 0,0,0");

is($atom0->count_dihedrals, 1, "atom0 knows it is in 1 dihedrals");
is($atom1->count_dihedrals, 1, "atom1 knows it is in 1 dihedrals");
is($atom2->count_dihedrals, 1, "atom2 knows it is in 1 dihedrals");
is($atom3->count_dihedrals, 1, "atom3 knows it is in 1 dihedrals");
is($atom0->get_dihedrals(0),$dihe, 'the atom1 is aware of dihedral');
is($atom1->get_dihedrals(0),$dihe, 'the atom1 is aware of dihedral');
is($atom2->get_dihedrals(0),$dihe, 'the atom1 is aware of dihedral');
is($atom3->get_dihedrals(0),$dihe, 'the atom1 is aware of dihedral');

$atom3->set_coords(0,V(-1.0,sqrt(2)/2,sqrt(2)/2));
cmp_ok($dihe->dihe,'==', -45.0, "-45 dihedral");

$atom3->set_coords(0,V(-1.0,-sqrt(2)/2,-sqrt(2)/2));
cmp_ok($dihe->dihe,'==', 135.0, "135 dihedral");


done_testing();
