use Modern::Perl;
use Test::Most;
use Test::Warnings;
use Test::Moose;
use Test::More;
use Math::Vector::Real;
use lib 'lib/HackaMol';
use Time::HiRes qw(time);
use Atom;
use Bond;


my @attributes = qw(
atoms bond_length bond_vector bond_order
);
my @methods = qw(
_build_bond_length _build_bond_vector 
clear_bond_length clear_bond_vector clear_bond_order
);

my @roles = qw(AtomsGroupRole);

map has_attribute_ok( 'Bond', $_ ), @attributes;
map can_ok( 'Bond', $_ ), @methods;
map does_ok( 'Bond', $_ ), @roles;

my $atom1 = Atom->new(
    name    => 'Hg',
    charges => [2,2,2,2,2,2,2,2,2,2],
    coords  => [ 
                V( 0.0, 0.0, 0.0 ), 
                V( 0.0, 1.0, 0.0 ), 
                V( 0.0, 2.0, 0.0 ), 
                V( 0.0, 3.0, 0.0 ), 
                V( 0.0, 4.0, 0.0 ), 
                V( 0.0, 5.0, 0.0 ), 
                V( 0.0, 6.0, 0.0 ), 
                V( 0.0, 7.0, 0.0 ), 
                V( 0.0, 8.0, 0.0 ), 
                V( 0.0, 9.0, 0.0 ), 
               ],
    symbol  => 'HG'
);

my $atom2 = Atom->new(
    name    => 'C1',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    coords  => [ 
                V( 0.0, 0.0, 0.0 ), 
                V( 1.0, 1.0, 0.0 ), 
                V( 2.0, 2.0, 0.0 ), 
                V( 3.0, 3.0, 0.0 ), 
                V( 4.0, 4.0, 0.0 ), 
                V( 5.0, 5.0, 0.0 ), 
                V( 6.0, 6.0, 0.0 ), 
                V( 7.0, 7.0, 0.0 ), 
                V( 8.0, 8.0, 0.0 ), 
                V( 9.0, 9.0, 0.0 ), 
               ],
    Z       => 6
);

my $atom3 = Atom->new(
    name    => 'C2',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    coords  => [
                V( -1.0, 0.0, 0.0 ),
                V( -1.0, 1.0, 0.0 ),
                V( -2.0, 2.0, 0.0 ),
                V( -3.0, 3.0, 0.0 ),
                V( -4.0, 4.0, 0.0 ),
                V( -5.0, 5.0, 0.0 ),
                V( -6.0, 6.0, 0.0 ),
                V( -7.0, 7.0, 0.0 ),
                V( -8.0, 8.0, 0.0 ),
                V( -9.0, 9.0, 0.0 ),
               ],
    Z => 6,
);


my $bond1 = Bond->new(atoms => [$atom1,$atom2]);
my $bond2 = Bond->new(atoms => [$atom1,$atom3]);

foreach my $t (0 .. 9){
  $bond1->gt($t);
  cmp_ok($bond1->bond_length,'==', $t, "t dependent bond length: $t");
  is_deeply($bond1->bond_vector, V($t,0,0), "t dependent bond vector: V ($t, 0, 0)");
  is($bond1->bond_order, 1, "bond order default");
  $bond1->bond_order(1.5);
}

is($bond1->bond_order, 1.5, "bond order set to num");
is($atom1->count_bonds, 2, "atom1 knows it has 2 bonds");
is($atom2->count_bonds, 1, "atom2 knows it has 1 bonds");
is($atom3->count_bonds, 1, "atom3 knows it has 1 bonds");
is($atom1->get_bonds(0),$bond1, 'the atom is aware of its bond');
is($atom2->get_bonds(0),$bond1, 'the atom is aware of its bond');
is($atom1->get_bonds(1),$bond2, 'the atom is aware of its bond');
is($atom3->get_bonds(0),$bond2, 'the atom is aware of its other bond');

print "use Molecule class to really delete bonds\n";
#$bond1 = undef; # = {};
#print $atom2->get_bonds(0)->dump;
#is($atom1->get_bonds(0),$bond1, 'the atom is aware of its bond');


#my $t1 = time;
#$atom1->distance($atom2) foreach 0 .. 10000;  
#my $t2 = time;
#$bond->bond_length foreach 0 .. 10000;  
#my $t3 = time;
#printf ("Calculate distance>   %10.3f per s \n", 10000/($t2-$t1));
#printf ("retrieve bond_length> %10.3f per s \n", 10000/($t3-$t2));
#print $bond->dump;

done_testing();

