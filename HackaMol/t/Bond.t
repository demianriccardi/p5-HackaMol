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
atoms dipole COM COZ dipole_moment total_charge atoms_bin
);
my @methods = qw(
_build_dipole _build_COM _build_COZ _build_total_charge _build_dipole_moment 
clear_atoms_bin clear_dipole clear_dipole_moment clear_COM 
clear_COZ clear_total_charge clear_dipole_moment
set_atoms_bin get_atoms_bin has_empty_bin count_unique_atoms all_unique_atoms 
atom_counts canonical_name Rg all_atoms push_atoms get_atoms delete_atoms count_atoms
clear_atoms
);

my @roles = qw(AtomsGroup);

map has_attribute_ok( 'Bond', $_ ), @attributes;
map can_ok( 'Bond', $_ ), @methods;
map does_ok( 'Bond', $_ ), @roles;

my $atom1 = Atom->new(
    name    => 'C',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
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
    Z       => 6
);
my $atom2 = Atom->new(
    name    => 'Hg',
    charges => [2,2,2,2,2,2,2,2,2,2],
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
    symbol  => 'HG'
);

my $bond = Bond->new(atoms => [$atom1,$atom2]);

foreach my $t (0 .. 9){
  $bond->t($t);
  cmp_ok($bond->bond_length,'==', $t, "t dependent bond length: $t");
  is_deeply($bond->bond_vector, V($t,0,0), "t dependent bond vector: V ($t, 0, 0)");
  is($bond->bond_order, 1, "bond order default");
  $bond->bond_order(1.5);
}

is($bond->bond_order, 1.5, "bond order set to num");

#my $t1 = time;
#$atom1->distance($atom2) foreach 0 .. 10000;  
#my $t2 = time;
#$bond->bond_length foreach 0 .. 10000;  
#my $t3 = time;
#printf ("Calculate distance>   %10.3f per s \n", 10000/($t2-$t1));
#printf ("retrieve bond_length> %10.3f per s \n", 10000/($t3-$t2));
#print $bond->dump;

done_testing();

