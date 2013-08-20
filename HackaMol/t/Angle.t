use Modern::Perl;
use Test::Most;
use Test::Warnings;
use Test::Moose;
use Test::More;
use Math::Vector::Real;
use lib 'lib/HackaMol';
use Time::HiRes qw(time);
use Atom;
use Angle;


my @attributes = qw(
atoms ang_eq ang_fc  
);
my @methods = qw(
ang ang_normvec angle_energy clear_ang_eq clear_ang_fc has_ang_eq has_ang_fc
);

my @roles = qw(AtomsGroupRole);

map has_attribute_ok( 'Angle', $_ ), @attributes;
map can_ok ( 'Angle', $_ ), @methods;
map does_ok( 'Angle', $_ ), @roles;

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
    Z       => 80
);
my $atom2 = Atom->new(
    name    => 'C1',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    coords  => [ 
                V( 1.0, 0.0, 0.0 ), 
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
    Z => 6,
);


my $atom3 = Atom->new(
    name    => 'C3',
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

my $atom4 = Atom->new(
    name    => 'C2',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    coords  => [
                V(  0.0, 0.0, 1.0 ),
                V(  0.0, 1.0, 1.0 ),
                V(  0.0, 2.0, 2.0 ),
                V(  0.0, 3.0, 3.0 ),
                V(  0.0, 4.0, 4.0 ),
                V(  0.0, 5.0, 5.0 ),
                V(  0.0, 6.0, 6.0 ),
                V(  0.0, 7.0, 7.0 ),
                V(  0.0, 8.0, 8.0 ),
                V(  0.0, 9.0, 9.0 ),
               ],
    Z => 6,
);

my $atom5 = Atom->new(
    name    => 'C3',
    charges => [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1],
    coords  => [
                V(  0.0, 0.0, 0.0 ),
                V( -1.0, 1.0, 1.0 ),
                V( -2.0, 2.0, 2.0 ),
                V( -3.0, 3.0, 3.0 ),
                V( -4.0, 4.0, 4.0 ),
                V( -5.0, 5.0, 5.0 ),
                V( -6.0, 6.0, 6.0 ),
                V( -7.0, 7.0, 7.0 ),
                V( -8.0, 8.0, 8.0 ),
                V( -9.0, 9.0, 9.0 ),
               ],
    Z => 6,
);

my $angle1 = Angle->new(atoms => [$atom2,$atom1,$atom3]);
my $angle2 = Angle->new(atoms => [$atom2,$atom1,$atom4]);
my $angle3 = Angle->new(atoms => [$atom2,$atom1,$atom5]);

foreach my $t (0 .. 9){
  $angle1->do_forall('t',$t);
  $angle2->do_forall('t',$t);
  cmp_ok($angle1->ang,'==', 180.0, "antiparallel t dependent angle: 180");
  cmp_ok($angle2->ang,'==', 90.0, "xz t dependent ang: 90");
  is_deeply($angle1->ang_normvec, V(0,0,0), "antiparallel t dependent ang_normvec: V (0, 0, 0)");
  is_deeply($angle1->COM, $atom1->get_coords($t), "antiparallel COM at Hg");
  is_deeply($angle2->ang_normvec, V(0,-1,0), "xz t dependent ang_normvec: V (0, 1, 0)");
}

is($atom1->count_angles, 3, "atom1 knows it is in 3 angles");
is($atom2->count_angles, 3, "atom2 knows it is in 3 angles");
is($atom3->count_angles, 1, "atom3 knows it is in 1 angle");
is($atom4->count_angles, 1, "atom4 knows it is in 1 angle");
is($atom5->count_angles, 1, "atom5 knows it is in 1 angle");
is($atom1->get_angles(0),$angle1, 'the atom1 is aware of angle1');
is($atom1->get_angles(1),$angle2, 'the atom1 is aware of angle2');
is($atom1->get_angles(2),$angle3, 'the atom1 is aware of angle3');
is($atom2->get_angles(0),$angle1, 'the atom2 is aware of angle1');
is($atom2->get_angles(1),$angle2, 'the atom2 is aware of angle2');
is($atom2->get_angles(2),$angle3, 'the atom2 is aware of angle3');
is($atom3->get_angles(0),$angle1, 'the atom1 is aware of angle1');
is($atom4->get_angles(0),$angle2, 'the atom1 is aware of angle2');
is($atom5->get_angles(0),$angle3, 'the atom1 is aware of angle3');

#my $t1 = time;
#$atom1->distance($atom2) foreach 0 .. 10000;  
#my $t2 = time;
#$angle->angle_length foreach 0 .. 10000;  
#my $t3 = time;
#printf ("Calculate distance>   %10.3f per s \n", 10000/($t2-$t1));
#printf ("retrieve angle_length> %10.3f per s \n", 10000/($t3-$t2));
#print $angle->dump;

done_testing();
