use Modern::Perl;
use Test::Most;
use Test::Warnings;
use Test::Moose;
use Test::More;
use Math::Vector::Real;
use Time::HiRes qw(time);
use HackaMol::Atom;
use HackaMol::Angle;

my @attributes = qw(
atoms ang_eq ang_fc  
);
my @methods = qw(
ang ang_normvec angle_energy clear_ang_eq clear_ang_fc has_ang_eq has_ang_fc
);

my @roles = qw(AtomGroupRole);

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

$angle1->ang_fc(1.0);
$angle1->ang_eq($angle1->ang - 0.5);

cmp_ok (abs(0.25-$angle1->angle_energy),'<',1E-7, 'simple angle energy test') ;

$angle1->angle_energy_func(
                          sub {
                               my $a = shift;
                               my $sum = 0;
                               $sum += $_*$a->ang foreach (@_);
                               return($sum);
                              }
                        );

cmp_ok (
        abs($angle1->angle_energy(1,2,3,4) - 10*$angle1->ang),
        '<', 1E-7, 'new nonsense energy'
       );

done_testing();
