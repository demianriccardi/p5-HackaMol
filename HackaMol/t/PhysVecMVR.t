{
  # see node_id=1049328 on perlmonks. This hack allows the required subroutine
  # to be included in the ClassCompositor
  package MooseX::ClassCompositor::ReqRole;
  use Moose;
  extends qw( MooseX::ClassCompositor );
  around class_for => sub {
    my $orig = shift;
    my $self = shift;
    my @roles = map {
      ref($_) eq q(HASH) ? 'Moose::Meta::Role'->create_anon_role(methods => $_)
: $_
    } @_;
    $self->$orig(@roles);
  };
}

use Test::Most;
use Test::Warnings;
use Test::Moose;
#use MooseX::ClassCompositor;    #use this for testing roles
use lib 'lib/roles', 't/lib';
use PhysVecMVR;                # v0.001;#To test for version availability
use Math::Vector::Real;           

my @attributes = qw( name t mass xyzfree is_fixed);
my @methods = qw(
  push_charges get_charges set_charges all_charges clear_charges
  push_coords get_coords set_coords all_coords clear_coords
  push_forces get_forces set_forces all_forces clear_forces
  distance 
  intra_dcoords intra_dforces intra_dcharges
  inter_dcoords inter_dforces inter_dcharges
  mean_coords mean_forces mean_charges
  msd_coords msd_forces msd_charges
  charge xyz force
);
my ( $obj1, $obj2, $obj3 );

my %methods = ('_build_mass' => sub{return 0}); # this is from the above

my $class = MooseX::ClassCompositor::ReqRole->new( { 
                                            class_basename => 'Test', 
                                          } )->class_for('PhysVecMVR',\%methods);

map has_attribute_ok( $class, $_ ), @attributes;
map can_ok( $class, $_ ), @methods;

lives_ok {
    $obj1 = $class->new( name => 'somephysvec', t => 1 );
}
'Test creation of an obj1';

$obj1->push_charges($_) foreach ( 0.3, 0.2, 0.1, -0.1, -0.2, -0.36 );
cmp_ok( $obj1->count_charges, '==', 6, 'pushed 6 _tcharges and we have 6' );

my $sum_charges = 0;
$sum_charges += $_ foreach $obj1->all_charges;
cmp_ok( $sum_charges, '==', -0.06, 'sum of charges as expected' );

cmp_ok( $obj1->get_charges(4), '==', -0.2, '5th _tcharges as expected' );

$obj1->set_charges( 4, 1.0 );
cmp_ok( $obj1->get_charges(4),
    '==', 1.0, '5th _tcharges set to 1.0 and get_charges as expected' );

lives_ok {
    $obj2 = $class->new( name => 'someotherphysvec', t => 0 );
}
'Test creation of an obj2';

$obj1->set_coords( $obj1->t, V( 1, 0, 0 ));
$obj2->set_coords( $obj2->t, V( 0, 1, 0 ) );
warning_is { $obj1->distance($obj2) }
"comparing objects with different times",
  "carp warning about obj1 and obj2 distance with different t";

$obj1->t(0);
$obj1->set_coords( $obj1->t, V( 1, 0, 0 ) );

my $sqrt_2 =
  sprintf( "%.16f", 1.4142135623730951 );    # double precision from fortran 90

cmp_ok( sprintf( "%.16f", $obj1->distance($obj2) ),
    '==',
    $sqrt_2, "distance between 100 and 010 to 16 decimal places: $sqrt_2" );

cmp_ok( sprintf( "%.16f", $obj2->distance($obj1) ),
    '==',
    $sqrt_2,
    "default obj2 distance between 100 and 010 to 16 decimal places: $sqrt_2" );

my $dq = $obj1->get_charges(1) - $obj1->get_charges(0);

cmp_ok( $obj1->intra_dcharges( 0, 1 ),
    '==', $dq, "change in charges between t = 0 and t = 1: $dq" );

$obj1->set_coords( 0, V( 0.0011,     -0.98458734, 1.0003984 ) );
$obj1->set_coords( 1, V( 1.00130011, 1.1,         2.0342 ) );
my $vec0 = $obj1->get_coords(0);
my $vec1 = $obj1->get_coords(1);
my @vec  = map { $vec1->[$_] - $vec0->[$_] } 0 .. 2;

is_deeply( $obj1->intra_dcoords( 0, 1 ), V(@vec), 'intra_dcoords' );

$obj1->push_forces( V( 0.0011,     -0.98458734, 1.0003984 ) );
$obj1->push_forces( V( 1.00130011, 1.1,         2.0342 ) );

$vec0 = $obj1->get_coords(0);
$vec1 = $obj1->get_coords(1);
@vec  = map { $vec1->[$_] - $vec0->[$_] } 0 .. 2;

is_deeply( $obj1->intra_dforces( 0, 1 ), V(@vec), 'intra_dforces' );

lives_ok {
    $obj3 = $class->new(
        name    => 'somephysvec3',
        t       => 0,
        charges => [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ],
        coords  => [
           V( 0, 0, 0 ),
           V( 1, 1, 1 ),
           V( 2, 2, 2 ),
           V( 3, 3, 3 ),
           V( 4, 4, 4 ),
           V( 5, 5, 5 ),
           V( 6, 6, 6 ),
           V( 7, 7, 7 ),
           V( 8, 8, 8 ),
           V( 9, 9, 9 ),
        ],
        forces => [
           V( 0, 0, 0 ),
           V( 1, 1, 1 ),
           V( 2, 2, 2 ),
           V( 3, 3, 3 ),
           V( 4, 4, 4 ),
           V( 5, 5, 5 ),
           V( 6, 6, 6 ),
           V( 7, 7, 7 ),
           V( 8, 8, 8 ),
           V( 9, 9, 9 ),
        ]
    );
}
'Test creation of an obj3';

is_deeply( $obj3->intra_dcoords( 0, 1 ), V( 1, 1, 1 ), 'intra_dcoords' );
is_deeply( $obj3->intra_dcharges( 0, 9 ), 9, 'intra_dcharges' );

for my $t ( 0 .. $obj3->count_coords - 1 ) {
    $obj3->t($t);
    is(
        $obj3->charge,
        $obj3->get_charges($t),
        "obj->charge returns get_charges($t)"
    );
    is_deeply(
        $obj3->xyz,
        $obj3->get_coords($t),
        "obj->xyz    returns get_coords($t)"
    );
    is_deeply(
        $obj3->force,
        $obj3->get_forces($t),
        "obj->force  returns get_forces($t)"
    );
}

is_deeply($obj3->mean_coords, V(4.5,4.5,4.5), "average coordinates");
is($obj3->msd_coords, 24.75, "mean square deviation forces");
is_deeply($obj3->mean_forces, V(4.5,4.5,4.5), "average coordinates");
is($obj3->msd_forces, 24.75, "mean square deviation forces");
is($obj3->mean_charges, 4.5, "average charges");
is($obj3->msd_charges, 8.25, "mean square deviation charges");

done_testing();

