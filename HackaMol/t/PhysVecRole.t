use Test::Most;
use Test::Warnings;
use Test::Moose;
use MooseX::ClassCompositor;    #use this for testing roles
use lib 'lib/roles', 't/lib';
use PhysVec;                # v0.001;#To test for version availability
use Math::VectorReal;           # just for example test of swapping coderefs

my @attributes = qw( name t mass xyzfree is_fixed
  distance_coderef inter_vector_coderef);
my @methods = qw(
  push_charges get_charges set_charges all_charges clear_charges
  push_coords get_coords set_coords all_coords clear_coords
  push_forces get_forces set_forces all_forces clear_forces
  distance delta_coords delta_forces delta_charges
);
my ( $obj1, $obj2, $obj3 );

my %methods = ('_build_mass' => sub{return 0});

my $class = MooseX::ClassCompositor->new( { 
                                            class_basename => 'Test', 
                                          } )->class_for('PhysVec');

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

$obj1->set_coords( $obj1->t, [ 1, 0, 0 ] );
$obj2->set_coords( $obj2->t, [ 0, 1, 0 ] );
warning_is { $obj1->distance($obj2) }
"you are comparing objects with different times",
  "carp warning about obj1 and obj2 distance with different t";

$obj1->t(0);
$obj1->set_coords( $obj1->t, [ 1, 0, 0 ] );

my $sqrt_2 =
  sprintf( "%.16f", 1.4142135623730951 );    # double precision from fortran 90

cmp_ok( sprintf( "%.16f", $obj1->distance($obj2) ),
    '==',
    $sqrt_2, "distance between 100 and 010 to 16 decimal places: $sqrt_2" );

$obj1->distance_coderef(
    sub {
        my $self  = shift;
        my $pvec  = shift;
        my $tself = $self->t;
        my $tpvec = $pvec->t;
        my $vec1  = $self->get_coords($tself);
        my $vec2  = $pvec->get_coords($tpvec);
        my $dist  = 0;
        $dist += ( $vec1->[$_] - $vec2->[$_] )**2 foreach 0 .. 2;

        #return (sqrt($dist));
        return ($dist);
    }
);

cmp_ok(
    sprintf( "%.16f", $obj1->distance($obj2) ),
    '==',
    sprintf( "%.16f", 2.0 ),
"new coderef obj1 distance**2 between 100 and 010 to 16 decimal places: 2.0000000000000000"
);

cmp_ok( sprintf( "%.16f", $obj2->distance($obj1) ),
    '==',
    $sqrt_2,
    "default obj2 distance between 100 and 010 to 16 decimal places: $sqrt_2" );

$obj1->distance_coderef(
    sub {
        my $self  = shift;
        my $pvec  = shift;
        my $tself = $self->t;
        my $tpvec = $pvec->t;
        my $vec1  = vector( @{ $self->get_coords($tself) } );
        my $vec2  = vector( @{ $pvec->get_coords($tpvec) } );
        my $dist  = ( $vec1 - $vec2 )->length;
        return ($dist);
    }
);
cmp_ok( sprintf( "%.16f", $obj1->distance($obj2) ), '==', $sqrt_2,
"swap VectorReal into obj1->distance_coderef between 100 and 010 to 16 decimal places: $sqrt_2"
);

my $dq = $obj1->get_charges(1) - $obj1->get_charges(0);

cmp_ok( $obj1->delta_charges( 0, 1 ),
    '==', $dq, "change in charges between t = 0 and t = 1: $dq" );

$obj1->set_coords( 0, [ 0.0011,     -0.98458734, 1.0003984 ] );
$obj1->set_coords( 1, [ 1.00130011, 1.1,         2.0342 ] );
my $vec0 = $obj1->get_coords(0);
my $vec1 = $obj1->get_coords(1);
my @vec  = map { $vec1->[$_] - $vec0->[$_] } 0 .. 2;

is_deeply( $obj1->delta_coords( 0, 1 ), \@vec, 'delta_coords' );

$obj1->push_forces( [ 0.0011,     -0.98458734, 1.0003984 ] );
$obj1->push_forces( [ 1.00130011, 1.1,         2.0342 ] );

$vec0 = $obj1->get_coords(0);
$vec1 = $obj1->get_coords(1);
@vec  = map { $vec1->[$_] - $vec0->[$_] } 0 .. 2;

is_deeply( $obj1->delta_forces( 0, 1 ), \@vec, 'delta_forces' );

lives_ok {
    $obj3 = $class->new(
        name    => 'somephysvec3',
        t       => 0,
        charges => [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ],
        coords  => [
            [ 0, 0, 0 ],
            [ 1, 1, 1 ],
            [ 2, 2, 2 ],
            [ 3, 3, 3 ],
            [ 4, 4, 4 ],
            [ 5, 5, 5 ],
            [ 6, 6, 6 ],
            [ 7, 7, 7 ],
            [ 8, 8, 8 ],
            [ 9, 9, 9 ],
        ],
        forces => [
            [ 0, 2, 3 ],
            [ 0, 3, 4 ],
            [ 0, 2, 3 ],
            [ 0, 3, 4 ],
            [ 0, 2, 3 ],
            [ 0, 3, 4 ],
            [ 0, 2, 3 ],
            [ 0, 3, 4 ],
            [ 0, 2, 3 ],
            [ 0, 3, 4 ],
        ]
    );
}
'Test creation of an obj3';

is_deeply( $obj3->delta_coords( 0, 1 ), [ 1, 1, 1 ], 'delta_coords' );
is_deeply( $obj3->delta_charges( 0, 9 ), 9, 'delta_charges' );

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

#Math::VectorReal example from the synopsis

$obj3->set_coords( $_, vector( @{ $obj3->get_coords($_) } ) )
  foreach ( 0 .. $obj3->count_coords - 1 );

my $sum_vec = vector( 0, 0, 0 );
$sum_vec += $_ foreach $obj3->all_coords;
my $avg = $sum_vec / $obj3->count_coords;

is_deeply(
    [ $avg->array ],
    [ 4.5, 4.5, 4.5 ],
    "average coordinates via Math::VectorReal"
);

done_testing();

