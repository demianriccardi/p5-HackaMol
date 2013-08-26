use Test::Most;
use Test::Warnings;
use Test::Moose;
use MooseX::ClassCompositor;    #use this for testing roles
use roles::QmRole;                # v0.001;#To test for version availability

my @attributes = qw(
basis
ecp
multiplicity
basis_geom
dummy
);
my @methods = qw(
);

my $class = MooseX::ClassCompositor->new( { 
                                            class_basename => 'Test', 
                                          } )->class_for('QmRole');

map has_attribute_ok( $class, $_ ), @attributes;
map can_ok( $class, $_ ), @methods;
my $obj;
lives_ok {
    $obj = $class->new();
}
'Test creation of an obj';

is($obj->basis, '6-31+G*', 'basis default');
is($obj->multiplicity        , 1     , 'multiplicity default');

done_testing();
