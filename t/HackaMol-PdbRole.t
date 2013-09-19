#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Warn;
use Test::Fatal qw(lives_ok);
use MooseX::ClassCompositor;    #use this for testing roles
use HackaMol::PdbRole;                # v0.001;#To test for version availability

my @attributes = qw(
record_name
serial    
occ       
bfact     
resname   
chain     
altloc    
resid     
iatom     
icode     
pdbid     
segid    
);
my @methods = qw(
);

my $class = MooseX::ClassCompositor->new( { 
                                            class_basename => 'Test', 
                                          } )->class_for('HackaMol::PdbRole');

map has_attribute_ok( $class, $_ ), @attributes;
map can_ok( $class, $_ ), @methods;
my $obj;
lives_ok {
    $obj = $class->new();
}
'Test creation of an obj';

is($obj->record_name, 'HETATM', 'record_name default');
is($obj->serial     , 0       , 'serial      default');
is($obj->occ        , 1.0     , 'occ         default');
is($obj->bfact      , 20.0    , 'bfact       default');
is($obj->resname    , 'ALA'   , 'resname     default');
is($obj->chain      , 'AA'    , 'chain       default');
is($obj->altloc     , ' '     , 'altloc      default');
is($obj->resid      , 64      , 'resid       default');
is($obj->iatom      , 0       , 'iatom       default');
is($obj->icode      , 'X'     , 'icode       default');
is($obj->pdbid      , '2CBA'  , 'pdbid       default');
is($obj->segid      , 'TIP3'  , 'segid       default');

done_testing();
