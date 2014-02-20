#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Warn;
use Test::Fatal qw(lives_ok);
use HackaMol::Atom;                # v0.001;#To test for version availability

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

map has_attribute_ok( 'HackaMol::Atom', $_ ), @attributes;
map can_ok( 'HackaMol::Atom', $_ ), @methods;
my $obj;
lives_ok {
    $obj = HackaMol::Atom->new(Z=>1);
}
'Test creation of an obj';

is($obj->record_name, 'HETATM', 'record_name default');
is($obj->occ        , 1.0     , 'occ         default');
is($obj->bfact      , 20.0    , 'bfact       default');
is($obj->resname    , ' '     , 'resname     default');
is($obj->chain      , ' '     , 'chain       default');
is($obj->altloc     , ' '     , 'altloc      default');
is($obj->resid      , 0       , 'resid       default');
is($obj->iatom      , 0       , 'iatom       default');
is($obj->serial     , 0       , 'serial      default');
is($obj->icode      , ' '     , 'icode       default');
is($obj->pdbid      , ' '     , 'pdbid       default');
is($obj->segid      , ' '     , 'segid       default');

done_testing();
