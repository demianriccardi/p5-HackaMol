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
is($obj->resname    , 'UNK'   , 'resname     default');
is($obj->chain      , ' '     , 'chain       default');
is($obj->altloc     , ' '     , 'altloc      default');
is($obj->resid      , 1       , 'resid       default');
is($obj->iatom      , 0       , 'iatom       default');
is($obj->serial     , 1       , 'serial      default');
is($obj->icode      , ' '     , 'icode       default');
is($obj->pdbid      , ' '     , 'pdbid       default');
is($obj->segid      , ' '     , 'segid       default');

# aa321 tests
{
  my $atom = HackaMol::Atom->new(Z=>1);
  my %aa321 = (
    ALA=>'A', ARG=>'R', ASN=>'N', ASP=>'D', CYS=>'C',
    GLU=>'E', GLN=>'Q', GLY=>'G', HIS=>'H', ILE=>'I',
    LEU=>'L', LYS=>'K', MET=>'M', PHE=>'F', PRO=>'P',
    SER=>'S', THR=>'T', TRP=>'W', TYR=>'Y', VAL=>'V',
    SEC=>'U',
  );

  foreach my $res ((keys %aa321), (map{lc($_)} keys %aa321)){  
    $atom->resname($res);
    my $onelett = $aa321{uc($res)};
    is($atom->aa321, $onelett, "aa321 $res -> $onelett");
  }
  $atom->resname("TRD");
  my $x;
  warning_is { $x = $atom->aa321}
  "PDBRole> residue TRD name has no 1 letter code; return X",
    "warning for unrecognized residue name";
  is($x, 'X', "aa321 returns X for unrecognized residue name");

}

done_testing();
