package PdbRole;
use Moose::Role;

has 'record_name'  , is => 'ro', isa => 'Str';
has 'serial'       , is => 'rw', isa => 'Int';
has 'occ'          , is => 'ro', isa => 'Num';
has 'bfact'        , is => 'ro', isa => 'Num';
has 'resname'      , is => 'ro', isa => 'Str';
has 'chain'        , is => 'ro', isa => 'Str';
has 'altloc'       , is => 'ro', isa => 'Str';
has 'resid'        , is => 'rw', isa => 'Int';
has 'segid'        , is => 'rw', isa => 'Int';
has 'iatom'        , is => 'rw', isa => 'Int';

no Moose::Role;
1;
