#!/usr/bin/env perl

use strict;
use warnings;
use Test::More;
use HackaMol;
use Path::Tiny;
use Time::HiRes qw(time);

my $ethane_json = path("t/lib/ethane.json");

my $mol = HackaMol::Molecule->new()->from_json($ethane_json->slurp);

$mol->print_xyz;
my $t0 = time;
$mol = HackaMol->new->pdbid_mol("2cba");
my $t1 = time;
my $json = $mol->to_json(1);
my $t2 = time;
my $mol2 = HackaMol::Molecule->new()->from_json($json);
my $t3 = time;

printf("pdbparse %.6f\n", $t1-$t0);
printf("json write%.6f\n", $t2-$t1);
printf("json read %.6f\n", $t3-$t2);

$mol2->print_pdb;

done_testing();
