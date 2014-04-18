#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Warn;
use Test::Fatal qw(lives_ok);
use HackaMol;                # v0.001;#To test for version availability

my @attributes = qw(
data
scratch
homedir
);
my @methods = qw(
);

map has_attribute_ok( 'HackaMol', $_ ), @attributes;
map can_ok( 'HackaMol', $_ ), @methods;
my $obj;
lives_ok {
    $obj = HackaMol->new(name=>'scratch');
}
'Test creation of an obj';

done_testing();
