#!/usr/bin/env perl

use strict;
use warnings;
use Test::Moose;
use Test::More;
use Test::Warn;
use Test::Dir;
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
    $obj = HackaMol->new();
}
'Test creation of an obj with nothing';

lives_ok {
    $obj = HackaMol->new(homedir=>"t", scratch => "t/tmp", data=> "t/lib");
}
'Test creation of an obj with directories';

is($obj->scratch, 't/tmp', "scratch set ok");
is($obj->data   , 't/lib', "data    set ok");
is($obj->homedir, 't',     "homedir set ok");
dir_exists_ok($obj->data,        'data directory exists');
dir_exists_ok($obj->homedir,     'homedir directory exists');
dir_not_exists_ok($obj->scratch, 'scratch directory does not exist');
$obj->scratch->mkpath;
dir_exists_ok($obj->scratch, 'scratch directory does exist after mkpath');
$obj->scratch->remove;
dir_not_exists_ok($obj->scratch, 'scratch directory deleted!');
$obj->scratch->mkpath;
dir_exists_ok($obj->scratch, 'scratch directory recreated');
my @storeit;
foreach (0 .. 2){
  my $fh = $obj->scratch->tempfile;
  print $fh $_ ;
  seek($fh, 0,0);
  while (my $line = <$fh>){
    is ($line, $_, "tempfile print/read")
  }
}

is (scalar($obj->scratch->children), 0,     "tempfiles were destroyed");

foreach (0 .. 2){
  my $file = $obj->scratch."/$_.txt";
  system ("touch $file");
}

is (scalar($obj->scratch->children), 3,     "touched 3 files");
my @txts = sort $obj->scratch->children;
my @list = sort qw(t/tmp/0.txt t/tmp/1.txt t/tmp/2.txt);
is_deeply(\@txts,\@list, "return contents of scratch");

$obj->scratch->rmtree;
$obj->scratch->remove;
dir_not_exists_ok($obj->scratch, 'scratch directory deleted!');

done_testing();
