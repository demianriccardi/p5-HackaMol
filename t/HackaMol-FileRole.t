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
input_fn
output_fn
log_fn
);
my @methods = qw(
);

map has_attribute_ok( 'HackaMol', $_ ), @attributes;
map can_ok( 'HackaMol', $_ ), @methods;
my $obj;

lives_ok {
    $obj = HackaMol->new(
                         scratch   => "t/tmp",
                         input_fn  => "t/tmp/blah.inp",
                         output_fn => "t/tmp/blah.out",
                         log_fn    => "t/tmp/blah.log",
                        );
}
'Test creation of an obj with files';

is($obj->input_fn,  't/tmp/blah.inp', "blah.inp name");
is($obj->output_fn, 't/tmp/blah.out', "blah.out name");
is($obj->log_fn,    't/tmp/blah.log', "blah.log name");

my $fhlog = $obj->log_fn->openw;

$obj->scratch->mkpath;
dir_exists_ok($obj->scratch, 'scratch directory does exist after mkpath');
print $fhlog "test 1\n";

$obj->input_fn->spew(join("\n", map{sprintf("testing %i",$_)} 1 .. 10 ));
$obj->input_fn->copy_to($obj->output_fn);
my @outlines = $obj->output_fn->slurp;
my @inlines  = $obj->input_fn->slurp;
is_deeply(\@inlines,\@outlines,"input written, copied to output, read back in");
print $fhlog "test 2\n";

my $string = join("\n", map{sprintf("testing %i",$_)} 11 .. 20 );
$obj->input_fn->spew($string);
$obj->output_fn->spew($string);
my $lines = $obj->input_fn->slurp;
is($lines,$string,"input written anew and slurped up again");

print $fhlog "test 3";
close($fhlog);

my $loglines = $obj->log_fn->slurp;
my $logstring = join("\n", map{sprintf("test %i",$_)} 1 .. 3);

is($loglines,$logstring,"log written 3 times and slurped");

{ #test two open fh
  my $fhi = $obj->input_fn->openr;
  my $fho = $obj->output_fn->openr;
  my @li = <$fhi>;
  my @lo = <$fho>;
  is_deeply(\@li,\@lo,"input/output filehandles opened and read");
}

done_testing();
