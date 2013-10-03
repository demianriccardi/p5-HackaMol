#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;
use Math::Vector::Real;
use Time::HiRes qw(time);

my $t1 = time;
my $hack = HackaMol->new( name => "hackitup" );

my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");
my $max_t = $atoms[0]->count_coords - 1;
my $mol   = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms] );

print_xyz($mol);

$mol->translate( -$mol->COM );
$mol->translate( V( 10, 10, 10 ) );

print_xyz($mol);

$mol->rotate( V( 1, 0, 0 ), 180, V( 10, 10, 10 ) );

$mol->print_xyz;

my @groups = $hack->group_by_atom_attr( 'resid', $mol->all_atoms );
$mol->push_groups(@groups);

$_->rotate( V( 1, 1, 1 ), 60, $_->COM, 1 ) foreach $mol->all_groups;

$mol->gt(1);
$mol->print_xyz;

