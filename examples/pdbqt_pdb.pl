#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;

my $hack = HackaMol->new( name => "hackitup" );
my @atoms = $hack->read_file_atoms("t/lib/test.pdbqt");

# all coordinates from NMR ensemble are loaded into atoms
my $mol = HackaMol::Molecule->new(
    name  => 'somedrug',
    atoms => [@atoms]
);

foreach my $t (0 .. $mol->tmax){
  $mol->gt($t);
  $mol->print_pdb;
}

