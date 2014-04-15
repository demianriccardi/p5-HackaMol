#!/usr/bin/env perl
#
# DMR 4-15-2014
# simple script to help get started analyzing cofactors/metals in proteins 
#
# wget http://pdb.org/pdb/files/1C7D.pdb
#
# run:
# perl feheme.pl 1C7D.pdb 
#
use Modern::Perl;
use HackaMol;

my $hack  = HackaMol->new ;
# load all atoms into an array
my @atoms = $hack->read_file_atoms(shift);

#pull out iron atoms
my @Fes   = grep {$_->symbol eq "Fe"} @atoms;

#pull out hemes, no iron
my @hemes = grep {
                  $_->resname eq "HEM" and
                  $_->symbol  ne "Fe"
                 } @atoms;

# find the bonds to fe foreach wrt hemes, print out lengths
foreach my $fe (@Fes){
  my @fe_bonds  = $hack ->find_bonds_brute(
                          bond_atoms => [$fe],
                          candidates => [@hemes],
#
#                         candidates => [grep {$_->symbol eq "N"} @hemes]
#
                          fudge      => 0.45,
                          max_bonds  => 6,
  );
  foreach my $bond (@fe_bonds){
    my @sym = map{ $_->symbol } $bond->all_atoms;
    printf ("%4s %4s %10.3f\n", @sym, $bond->bond_length); 
  }
}


