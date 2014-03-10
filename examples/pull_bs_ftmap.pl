#!/usr/bin/env perl
# Demian Riccardi 03/09/2014
# print out COM of the binding sites predicted by FTMap
#
# The clusters of small molecule probes begin with HEADER: 
#   HEADER crosscluster.\d+.\d+.pdb 
#   
# This dum script will: 
#
#   1. write them all out into temporary xyz files 
#   2. read them in and print out the COM of the cluster
#   3. unlink the temporary xyz file
#
# need to use xyz file, because FTMap does not write a proper PDB 
# may need a dirty pdb reader
#
use Modern::Perl;
use HackaMol;
use File::Temp qw/ tempfile tempdir mkstemps/;
use File::chdir;
use FileHandle;
use Data::Dumper;

my $file = shift;
my $fh   = FileHandle->new("< $file");

my @CLUSTERS;

my $i = 0;
while (<$fh>){
  if (/HEADER crosscluster/../REMARK/ ){
    if (/ATOM\s+\d+/){
      my ($symbol, $xyz) = unpack "x13A1x16A24", $_;
      push @{$CLUSTERS[$i]}, "$symbol $xyz\n";
    }
    $i++ if (/REMARK/); 
  }
}

my $dir = tempdir( CLEANUP=>1 );

my $hack = HackaMol->new();
foreach my $clust (@CLUSTERS){
  my ($fh, $filename) = mkstemps( "hackXXXX",'.xyz');
  #write xyz file
  print $fh scalar(@{$clust})."\n\n";
  print $fh $_ foreach @{$clust};
  close($fh);
  #read in xyz and print COM
  my $mol = $hack->read_file_mol($filename);
  printf("%10.3f %10.3f %10.3f\n", @{$mol->COM});
  unlink($filename);
}


