#!/usr/bin/env perl
use Modern::Perl;
use HackaMol;
use Time::HiRes qw(time);

my $t1 = time;

my $hack = HackaMol->new( name => "hackitup" );

#my @atoms = readinto_atoms("t/lib/2CBA.pdb");
my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");
my $max_t = $atoms[0]->count_coords - 1;
my $mol   = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms] );

#backbone
my @N_CA_C =
  grep { $_->name eq 'N' or $_->name eq 'CA' or $_->name eq 'C' } @atoms;

my @dihedrals = $hack->build_dihedrals(@N_CA_C);

foreach my $dihe (@dihedrals) {
    printf( "%20s ", $dihe->name );
    do { $dihe->gt($_); printf( "%7.2f ", $dihe->dihe_deg ) }
      foreach 0 .. $max_t / 4;
    print "\n";

    # %10.2f (%.2f)\n", $dihe->name , avg_rmsd(@vals));
}

my $t2 = time;

printf( "time: %10.6f\n", $t2 - $t1 );

sub avg {
    my $sum = 0;
    $sum += $_ foreach @_;
    return ( $sum / scalar(@_) );
}

sub avg_rmsd {
    my $avg = avg(@_);
    my $sum = 0;
    $sum += ( $_ - $avg )**2 foreach @_;
    return ( $avg, sqrt( $sum / scalar(@_) ) );
}
