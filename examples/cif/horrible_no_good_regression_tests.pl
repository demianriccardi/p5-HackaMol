use Modern::Perl;
use HackaMol;
use Data::Dumper;

#die "run this as much as you want! but it will take a long time to sync an empty dir";
#
my $bldr = HackaMol->new();
#my $file = $bldr->pdbid_local_path('1gna','cif');
foreach my $pdbid ( 
    #grep {m/3gtk/} 
    qw/
    1a2c 3gtk 2ahn 2ah8 3oe0 4gpg 2A0Z 207l 3gtk 2a2x 
    5a71 3j7p 1aq5 5aiy 3A34 2l9h
    /){
    my $file = $bldr->pdbid_local_path($pdbid,'cif');
    say $file->stringify;
    my $fh = $file->openr_raw;
    my $info = $bldr->read_cif_info($fh);
    print Dumper $info;
    say $info->{exp_method};
}
#print $mols->[0]->get_atoms(0)->dump;

