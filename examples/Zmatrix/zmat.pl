use HackaMol;
my $mol = HackaMol->new->read_file_mol(shift);
$mol->print_pdb;
