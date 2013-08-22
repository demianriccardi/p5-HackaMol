use Modern::Perl;
use lib 'lib/HackaMol','t/lib';
use Molecule;
use Dihedral;
use PDBintoAtoms qw(readinto_atoms);
use Time::HiRes qw(time);

my $t1 = time; 
my @atoms = readinto_atoms("t/lib/1L2Y.pdb");
my $max_t = $atoms[0]->count_coords -1;
my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);
my @N_CA_C = grep {
                   $_->name eq 'N'  or 
                   $_->name eq 'CA' or 
                   $_->name eq 'C'   
                  } @atoms;
my @dihedrals ;

# abcdefgh
my $k = 0;
while ($k+3 <= $#N_CA_C){
  my $name; 
  $name .= $_->name.$_->resid foreach (@N_CA_C[$k .. $k+3]);
  push @dihedrals, Dihedral->new(name=>$name, atoms=>[ @N_CA_C[$k .. $k+3] ]);
  $k++;
}

printf("%20s %10.2f\n", $_->name , $_->dihe) foreach @dihedrals;
my $t2 = time;

printf("time: %10.6f\n", $t2-$t1);

