# Demian Riccardi August, 22, 2013
#
# Description
# grep out the backbone atoms and rotate the dihedrals to the angle read in
# adding the sidechains shouldn't be too difficult.  Just have to identify which
# atoms are moving
use Modern::Perl;
use lib 'lib','t/lib';
use HackaMol::Molecule;
use HackaMol::Dihedral;
use PDBintoAtoms qw(readinto_atoms);
use Time::HiRes qw(time);
use Scalar::Util qw(refaddr);

my $t1 = time; 
my $angle = shift ;
$angle = 180 unless (defined($angle));

my @all_atoms = readinto_atoms("t/lib/1L2Y.pdb");
#my @all_atoms = readinto_atoms("t/lib/2LL5.pdb");
#to keep example simple, keep only the backbone
my @atoms = grep {
               $_->name eq 'N'  or
               $_->name eq 'CA' or
               $_->name eq 'C'
              } @all_atoms;
#reset iatom
$atoms[$_]->iatom($_) foreach 0 .. $#atoms;

my $max_t = $atoms[0]->count_coords -1;
my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

my @dihedrals ; 

# build the dihedrals 
my $k = 0;
while ($k+3 <= $#atoms){
  my $name; 
  $name .= $_->name.$_->resid foreach (@atoms[$k .. $k+3]);
  push @dihedrals, Dihedral->new(name=>$name, atoms=>[ @atoms[$k .. $k+3] ]);
  $k++;
}

my $natoms = $mol->count_atoms;
my $t = 0;
$mol->t($t);

my @iatoms = map{$_->iatom} @atoms;

foreach my $dihe (@dihedrals){

  my $ratom1 = $dihe->get_atoms(1);
  my $ratom2 = $dihe->get_atoms(2);

  # atoms from nterm to ratom1 and from ratom2 to cterm
  my @nterm = 0 .. $ratom1->iatom - 1;
  my @cterm = $ratom2->iatom +1 .. $natoms-1;

  # use the smaller list for rotation
  my $r_these = \@nterm;
  $r_these = \@cterm if (@nterm > @cterm);
 
  #set angle to rotate
  my $rang = -1*($dihe->dihe + $angle) ;
  #switch nterm to cterm switches sign on angle
  $rang *= -1 if (@nterm>@cterm); 
  my @slice = @atoms[@{ $r_these}]; 
  #ready to rotate!
  $mol->dihedral_rotate_atoms($dihe,$rang,\@slice);

}
 
print "$natoms \n\n"; 
printf("%5s %8.3f %8.3f %8.3f\n", $_->symbol, @{$_->get_coords($t)}) foreach @atoms;

my $t2 = time;

printf("time: %10.6f\n", $t2-$t1);
