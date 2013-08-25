# Demian Riccardi August, 22, 2013
#
# Description
# grep out the backbone atoms and rotate the dihedrals to the angle read in
# adding the sidechains shouldn't be too difficult.  Just have to identify which
# atoms are moving
use Modern::Perl;
use lib 'lib/HackaMol','t/lib';
use Molecule;
use Dihedral;
use PDBintoAtoms qw(readinto_atoms);
use Math::Vector::Real;
use Time::HiRes qw(time);
use Scalar::Util qw(refaddr);

my $t1 = time; 
my $angle = shift ;
$angle = 180 unless (defined($angle));

my @atoms = readinto_atoms("t/lib/1L2Y.pdb");
my $max_t = $atoms[0]->count_coords -1;
my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

$mol->push_groups_by_atom_attr('resid');
$_->translate(-$_->COM,1) foreach $mol->all_groups; 

$mol->print_xyz(0);
$mol->print_xyz(1);

my %Z;
foreach my $atom ($mol->all_atoms){
  my $z = $atom->Z;
  $Z{$z}++;
  $atom->push_coords(V(3*$z,$Z{$z}/4,0));
}

$mol->print_xyz($max_t+1);
$_->push_coords(V(3*$_->Z,0,0)) foreach $mol->all_atoms;
$mol->print_xyz($max_t+2);
$_->push_coords(V(0,0,0)) foreach $mol->all_atoms;
$mol->print_xyz($max_t+3);


