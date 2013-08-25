use Modern::Perl;
use lib 'lib/HackaMol','t/lib';
use Molecule;
use Dihedral;
use PDBintoAtoms qw(readinto_atoms);
use Math::Vector::Real;
use Time::HiRes qw(time);


my $t1 = time; 
my @atoms = readinto_atoms("t/lib/1L2Y.pdb");
my $max_t = $atoms[0]->count_coords -1;
my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

say $mol->count_atoms;
print "\n";
printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

$mol->translate(-$mol->COM);
$mol->translate(V(10,10,10));

say $mol->count_atoms;
print "\n";
printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

$mol->rotate(V(1,0,0), 180, V(10,10,10));

say $mol->count_atoms;
print "\n";
printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

$mol->push_groups_by_atom_attr('resid');
$_->rotate(V(1,1,1),60,$_->COM,1) foreach $mol->all_groups;

$mol->gt(1);

say $mol->count_atoms;
print "\n";
printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;


sub print_xyz {
  my $mol = shift;
  print $mol->count_atoms;
  print "\n\n";
  printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;
}
