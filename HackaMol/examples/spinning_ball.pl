use Modern::Perl;
use lib 'lib','t/lib';
use HackaMol::Molecule;
use PDBintoAtoms qw(readinto_atoms);
use Math::Vector::Real;


my @atoms = readinto_atoms("t/lib/1L2Y.pdb");
my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

$mol->translate(-$mol->COM);
$mol->translate(V(90,0,0));

foreach (1 .. 360){
  $mol->translate(V(-0.5,0,0));
  $mol->rotate(V(1,1,1), 10, $mol->COM);
  print_xyz($mol);
}  

sub print_xyz {
  my $mol = shift;
  print $mol->count_atoms;
  print "\n\n";
  printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;
}
