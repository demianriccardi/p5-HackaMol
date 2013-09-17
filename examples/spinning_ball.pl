use Modern::Perl;
use HackaMol;
use Math::Vector::Real;

my $hack  = HackaMol->new( name=>"hackitup" );
my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

my $mol = HackaMol::Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

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
  foreach my $atom ($mol->all_atoms) {
    printf("%5s %12.6f %12.6f %12.6f\n", $atom->Z, @{$atom->xyz});
  }
}
