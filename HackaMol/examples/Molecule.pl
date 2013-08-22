use Modern::Perl;
use Math::Vector::Real;
use lib 'lib/HackaMol';
use Time::HiRes qw(time);
use Benchmark qw(cmpthese);
use Scalar::Util qw(refaddr);
use Atom;
use Molecule;

my $natoms = 100000;
print
"Atom-timer will time the construction of and array of  $natoms atoms to give idea about speed\n";

my $t1 = time;
my @atoms = map { Atom->new( Z => 80 ), coords =>[V(1,2,3)] } 1 .. $natoms;
my $molecule = Molecule->new(name => "molecule");

my $t2 = time;
printf( "time to Atom->new(Z => 80) for $natoms atoms: %10.3f\n", $t2 - $t1 );

exit;
