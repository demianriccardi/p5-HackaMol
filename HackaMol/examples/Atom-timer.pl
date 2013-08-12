use Modern::Perl;
use lib 'lib/HackaMol';
use Time::HiRes qw(time);
use Benchmark qw(cmpthese);
use Atom;

my $natoms = 100000; 
print "Atom-timer will time the construction of and array of  $natoms atoms to give idea about speed\n";

my $t1 = time;
my @atoms = map{
                 Atom->new(name => 'Mercury', coords => [[0,0,0]], symbol => "HG")
               } 1 .. $natoms;

my $t2 = time;
printf ("time to load $natoms atoms: %10.3f\n", $t2-$t1);

$_->push_coords([2.000 ,2.000,2.000]) foreach @atoms;
my $t3 = time;
printf ("time to push_coords( [2.000,2.000,2.000] ) for $natoms atoms: %10.3f\n", $t3-$t2);

$_->set_coords(0,[1.000 ,1.000,1.000]) foreach @atoms;
my $t4 = time;
printf ("time to set_coords(0, [1.000 ,1.000,1.000]) for $natoms atoms: %10.3f\n", $t4-$t3);

print "dump the last atom: ", $#atoms, "\n";
print $atoms[$#atoms]->dump;

print "construction comparison benchmarks Z vs Symbol\n";
cmpthese (50000, {
                  'symb construct' => sub{my $atom = Atom->new(name=>'shit', coords => [[0,0,0]], symbol => "HG" ) },
                  'Z    construct' => sub{my $atom = Atom->new(name=>'shit', coords => [[0,0,0]], Z      =>  80  ) },
                 }
          );




