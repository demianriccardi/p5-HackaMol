use Modern::Perl;
use lib 'lib/HackaMol';
use Test::More;
use Time::HiRes qw(time);
use Atom;

my $atom1 = Atom->new(name=>'shit', coord =>[1,0,0], Z=>6);

$atom1->push_coords([1,1,1]);

print $atom1->dump;
say $atom1->symbol;
print $atom1->dump;

my $atom2 = Atom->new(name=>'shit', coord =>[1,0,0], symbol=>"HG");

print $atom2->dump;
say $atom2->Z;
print $atom2->dump;
=shit
my $t1 = time;
my @atoms = map{
                 Atom->new(name=>'shit', coord => [0,0,0], symbol => "HG")
               } 1 .. 100000;
my $t2 = time;
printf ("time to load 100000 atoms: %10.3f\n", $t2-$t1);

$_->push_coords([2,2,2]) foreach @atoms;
my $t3 = time;
printf ("time to push coords on to 100000 atoms: %10.3f\n", $t3-$t2);

$_->set_coords(0,[1,1,1]) foreach @atoms;
my $t4 = time;
printf ("time to set coords on to 100000 atoms: %10.3f\n", $t4-$t3);
=cut


use Benchmark qw(cmpthese);

cmpthese (50000, {
                  'symb construct' => sub{my $atom = Atom->new(name=>'shit', coord => [0,0,0], symbol => "HG" ) },
                  'Z    construct' => sub{my $atom = Atom->new(name=>'shit', coord => [0,0,0], Z      =>  80  ) },
                 }
          );
