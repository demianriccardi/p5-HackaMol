use Modern::Perl;
use WWW::Search;
use Data::Dumper;
use lib 'lib';
use HackaMol::Atom;
use Math::Vector::Real;

my $search = new WWW::Search('PubChem');

my @ids = qw/6037/;
$search->native_query( \@ids );
my $chem = $search->next_result;
my $smiles = $chem->{smiles};
my @a = `obabel -:"$smiles" -oxyz --gen3D`;
chomp @a;

my @atoms = map  {
                  HackaMol::Atom->new(symbol=>$_->[0],coords=>[V($_->[1],$_->[2],$_->[3])])
                 } 
            map  { [split] } 
            grep { m/\w+\s+-*\d+.\d+/ } @a;

print Dumper \@atoms;

foreach (my $i = 0 ; $i < @atoms; $i++) {
  foreach (my $j = $i+1 ; $j < @atoms; $j++){
    
  }
}


