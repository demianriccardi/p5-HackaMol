use Modern::Perl;
use Math::Vector::Real;
use Math::Vector::Real::Neighbors;
use Math::Vector::Real::Random;
use Data::Dumper;
 
my @v = map Math::Vector::Real->random_normal(3), 0..1000;
 
my @ixs = Math::Vector::Real::Neighbors->neighbors(@v);

foreach my $ivec (0 .. $#v){
  printf ("C%i %10.3f %10.3f %10.3f\n",$ivec, @{$v[$ivec]});
  printf ("N%i %10.3f %10.3f %10.3f\n",$ivec, @{$v[$ixs[$ivec]]});
}


