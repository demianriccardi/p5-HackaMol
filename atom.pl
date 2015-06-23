use Modern::Perl;
use Math::Vector::Real;
use HackaMol;

my $atom = HackaMol::Atom->new(Z=>1, coords=>[V(1,1,2)]);
$atom->store('atom.json');
print $atom->dump;

my $atom_again = HackaMol::Atom->load('atom.json');
print $atom_again->dump;

