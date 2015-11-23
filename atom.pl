use Modern::Perl;
use Math::Vector::Real;
use HackaMol;

#my $mol = HackaMol->new->pdbid_mol("2cba");

my $mol = HackaMol::Atom->new(name=>'shit', Z=>1, coords =>[
V(1,2,3), 
V(4,5,6),
V(4,5,6),
V(4,5,6),
V(4,5,6),
V(4,5,6),
V(4,5,6),
]
);#new->pdbid_mol("2cba");
$mol->store('shit.json');

#exit;
my $mol2 = HackaMol::Atom->load('shit.json');

print $mol2->dump;

print abs($mol2->origin);

