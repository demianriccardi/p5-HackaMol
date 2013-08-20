package HackaMol;
# ABSTRACT: object-oriented library for multi-scale molecular hacking 
use lib 'lib';
use Math::Vector::Real qw(V);
require Exporter;
our @ISA=qw(Exporter);
our @EXPORT = qw(V);
use HackaMol::Atom;
use HackaMol::Molecule;
use HackaMol::Bond;
use HackaMol::Angle;
use HackaMol::Dihedral;


1;

__END__

