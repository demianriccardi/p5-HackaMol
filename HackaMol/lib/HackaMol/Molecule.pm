package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/roles';
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecMVRRole';


1;

__END__

=pod

=head1 NAME

Molecule

=cut
