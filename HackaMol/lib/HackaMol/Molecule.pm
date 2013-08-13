package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/roles';
use MooseX::Storage;
with Storage('io' => 'StorableFile');#, 'PhysVecRole';


1;

__END__

=pod

=head1 NAME

Molecule

=cut
