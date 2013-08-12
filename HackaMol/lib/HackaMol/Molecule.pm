package Molecule;
use Moose;
use lib 'lib/roles';
use MooseX::Storage;
with Storage('io' => 'StorableFile');

with 'PhysVecRole';

1;

__END__

=pod

=head1 NAME

Molecule

=cut
