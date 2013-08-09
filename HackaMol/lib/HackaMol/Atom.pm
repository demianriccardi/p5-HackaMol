package Atom;
use Moose;
use lib 'roles'; 
use MooseX::Storage;
with Storage('io' => 'StorableFile');

with 'PhysVecRole';

1;

__END__

=pod

=head1 NAME

Atom

=cut
