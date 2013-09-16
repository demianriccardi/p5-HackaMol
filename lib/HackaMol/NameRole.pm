package HackaMol::NameRole;
#ABSTRACT: provides name attribute 
use 5.008;
use Moose::Role;

has 'name', is => 'rw', isa => 'Str';

no Moose::Role;

1;

__END__
=head1 DESCRIPTION

simple role for the shared attribute 'name'. isa Str that is rw. useful for labeling, 
bookkeeping...
