package HackaMol::ScratchRole;

#ABSTRACT:  
use 5.008;
use Moose::Role;
use MooseX::Types::Path::Tiny qw/Path AbsPath AbsPaths/;

has 'homedir' => (
    is       => 'ro',
    isa      => AbsPath,
    coerce   => 1,
);
 
has 'scratch' => (
    is       => 'ro',
    isa      => AbsPath,
    coerce   => 1,
);

has 'data' => (
    is       => 'ro',
    isa      => AbsPath,
    coerce   => 1,
); 

has 'dirs' => (
    is       => 'ro',
    isa      => AbsPaths,
    coerce   => 1,
); 

no Moose::Role;

1;

__END__
=head1 DESCRIPTION

This role adds directories to a class. This is still a work in progress, and it
will probably change (suggestions welcome).  The goal is to reduce the amount
code required for manipulating several paths used for work, and to allow 
scripts to be more platform independent. MooseX::Types::Path::Class is used to
coerce the attributes into Path::Class::Dir objects. See Path::Class for 
associated methods. 

=attr scratch 

isa Path::Class::Dir that is 'ro'  

Intended to be temporary, but there's nothing enforcing that for now.

=attr homedir

isa Path::Class::Dir that is 'ro'

Intended to be the mother ship, but there's nothing enforcing that for now.

=attr homedir

isa Path::Class::Dir that is 'ro'

Intended to be a place that data is retrieved from, but there's nothing 
enforcing that for now.

