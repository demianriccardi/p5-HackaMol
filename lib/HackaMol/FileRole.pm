package HackaMol::FileRole;

#ABSTRACT:  
use 5.008;
use Moose::Role;
use MooseX::Types::Path::Class;

has 'input_fn' => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    coerce   => 1,
);
 
has 'output_fn' => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    coerce   => 1,
);

has 'log_fn' => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    coerce   => 1,
); 

no Moose::Role;

1;

__END__
=head1 DESCRIPTION

This role adds files (log_fn,input_fn,output_fn) to a class. This is still a work in progress, 
and it will probably change (suggestions welcome).  The goal is to reduce the amount
code required for creating inputs, processing outputs, and monitoring it all
in a platform independent way. MooseX::Types::Path::Class is used to
coerce the attributes into Path::Class::File objects. See Path::Class for 
associated methods. 

=attr log_fn 

isa Path::Class::File that is 'ro'  

Intended for logging, but there's nothing enforcing that for now.

=attr input_fn

isa Path::Class::File that is 'ro'

writing input, but there's nothing enforcing that for now.

=attr output_fn

isa Path::Class::File that is 'ro'

reading output, but there's nothing enforcing that for now.

