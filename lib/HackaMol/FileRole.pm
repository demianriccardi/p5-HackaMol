package HackaMol::FileRole;

#ABSTRACT:  
use 5.008;
use Moose::Role;
use MooseX::Types::Path::Tiny qw/Path Paths AbsPath/;

has 'in_fn' => (
    is       => 'ro',
    isa      => Path,
    coerce   => 1,
);
 
has 'out_fn' => (
    is       => 'ro',
    isa      => Path,
    coerce   => 1,
);

has 'log_fn' => (
    is       => 'ro',
    isa      => Path,
    coerce   => 1,
); 

has forts => (
    is       => 'ro',
    isa      => Paths,
    coerce   => 1,
);

no Moose::Role;

1;

__END__
=head1 DESCRIPTION

This role adds files (log_fn,in_fn,out_fn) to a class. This is still a work in progress, 
and it will probably change (suggestions welcome).  The goal is to reduce the amount
code required for creating inputs, processing outputs, and monitoring it all
in a platform independent way. MooseX::Types::Path::Class is used to
coerce the attributes into Path::Class::File objects. See Path::Class for 
associated methods. 

=attr log_fn 

isa Path::Class::File that is 'ro'  

Intended for logging, but there's nothing enforcing that for now.

=attr in_fn

isa Path::Class::File that is 'ro'

writing input, but there's nothing enforcing that for now.

=attr out_fn

isa Path::Class::File that is 'ro'

reading output, but there's nothing enforcing that for now.

=attr fort1_fn fort2_fn fort3_fn fort4_fn fort5_fn

isa Path::Class::File that is 'rw'

a place for those extra, annoying files


