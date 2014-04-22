package HackaMol::ExeRole;

#ABSTRACT:  
use 5.008;
use Moose::Role;
use Carp;

has 'exe'       => (
                      is  => 'rw',
                      isa => 'Str',
                      predicate => 'has_exe',
                      clearer   => 'clear_exe',
                   );


#options after input filename... bleh
# ie. dftd3 shit.xyz -func b3pw91 -bj              
has 'exe_endops'       => (
                      is  => 'rw',
                      isa => 'Str',
                      predicate => 'has_exe_endops',
                      clearer   => 'clear_exe_endops',
                   );

# this will be called by capture::tiny
has 'command'       => (
                        is  => 'rw',
                        isa => 'Str',
                        predicate => 'has_command',
                        clearer   => 'clear_command',
                       );


sub exists_exe {
  my $self = shift;
  if (-e $self->exe){
    return 1;
  }
  else {
    carp $self->exe . " does not exist";
    return 0;
  } 
}

no Moose::Role;

1;

__END__
=head1 DESCRIPTION

This role adds executables/commands for running external programs. This is still a work in progress, and it will 
probably change (suggestions and help very much welcome!).  The goal is to reduce the amount code required for 
building interfaces to external programs to be run on inputs to generate output in some directory that may be 
temporary... or not.  Of course, exes do all sorts of things where other files may be written. Requirements (e.g. a 
method that tests functionality) for interfaces are still under development. Considering the trickiness of this 
sort of abstraction, it will cowardly left to the extensions to figure out. Recommendation: Capture::Tiny 

=attr command 

isa Str that is rw   

to be constructed from exe, exe_endops, in_fn, out_fn, etc. Then run and 
captured, which is left to scripts/interfaces

=attr exe 

isa Str that is rw

the program to be run.  $self->command($self->exe . " < " . $self->in_fn . " > " . $self->out_fn); 

=attr exe_endops 

isa Str that is rw

options to be catenated to the end of the exe.  For those command line tools that use options after input filename

=method exists_exe

return 1 if exe exists, carp warning and return 0 if exe does not exist


