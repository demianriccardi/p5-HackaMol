package PhysicalVec;
use Moose::Role;
use Carp;

my @t_dep = qw(coords forces charges); 

has 'name'   ,    is => 'ro', isa => 'Str' ;
has 'mass'   ,    is => 'rw', isa => 'Num' ;  
# t ScalarRef allows one to set everything to the same t... but that would be 
# more useful if there was a fast function that would return the value
# regardless of whether it was a Scalar or ScalarRef
has 't'      ,    is => 'rw', isa => 'Int|ScalarRef' ,  default => 0;
has 'origin' => (
                 is      => 'rw', 
                 isa     => 'ArrayRef', 
                 default => sub{[0,0,0]},
                 lazy    => 1,
                );

#charge was here once. we have tcharges now, charge into Molecule
has $_ =>   (
              is  => 'rw', 
              isa =>'Num', 
              predicate => 'has_charge', 
              lazy=> 1, 
            ) for qw(mass); 

has "_t$_"  => (
                traits   => [ 'Array' ],
                isa      => 'ArrayRef',
                default  => sub { [] },
                handles  =>
                {
                  "add_$_" => 'push'    ,
                  "get_$_" => 'get'     ,
                  "set_$_" => 'set'     ,
                  "all_$_" => 'elements',
                "count_$_" => 'count'   ,
                },
                lazy   => 1,
               ) for @t_dep;

has 'units'  ,    is => 'rw', isa => 'Str'  ; #flag for future use [SI]
has 'xyzfree' => (  
                  is      => 'rw', 
                  isa     => 'ArrayRef[Int]', 
                  default => sub{[1,1,1]},
                  lazy    => 1,
                 );

do need buildargs


sub charge {
  my $self = shift;
  if (@_){
    $self->set_charges($self->t,$_[0]);
  }
  else {
    return ($self->get_charges($self->t);
  }
}
1;


