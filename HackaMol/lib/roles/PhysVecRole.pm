package PhysVecRole;
# ABSTRACT: PhysVecRole is at the heart of HackaMol Atoms and Molecules.  As documented below,
# blah
use Moose::Role;
use Carp;

my @t_dep = qw(coords forces charges); 

has 'name'   ,    is => 'ro', isa => 'Str' ;
has 'mass'   ,    is => 'rw', isa => 'Num' ;  

=attribute
t

isa Int or ScalarRef that is rw with default of 0

t is the index of the PhysVecRol attributes that have ARRAY traits. It is important to accessing
_tcoords, _tforces, and _tcharges.   Setting t to a 
ScalarRef would allow one to set all objects to the same t.

my $t  = 0;
my $rt = \$t;  
$_->t($rt) for (@objects);
$t = 1; # change t for all objects.

the physical vectorScalarRef allows one to set everything to the same t... but that would be 
 more useful if there was a fast function that would return the value
 regardless of whether it was a Scalar or ScalarRef
=cut
has 't'      ,    is => 'rw', isa => 'Int|ScalarRef' ,  default => 0;
has 'origin' => (
                 is      => 'rw', 
                 isa     => 'ArrayRef', 
                 default => sub{[0,0,0]},
                 lazy    => 1,
                );

has $_ =>   (
              is  => 'rw', 
              isa =>'Num', 
              predicate => 'has_charge', 
              lazy=> 1, 
            ) for qw(mass); 

=private_attribute
_tcharges

isa ArrayRef that is lazy with public ARRAY traits described in ARRAY_METHODS
 
gives atoms and molecules t-dependent arrays of charges for the purpose of 
analysis.  e.g. store and analyze atomic charges from a quantum 
mechanical molecule in several intramolecular configurations or a fixed
configuration in varied potentials.   
=cut

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
=attribute
xyzfree

lazy, default value [1,1,1]

=cut
has 'xyzfree' => (  
                  is      => 'rw', 
                  isa     => 'ArrayRef[Int]', 
                  default => sub{[1,1,1]},
                  lazy    => 1,
                 );

=attribute
charge

lazy, default: (build from t-dependent charges) or 0
 
As decribed in ARRAY_ATTRIBUTES, atoms and molecules have t-dependent arrays 
of charges for the purpose of analysis.  e.g. one could store and analyze 
atomic charges 
from a quantum mechanical molecule in several intramolecular configurations 
or under varying environmental influences.  

Often thinking in terms of a given atom/molecule having a "charge" is more
intuitive and convenient.  The default value of the "charge" attribute is 
the t index of "charges" or 0.0 if "charges" have been cleared. 

=cut

has 'charge' => (
     is         => 'rw',
     isa        => 'Num',
     builder    => '_build_charge',
     lazy       => 1,
                ); 

sub _build_charge {
  my $self = shift;
  if ($self->count_charges){
    return ( $self->get_charges($self->t) );
  }
  else {
    return (0);
  }
}


1;


