package PhysVecRole;
# ABSTRACT: Provides the core of HackaMol Atom and Molecule classes. 
use Moose::Role;
use Carp;

has 'name'     , is => 'rw', isa => 'Str' ;

has 'mass'     , is => 'rw', isa => 'Num' , lazy => 1, default => 0;  

has 't'        , is => 'rw', isa => 'Int|ScalarRef' ,  default => 0;

has 'is_fixed' , is => 'rw', isa => 'Bool', lazy => 1, default   => 0;


my @t_dep = qw(coords forces charges); 

has "_t$_"  => (
                traits   => [ 'Array' ],
                isa      => 'ArrayRef',
                default  => sub { [] },
                handles  =>
                {
                  "add_$_"   => 'push'    ,
                  "get_$_"   => 'get'     ,
                  "set_$_"   => 'set'     ,
                  "all_$_"   => 'elements',
                  "clear_$_" => 'clear'   ,
                  "count_$_" => 'count'   ,
                },
                lazy   => 1,
               ) for @t_dep;


has 'units'  ,    is => 'rw', isa => 'Str'  ; #flag for future use [SI]

has 'origin' => (
                 is      => 'rw', 
                 isa     => 'ArrayRef', 
                 default => sub{[0,0,0]},
                 lazy    => 1,
                );

has 'xyzfree' => (  
                  is      => 'rw', 
                  isa     => 'ArrayRef[Int]', 
                  default => sub{[1,1,1]},
                  lazy    => 1,
                 );

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

__END__

=pod

=encoding utf8

=head1 SYNOPSIS

# instance of class that consumes the PhysVecRol

my $obj = Class_with_PhysVecRole->new( name => 'foo', t => 0 ); 

# add some charges

$obj->add_charges($_) foreach ( 0.3, 0.2, 0.1, -0.1, -0.2, -0.36 );

my $sum_charges = 0;

$sum_charges += $_ foreach $obj->all_charges;

print $obj->charge . "\n";

#set charge to the average charge

$obj->charge( $sum_charges / $obj->count_charges ); 

print $obj->charge . "\n";

#charge has no effect on array of charges
print $_ foreach $obj->all_charges;

# add some coordinates

$obj->add_coords($_) foreach ( [0,0,0], [1,1,1], [-1.0,2.0,-4.0] );

#example of flexibility: 

#Inplace conversion of coordinates to Math::VectorReal objects

use Math::VectorReal; # exports 'vector' method

$obj->set_coords( 
                  $_, vector( @{ $obj->get_coords($_) } ) 
                ) foreach (0 .. $obj->count_coords -1);

# now ready for Math::VectorReal fun

my $sum_vec = vector(0, 0, 0);

$sum_vec += $_ foreach $obj->all_coords;

print $sum_vec / $obj->count_coords; # average x,y,z coordinates

=head1 DESCRIPTION

The PhysVecRole provides the core attributes and methods shared between Atom 
and Molecule classes. Consuming this role gives Classes a place to store 
coordinates, forces, and charges, perhaps, over the course of a simulation 
or for a collection 
configurations for which all other object metadata (name, mass, etc) remains 
fixed. As such, the 't' attribute, described below, is important to understand. 
The PhysVecRole is very flexible. Several attribute types are left agnostic and all
are rw so that they may be filled with whatever the user defines them to be... on 
the fly.  This seems most intuitive from the perspective of carrying out computational 
work on molecules.  Thus, HackaMol bravely ignores Moose recommendations to use mostly
'ro' attributes and to generate objects on the fly.  HackaMol may be coerced to be more
rigid in future releases.    

Comparing the PhysVecRole within Atom and Molecule may be helpful. For both, the PhysVecRol 
generates a little metadata (mass, name, etc.) and an array of coordinates, forces, and 
charges.  For an atom, the array of coordinates gives an atom (with fixed metadata) the ability 
to store multiple [x,y,z] positions (as a function of time, symmetry, distribution, etc.). What 
is the array of coordinates for Molecule? Usually, the coordinates for a molecule will likely 
remain empty (because the atoms that Molecule contains have the more useful coordinates), but we
can imagine using the coordinates array to track the center of mass of the molecule if needed. 
But for much larger systems, the atoms may be ignored while the Molecule coordinates array could 
be filled with PDLs from Perl Data Language for much faster analyses. I.e. flexible arrays of
coordinates are incredibly powerful. 

=array_method add_$_, all_$_, get_$_, set_$_, count_$_, clear_$_ foreach qw(charges coords forces)

  ARRAY traits, respectively: push, get, set, all, elements, clear
  Descriptions for charges and coords follows.  forces analogous to coords.
  
=array_method add_charges

  push value on to _tcharges array

  $obj->add_charges($_) foreach (0.15, 0.17, 0.14, 0.13);

=array_method all_charges

  returns array of all elements in _tcharges array

    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 0.14 0.13

=array_method get_charges

  return element by index from _tcharges array

    print $obj->get_charges(1); # prints 0.17

=array_method set_charges

  set value of element by index from _tcharges array

    $obj->set_charges(2, 1.00);
    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 1.00 0.13

=array_method count_charges

  return number of elements in _tcharges array
  
    print $obj->count_charges; # prints 4 

=array_method clear_charges
  
  clears _tcharges array
    
    $obj->clear_charges;
    print $_ . " " foreach $obj->all_charges; # does nothing 
    print $obj->count_charges; # prints 0

=array_method add_coords

  push value on to _tcoords array

  $obj->add_coords($_) foreach ([0,0,0],[1,1,1],[-1.0,2.0,-4.0], [3,3,3]);

=array_method all_coords

  returns array of all elements in _tcoords array

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

    my @new_coords = map {[$_->[0]+1,$_->[1]+1,$_->[2]+1]} $obj->all_coords;

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach @new_coords; 

=array_method get_coords

  return element by index from _tcoords array

    printf ("%10.3f %10.3f %10.3f \n", @{$obj->get_coords(1)}); # 

=array_method set_coords

  set value of element by index from _tcoords array

    $obj->set_coords(2, [100,100,100]);
    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

=array_method count_coords

  return number of elements in _tcoords array
  
    print $obj->count_coords; # prints 4 

=array_method clear_coords
  
  clears _tcoords array
    
    $obj->clear_coords;
    print $_ . " " foreach $obj->all_coords; # does nothing 
    print $obj->count_coords # prints 0

=attr name

isa Str that is rw. useful for labeling, bookkeeping...

=attr t 

isa Int or ScalarRef that is rw with default of 0

t is intended to describe the current "setting" of the object. Objects that 
consume the PhysVecRol have  arrays of coordinates, forces, and charges to 
allow storage with the passing of time (hence t) or the generation alternative 
configurations.  For example, a crystal lattice can be stored as a single 
object that consumes PhysVecRole (with associated metadata) along with the 
array of 3d-coordinates resulting from lattice vector translations. 

Experimental: Setting t to a ScalarRef allows all objects to share the same t.  
Although, to use this, a $self->t accessor that dereferences the value would 
seem to be required.  There is nothing, currently, in the core to do so. Not 
sure yet if it is a good or bad idea to do so.

    my $t  = 0;
    my $rt = \$t;  
    $_->t($rt) for (@objects);
    $t = 1; # change t for all objects.

=attr mass 

isa Num that is rw and lazy with a default of 0

=attr xyzfree

isa ArrayRef that is rw and lazy with a default value [1,1,1]. Using this array
allows the object to be fixed for calculations that support it.  For example, to
fix the X and Y coordinates:

  $obj->xyzfree([0,0,1]);

=attr is_fixed

isa Bool that is rw and lazy with a default of 0 (false)

=attr 
charge

is Num that is rw and lazy with a default of (build from t-dependent charges) or 0

why have a charge attribute too? As decribed in ARRAY_ATTRIBUTES, atoms and molecules 
have t-dependent arrays of charges for the purpose of analysis. Often, thinking in 
terms of a given atom/molecule having a "charge" is more intuitive and convenient.  The 
default value of the "charge" attribute is the t index of "charges" or 0.0 if "charges" 
have been cleared before calling on charge.  Being 'rw' charge can be set to whatever,
whenever without effect on the _tcharges array. 

=private_attr _tcharges _tcoords _tforces

isa ArrayRef that is lazy with public ARRAY traits described in ARRAY_METHODS

Gives atoms and molecules t-dependent arrays of charges, coordinates, and forces,
for the purpose of analysis.  e.g. store and analyze atomic charges from a 
quantum mechanical molecule in several intramolecular configurations or a fixed 
configuration in varied external potentials.    

The types of the attributes are left agnostic so that they may be filled with 
whatever the user defines them to be.  Perhaps this is too flexible? One example 
to argue for the flexibility: A molecule consumes PhysVecRole so it has an array 
of coordinates for itself that is likely to remain empty (because the atoms that 
Molecule contains have the more useful coordinates).  For much larger systems, 
the atoms may be ignored while the Molecule t-dependent coordinate array could 
be filled with PDLs from Perl Data Language.  

=cut
