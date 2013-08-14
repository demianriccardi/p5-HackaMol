package PhysVecMVR;
# ABSTRACT: Provides the core of HackaMol Atom and Molecule classes.
use Math::Vector::Real;
use Moose::Role;
use Carp;

requires '_build_mass';

has 'name', is => 'rw', isa => 'Str';

has 'mass', is => 'rw', isa => 'Num', lazy => 1, builder => '_build_mass';

has 't', is => 'rw', isa => 'Int|ScalarRef', default => 0;

has 'is_fixed', is => 'rw', isa => 'Bool', lazy => 1, default => 0;

my @t_dep = qw(coords forces);

#DMR notes: tried using coercions where we'd push a number or whatever into 
#       [] to clean up the construction, but there's a pretty big (2x) time hit
#   isa => ArrayRef[Num|ArrayRef|Object]  also slows it down.  although, the benefits may be worth it... 

has "$_" => (
    traits  => ['Array'],
    isa     => 'ArrayRef[Math::Vector::Real]',
    default => sub { [] },
    handles => {
        "push_$_"  => 'push',
        "get_$_"   => 'get',
        "set_$_"   => 'set',
        "all_$_"   => 'elements',
        "clear_$_" => 'clear',
        "count_$_" => 'count',
    },
    lazy   => 1,
) for @t_dep;

has "$_" => (
    traits  => ['Array'],
    isa     => 'ArrayRef[Num]',
    default => sub { [] },
    handles => {
        "push_$_"  => 'push',
        "get_$_"   => 'get',
        "set_$_"   => 'set',
        "all_$_"   => 'elements',
        "clear_$_" => 'clear',
        "count_$_" => 'count',
    },
    lazy   => 1,
) for qw(charges);


has 'units', is => 'rw', isa => 'Str';    #flag for future use [SI]

has 'origin' => (
    is      => 'rw',
    isa     => 'Math::Vector::Real',
    default => sub { V(0,0,0) },
    lazy    => 1,
);

has 'xyzfree' => (
    is      => 'rw',
    isa     => 'ArrayRef[Int]',
    default => sub { [ 1, 1, 1 ] },
    lazy    => 1,
);

#intraobject methods
sub intra_dcharges {
    my $self = shift;
    croak "delta_charges> pass initial and final time" unless ( @_ == 2 );
    my $ti = shift;
    my $tf = shift;
    return($self->get_charges($tf) - $self->get_charges($ti));
}

sub mean_charges {
    my $self = shift;
    my @tcharges = $self->all_charges;
    my $sum = 0;
    $sum += $_ foreach @tcharges;
    return($sum / $self->count_charges);
}

sub msd_charges {
    my $self = shift;
    my $avg  = $self->mean_charges;
    my @tcharges = $self->all_charges;
    my $sum  = 0;
    $sum += ($_ - $avg)**2 foreach @tcharges;
    return ($sum / $self->count_charges);
}

sub intra_dcoords {
#M::V::R makes much simpler
    my $self = shift;
    croak "delta_coords> pass initial and final time" unless ( @_ == 2 );
    my $ti = shift;
    my $tf = shift;
    return( $self->get_coords($tf) - $self->get_coords($ti) );
}

sub mean_coords {
    my $self    = shift;
    my @tcoords = $self->all_coords;
    my $sum = V(0, 0, 0);
    $sum += $_ foreach @tcoords;
    return ($sum / $self->count_coords); 
}

sub msd_coords {
# returns scalar
    my $self    = shift;
    my $avg     = $self->mean_coords;
    my @tcoords = $self->all_coords;
    my $sum = 0;
    foreach my $c (@tcoords){
      my $dc = $c - $avg;
      $sum += $dc*$dc;
    }
    return ($sum / $self->count_coords);
}

sub intra_dforces {
#M::V::R makes much simpler
    my $self = shift;
    croak "delta_forces> pass initial and final time" unless ( @_ == 2 );
    my $ti = shift;
    my $tf = shift;
    return( $self->get_forces($tf) - $self->get_forces($ti) );
}

sub mean_forces {
    my $self    = shift;
    my @tforces = $self->all_forces;
    my $sum = V(0, 0, 0);
    $sum += $_ foreach @tforces;
    return ($sum / $self->count_forces);
}

sub msd_forces {
# returns scalar
    my $self    = shift;
    my $avg     = $self->mean_forces;
    my @tforces = $self->all_forces;
    my $sum = 0;
    foreach my $c (@tforces){
      my $dc = $c - $avg;
      $sum += $dc*$dc;
    }
    return ($sum / $self->count_forces);
}

#interobject methods
sub distance {
    my $self = shift;
    my $obj2 = shift or croak "need to pass another obj that does PhysVec";
    my ($ts, $t2) = ($self->t, $obj2->t);
    carp "comparing objects with different times" unless($ts == $t2);
    my $vs = $self->get_coords($ts); 
    my $v2 = $obj2->get_coords($t2); 
    return ( $vs->dist( $v2 ) );
}

sub inter_dcharges {      
    my $self = shift;
    my $obj2 = shift or croak "need to pass another obj that does PhysVec";
    my ($ts, $t2) = ($self->t, $obj2->t);
    carp "comparing objects with different times" unless($ts == $t2);
    my $dvec = $obj2->get_charges($t2) - $self->get_charges($ts);
    return ($dvec);
}

sub inter_dcoords { 
    my $self = shift;
    my $obj2 = shift or croak "need to pass another obj that does PhysVec";
    my ($ts, $t2) = ($self->t, $obj2->t);
    carp "comparing objects with different times" unless($ts == $t2);
    my $dvec = $obj2->get_coords($t2) - $self->get_coords($ts);
    return ($dvec);
}

sub inter_dforces {      
    my $self = shift;
    my $obj2 = shift or croak "need to pass another obj that does PhysVec";
    my ($ts, $t2) = ($self->t, $obj2->t);
    carp "comparing objects with different times" unless($ts == $t2);
    my $dvec = $obj2->get_forces($t2) - $self->get_forces($ts);
    return ($dvec);
}

sub charge {
  my $self = shift;
  carp "charge> takes no arguments. returns get_charges(t)\n" if (@_);
  return($self->get_charges($self->t));
}

sub xyz {
  my $self = shift;
  carp "xyz> takes no arguments. returns get_coords(t)\n" if (@_);
  return($self->get_coords($self->t));
}

sub force {
  my $self = shift;
  carp "force> takes no arguments. returns get_forces(t)\n" if (@_);
  return($self->get_forces($self->t));
}

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

# instance of class that consumes the PhysVecRol

my $obj = Class_with_PhysVec->new( name => 'foo', t => 0 , charges => [0.1], coords => [ V(0,1,2) ]); 

# add some charges

$obj->push_charges($_) foreach ( 0.3, 0.2, 0.1, -0.1, -0.2, -0.36 );

my $sum_charges = 0;

$sum_charges += $_ foreach $obj->all_charges;
print "average charge: ", $sum_charges / $obj->count_charges;

# add some coordinates

$obj->push_coords($_) foreach ( map{V(@{$_})} foreach ([0,0,0], [1,1,1], [-1.0,2.0,-4.0]) );

print $obj->mean_charges . "\n";
print $obj->msd_charges . "\n";
printf ("%10.3f %10.3f %10.3f \n", @{$obj->mean_coords};
print $obj->msd_coords . "\n";


=head1 DESCRIPTION

PhysVec provides the core attributes and methods shared between Atom 
and Molecule classes. Consuming this role gives Classes a place to store 
coordinates, forces, and charges, perhaps, over the course of a simulation 
or for a collection 
configurations for which all other object metadata (name, mass, etc) remains 
fixed. As such, the 't' attribute, described below, is important to understand. 
The PhysVec is very flexible. Several attribute types are left agnostic and all
are rw so that they may be filled with whatever the user defines them to be... on 
the fly.  This seems most intuitive from the perspective of carrying out computational 
work on molecules.  Thus, HackaMol bravely ignores Moose recommendations to use mostly
'ro' attributes and to generate objects on the fly.  HackaMol may be coerced to be more
rigid in future releases.    

Comparing the PhysVec within Atom and Molecule may be helpful. For both, the PhysVecRol 
generates a little metadata (mass, name, etc.) and an array of coordinates, forces, and 
charges.  For an atom, the array of coordinates gives an atom (with fixed metadata) the ability 
to store multiple [x,y,z] positions (as a function of time, symmetry, distribution, etc.). What 
is the array of coordinates for Molecule? Usually, the coordinates for a molecule will likely 
remain empty (because the atoms that Molecule contains have the more useful coordinates), but we
can imagine using the coordinates array to track the center of mass of the molecule if needed. 
But for much larger systems, the atoms may be ignored while the Molecule coordinates array could 
be filled with PDLs from Perl Data Language for much faster analyses. I.e. flexible arrays of
coordinates are incredibly powerful. 

=array_method push_$_, all_$_, get_$_, set_$_, count_$_, clear_$_ foreach qw(charges coords forces)

  ARRAY traits, respectively: push, get, set, all, elements, clear
  Descriptions for charges and coords follows.  forces analogous to coords.
  
=array_method push_charges

  push value on to charges array

  $obj->push_charges($_) foreach (0.15, 0.17, 0.14, 0.13);

=array_method all_charges

  returns array of all elements in charges array

    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 0.14 0.13

=array_method get_charges

  return element by index from charges array

    print $obj->get_charges(1); # prints 0.17

=array_method set_charges

  set value of element by index from charges array

    $obj->set_charges(2, 1.00);
    print $_ . " " foreach $obj->all_charges; # prints 0.15 0.17 1.00 0.13

=array_method count_charges

  return number of elements in charges array
  
    print $obj->count_charges; # prints 4 

=array_method clear_charges
  
  clears charges array
    
    $obj->clear_charges;
    print $_ . " " foreach $obj->all_charges; # does nothing 
    print $obj->count_charges; # prints 0

=array_method push_coords

  push value on to coords array

  $obj->push_coords($_) foreach ([0,0,0],[1,1,1],[-1.0,2.0,-4.0], [3,3,3]);

=array_method all_coords

  returns array of all elements in coords array

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

    my @new_coords = map {[$_->[0]+1,$_->[1]+1,$_->[2]+1]} $obj->all_coords;

    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach @new_coords; 

=array_method get_coords

  return element by index from coords array

    printf ("%10.3f %10.3f %10.3f \n", @{$obj->get_coords(1)}); # 

=array_method set_coords

  set value of element by index from coords array

    $obj->set_coords(2, [100,100,100]);
    printf ("%10.3f %10.3f %10.3f \n", @{$_}) foreach $obj->all_coords; 

=array_method count_coords

  return number of elements in coords array
  
    print $obj->count_coords; # prints 4 

=array_method clear_coords
  
  clears coords array
    
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
object that consumes PhysVec (with associated metadata) along with the 
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

=attr charges coords forces

isa ArrayRef that is lazy with public ARRAY traits described in ARRAY_METHODS

Gives atoms and molecules t-dependent arrays of charges, coordinates, and forces,
for the purpose of analysis.  e.g. store and analyze atomic charges from a 
quantum mechanical molecule in several intramolecular configurations or a fixed 
configuration in varied external potentials.    

The types of the attributes are left agnostic so that they may be filled with 
whatever the user defines them to be.  Perhaps this is too flexible? One example 
to argue for the flexibility: A molecule consumes PhysVec so it has an array 
of coordinates for itself that is likely to remain empty (because the atoms that 
Molecule contains have the more useful coordinates).  For much larger systems, 
the atoms may be ignored while the Molecule t-dependent coordinate array could 
be filled with PDLs from Perl Data Language.  

=attr distance_coderef

isa CodeRef that is rw and lazy with default provided by
_builder_distance_coderef. Default takes two arguments:

  &{$self->distance_coderef}($obj1,$obj2);

Wrapped in the method "distance" for cleaner interface. 
see _build_distance_coderef that is default

=method distance

Takes one argument ($obj) and passes $self and 
$obj to the "distance_coderef". I.e.:

  $obj1->distance($obj2);

  which does this inside: &{$self->distance_coderef}($self,$obj2)

define distance_coderef as you wish.  See default generated from _build_distance_coderef for an example.

=private_attr _build_distance_coderef

builds the default distance_coderef attribute that is wrapped in the distance
method.

=attr delta_charges_coderef

isa CodeRef that is rw and lazy with default provided by
_builder_delta_charges_coderef. Default takes three arguments:

  &{$self->delta_charges_coderef}($self,$ti,$tf);

The delta_charges_coderef, delta_coords_coderef, and delta_forces_coderef are
all analogous (they can probably be abstracted into other roles), 
and are similar to distance_coderef.  In contrast to distance_coderef, which
takes two objects as arguments, delta_charges_coderef subtracts the initial
charges at $ti from that at the final $tf.  Wrapped in "delta_charges" method
for a cleaner interface. The delta_charges_coderef default behaviour is built
from _build_charges_coderef, and it can be adjusted on the fly
to redefine "delta_charges".  You are free to define new Charge classes and work
them into your objects. 

=private_attr _build_delta_charges_coderef

builds the default delta_charges_coderef attribute that is wrapped in the
delta_charges method

=method delta_charges 

Takes initial t ($ti) and final t ($tf) arguments, and passes $self,
$ti, $tf to the delta_charges_coderef. I.e.:

  $obj1->delta_charges($ti,$tf);

  which does this inside: &{$self->delta_charges_coderef}($self,$ti,$tf)

  see delta_charges_coderef and delta_distance_coderef and the builders for
  further discussions

=attr delta_coords_coderef delta_forces_coderef

isa CodeRef that is rw and lazy with default provided by
default behaviour of delta_coords_coderef and delta_forces_coderef are exactly
the same.  Default takes three arguments:

  &{$self->delta_coords_coderef}($self,$ti,$tf);
  &{$self->delta_forces_coderef}($self,$ti,$tf);

See delta_charges_coderef for additional, analogous discussions.
_build_delta_coords_coderef and _build_delta_forces_coderef are default builders

=private_attr _build_delta_coords_coderef _build_delta_forces_coderef

builds the default delta_coords_coderef and delta_forces_coderef attribute that is wrapped in the
delta_coords and delta_forces methods.

_build_delta_coords is the same with s/forces/coords/ .

=method delta_coords delta_forces 

Takes initial t ($ti) and final t ($tf) arguments, and passes $self,
$ti, $tf to the delta_charges_coderef. I.e.:

  $obj1->delta_coords($ti,$tf);
  
  which does this inside: &{$self->delta_coords_coderef}($self,$ti,$tf)

  $obj1->delta_forces($ti,$tf);

  which does this inside: &{$self->delta_forces_coderef}($self,$ti,$tf)

  see delta_charges_coderef, delta_coords_coderef and delta_distance_coderef and the builders for
  further discussions.

=method charge

called with no arguments.  returns $self->get_charges($self->t);

=method xyz

called with no arguments.  returns $self->get_coords($self->t);

=method force

called with no arguments.  returns $self->get_forces($self->t);

=cut
