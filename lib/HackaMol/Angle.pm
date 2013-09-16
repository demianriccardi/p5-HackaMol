package HackaMol::Angle;
#ABSTRACT: Angle class for HackaMol
use 5.008;
use Moose;
use namespace::autoclean;
use Carp;
use Math::Vector::Real;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' ),'HackaMol::NameRole','HackaMol::AtomGroupRole';

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 0 ,
            lazy    => 1 ,
            clearer => "clear_$_",
            predicate => "has_$_",
          ) foreach qw(ang_fc ang_eq);

sub ang_normvec{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  my $ang  = $self->ang_deg;
  return V(0,0,0) if ($ang == 0 or $ang == 180);
  my $vec1 = $atoms[1]->inter_dcoords($atoms[0]);
  my $vec2 = $atoms[1]->inter_dcoords($atoms[2]);
  my $v1xv2 = $vec1 x $vec2;
  return ($v1xv2->versor);
}

sub ang_deg{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[1]->angle_deg($atoms[0],$atoms[2]));
}

sub ang_rad{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[1]->angle_rad($atoms[0],$atoms[2]));
}

has 'angle_energy_func' => (
    is      => 'rw',
    isa     => 'CodeRef',
    builder => "_build_angle_energy_func",
    lazy    => 1,
);


sub _build_angle_energy_func {
    #my $self = shift; #self is passed by moose, but we don't use it here
    my $subref = sub {
        my $angle = shift;
        my $val = ( $angle->ang_deg - $angle->ang_eq )**2;
        return ($angle->ang_fc*$val);
    };
    return ($subref);
}

sub angle_energy {
    my $self  = shift;
    return (0) unless ($self->ang_fc > 0);
    my $energy = &{$self->angle_energy_func}($self,@_);
    return ($energy);
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

use HackaMol::Atom;
use HackaMol::Angle;

my $atom1 = HackaMol::Atom->new(
    name    => 'O1',
    coords  => [ V( 2.05274, 0.01959, -0.07701 ) ],
    Z       => 8,
);

my $atom2 = HackaMol::Atom->new(
    name    => 'H1',
    coords  => [ V( 1.08388, 0.02164, -0.12303 ) ],
    Z       => 1,
);

my $atom3 = HackaMol::Atom->new(
    name    => 'H2',
    coords  => [ V( 2.33092, 0.06098, -1.00332 ) ],
    Z       => 1,
);

my $angle1 = HackaMol::Angle->new(name=>'OH2', atoms=>[$atom1,$atom2,$atom3]);
my $angle2 = HackaMol::Angle->new(name=>'OH2', atoms=>[$atom2,$atom1,$atom3]);

foreach my $angle ($angle1, $angle2){
  my $pangle = sprintf(
                "Angle: %s, angle: %.2f, vector normal to angle plane: %.5 %.5 %.5 \n",
                $angle->name, 
                $angle->ang, 
                @{$angle->ang_normvec},
                     );

  print $pangle;
}

my @COM_ats = map {HackaMol::Atom->new(
                      name    => "X".$_->name."X",
                      coords  => [ $_->COM ],
                      Z       => 1,
                            )
                  } ($angle1, $angle2);

my @ang_w_HH = grep { $_->get_atoms(0)->Z == 1 and
                      $_->get_atoms(1)->Z == 1} ($angle1, $angle2);


=head1 DESCRIPTION

The HackaMol Angle class provides a set of methods and attributes for working with 
two connections between three atoms.  Like the Bond, the Angle class consumes the 
AtomGroupRole providing methods to determine the center of mass, total charge, etc (see 
AtomGroupRole). An Angle containing (atom1,atom2,atom3) 
produce angles between the atom21 atom23 interatomic vectors.

The Angle class also provides attributes and methods to set force_constants and 
measure energy.  The angle_energy method calls on a CodeRef attribute that the 
user may define.  See descriptions below.  

=attr atoms

isa ArrayRef[Atom] that is lazy with public ARRAY traits provided by the AtomGroupRole (see documentation
for more details).

=attr name

isa Str that is lazy and rw. useful for labeling, bookkeeping...

=attr angle_fc

isa Num that is lazy and rw. default = 0.  force constant for harmonic potentials.

=attr ang_eq

isa Num that is lazy and rw. default = 0.  Equilibrium angle.  The ang method returns angle
in degrees.

=attr angle_energy_func

isa CodeRef that is lazy and rw. default uses builder to generate harmonic potential 
from the angle_fc, ang_eq, and ang.  See the _build_angle_energy_func,
if interested in changing the function form.

=method ang_normvec

no arguments. returns Math::Vector::Real object from the normalized cross product of the atom21 and 
atom23 interatomic vectors.

=method ang_deg 

no arguments. returns the angle (degrees) between the atom21 and atom23 vectors. 

=method ang_rad 

no arguments. returns the angle (radians) between the atom21 and atom23 vectors. 

=method angle_energy

arguments, as many as you want. Calculates energy using the angle_energy_func 
described below, if the attribute, angle_fc > 0.  The angle_energy method calls 
the angle_energy_func as follows: 

my $energy = &{$self->angle_energy_func}($self,@_);

which will pass $self and that in @_ array to angle_energy_func, which, similar to the Bond class, can be redefined.


