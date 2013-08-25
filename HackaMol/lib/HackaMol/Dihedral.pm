package Dihedral;
#ABSTRACT: Dihedral Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Trig;
use Scalar::Util qw(weaken);
with Storage( 'io' => 'StorableFile' ),'AtomGroupRole';

has 'name' => (
    is   => 'rw',
    isa  => 'Str',
); 

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 0 ,
            lazy    => 1 ,
            clearer => "clear_$_",
          ) foreach qw(dihe_fc dihe_dphase dihe_eq dihe_mult);

sub dihe{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->dihedral($atoms[1],$atoms[2],$atoms[3]));
}

sub dihe_rad{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->dihedral_rad($atoms[1],$atoms[2],$atoms[3]));
}

has 'improper_dihe_energy_func' => (
    is      => 'rw',
    isa     => 'CodeRef',
    builder => "_build_improper_dihe_energy_func",
    lazy    => 1,
);


sub _build_improper_dihe_energy_func {
    my $subref = sub {
        my $dihedral = shift;
        my $val = ( $dihedral->dihe - $dihedral->dihe_eq )**2;
        return ($dihedral->dihe_fc*$val);
    };
    return ($subref);
}

has 'torsion_energy_func' => (
    is      => 'rw',
    isa     => 'CodeRef',
    builder => "_build_torsion_energy_func",
    lazy    => 1,
);


sub _build_torsion_energy_func {
    my $subref = sub {
        my $dihedral = shift;
        my $val = 1 + cos($dihedral->dihe_mult*$dihedral->dihe_rad 
                          - $dihedral->dihe_dphase);
        return ($dihedral->dihe_fc*$val);
    };
    return ($subref);
}

sub torsion_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $energy = &{$self->torsion_energy_func}($self,@_);
    return ($energy);
}

sub improper_dihe_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $energy = &{$self->improper_dihe_energy_func}($self,@_);
    return ($energy);
}

sub BUILD {
    my $self = shift;
    weaken($self);
    $_->push_dihedrals($self) foreach $self->all_atoms;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

use HackaMol::Atom;
use HackaMol::Dihedral;

my ($atom1,$atom4) = map {
                       Atom->new(
                          name    => "C".($_+1),
                          charges => [0],
                          coords  => [ V( $_, $_, 0) ],
                          Z       => 6, 
                       )} (-1, 1);

my ($atom2,$atom3) = map {
                       Atom->new(
                          name    => "S".($_+1),
                          charges => [0],
                          coords  => [ V( $_, 0, 0) ],
                          Z       => 16, 
                       )} (-1, 1);

my $dihe = HackaMol::Dihedral->new(name=>'disulfide', 
                                  atoms=>[$atom1,$atom2,$atom3,$atom4]);

my $pdihe = sprintf(
                "Dihedral: %s, angle: %.2f\n"
                $dihe->name, 
                $dihe->ang, 
                     );
print $pdihe;

my $COM_atom = HackaMol::Atom->new(
                                name    => "X".$_->name."X",
                                coords  => [ $dihe->COM ],
                                Z       => 1,
                                  );


=head1 DESCRIPTION

The HackaMol Dihedral class provides a set of methods and attributes for working 
with three connections between four atoms.  Like the Bond and Angle classes, the 
Dihedral class consumes the AtomGroupRole providing methods to determine the 
center of mass, total charge, etc. (see AtomGroupRole). The Dihedral class is 
flexible.   Instantiation of a Dihedral object also adds that Dihedral to the 
atoms in the dihedral (during the BUILD phase). In contrast, pushing or resetting 
atoms for a Dihedral instance will not add that Dihedral object to the 
atoms. A $dihedral containing (atom1,atom2,atom3,atom4) produces the angle 
($dihedral->dihe) between the planes containing (atom1, atom2, atom3) and 
(atom2, atom3, atom4).

The Dihedral class also provides attributes and methods to set parameters and  
functions to measure energy.  The energy methods call on CodeRef attributes that 
the user may define.  See descriptions below.  

=attr atoms

isa ArrayRef[Atom] that is lazy with public ARRAY traits provided by the 
AtomGroupRole (see documentation for more details).

=attr name

isa Str that is lazy and rw. useful for labeling, bookkeeping...

=attr dihe_dphase

isa Num that is lazy and rw. default = 0.  phase shift for torsion potentials.

=attr dihe_mult

isa Num that is lazy and rw. default = 0.  multiplicity for torsion potentials.

=attr dihe_fc

isa Num that is lazy and rw. default = 0.  force constant for harmonic bond potentials.

=attr dihe_eq

isa Num that is lazy and rw. default = 0.  Equilibrium dihedral angle.  The dihe 
method returns dihedral angle in degrees.

=attr improper_dihe_energy_func torsion_energy_func

isa CodeRef that is lazy and rw. default uses builder to generate a 
harmonic potential for the improper_dihedral and a torsion potential. 

=method dihe 

no arguments. returns the angle between the planes containing (atom1,atom2,atom3) 
and (atom2, atom3, atom4). 

=method improper_dihe_energy torsion_energy

arguments, as many as you want. Calculates energy using the 
improper_dihe_energy_func described below, if the attribute, dihe_fc > 0.  
The improper_dihe_energy method calls the improper_dihe_energy_func as follows: 

my $energy = &{$self->improper_dihe_energy_func}($self,@_);

which will pass $self and that in @_ array to improper_dihe_energy_func, which, similar to the Bond and Angle classes, can be redefined. torsion_energy is analogous.

