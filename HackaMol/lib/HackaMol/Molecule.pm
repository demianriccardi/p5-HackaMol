package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/HackaMol','lib/roles';
use AtomGroup;
use Carp;
use Math::Trig;
use Scalar::Util qw(refaddr);
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecMVRRole',
'BondsAnglesDihedralsRole','QmRole';

extends 'AtomGroup';

has 'atomgroups'  => (
    traits   => ['Array'],
    is       => 'ro',
    isa      => 'ArrayRef[AtomGroup]',
    default  => sub { [] },
    lazy     => 1,  
    handles  => {
        "has_groups"   => 'count',
        "push_groups"   => 'push',
        "get_groups"   => 'get',
        "set_groups"   => 'set',
        "all_groups"   => 'elements',
        "count_groups" => 'count',
        "delete_groups" => 'delete',
        "clear_groups" => 'clear',
    },
);

after 't' => sub {
  my $self = shift;
  $self->gt(@_) if (@_); # set t for all in group
};

sub _build_mass {
  my $self = shift;
  my $mass = 0;
  $mass += $_->mass foreach $self->all_atoms;
  return ($mass); 
}

sub push_groups_by_atom_attr {

  my $self = shift;
  my $attr = shift;

  my %group;
  foreach my $atom ($self->all_atoms){
    push @{$group{$atom->$attr}},$atom;
  }

  my @atomsgroups = map{
                        AtomGroup->new(atoms=>$group{$_})
                       } sort keys (%group);
                      # } sort{$a<=>$b} keys (%group);

  $self->push_groups(@atomsgroups);
  
}

sub dihedral_rotate {
  my $self = shift;
  croak "pass Dihedral, rotation angle (deg), atoms to rotate" unless @_ == 3;
  my $t = $self->t;
  my ($dihe,$ang,$atoms) = @_;
  my ($atom0, $ratom1, $ratom2, $atom3) = $dihe->all_atoms;
  if (2 == scalar( grep{ refaddr($atom0) == refaddr($_) or 
                         refaddr($atom3) == refaddr($_) } @{$atoms} )){ 
    croak "will not rotate atoms on both sides of dihedral " ;
  } 
  my $rvec = ($ratom2->inter_dcoords($ratom1))->versor;
  my @cor  = map{$_->get_coords($t)-$ratom1->xyz} @$atoms; #shift origin too
  my @rcor = $rvec->rotate_3d( deg2rad($ang), @cor );
#shift origin back
  $atoms->[$_]->set_coords($t, $rcor[$_]+$ratom1->xyz) foreach 0 .. $#rcor; 

}


1;

__END__

=head1 SYNOPSIS

use HackaMol::Molecule;

use Math::Vector::Real;

use PDBintoAtoms qw(readinto_atoms); # barebones PDB reader 

my @atoms = readinto_atoms("t/lib/1L2Y.pdb"); 

my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

$mol->translate(-$mol->COM);

$mol->rotate(V(1,0,0), 180, V(10,10,10));

say $mol->count_atoms;

print "\n";

printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

$mol->push_groups_by_atom_attr('resid'); #populate groups by atom resid attr

$_->rotate(V(1,1,1),60,$_->COM,1) foreach $mol->all_groups; # mess up all the amino acids

say $mol->count_atoms;

print "\n";

printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms; 

=head1 DESCRIPTION

The Molecule class provides methods and attributes for collections of atoms that may be divided
into groups, placed into bonds, angles, and dihedrals. The Molecule class extends the AtomGroup 
parent class, which consumes the AtomGroupRole, and consumes PhysVecMVRRole, QmRole, and 
BondsAnglesDihedralsRole. See the documentation of those classes and roles for details.  

In addition to Bonds, Angles, and Dihedrals, which also consume the AtomGroupRole, the Molecule
class has the atomgroups attr.  The atomgroups attr is an ArrayRef[AtomGroup] with native array
traits that allows all the atoms in the Molecule to be grouped and regroup at will. Thus, the 
Molecule class provides a suite of methods and attributes that is very powerful. For example,
a HackaMolX extension for proteins could group the atoms by sidechains and backbones, populate bonds,
and then use Math::Vector::Real objects to sample alternative conformations of the sidechains and 
backbone. 

=array_method push_groups, get_groups, set_groups, all_groups, count_groups, delete_groups, clear_groups

ARRAY traits for the groups attribute, respectively: push, get, set, elements, count, delete, clear

=array_method push_groups

push bond on to groups array

$group->push_groups($bond1, $bond2, @othergroups);

=array_method all_groups

returns array of all elements in groups array

print $_->bond_order, "\n" foreach $group->all_groups; 

=array_method get_groups

return element by index from groups array

print $group->get_groups(1); # returns $bond2 from that pushed above

=array_method set_groups

set groups array by index

$group->set_groups(1, $bond1);

=array_method count_groups

return number of groups in the array  
  
print $group->count_groups; 

=array_method has_groups

same as count_groups, allows clearer conditional code. i.e.  doing something if $mol->has_groups;

=method t 

t is the same attr as before.  Molecule modifies t.  the $mol->t accessor behaves as before.  The $mol->(1)
setter $self->gt(1) to set t for all atoms in the molecule.

=method push_groups_by_atom_attr

takes atom attribute as argument.  pushes the atoms into the atomgroup array by attribute

=method dihedral_rotate

takes Dihedral object, an angle (degress), and active atoms as arguments. rotates the active atoms
about the dihedral and stores rotated coordinates in place ($atom->set_coords($mol->t,$rotated_coor).


=head1 SEE ALSO

=for :list
* L<PhysVecMVRRole>
* L<BondsAnglesDihedralsRole>
* L<PdbRole>
* L<QmRole>
* L<PerlMol>

