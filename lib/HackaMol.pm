package HackaMol;

#ABSTRACT: HackaMol: Object-Oriented Library for Molecular Hacking
use 5.008;
use Moose;
use HackaMol::AtomGroup;
use HackaMol::Molecule;
use HackaMol::Atom;
use HackaMol::Bond;
use HackaMol::Angle;
use HackaMol::Dihedral;
use Carp;

with 'HackaMol::NameRole','HackaMol::MolReadRole';

sub build_bonds {
#take a list of n, atoms; walk down list and generate bonds 
  my $self  = shift;
  my @atoms = @_;
  croak "<2 atoms passed to build_dihedrals" unless (@atoms > 1);
  my @bonds ;

  # build the bonds 
  my $k = 0;
  while ($k+1 <= $#atoms){
    my $name;
    $name .= $_->name.$_->resid foreach (@atoms[$k, $k+1]);
    push @bonds, HackaMol::Bond->new(
                        name => $name,
                        atoms=>[ @atoms[$k, $k+1] ] );
    $k++;
  }
  return (@bonds);
}

sub build_angles {
#take a list of n, atoms; walk down list and generate angles 
  my $self  = shift;
  my @atoms = @_;
  croak "<3 atoms passed to build_dihedrals" unless (@atoms > 2);
  my @angles ;

  # build the angles 
  my $k = 0;
  while ($k+2 <= $#atoms){
    my $name;
    $name .= $_->name.$_->resid foreach (@atoms[$k .. $k+2]);
    push @angles, HackaMol::Angle->new(
                        name => $name,
                        atoms=>[ @atoms[$k .. $k+2] ] );
    $k++;
  }
  return (@angles);
}

sub build_dihedrals {
#take a list of n, atoms; walk down list and generate dihedrals 
  my $self  = shift;
  my @atoms = @_;
  croak "<4 atoms passed to build_dihedrals" unless (@atoms > 3);
  my @dihedrals ;

  # build the dihedrals 
  my $k = 0;
  while ($k+3 <= $#atoms){
    my $name;
    $name .= $_->name.$_->resid foreach (@atoms[$k .. $k+3]);
    push @dihedrals, HackaMol::Dihedral->new(
                        name => $name, 
                        atoms=>[ @atoms[$k .. $k+3] ] );
    $k++;
  }
  return (@dihedrals);
}

sub group_by_atom_attr {
# group atoms by attribute
# Z, name, bond_count, etc. 
  my $self = shift;
  my $attr = shift;
  
  my %group;
  foreach my $atom ( $self->all_atoms ) {
      push @{ $group{ $atom->$attr } }, $atom;
  }

  my @atomgroups =
    map { HackaMol::AtomGroup->new( atoms => $group{$_} ) } sort
    keys(%group);

  return(@atomgroups);

}


__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

use HackaMol;

my $hack  = HackaMol->new(name => "nofunnystuff");

@atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

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

=array_method push_bonds, set_bonds, delete_bonds, clear_bonds

MODIFIED ARRAY traits for the bonds attribute provided by BondsAnglesDihedralsRole

=array_method push_bonds

before push_bonds, bond_count is incremented for all atoms in all bonds to be pushed.

=array_method set_bonds

around set_bonds, bound_count decremented for all atoms in bond being replaced. Then, bond_count is 
incremented for all atoms in new bond

=array_method delete_bonds

before deleting bond, bond_count decremented for all atoms in bond.

=array_method clear_bonds

before clearing bonds, bond_count decremented for all atoms in all bonds.

=method t 

t is the same attr as before.  Molecule modifies t.  the $mol->t accessor behaves as before.  The $mol->(1)
setter $self->gt(1) to set t for all atoms in the molecule.

=method push_groups_by_atom_attr

takes atom attribute as argument.  pushes the atoms into the atomgroup array by attribute

=method all_bonds_atoms  

takes array of atoms as argument, returns array of bonds that includes 1 or more of those atoms

=method all_angles_atoms  

takes array of atoms as argument, returns array of angles that includes 1 or 
more of those atoms

=method all_dihedrals_atoms  

takes array of atoms as argument, returns array of dihedrals that includes 1 or 
more of those atoms 

=method bond_stretch_atoms

takes Bond object, a distance (angstroms, typically), and active atoms as arguments. 
translates the active atoms along the bond_vector by the distance and stores coordinates 
in place ($atom->set_coords($mol->t,$translated_coors)).

=method bond_stretch_groups

takes Bond object, a distance (angstroms, typically), and active groups as arguments. 
translates the atoms in the active groups along the bond_vector by the distance and 
stores coordinates in place.

=method angle_bend_atoms

takes Angle object, an angle (degress), and active atoms as arguments. rotates the active atoms
about the vector normal to be angle and stores rotated coordinates in place 
($atom->set_coords($mol->t,$rotated_coor)).

=method angle_bend_groups

takes Angle object, an angle (degress), and active groups as arguments. rotates the atoms
in the active groups about the vector normal to be angle and stores rotated coordinates 
in place ($atom->set_coords($mol->t,$rotated_coor)).

=method dihedral_rotate_atoms

takes Dihedral object, an angle (degress), and active atoms as arguments. rotates the active atoms
about the dihedral and stores rotated coordinates in place 
($atom->set_coords($mol->t,$rotated_coor)).

=method dihedral_rotate_groups

takes Dihedral object, an angle (degress), and active groups as arguments. rotates atoms in 
groups about the dihedral and stores rotated coordinates in place 
($atom->set_coords($mol->t,$rotated_coor)).

=head1 SEE ALSO

=for :list
* L<PhysVecMVRRole>
* L<BondsAnglesDihedralsRole>
* L<PdbRole>
* L<QmRole>
* L<PerlMol>

