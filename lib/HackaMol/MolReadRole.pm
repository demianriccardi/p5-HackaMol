package HackaMol::MolReadRole;
# ABSTRACT: Read XYZ and PDB files 
use Moose::Role;
use Carp;
use Math::Vector::Real;
use FileHandle;

sub read_file_atoms {
  my $self = shift;
  my $file = shift;
  my @atoms;

  if ($file =~ m/\.pdb$/){
    @atoms = $self->read_pdb_atoms($file);
  }
  elsif ($file =~ m/\.xyz$/){
    @atoms = $self->read_xyz_atoms($file);
  }
  else {
    croak "$file format not supported";
  }
  return (@atoms); 
}

sub read_pdb_atoms {
#read pdb file and generate list of Atom objects
  my $self  = shift;
  my $file  = shift;
  my $segid = $file;
  $segid =~ s/\.pdb//;
  $segid =~ s/t\/lib\///;
  my $fh = FileHandle->new("<$file");

  my @atoms;
  my ( $n, $t ) = ( 0, 0 );

  while (<$fh>) {

      if (/^(?:MODEL\s+(\d+))/) {
          $t = $1 - 1;
          $n = 0;
      }
      elsif (/^(?:HETATM|ATOM)/) {

          my (
              $record_name, $serial, $name, $altloc,  $resName,
              $chainID,     $resSeq, $icod, $x,       $y,
              $z,           $occ,    $B,    $element, $charge
          ) = unpack "A6A5x1A4A1A3x1A1A4A1x3A8A8A8A6A6x10A2A2", $_;

          if   ( $charge =~ m/\d/ ) { $charge = _qstring_num($charge) }
          else                      { $charge = 0 }

          if   ( $chainID =~ m/\w/ ) { $chainID = uc( _trim($chainID) ) }
          else                       { $chainID = 'AA' }

          $element = ucfirst( lc( _trim($element) ) );
          $name    = _trim($name);
          $resName = _trim($resName);
          $resSeq  = _trim($resSeq);
          $resSeq  = 0 if ( $resSeq < 0 );
          $serial  = _trim($serial);
          my $xyz = V( $x, $y, $z );

          if ( $t == 0 ) {
              $atoms[$n] = HackaMol::Atom->new(
                  name        => $name,
                  record_name => $record_name,
                  serial      => $serial,
                  chain       => $chainID,
                  symbol      => $element,
                  charges     => [$charge],
                  coords      => [$xyz],
                  occ         => $occ * 1,
                  bfact       => $B * 1,
                  resname     => $resName,
                  resid       => $resSeq,
                  altloc      => $altloc,
              );
          }
          else {
              croak "atoms have changed from last model to current: $t\n"
                if ( 
                    $name    ne $atoms[$n]->name   or 
                    $element ne $atoms[$n]->symbol
                    );

              $atoms[$n]->set_coords( $t, $xyz );
          }
          $n++;
      }
  }

  # set iatom to track the array.  diff from serial which refers to pdb
  $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
  return (@atoms);
}

sub read_xyz_atoms {
#read pdb file and generate list of Atom objects
  my $self  = shift;
  my $file  = shift;
  my $segid = $file;
  $segid =~ s/\.xyz//;
  $segid =~ s/t\/lib\///;
  my $fh = FileHandle->new("<$file");
use Data::Dumper;
  my @atoms;
  my ( $n, $t ) = ( 0, 0 );

  my $nat  = undef;
  while (<$fh>) {

    if (/^(\s*\d+\s*)$/){
      $n = (split)[0];
      if (defined($nat)){
        croak "number of atoms has changed\n" unless ($nat == $n);
        $t++;
      }
      $nat  = $n;
      $n = 0;
    }
    elsif (/(\w+|\d+)(\s+-*\d+\.\d+){3}/) {
      my @stuff = split;
      my $sym = $stuff[0];
      my $xyz = V( @stuff[1,2,3] );
      if ( $t == 0 ) {
        if ($sym =~ /\d/)
             {
         $atoms[$n] = HackaMol::Atom->new( Z      => $sym,coords => [$xyz])
        }
        else {
         $atoms[$n] = HackaMol::Atom->new( symbol => $sym,coords => [$xyz])
        }
      }
      else {
        if ($sym =~ /\d/){
          croak "atoms have changed from last model to current: $t\n"
                if ( $sym != $atoms[$n]->Z );
        }
        else {
          croak "atoms have changed from last model to current: $t\n"
                if ( $sym ne $atoms[$n]->symbol );
        }
        $atoms[$n]->set_coords( $t, $xyz );

      }
      $n++;
    }
  }
  # set iatom to track the array.  diff from serial which refers to pdb
  $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
  return (@atoms);
}                                                    

sub _trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub _qstring_num {

    # _qstring something like 2+  or 2-
    my $string = shift;
    $string =~ s/\+//;
    $string =~ s/(.*?)(\-)/$2$1/;
    $string = sprintf( "%g", $string );
    return $string;

}

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

use HackaMol;

my @atoms = read_file_atoms("t/lib/1L2Y.pdb"); 

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

