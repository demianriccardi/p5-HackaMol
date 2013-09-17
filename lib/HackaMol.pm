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

my $hack  = HackaMol->new(name => "hackitup");

@atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

$mol->translate(-$mol->COM);

$mol->rotate(V(1,0,0), 180, V(10,10,10));

say $mol->count_atoms;

print "\n";

printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

my @groups = $hack->group_by_atom_attr('resid'); #populate groups by atom resid attr
$mol->push_groups(@groups); 

$_->rotate(V(1,1,1),60,$_->COM,1) foreach $mol->all_groups; # mess up all the amino acids

print $mol->count_atoms . "\n";

print "\n";

printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms; 

=head1 DESCRIPTION

The HackaMol library enables users to build simple, yet powerful scripts for carrying out 
computational work on molecules at multiple scales. The molecular object system organizes 
atoms within molecules using groups, bonds, angles, and dihedrals.  HackaMol seeks to provide 
intuitive attributes and methods that may be harnessed to coerce computational chemistry 
through a common core. 

The HackaMol class uses the core classes to provide some object building
utilities described below.  This class consumes HackaMol::MolReadRole to provide
structure readers for xyz and pdb coordinates.  Additional formats are pretty
easy to add, but using open babel to do so may be a more robust approach.

=attr name 

name is a rw str provided by HackaMol::NameRole.

=method build_bonds

takes a list of atoms and returns a list of bonds.  The bonds are generated for
"list neighbors" by simply stepping through the atom list one at a time. e.g.

  my @bonds = $hack->build_bonds(@atoms[1,3,5]);

  will return two bonds: B13 and B35 

=method build_angles

takes a list of atoms and returns a list of angles. The angles are generated for
"list neighbors" by simply stepping through the atom list one at a time. e.g.

  my @angles = $hack->build_angles(@atoms[1,3,5]);

  will return one angle: A135

=method build_dihedrals

takes a list of atoms and returns a list of dihedrals. The dihedrals are generated for
"list neighbors" by simply stepping through the atom list one at a time. e.g.

  my @dihedral = $hack->build_dihedrals(@atoms[1,3,5]);

  will croak!  you need atleast four atoms.

  my @dihedral = $hack->build_dihedrals(@atoms[1,3,5,6,9]);

  will return two dihedrals: D1356 and D3569

=method group_by_atom_attr

takes atom attribute as argument and builds AtomGroup objects by attribute

=head1 SEE ALSO

=for :list
* L<HackaMol::Atom>
* L<HackaMol::Bond>
* L<HackaMol::Angle>
* L<HackaMol::Dihedral>
* L<HackaMol::AtomGroup>
* L<HackaMol::Molecule>
* L<PerlMol>

