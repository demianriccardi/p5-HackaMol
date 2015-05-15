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
use namespace::autoclean;
use MooseX::StrictConstructor;
use Scalar::Util qw(refaddr);
use Carp;

with 'HackaMol::Roles::NameRole', 
     'HackaMol::Roles::MolReadRole', 
     'HackaMol::Roles::PathRole',
     'HackaMol::Roles::ExeRole', 
     'HackaMol::Roles::FileFetchRole';

sub pdbid_mol {
  my $self  = shift;
  my $pdbid = shift || croak "pass pdbid";
  my ($file,$rc) = $self->getstore_pdbid($pdbid);
  return($self->read_file_mol($file));
}

sub read_file_push_coords_mol{
    my $self = shift;
    my $file = shift;
    my $mol  = shift or croak "must pass molecule to add coords to";

    my @atoms = $self->read_file_atoms($file);
    my @matoms= $mol->all_atoms;

    if (scalar(@matoms) != scalar(@atoms) ){
      croak "file_push_coords_mol> number of atoms not same";
    }
    foreach my $i (0 .. $#atoms) {
      if ($matoms[$i]->Z != $atoms[$i]->Z){
        croak "file_push_coords_mol> atom mismatch";
      }
      $matoms[$i]->push_coords($_) foreach ($atoms[$i]->all_coords);
    }
}

sub read_file_mol{
    my $self = shift;
    my $file = shift;

    my @atoms = $self->read_file_atoms($file);
    my $name = $file . ".mol";
    return (HackaMol::Molecule->new(name=>$name, atoms=>[@atoms]));
}

sub build_bonds {

    #take a list of n, atoms; walk down list and generate bonds
    my $self  = shift;
    my @atoms = @_;
    croak "<2 atoms passed to build_bonds" unless ( @atoms > 1 );
    my @bonds;
    # build the bonds
    my $k = 0;
    while ( $k + 1 <= $#atoms ) {
        my $name =
              join( "_", map { _name_resid($_,'B')} @atoms[ $k .. $k + 1 ] );
        push @bonds,
          HackaMol::Bond->new(
            name  => $name,
            atoms => [ @atoms[ $k, $k + 1 ] ]
          );
        $k++;
    }
    return (@bonds);
}

sub build_angles {

    #take a list of n, atoms; walk down list and generate angles
    my $self  = shift;
    my @atoms = @_;
    croak "<3 atoms passed to build_angles" unless ( @atoms > 2 );
    my @angles;

    # build the angles
    my $k = 0;

    while ( $k + 2 <= $#atoms ) {
        my $name =
              join( "_", map { _name_resid($_,'A')} @atoms[ $k .. $k + 2 ] );
        push @angles,
          HackaMol::Angle->new(
            name  => $name,
            atoms => [ @atoms[ $k .. $k + 2 ] ]
          );
        $k++;
    }
    return (@angles);
}

sub _name_resid {
  my $atom    = shift;
  my $default = shift;
  return ($default    . $atom->resid) unless $atom->has_name; # this comes up when undefined atoms are passed 
  return ($atom->name . $atom->resid);
}

sub build_dihedrals {

    #take a list of n, atoms; walk down list and generate dihedrals
    my $self  = shift;
    my @atoms = @_;
    croak "<4 atoms passed to build_dihedrals" unless ( @atoms > 3 );
    my @dihedrals;

    # build the dihedrals
    my $k = 0;
    while ( $k + 3 <= $#atoms ) {
        my $name =
              join( "_", map { _name_resid($_,'D')} @atoms[ $k .. $k + 3 ] );
        push @dihedrals,
          HackaMol::Dihedral->new(
            name  => $name,
            atoms => [ @atoms[ $k .. $k + 3 ] ]
          );
        $k++;
    }
    return (@dihedrals);
}

sub group_by_atom_attr {

    # group atoms by attribute
    # Z, name, bond_count, etc.
    my $self  = shift;
    my $attr  = shift;
    my @atoms = @_;

    my %group;
    foreach my $atom (@atoms) {
        push @{ $group{ $atom->$attr } }, $atom;
    }

    my @atomgroups =
      map { HackaMol::AtomGroup->new( atoms => $group{$_} ) } sort
      keys(%group);

    return (@atomgroups);

}

sub find_disulfide_bonds {
    my $self  = shift;
    my @sulf  = grep {$_->Z == 16} @_;
    my @ss = $self->find_bonds_brute(
                                     bond_atoms => [@sulf],
                                     candidates => [@sulf],
                                     fudge      => 0.45,
                                     max_bonds  => 1,
                                    );
    return @ss;
}

sub find_bonds_brute {
    my $self       = shift;
    my %args       = @_;
    my @bond_atoms = @{ $args{bond_atoms} };
    my @atoms      = @{ $args{candidates} };

    my $fudge = 0.45;
    my $max_bonds = 99;

    $fudge     = $args{fudge}     if ( exists( $args{fudge} ) );
    $max_bonds = $args{max_bonds} if ( exists( $args{max_bonds} ) );

    my @init_bond_counts = map {$_->bond_count} (@bond_atoms,@atoms);

    my @bonds;
    my %name;

    foreach my $at_i (@bond_atoms) {
        next if ($at_i->bond_count >= $max_bonds);
        my $cov_i = $at_i->covalent_radius;
        my $xyz_i = $at_i->xyz;

        foreach my $at_j (@atoms) {
            next if ( refaddr($at_i) == refaddr($at_j) );
            next if ($at_j->bond_count >= $max_bonds);
            my $cov_j = $at_j->covalent_radius;
            my $dist  = $at_j->distance($at_i);

            if ( $dist <= $cov_i + $cov_j + $fudge ) {
                my $nm = $at_i->symbol . "-" . $at_j->symbol;
                $name{$nm}++;
                push @bonds,
                  HackaMol::Bond->new(
                    name  => "$nm\_" . $name{$nm},
                    atoms => [ $at_i, $at_j ],
                  );
                $at_i->inc_bond_count;
                $at_j->inc_bond_count;
            }

        }
    }
  
    my $i = 0;
    foreach my $at (@bond_atoms,@atoms){
      $at->reset_bond_count;
      $at->inc_bond_count($init_bond_counts[$i]);
      $i++;
    }
    return (@bonds);

}

__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS 

       # simple example: load pdb file and extract the disulfide bonds

       use HackaMol;

       my $bldr = HackaMol->new( name => 'builder');
       my $mol  = $bldr->pdbid_mol('1kni');

       my @disulfide_bonds = $bldr->find_disulfide_bonds( $mol->all_atoms );

       print $_->dump foreach @disulfide_bonds;

See the above executed in this L<linked notebook | https://github.com/demianriccardi/p5-IPerl-Notebooks/blob/master/HackaMol/HackaMol-Synopsis.ipynb>

=head1 DESCRIPTION

The L<HackaMol publication | http://pubs.acs.org/doi/abs/10.1021/ci500359e> has
a more complete description of the library (L<pdf available from researchgate| http://www.researchgate.net/profile/Demian_Riccardi/publication/273778191_HackaMol_an_object-oriented_Modern_Perl_library_for_molecular_hacking_on_multiple_scales/links/550ebec60cf27526109e6ade.pdf>). 

Citation: J. Chem. Inf. Model., 2015, 55, 721 

Loading the HackaMol library in a script with 
 
       use HackaMol;

provides attributes and methods of a builder class. It also loads all the 
classes provided by the core so including them is not necessary, e.g.:
 
       use HackaMol::Atom;
       use HackaMol::Bond;
       use HackaMol::Angle;
       use HackaMol::Dihedral;
       use HackaMol::AtomGroup;
       use HackaMol::Molecule;

The methods, described below, facilitate the creation of objects from files and
other objects.

=attr name 

name is a rw str provided by HackaMol::NameRole.

=method pdbid_mol

one argument: pdbid

This method will download the pdb, unless it exists, and load it into a 
HackaMol::Molecule object. For example,

      my $mol = HackaMol->new->pdbid_mol('2cba');

=method read_file_atoms

one argument: filename.

This method parses the file (e.g. file.xyz, file.pdb) and returns an 
array of HackaMol::Atom objects.

=method read_file_mol

one argument: filename.

This method parses the file (e.g. file.xyz, file.pdb) and returns a 
HackaMol::Molecule object.

=method read_file_push_coords_mol

two arguments: filename and a HackaMol::Molecule object. 

This method reads the coordinates from a file and pushes them into the atoms 
contained in the molecule. Thus, the atoms in the molecule and the atoms in 
the file must be the same.

=method build_bonds

takes a list of atoms and returns a list of bonds.  The bonds are generated for
"list neighbors" by simply stepping through the atom list one at a time. e.g.

  my @bonds = $hack->build_bonds(@atoms[1,3,5]);

will return two bonds: B13 and B35 

=method build_angles

takes a list of atoms and returns a list of angles. The angles are generated 
analagously to build_bonds, e.g.

  my @angles = $hack->build_angles(@atoms[1,3,5]);

will return one angle: A135

=method build_dihedrals

takes a list of atoms and returns a list of dihedrals. The dihedrals are generated 
analagously to build_bonds, e.g.

  my @dihedral = $hack->build_dihedrals(@atoms[1,3,5]);

will croak!  you need atleast four atoms.

  my @dihedral = $hack->build_dihedrals(@atoms[1,3,5,6,9]);

will return two dihedrals: D1356 and D3569

=method group_by_atom_attr

arguments are an atom attribute and then a list of atoms. 

This method builds AtomGroup objects that are grouped by attribute.

=method find_bonds_brute 

The arguments are key_value pairs of bonding criteria (see example below). 

This method returns bonds between bond_atoms and the candidates using the 
criteria (many of wich have defaults).

  my @oxy_bonds = $hack->find_bonds_brute(
                                    bond_atoms => [$hg],
                                    candidates => [$mol->all_atoms],
                                    fudge      => 0.45,
                                    max_bonds  => 6,
  );

fudge is optional with Default is 0.45 (open babel uses same default); 
max_bonds is optional with default of 99. max_bonds is compared against
the atom bond count, which are incremented during the search. Before returning
the bonds, the bond_count are returned the values before the search.  For now,
molecules are responsible for setting the number of bonds in atoms. 
find_bonds_brute uses a bruteforce algorithm that tests the interatomic 
separation against the sum of the covalent radii + fudge. It will not test
for bond between atoms if either atom has >= max_bonds. It does not return 
a self bond for an atom (C< next if refaddr($ati) == refaddr($atj) >).

=method find_disulfide_bonds

the argument is a list of atoms, e.g. '($mol->all_atoms)'. 

this method returns disulfide bonds as bond objects.

=head1 SEE ALSO

=for :list
* L<HackaMol::Atom>
* L<HackaMol::Bond>
* L<HackaMol::Angle>
* L<HackaMol::Dihedral>
* L<HackaMol::AtomGroup>
* L<HackaMol::Molecule>
* L<HackaMol::X::Calculator>
* L<HackaMol::X::Vina>
* L<Protein Data Bank | http://pdb.org>
* L<VMD | http://www.ks.uiuc.edu/Research/vmd/>

