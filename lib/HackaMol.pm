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
use Scalar::Util qw(refaddr);
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
    my $name = join("_", map{ 
                           $_->name.$_->resid  
                            } @atoms[$k , $k+1] 
    );
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
    my $name = join("_", map{ 
                           $_->name.$_->resid  
                            } @atoms[$k .. $k+2] 
    );
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
    my $name = join("_", map{ 
                           $_->name.$_->resid  
                            } @atoms[$k .. $k+3] 
    );
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
  my $self  = shift;
  my $attr  = shift;
  my @atoms = @_; 
 
  my %group;
  foreach my $atom ( @atoms ) {
      push @{ $group{ $atom->$attr } }, $atom;
  }

  my @atomgroups =
    map { HackaMol::AtomGroup->new( atoms => $group{$_} ) } sort
    keys(%group);

  return(@atomgroups);

}

sub find_bonds_brute {
  my $self       = shift;
  my %args       = @_;
  my @bond_atoms = @{ $args{bond_atoms} };
  my @atoms      = @{ $args{candidates} };

  my $fudge      = 0.45; 

  $fudge = $args{fudge} if ( exists($args{fudge}) );

  my @bonds;
  my %name;

  foreach my $at_i (@bond_atoms){
    my $cov_i = $at_i->covalent_radius;
    my $xyz_i = $at_i->xyz;

    foreach my $at_j (@atoms){
        next if (refaddr($at_i) == refaddr($at_j));
        my $cov_j = $at_j->covalent_radius;
        my $dist  = $at_j->distance( $at_i );

        if ($dist <= $cov_i + $cov_j + $fudge){
          my $nm=$at_i->symbol."-".$at_j->symbol;
          $name{$nm}++;
          push @bonds,
            HackaMol::Bond->new(
              name  => "$nm\_".$name{$nm},
              atoms => [ $at_i, $at_j ],
            );
        }

    }
  }
  return (@bonds);
}
 
__PACKAGE__->meta->make_immutable;

1;

__END__

=head1 SYNOPSIS

use HackaMol;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;

my $hack = HackaMol->new( name => "hackitup" );

my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

my $mol = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms] );

# all coordinates from NMR ensemble are loaded

my $fh = $mol->print_xyz( $mol->name . ".xyz" );
foreach my $t ( 1 .. 4 ) {
    $mol->t($t);
    $mol->print_xyz($fh);
}

$mol->t(0);

$mol->translate( -$mol->COM );

$mol->rotate( V( 1, 0, 0 ), 180, V( 10, 10, 10 ) );

$mol->print_xyz($fh);

$mol->translate( -$mol->COM );

$mol->print_xyz($fh);

# translate/rotate method is provided by AtomGroupRole
#populate groups byatom resid attr
my @groups = $hack->group_by_atom_attr( 'resid', $mol->all_atoms );
$mol->push_groups(@groups);

foreach my $ang ( 1 .. 36 ) {
    $_->rotate( V( 1, 1, 1 ), 10, $_->COM ) foreach $mol->all_groups;

    #  $mol->get_groups(1)->print_xyz;
    $mol->print_xyz($fh);
}

$fh->close;

my $radius = 20;
my $natoms = int( 0.0334 * ( $radius**3 ) * 4 * pi / 3 );

my @sphatoms =
  map { HackaMol::Atom->new( Z => 8, charges => [0], coords => [$_] ) }
  map { Math::Vector::Real->random_in_sphere( 3, $radius ) } 1 .. $natoms;

my $sphere = HackaMol::Molecule->new(
    name  => "ball",
    atoms => [@sphatoms]
);

my $bigmol = HackaMol::Molecule->new(
    name  => "bigoverlap",
    atoms => [ $mol->all_atoms, $sphere->all_atoms ],
);

$fh = $bigmol->print_xyz( $bigmol->name . ".xyz" );

foreach my $ang ( 1 .. 36 ) {
    $sphere->rotate( V( 1, 1, 1 ), 20, $sphere->COM );
    $bigmol->print_xyz($fh);
}


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

=method find_bonds_brute 

takes hash argument list and returns bonds.  Find bonds between bond_atoms and the candidates.

  my @oxy_bonds = $hack->find_bonds_brute(
                                    bond_atoms => [$hg],
                                    candidates => [$mol->all_atoms],
                                    fudge      => 0.45,
                  );

fudge is an optional argument. Default is 0.45 (open babel uses same default). find_bonds_brute
uses a bruteforce algorithm that tests the interatomic separation against the sum of the 
covalent radii + fudge. It does not return a self bond for an atom 
( next if refaddr($ati) == refaddr($atj) ).

=head1 SEE ALSO

=for :list
* L<HackaMol::Atom>
* L<HackaMol::Bond>
* L<HackaMol::Angle>
* L<HackaMol::Dihedral>
* L<HackaMol::AtomGroup>
* L<HackaMol::Molecule>
* L<http://perlmol.org>

