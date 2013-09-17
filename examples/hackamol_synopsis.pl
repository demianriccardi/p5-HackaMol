use Modern::Perl;
use HackaMol;
use Math::Vector::Real;
use Math::Vector::Real::Random;
use Math::Trig;

my $hack = HackaMol->new( name => "hackitup" );

my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

my $mol = HackaMol::Molecule->new( name => 'trp-cage', atoms => [@atoms] );

$mol->translate( -$mol->COM );

$mol->rotate( V( 1, 0, 0 ), 180, V( 10, 10, 10 ) );

$mol->print_xyz;

#populate groups byatom resid attr 
my @groups = $hack->group_by_atom_attr('resid');    
$mol->push_groups(@groups);

# mess upall the amino acids
$_->rotate( V( 1, 1, 1 ),  60, $_->COM, 1 ) foreach $mol->all_groups;   

$mol->gt(1);
$mol->print_xyz;

$_->rotate( V( 1, 1, 1 ), -60, $_->COM, 1 ) foreach $mol->all_groups;   
$mol->print_xyz;

my $radius = 30;
my $natoms = int( 0.0334 * ( $radius**3 ) * 4 * pi / 3 );

my @sphatoms =
  map { HackaMol::Atom->new( Z => 8, charges => [0], coords => [$_] ) }
  map { Math::Vector::Real->random_in_sphere( 3, $radius ) } 1 .. $natoms;

my $sphere = HackaMol::Molecule->new(
    name  => "ball",
    atoms => [@sphatoms]
);

$sphere->print_xyz;

my $bigmol = HackaMol::Molecule->new(
    name  => "bigoverlap",
    atoms => [ $mol->all_atoms, $sphere->all_atoms ],
);

$bigmol->print_xyz;

