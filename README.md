HackaMol
========
Object-Oriented Perl 5, Moose Library for Molecular Hacking on multiple scales

VERSION
========
       version 0.00_01

SYNOPSIS
========
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


DESCRIPTION
============

       The HackaMol library aims to enable users to build simple, yet powerful
       scripts for carrying out computational work on molecules at multiple
       scales. The molecular object system organizes atoms within molecules
       using groups, bonds, angles, and dihedrals.  HackaMol seeks to provide
       intuitive attributes and methods that may be harnessed to coerce
       computational chemistry through a common core.

The library is inspired by PerlMol, BioPerl, MMTSB (a Perl library developed
by Michael Feig and coworkers), and my own experiences as a researcher.
