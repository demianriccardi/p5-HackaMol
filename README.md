HackaMol
========
Object-oriented Perl 5, Moose library for molecular hacking on multiple scales

VERSION
========
developer version 0.00_12 
Available for testing from cpan.org:
       
please see *[HackaMol on MetaCPAN](https://metacpan.org/release/DEMIAN/HackaMol-0.00_12) for formatted documentation.  
       
SYNOPSIS
========
       use HackaMol;
       use Math::Vector::Real;
       use Math::Vector::Real::Random;
       use Math::Trig;
       
       my $hack = HackaMol->new( name => "hackitup" );
       my @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");
       
       # all coordinates from NMR ensemble are loaded into atoms
       my $mol = HackaMol::Molecule->new(
           name  => 'trp-cage',
           atoms => [@atoms]
       );
       
       #recenter all coordinates to center of mass
       foreach my $t ( 0 .. $atoms[0]->count_coords - 1 ) {
           $mol->t($t);
           $mol->translate( -$mol->COM );
       }
       
       #create array of all alanine, CA atoms

       my @CAs = grep {$_->name eq 'CA'} 
                 grep {$_->resname eq 'ALA'} $mol->all_atoms;
       
       # print coordinates from t=0 to trp-cage.xyz and return filehandle
       my $fh = $mol->print_xyz( $mol->name . ".xyz" );
       
       # print coordinates for @t=(1..4) to same filehandle
       foreach my $t ( 1 .. 4 ) {
           $mol->t($t);
           $mol->print_xyz($fh);
       }
       
       $mol->t(0);
       foreach ( 1 .. 10 ) {
           $mol->rotate(
               V( 0, 0, 1 ),    # rotation vector
               36,              # rotate by 180 degrees
               V( 5, 0, 0 )     # origin of rotation
           );
           $mol->print_xyz($fh);
       }
       
       # translate/rotate method is provided by AtomGroupRole
       # populate groups byatom resid attr
       my @groups = $hack->group_by_atom_attr( 'resid', $mol->all_atoms );
       $mol->push_groups(@groups);
       
       foreach my $ang ( 1 .. 10 ) {
           $_->rotate( V( 1, 1, 1 ), 36, $_->COM ) foreach $mol->all_groups;
           $mol->print_xyz($fh);
       }
       
 

DESCRIPTION
============
The HackaMol library enables users to build simple, yet powerful scripts 
for carrying out computational work on molecules at multiple scales. The 
molecular object system organizes atoms within molecules using groups, bonds, 
angles, and dihedrals.  HackaMol seeks to provide intuitive attributes and 
methods that may be harnessed to coerce computational chemistry through a 
common core. The library is inspired by 
*[PerlMol](http://www.perlmol.org)*, *[BioPerl](http://bioperl.org)*, *[MMTSB](http://www.mmtsb.org)*, 
and my own experiences as a researcher.  A goal of this library is to reduce
the "viscosity" of setting up computations and managing data.
       
The library is organized into two regions: HackaMol, the core (contained 
here) that has classes for atoms and molecules, and HackaMolX, the 
extensions, such as HackaMolX::PDB, a parser for protein databank files, 
and HackaMolX::Calculator, an abstract calculator for coercing 
computational chemistry, that use the core. The three major goals of the 
core are for it to be well-tested, well-documented, and easy to install. 
The goal of the extensions is to provide a more flexible space for 
researchers to develop and share new methods that use the core. 
Extensions are in the works, but the HackaMolX namespace has not been 
established yet! 
       
HackaMol uses *[Math::Vector::Real](https://metacpan.org/module/Math::Vector::Real)* (MVR) for all the vector operations. 
MVR is a lightweight solution with a fast XS dropin that overlaps very 
well with the desirables for working with atomic coordinates. Extensions 
that treat much larger systems will definitely benefit from the 
capabilities *[Perl Data Language](http://pdl.perl.org)* (PDL) or *[Math::GSL](https://metacpan.org/module/Math::GSL)*.
       
INSTALLATION
============
To install HackaMol, I recommend using cpanminus and local::lib. This approach avoids the need for root privileges and uses the system Perl 
(available on most systems; see Strawberry Perl for Windows). 
I use *[Perlbrew](http://perlbrew.pl)* to manage local versions of Perl, but Perlbrew is overkill unless your system Perl is very old (if you do use Perlbrew, you won't need local::lib).

Here is a quick summary of a *[step by step post on installing Perl modules] (http://perlmaven.com/install-perl-modules-without-root-rights-on-linux-ubuntu-13-10)*. Below, prompt> denotes your terminal prompt.  

  1. install cpanminus 

       prompt> curl -L http://cpanmin.us | perl - App::cpanminus

    Execution without superuser privileges will complain about not being to write to a path.  This is ok! It will create a perl5 directory in your home directory, in which all modules will be installed. Thus, an uninstall involves deleting the ~/perl5 directory.
    
  2. install local::lib

       prompt> ~/perl5/bin/cpanm --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)  
    
    Execution should install local::lib and configure the cpanm command to know about the ~/perl5 directory in the local shell.  

  3. configure .bash\_profile

       prompt> echo 'eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)' >> ~/.bash\_profile

    Execution will add the stuff between the quotes to your bash_profile, which enables automatic configuration of cpanminus and local::lib with every new terminal shell.
   
  4. install HackaMol

       prompt> cpanm HackaMol
       
    Execution will fail to find HackaMol until it is officially released. For now, you can install the developer release:
       
       prompt> cpanm DEMIAN/HackaMol-0.00_12.tar.gz
       
Feedback and contributions welcome!
