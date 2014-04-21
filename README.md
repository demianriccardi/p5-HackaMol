HackaMol
========
Object-oriented Perl 5, Moose library for molecular hacking on multiple scales

VERSION
========
developer version 0.00_15 
Available for testing from cpan.org:
       
please see *[HackaMol on MetaCPAN](https://metacpan.org/release/DEMIAN/HackaMol-0.00_15) for formatted documentation.  
       
SYNOPSIS
========
       use HackaMol;
       use Math::Vector::Real;
       
       my $hack = HackaMol->new( name => "hackitup" );

       # all coordinates from NMR ensemble are loaded into atoms
       my $mol  = $hack->read_file_mol("1L2Y.pdb");
       
       #recenter all coordinates to center of mass
       foreach my $t ( 0 .. $mol->tmax) {
           $mol->t($t);
           $mol->translate( -$mol->COM );
       }
       
       #create array of CA atoms with full occupancy 

       my @CAs = grep {
                        $_->name    eq 'CA'  and
                        $_->occ == 1 
                      } $mol->all_atoms;
      
       #print out the pdb with CA for several models from the NMR 
       HackaMol::Molecule->new( 
                                atoms=>[@CAs] 
                              )-> print_pdb_ts([8,2,4,6,8,0], 'some.pdb');

       # print coordinates to trp-cage.xyz and return filehandle for future
       # writing
       my $fh = $mol->print_xyz( $mol->name . ".xyz" );
       foreach ( 1 .. 10 ) {
           $mol->rotate(
               V( 0, 0, 1 ),    # rotation vector
               36,              # rotate by 36 degrees
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

    Execution will add the stuff between the quotes to your bash_profile, which enables automatic configuration of cpanminus and local::lib with every new terminal shell. To uninstall, delete this line from .bash\_profile
   
    You are now ready to install anything installable from CPAN!!!

  4. install HackaMol

       prompt> cpanm HackaMol
       
    Execution will fail to find HackaMol until it is officially released. For now, you can install the developer release:
       
       prompt> cpanm DEMIAN/HackaMol-0.00_15.tar.gz
       
HackaMol is too fun to not be experimental! Feedback and contributions welcome!
