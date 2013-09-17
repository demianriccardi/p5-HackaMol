HackaMol
========
Object-Oriented Perl 5, Moose Library for Molecular Hacking on multiple scales

VERSION
========
       version 0.00_01

SYNOPSIS
========
       use HackaMol;

       my $hack  = HackaMol->new(name => "hackitup");

       @atoms = $hack->read_file_atoms("t/lib/1L2Y.pdb");

       my $mol = Molecule->new(name=> 'trp-cage', atoms=>[@atoms]);

       $mol->translate(-$mol->COM);

       $mol->rotate(V(1,0,0), 180, V(10,10,10));

       print $mol->count_atoms;

       print "\n";

       printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach
       $mol->all_atoms;

       my @groups = $hack->group_by_atom_attr('resid'); #group by atom resid attr 
       $mol->push_groups(@groups);

       $_->rotate(V(1,1,1),60,$_->COM,1) foreach $mol->all_groups; # mess up
       all the amino acids

       print $mol->count_atoms . "\n";

       print "\n";

       printf("%5s %8.3f %8.3f %8.3f\n", $_->Z, @{$_->xyz}) foreach $mol->all_atoms;

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
