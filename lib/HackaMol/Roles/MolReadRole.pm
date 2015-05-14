package HackaMol::Roles::MolReadRole;

# ABSTRACT: Read files with molecular information
use Moose::Role;
use Carp;
use Math::Vector::Real;
use HackaMol::PeriodicTable qw(%KNOWN_NAMES);
use FileHandle;
use HackaMol::Atom;    #add the code,the Role may better to understand
use List::MoreUtils qw(singleton);

with qw(
        HackaMol::Roles::ReadZmatRole
        HackaMol::Roles::ReadPdbRole
        HackaMol::Roles::ReadPdbqtRole
        HackaMol::Roles::ReadXyzRole
);

has 'hush_read' => (
    is      => 'rw',
    isa     => 'Bool',
    lazy    => 1,
    default => 0,
);

sub read_file_atoms {
    my $self = shift;
    my $file = shift;
    my @atoms;

    if ( $file =~ m/\.pdb$/ ) {
        @atoms = $self->read_pdb_atoms($file);
    }
    elsif ( $file =~ m/\.pdbqt$/ ) {
        @atoms = $self->read_pdbqt_atoms($file);
    }
    elsif ( $file =~ m/\.xyz$/ ) {
        @atoms = $self->read_xyz_atoms($file);
    }
    elsif ( $file =~ m/\.zmat$/ ) {
        @atoms = $self->read_zmat_atoms($file);
    }
    else {
        croak "$file format not supported";
    }
    return (@atoms);
}

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

   use HackaMol;

   my $hack   = HackaMol->new( name => "hackitup" );

   # build array of carbon atoms from pdb [xyz,pdbqt] file
   my @carbons  = grep {
                        $_->symbol eq "C"
                       } $hack->read_file_atoms("t/lib/1L2Y.pdb"); 

   my $Cmol     = HackaMol::Molecule->new(
                        name => "carbonprotein", 
                        atoms => [ @carbons ]
                  );

   $Cmol->print_pdb;   
   $Cmol->print_xyz;     

   # build molecule from xyz [pdb,pdbqt] file
   my $mol    = $hack->read_file_mol("some.xyz");
   $mol->print_pdb; # not so easy from xyz to pdb! 

=head1 DESCRIPTION

The HackaMol::MolReadRole role provided methods for reading common structural files.  Currently,
pdb and xyz are provided in the core, but others will be likely added.  

=attr hush_read

isa Bool that is lazy. $hack->hush_read(1) will quiet some warnings that may be ignored under some instances.

=method read_file_atoms

takes the name of the file as input, parses the file, builds Atom objects, and returns them.
Matches the filename extension and calls on either read_pdb_atoms or read_xyz_atoms

=method read_pdb_atoms

takes the name of the file as input, parses the pdb file to return the list of built 
Atom objects. This is a barebones parser.  A more advanced PDB parser will be released 
soon as an extension. 

According to the PDB specification, the element symbol should be present in columns 77-78.  
The element is often ommitted by programs, such as charmm, that can write pdbs because it makes the
file larger, and the information is accessible somewhere else. Unfortunately, other programs require
the information.  HackaMol::MolReadRole, loads a hash (KNOWN_NAMES) from HackaMol::PeriodicTable 
that maps common names to the element (e.g. POT => 'K'). read_pdb_atoms will carp if the name is 
not in the hash, and then set the element to the first letter of the name. This will be improved when
HackaMol::PeriodicTable is improved. See TODO.

=method read_xyz_atoms

takes the name of the file as input, parses the xyz file to return the list of built 
Atom objects.  

=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<Protein Data Bank | http://pdb.org>

