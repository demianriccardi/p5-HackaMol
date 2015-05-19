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
   $mol->print_pdb; # 

=head1 DESCRIPTION

The HackaMol::Role::MolReadRole role provides methods for reading common structural files.  Currently,
pdb, pdbqt, Z-matrix, and xyz are provided.  The methods are all provided in separate roles.  Adding 
additional formats is straightforward: 

    1. Add a Role that parses the file and returns a list of HackaMol::Atoms. 

    2. Add the code here to consume the role and call the method based on the file ending.

=attr hush_read

isa Bool that is lazy. $hack->hush_read(1) will quiet some warnings that may be ignored under some instances.

=method read_file_atoms

one argument: the name of a file (.xyz, .pdb, .pdbqt, .zmat)

returns a list of HackaMol::Atom objects

=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<HackaMol::Roles::ReadXyzRole>
* L<HackaMol::Roles::ReadPdbRole>
* L<HackaMol::Roles::ReadPdbqtRole>
* L<HackaMol::Roles::ReadZmatRole>
* L<HackaMol::Atom>
* L<HackaMol::Molecule>

