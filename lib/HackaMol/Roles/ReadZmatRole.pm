package HackaMol::Roles::ReadZmatRole;

# ABSTRACT: Read files with molecular information
use Moo::Role;
use strictures 2;
use HackaMol::PeriodicTable qw(%KNOWN_NAMES);
use Math::Vector::Real;
use Carp;
use List::MoreUtils qw(singleton);

with qw(
        HackaMol::Roles::NERFRole
);

sub read_zmat_atoms {

    #xyz file and generate list of Atom object
    my $self = shift;
    my $file = shift;
    my $fh   = FileHandle->new("<$file") or croak "unable to open $file";

    my @atoms;
    my ( $n, $t ) = ( 0, 0 );

    my $nat  = undef;
    my @zmat = <$fh>;
    chomp @zmat;
    #use Data::Dumper; 
    #print Dumper \@zmat; exit;

    # we have 5 types of extensions
    # A. SYM 0 x y z
    # B. SYM
    # C. SYM i R
    # D. SYM i R j Ang
    # E. SYM i R j Ang k Tors
    # we need to filter the indices (can't lose the location)

    #type A
    my @iA = grep { $zmat[$_] =~ m/^\s*\w+\s+0(\s+\d+\.*\d*){3}/ } 0 .. $#zmat;
    my @inA = singleton( 0 .. $#zmat, @iA );
    #type B
    my @iB = grep { $zmat[$_] =~ m/^\s*\w+\s*$/ } @inA;
    #type C
    my @iC = grep { $zmat[$_] =~ m/^\s*\w+(\s+\d+\s+\d+\.*\d*)\s*$/ } @inA;
    #type D
    my @iD = grep { $zmat[$_] =~ m/^\s*\w+(\s+\d+\s+\d+\.*\d*){2}\s*$/ } @inA;
    #type E
    my @iE = grep {
        $zmat[$_] =~ m/^\s*\w+(\s+\d+\s+\d+.*\d*){2}\s+\d+\s+-*\d+.*\d*\s*$/
    } @inA;

    foreach my $ia (@iA) {
        my ( $sym, $iat1, @xyz ) = split( / /, $zmat[$ia] );
        $atoms[$ia] = HackaMol::Atom->new(
                        symbol => $sym,
                        coords => [ V(@xyz) ]
        );
    }

    #print Dump 'A', \%mol;

    foreach my $ib (@iB) {
        my $sym = $zmat[$ib];
        my $a   = $self->init;
        $atoms[$ib] = HackaMol::Atom->new(
            symbol => $sym,
            coords => [$a]
        );
    }

    #print Dump 'B', \%mol;

    foreach my $ic (@iC) {
        my ( $sym, $iat1, $R ) = split( / /, $zmat[$ic] );
        my $a = $atoms[ $iat1 - 1 ]->xyz;
        my $b = $self->extend_a( $a, $R );
        $atoms[$ic] = HackaMol::Atom->new(
            symbol => $sym,
            coords => [$b]
        );
    }

    #print Dump 'C', \%mol;

    foreach my $id (@iD) {
        my ( $sym, $iat1, $R, $iat2, $ang ) = split( / /, $zmat[$id] );
        my $a = $atoms[ $iat1 - 1 ]->xyz;
        my $b = $atoms[ $iat2 - 1 ]->xyz;
        my $c = $self->extend_ab( $b, $a, $R, $ang );
        $atoms[$id] = HackaMol::Atom->new(
            symbol => $sym,
            coords => [$c]
        );
    }

    #print Dump 'D', \%mol;

    foreach my $ie (@iE) {
        my ( $sym, $iat1, $R, $iat2, $ang, $iat3, $tor ) =
          split( / /, $zmat[$ie] );
        my $a = $atoms[ $iat1 - 1 ]->xyz;
        my $b = $atoms[ $iat2 - 1 ]->xyz;
        my $c = $atoms[ $iat3 - 1 ]->xyz;
        my $d = $self->extend_abc( $c, $b, $a, $R, $ang, $tor );
        $atoms[$ie] = HackaMol::Atom->new(
            symbol => $sym,
            coords => [$d]
        );
    }
    return (@atoms);

}

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

