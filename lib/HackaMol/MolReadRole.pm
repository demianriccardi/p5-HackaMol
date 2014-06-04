package HackaMol::MolReadRole;

# ABSTRACT: Read files with molecular information
use Moose::Role;
use Carp;
use Math::Vector::Real;
use HackaMol::PeriodicTable qw(%KNOWN_NAMES);
use FileHandle;

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
    else {
        croak "$file format not supported";
    }
    return (@atoms);
}



sub read_pdb_atoms {

    #read pdb file and generate list of Atom objects
    my $self  = shift;
    my $file  = shift;
    my $fh = FileHandle->new("<$file") or croak "unable to open $file";

    my @atoms;
    my ( $n, $t ) = ( 0, 0 );
    my $q_tbad = 0;
    my $something_dirty = 0;

    while (<$fh>) {

        if (/^(?:MODEL\s+(\d+))/) {
            #$t = $1 - 1; # I don't like this!!  set increment t instead.. below
            $n = 0;
            $q_tbad = 0; # flag a bad model and never read again!
        }
        elsif (/^(?:ENDMDL)/){
            $t++;
        }
        elsif (/^(?:HETATM|ATOM)/) {
            next if $q_tbad;
            my (
                $record_name, $serial, $name, $altloc,  $resName,
                $chainID,     $resSeq, $icod, $x,       $y,
                $z,           $occ,    $B,    $segID,   $element, $charge
            ) = unpack "A6A5x1A4A1A3x1A1A4A1x3A8A8A8A6A6x6A4A2A2", $_;

            if   ( $charge =~ m/\d/ ) { $charge = _qstring_num($charge) }
            else                      { $charge = 0 }

            if   ( $chainID =~ m/\w/ ) { $chainID = uc( _trim($chainID) ) }
            else                       { $chainID = ' ' }

            
            $name    = _trim($name);
            $resName = _trim($resName);
            $resSeq  = _trim($resSeq);
            #$resSeq  = 0 if ( $resSeq < 0 );
            $serial  = _trim($serial);
            $segID   = _trim($segID);
  
            $element = ucfirst( lc( _trim($element) ) );
            my $qdirt = 0;
            ($element,$qdirt) = _element_name($name) unless ($element =~ /\w+/);
            $something_dirty++ if ($qdirt); 
            my $xyz = V( $x, $y, $z );

            if ( $t == 0 ) {
                $atoms[$n] = HackaMol::Atom->new(
                    name        => $name,
                    record_name => $record_name,
                    serial      => $serial,
                    chain       => $chainID,
                    symbol      => $element,
                    charges     => [$charge],
                    coords      => [$xyz],
                    occ         => $occ * 1,
                    bfact       => $B * 1,
                    resname     => $resName,
                    resid       => $resSeq,
                    segid       => $segID ,
                    altloc      => $altloc,
                );
                $atoms[$n]->is_dirty($qdirt) unless $atoms[$n]->is_dirty; 
            }
            else {
                #croak condition if atom changes between models
                if ( $name ne $atoms[$n]->name
                    or $element ne $atoms[$n]->symbol ) {
                  my $carp_message = "BAD t->$t PDB Atom $n "
                                    ."serial $serial resname $resName "
                                    ."has changed";
                  carp $carp_message;
                  $q_tbad = $t; # this is a bad model!
                  #wipe out all the coords prior
                  $atoms[$_]->delete_coords($t) foreach 0 .. $n-1; 
                  $t--;
                  next;
                }
                $atoms[$n]->set_coords( $t, $xyz );
            }
            $n++;
        }
    }

    # set iatom to track the array.  diff from serial which refers to pdb
    $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
    if ($something_dirty){
      unless($self->hush_read){
        my $message = "MolReadRole> found $something_dirty dirty atoms. ";
        $message .= "Check symbols and lookup names";
        carp $message;
      }
    }
    return (@atoms);
}

sub read_pdbqt_atoms {
    # this is too similar to reading pdb for it too exist separately... think about
    my $self  = shift;
    my $file  = shift;
    
    my $fh = FileHandle->new("<$file") or croak "unable to open $file";

    # $RtBrnch{ROOT}{iatoms}    = LIST integers of atoms in root 
    # $RtBrnch{BRNCH1}{iatoms}  = LIST integers of atoms in BRNCH1 
    # $RtBrnch{BRNCH1}{ROOT}    = LIST two integers of rotable bond for BRNCH1
    # $RtBrnch{BRNCH2}{iatoms}  = LIST integers of atoms in BRNCH2
    # $RtBrnch{BRNCH2}{SBRNCH1} iatoms}  = LIST integers of atoms in BRNCH2
    # each branch is rigid block of atoms
    #  NONONONONO!!!!  for now, just read in atoms.  We'll figure out how to save the tree later 
    my %RtBrnch; 

    my @atoms;
    my ( $n, $t ) = ( 0, 0 );
    my $q_tbad = 0;
    my $something_dirty = 0;

    while (<$fh>) {

        if (/^(?:MODEL\s+(\d+))/) {
            $n = 0;
            $q_tbad = 0; # flag a bad model and never read again!
        }
        elsif (/^(?:ENDMDL)/){
            $t++;
        }
        elsif (/^(?:HETATM|ATOM)/) {
            next if $q_tbad;
            my (
                $record_name, $serial, $name, $altloc,  $resName,
                $chainID,     $resSeq, $icod, $x,       $y,
                $z,           $occ,    $B,    $charge,  $ADTtype
            ) = unpack "A6A5x1A4A1A3x1A1A4A1x3A8A8A8A6A6x4A6x1A2", $_;
#ATOM      1  O   LIG d   1       8.299   4.799  79.371  0.00  0.00    -0.292 OA
#-----|----|x---||--|x|---||xxx-------|-------|-------|-----|-----|xxxx-----|x-|

            if   ( $chainID =~ m/\w/ ) { $chainID = uc( _trim($chainID) ) }
            else                       { $chainID = ' ' }

            
            $name    = _trim($name);
            $resName = _trim($resName);
            $resSeq  = _trim($resSeq);
            #$resSeq  = 0 if ( $resSeq < 0 );
            $serial  = _trim($serial);
            $charge  = _trim($charge);
            $ADTtype = _trim($ADTtype);
  
            my ($element,$qdirt) = _element_name($ADTtype);
            $element = 'C' if ($element eq 'A'); #aromatic, is this dirty?
            $something_dirty++ if ($qdirt); 
            my $xyz = V( $x, $y, $z );

            if ( $t == 0 ) {
                $atoms[$n] = HackaMol::Atom->new(
                    name        => $name,
                    record_name => $record_name,
                    serial      => $serial,
                    chain       => $chainID,
                    symbol      => $element,
                    charges     => [$charge],
                    coords      => [$xyz],
                    occ         => $occ * 1,
                    bfact       => $B * 1,
                    resname     => $resName,
                    resid       => $resSeq,
                    segid       => $ADTtype ,
                    altloc      => $altloc,
                );
                $atoms[$n]->is_dirty($qdirt) unless $atoms[$n]->is_dirty; 
            }
            else {
                #croak condition if atom changes between models
                if ( $name ne $atoms[$n]->name
                    or $element ne $atoms[$n]->symbol ) {
                  my $carp_message = "BAD t->$t PDB Atom $n "
                                    ."serial $serial resname $resName "
                                    ."has changed";
                  carp $carp_message;
                  $q_tbad = $t; # this is a bad model!
                  #wipe out all the coords prior
                  $atoms[$_]->delete_coords($t) foreach 0 .. $n-1; 
                  $t--;
                  next;
                }
                $atoms[$n]->set_coords( $t, $xyz );
            }
            $n++;
        }
    }

    # set iatom to track the array.  diff from serial which refers to pdb
    $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
    if ($something_dirty){
      unless($self->hush_read){
        my $message = "MolReadRole> found $something_dirty dirty atoms. ";
        $message .= "Check symbols and lookup names";
        carp $message; 
      }
    }  
    return (@atoms);
}

sub read_xyz_atoms {

    #read pdb file and generate list of Atom objects
    my $self  = shift;
    my $file  = shift;
    my $fh = FileHandle->new("<$file") or croak "unable to open $file";

    my @atoms;
    my ( $n, $t ) = ( 0, 0 );

    my $nat = undef;
    while (<$fh>) {

        if (/^(\s*\d+\s*)$/) {
            $n = (split)[0];
            if ( defined($nat) ) {
                croak "number of atoms has changed\n" unless ( $nat == $n );
                $t++;
            }
            $nat = $n;
            $n   = 0;
        }
        elsif (/(\w+|\d+)(\s+-*\d+\.\d+){3}/) {
            my @stuff = split;
            my $sym   = $stuff[0];
            my $xyz   = V( @stuff[ 1, 2, 3 ] );
            if ( $t == 0 ) {
                if ( $sym =~ /\d/ ) {
                    $atoms[$n] =
                      HackaMol::Atom->new(name=> "at$n", Z => $sym, coords => [$xyz] );
                }
                else {
                    $atoms[$n] =
                      HackaMol::Atom->new(name=> "at$n", symbol => $sym, coords => [$xyz] );
                }
            }
            else {
                if ( $sym =~ /\d/ ) {
                    croak "atoms have changed from last model to current: $t\n"
                      if ( $sym != $atoms[$n]->Z );
                }
                else {
                    croak "atoms have changed from last model to current: $t\n"
                      if ( $sym ne $atoms[$n]->symbol );
                }
                $atoms[$n]->set_coords( $t, $xyz );

            }
            $n++;
        }
    }

    # set iatom to track the array.  diff from serial which refers to pdb
    $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
    return (@atoms);
}

sub _trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}


sub _qstring_num {

    # _qstring something like 2+  or 2-
    my $string = shift;
    $string =~ s/\+//;
    $string =~ s/(.*?)(\-)/$2$1/;
    $string = sprintf( "%g", $string );
    return $string;

}

sub _element_name{
# guess the element using the atom name
  my $name  = uc(shift);
  my $dirt  = 0;
  unless (exists($KNOWN_NAMES{$name})){
    #carp "$name doesn not exist in HackaMol::PeriodicTable, if common please add to KNOWN_NAMES";
    $dirt = 1;
    my $symbol = substr $name, 0,1;  #doesn't work if two letters for symbol
    return ($symbol,$dirt);
  }
  return ($KNOWN_NAMES{$name},$dirt);
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

