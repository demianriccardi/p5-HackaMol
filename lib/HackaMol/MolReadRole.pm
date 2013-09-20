package HackaMol::MolReadRole;
# ABSTRACT: Read XYZ and PDB files 
use Moose::Role;
use Carp;
use Math::Vector::Real;
use FileHandle;

sub read_file_atoms {
  my $self = shift;
  my $file = shift;
  my @atoms;

  if ($file =~ m/\.pdb$/){
    @atoms = $self->read_pdb_atoms($file);
  }
  elsif ($file =~ m/\.xyz$/){
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
  my $segid = $file;
  $segid =~ s/\.pdb//;
  $segid =~ s/t\/lib\///;
  my $fh = FileHandle->new("<$file") or croak "unable to open $file";;

  my @atoms;
  my ( $n, $t ) = ( 0, 0 );

  while (<$fh>) {

      if (/^(?:MODEL\s+(\d+))/) {
          $t = $1 - 1;
          $n = 0;
      }
      elsif (/^(?:HETATM|ATOM)/) {

          my (
              $record_name, $serial, $name, $altloc,  $resName,
              $chainID,     $resSeq, $icod, $x,       $y,
              $z,           $occ,    $B,    $element, $charge
          ) = unpack "A6A5x1A4A1A3x1A1A4A1x3A8A8A8A6A6x10A2A2", $_;

          if   ( $charge =~ m/\d/ ) { $charge = _qstring_num($charge) }
          else                      { $charge = 0 }

          if   ( $chainID =~ m/\w/ ) { $chainID = uc( _trim($chainID) ) }
          else                       { $chainID = 'AA' }

          $element = ucfirst( lc( _trim($element) ) );
          $name    = _trim($name);
          $resName = _trim($resName);
          $resSeq  = _trim($resSeq);
          $resSeq  = 0 if ( $resSeq < 0 );
          $serial  = _trim($serial);
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
                  altloc      => $altloc,
              );
          }
          else {
              croak "atoms have changed from last model to current: $t\n"
                if ( 
                    $name    ne $atoms[$n]->name   or 
                    $element ne $atoms[$n]->symbol
                    );

              $atoms[$n]->set_coords( $t, $xyz );
          }
          $n++;
      }
  }

  # set iatom to track the array.  diff from serial which refers to pdb
  $atoms[$_]->iatom($_) foreach ( 0 .. $#atoms );
  return (@atoms);
}

sub read_xyz_atoms {
#read pdb file and generate list of Atom objects
  my $self  = shift;
  my $file  = shift;
  my $segid = $file;
  $segid =~ s/\.xyz//;
  $segid =~ s/t\/lib\///;
  my $fh = FileHandle->new("<$file") or croak "unable to open $file";;

  my @atoms;
  my ( $n, $t ) = ( 0, 0 );

  my $nat  = undef;
  while (<$fh>) {

    if (/^(\s*\d+\s*)$/){
      $n = (split)[0];
      if (defined($nat)){
        croak "number of atoms has changed\n" unless ($nat == $n);
        $t++;
      }
      $nat  = $n;
      $n = 0;
    }
    elsif (/(\w+|\d+)(\s+-*\d+\.\d+){3}/) {
      my @stuff = split;
      my $sym = $stuff[0];
      my $xyz = V( @stuff[1,2,3] );
      if ( $t == 0 ) {
        if ($sym =~ /\d/)
             {
         $atoms[$n] = HackaMol::Atom->new( Z      => $sym,coords => [$xyz])
        }
        else {
         $atoms[$n] = HackaMol::Atom->new( symbol => $sym,coords => [$xyz])
        }
      }
      else {
        if ($sym =~ /\d/){
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

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

use HackaMol;

my $hack   = HackaMol->new( name => "hackitup" );
my @atoms1 = $hack->read_file_atoms("t/lib/1L2Y.pdb"); 
my @atoms2 = $hack->read_file_atoms("t/lib/something.xyz"); 

=head1 DESCRIPTION

The HackaMol::MolReadRole role provided methods for reading common structural files.  Currently,
pdb and xyz are provided in the core, but others will be likely added.  

=method read_file_atoms

takes the name of the file as input, parses the file, builds Atom objects, and returns them.
Matches the filename extension and calls on either read_pdb_atoms or read_xyz_atoms

=method read_pdb_atoms

takes the name of the file as input, parses the pdb file to return the list of built 
Atom objects. This is a barebones parser.  A more advanced PDB parser will be released 
soon as an extension.

=method read_xyz_atoms

takes the name of the file as input, parses the xyz file to return the list of built 
Atom objects.  

=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<Protein Data Bank | http://pdb.org>

