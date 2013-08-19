package PDBintoAtoms;
use Carp;
use Math::Vector::Real;
require Exporter;
our @ISA=qw(Exporter);
our @EXPORT_OK = qw(readinto_atoms);
use lib 'lib/HackaMol';
use Atom;
use FileHandle;

sub readinto_atoms{
    my $file = shift;
    my $segid = $file;
    $segid =~ s/\.pdb//;
    $segid =~ s/t\/lib\///;
    my $fh = FileHandle->new("<$file");

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

            if ($charge =~ m/\d/){$charge  = _qstring_num($charge)}
                            else {$charge  = 0}

            if ($chainID =~ m/\w/){$chainID = uc(_trim($chainID) ) }
                              else{$chainID = 'AA'}

            $element = ucfirst( lc( _trim($element) ) );
            $name = _trim($name);
            $resName = _trim($resName);
            $resSeq = _trim($resSeq);
            $resSeq = 0 if ($resSeq<0);
            $serial = _trim($serial);
  
            my $xyz = V($x, $y, $z);

            if ( $t == 0 ) {
                $atoms[$n] = Atom->new(
                    name        => $name,
                    record_name => $record_name,
                    serial      => $serial,
                    chain       => $chainID,
                    symbol      => $element,
                    charges     => [$charge],
                    coords      => [$xyz],
                    occ         => $occ*1,
                    bfact       => $B*1,
                    resname     => $resName,
                    resid       => $resSeq,
                    altloc      => $altloc,
                );
            }
            else {
                croak "atoms have changed from last model to current: $t\n"
                  if ( $name ne $atoms[$n]->name );
 
                $atoms[$n]->set_coords( $t, $xyz );
            }
            $n++;
        }
    }
    # set iatom to track the array.  diff from serial which refers to pdb
    $atoms[$_]->iatom($_) foreach (0 .. $#atoms); 
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

1;

