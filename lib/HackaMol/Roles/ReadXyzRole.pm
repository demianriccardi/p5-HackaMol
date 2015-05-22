package HackaMol::Roles::ReadXyzRole;

# ABSTRACT: Read files with molecular information
use Moo::Role;
use strictures 2; 
use Carp;
use Math::Vector::Real;
use FileHandle;

sub read_xyz_atoms {

    #read xyz file and generate list of Atom objects
    my $self = shift;
    my $file = shift;
    my $fh   = FileHandle->new("<$file") or croak "unable to open $file";

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
                    $atoms[$n] = HackaMol::Atom->new(
                        name   => "at$n",
                        Z      => $sym,
                        coords => [$xyz]
                    );
                }
                else {
                    $atoms[$n] = HackaMol::Atom->new(
                        name   => "at$n",
                        symbol => $sym,
                        coords => [$xyz]
                    );
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

1;

__END__

=head1 SYNOPSIS

   my @atoms = HackaMol->new
                       ->read_xyz_atoms("some.xyz");

=head1 DESCRIPTION

The HackaMol::Roles::ReadXyzRole provides read_xyz_atoms reading xyz files.

=method read_xyz_atoms

One argument: the filename
Returns a list of HackaMol::Atom objects.

=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<HackaMol::Atom>
* L<HackaMol::Roles::MolReadRole>

