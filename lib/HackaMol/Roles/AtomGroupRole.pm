package HackaMol::Roles::AtomGroupRole;

#ABSTRACT: Role for a group of atoms
use Moose::Role;
use Carp;
use Math::Trig;
use Math::Vector::Real;
use FileHandle;
use Scalar::Util 'reftype';
use List::Util qw(sum);

#use MooseX::Storage;
#with Storage( 'format' => 'JSON', 'io' => 'File', traits => ['OnlyWhenBuilt'] );

my $angste_debye = 4.80320;

has 'atoms' => (
    traits  => ['Array'],
    is      => 'ro',
    isa     => 'ArrayRef[HackaMol::Atom]',
    default => sub { [] },
    handles => {
        unshift_atoms => 'unshift',
        push_atoms    => 'push',
        select_atoms  => 'grep',
        map_atoms     => 'map',
        sort_atoms    => 'sort',
        get_atoms     => 'get',
        set_atoms     => 'set',
        insert_atoms  => 'insert',
        delete_atoms  => 'delete',
        all_atoms     => 'elements',
        count_atoms   => 'count',
        natoms        => 'count',
        clear_atoms   => 'clear',
        has_atoms     => 'count',
    },
    lazy => 1,
);

has 'is_constrained' => (
    is      => 'rw',
    isa     => 'Bool',
    lazy    => 1,
    default => 0,
);

has 'qcat_print' => (
    is      => 'rw',
    isa     => 'Bool',
    lazy    => 1,
    default => 0,
);

has 'info' => (
    is      => 'rw',
    isa     => 'Str',
    lazy    => 1,
    default => "",
);

sub _lsq_slope {
    # private function used for each coordinate 
    # translation of tcl function written by Justin Gullingsrud @ uiuc.edu 
    # algorithm reference: Bevington
    # Fit the points x to x = ai + b, i=0...N-1, and return the value of a 
    # a = 12/( (N(N^2 - 1)) ) sum[ (i-(N-1)/2) * xi]    
    my $xis = shift || die "expecting array_ref of cartesian coordinate [x y or z]";
    my $n = @$xis;
    my $sum = 0;
    my $d = ($n-1)/2;
    my $i = 0;
    $sum += ($i++ - $d)*$_ foreach @$xis; 
    return $sum * 12 / ($n * ( $n * $n - 1)) ;
} 

sub centered_vector {
    my $self = shift;
    my $centered_flag = shift;
    my @mvrs = map {$_->xyz} $self->all_atoms;
    die "2 atoms needed for a centered_vector" unless @mvrs > 1;
    my @x = map { $_->[0] } @mvrs;
    my @y = map { $_->[1] } @mvrs;
    my @z = map { $_->[2] } @mvrs;
    my $mvr = V( map { _lsq_slope($_) } \@x,\@y,\@z);
    return $mvr->versor;
}

sub calc_bfps {

    # this should be rerun for each selection
    #10.1016/j.jmb.2015.09.024
    my $self = shift;
    unless ( $self->count_atoms > 1 ) {
        warn "calc_bfps> group not large enough\n";
        return;
    }
    my @atoms      = $self->all_atoms;
    my @bfacts     = map { $_->bfact } @atoms;
    my $bfact_mean = sum(@bfacts) / @bfacts;
    my $sd         = 0;
    $sd += ( $_ - $bfact_mean )**2 foreach @bfacts;
    unless ( $sd > 0 ) {
        warn "calc_bfps> no variance in the group bfactors\n";
        return;
    }
    my $bfact_std = sqrt( $sd / ( @bfacts - 1 ) );
    foreach my $atom (@atoms) {
        my $bfp = ( $atom->bfact - $bfact_mean ) / $bfact_std;
        $atom->bfp($bfp);
    }
    return map { $_->bfp } @atoms;
}

sub dipole {
    my $self = shift;
    return ( V(0) ) unless ( $self->count_atoms );
    my @atoms   = $self->all_atoms;
    my @vectors = grep { defined } map { $_->get_coords( $_->t ) } @atoms;
    my @charges = grep { defined } map { $_->get_charges( $_->t ) } @atoms;
    my $dipole  = V( 0, 0, 0 );
    if ( $#vectors != $#charges ) {
        carp
          "build_dipole> mismatch number of coords and charges. all defined?";
        return $dipole;
    }
    $dipole += $vectors[$_] * $charges[$_] foreach 0 .. $#charges;
    return ($dipole);
}

sub COM {
    my $self = shift;
    return ( V(0) ) unless ( $self->count_atoms );
    my @atoms     = $self->all_atoms;
    my @m_vectors = map { $_->mass * $_->get_coords( $_->t ) } @atoms;
    my $com       = V( 0, 0, 0 );
    $com += $_ foreach @m_vectors;
    return ( $com / $self->total_mass );
}

sub center {
    my $self = shift;
    return ( V(0) ) unless ( $self->count_atoms );
    my @atoms   = $self->all_atoms;
    my @vectors = map { $_->get_coords( $_->t ) } @atoms;
    my $center  = V( 0, 0, 0 );
    $center += $_ foreach @vectors;
    return ( $center / $self->count_atoms );
}

sub COZ {
    my $self = shift;
    return ( V(0) ) unless ( $self->count_atoms );
    my @atoms     = $self->all_atoms;
    my @z_vectors = map { $_->Z * $_->get_coords( $_->t ) } @atoms;
    my $coz       = V( 0, 0, 0 );
    $coz += $_ foreach @z_vectors;
    return ( $coz / $self->total_Z );
}

sub gt {

    #set group time
    my $self = shift;
    my $t    = shift;
    $self->do_forall( 't', $t );
}

sub do_forall {
    my $self   = shift;
    my $method = shift;
    do { carp "doing nothing for all"; return } unless (@_);
    my @atoms = $self->all_atoms;
    $_->$method(@_) foreach @atoms;
}

sub total_charge {
    my $self = shift;
    return (0) unless ( $self->count_atoms );
    my @atoms   = $self->all_atoms;
    my @charges = map { $_->get_charges( $_->t ) } @atoms;
    my $sum     = 0;
    $sum += $_ foreach @charges;
    return ($sum);
}

sub total_mass {
    my $self = shift;
    return (0) unless ( $self->count_atoms );
    my @masses = map { $_->mass } $self->all_atoms;
    my $sum = 0;
    $sum += $_ foreach @masses;
    return ($sum);
}

sub total_Z {
    my $self = shift;
    return (0) unless ( $self->count_atoms );
    my @Zs = map { $_->Z } $self->all_atoms;
    my $sum = 0;
    $sum += $_ foreach @Zs;
    return ($sum);
}

sub dipole_moment {
    my $self = shift;
    return ( abs( $self->dipole ) * $angste_debye );
}

sub bin_atoms {

    # Called with no arguments.
    # Returns a hash with a count of unique atom symbols
    my $self   = shift;
    my $bin_hr = $self->bin_this('symbol');
    return ($bin_hr);
}

sub count_unique_atoms {
    my $self   = shift;
    my $bin_hr = $self->bin_atoms;
    return ( scalar( keys %{$bin_hr} ) );
}

sub bin_atoms_name {

    # return something like C4H10 sort in order of descending Z
    my $self   = shift;
    my $bin_hr = $self->bin_atoms;
    my $z_hr;
    $z_hr->{ $_->symbol } = $_->Z foreach $self->all_atoms;

    my @names = map {
        my $name = $_ . $bin_hr->{$_};
        $name =~ s/(\w+)1$/$1/;
        $name;    # substitue 1 away?
      }
      sort {
        $z_hr->{$b} <=> $z_hr->{$a}    # sort by Z!  see above...
      } keys %{$bin_hr};
    return join( '', @names );
}

sub translate {
    my $self = shift;
    my $tvec = shift or croak "pass MVR translation vector";
    my $tf   = shift;

    my @atoms = $self->all_atoms;
    do { carp "no atoms to translate"; return } unless (@atoms);
    $tf = $atoms[0]->t unless ( defined($tf) );

    foreach my $at (@atoms) {
        my $v = $at->xyz + $tvec;
        $at->set_coords( $tf, $v );
    }
}

sub rotate {

    #rotate about origin. having origin allows rotation of subgroup
    #without having to translate everything.
    my $self = shift;
    my $rvec = shift or croak "pass MVR rotation vector";
    my $ang  = shift or croak "pass rotation angle";
    my $orig = shift or croak "pass MVR origin";
    my $tf   = shift;

    my @atoms = $self->all_atoms;
    my $t     = $atoms[0]->t;
    $tf = $t unless ( defined($tf) );
    $rvec = $rvec->versor;    #unit vector

    my @cor = map { $_->get_coords($t) - $orig } @atoms;
    my @rcor = $rvec->rotate_3d( deg2rad($ang), @cor );

    $atoms[$_]->set_coords( $tf, $rcor[$_] + $orig ) foreach 0 .. $#rcor;
}

sub rotate_translate {

# args:
#      1. rotation matrix (3x3): ar, each column (cx,cy,cz, below) is a Math::Vector::Real
#      2. translate vector, MVR
#               r' = x*cx + y*cy + z*cz + v_trans
#      3. t final, the final t location for transformed coordinates [defaults to current t]
    my $self  = shift;
    my $rmat  = shift;
    my $trns  = shift;
    my $tf    = shift;
    my @atoms = $self->all_atoms;
    my ( $cx, $cy, $cz ) = @{$rmat};

    my $t = $atoms[0]->t;
    $tf = $t unless ( defined($tf) );

    foreach my $atom (@atoms) {
        my $xyz = $atom->xyz;

        #my ($x,$y,$z) = @{$xyz};
        my $xr = $cx * $xyz;
        my $yr = $cy * $xyz;
        my $zr = $cz * $xyz;

        #my $xyz_new = $x*$cx + $y*$cy + $z*$cz + $trns;
        my $xyz_new = V( $xr, $yr, $zr ) + $trns;
        $atom->set_coords( $tf, $xyz_new );
    }

}

sub fix_serial {
    my @atoms  = shift->all_atoms;
    my $offset = shift;
    $offset = 1 unless defined($offset);
    $atoms[$_]->{serial} = $_ + $offset foreach ( 0 .. $#atoms );
    return $offset;
}

sub print_xyz_ts {
    _print_ts( 'print_xyz', @_ );
}

sub print_pdb_ts {
    _print_ts( 'print_pdb', @_ );
}

sub _print_ts {

    #use one sub for xyz_ts and pdb_ts writing
    my $print_method = shift;

    # two args: \@ts, optional filename
    my $self = shift;
    my $ts   = shift;
    unless ( defined($ts) ) {
        croak "must pass arrayref containing ts";
    }
    my @ts = @$ts;
    unless ( scalar(@ts) ) {
        croak "must pass array with atleast one t";
    }
    my $tmax = $self->tmax;
    my $nt = grep { $_ > $tmax } @ts;
    croak "$nt ts out of bounds" if $nt;

    my $tnow = $self->what_time;

    # take the first out of the loop to setup fh
    $self->gt( shift @ts );
    my $fh = $self->$print_method(@_);

    foreach my $t (@ts) {
        $self->gt($t);
        $fh = $self->$print_method($fh);
    }

    # return to original t
    $self->gt($tnow);
}

sub bin_this {

    #return hash{$_->method}++
    my $self   = shift;
    my $method = shift;

    return ( {} ) unless $self->count_atoms;

    my @atoms = $self->all_atoms;

    # just test the first one...
    croak "Atom does not do $method" unless $atoms[0]->can($method);

    my $bin;
    $bin->{$_}++ foreach ( map { $_->$method } @atoms );
    return ($bin);

}

sub tmax {

    # still not the best implementation! what about atoms without coords?
    my $self = shift;
    my $tbin = $self->bin_this('count_coords');
    my @ts   = keys(%$tbin);
    croak "tmax differences within group" if ( scalar(@ts) > 1 );
    $ts[0] ? return $ts[0] - 1 : return 0;
}

sub what_time {
    my $self = shift;
    my $tbin = $self->bin_this('t');
    my @ts   = keys(%$tbin);
    carp "what_time> t differences within group!!" if ( scalar(@ts) > 1 );
    return $ts[0];
}

sub string_xyz {
    my $self              = shift;
    my $add_info_to_blank = shift;

    my $string;
    $string .= $self->count_atoms . "\n" unless $self->qcat_print;
    $string .= $add_info_to_blank if ( defined($add_info_to_blank) );
    $string .= "\n";

    foreach my $at ( $self->all_atoms ) {
        $string .= sprintf( "%3s %10.6f %10.6f %10.6f\n",
            $at->symbol, @{ $at->get_coords( $at->t ) } );
    }
    return $string;
}

sub print_xyz {
    my $self = shift;
    my $fh   = _open_file_unless_fh(shift);

    print $fh $self->string_xyz;

    # my @atoms = $self->all_atoms;
    #print $fh $self->count_atoms . "\n\n" unless $self->qcat_print;
    #foreach my $at (@atoms) {
    #    printf $fh (
    #        "%3s %10.6f %10.6f %10.6f\n",
    #        $at->symbol, @{ $at->get_coords( $at->t ) }
    #    );
    #}

    return ($fh);    # returns filehandle for future writing

}

sub string_pdb {
    my $self = shift;

    my $t     = $self->what_time;
    my @atoms = $self->all_atoms;

    my $string;
    $string .= sprintf( "MODEL       %2i\n", $t + 1 ) unless $self->qcat_print;

    my $atform =
      "%-6s%5i  %-3s%1s%3s%2s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n";

    foreach my $at (@atoms) {

        # front pad one space if name length is < 4
        my $form = $atform;
        if ( length $at->name > 3 ) {
            $form =
      "%-6s%5i %4s%1s%3s%2s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n";
        }
        $string .= sprintf(
            $form,
            (
                map { $at->$_ }
                  qw (
                  record_name
                  serial
                  name
                  altloc
                  resname
                  chain
                  resid
                  icode
                  )
            ),
            @{ $at->get_coords( $at->t ) },
            $at->occ,
            $at->bfact,
            $at->segid,
            $at->symbol,    # $at->charge
        );
    }
    $string .= "ENDMDL\n" unless $self->qcat_print;
    return $string;
}

sub print_pdb {
    my $self = shift;
    my $fh   = _open_file_unless_fh(shift);

    print $fh $self->string_pdb;

#   my @atoms = $self->all_atoms;
#   printf $fh ( "MODEL       %2i\n", $atoms[0]->t + 1 ) unless $self->qcat_print;
#   my $atform = "%-6s%5i  %-3s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n";

#   foreach my $at (@atoms) {
#       # front pad one space if name length is < 4
#       my $form = $atform;
#       if (length $at->name > 3){
#         $form = "%-6s%5i %4s%1s%3s %1s%4i%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s\n"
#       }
#       printf $fh (
#           $form,
#           (
#               map { $at->$_ }
#                 qw (
#                 record_name
#                 serial
#                 name
#                 altloc
#                 resname
#                 chain
#                 resid
#                 icode
#                 )
#           ),
#           @{ $at->get_coords( $at->t ) },
#           $at->occ,
#           $at->bfact,
#           $at->segid,
#           $at->symbol,    # $at->charge
#       );

    #   }
    #   print $fh "ENDMDL\n" unless $self->qcat_print;

    return ($fh);    # returns filehandle for future writing
}

sub _open_file_unless_fh {

    my $file = shift;    # could be file or filehandle

    my $fh = \*STDOUT;   # default to standard out
                         # if argument is passed, check if filehandle
    if ( defined($file) ) {
        if ( ref($file) ) {
            if ( reftype($file) eq "GLOB" ) {
                $fh = $file;
            }
            else {
                croak "trying write to reference that is not a GLOB";
            }
        }
        else {
            carp "overwrite $file" if ( -e $file );
            $fh = FileHandle->new(">$file");
        }
    }

    return ($fh);
}

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

    use HackaMol::AtomGroup;
    use HackaMol::Atom;

    my $atom1 = HackaMol::Atom->new(
        name    => 'O1',
        coords  => [ V( 2.05274, 0.01959, -0.07701 ) ],
        Z       => 8,
    );
    
    my $atom2 = HackaMol::Atom->new(
        name    => 'H1',
        coords  => [ V( 1.08388, 0.02164, -0.12303 ) ],
        Z       => 1,
    );
    
    my $atom3 = HackaMol::Atom->new(
        name    => 'H2',
        coords  => [ V( 2.33092, 0.06098, -1.00332 ) ],
        Z       => 1,
    );
    
    $atom1->push_charges(-0.834);
    $_->push_charges(0.417) foreach ($atom1, $atom2);
    
    # instance of class that consumes the AtomGroupRole 
    
    my $group = HackaMol::AtomGroup->new(atoms=> [$atom1,$atom2,$atom3]);
    
    print $group->count_atoms . "\n"; #3
    print $group->total_charge . "\n"; # 0
    print $group->total_mass . "\n";  
    
    my @atoms = $group->all_atoms;
    
    print $group->dipole_moment . "\n";
    
    $group->do_forall('push_charges',0);
    $group->do_forall('push_coords',$group->COM);

    $group->gt(1); # same as $group->do_forall('t',1);
    
    print $group->dipole_moment . "\n";
    print $group->bin_atoms_name . "\n";
    print $group->unique_atoms . "\n";

    $group->translate(V(10,0,0));

    $group->rotate( V(1,0,0),
                         180,
                    V(0,0,0));
  
    $group->print_xyz ; #STDOUT

    my $fh = $group->print_xyz("hackagroup.xyz"); #returns filehandle
    $group->print_xyz($fh) foreach (1 .. 9);     # boring VMD movie with 10 frames

=head1 DESCRIPTION

The HackaMol AtomGroupRole class provides core methods and attributes for 
consuming classes that use groups of atoms. The original implementation of 
this role relied heavily on attributes, builders, and clearers.  Such an approach
naturally gives fast lookup tables, but the ability to change atoms and coordinates
made the role to difficult.  Such an approach may be pursued again (without changing
the API) in the future after the API has matured.  The AtomGroupRole calculates all
values for atoms using their own t attributes.

=array_method push_atoms, get_atoms, set_atoms, all_atoms, count_atoms, clear_atoms

ARRAY traits for the atoms attribute, respectively: push, get, set, elements, count, clear

=array_method push_atoms, unshift_atoms

push atom on to the end of the atoms array 
  or 
unshift_atoms on to the front of the array

  $group->push_atoms($atom1, $atom2, @otheratoms);
  $group->unshift_atoms($atom1, $atom2, @otheratoms); # maybe in reverse

=array_method all_atoms

returns array of all elements in atoms array

  print $_->symbol, "\n" foreach $group->all_atoms; 

=array_method get_atoms

return element by index from atoms array

  print $group->get_atoms(1); # returns $atom2 from above

=array_method set_atoms

set atoms array by index

  $group->set_atoms(1, $atom1);

=array_method count_atoms

return number of atoms in group 
  
  print $group->count_atoms; 

=array_method clear_atoms

clears atoms array

=method do_for_all

pass method and arguments down to atoms in group

  $group->do_for_all('t',1); #sets t to 1 for all atoms

=method gt
  
integer argument. wraps do_for_all for setting time within group

  $group->gt(1);

=method dipole

no arguments. return dipole calculated from charges and coordinates as Math::Vector::Real object  

=method COM

no arguments. return center of mass calculated from masses and coordinates as Math::Vector::Real object  
 
=method COZ

no arguments. return center of nuclear charge calculated from Zs and coordinates as Math::Vector::Real object  

=method total_charge

no arguments. return sum of atom charges.

=method total_mass

no arguments. return sum of atom masses.

=method total_Z

no arguments. return sum of Zs.

=method dipole_moment

no arguments. returns the norm of the dipole in debye (assuming charges in electrons, AKMA)

=method bin_atoms

Called with no arguments. Returns a hash with a count for each unique atom symbol.

=method count_unique_atoms

no arguments. returns the number of unique atoms 

=method bin_atoms_name

no arguments. returns a string summary of the atoms in the group sorted by decreasing atomic number. For example; OH2 for water or O2H2 for peroxide.

=attr atoms

isa ArrayRef[Atom] that is lazy with public ARRAY traits described in ARRAY_METHODS

=attr qcat_print

isa Bool that has a lazy default value of 0.  if qcat_print, print all atoms coordinates in one go (no model breaks)

=method tmax

return (count_coords-1) if > 0; return 0 otherwise; croaks if not all atoms share the same tmax.

=method translate

requires L<Math::Vector::Real> vector argument. Optional argument: integer tf.  

Translates all atoms in group by the MVR vector.  Pass tf to the translate 
method to store new coordinates in tf rather than atom->t.

=method rotate

requires Math::Vector::Real vector, an angle (in degrees), and a MVR vector 
origin as arguments. Optional argument: integer tf.  

Rotates all atoms in the group around the MVR vector. Pass tf to the translate 
method to store new coordinates in tf rather than atom->t.

=method print_xyz

optional argument: filename or filehandle.  with no argument, prints xyz formatted output to STDOUT. pass 
a filename and an xyz file with that name will be written or overwritten (with warning). pass filehandle 
for continuous writing to an open filehandle.

=method print_xyz_ts

argument: array_ref containing the values of t to be used for printing.  
optional argument: filename or filehandle for writing out to file. For example,

  $mol->print_xyz_ts([0 .. 3, 8, 4], 'fun.xyz');

will write the coordinates for all group atoms at t=0,1,2,3,8,4 to a file, in
that order.

=method print_pdb

same as print_xyz, but for pdb formatted output

=method print_pdb_ts

same as print_xyz_ts, but for pdb formatted output

=method bin_this

argument: Str , return hash_ref of binned $self->Str.

  $hash_ref{$_}++ foreach ( map {$_->$Str} $self->all_atoms );

=method what_time

returns the current setting of t by checking against all members of group.

=method fix_serial

argument, optional: Int, offset for resetting the serial number of atoms.  
Returns the offset. 

  $group->fix_serial(0); # serial starts from zero

=method centered_vector 

calculates least squares fitted vector for the AtomGroup. Returns normalized Math::Vector::Real 
object with origin V(0,0,0). 

  $mvr = $group->centered_vector; # unit vector origin 0,0,0 
  # place two mercury atoms along the vector to visualize the fit 
  my $hg_1 = HackaMol::Atom->new(Z => 80, coords => [$group->center]);
  my $hg_2 = HackaMol::Atom->new(Z => 80, coords => [$group->center + $mvr]);

=head1 SEE ALSO

=for :list
* L<HackaMol::AtomGroup>
* L<HackaMol::Bond>
* L<HackaMol::Angle>
* L<HackaMol::Dihedral>
