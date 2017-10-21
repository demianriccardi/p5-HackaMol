package HackaMol::Roles::JSONRole;

#ABSTRACT: Role for a group of atoms
use Moose::Role;
use Math::Vector::Real;
use JSON::XS;

sub to_json {
  my $self = shift;
  my $pdb  = shift;

  my $grp = {
              atoms => {},
  };

  my $atoms = $grp->{atoms};
  # is there MOP access to attrs from a Role?
  foreach my $atom ($self->all_atoms){
    push @{$atoms->{elements}{Z}},       $atom->Z;
    if ($pdb){
      push @{$atoms->{elements}{resname}}, $atom->resname;
      push @{$atoms->{elements}{occ}}, $atom->occ;
      push @{$atoms->{elements}{bfact}}, $atom->bfact;
      push @{$atoms->{elements}{iatom}}, $atom->iatom;
      push @{$atoms->{elements}{record_name}}, $atom->record_name;
      push @{$atoms->{elements}{serial}}, $atom->serial;
    }

    foreach my $t (0 .. $self->tmax){
      push @{$atoms->{coords}{$t}}, @{$atom->get_coords($t)} 
    }
  }  

  return JSON::XS::encode_json($grp);   
}


sub from_json{
  my $self = shift;
  my $json = shift;

  my $grp = JSON::XS::decode_json($json);  
  my $atoms = delete $grp->{atoms};
  my @atoms = _build_atoms($atoms);
  return HackaMol::Molecule->new(atoms=>\@atoms);
}

sub _parse_elements{
  my $elements = shift;
  #number via Chemical JSON or Z
  my $Zs = exists($elements->{number}) ? $elements->{number} :
               exists($elements->{Z}) ? $elements->{Z} :
               die "number or Z not defined";

  return ( map{HackaMol::Atom->new( Z=> $_ ) } @$Zs);
}

sub _parse_coords{
  my $coords = shift;
  # chemical json, "3d" or "3d fractional" 
  die "fractional coordinates not supported yet" if exists($coords->{"3d fractional"});  
  my @mvrs; # i -> atoms, j->t
  
  if (exists($coords->{"3d"})){
    @mvrs = _3ds_mvrs(delete($coords->{"3d"}));
  }

  foreach my $t (sort{$a <=> $b} keys %$coords  ) {
    #n, 3-dimensional mvrs
    my $xyz = $coords->{$t};
    my @l_mvrs = _3ds_mvrs($coords->{$t});
    push @{$mvrs[$_]},$l_mvrs[$_] foreach 0 .. $#l_mvrs; 
  }
  return @mvrs;
}

sub _3ds_mvrs {
  my $ar = shift;
  die "dimension problem " if scalar(@$ar) % 3;
  my $natoms = scalar(@$ar)/ 3;
  my @mvrs;
  foreach my $n ( 0 .. $natoms -1){
    my $ind = 3*$n;
    push @mvrs, [ V( map{ $ar->[$ind+$_] } 0, 1, 2 )];
  }
  return @mvrs;  
}


sub _build_atoms{
  my $atoms = shift;
  my @atoms = _parse_elements( delete $atoms->{elements} );
  my @atoms_mvrs = _parse_coords( delete $atoms->{coords} );
  foreach my $i (0 .. $#atoms){
    $atoms[$i]->{ coords } = $atoms_mvrs[$i];
  }
  return @atoms;
 
}

no Moose::Role;

1;

__END__
