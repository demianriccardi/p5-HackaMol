package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/HackaMol','lib/roles';
use AtomGroup;
use Carp;
use Math::Trig;
use Scalar::Util qw(refaddr);
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecMVRRole',
'BondsAnglesDihedralsRole','QmRole';

extends 'AtomGroup';

has 'atomgroups'  => (
    traits   => ['Array'],
    is       => 'ro',
    isa      => 'ArrayRef[AtomGroup]',
    default  => sub { [] },
    lazy     => 1,  
    handles  => {
        "has_groups"   => 'count',
        "push_groups"   => 'push',
        "get_groups"   => 'get',
        "set_groups"   => 'set',
        "all_groups"   => 'elements',
        "count_groups" => 'count',
        "delete_groups" => 'delete',
        "clear_groups" => 'clear',
    },
);

after 't' => sub {
  my $self = shift;
  $self->gt(@_) if (@_); # set t for all in group
};

sub _build_mass {
  my $self = shift;
  my $mass = 0;
  $mass += $_->mass foreach $self->all_atoms;
  return ($mass); 
}

sub push_groups_by_atom_attr {

  my $self = shift;
  my $attr = shift;

  my %group;
  foreach my $atom ($self->all_atoms){
    push @{$group{$atom->$attr}},$atom;
  }

  my @atomsgroups = map{
                        AtomGroup->new(atoms=>$group{$_})
                       } sort keys (%group);
                      # } sort{$a<=>$b} keys (%group);

  $self->push_groups(@atomsgroups);
  
}

sub dihedral_rotate {
  my $self = shift;
  croak "pass Dihedral, rotation angle (deg), atoms to rotate" unless @_ == 3;
  my $t = $self->t;
  my ($dihe,$ang,$atoms) = @_;
  my ($atom0, $ratom1, $ratom2, $atom3) = $dihe->all_atoms;
  if (2 == scalar( grep{ refaddr($atom0) == refaddr($_) or 
                         refaddr($atom3) == refaddr($_) } @{$atoms} )){ 
    croak "will not rotate atoms on both sides of dihedral " ;
  } 
  my $rvec = ($ratom2->inter_dcoords($ratom1))->versor;
  my @cor  = map{$_->get_coords($t)-$ratom1->xyz} @$atoms; #shift origin too
  my @rcor = $rvec->rotate_3d( deg2rad($ang), @cor );
#shift origin back
  $atoms->[$_]->set_coords($t, $rcor[$_]+$ratom1->xyz) foreach 0 .. $#rcor; 

}


1;

__END__

=pod

=head1 NAME

Molecule

=cut
