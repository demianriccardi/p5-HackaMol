package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/HackaMol','lib/roles';
use AtomsGroup;
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecMVRRole',
'BondsAnglesDihedralsRole';

extends 'AtomsGroup';

has 'atomgroups'  => (
    traits   => ['Array'],
    is       => 'ro',
    isa      => 'ArrayRef[AtomsGroup]',
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
                        AtomsGroup->new(atoms=>$group{$_})
                       } sort keys (%group);
                      # } sort{$a<=>$b} keys (%group);

  $self->push_groups(@atomsgroups);
  
}

1;

__END__

=pod

=head1 NAME

Molecule

=cut
