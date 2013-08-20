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
        "break_groups" => 'delete',
        "clear_groups" => 'clear',
    },
);

after 'gt' => sub {
  my $self = shift;
  $_->_clear_group_attrs foreach $self->all_groups;
  $_->_clear_group_attrs foreach $self->all_bonds;
  $_->_clear_group_attrs foreach $self->all_angles;
  $_->_clear_group_attrs foreach $self->all_dihedrals;
}; 

after 'push_bonds' => sub {shift->gt};

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
