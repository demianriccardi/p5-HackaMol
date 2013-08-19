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

sub _build_mass {
  my $self = shift;
  my $mass = 0;
  $mass += $_->mass foreach $self->all_atoms;
  return ($mass); 
}

1;

__END__

=pod

=head1 NAME

Molecule

=cut
