package Molecule;
#ABSTRACT: Molecule class for HackaMol
use Moose;
use lib 'lib/roles';
use MooseX::Storage;
with Storage('io' => 'StorableFile'), 'PhysVecMVRRole',
'BondsAnglesDihedralsRole';

has 'atomgroups'  => (
    traits   => ['Array'],
    isa      => 'ArrayRef[AtomsGroup]',
    default  => sub { [] },
    lazy     => 1,
    handles  => {
        "has_atomgroups"   => 'count',
        "push_atomgroups"   => 'push',
        "get_atomgroups"   => 'get',
        "set_atomgroups"   => 'set',
        "all_atomgroups"   => 'elements',
        "count_atomgroups" => 'count',
        "break_atomgroups" => 'delete',
        "clear_atomgroups" => 'clear',
    },
);

sub atoms {
  my $self = shift;
  my @atoms;
  push @atoms,$_->all_atoms  foreach $self->all_atomgroups;
  return @atoms;
}

sub _build_mass {
  my $self = shift;
  my $mass = 0;
  $mass += $_->mass foreach $self->atoms;
  return ($mass); 
}

1;

__END__

=pod

=head1 NAME

Molecule

=cut
