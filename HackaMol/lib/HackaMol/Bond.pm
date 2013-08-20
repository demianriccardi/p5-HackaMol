package Bond;
#ABSTRACT: Bond class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' ),'AtomsGroupRole';

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 1 ,
            lazy    => 1 ,
            clearer => 'clear_bond_order',
          ) foreach qw(bond_order);

sub BUILD {
    my $self = shift;
    # atoms know about bonds they have
    $_->push_bonds($self) foreach $self->all_atoms;
}

sub bond_vector{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->inter_dcoords($atoms[1]));
}

sub bond_length{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->distance($atoms[1]));
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
