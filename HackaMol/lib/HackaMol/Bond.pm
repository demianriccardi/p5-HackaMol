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

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 0 ,
            lazy    => 1 ,
            clearer => "clear_$_",
            predicate => "has_$_",
          ) foreach qw(bond_fc bond_length_eq);

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

sub bond_energy {
    my $self  = shift;
    return (0) unless ($self->bond_fc > 0);
    my $bndsd = ( $self->bond_length - $self->bond_length_eq )**2;
    return ($self->bond_fc*$bndsd);
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
