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

has 'bond_energy_func' => (
    is      => 'ro',
    isa     => 'CodeRef',
    builder => "_build_bond_energy_func",
    lazy    => 1,
);


sub _build_bond_energy_func {
    #my $self = shift; #self is passed by moose, but we don't use it here
    my $subref = sub {
        my $bond = shift;
        my $val = ($bond->bond_length - $bond->bond_length_eq )**2;
        return ($bond->bond_fc*$val);
    };
    return ($subref);
}

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
    my $energy = &{$self->bond_energy_func}($self);
    return ($energy);
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
