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
            is      => 'rw'  ,
            isa     => 'Math::Vector::Real' ,
            clearer => 'clear_bond_vector',
            builder => '_build_bond_vector',
            lazy    => 1,
          ) foreach qw(bond_vector);

sub BUILD {
    my $self = shift;
    $_->push_bonds($self) foreach $self->all_atoms;
}

sub _build_bond_vector{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->inter_dcoords($atoms[1]));
}

has $_ => (
            is      => 'rw'  ,
            isa     => 'Num' ,
            clearer => 'clear_bond_length',
            builder => '_build_bond_length',
            lazy    => 1,
          ) foreach qw(bond_length);

sub _build_bond_length{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->distance($atoms[1]));
}

before 'bond_length' => sub {
    my $self = shift;
    if (grep {$_->is_dirty} $self->all_atoms){
      $self->clear_bond_length;
    }
};

before 'bond_vector' => sub {
    my $self = shift;
    if (grep {$_->is_dirty} $self->all_atoms){
      $self->clear_bond_vector;
    }
};


sub _clear_group_attrs {
    my $self = shift;
    foreach my $clearthis (qw(clear_dipole clear_COM clear_COZ
                              clear_dipole_moment clear_total_charge
                              clear_total_mass clear_total_Z 
                              clear_atoms_bin clear_bond_length
                              clear_bond_vector clear_bond_order)){
      $self->$clearthis;
    }
}


__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
