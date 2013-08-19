package Dihedral;
#ABSTRACT: Dihedral Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Vector::Real;
with Storage( 'io' => 'StorableFile' ),'AtomsGroupRole';

has $_ => (
            is      => 'rw'  ,
            isa     => 'Num' ,
            clearer => "clear_$_",
            builder => "_build_$_",
            lazy    => 1,
          ) foreach qw(dihe);

before 'dihe' => sub {
    my $self = shift;
    if (grep {$_->is_dirty} $self->all_atoms){
      $self->clear_dihe;
    }
};

sub _build_dihe{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->dihedral($atoms[1],$atoms[2],$atoms[3]));
}

sub _clear_group_attrs {
    my $self = shift;
    foreach my $clearthis (qw(clear_dipole clear_COM clear_COZ
                              clear_dipole_moment clear_total_charge
                              clear_total_mass clear_total_Z 
                              clear_atoms_bin 
                              clear_dihe)) {
      $self->$clearthis;
    }
}


sub BUILD {
    my $self = shift;
    $_->push_dihedrals($self) foreach $self->all_atoms;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
