package Angle;
#ABSTRACT: Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Vector::Real;
with Storage( 'io' => 'StorableFile' ),'AtomsGroupRole';

has $_ => (
            is      => 'rw'  ,
            isa     => 'Math::Vector::Real' ,
            clearer => "clear_$_",
            builder => "_build_$_",
            lazy    => 1,
          ) foreach qw(ang_normvec);

sub _build_ang_normvec{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  my $ang  = $self->ang;
  return V(0,0,0) if ($ang == 0 or $ang == 180);
  my $vec1 = $atoms[1]->inter_dcoords($atoms[0]);
  my $vec2 = $atoms[1]->inter_dcoords($atoms[2]);
  my $v1xv2 = $vec1 x $vec2;
  return ($v1xv2->versor);
}

has $_ => (
            is      => 'rw'  ,
            isa     => 'Num' ,
            clearer => "clear_$_",
            builder => "_build_$_",
            lazy    => 1,
          ) foreach qw(ang);

sub _build_ang{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[1]->angle($atoms[0],$atoms[2]));
}

before $_ => sub {
    my $self = shift;
    if (grep {$_->is_dirty} $self->all_atoms){
      my $method = "clear_$_";
      $self->$method;

    }
} foreach (qw(ang ang_normvec));

sub _clear_group_attrs {
    my $self = shift;
    foreach my $clearthis (qw(clear_dipole clear_COM clear_COZ
                              clear_dipole_moment clear_total_charge
                              clear_total_mass clear_total_Z 
                              clear_atoms_bin 
                              clear_ang clear_ang_normvec )){
      $self->$clearthis;
    }
}


sub BUILD {
    my $self = shift;
    $_->push_angles($self) foreach $self->all_atoms;
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
