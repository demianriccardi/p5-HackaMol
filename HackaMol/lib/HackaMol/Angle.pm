package Angle;
#ABSTRACT: Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Vector::Real;
with Storage( 'io' => 'StorableFile' ),'AtomsGroupRole';

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 0 ,
            lazy    => 1 ,
            clearer => "clear_$_",
            predicate => "has_$_",
          ) foreach qw(ang_fc ang_eq);

sub ang_normvec{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  my $ang  = $self->ang;
  return V(0,0,0) if ($ang == 0 or $ang == 180);
  my $vec1 = $atoms[1]->inter_dcoords($atoms[0]);
  my $vec2 = $atoms[1]->inter_dcoords($atoms[2]);
  my $v1xv2 = $vec1 x $vec2;
  return ($v1xv2->versor);
}

sub ang{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[1]->angle($atoms[0],$atoms[2]));
}

sub BUILD {
    my $self = shift;
    #atoms know about the angles they are in
    $_->push_angles($self) foreach $self->all_atoms;
}

has 'angle_energy_func' => (
    is      => 'ro',
    isa     => 'CodeRef',
    builder => "_build_angle_energy_func",
    lazy    => 1,
);


sub _build_angle_energy_func {
    #my $self = shift; #self is passed by moose, but we don't use it here
    my $subref = sub {
        my $angle = shift;
        my $val = ( $angle->ang - $angle->ang_eq )**2;
        return ($angle->ang_fc*$val);
    };
    return ($subref);
}

sub angle_energy {
    my $self  = shift;
    return (0) unless ($self->ang_fc > 0);
    my $energy = &{$self->angle_energy_func}($self,@_);
    return ($energy);
}

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
