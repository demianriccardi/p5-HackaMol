package Dihedral;
#ABSTRACT: Dihedral Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Trig;
with Storage( 'io' => 'StorableFile' ),'AtomsGroupRole';

has $_ => (
            is  => 'rw'  ,
            isa => 'Num' ,
            default => 0 ,
            lazy    => 1 ,
            clearer => "clear_$_",
          ) foreach qw(dihe_fc dihe_dphase dihe_eq dihe_multi);

sub dihe{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->dihedral($atoms[1],$atoms[2],$atoms[3]));
}

sub torsion_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $tfunc = 1 + cos($self->dihe_multi*$self->dihe - $self->dihe_dphase);
    return ($self->dihe_fc*$tfunc);
}

sub improper_dihe_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $angsd = ($self->dihe - $self->dihe_eq)**2 ;
    return ($self->dihe_fc*$angsd);
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
