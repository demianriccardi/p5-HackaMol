package Dihedral;
#ABSTRACT: Dihedral Angle class for HackaMol
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
use Math::Trig;
with Storage( 'io' => 'StorableFile' ),'AtomGroupRole';

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

sub dihe_rad{
  my $self  = shift;
  my @atoms = $self->all_atoms;
  return ($atoms[0]->dihedral_rad($atoms[1],$atoms[2],$atoms[3]));
}

has 'improper_dihe_energy_func' => (
    is      => 'rw',
    isa     => 'CodeRef',
    builder => "_build_improper_dihe_energy_func",
    lazy    => 1,
);


sub _build_improper_dihe_energy_func {
    my $subref = sub {
        my $dihedral = shift;
        my $val = ( $dihedral->dihe - $dihedral->dihe_eq )**2;
        return ($dihedral->dihe_fc*$val);
    };
    return ($subref);
}

has 'torsion_energy_func' => (
    is      => 'rw',
    isa     => 'CodeRef',
    builder => "_build_torsion_energy_func",
    lazy    => 1,
);


sub _build_torsion_energy_func {
    my $subref = sub {
        my $dihedral = shift;
        my $val = 1 + cos($dihedral->dihe_multi*$dihedral->dihe_rad 
                          - $dihedral->dihe_dphase);
        return ($dihedral->dihe_fc*$val);
    };
    return ($subref);
}

sub torsion_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $energy = &{$self->torsion_energy_func}($self,@_);
    return ($energy);
}

sub improper_dihe_energy {
    my $self  = shift;
    return (0) unless ($self->dihe_fc > 0 );
    my $energy = &{$self->improper_dihe_energy_func}($self,@_);
    return ($energy);
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
