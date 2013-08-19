package AtomsGroup;
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

sub Rg {
 #radius of gyration. 
    my $self         = shift;
    return(0) unless ($self->count_atoms);
    my @atoms        = $self->all_atoms;
    my $com          = $self->COM;
    my $total_mass   = $self->total_mass;
    my @masses = map { $_->mass} @atoms;
    my @dvec2  = map{$_*$_} map { $_->get_coords($_->t) - $com } @atoms;
    my $sum    = 0;
    $sum      += $masses[$_]*$dvec2[$_] foreach 0 .. $#dvec2;
    return( sqrt($sum/$total_mass) );
}

sub _clear_group_stuff => sub {
    my $self = shift;
    foreach my $clearthis (qw(
                              clear_bond_length
                              clear_bond_vector
                              clear_bond_order
                             )){
    $self->$clearthis;
    }
};


__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

Bond

=cut
