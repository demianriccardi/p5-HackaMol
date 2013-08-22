package AtomGroup;
#ABSTRACT: HackaMol AtomGroup class 
use Moose;
use lib 'lib/roles';
use Carp;
use MooseX::Storage;
with Storage( 'io' => 'StorableFile' ),'AtomGroupRole';

has 'gname' => (
    is  => 'rw',
    isa => 'Str',
);

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

__PACKAGE__->meta->make_immutable;

1;

__END__

=pod

=head1 NAME

AtomGroup - HackaMol AtomGroup class 

=head1 NAME

Bond

=cut
