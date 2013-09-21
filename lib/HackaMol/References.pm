package references;

# ABSTRACT: some useful references
use Moose;
use namespace::autoclean;

has '_refs' => (
    traits  => ['Array'],
    isa     => 'ArrayRef',
    default => sub { [] },
    handles => {
        "add_refs"   => 'push',
        "get_refs"   => 'get',
        "all_refs"   => 'elements',
        "count_refs" => 'count',
    },
    lazy => 1,
);

sub build_references {
    my $self = shift;
    $self->add_refs($_) foreach qw(http://perlmol.org);
    $self->add_refs(
        '@article{o2011open,
                    title={Open Babel: An open chemical toolbox},
                    author={O’Boyle, Noel M and Banck, Michael and James, Craig A and Morley, Chris and Vandermeersch, Tim and Hutchison, Geoffrey R},
                    journal={Journal of cheminformatics},
                    volume={3},
                    number={1},
                    pages={1--14},
                    year={2011},
                    publisher={Springer}
                   }'
    );
    $self->add_refs(
        '@article{hanwell2012avogadro,
                    title={Avogadro: an advanced semantic chemical editor, visualization, and analysis platform},
                    author={Hanwell, Marcus D and Curtis, Donald E and Lonie, David C and Vandermeersch, Tim and Zurek, Eva and Hutchison, Geoffrey R},
                    journal={Journal of cheminformatics},
                    volume={4},
                    number={1},
                    pages={1--17},
                    year={2012},
                    publisher={Springer}
                   }'
    );
}

__PACKAGE__->meta->make_immutable;

1;

__END__
need to expand...  maybe drop refs into YAML, JSON... actually that would be a
good module to have...

