package BondsAnglesDihedrals;
# ABSTRACT: Array traits for containers of HackaMol Bonds, Angles, Dihedrals.
use Moose::Role;
use Carp;

has 'bonds'  => (
    traits   => ['Array'],
    isa      => 'ArrayRef[Bond]',
    default  => sub { [] },
    lazy     => 1,
    weak_ref => 1,
    handles  => {
        "has_bonds"   => 'count',
        "add_bonds"   => 'push',
        "get_bonds"   => 'get',
        "set_bonds"   => 'set',
        "all_bonds"   => 'elements',
        "count_bonds" => 'count',
        "break_bonds" => 'delete',
        "clear_bonds" => 'clear',
    },
); 

has 'angles'  => (
    traits   => ['Array'],
    isa      => 'ArrayRef[Angle]',
    default  => sub { [] },
    lazy     => 1,
    weak_ref => 1,
    handles  => {
        "has_angles"   => 'count',
        "add_angles"   => 'push',
        "get_angles"   => 'get',
        "set_angles"   => 'set',
        "all_angles"   => 'elements',
        "count_angles" => 'count',
        "break_angles" => 'delete',
        "clear_angles" => 'clear',
    },
);

has 'dihedrals'  => (
    traits   => ['Array'],
    isa      => 'ArrayRef[Dihedral]',
    default  => sub { [] },
    lazy     => 1,
    weak_ref => 1,
    handles  => {
        "has_dihedrals"   => 'count',
        "add_dihedrals"   => 'push',
        "get_dihedrals"   => 'get',
        "set_dihedrals"   => 'set',
        "all_dihedrals"   => 'elements',
        "count_dihedrals" => 'count',
        "break_dihedrals" => 'delete',
        "clear_dihedrals" => 'clear',
    },
);

no Moose::Role;

1;

__END__

=pod

