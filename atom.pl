use Modern::Perl;
use Math::Vector::Real;
use HackaMol;
use KiokuDB;
use KiokuDB::Backend::DBI;
use KiokuDB::TypeMap::Entry::Callback;

my $mvr_tp = KiokuDB::TypeMap->new(
    isa_entries => {
        'Math::Vector::Real' => KiokuDB::TypeMap::Entry::Callback->new(
          collapse => sub {my $v   = shift; Math::Vector::Real->new(@{$v})},
          expand   => sub {my $mvr = shift; return (unbless($mvr)) }, 
        ),
    },
);

my $dir = KiokuDB->new(
    backend => KiokuDB::Backend::DBI->new, #"dbi:SQLite:dbname=kiokudb_tutorial.db",
    typemap => $mvr_tp, 
);



$dir->connect(
    "dbi:SQLite:dbname=kiokudb_tutorial.db",
    create => 1, # this causes the tables to be created
);

#'Math::Vector::Real' => KiokuDB::TypeMap::Entry::Callback->new(
#    collapse => sub {my $v   = shift; Math::Vector::Real->new(@{$v})},
#    expand   => sub {my $mvr = shift; return (unbless($mvr)) }, 
#);




my $atom = HackaMol::Atom->new(Z=>1, coords=>[V(1,1,2)]);
my $scope = $dir->new_scope;
$dir->store( shit => $atom );

exit;
$atom->store('atom.json');
print $atom->dump;

my $atom_again = HackaMol::Atom->load('atom.json');
print $atom_again->dump;

