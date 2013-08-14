package AtomsGroup;
use Moose::Role;
use Carp;
use MooseX::Storage;
with Storage('io' => 'StorableFile'); 

has 'atoms'   =>  (   
                    traits   => [ 'Array' ]      ,
                    isa      => 'ArrayRef[Atom]' ,
                    default  => sub{[]}          ,
                    handles  =>
                    {
                         add_atoms => 'push'     ,
                         get_atoms => 'get'      ,
                         set_atoms => 'set'      ,
                      delete_atoms => 'delete'   ,
                         all_atoms => 'elements' ,
                       count_atoms => 'count'    ,
                       clear_atoms => 'clear'    ,
                    },
                   );

has 'dipole_moment'  , is => 'rw', isa => 'Num';
has 'dipole'         , is => 'rw', isa => 'ArrayRef|Object';
has 'com'            , 

has 'atomtype_bin'  => (
                         is        => 'rw', 
                         isa       => 'HashRef', 
                         default   => sub{{}},
                         clearer   => 'clear_atomtype_bin',
                         predicate => 'has_atomtype_bin',
                       );

after 'add_groups' => sub {
    my $self = shift;
    $self->countbin_atoms;
    my $atoms = $self->atoms; #this need to be cleaned up, maybe an attribute with all atoms?      
do i need this!!!     
    $atoms->[$_]->iatom($_) foreach (0 .. $#{$atoms});
};  #should be run afterward

after 'set_groups' => sub {
    my $self = shift;
    $self->atomtype_bin({});
    $self->countbin_atoms;
    my $atoms = $self->atoms; #this need to be cleaned up, maybe an attribute with all atoms?      
    $atoms->[$_]->iatom($_) foreach (0 .. $#{$atoms});
};

sub countbin_atoms  
{
  my $self=shift;
  my $syms = $self->atoms_attr('symbol');
  $self->atomtype_bin({}) unless $self->has_atomtype_bin;
  foreach my $sym (@$syms){
    $self->atomtype_bin->{$sym}++;
  }
  $self->natoms( scalar(@$syms) );
}

sub add_to_basis_atoms
{
  my ($self,$atom) = @_;
  push @{$self->basis_atoms}, $atom if (exists $self->atomtype_bin->{$atom->symbol}); 
}

sub set_basis_atoms{
  my ($self,$atoms) = @_;
  my @ats  =  grep{ exists ($self->atomtype_bin->{$_->symbol})} @$atoms;
  $self->basis_atoms(\@ats);
}

sub canonical_name {
use Atom qw(_symbol_to_Z);
#H2O to OH2  CH3CH2CH2CH3-> C4H10
# should this be C4_H10?  HgH or HHg... the question is whether there is a regex way to back out.
# underscore seems convenient
# sort by number seen ... or  sort by Z!
  my $self = shift;
  my %atombin = %{$self->atomtype_bin};

  #my @a = map{$_.$atombin{$_}} 
  #           sort{ 
  #                if   ($atombin{$a} == $atombin{$b}){return $a cmp $b}
  #                else {return $atombin{$a} <=> $atombin{$b}}
  #               } keys(%atombin);

  my @a = map{$_.$atombin{$_}} sort{Atom::_symbol_to_Z($b) <=> Atom::_symbol_to_Z($a)} keys(%atombin);


  return join("_",@a);

}

 __PACKAGE__->meta->make_immutable;

1;
