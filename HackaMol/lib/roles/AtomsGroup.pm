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

has 'dipole'         , is => 'rw', isa => 'Math::Vector::Real';
has 'com'            , is => 'rw', isa => 'Math::Vector::Real';

has $_ => (
           is      => 'rw',
           isa     => 'Math::Vector::Real',
           builder => "_build_$_",
           clearer => "_clear_$_",
          ) foreach (qw(dipole COM COZ));

sub _build_dipole {
    my $self    = shift;
    my @atoms   = $self->all_atoms;
    my @vectors = map { $_->get_coords($_->t)} @atoms;
    my @charges = map { $_->get_charges($_->t)} @atoms;
    croak "mismatch number of coords and charges" if ($#vectors != $#charges);
    my $dipole  = V(0,0,0);
    $dipole += $vectors[$_]*$charges[$_] foreach 0 .. $#charges;
    return ($dipole);
}

sub _build_COM {
    my $self = shift;

}

sub _build_COZ {
    my $self = shift;
    
}


has $_ => (
           is      => 'rw',
           isa     => 'Num',
           builder => "_build_$_",
           clearer => "_clear_$_",
          ) foreach (qw(dipole_moment total_charge));

sub _build_total_charge {
  my $self    = shift;
  my @atoms   = $self->all_atoms;
  my @charges = map{$_->get_charges($_->t)} @atoms;
  my $sum = 0;
  $sum += $_ foreach @charges;
  return $sum; 
}

sub _build_dipole_moment {
  my $self    = shift;
  my @atoms   = $self->all_atoms;
  my @charges = map{$_->get_charges($_->t)} @atoms;
  my $sum = 0;
  $sum += $_ foreach @charges;
  return $sum;
}



has 'atomtype_bin'  => (
                         is        => 'rw', 
                         isa       => 'HashRef', 
                         default   => sub{{}},
                         clearer   => 'clear_atomtype_bin',
                         predicate => 'has_atomtype_bin',
                       );

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

no Moose::Role;

1;
