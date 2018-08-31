package HackaMol::Roles::ReadPdbxRole;

# ABSTRACT: parse PDBx/mmCIF files 
use Moose::Role;
use HackaMol::PeriodicTable qw(_element_name _trim _qstring_num);
use Math::Vector::Real;
use Carp;

requires 'readline_func';

sub read_cif_atoms {
    # merge atoms here??? 
    my $self = shift;
    my  ($atoms, $models, $count_check) = $self->_read_cif_atoms(@_);
    my @sets;
    foreach my $model_num (@{$models}){
        push @sets, [grep {$_->model_num == $model_num} @$atoms]
    }
    return @sets;
}

sub _read_cif_atoms {

    # read pdb file and generate list of Atom objects
    # no model branching
    my $self = shift;
    my $fh   = shift;
    my @atoms;
    my %track_models;

    while (<$fh>) {
        if($self->has_readline_func){
            next if $self->readline_func->($_) eq 'PDB_SKIP';
        }
# assuming order of columns, e.g.
#loop_
#_atom_site.group_PDB            (record_name)
#_atom_site.id                   (serial)
#_atom_site.type_symbol          (element)
#_atom_site.label_atom_id        (name)
#_atom_site.label_alt_id         (altloc)
#_atom_site.label_comp_id        (res_name)
#_atom_site.label_asym_id        (chain_id)
#_atom_site.label_entity_id      (entity_id)
#_atom_site.label_seq_id         (res_id)
#_atom_site.pdbx_PDB_ins_code    (icod)
#_atom_site.Cartn_x              (x)
#_atom_site.Cartn_y              (y)
#_atom_site.Cartn_z              (z)
#_atom_site.occupancy            (occ)
#_atom_site.B_iso_or_equiv       (B)
#_atom_site.pdbx_formal_charge   (charge)
#_atom_site.auth_seq_id          (auth_seq_id)
#_atom_site.auth_comp_id         (auth_comp_id)
#_atom_site.auth_asym_id         (auth_asym_id)
#_atom_site.auth_atom_id         (auth_atom_id)
#_atom_site.pdbx_PDB_model_num   (model_num)
#ATOM 1       N N   . PRO A   1 1   ? 393.230 1016.300 385.017 1.0 0.00 ? 1   PRO g8 N   1
        if (/^(?:HETATM|ATOM)/) { 

            my (
                $record_name, $serial,  $element, $name,  $altloc,
                $res_name,     $chain_id, $entity_id, $res_id,  $icod,
                $x,           $y,       $z,       $occ,
                $B,           $charge, $auth_seq_id,   $auth_comp_id, 
                $auth_asym_id, $auth_atom_id, $model_num
            ) = split ; 

            # ignore charge for now...

            $element = ucfirst( lc( $element ) );
            my $xyz = V( $x, $y, $z );

            my $atom = HackaMol::Atom->new(
                name         => $name,
                record_name  => $record_name,
                serial       => $serial,
                chain        => $chain_id,
                symbol       => $element,
                coords       => [$xyz],
                occ          => $occ * 1,
                bfact        => $B * 1,
                resname      => $res_name,
                resid        => $res_id,
                auth_seq_id  => $auth_seq_id,
                auth_comp_id => $auth_comp_id,
                auth_asym_id => $auth_asym_id,
                auth_atom_id => $auth_atom_id,
                entity_id    => $entity_id,
                model_num    => $model_num,
            );
            $atom->altloc($altloc) if ($altloc ne '.');
            $atom->icode($icod) if ($icod ne '?');
            $track_models{$model_num}++;
            push @atoms, $atom;
        }
    }

    my %model_atom_count = reverse %track_models;

    if (keys (%model_atom_count) > 1){
        warn "atom number inconsistencies in models. use third return value\n";
    }

    return (\@atoms, [ sort {$a <=> $b} keys %track_models ], \%model_atom_count);
}

no Moose::Role;

1;

__END__

=head1 SYNOPSIS

   my @atoms = HackaMol->new
                       ->read_cif_atoms("some.cif");

=head1 DESCRIPTION

The HackaMol::Roles::ReadPdbxRole provides methods (such as read_cif_atoms for building HM atoms) for pulling information from PDBx/mmcif.
More info available in cif and more attributes are populated than via ReadPdbRole.

=method read_pdb_atoms

One argument: the filehandle
Returns a list of HackaMol::Atom objects.

The implementation differs from the PDB parser in how models are handled.  Each atom tracks the model number through the model_num attribute;
the read_pdb_atoms adds configurations to each subsequent atom into t (with identical metadata) as model number increases. Still considering the
best way to handle this.

=head1 SEE ALSO

=for :list
* L<HackaMol>
* L<HackaMol::Atom>
* L<HackaMol::Roles::MolReadRole>
* L<Protein Data Bank|http://pdb.org>

