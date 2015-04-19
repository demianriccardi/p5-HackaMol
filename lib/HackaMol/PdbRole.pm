package HackaMol::PdbRole;

#ABSTRACT: PdbRole of lazy attributes for HackaMol atoms
use Moose::Role;
use LWP::Simple;
use Carp;


#something like this: $mol->fetch_pdbid(); PDBRole may be the best place for it 
sub fetch_pdbid
{
	my $self=shift;
	my $pdbid=shift;
	my $file;
	if(length($pdbid)==4)
	{
		$file=$pdbid;
	}
	else
	{
		$file=$pdbid.'.pdb';
	}
	if(-f $file)
	{
		return 1;
	}
	else
	{
		getstore("http://pdb.org/pdb/files/$file","$file")  or carp "please check the pdbid:$pdbid";
	}
	
	
}

sub fetch_pdbids
{
	my $self=shift;
	my $pdbids=shift;
	foreach my $pdbid(@pdbids)
	{
	$self->fetch_pdbid($pdbid);
  }
	
}


my %aa321 = (
  ALA=>'A', ARG=>'R', ASN=>'N', ASP=>'D', CYS=>'C',
  GLU=>'E', GLN=>'Q', GLY=>'G', HIS=>'H', ILE=>'I',
  LEU=>'L', LYS=>'K', MET=>'M', PHE=>'F', PRO=>'P',
  SER=>'S', THR=>'T', TRP=>'W', TYR=>'Y', VAL=>'V',
  SEC=>'U',
);

sub aa321 {
  my $resname = uc($_[0]->resname);
  unless ( exists($aa321{$resname}) ){
    carp "PDBRole> residue $resname name has no 1 letter code; return X";
    return ('X');
  }
  return ($aa321{$resname});
}

has 'record_name', is => 'rw', isa => 'Str', lazy => 1, default => 'HETATM';
has 'occ',         is => 'rw', isa => 'Num', lazy => 1, default => 1.0;
has 'bfact',       is => 'rw', isa => 'Num', lazy => 1, default => 20.0;
has 'resname',     is => 'rw', isa => 'Str', lazy => 1, default => 'UNK';
has 'chain',       is => 'rw', isa => 'Str', lazy => 1, default => ' ';
has 'altloc',      is => 'rw', isa => 'Str', lazy => 1, default => ' ';
has 'icode',       is => 'rw', isa => 'Str', lazy => 1, default => ' ';
has 'pdbid',       is => 'rw', isa => 'Str', lazy => 1, default => ' ';
has 'segid',       is => 'rw', isa => 'Str', lazy => 1, default => ' ';
has 'iatom',       is => 'rw', isa => 'Int', lazy => 1, default => 0;
has 'resid',       is => 'rw', isa => 'Int', lazy => 1, default => 1;
has 'serial',      is => 'rw', isa => 'Int', lazy => 1, default => 1;

no Moose::Role;
1;

# added attributes for parsing a PDB file
# Histidine 64 for Human Carbonic Anhydrase II from pdbid: 2CBA
#some lines for reference that did not translate well into POD
#my $lines_H64_2cba = '
#         1         2         3         4         5         6         7         8
#12345678901234567890123456789012345678901234567890123456789012345678901234567890
#ATOM    504  N   HIS A  64      -0.822  -1.995   6.439  1.00 10.63           N
#ATOM    505  CA  HIS A  64      -0.847  -2.773   7.721  1.00 10.22           C
#ATOM    506  C   HIS A  64      -2.270  -3.278   7.949  1.00  9.51           C
#ATOM    507  O   HIS A  64      -2.449  -4.225   8.736  1.00 11.94           O
#ATOM    508  CB AHIS A  64      -0.335  -2.173   8.991  0.70 10.98           C
#ATOM    509  CB BHIS A  64      -0.409  -1.888   8.917  0.20  9.21           C
#ATOM    510  CG AHIS A  64      -0.736  -0.782   9.258  0.70 12.16           C
#ATOM    511  CG BHIS A  64       0.714  -0.964   8.521  0.20  8.70           C
#ATOM    512  ND1AHIS A  64      -1.795  -0.353  10.014  0.70 14.34           N
#ATOM    513  ND1BHIS A  64       0.513   0.338   8.162  0.20  8.57           N
#ATOM    514  CD2AHIS A  64      -0.117   0.344   8.805  0.70 12.79           C
#ATOM    515  CD2BHIS A  64       2.037  -1.198   8.364  0.20  8.90           C
#ATOM    516  CE1AHIS A  64      -1.809   0.971  10.061  0.70 13.80           C
#ATOM    517  CE1BHIS A  64       1.691   0.868   7.814  0.20  8.65           C
#ATOM    518  NE2AHIS A  64      -0.844   1.403   9.281  0.70 15.45           N
#ATOM    519  NE2BHIS A  64       2.615  -0.026   7.936  0.20  6.94           N
#';
#
#my $ATOM_record_format = '
#Record Format:  http://www.wwpdb.org/documentation/format23/sect9.html
#COLUMNS      DATA TYPE        FIELD      DEFINITION
#------------------------------------------------------
# 1 -  6      Record name      "ATOM    "
# 7 - 11      Integer          serial     Atom serial number.
#13 - 16      Atom             name       Atom name.
#17           Character        altLoc     Alternate location indicator.
#18 - 20      Residue name     resName    Residue name.
#22           Character        chainID    Chain identifier.
#23 - 26      Integer          resSeq     Residue sequence number.
#27           AChar            iCode      Code for insertion of residues.
#31 - 38      Real(8.3)        x          Orthogonal coordinates for X in
#                                         Angstroms
#39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in
#                                         Angstroms
#47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in
#                                         Angstroms
#55 - 60      Real(6.2)        occupancy  Occupancy.
#61 - 66      Real(6.2)        tempFactor Temperature factor.
#77 - 78      LString(2)       element    Element symbol, right-justified.
#79 - 80      LString(2)       charge     Charge on the atom.
#';
#

__END__

=head1 SYNOPSIS

   use HackaMol::Atom;

   my $atom = Atom->new(Z=>6);

   print $atom->$_ foreach ( qw( 
                              record_name  
                              serial       
                              occ          
                              bfact        
                              resname      
                              chain        
                              altloc       
                              resid        
                              iatom        
                              pdbid        
                              segid  
                            )
   );
    
   foreach my $resn (qw(TRP TYR GLN ASN ETC)){
     $atom->resname($resn);
     print $atom->aa321;
   }

=head1 DESCRIPTION

PdbRole provides atom attributes for PDB parsing.  All attributes are 'rw' and
lazy, so they will not contaminate the namespace unless called upon. The
functionality of the PdbRole may be extended in the future.  An extension 
(HackaMolX::PDB or HackaMol::X::PDB) will be released soon.

=method aa321

returns 1 letter code for amino acid. Generated from resname attribute. throws a warning and returns 'X' if residue name is unrecognized.
     
=attr  record_name  

pdb cols 1-6:  ATOM|HETATM.  default = 'HETATM'

=attr  serial       

pdb cols 7-11: index. default = 0

=attr  occ          

pdb cols 55-60: occupancy.  range 0 to 1. default = 1.0  (all there)

=attr  bfact        

pdb cols 61-66: temperature factor. see Willis and Pryor. default = 20.0

=attr  altloc       

pdb col 17. is AHIS and BHIS above. default = ' '

=attr  resname      

pdb cols 18-20: residue name. defaults to 'ALA' 

=attr  chain        

pdb cols 22. protein chain.  default = '  '

=attr  resid        

pdb cols 23-26. residue index. PDB calls is resSeq, but resid is more familiar
(to me). default = 64  

=attr  icode

pdb cols 27. see code comments or
L<http://www.wwpdb.org/documentation/format23/sect9.html>
                            
=attr  pdbid        

a place to store which PDB it came from

=attr  segid 

a CHARMMish parameter that may not belong here. default = 'TIP3'
                              
=attr  iatom        
                              
this is a place for the actual index.  Segid does not have to start from 0 (cut
and paste as above, etc). 
                              
=head1 SEE ALSO

=for :list
* L<http://www.pdb.org>
* L<Bio::PDB::Structure::Atom>
                              
