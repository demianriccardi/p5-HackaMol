package HackaMol::FileFetchRole;

#ABSTRACT: Role for using LWP::Simple to fetch files from www 
use Moose::Role;
use Carp;

has 'pdbserver',   is => 'rw', isa => 'Str', lazy => 1, default => 'http://pdb.org/pdb/files/';

sub fetch_pdbid{
  #return array of lines from pdb downloaded from pdb.org
  use LWP::Simple;
  my $self = shift;
  my $pdbid = shift;
  $pdbid =~ s/\.pdb//; #just in case
  $pdbid .= '.pdb';
  my $pdb = get($self->pdbserver.$pdbid);
  return ( $pdb );

}

no Moose::Role;
1;

__END__

=head1 SYNOPSIS

   use HackaMol;

   my $pdb = $HackaMol->new->fetch_pdbid("2cba");
   print $pdb;

=head1 DESCRIPTION

FileFetchRole provides attributes and methods for pulling files from the internet.
Currently, the Role has one method and one attribute for interacting with the Protein Database.

=method fetch_pdbid 

fetches a pdb from pdb.org and returns the file in a string.
     
=attr  pdbserver  

isa lazy rw Str that defaults to http://pdb.org/pdb/files/

=head1 SEE ALSO

=for :list
* L<http://www.pdb.org>
* L<LWP::Simple>
                              
