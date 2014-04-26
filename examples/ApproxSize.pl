# upper limit size of a collection of atoms
# 
use Modern::Perl;
use HackaMol;
use Math::Vector::Real;

my $hack = HackaMol->new(data=>'t/lib');
foreach my $fxyz (grep {!m/bad/} grep {m/\.xyz/} $hack->data->children){
  say $fxyz;
  my @mvr  = map{$_->xyz} $hack->read_file_atoms($fxyz);
  my ($bot,$top) = Math::Vector::Real->box(@mvr);
  printf ("$fxyz %10.2f\n", abs($bot-$top));
  
}

