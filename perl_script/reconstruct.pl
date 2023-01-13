#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 4){
 print "\nFonctionnement:\n> perl process_Genlib.pl <fichier Proband_Haplotypes.txt> <fichier *.hap> <fichier *.map> <taille de la sequence en BP>\n\n";
 exit;
}

my $fichier = $ARGV[0];
my $anc     = $ARGV[1];
my $posFile = $ARGV[2];
my $bpTot   = $ARGV[3]; 


## FICHIER ANC
open ANC, $anc || die "peut pas ouvrir $anc\n";
print "\nProcessing $anc \n";

my %seqAnc=();
while(my $line = <ANC>)
{
 chomp($line);
 my @tmp = split(" ", $line);
 $seqAnc{$tmp[0]} = $tmp[1];
}
print keys( %seqAnc )." seqAnc\n";


## FICHIER POS
print "\nProcessing $posFile\n";
open POS, $posFile || die "peut pas ouvrir $posFile\n";

# On assume que les positions des SNPS sont ordonnees.
my @pos=();
while(my $line = <POS>)
{
 chomp($line);
 push @pos, $line;
}
print @pos." positions\n";

## FICHIER PROBAND
print "\nProcessing $fichier\n";
my $fichierOut = "Probands.hap";

open(FILE, $fichier)  || die "peut pas ouvrir $fichier\n";
open(OUT, ">", $fichierOut)  || die "peut pas ouvrir $fichierOut\n";

my $num=0; my $nbSim; my $nbPro; my $sim=0;

while(my $line = <FILE>)
{
 chomp($line);
 
 if($num==0){
   my @tmp = split(";", $line);
   $nbSim = $tmp[0];
   $nbPro = $tmp[1];
 }
 else{
   
   # separation des 3 parties et on enleve les crochets de debut et fin :
   my @tmp  = split("\\}\\{", $line);
   $tmp[0] =~ s/{//g;
   $tmp[2] =~ s/}//g;
   
   # separe les composantes de chaque partie:
   my @specs = split(";", $tmp[0]);
   my @hap1  = split(";", $tmp[1]);
   my @hap2  = split(";", $tmp[2]);

   # on vient de passer à une nouvelle simulation
   if($sim != $specs[0]) {
    $sim = $specs[0];
    if($sim%50==0){ print "simulation #$sim / $nbSim\n"; }
   }
   
   my $proID = $specs[1];
   
   #Enlever les 0 de depart.
   shift @hap1;
   shift @hap2;
   # ecrire les sequences pour chaque proband, un haplotype par ligne, donc 2 lignes par proband.
   print OUT $sim." ".$proID." ";
   getSequence(@hap1); 
   
   print OUT $sim." ".$proID." ";
   getSequence(@hap2);
 }
 $num++;
}

print "Number of :\n  simulations: $nbSim\n  probands: $nbPro\n";
close(FILE);

sub getSequence{

  my @hap = @_;
  my $deb = 0;
  my $fin;
  my $seqTester = 0;
  my $nbRec = 0;
  
  my $index = @hap/2 - 1;
  if($index == 0) {
    print OUT $seqAnc{$hap[0]}."\n";
  }
  else {
   foreach my $i (0..$index) {
    my $ind1 = $hap[2*$i];
    my $pos1 = $hap[2*$i+1];
    
    if($pos1 == $bpTot) {
     $fin = @pos;
    }
    else {
     my $posRec = $pos1 ;
     # on cherche la position de recombinaison.
     $fin=$deb;
     #la recombinaison peut arriver entre le dernier SNP et la fin de la sequence.
     while($pos[$fin] - $posRec < 0 ){
       $fin++;
       if($fin == (@pos-1)){
         if($pos[$fin] - $posRec < 0){ $fin++; }
         last;
       }
     } 
    }
    $seqTester = $seqTester + length(substr($seqAnc{$ind1}, $deb, ($fin-$deb) ));
    print OUT substr($seqAnc{$ind1}, $deb, ($fin-$deb) );
    $deb=$fin;
   }
   if($seqTester != @pos){ print "Probleme: $seqTester vs ".@pos." @hap\n"; exit}
   print OUT "\n";
  }
}



