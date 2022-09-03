#!/usr/bin/perl
use strict;
use warnings;

if(@ARGV != 4){
 print "\nFonctionnement:\n> perl process_Genlib.pl <fichier Proband_Haplotypes.txt> <fichier *.hap> <fichier *.map> <taille de la sequence en Mb>\n\n";
 exit;
}

my $fichier = $ARGV[0];
my $anc     = $ARGV[1];
my $posFile = $ARGV[2];
my $bpTot   = $ARGV[3]; # donne en Mb donc on va *1000000

$bpTot = $bpTot*1000000;


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
my $nbRec; my $nbMeioses;# my $nbMatch;

while(my $line = <FILE>)
{
 chomp($line);
 
 if($num==0){
   my @tmp = split(";", $line);
   $nbSim = $tmp[0];
   $nbPro = $tmp[1];
 }
 else{
   # Le format: {1;342;540;5320}{0;69362.2;1}{0;18433.2;0.241942;...;18351.2;1}  :
   #		{#sim;proID;nbRec;nbMeioses}
   #		{composition de l'haplotype1}  : {0;individu1.haplotype;positionFin1;individu2.haplotype;positionFin2...;1}
   #		{composition de l'haplotype2}  : {0;individu1.haplotype;positionFin1;individu2.haplotype;positionFin2...;1}
   
   # separation des 3 parties et on enleve les crochets de debut et fin :
   my @tmp  = split("\\}\\{", $line);
   $tmp[0] =~ s/{//g;
   $tmp[2] =~ s/}//g;
   
   # separe les composantes de chaque partie:
   my @specs = split(";", $tmp[0]);
   my @hap1  = split(";", $tmp[1]);
   my @hap2  = split(";", $tmp[2]);

   # on commence les donnees d'haplotypes. Comme on travaille sur une genealogie a la fois, nbRec et nbMeiose ne change pas.
   if($num==1) { 
    $nbRec    = $specs[2];
    $nbMeioses = $specs[3];
   }
   
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
#   if($proID==408828){ print "@hap1\n"; }
   getSequence(@hap1); #, $sim, $proID);
   
   print OUT $sim." ".$proID." ";
   getSequence(@hap2); #, $sim, $proID);
#   if($sim == 736){ exit; }
#   exit;
 }
 $num++;
}

print "Number of :\n  simulations: $nbSim\n  probands: $nbPro\n";#  matches: $nbMatch\n";
close(FILE);

sub getSequence{
#  my $proID = pop(@_);
#  my $sim = pop(@_);
#  print @_."\n";
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
    
    if($pos1 == 1) {
     $fin = @pos;
    }
    else {
     # on traduit la pos relative en pos absolue et on va chercher l'indice auquel ca correspond.
     my $posRec = $pos1 * $bpTot;
     # on cherche la position de recombinaison.
     $fin=$deb;
     #la recombinaison peut arriver entre le dernier SNP et la fin de la sequence.
     while($pos[$fin] - $posRec < 0 ){
       $fin++;
       if($fin == (@pos-1)){
         if($pos[$fin] - $posRec < 0){ $fin++; }
         last;
       }
     } #$fin < (@pos-1) & 
#    if($sim == 735 & $proID == "408366"){ print $deb." - ".$fin." ; ".$pos1." - ".$posRec." ; ".$pos[$fin-1]." de ".$pos[@pos-1]."(".@pos.")\n"; }
    }
    $seqTester = $seqTester + length(substr($seqAnc{$ind1}, $deb, ($fin-$deb) ));
    print OUT substr($seqAnc{$ind1}, $deb, ($fin-$deb) );
    $deb=$fin;
   }
   if($seqTester != @pos){ print "Probleme: $seqTester vs ".@pos." @hap\n"; exit}
   print OUT "\n";
  }
}



