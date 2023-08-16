
#!/usr/bin/perl -w

use strict;
use Carp;
use Getopt::Long;

my($file, $outfile, $line, $taux_diff, $tot1, $M1, $NM1, $NC1, $U1, $nb_pat1, $tot2, $M2, $NM2, $U2, $NC2, $nb_pat2, $seuil_min_couv, $contexte, $sonde, $key, $seuil_diff_taux_min);
my(@ligne);
my(%good_hash);

$seuil_min_couv=10;
$seuil_diff_taux_min=0.20;

GetOptions(
"f=s" => \$file

);

$outfile="$file.best.tsv";
open(FILE, "$file") or die "Can't open $file !\n";
open(OUTFILE, ">$outfile") or die "Can't open $outfile !\n";

while ($line = <FILE>)
{
	chomp($line);
	if($line !~ m/Sonde/)
	{
		@ligne=split("\t", $line);
		$sonde=$ligne[0];
		$contexte=$ligne[2];
		$taux_diff=$ligne[4];
		$tot1=$ligne[5];
		$M1=$ligne[8];
		$NM1=$ligne[9];
		$U1=$ligne[10];
		$NC1=$ligne[11];
		$nb_pat1=$ligne[12];

		$tot2=$ligne[13];
		$M2=$ligne[16];
		$NM2=$ligne[17];
		$U2=$ligne[18];
		$NC2=$ligne[19];
		$nb_pat2=$ligne[20];
		#print "$taux_diff\t$tot1\t$M1\t$NM1\t$U1\t$NC1\t$nb_pat1\t$tot2\t$M2\t$NM2\t$U2\t$NC2\t$nb_pat2\n";

		#si la couverture moyenne est inférieur au seuil d'un côté ou de l'autre
		if( $tot1/$nb_pat1 > $seuil_min_couv and $tot2/$nb_pat2 > $seuil_min_couv )
		{
			#si tout n'est pas NM et tout n'est pas M
			if($NM1+$NM2 != $nb_pat1+$nb_pat2 and $M1+$M2 != $nb_pat1+$nb_pat2)
			{
				#si le nombre de non couvert est inférieur à la moitié des patients
				if($NC1 < $nb_pat1/2 and $NC2 < $nb_pat2/2)
				{
					#si la différence de taux est supérieur au seuil d'interet
					if($taux_diff gt $seuil_diff_taux_min)
					{
						#s'il y a plus de la moitié de U,on définit le contexte
						if($U1+$U2 > ($nb_pat1+$nb_pat2)/2)
						{
							$ligne[2]="U50";
						}
						else
						{
							$ligne[2]="D50";
						}
						$good_hash{$sonde}=[@ligne];
					}
					 
				}
			}	
		}
		
	}
	else
	{
		print OUTFILE "$line\n";
	}
}

for $key ( sort { $good_hash{$b}[4] <=> $good_hash{$a}[4] } keys %good_hash )
{
	print OUTFILE join("\t", @{ $good_hash{$key}}, "\n");
}



