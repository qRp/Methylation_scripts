
#!/usr/bin/perl -w

use strict;
use Carp;
use Getopt::Long;

my($line, $chr, $pos, $percentage, $number_of_line, $cpt, $chemin, $ID, $statut, $read, $verbose, $outfile, $element, $methyl, $non_methyl, $total, $percent_methyl, $data_point, $id_patient, $id_list, $patient, $ligne_to_print, $total_ligne, $cpt_patient, $seuil, $id, $position, $data_file, $ext, $id1_methyl, $id1_total, $id1_taux, $id2_methyl, $id2_total, $id2_taux, $id1, $id2, $id1_statut, $id2_statut, $seuil_couv, $seuil_haut, $seuil_bas, $bed_file, $region, $stop, $start, $nearest_region, $nearest_start, $header, $temp, $key, $offset, $methyl_index, $total_index, $taux_index, $taux_global, $mode );
my(@ligne, @listid, @args);
my(%methyl_hash, %ref_hash, %patients, %bed_hash, %ref_inv_hash, %all_methyl);

%methyl_hash=();
%bed_hash=();
%ref_hash=();
$cpt_patient=0;
$seuil_couv=10;
$seuil_bas=0.20;
$seuil_haut=0.80;

GetOptions(
"f=s" => \$data_file,

"b=s" => \$bed_file,

"m=s" => \$mode,

"v" =>\$verbose

);

if($mode ne "S" && $mode ne "C")
{
	print "le mode doit être défini sur 'S' ou 'C' !";
	die;
}


sub load_methyl {
	open(METHYL, "$data_file") or die "Can't open Methyl file : $data_file !";
	if($verbose) #Affiche le pourcentage, mais coute jusqu'à 1 minute initiale + temps affichage et surcharge le nohup"
	{
		$cpt=0;
		$number_of_line = `wc -l < $chemin`;
		print "Nombre de ligne : $number_of_line\n";
		chomp($number_of_line);
	}
	while ($line = <METHYL>)
	{
		chomp($line);
		if($verbose)
		{
			$cpt=$cpt+1;
			if($cpt%10000 == 1 )
			{
				$percentage=($cpt/$number_of_line)*100;
				$percentage=int($percentage);
				print "\r $percentage pourcent de $chemin chargés";
			}
		}
		if($line =~ "Position")
		{
			@ligne=split("\t",$line);
			$position=$ligne[0];
			for (my $i=1; $i<=$#ligne ; $i++)
			{
				$id=$ligne[$i];
				$ext=$ligne[$i];
				$id=~s/\..*//g;
				$ext=~s/.*\.//g;
#				print "$id $ext \n" ; 
				$ref_hash{$i}=();
				$ref_hash{$i}{"ext"}=$ext; 
				$ref_hash{$i}{"id"}=$id;
			}
		}
		else
		{
			@ligne=split("\t",$line);
			$position=$ligne[0];
			$chr=$position;
			$chr=~s/chr//;
			$chr=~s/_.*//g;
			$position=~s/.*_//g;
#			print "$chr $position \n";
			if (! exists($methyl_hash{$chr}))
			{
				$methyl_hash{$chr}={};
			}
			$methyl_hash{$chr}{$position}={};
			for (my $i=1; $i<=$#ligne; $i++)
			{
				$methyl_hash{$chr}{$position}{$i}=$ligne[$i];
			}
			
		}
	}
	open(BED, "$bed_file") or die "Can't open bed file : $bed_file !";
	while($line=<BED>)
	{
		@ligne=split("\t",$line);
		$region=$ligne[0];
		$chr=$ligne[1];
		$start=$ligne[2];
		$stop=$ligne[3];
		if (! exists($bed_hash{$chr}))
		{
			$bed_hash{$chr}={};
		}
		$bed_hash{$chr}{$region}={start=>$start, stop=>$stop};
	}
	
}

sub appariement{
	#on recrée un hash dans le bon sens
	%ref_inv_hash={};
	for $key (keys %ref_hash)
	{
		if(!exists($ref_inv_hash{$ref_hash{$key}{"id"}}))
		{
			$ref_inv_hash{$ref_hash{"id"}}=();
		}
		$ref_inv_hash{$ref_hash{$key}{"id"}}{$ref_hash{$key}{"ext"}}=$key;
#		print "$key $ref_hash{$key}{'id'} $ref_hash{$key}{'ext'}\n";
	}	
}

sub find_region {
	@args=@_;
	$chr=$args[0];
	$position=$args[1];
	$nearest_region="";
	$nearest_start=1000000;
	if(exists($bed_hash{$chr}))
	{
		for $region (keys %{ $bed_hash{$chr}})
		{
			$start=$bed_hash{$chr}{$region}{"start"};
			$stop=$bed_hash{$chr}{$region}{"stop"};
			if($start > $stop)
			{
				$temp=$stop;
				$stop=$start;
				$start=$temp;
				print "$region inversée : $stop $start\n";
			}
			if($position >= $start and $position <=$stop)
			{
				return $region;
			}
			elsif( abs($position-$start) < abs($position-$nearest_start))
			{
				$nearest_start=$start;
				$nearest_region=$region.".isclosest";
			}
		}
		return $nearest_region;
	}
	else
	{
		return "chr_not_in_bed";
	}

}

	

sub comp_methyl_all {
	if($mode eq "S")
	{
		open(OUTFILE, ">Glob_S.methyl_stat.tsv");
	}
	elsif($mode eq "C")
	{
		open(OUTFILE, ">Glob_C.methyl_stat.tsv");
	}

	$header="chr\tposition\tregion\tsomme.methyl\tsomme.total\ttaux.global\tnb_M\tnb_NM\tnb_U\tnb_NC\ttotal\n";
	print OUTFILE "$header";
	
	for $chr (keys %methyl_hash)
	{
		for $position (keys %{ $methyl_hash{$chr}})
		{
			for $id (keys %ref_inv_hash)
			{
			if($id =~m/[0-9][0-9][0-9]$mode/)
			{
			#calcul si M, NM, NC ou U (TODO adapter id1)
			#si moins de 20 pourcent methyl, statu NM, si plus de 80, statu M, si entre les deux ou si pas couvert, statut U
			$methyl_index=$ref_inv_hash{$id}{"nb_methyl"};
			$total_index=$ref_inv_hash{$id}{"total_couv"};
			$taux_index=$ref_inv_hash{$id}{"pourcentage_methyl"};

			#initiaisation
			if(!exists($all_methyl{$chr}{$position}))
			{
				$all_methyl{$chr}{$position}{"cpt_all"}=0;
				$all_methyl{$chr}{$position}{"cpt_M"}=0;
				$all_methyl{$chr}{$position}{"cpt_NM"}=0;
				$all_methyl{$chr}{$position}{"cpt_U"}=0;
				$all_methyl{$chr}{$position}{"cpt_NC"}=0;
				$all_methyl{$chr}{$position}{"methyl_somme"}=0;
				$all_methyl{$chr}{$position}{"couv_somme"}=0;
			}

			if($methyl_hash{$chr}{$position}{$taux_index} < $seuil_bas){$statut="NM";}
			elsif($methyl_hash{$chr}{$position}{$taux_index} > $seuil_haut){$statut="M";}
			else{$statut="U";};
			if($methyl_hash{$chr}{$position}{$total_index} < $seuil_couv){$statut="NC";} 
			
			$all_methyl{$chr}{$position}{"cpt_all"}=$all_methyl{$chr}{$position}{"cpt_all"}+1;
			$all_methyl{$chr}{$position}{"methyl_somme"}=$all_methyl{$chr}{$position}{"methyl_somme"}+$methyl_hash{$chr}{$position}{$methyl_index};
			$all_methyl{$chr}{$position}{"couv_somme"}=$all_methyl{$chr}{$position}{"couv_somme"}+$methyl_hash{$chr}{$position}{$total_index};
			if($statut eq "M"){$all_methyl{$chr}{$position}{"cpt_M"}=$all_methyl{$chr}{$position}{"cpt_M"}+1;}
			if($statut eq "NM"){$all_methyl{$chr}{$position}{"cpt_NM"}=$all_methyl{$chr}{$position}{"cpt_NM"}+1;}
			if($statut eq "U"){$all_methyl{$chr}{$position}{"cpt_U"}=$all_methyl{$chr}{$position}{"cpt_U"}+1;}
			if($statut eq "NC"){$all_methyl{$chr}{$position}{"cpt_NC"}=$all_methyl{$chr}{$position}{"cpt_NC"}+1;}
			
			}#if S
			}#boucle id
		}#boucle pos
		print "fin du chromosome $chr\n";
	}#boucle chr
	print "Fin du Load, écriture des résultats\n";


	for $chr (keys %all_methyl)
	{
		for $position (keys %{ $all_methyl{$chr}})
		{
			#recherche dans le bed de la région correspondante
			$region=&find_region($chr, $position);
	$header="chr\tposition\tregion\tsomme.methyl\tsomme.total\ttaux.global\tnb_methyl\tnb_non_methyl\tnb_U\tnb_NC\n";
			$taux_global=$all_methyl{$chr}{$position}{"methyl_somme"}/$all_methyl{$chr}{$position}{"couv_somme"};
			print OUTFILE "$chr\t$position\t$region\t$all_methyl{$chr}{$position}{couv_somme}\t$all_methyl{$chr}{$position}{methyl_somme}\t$taux_global\t$all_methyl{$chr}{$position}{cpt_M}\t$all_methyl{$chr}{$position}{cpt_NM}\t$all_methyl{$chr}{$position}{cpt_U}\t$all_methyl{$chr}{$position}{cpt_NC}\t$all_methyl{$chr}{$position}{cpt_all}\n"; 
		}
	}
}



&load_methyl();
&appariement();
$offset=$offset*6;
&comp_methyl_all();
#qprint "$methyl_hash{7}{16400206}{6} $ref_hash{6}{'ext'} $ref_hash{6}{'id'}\n";

