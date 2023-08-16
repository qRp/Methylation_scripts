
#!/usr/bin/perl -w

use strict;
use Carp;
use Getopt::Long;

my($line, $file1, $file2, $chr, $pos, $sonde, $score, $nb_M, $nb_NM, $nb_U, $nb_NC, $ScoreM, $ScoreNM, $ScoreU, $ScoreNC, $ScoreDM, $ScoreDNM, $ScoreDU, $ScoreDNC, $cpt_non_matching, $outfile, $ligne2, $force, $contexte); 
my(@ligne);
my(%comp_hash);

$ScoreM=10;
$ScoreNM=10;
$ScoreU=3;
$ScoreNC=1;
$ScoreDM=-10;
$ScoreDNM=-10;
$ScoreDU=-3;
$ScoreDNC=-1;

$cpt_non_matching=0;
$contexte="UNK";

GetOptions(
"f2=s" => \$file2,

"f1=s" => \$file1,

"o=s" => \$outfile,

"f" => \$force,

"c=s" => \$contexte

);


if ($contexte ne "CpG" and $contexte ne "CHG" and $contexte ne "CHH")
{
	$contexte="UNK";
}


open(FILE1, "$file1") or die "Can't open $file1 !\n";
open(FILE2, "$file2") or die "Can't open $file2 !\n";
open(OUTFILE, ">$outfile") or die "Can't open $outfile !\n";
print OUTFILE "Sonde\tRegion\tContexte\tSimilarity_Score\tsum_methyl.$file1\tsum_total.$file1\ttaux.glob.$file1\tnb_M.$file1\tnb_NM.$file1\tnb_U.$file1\tnb_NC.$file1\tnb_all.$file1\tum_methyl.$file2\tsum_total.$file2\ttaux.glob.$file2\tnb_M.$file2\tnb_NM.$file2\tnb_U.$file2\tnb_NC.$file2\tnb_all.$file2\n";

while ($line = <FILE1>)
{
	chomp($line);
	if($line !~ m/chr/)
	{
		@ligne=split("\t", $line);
		$chr=$ligne[0];
		$pos=$ligne[1];
		$sonde=$chr."_".$pos;
		$comp_hash{$sonde}{"region"}=$ligne[2];
		$comp_hash{$sonde}{"nb_M"}=$ligne[6];
		$comp_hash{$sonde}{"nb_NM"}=$ligne[7];
		$comp_hash{$sonde}{"nb_U"}=$ligne[8];
		$comp_hash{$sonde}{"nb_NC"}=$ligne[9];
		$comp_hash{$sonde}{"nb_all"}=$ligne[10];
		$comp_hash{$sonde}{"line"}=$ligne[3]."\t".$ligne[4]."\t".$ligne[5]."\t".$ligne[6]."\t".$ligne[7]."\t".$ligne[8]."\t".$ligne[9]."\t".$ligne[10];
	}
}


while ($line = <FILE2>)
{
	chomp($line);
	if($line !~ m/chr/)
	{
		@ligne=split("\t", $line);
		$chr=$ligne[0];
		$pos=$ligne[1];
		$sonde=$chr."_".$pos;
		if(exists($comp_hash{$sonde}) and ($comp_hash{$sonde}{"nb_all"} eq $ligne[10] or $force eq 1))
		{
			$score=0;
			$nb_M=$ligne[6];
			$nb_NM=$ligne[7];
			$nb_U=$ligne[8];
			$nb_NC=$ligne[9];
			$ligne2=$ligne[3]."\t".$ligne[4]."\t".$ligne[5]."\t".$ligne[6]."\t".$ligne[7]."\t".$ligne[8]."\t".$ligne[9]."\t".$ligne[10];
			$score=$score+&same_score($nb_M, $comp_hash{$sonde}{"nb_M"}, $ScoreM);
			$score=$score+&same_score($nb_NM, $comp_hash{$sonde}{"nb_NM"}, $ScoreNM);
			$score=$score+&same_score($nb_U, $comp_hash{$sonde}{"nb_U"}, $ScoreU);
			$score=$score+&same_score($nb_NC, $comp_hash{$sonde}{"nb_NC"}, $ScoreNC);
			$score=$score+($ScoreDM*abs($nb_M-$comp_hash{$sonde}{"nb_M"}));
			$score=$score+($ScoreDNM*abs($nb_NM-$comp_hash{$sonde}{"nb_NM"}));
			$score=$score+($ScoreDU*abs($nb_U-$comp_hash{$sonde}{"nb_U"}));
			$score=$score+($ScoreDNC*abs($nb_NC-$comp_hash{$sonde}{"nb_NC"}));
			print OUTFILE "$sonde\t$comp_hash{$sonde}{region}\t$contexte\t$score\t$comp_hash{$sonde}{line}\t$ligne2\n";

		}	
		else
		{
			$cpt_non_matching=$cpt_non_matching+1;
		}
	}
}

sub same_score {
	my($val1, $val2, $factor);
	$val1=$_[0];
	$val2=$_[1];
	$factor=$_[2];
	if($val1 <= $val2)
	{
		return $val1*$factor;
	}
	else
	{
		return $val2*$factor;
	}
}

print "Total non matching : $cpt_non_matching\n";

