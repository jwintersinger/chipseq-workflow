#!/usr/bin/perl

=head2 overlappercent_outputsimple_sortedonly.pl

PROGRAM: overlappercent_outputsimple_sortedonly.pl
DESIGNED TO BE MUCH FASTER!
GOAL: calculate the overlap in two gff files, allow near overlaps... and output the differences into two files...
INPUT: two SORTED AND COMMENTS KILLED GFF only filenames then a slopvalue (set to zero for strict overlaps)
ALL must be gff format and SORTED AND COMMENTS KILLED!

OUTPUT:
(1) two lines of text quantifying overlaps for each
(2) the files diffs_firstfile, diffs_secondfile
 (3) the files same_firstfile, same_secondfile; these files have the ones in common


EXAMPLE:
overlappercent_outputsimple_sortedonly.pl  E2F1_MARKMAX_H_T05P05_S50_CHR_G2.gff TTTSSCGC_CHR.gff 0


EXAMPLE OF OUTPUT (from above example):
E2F1_MARKMAX_H_T05P05_S50_CHR_G2.gff length is 583 num_overlapped: 36 percent: 6.17495711835334
TTTSSCGC_CHR.gff length is 511 num_overlapped: 39 percent: 7.6320939334638

NOTES:
(0) note that testing notes apply to previous diffs version; this one is untested
(1) some  testing, verification of output...
(2) will check for the first field being the same - will only match chr17 to chr17 , for example
(3) runs pretty fast - above example ran in just a few seconds (may be 3 seconds or so)

=cut


sub isoverlappingslop
{

	my $extrabit=$_[0]; #allows for sloppiness - if this is set to zero, then there is no sloppiness
	my $acor=$_[1]-$extrabit;
	my $bcor=$_[2]+$extrabit;
	my $ccor=$_[3]-$extrabit;
	my $dcor=$_[4]+$extrabit;
	
	if ($acor>=$ccor && $acor<=$dcor)
	{return 1;}
	if ($bcor>=$ccor && $bcor<=$dcor)
	{return 1;}
	if ($ccor>=$acor && $ccor<=$bcor)
	{return 1;}
	if ($dcor>=$acor && $dcor<=$bcor)
	{return 1;}
	return 0;
}


$fname1=$ARGV[0]; #input filenam
$fname2=$ARGV[1]; #second filename
$slopval=$ARGV[2]; #extra bit to add to each piece to allow near overlaps

$fname1_out="diffs_firstfile";
$fname2_out="diffs_secondfile";
$fname1_same="same_firstfile";
$fname2_same="same_secondfile";

open (DOG,"<$fname1");
@first=<DOG>;
close(DOG);



open (DOG2,"<$fname2");
@second=<DOG2>;
close(DOG2);

$first_numlines=$#first +1;
$second_numlines=$#second+1;

#print "first numlines = $first_numlines second numlines = $second_numlines\n";

for ($i=0; $i<$first_numlines; $i++)
{
	$first_overlaps[$i]=0;
}

for ($i=0; $i<$second_numlines; $i++)
{
	$second_overlaps[$i]=0;
}

$posinsecond=0;
$inchr=0;
for ($i=0; $i < $first_numlines; $i++)
{
	$inchr=0;
	@pieces = split "\t", $first[$i];
	$fstart=$pieces[3];
	$fend=$pieces[4];
	#MB note: PRINT BELOW!
	#print "posinsecond=$posinsecond\n";
	MYSTEP: for ($j=$posinsecond; $j < $second_numlines; $j++)
	{
		
		@spieces = split "\t", $second[$j];
		$sstart=$spieces[3];
		$send=$spieces[4];
		#print "fstart is $fstart fend is $fend sstart is $sstart send is $send\n";
		#!!!!!!MUST FIX NEXT STATEMENT TO BREAK LOOP!!!!!
		if (($pieces[0] ne $spieces[0]) and ($inchr==1)) {$inchr=0; last MYSTEP;}
		if (($pieces[0] eq $spieces[0])) 
			{ if ($inchr==0) {$posinsecond=$j;$inchr=1};
				if (isoverlappingslop($slopval, $fstart,$fend,$sstart,$send))
				{
					$first_overlaps[$i]=$first_overlaps[$i] +1;
					$second_overlaps[$j] = $second_overlaps[$j] + 1;
					$posinsecond=$j;
					
			#print $first[$i] . $second[$j] . "******\n";
				}
		#print "second overlaps 0 is $second_overlaps[0]\n";
		#print "i is $i j is $j\n";
		
			}
		else
			{$inchr=0;}
	}

}
#print "second overlaps[9] is " .  $second_overlaps[9] . "\n";
$fcount=0;
$scount=0;

open(OUTSAME,">$fname1_same");
open (OUTFIRST,">$fname1_out");
for ($i=0; $i < $first_numlines; $i++)
{
	if ($first_overlaps[$i]>0)
	{
		print OUTSAME $first[$i];
		$fcount++;
	}
	else
	{
		print OUTFIRST $first[$i];
	}
}
close(OUTFIRST);
close(OUTSAME);

$fper=100*$fcount/$first_numlines;

open(OUTSAME,">$fname2_same");
open (OUTSECOND ,">$fname2_out");
for ($i=0; $i < $second_numlines; $i++)
{
	if ($second_overlaps[$i]>0)
	{
		print OUTSAME $second[$i];
		$scount++;
	}
	else
	{
		print OUTSECOND $second[$i];
	}
}

close(OUTSECOND);
close(OUTSAME);

$sper=100*$scount/$second_numlines;

print "$fname1_out was created and $fname2_out was created\n";
print "overlapped: $fname1 $fcount of $first_numlines were overlapped for $fper %\n";
print "overlapped: $fname2 $scount of $second_numlines were overlapped for $sper %\n";
#print "$fname2 length is $second_numlines num_overlapped: $scount percent: $sper\n";

