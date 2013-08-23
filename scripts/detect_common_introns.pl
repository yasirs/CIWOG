#!/usr/bin/perl

#############################
# Title:
# Description:
# Date: Original Inception.  
# Current Version:
# License: 
#############################

require 'conf.pl';
my $cid = shift; #alignment and intron file prefix.

###
# parameters
###
	my $maxSlide = shift; # number of amino acids that can be associated with an intron's position
	my $maxGap = shift; # number of gaps that can be claimed for an intron's position
	my $flank = shift; # simRegion length of flanking regions of introns that are assesed for sequence similarity.  
	my $minSeqSim = shift; # minimum percentage of identical amino acids in flanking and intronSeg regions required for two sequences to be similar
	my $minAbsentSim = shift;
	my $freeGaps = shift; # number of gaps that are allowed in assessing sequence similarity

	#not used
	my $repName = shift;
	$repName = ($repName ne $cVer) ? $repName : "";

print "# $cid: maxSlide $maxSlide, maxGap $maxGap, flank $flank, minSeqSim $minSeqSim, freeGaps $freeGaps\n";

my $clusterFile = "${INPUT_DIR}$cid/$cid.aln";
my $intronsFile = "${OUTPUT_DIR}$cid/introns.csv";

$pv = 1; #print verbose
my %g;
my @introns = ();
my %bounds = ();
my %membersWintron;

###
# Main
###

    loadAln($clusterFile);

    fetchIntrons($cid);
    setIntronBoundaries();
    my $g = groupIntrons();
    if ($pv){print "#add absent\n";}
    $g = addAbsentIntrons($g);
    
    printSQL($g,$cid);

###
# End Main
###


###
# Functions
###


sub printSQL{
  # tables: Clusters, ClusterIntron2Gene, ClusterIntrons
  my ($g,$cid) = @_;
  my $icnt=1;

  #clear tables
  my $mc=0;
  foreach my $k (keys(%aln)){
      foreach my $j (keys (%{$aln{$k}})){
	  $mc++;
      }
  }

  print "insert into Clusters$cVer$repName (clusterUID,distinctIntronCount,geneCount) VALUES ('$cid',$#$g+1,$mc);\n";

  $g = sortIntronPos($g);

  foreach my $i (@$g){
	#print "group $i->[0] to $i->[1]\n";

	# ClusterIntrons
        print "#ClusterIntron UID $icnt\n";
	print "insert into ClusterIntrons$cVer$repName (clusterUID,clusterIntronUID,alnPosStart,alnPosEnd) VALUES ('$cid',$icnt,'$i->[0]','$i->[1]');\n";


	foreach my $j (@{$i}[2..$#$i]){
	    # ClusterIntron2Gene
	    my ($org,$annotationUID,$intronUID,$intronNum,$proteinPos,$phase,$spliceType,$left,$alnPos,$right) = @{$introns[$j]};
	    print "insert into ClusterIntron2Gene$cVer$repName (clusterUID,clusterIntronUID,org,annotationUID,intronUID,intronNum,intronType,phase,alnPos,intronSeq,start,stop,segment,proteinPos,strand) VALUES";
	    if ($spliceType eq "absent"){
		print "('$cid',$icnt,'$org','$annotationUID',NULL,NULL,'$spliceType',NULL,'$i->[0]','$intronDat{$j}[0]','$intronDat{$j}[1]','$intronDat{$j}[2]','$intronDat{$j}[3]',NULL,'$intronDat{$j}[4]');\n";
	    }else{
		print "('$cid',$icnt,'$org','$annotationUID',$intronUID,$intronNum,'$spliceType','$phase','$alnPos','$intronDat{$j}[0]','$intronDat{$j}[1]','$intronDat{$j}[2]','$intronDat{$j}[3]','$proteinPos','$intronDat{$j}[4]');\n";
	    }
	}
	$icnt++;
  }

  return;
}

sub printIntronGroups{
  my $g = shift;
  foreach my $i (sort sort_lr @$g){
	print "#group $i->[0] to $i->[1]\n";
	foreach my $j (@{$i}[2..$#$i]){
		print "#\t".join(",",@{$introns[$j]})."\n";
	}
  }
}


sub addAbsentIntrons{
    my $g = shift;
    # for each sequence,for each intron group, check boundaries, check membership in group, check if similar to any of intron-containing sequences, add absent, im, em.
      foreach my $m (keys(%aln)){
	  foreach my $n (keys (%{$aln{$m}})){
	      INTL: for(my $i=0;$i<scalar(@$g);$i++){
		  # skip utr
		  if ($g->[$i]->[0] =~ /utr/){
		      next;
		  }
		  # determine if similar to other sequences with intron
		  my $haveSim=0;
          for (my $j=2;$j<scalar(@{$g->[$i]});$j++){
		      my $t = $introns[$g->[$i]->[$j]];
              if ($m eq $t->[0] && $n eq $t->[1]){
                         # sequence is at this position, skip
                         if($pv){print "# at $g->[$i]->[1] position $m $n $t->[6]\n";}
                          next INTL;
              }
          }
		 # if intron is not within boundaries of gene, skip
                  if (!($bounds{$m}{$n}{left} <= $g->[$i]->[1] && $bounds{$m}{$n}{right} >= $g->[$i]->[0])){
                       $introns[++$#introns] = [$m,$n,NULL,NULL,NULL,NULL,'em'];
                       push @{$g->[$i]},$#introns;
                      next INTL;
                  }
		  my $score;
		  for (my $j=2;$j<scalar(@{$g->[$i]});$j++){
		      my $t = $introns[$g->[$i]->[$j]];
		      if ($t->[6] eq 'im' or $t->[6] eq 'em'){
				  # internal and external missing states cannot link more sequences into the group
				  next;
		      }
		      my ($intronSeg,$intronSegsize,$flankL,$flankR,$flankLsize,$flankRsize,$gapCount,$simSeg,$simL,$simR,$sGap) = cFlanks($t->[0],$t->[1],$m,$n,$g->[$i]->[0],$g->[$i]->[1]);
		      my $simCond=0;


                      $score = ($intronSegSize+$flankLsize+$flankRsize > 0) ? ($intronSeg+$flankL+$flankR)/($intronSegSize+$flankLsize+$flankRsize) : 0;
		      #if ($score >= $minAbsentSim){ #same score as introns
			#$simCond=2;
                      #}

		      if ($flankL > 0 && $flankL/$flankLsize >= $minAbsentSim){
		      			$simCond++;
		       }
		      if ($flankR > 0 && $flankR/$flankRsize >= $minAbsentSim){
		      			$simCond++;
		      }


		      if ($simCond==2 && ($sGap==0) && $introns[ $g->[$i]->[2] ]->[7] != -2 ){	
		      #if ($simCond==2 && $introns[ $g->[$i]->[2] ]->[7] != -2){
			  	$best = $j;
			  	$haveSim=1;
			  	last;
		      }elsif($gapCount > ((1-$minAbsentSim)*($flankLsize+$flankRsize+$intronSegsize)) ){
			  	# if region is full of gaps more than the minSeq requirement, it is not able to match to any other seq
			  	if ($pv){print "#all gap\n";}
			  	$introns[++$#introns] = [$m,$n,NULL,NULL,NULL,NULL,'im'];
			  	push @{$g->[$i]},$#introns;
			  	next INTL;  ## slight speed improvement when freeGaps or flank is large
		      }else{
			  	if($pv){print "#sim too low $t->[0],$t->[1],$m,$n : $score\n";}
		      }
		  }
		  if ($haveSim == 1){
 		       if ($pv){print "#absent add ($score);\n";}
		       $introns[++$#introns] = [$m,$n,NULL,NULL,NULL,NULL,'absent'];
		       push @{$g->[$i]},$#introns;
		  }else{
		       $introns[++$#introns] = [$m,$n,NULL,NULL,NULL,NULL,'im'];
			   push @{$g->[$i]},$#introns;
		  }
	      }
	  } # end gene
      } # end org
	return $g;
}


sub groupIntrons{
    # groups introns into single linkage clusters based on leftAlnPos,rightAlnPos,and sequence similarity
    # additional rule: one intron per gene per cluster
    # group  list
    my $g = initGroups();
    $g = makeGroups($g);
    return $g;
}


sub initGroups{
    my $g = [];
    for (my $i=0;$i<scalar(@introns);$i++){
	$g->[$i] = [ $introns[$i]->[7], $introns[$i]->[9], $i] ;
    }
    return $g;
}

sub makeGroups{
  my ($g) = @_;
  for (my $e1=0;$e1<scalar(@$g);$e1++){
        if ($g->[$e1]->[0] ne "-1"){
          #test print "search\n";
          $g = groupIdent($g,$e1);
        }
  }
  for (my $e1=0;$e1<scalar(@$g);$e1++){
	if ($g->[$e1]->[0] ne "-1" ){
	  $g = groupOne($g,$e1);
	}
  }
  my $g2 = [];
  for (my $e1=0;$e1<scalar(@$g);$e1++){

	# ignore utr groups
	if ($g->[$e1]->[0] =~ /utr/){
		$g2->[++$#$g2] = $g->[$e1];
		next;
	}
	
	if ($g->[$e1]->[0] != -1){ # remove empty groups
		my ($l,$r) = (99999,-99999);
		# set l and r to l and r of intron group without slide/gap
		for (my $i=2;$i<scalar(@{$g->[$e1]});$i++){
			$l = min($l,$introns[$g->[$e1]->[$i]]->[8]);
			$r = max($r,$introns[$g->[$e1]->[$i]]->[8]);
		}
		$g2->[++$#$g2] = [$l,$r, @{$g->[$e1]}[2..$#{$g->[$e1]}] ];
	}
  }

  return $g2;
}



sub groupIdent{
    my ($g,$one) = @_;
    my ($best,$overlap) = (-9999,0);

    for (my $e1=0;$e1<scalar(@$g);$e1++){
       if ($e1 == $one){ next;} # if records are the same
        if ($g->[$one]->[0] == -1){next;} # skip empty groups

        # utr introns
        if ($g->[$one]->[0] eq $g->[$e1]->[0] && $g->[$one]->[0] =~ /utr/){
                @{$g->[$e1]} = (@{$g->[$one]}[0..1] , @{$g->[$e1]}[2..$#{$g->[$e1]}], @{$g->[$one]}[2..$#{$g->[$one]}] );
                $g->[$one] = [-1,-1,-1];
                return $g;
        }        
        if ($g->[$one]->[0] =~ /utr/ || $g->[$e1]->[0] =~ /utr/){ next;} #skip utr.

        	# protein coding introns
            my $l = min($g->[$one]->[0],$g->[$e1]->[0]);
            my $r = max($g->[$e1]->[1],$g->[$one]->[1]);
            if (
		            $introns[ $g->[$one]->[2] ]->[8] == $introns[$g->[$e1]->[2] ]->[8] #identical position
					&&
					compareSeq($g->[$e1], $g->[$one], $l, $r) #sequence matching condition
					&&
					intronSetOverlap($g->[$e1],$g->[$one]) # one intron per group per sequence
                   ){
                        if ($pv){ print "***** $introns[ $g->[$one]->[2] ]->[8],$introns[$g->[$e1]->[2] ]->[8],$introns[$g->[$e1]->[2] ]->[2],$introns[$g->[$one]->[2] ]->[2]\n"; }
						$best = $e1;
                }		
     }
    if ($best != -9999){
               my $l = min($g->[$one]->[0],$g->[$best]->[0]);
                my $r = max($g->[$best]->[1],$g->[$one]->[1]);
                @{$g->[$best]} = ($l,$r , @{$g->[$best]}[2..$#{$g->[$best]}], @{$g->[$one]}[2..$#{$g->[$one]}] );
                $g->[$one] = [-1,-1,-1];
    }
    return $g;
}



sub groupOne{
    my ($g,$one) = @_;
    my ($best,$overlap) = (-9999,0);
    for (my $e1=0;$e1<scalar(@$g);$e1++){
       if ($e1 == $one){ next;} # if records are the same

	if ($g->[$one]->[0] == -1){next;} # skip empty groups

	# utr introns
	if ($g->[$one]->[0] eq $g->[$e1]->[0] && $g->[$one]->[0] =~ /utr/){
		@{$g->[$e1]} = (@{$g->[$one]}[0..1] , @{$g->[$e1]}[2..$#{$g->[$e1]}], @{$g->[$one]}[2..$#{$g->[$one]}] );
		$g->[$one] = [-1,-1,-1];
		return $g;
	}
	if ($g->[$one]->[0] =~ /utr/ || $g->[$e1]->[0] =~ /utr/){ next;} #remove utr from consideration

	# protein coding introns
	if ($g->[$one]->[0] <= $g->[$e1]->[1] && $g->[$one]->[1] >= $g->[$e1]->[0]  ){
	    my $l = min($g->[$one]->[0],$g->[$e1]->[0]);
	    my $r = max($g->[$e1]->[1],$g->[$one]->[1]);
	    if( 
		compareSeq($g->[$e1], $g->[$one], $l, $r)
		&&
		intronSetOverlap($g->[$e1],$g->[$one]) ){

		my $tmp_overlap = 0;
		for (my $i=$l;$i<=$r;$i++){
		    if ($i <= $g->[$e1]->[1] && $i <= $g->[$one]->[1] && $i >= $g->[$e1]->[0] && $i >= $g->[$one]->[0]){
			$tmp_overlap++;
		    }
		}
		$best = ($tmp_overlap > $best) ? $e1 : $best;

	    }
	}
    }
    if ($best != -9999){
	        my $l = min($g->[$one]->[0],$g->[$best]->[0]);
                my $r = max($g->[$best]->[1],$g->[$one]->[1]);
		@{$g->[$best]} = ($l,$r , @{$g->[$best]}[2..$#{$g->[$best]}], @{$g->[$one]}[2..$#{$g->[$one]}] );
		$g->[$one] = [-1,-1,-1];
    }

    return $g;
}

sub intronSetOverlap{
    # this determines if there any common between the two sets
    my ($a,$b) = @_;
    foreach my $i (@{$a}[2..$#$a]){
	foreach my $j (@{$b}[2..$#$b]){
	    if ($pv){print "#check\n";}
	    if ($introns[$i]->[0] eq $introns[$j]->[0] && $introns[$i]->[1] eq $introns[$j]->[1]){
		print "$introns[$i]->[0] eq $introns[$j]->[0] && $introns[$i]->[1] == $introns[$j]->[1] \n";
		return 0;
	    }
	}
    }
    return 1;
}

sub compareSeq{
    # this compares one set of sequences with a different set and determines if there is at least one pair of sequences between the two sets that has minimum sequence similarity
    my ($a,$b,$l,$r) = @_;
    foreach my $i (@$a[2..$#$a]){
	foreach my $j (@$b[2..$#$b]){
	    #test print "data: @{$introns[$i]}[0..1],@{$introns[$j]}[0..1],$l,$r\n";
	    my ($intronSeg,$intronSegsize,$flankL,$flankR,$flankLsize,$flankRsize,$gapCount,$simCount) = cFlanks(@{$introns[$i]}[0..1],@{$introns[$j]}[0..1],$l,$r);
	    my $score = ($intronSegSize+$flankLsize+$flankRsize > 0) ? ($intronSeg+$flankL+$flankR)/($intronSegSize+$flankLsize+$flankRsize) : 0;

	    if ($score >= $minSeqSim){
		return 1;
	    }

	}
    }
    return 1;
    return 0;
}


sub setIntronBoundaries{
    # description.
    # for each intron within alignment
    #this sets leftAlnPos,AlnPos,rightAlnPos which are defined by by maxSlide and maxGap
   
    my @last = (0,0,0,0);
    my $initGap = 0;
   
    if($pv){print "# print # introns ".scalar(@introns)." \n";}
    
    for (my $i=0;$i<scalar(@introns);$i++){
	my ($org,$annotationUID,$intronUID,$intronNum,$proteinPos,$phase,$spliceType,$ncFail) = @{$introns[$i]};

	my $pos;
	my $aapos;

	#print "$org,$annotationUID,$intronUID,$intronNum,$proteinPos,$phase,$spliceType,$ncFail\n";

	if ($proteinPos eq "5utr"){
		@{$introns[$i]} = (@{$introns[$i]}[0..6], "5utr","5utr","5utr");
		@last = ($org,$annotationUID,-1,0,$intronNum);
		next;
	}
	
	if ($proteinPos eq "3utr"){
		@{$introns[$i]} = (@{$introns[$i]}[0..6], "3utr","3utr","3utr");
		@last = ($org,$annotationUID,-1,0,$intronNum);
		next;
	}

	if ($last[0] eq $org && $last[1] eq $annotationUID and $intronNum > $last[4]){
	    $pos = $last[2]+1;
	    $aapos = $last[3];
	    $initGap = 1;
	}else{
	    $initGap = 0;	
	    $aapos = 0;  #zero position introns
	    $pos = 0;
	    $utr3cnt=0;
	    @last = ();
	}
        
	my $seq = $aln{$org}{$annotationUID};
	for(my $j=$pos;$j<length($seq);$j++){ # by characters in alignment
	    if (substr($seq,$j,1) eq "-"){
			next;
	    }else{
			$initGap = 1;  # inital gap passed
	    }	
	    $aapos++;
	    
	    if (($aapos == $proteinPos) || ($proteinPos ==0 && $aapos==1) ){ # found alignment position for intron by protein position.
		
			my ($leftAlnPos,$rightAlnPos,$slideCount,$gapCount);
			$leftAlnPos = $rightAlnPos = $j;
				
			#### find left and right boundaries.
				$slideCount = 0;
				$gapCount = 0;	
				### left boundary
				for(my $k=$j-1;$k>0;$k--){
				    if (substr($seq,$k,1) eq "-"){
						$gapCount++;
				    }else{
						$slideCount++;
				    }
				    if ($slideCount <= $maxSlide && $gapCount <= $maxGap){
						$leftAlnPos = $k;
				    }else{
						last;
				    }
				}
	
			### right boundary
				$slideCount = 0;
				$gapCount = 0;
				for(my $k=$j+1;$k<length($seq);$k++){
				    if (substr($seq,$k,1) eq "-"){
					$gapCount++;
				    }else{
					$slideCount++;
				    }
				    if ($slideCount <= $maxSlide && $gapCount <= $maxGap){
						$rightAlnPos = $k;
				    }else{
						last;
				    }
				}
			@{$introns[$i]} = (@{$introns[$i]}[0..6], ($leftAlnPos,$j,$rightAlnPos));
			@last = ($org,$annotationUID,$j,$aapos,$intronNum);
			last;
		    } #end intron found.
		} # end alignment search
    } # end intron for loop

    return;
}


sub fetchIntrons{
    my $cid = shift;
    @introns = ();
    %intronDat = ();

	open(IF,"$intronsFile") || die "cannot open $cid";
	while(<IF>){	

	chomp;
    my @f = split /\,/, $_;
    my ($intronNum,$Dsite,$Asite, $intronStart, $intronLength,$intronSeq,$phase,$proteinPos,$annotationUID,$chromosome,$organism,$strand) = @f;

	my $intronUID = "NULL";
	
	$introns[++$#introns] = [$organism,$annotationUID,$intronUID,$intronNum,$proteinPos,$phase,'U2',0];
	$intronDat{$#introns} = [$intronSeq,$Dsite,$Asite,$chromosome,$strand];
	$membersWintron{$org.$annotationUID}++;
    }
}

sub cFlanks{
    # this function returns similar character counts in flanking sequences
    my ($org1,$annotationUID1,$org2,$annotationUID2,$alnPosStart,$alnPosEnd) = @_;

    # alnPos is 1based
    #removed 9.24 my $c = 2;
    my $c=0;

    # intron Seg
    my $is1 = substr($aln{$org1}{$annotationUID1},$alnPosStart-$c,$alnPosEnd-$alnPosStart+1);
    my $is2 = substr($aln{$org2}{$annotationUID2},$alnPosStart-$c,$alnPosEnd-$alnPosStart+1);

    # flank Left
    my $fl1 = reverse substr($aln{$org1}{$annotationUID1},0,$alnPosStart-$c);
    my $fl2 = reverse substr($aln{$org2}{$annotationUID2},0,$alnPosStart-$c);

    # flank Right
    my $fr1 = substr($aln{$org1}{$annotationUID1},$alnPosEnd-$c+1);
    my $fr2 = substr($aln{$org2}{$annotationUID2},$alnPosEnd-$c+1);


    my ($intronSeg,$intronSegsize,$iLen,$igapcount,$simSeg,$iSingleGap) = compareString($is1,$is2,0,$alnPosEnd-$alnPosStart+1);
    my ($flankL,$flankLsize,$lLen,$lgapcount,$simL,$LSingleGap) = compareString( $fl1, $fl2,0,$flank);
    my ($flankR,$flankRsize,$rLen,$rgapcount,$simR,$RSingleGap) = compareString($fr1,$fr2,0,$flank);

    if ($pv){
    print "#cflanks $alnPosStart $alnPosEnd - $org1 $annotationUID1 vs. $org2 $annotationUID2 :\n";
    print "#$is1:$is2 $intronSeg/$intronSegsize ".substr($fl1,0,$lLen).":".substr($fl2,0,$lLen)." $flankL/$flankLsize ".substr($fr1,0,$rLen).":".substr($fr2,0,$rLen)." $flankR/$flankRsize  ($igapcount+$lgapcount+$rgapcount) $iSingleGap $LSingleGap $RSingleGap\n";
    }

    return ($intronSeg,$intronSegsize,$flankL,$flankR,$flankLsize,$flankRsize, ($igapcount+$lgapcount+$rgapcount),$simSeg,$simL,$simR,$iSingleGap+$LSingleGap+$RSingleGap );

}

sub compareString{
    # this function returns the amount of position specific characters in common between two string
    # if full == 1, the entire sequences are compared
    #    maxsearchlength = flank
    # otherwise, a segment of aligned residues of size flank including at most freeGaps gaps is compared
    #    maxsearchlength = flank + freeGaps
    #    max
 
    my ($s1,$s2,$full,$flank) = @_;
    my $count = 0; # number of identical non-gap characters
    my $nonGapPos = 0; # number of aligned characters
    my $gapCount = 0; # number of gaps in either sequence
    my $sGapCount=0; #number of gaps in one sequence but not both
    my $compGap = 0;  # number of gaps in second sequence
    my $tailGap = 0;  # number of ending gaps in sequence
    my $sim = 0; # number of similar residues between sequences
    my @s1R = $s1 =~ /./g;
    my @s2R = $s2 =~ /./g;
    my $j;
    for($j=0;$j<scalar(@s1R);$j++){

	    if ($s2R[$j] eq "-"){
		$compGap++;
		if ($j >= ($flank)){
		   $tailGap++;
	        }
	    }

	    if (($s1R[$j] eq "-" || $s2R[$j] eq "-")&&!($s1R[$j] eq "-" && $s2R[$j] eq "-")){
		$sGapCount++;
	    }

	    # if both characters are not gaps.  two gaps are not identical amino acids.
            if ($s1R[$j] ne "-" && $s2R[$j] ne "-"){

                $nonGapPos++;


                if ($s1R[$j] eq $s2R[$j]){
                    $count++;
                }
		if ($proteinMatrix->{$s1R[$j]}->{$s2R[$j]} > 0){
		    $sim++;
		}

	        if ( ($nonGapPos >= $flank)   && $full==0){
                    last;
                }

	    }else{

		#changed 9.24 if ( ($j >= ($flank+$freeGaps) || ($gapCount >= $freeGaps && $j >= $flank) )&& $full==0){

		$gapCount++;
	    }

		if ( ($j >= ($flank+$freeGaps) || ($gapCount+$nonGapPos >= $freeGaps + $flank) )&& $full==0){
		    last;
		}

    }
#    my $fragLength = ($full==0) ? max($nonGapPos,$flank) : $nonGapPos;
     my $fragLength = $flank;  #always flank positions, even if near start or end of gene.
    return ($count,$fragLength,$j+1,$compGap-$tailGap,$sim,$sGapCount);
}

sub sort_utr3{
#  my ($a,$b) = @_;
  my ($an) = $a->[0] =~ /utr\.(\d+)/;
  my ($bn) = $b->[0] =~ /utr\.(\d+)/;
  $an <=> $bn
}

sub sort_utr5{
#  my ($a,$b) = @_;
  my ($an) = $a->[0] =~ /utr\.(\d+)/;
  my ($bn) = $b->[0] =~ /utr\.(\d+)/;
  $bn = int($bn);
  $an = int($an);
  $bn <=> $an
}

sub sort_lr{
  $a->[0] <=> $b->[0]
	||
  $a->[1] <=> $a->[1]
}

sub sortIntronPos{
  my $g = shift;
  my $utr3;
  my $utr5;
  my $pc;
  for (my $i=0;$i<scalar(@$g);$i++){
    if ($g->[$i]->[0] =~ /5utr/){
	$utr5->[++$#$utr5] = $g->[$i];
    }elsif($g->[$i]->[0] =~ /3utr/){
	$utr3->[++$#$utr3] = $g->[$i];
    }else{
	$pc->[++$#$pc] = $g->[$i];
    }
  }
  @$utr5 = sort sort_utr5 @$utr5;
  @$utr3 = sort sort_utr3 @$utr3;
  @$pc = sort sort_lr @$pc;
  $g = [@$utr5,@$pc,@$utr3];
  return $g;
}

sub loadAln{
my $cid = shift;

%aln = ();
my $fail=0;
open(IF,"$clusterFile") || die "no file $cid $infile";

my $tmp = <IF>; # clustal title line
# pass on single gene clusters for now
if ($tmp !~ /MUSCLE/ && $tmp !~ /CLUSTAL/){
  die "bad alignment parse $cid\n";
  next;
}
while (<IF>){
  if ($_ =~ /^\w/){ # skip space

	my ($org,$id,$str) = parseId($_);
    print "algn $org $id $str\n";
    $aln{$org}{$id} .= $str;

  }
}

if (keys(%aln)==0){die "bad alignment parse\n";}

%bounds = ();

foreach my $j (keys(%aln)){
    foreach my $k (keys (%{$aln{$j}})){
	my ($l,$m,$r) = $aln{$j}{$k} =~ /^(-*)([^-].+?[^-])(-*)$/;
	$r = length($l.$m);
	$l = length($l);
	#test print "$j $k $l $r\n";
	$bounds{$j}{$k}{left} = $l;
	$bounds{$j}{$k}{right} = $r;

    }
}


close IF;

}
