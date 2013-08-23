#!/usr/local/bin/perl
use GD;
use GD::Text::Align;

#($rpath) = $tdir =~ /(.+?)cgi/;
require "conf.pl";


#aln_data_graph("cluster1.txt","/home/mwilkers/workspace/ciwogDist/output/cluster1.txt/cluster1.txt.clw");

 
sub aln_data_graph{
my $cid = shift; #clusterUID.
my $infile = shift; #alignment

# parse aln for seq position and gaps

my %aln;
my $fail=0;
my $alnLen = 0;
my @seqOrder;
my %seqNum;
open(IF,"$infile") || die "error with mus file $infile";
my $seqNumCnt = 0;
my $tmp = <IF>; # clustal title line
while (<IF>){
  if ($_ =~ /^\w/){ # skip space

    my ($org,$id,$str) = parseId($_);

    #my $id = "$source-$id";
    $aln{$org}{$id} .= $str;
    if (!defined($seqNum{$org.$id})){
      $seqNum{$org.$id} = $seqNumCnt;
      $seqNumCnt++;
      $seqOrder[++$#seqOrder] = [$org,$id];
    }
  }
}
close IF;



my @seq;
my @gaps;
#my $seqNumCnt = 0;
my @labels;
 @orgLabels;
 @geneLabels;
#foreach my $k1 (sort keys(%aln)){
#  foreach my $k2 (sort keys(%{$aln{$k1}}) ){
#    $seqNum{$k1.$k2} = $seqNumCnt;

foreach my $i (@seqOrder){
    my ($k1,$k2) = @$i;

    my $seqNumCnt = $#seq + 1;
    $seqNum{$k1.$k2} = $seqNumCnt;
    my ($gap1,$mid,$gap2) = $aln{$k1}{$k2} =~ /^(-*)([^-].+?)(-*)$/;
    my $s = length($gap1);
    my $e = $s + length($mid);
    my $len = $e-$s;
    $alnLen = max( $alnLen,length($aln{$k1}{$k2}) );

    $seq[$seqNumCnt] = [$len,$s,$e];
    $labels[$seqNumCnt] = $k1.$k2;
    $orgLabels[$seqNumCnt] = $k1;
    $geneLabels[$seqNumCnt] = $k2;
    
    my @gArr;
    my $mr = $s;
    my $sp;
    my $sg;
    while ( ($sp,$ig) = $mid =~ /^([^-]+?)(-+)/){
     # print "$sp\n\t$ig\n";
      $gArr[++$#gArr] = $mr + length($sp);
      $gArr[++$#gArr] = $mr + length($sp)+length($ig);
      $mr  = $mr + length($sp)+length($ig);
      $mid = substr($mid,length($sp)+length($ig));
     # last;
    }
    $gaps[$seqNumCnt] = [@gArr];
     
#    $seqNumCnt++; 
#  }
}


# fetch cluster introns
$sql = "select t2.org,annotationUID,alnPosStart,t1.clusterIntronUID,t2.intronType,intronType,intronNum from ${queryDb}ClusterIntrons$cVer$trialRep as t1 inner join ${queryDb}ClusterIntron2Gene$cVer$trialRep as t2 on (t1.clusterUID = t2.clusterUID and t1.clusterIntronUID = t2.clusterIntronUID) where t1.clusterUID = '$cid' and intronType not in ('absent','im','em') order by t1.clusterIntronUID";
my $ref = $NRIdbh->selectall_arrayref($sql);
my @ci;
my @u5 = (0);
my @u3 = (0);
my @u5_u12 = (0);
my @u5_dual = (0);
my @u3_u12 = (0);
my @u3_dual = (0);

foreach my $i (@$ref){
  my ($org,$annotationUID,$alnPos,$clusterIntronUID,$intronType,$mod_intronType,$intronNum) = @$i;
#  print STDERR join(",",@$i)."\n";
  if ($alnPos =~ /5utr/){
	$u5[$seqNum{$org.$annotationUID}+1] += 1;
	if ($mod_intronType eq "u12"){
		$u5_u12[$seqNum{$org.$annotationUID}+1] += 1;
	}elsif($mod_intronType eq "dual"){
		$u5_dual[$seqNum{$org.$annotationUID}+1] += 1;
	}

  }elsif($alnPos =~ /3utr/){
	$u3[$seqNum{$org.$annotationUID}+1] += 1;
	if ($mod_intronType eq "u12"){
		$u3_u12[$seqNum{$org.$annotationUID}+1] += 1;
	}elsif($mod_intronType eq "dual"){
		$u3_dual[$seqNum{$org.$annotationUID}+1] += 1;
	}
  }else{
	$ci[$clusterIntronUID-1][0] = $alnPos;
	$ci[$clusterIntronUID-1][$seqNum{$org.$annotationUID}+1] = [$intronType,$mod_intronType,$intronNum];
  }
}

#test 
print STDERR "$alnLen $#seq $#gaps $#ci\n";
draw_aln_graph($alnLen,\@seq,\@gaps,\@ci,\@labels,\@u5,\@u3,\@u5_u12,\@u5_dual,\@u3_u12,\@u3_dual);
}

sub draw_aln_graph{

my ($alnLen,$seqRef,$gapsRef,$ciRef,$labelRef,$u5Ref,$u3Ref,$u5_u12Ref,$u5_dualRef,$u3_u12Ref,$u3_dualRef) = @_;

my @seq = @$seqRef;
my @gaps = @$gapsRef;
my @ci = @$ciRef;
my @labels = @$labelRef;
my @u5 = @$u5Ref;
my @u3 = @$u3Ref;
my @u5_u12 = @$u5_u12Ref;
my @u5_dual = @$u5_dualRef;
my @u3_u12  = @$u3_u12Ref;
my @u3_dual = @$u3_dualRef;

$path = "${HTMLpath}/tmp/";
$iname = "aln$cid.$cVer$trialRep";


my $map = "<map name='alnMap$cid'>\n";

#example input

# length, start, end
@xseq = ( [3000,0,3000],
	 [2500,500,3000],
	 [2000,100,2100],
	 [1000,2000,3000]
);

# gaps, row is seq, 0 is start, 1 is end, ...
@xgaps = ( [10,200,500,1000],
	[],
	[],
	[]
	);


# row is intron group, column is alnStart, then sequence possession or not
@xci = ( [2100,1,1,1,1],
	[2200,0,0,0,1],
	[2700,1,1,0,0]
);


$pwidth = 1200;
$pheight = 100;
$image = new GD::Image($pwidth,$pheight);

# colors
$white = $image->colorAllocate(255,255,255);
$black = $image->colorAllocate(0,0,0);       
$red = $image->colorAllocate(255,0,0);      
$blue = $image->colorAllocate(0,0,255);
$orange = $image->colorAllocate(249,166,5);
$yellow = $image->colorAllocate(246,249,5);
$magenta = $image->colorAllocate(223,5,249);
$purple = $image->colorAllocate(133,31,236);
$green = $image->colorAllocate(10,136,10);
$lblue = $image->colorAllocate(6,241,252);
$lgreen = $image->colorAllocate(6,251,17);
$grey = $image->colorAllocate(161,161,161);
$brown = $image->colorAllocate(148,98,10);
$lgrey = $image->colorAllocate(207,207,207);

@cycle = ($red,$orange,$yellow,$magenta,$purple,$green,$lblue,$lgreen,$brown);


# org colors
%colors = ( 
	    "at" => $blue, 
	    "os" => $red,
	    "hs" => $magenta,
	    "mm" => $purple,
	    "rn" => $green,
	    "fr" => $orange,
	    "dr" => $lgreen,
	    "gg" => $lblue,
	    "dm" => $yellow
	    );
# sort by org


$pad = 175; # label padding
$vpad = 15;
#$sHeight = ($pheight-$vpad*2) / scalar(@seq) / 2;
$sHeight = 4;
$vSpace = $sHeight+15;

$utrPad = 50;
$hScale = $alnLen / ($pwidth - 1 *$pad - $utrPad*2);

$intronWidth = 2;
$intronB = 5; #intron vertical height


# draw sequences
my $r=$vpad;
for( my $i=0;$i<scalar(@seq);$i++){
  my ($length,$start,$end) = @{$seq[$i]};
  
  # add utr number;
  $image->string(gdSmallFont,$pad,$r,"$u5[$i+1]",$black);
  $image->string(gdSmallFont,$pwidth-$utrPad+2,$r,"$u3[$i+1]",$black);

  # add utr markers
  if ($u5_u12[$i+1] >= 1){
    $image->string(gdSmallFont,$pad,$r-$vSpace/2-$intronB,"*",$black);
  }elsif($u5_dual[$i+1] >= 1){
    $image->string(gdSmallFont,$pad,$r-$vSpace/2-$intronB,"*",$brown);
  }

  if ($u3_u12[$i+1] >= 1){
    $image->string(gdSmallFont,$pwidth-$utrPad+2,$r-$vSpace/2-$intronB,"*",$black);
  }elsif($u3_dual[$i+1] >= 1){
    $image->string(gdSmallFont,$pwidth-$utrPad+2,$r-$vSpace/2-$intronB,"*",$brown);
  }

  # draw sequence
  $length = $length / $hScale;
  $start = $start / $hScale;
  $end = $end / $hScale;
  $image->filledRectangle($pad+$start+$utrPad,$r,$pad+$utrPad+$end,$r+$sHeight,$black);  

  # add gaps
  for (my $j=0;$j<scalar(@{$gaps[$i]});$j=$j+2){
    if (($gaps[$i][$j+1] - $gaps[$i][$j]) > 10){ 
      my $gs = $pad  + $utrPad + $gaps[$i][$j]/$hScale;
      my $ge = $pad  + $utrPad + $gaps[$i][$j+1]/$hScale;
      $image->filledRectangle($gs,$r,$ge,$r+$sHeight,$lgrey);
    }
  }

  # add introns
  for (my $k=0;$k<scalar(@ci);$k++){
#    my $icolor = ($k % 2 == 0) ? $yellow : $orange;
    my $icolor = $cycle[$k%scalar(@cycle)];
#    if ($ci[$k][$i+1] ne ""){
     if (scalar(@{$ci[$k][$i+1]}) > 0){


      # glyphs are only used for mod_intronType
      if (0 && $ci[$k][$i+1][0] eq "u12" || $ci[$k][$i+1][0] eq "dual"){
  	  my $dotColor = ($ci[$k][$i+1][0] eq "u12") ? $black : $brown;
	  my $gdSmallFontWidth=5;
	  $image->string(gdSmallFont,$ci[$k][0]/$hScale + $pad + $utrPad -$gdSmallFontWidth/2,$r-$vSpace/2-$intronB,"*",$dotColor);
#	  $image->filledArc($ci[$k][0]/$hScale + $pad + $utrPad+$intronWidth/2,$r+$sHeight/2,$intronWidth*4,$intronWidth*4,0,360,$white);
#	  $image->filledRectangle($ci[$k][0]/$hScale + $pad + $utrPad-2*$intronWidth,$r+$sHeight/3,$ci[$k][0]/$hScale + $pad + $utrPad + 2*$intronWidth,$r+2*$sHeight/3,$dotColor);
      }

      #marker for mod_type
      if ($ci[$k][$i+1][0] eq "U12" || $ci[$k][$i+1][0] eq "GC" ){
          my $dotColor = ($ci[$k][$i+1][1] eq "u12") ? $black : $black;
          my $gdSmallFontWidth=2; #5
          #$image->string(gdTinyFont,$ci[$k][0]/$hScale + $pad + $utrPad -$gdSmallFontWidth/2,$r-$vSpace/2-$intronB,"*",$dotColor);
	  $image->filledEllipse($ci[$k][0]/$hScale + $pad + $utrPad + 1,$r-$intronB,8,8,$icolor);
      }


      $image->filledRectangle($ci[$k][0]/$hScale + $pad + $utrPad,$r-$intronB,$ci[$k][0]/$hScale + $pad + $utrPad + $intronWidth,$r+$sHeight+$intronB,$icolor);
      my $c1 =$ci[$k][0]/$hScale + $pad + $utrPad;
      my $c2=$r-$intronB;
      my $c3 = $ci[$k][0]/$hScale + $pad + $utrPad + $intronWidth;
      my $c4 = $r+$sHeight+$intronB;
      $map .= "<AREA SHAPE='RECT' COORDS='$c1,$c2,$c3,$c4' HREF='";
      $map .= "${CGIurl}display.pl?cid=$cid";
      $map .= "#i".$labels[$i].$ci[$k][$i+1][2]."'  onmouseover='' onmouseout=''>\n";
    }
  }
  # add label
  $image->string(gdTinyFont,0,$r - $vSpace/2,$orgLabels[$i],$black);
  $image->string(gdSmallFont,0,$r - $sHeight,$geneLabels[$i],$black);



  $r += $vSpace+$sHeight;
  if ($r +$vSpace + $sHeight > $pheight){
	my $nh = $r +$vSpace + $sHeight;
	my $tmpIm=new GD::Image($pwidth,$nh);
        define_colors($tmpIm);
        $tmpIm->copy($image,0,0,0,0,$pwidth,$nh);
	$image = $tmpIm;
	$pheight = $nh;
  }
}


open (OF, ">${path}$iname.png") || die "cannot open ${path}$iname.png";
binmode OF;
print OF $image->png;
close OF;
$map .= "</map>";

open (MF, ">${path}$iname.imap.html") || die "cannot open >${path}$iname.imap.html";
print MF "$map";
close MF;

return $map;

}


sub define_colors{
  my $image = shift;
# colors
$white = $image->colorAllocate(255,255,255);
$black = $image->colorAllocate(0,0,0);
$red = $image->colorAllocate(255,0,0);
$blue = $image->colorAllocate(0,0,255);
$orange = $image->colorAllocate(249,166,5);
$yellow = $image->colorAllocate(246,249,5);
$magenta = $image->colorAllocate(223,5,249);
$purple = $image->colorAllocate(133,31,236);
$green = $image->colorAllocate(10,136,10);
$lblue = $image->colorAllocate(6,241,252);
$lgreen = $image->colorAllocate(6,251,17);
$grey = $image->colorAllocate(161,161,161);
$brown = $image->colorAllocate(148,98,10);
$lgrey = $image->colorAllocate(207,207,207);
  
}

1;
