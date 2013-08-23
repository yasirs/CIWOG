#!/usr/bin/perl

use CGI ":all";
use Cwd;
my $tdir = cwd;


($rpath) = $tdir =~ /(.+?)cgi/;

#require "${rpath}conf.pl";
require 'conf.pl';

$cid = param('cid');


#$aln = "fsa.mus";
$infile = $ALN_DIR."/$cid/$cid.aln";


my @seqOrder; # for alignments intron table
my @seqOrder_alt;
my @altg;
my $cref; #cogs array reference
my $iref; #clusterintrons ref
my ($utr3,$utr5) = 0;


my $maxl;

print header();
print "<html>";
print "<head>
<style>
.itable {font-family:'courier new','courier';font-size:12px;border:1px solid black;}
.altable {font-family:'courier new','courier';font-size:12px;letter-spacing:1px;padding:15px}
.intron1 {color:red}
.intron2 {color:purple}
.btable{font-family:'verdana','arial';font-size:8pt;empty-cells:show;border:1px solid gray;}
.btable TH {text-decoration:none;font-weight:bold;font-size:8pt;}
.btable TD {border-top:gray 1px solid;}
.t1{font-family:'verdana','arial';font-size:18pt;font-weight:bold}
.t2{font-family:'verdana','arial';font-size:10pt}

</style>

<script>
    function hintron(uid,color){
	document.getElementById('r'+uid).style.backgroundColor = (color == '') ? 'yellow' : 'white';
	document.getElementById('i'+uid).style.backgroundColor = (color == '') ? 'yellow' : color;
	return;
    }

function hrow(uid,color){
    hintron(uid,color);
    return;

}

</script>

</head>";

@altColors = ("red","pink","purple","violet","green","blue","lightblue","brown","gray","beige");


my ($map) = writeAlnGraph(); # draws image
($idArr,$seqArr) = writeAln();
#$geneTable = writeGeneTable();
($itArr) = writeIntronsTable();

print "<body>";
my $titleT = "$globalTitle Cluster ";
print "<table border=0><tr><td width = 500 valign='top'><a href='$homeURL' ><img border=0 src='${HTMLurl}$imgTitle'> <font class=t2>home</font> </a> 
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;

<font class=t1 >$titleT $cid</font></td> ";

#cluster summary table

my $sql =  "select count(distinct annotationUID), count(distinct clusterIntronUID), group_concat(distinct org,' '), group_concat(distinct intronType,' ') from ClusterIntron2Gene where clusterUID = '$cid' ";
my $ref = $NRIdbh->selectall_arrayref($sql);

print "<td><table align='right' class=t2 bgcolor='lightgray' width=200>";
print "<tr><td>Organisms</td><td>$ref->[0]->[2]</td></tr>";
print "<tr><td>Gene count</td><td>$ref->[0]->[0]</td></tr>";
print "<tr><td>Cintron count</td><td>$ref->[0]->[1]</td></tr>";
print "<tr><td>Cintron classes</td><td>$ref->[0]->[3]</td></tr>";
print "</table></td></tr>";

print "<tr><td colspan=2>";
print "$map<img src='${HTMLurl}tmp/aln$cid.$cVer.png' border=0 useMap='#alnMap$cid'></td>";
#print "<td>
#<font class=t2>Cluster Summary</font>
#".writeClusterSummary()."<br>
#<font class=t2>Gene Summary</font>
#$geneTable</td>";

print "</tr></table><br><br>";


################
### Make Display
################

my $ind = 0;


my $scale = "    .    :" x 6;
my $ftable;
$ftable = "\n\n<table class=t2>";

## add legend
$ftable .= "<tr><td></td><td></td><td>$tableIntronHash->{'legend'}</td></tr>";

## add 3utr
$ftable .= "<tr><td></td><td></td><td ><br><br>";
$ftable .= ($utr5==1) ? "<a href='displayIntronGroup.pl?cid=$cid&ciUID=".(1)."' ><font class=t2>5' UTR introns</font></a>": "" ;
$ftable .= "</td></tr>";


### add legend
$ftable .= "<tr>
				<td class=altable></td>
				<td class=altable></td>
				<td>
				<table class=itable>

	    					<tr>
							<td>Cluster Intron #: Alignment Start-Stop<br></td>
						</tr><tr>
		    					<td !width=200>Gene ID, Intron #, Phase,Intron Type,Intron Seq
		    				</td></tr>
				</table>
				</td>
				</tr>";


for(my $i=0;$i<scalar(@$seqArr);$i++){
	$ftable .= "<tr class=itable>
				<td class=altable>$idArr->[0]</td>
				<td class=altable>$seqArr->[$i]</td>
				<td >$itArr->[$i]</td>
				</tr>";
}


## add 3utr
$ftable .= "<tr><td></td><td></td><td><br><br>";
$ftable .= ($utr3==1) ? "<a href='displayIntronGroup.pl?cid=$cid&ciUID=".(1+length(keys(%$tableIntronHash)))."' ><font class=t2>3' UTR introns</font></a>": "" ;
$ftable .= "</td></tr>";



$ftable .= "</table>";

print $ftable;
print "<br>"x20;
print "</body></html>";

sub writeClusterSummary{
  my $sql = "select geneCount,distinctOrgL,distinctIntronCount,avgIntronAmount,intronRange,u12GeneCount from Clusters$cVer where clusterUID = $cid";
  my @r = $NRIdbh->selectrow_array($sql);
  my $table = "<table class='btable'>";
  $table .= "<tr><td><b>Gene Count:</b> $r[0]</td></tr>";
  $table .= "<tr><td><b>Intron Sites:</b> $r[2]</td></tr>\n";
#  $table .= "<tr><td><b>Intron Range:</b> $r[4]</td></tr>\n";
#  $table .= "<tr><td><b>Organisms:</b> $r[1]</td></tr>";
  #<td><b>Avg. Intron #:</b> $r[3]</td></tr>\n";
  $table .= "<tr><td><b>Related Uniref Descriptions:</b><br>".writeUniref()."</td></tr>";

#  $table .= "<tr><td><b>Gene w/ u12 #:</b> $r[5]</td><td><b>Member w/ multiple u12:</b> $r[6]</td></tr>\n";
#  $table .= "<tr><td><b>Conserved # u12:</b> $r[8]</td><td><b>Conserved # u12 introns <br> among organisms:</b> $r[10]</td></tr>\n";
#  $table .= "<tr><td><b>Avg. % Identity:</b> $r[11]</td><td><b>Avg % Coverage:</b> $r[12]</td></tr>\n";
  $table .= "</table>"; 

  return $table;
}


sub writeAln{
	my $sql = "select distinctrow annotationUID,t1.org,phase,intronNum,segment,start,stop,alnPos,segment,intronSeq,t1.clusterIntronUID, intronType,alnPosStart,alnPosEnd from ClusterIntron2Gene$cVer as t1 inner join ClusterIntrons$cVer as t2 using(clusterUID,clusterIntronUID) where t1.clusterUID = '$cid' order by t1.org,t1.annotationUID,t1.intronNum";
	$iref = $NRIdbh->selectall_arrayref($sql);
	
	%geneIpos; #used in intron Table
	my %geneCpos;
	my %iType;
	my %geneDetail;
	my %intronType;
	my @idArr;
	for(my $k=0;$k<scalar(@$iref);$k++){
	  my ($annotationUID,$org,$phase,$intronNum,$segment,$start,$stop,$alnPos,$chr,$intronSeq,$clusterIntronUID, $intronType) = @{$iref->[$k]};
	  $geneIpos{$org}{$annotationUID}{$alnPos} = $k;
	  #print STDERR "$org $annotationUID $alnPos\n";	  
	}
	#parse alignment
	%aln = ();
	my $fail=0;
	open(IF,"$infile") || die "no file $cid $infile";
	my $tmp = <IF>; # clustal title line
	# pass on single gene clusters for now
	if ($tmp !~ /MUSCLE/ && $tmp !~ /CLUSTAL/){
	  die "bad alignment parse $cid\n";
	  next;
	}
	while (<IF>){
	  if ($_ =~ /^\w/){ # skip space
	
	    my ($org,$id,$str) = parseId($_);
		if (!defined($aln{$org}{$id}) ){
			$seqOrder[++$#seqOrder] = $org."~".$id;
	    }
	    $aln{$org}{$id} .= $str;	    
	    
	  }
	}
	close IF;

	# find maximum of alignment
	$maxl=0;
	foreach my $k (keys (%aln)){
	  foreach my $j (keys (%{$aln{$k}})){
	    $maxl = ($maxl < length($aln{$k}{$j})) ? length($aln{$k}{$j}) : $maxl;
	  }
	}
	
	my @seqarr;
	my $j;
	$cnt = -1;
	for ($l=0;$l<$maxl;$l=$l+60){ # for each block of seq AA
		$cnt++;
	  	foreach my $y (@seqOrder){  # in order of alignment 
		    my ($k,$j) = $y =~ /^(\S+?)~(\S+)/;
		    #print STDERR "se $k $j\n";
		    $idArr[$cnt] .= $y."<br>\n";
		    my $as="";
			for (my $i=$l;$i<$l+60;$i++){ # for each individual base
				
				#search to see if there is an intron for this gene.
				
				my $b = substr($aln{$k}{$j},$i,1);
				my $ih = $geneIpos{$k}{$j}{$i};
				if ($ih ne ""){
					my $type = $iref->[$geneIpos{$k}{$j}{$i}]->[11];
					my $col = $intronColors{$type};
					
					$jscript = "style='cursor:pointer' onmouseover=\"hintron('".$k.$j.$iref->[$ih]->[3]."','');\" onmouseout=\"hintron('".$k.$j.$iref->[$ih]->[3]."','$col');\" ";
					$b="<font style='background:$col' ><a id='i".$k.$j."".$iref->[$ih]->[3]."' $jscript>$b</a></font>";
				}
				$as .= "$b";
				#print STDERR "$k $j $geneIpos{$k}{$j}{$i}\n";
	  		}
			$seqarr[$cnt] .= "$as<br>\n";
	  	}
	}

	return (\@idArr,\@seqarr);

}



sub link_to_browser(){
  return;
}

sub writeIntronsTable{
	#for geneIpos by 60s,
	
	#find clusterintrons by alnPosStart
	my %ci;
	
	my $cnt=-1;
		
	my @it;
	
	for(my $j=0;$j<scalar(@$iref);$j++){
		#print STDERR "$iref->[$j]->[12]\n";
		$ci{$iref->[$j]->[12]} = [@{$ci{$iref->[$j]->[12]}}, $j];
		if($iref->[$j]->[12] eq "5utr" & !$utr5){
			$utr5=1;
    						
		}
		if($iref->[$j]->[12] eq "3utr"){
			$utr3=1;
		}
		
	}


	for ($l=0;$l<$maxl;$l=$l+60){ # for each block of seq AA
		$cnt++;
	    for (my $i=$l;$i<$l+60;$i++){ # for each individual base
	    		my %r={};
	    	if (exists($ci{$i})){
	    		for(my $cir=0;$cir<scalar(@{$ci{$i}});$cir++){
	    			my $x = $iref->[$ci{$i}->[$cir]];
	    			
	    			if($x->[11] eq "absent" || $x->[11] eq "im" || $x->[11] eq "em"){
	    				next;
	    			}
	    			my $col = $intronColors{$x->[11]};
	    			my $markup = "id='r$x->[1]$x->[0]$x->[3]' style='cursor:pointer' onmouseover=\"hrow('$x->[1]$x->[0]$x->[3]','');\" onmouseout=\"hrow('$x->[1]$x->[0]$x->[3]','$col');\" ";
	    			
		    		$r{$x->[1]."~".$x->[0]}= "
		    				<tr $markup>
		    					<td>".$x->[0]."</td>
		    					<td>".$x->[3]."</td>
		    					<td>".$x->[2]."</td>
		    					<td>".$x->[11]."</td>
		    					<td>".substr($x->[9],0,8)."..".substr($x->[9],length($x->[9])-4)."</td>	    					
		    					
		    				</td></tr>
		    				";
	    		}
	    		my $defspacing = "<tr>
		    					<td width=200></td>
		    					<td width=10></td>
		    					<td width=10></td>
		    					<td width=10></td>
		    					<td width=20></td>	    					
		    				</td></tr>";
		    				
				#order by seq order
				my $nr = "";
				foreach my $k (@seqOrder){
					$nr .= $r{$k};
				}
		    	my $x = $iref->[$ci{$i}->[0]]; #one for cluter intron data
		    				
	    		$it[$cnt] .= "<table class=itable>
	    						$defspacing
	    						<tr><td><a href='displayIntronGroup.pl?cid=$cid&ciUID=$x->[10]'>$x->[10]</a>: $x->[12]-$x->[13] </td></tr>
	    						$nr
	    						</table>";
	    		
	    	}
		}
	}


	

	return \@it;
}






sub writeAlnGraph{
   require "aln_graph.pl";
   my ($map) = aln_data_graph($cid,$infile); 
   my $map = `cat ../html/tmp/aln$cid.$cVer.imap.html`;
   return ($map);
}
sub by_range_sort{
    my @a1 = split /\-/, $a;
    my @b1 = split /\-/, $b;
    $a1[0] <=> $b1[0];
}

sub by_utr_sort{
    my @a1 = split /utr\./, $a;
    my @b1 = split /utr\./, $b;
    $a1[1] <=> $b1[1]
}


sub by_first_sort{
    # this sorts intron records by org, then altgroup
    my @a1 = @$a;
    my @b1 = @$b;
    $a1[0] cmp $b1[0]
	||
    $a1[1] <=> $b1[1]
}
