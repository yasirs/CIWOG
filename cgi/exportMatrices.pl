#!/usr/bin/perl

use CGI ":all";

require 'conf.pl';

my $ci = param("cid");
my $cii = param("cii");



print header;

if ($ci eq ""){
#	exit();
}

my $sql = "select clusterUID,clusterIntronUID,intronNum,org,annotationUID,intronType from ClusterIntron2Gene where alnPos not regexp 'utr' order by clusterUID, org,annotationUID,clusterIntronUID";
$sql .=  (!$ci > 0 )? "" : " and clusterUID = '$ci'";
$sql .=  (!$cii > 0 )? "" : " and clusterIntronUID = $cii";

my $ref = $NRIdbh->selectall_arrayref($sql);

print "<html><body><pre>#format:\n ##ClusterUID.\n# cintron classes where 0=absent,1=U2,2=GC,3=U12,4=im,5=em\n";

my %code = ("absent" => 0, "U2" => 1, "GC" => 2, "U12" => 3, "im" => 4, "em" => 5);

my $lastGene = -1;
my $lastc = -1;
my $txt = "";

foreach my $i (@$ref){
  ($clusterUID, $clusterIntronUID, $intronNum, $org, $annotationUID, $intronType) = @$i;
  if ($clusterUID eq $lastc){
    if ($annotationUID eq $lastGene){
	$txt .= "\t$code{$intronType}";
    }else{
	$txt .= "\n$org~$annotationUID\t$code{$intronType}";
	$lastGene = $annotationUID;
    }
  }else{
    if ($lastc != -1){
	    print "#$lastc\n$txt\n\n";
    }
    $txt = "";
    $lastc = $clusterUID;
    $lastGene = $annotationUID;
    $txt .= "\n$org~$annotationUID\t$code{$intronType}"
  }
}

print "#$clusterUID\n$txt\n\n";

print "</pre></body></html>";


