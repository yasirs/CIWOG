#!/usr/bin/perl

use CGI ":all";

require 'conf.pl';

my $ci = param("cid");
my $cii = param("cii");



print header;

if ($ci eq ""){
#	exit();
}

my $sql = "select clusterUID,clusterIntronUID,intronNum,org,annotationUID,intronSeq from ClusterIntron2Gene where intronType not in ('absent','im','em') ";
$sql .=  (!$ci > 0 )? "" : " and clusterUID = '$ci'";
$sql .=  (!$cii > 0 )? "" : " and clusterIntronUID = $cii";

my $ref = $NRIdbh->selectall_arrayref($sql);

print "<html><body><pre>#format is >cluster ID|cintron ID|organism|gene ID|intron number\n\n";

foreach my $i (@$ref){
  my ($clusterUID, $clusterIntronUID, $intronNum, $org, $annotationUID, $intronSeq) = @$i;
  print ">$clusterUID\|$clusterIntronUID\|$org\|$annotationUID\|$intronNum\n$intronSeq\n";
}

print "</pre></body></html>";
