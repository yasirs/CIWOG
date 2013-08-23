#!/usr/bin/perl
use DBI;
use Switch;
use CGI ":all";

$LOCAL_DIR="/ciwog/scripts/"; #location of ciwog scripts directory

$OUTPUT_DIR="/ciwog/output/"; #location of ciwog output directory
$INPUT_DIR="/ciwog/input/"; #location of ciwog input directory
$ALN_DIR="/ciwog/input/"; #location of ciwog alignments, usually the same as INPUT_DIR

$HTMLpath="/ciwog/html/"; #location of ciwog html directory
$CGIurl="http://localhost/ciwog-cgi/"; #URL configured in httpd.conf
$HTMLurl="http://localhost/ciwog/"; #URL configured in httpd.conf

$NRIdbh=DBI->connect("DBI:mysql:ciwog:localhost","ciwogUser","",{PrintError=>1,RaiseError=>1}); #MySQL database connection string

require "${LOCAL_DIR}ciwog.functions.pl";
$imgTitle = "ciwogPlants.logo.jpg"; #logo for website
$homeURL = "${CGIurl}home.pl";

%intronColors = ("U2" => 'orange', "U12" => 'red', "GC" => '#fb3deb', "absent" => "lightgray");


$min_u12_donor = 4;
$min_u12_branch = 4;



sub genomeLink{
  my ($org,$chr,$l,$r) = @_;
  my $link;
  switch ($org){
    case 'AT' { $link = "http://www.plantgdb.org/AtGDB/cgi-bin/getRegion.pl?dbid=3&chr=$org&l_pos=$l&r_pos=$r"; }  #URLs for genomic context.
    
  }
    
  return $link

}


