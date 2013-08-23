#!/usr/bin/perl

use CGI ":all";
use Cwd;
my $tdir = cwd;

#($rpath) = $tdir =~ /(.+?)cgi/;

require "conf.pl";


my $cUID = param("cid");
my $ciUID = param("ciUID");

print STDERR "$cUID $ciUID\n";

my $sql = "select org,annotationUID,segment,start,stop,intronType,phase,intronSeq from ClusterIntron2Gene$cVer where clusterUID = '$cUID' and
 clusterIntronUID = $ciUID and intronType != 'em' and intronType != 'absent' and intronType != 'im'  ";
my $ref = $NRIdbh->selectall_arrayref($sql);

my $t="";
for(my $i=0;$i<scalar(@$ref);$i++){
        #$t .= "<tr><td>".join("</td><td>",@{$ref->[$i]})."</td></tr>\n";
        my $x = $ref->[$i];

        my $is = (length($x->[7])>50) ? substr($x->[7],0,25)."...".substr($x->[7],length($x->[6])-25) : $x->[7];

        $t .= "<tr>
                        <td>$x->[0]</td>
                        <td>$x->[1]</td>
                        <td><a href='".genomeLink($x->[0],$x->[2],$x->[3],$x->[4])."'> $x->[2] </a></td>
                        <td>$x->[3]</td>
                        <td>$x->[4]</td>
                        <td>$x->[5]</td>
                        <td>$x->[6]</td>

                        <td>$is</td>
                        </tr>";
}

print header();

print "<html>
<style>
.header{font-weight:bold; background-color:lightgray;}
.tab{font-family:Courier;font-size:12px;border-inset:black 1px solid;}
</style>

";

print "<body>
CIWOG
<table border=0 class=tab>
        <tr>
                        <td class=header>Organism</td>
                        <td class=header>Gene ID</td>
                        <td class=header>Genome Segment</td>
                        <td class=header>Intron Start</td>
                        <td class=header>Intron Stop</td>
                        <td class=header >Intron Type</td>
                        <td class=header>Intron Phase</td>
                        <td class=header >Intron Sequence</td>
                        </tr>
$t
</table></body></html>";



