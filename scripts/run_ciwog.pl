#!/usr/bin/perl

require 'conf.pl';
require 'ciwog.functions.pl';

#parse introns and load intron and gene tables
my $res = readDirect($INPUT_DIR,$OUTPUT_DIR);


#clear database
print `cat ../ciwog.sql | mysql ciwog -u ciwogUser`;


###
# CIDA options
###
my $maxSlide=1;	
my $maxGap=1;
my $simRegion=10;
my $minSim=0.1;
my $minAbsentSim=0.3;
###
# end CIDA options
###

# detect_common_introns and load into database
  opendir(DIR,$INPUT_DIR) || die "cannot open directory $INPUT_DIR";
  my @f = grep /\w/,readdir(DIR);
  closedir DIR;
  
  foreach my $f (@f){
  	my $cmd = "perl -I../conf detect_common_introns.pl $f  $maxSlide $maxGap $simRegion $minSim $minAbsentSim 0  | grep ins | mysql ciwog -u ciwogUser";
	print "$cmd\n";
        `$cmd`;
  }

#set intronClass.
print `perl -I../conf U12_Predict.pl | mysql ciwog -u ciwogUser`;


  
