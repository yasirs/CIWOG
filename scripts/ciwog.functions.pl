#!/usr/bin/perl

use strict vars;
use vars qw{ %codon};

###############
# This file contains functions for ciwog
###############


sub readDirect(){
  my $directory = shift;
  my $opathOut = shift;
  opendir(DIR,$directory) || die "cannot open directory $directory";
  my @f = grep /\w/,readdir(DIR);
  closedir DIR;
  
  foreach my $f (@f){
  	my $pathOut = "$opathOut/$f";
  	system("rm -r $pathOut");
  	system("mkdir $pathOut");
  	$pathOut = "$pathOut/";
  	parseTextFile($f,$directory,$pathOut);
  }
}





sub parseTextFile{
	###
	# this is the functions script for parsing intron sequence, proteins, and intron positions
	###

	my $fileName = shift; #location of intronFile
	my $dirp=shift;
	my $pathOut = shift;
	print STDERR "processing $fileName\n";
	#read elements by file
	if ( ! -e "$dirp/$fileName/$fileName.aln" ){
		#skip if cluster/cluster does not exist.
		next;
	}

	open(IF, "$dirp/$fileName/$fileName") || die "cannot open $fileName";
	my $ind = -1;
	my @fields = ("organism","geneId","chromosome","geneStructure","cdsStart","cdsStop","genomeStart","genomeStop","genomeSequence");
	my $dat={};
	my $gc=1;
	while(<IF>){
		if (/^#/){
			if($ind==8){
				$dat->{$gc}->{genomeSequence} =~ s/\n//g; #remove line breaks if any
				
				#reformat gene ID HACK for Yuanbin's line
				($dat->{$gc}->{geneId}) = $dat->{$gc}->{geneId} =~ /([^\|]+?)$/;
				
				#new gene
				$ind=-1;
				$gc++;
			}
						
			$ind++;
			chomp;
			my ($tmp) = $_ =~ /#(.+)/;
			if($fields[$ind] eq "genomeSequence" && $tmp =~ /^>/){
				next;
			}
			$dat->{$gc}->{$fields[$ind]}=$tmp;
		}else{
			$dat->{$gc}->{$fields[$ind]}.=$_;
		}
		
	}
	#reformat genomic sequence
	$dat->{$gc}->{genomeSequence} =~ s/\n//g; #remove line breaks if any	
	#reformat gene ID HACK for Yuanbin's line
	($dat->{$gc}->{geneId}) = $dat->{$gc}->{geneId} =~ /([^\|]+?)$/;
	
	close IF;
	
	foreach my $gc (keys(%$dat)){
		#print "$gc $dat->{$gc}->{geneId} \n";
		get_protein_seq_w_introns($dat->{$gc}->{geneId},$dat->{$gc}->{geneStructure},$dat->{$gc}->{cdsStart},$dat->{$gc}->{cdsStop},$dat->{$gc}->{genomeSequence},$dat->{$gc}->{genomeStart},$dat->{$gc}->{genomeStop},$pathOut,$dat->{$gc}->{chromosome}, $dat->{$gc}->{organism} );
		
	}
	
	#print "$dat->{geneStructure},$dat->{cdsStart},$dat->{cdsStop},$dat->{genomeSequence},$dat->{genomeStart},$dat->{genomeStop},$pathOut ";

}



sub get_protein_seq_w_introns {

  	## need to test exons < 3
  	my ($geneId,$gs,$cdsStart,$cdsStop,$genomeSequence,$genomeStart,$genomeStop,$pathOut,$chromosome,$organism) = @_;

    
  	#get boundaries of $cds;
	my ($l,$r) = $gs =~ /^\D+?(\d+).+?(\d+)\D+?$/;
	#trim genomic sequence
		if($genomeStart < $l){
			my $lpad = $l - $genomeStart; 
			$genomeSequence = substr($genomeSequence,$lpad);
		}
		if($genomeStop > $r){
			my $rpad = $genomeStop - $r;
			$genomeSequence = substr($genomeSequence,0,length($genomeSequence)-$rpad);
		}
	
	if ($gs =~ /comp/){
		#reverse complement if needed
		$genomeSequence =~ tr/a-z/A-Z/; 
		$genomeSequence = reverse $genomeSequence;
      	$genomeSequence =~ tr/ACGT/TGCA/;
	}
	
	#extract introns from gene structure
	my @introns;
	my @if = $gs =~ /(\d+?\,\d+)/g;
	for(my $i=0;$i<scalar(@if);$i++){
		my @c = split /\,/, $if[$i];
		$c[0]++; #intron coord
		$c[1]--; #intron coord
		
		my $intronStart = $c[0]-$l; #relative to beginning of sequence reading direction start
		my $intronLength = $c[1]-$c[0]+1;
		
		if($gs =~ /comp/){
			$intronStart = $r - $c[1];
		}
		#print "$c[0],$c[1], $intronStart, $intronLength";
		$introns[$i] = [$i,$c[0],$c[1], $intronStart, $intronLength];
	}
	
	
	my $pcSeq = "";
	#scan to first intron
	
	my $cdsPad = ($gs =~/comp/) ? $r-$cdsStart : $cdsStart - $l;
	
	

	my $cdsLast = ($gs =~/comp/) ? $r-$cdsStop : $cdsStop - $l;
	
	my $lastgp=$cdsPad;
	
	if ($gs =~ /comp/){
		@introns = reverse @introns;
		for (my $i=0;$i<scalar(@introns);$i++){
			$introns[$i][0] = $i;
		}
	}
	my $lastCode=-1;
	
	#Function proceeds by intron, adding protein coding sequence
	
	#write introns to file
	open(intF,">>${pathOut}introns.csv") || die "cannot open introns.tsv";
	for(my $i=0;$i<scalar(@introns);$i++){
		
		my $pcNt = substr($genomeSequence,$lastgp,$introns[$i][3] - $lastgp );
		#print "$pcSeq\n";
		$lastgp = $introns[$i][3] + $introns[$i][4]  ; #intronStart + intronLength + 1	

		my ($phase,$proteinPos);
		#5' utr intron.	
		if(
			($introns[$i][3] <= $cdsPad)
		){
			#could add print stat
			$pcNt=$phase="";
			$proteinPos = "5utr";
			$lastgp=$cdsPad; #keep as original default.
		}
		
		#3' utr intron
		if(
			($introns[$i][3] > $cdsLast)
		){
			#could add print stat
			$pcNt=$phase="";
			$proteinPos = "3utr";
			#$lastgp=$cdsPad; #original.
		}
		$lastCode = ($introns[$i][3] > $cdsLast && $lastCode==-1) ? $i : $lastCode; #last coding intron number.  
		$pcSeq .= $pcNt;
		if(!$proteinPos){
			$phase = length($pcSeq) % 3;
			$proteinPos = int(length($pcSeq)/3); #proteinPos is number of prior complete codons.
		}
		
		#print "$introns[$i]->[1],$introns[$i]->[2],$introns[$i]->[3],$introns[$i]->[4]\n";
		my $intronSeq = substr($genomeSequence,$introns[$i][3],$introns[$i][4]);
		my $strand = ($gs =~ /comp/) ? "-" : "+";
		$introns[$i] = [@{$introns[$i]},$intronSeq,$phase,$proteinPos,$geneId,$chromosome,$organism,$strand];
		print intF join(",",@{$introns[$i]})."\n";
		
	}
	close intF;

	##
	#add last exon peptide contribution
	##
		my ($lastExonStart);
		if(scalar(@introns)>0){
			if ($lastCode != -1 ){ # -1 index for last intron.
				$lastCode = $lastCode-1; #last PC intron before UTR
			}
			$lastExonStart = $introns[$lastCode][3] + $introns[$lastCode][4] ;
			$lastExonStart = ($lastExonStart < $cdsPad) ? $cdsPad : $lastExonStart;		# CDS completely within last exon.

			#CDS within first exon
			$lastExonStart = ($lastExonStart > $cdsLast) ? $cdsPad : $lastExonStart;

		}else{
			$lastExonStart = $cdsPad; #single exon genes.
		}

		my $lastNT = substr($genomeSequence,$lastExonStart, $cdsLast - $lastExonStart+1 );
		$pcSeq .= $lastNT;
		#print "$lastCode $lastExonStart $cdsLast \n";	

	
	
	#calc protein
	my $pep="";
	for(my $i=0;$i<length($pcSeq);$i=$i+3){
		$pep .= $codon{substr($pcSeq,$i,3)};
	}

	#write ProteinSequence to file
	open(proteinF,">>${pathOut}protein.fsa") || die "cannot open protein.fsa";
	print proteinF ">$organism\|$geneId\n";
	print proteinF "$pep\n";
	#print ">$organism\|$geneId\n$pep\n";
	close proteinF;


	####
	# integrity check for proper protein sequence parsing..
	####
	if ($pep =~ /\W.*?.$/){
		print STDERR "\ninternal stop in $organism $geneId\n";
	}

	if ($cdsLast < $lastExonStart){
		print STDERR "\npeptide parsing error $organism $geneId\n";
	}

		
	###
	# write sql insert statements.
	###
	
}	



	
sub min{
  my ($a,$b) = @_;
  if ($a eq ""){
    return $b;
  }
  if ($b eq ""){
    return $a;
  }
  if ($a >$b){
    return $b;
  }else{
    return $a;
  }
}

sub max{
  my ($a,$b) = @_;
  if ($a eq ""){
    return $b;
  }
  if ($b eq ""){
    return $a;
  }

  if ($a < $b){
    return $b;
  }else{
    return $a;
  }
}


sub revcomp{
    my $f = shift;
    $f = reverse $f;
    $f =~ tr/ATCGatcg/TAGCtagc/;
    return $f;
}


sub parseId{
	my $str = shift;
	my ($org,$id,$str) = $str =~ /^(\S+?)\|.+?([^\|]+?)\s+(\S+)/;
	#my ($org,$id,$str) = $str =~ /^(\S+?)\|(\S+?)\s+?(\S+)/;
	return($org,$id,$str);
}

	
###########
#variables
###########

## standard vertebrate codon table
 %codon = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    '---' => '-',    # In-frame gap
    '...' => '_',    # In-frame gap
    'NNN' => 'N',	 # UNK
    '???' => '?',    # UNK
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );


#protein sequence
#proteinPos
#phase
#intron sequences

#my ($proteinIseq,$protein_seq,$iPos,$iPhase) = get_protein_seq_w_introns($cds,$cdsStart,$cdsEnd,$chrom,$org,0);






1;
