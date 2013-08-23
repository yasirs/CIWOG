#!/usr/bin/perl




my $SEQ_START=0; #relative to the donor site
my $SEQ_END=0;   #relative to the acceptor site
my @BPS_REG= (-40, -5); #relative to the acceptor site

#Weight Matrixes from Zhu, 2003.
my %naCode=('A'=>0,'a'=>0,'T'=>1,'t'=>1,'G'=>2,'g'=>2,'C'=>3,'c'=>3);
my $don_mat = 
#   A         T            G        C         Position
[
  [0.000000, 0.000000, 0.000000, 0.000000],   #  1
  [0.000000, 0.000000, 0.000000, 0.000000],
  [0.501196, -3.128329, -2.627466, -1.296697],
  [-4.820028, 1.832309, -1.395475, -2.880061],
  [-3.411230, -3.379540, -4.713417, 3.376372],
  [-3.593073, -4.734689, -2.492875, 2.673939],
  [-2.981455, 1.195474, -2.623661, -3.183019],
  [-3.796760, 1.031056, -1.306239, -1.827300],
  [-0.562288, 0.804143, -2.239227, -1.611923]
];

my $u12_acc_mat = 
[
  [-2.893089, -0.342887, -4.115485, -5.700422],
  [-4.115485, -2.530514, -5.700422, -0.415037],
  [-4.700422, -3.700441, -5.700422, -0.208586],
  [-5.700422, -0.085729, -5.700422, -5.700422],
  [-5.700422, -0.115477, -4.700422, -5.700422],
  [-0.842460, -4.700422, -1.452511, -4.700422],
  [-0.176878, -5.700422, -3.700441, -5.700422],  # BPS
  [-5.700422, -2.378509, -5.700422, -0.378511]
];


my $u2_acc_mat = 

[
  [-2.208693, -1.072045, -2.585847, -2.821525],
  [-2.245258, -1.063170, -2.570625, -2.814278],
  [-2.300486, -1.022225, -2.607229, -2.831000],
  [-2.306981, -1.006985, -2.620314, -2.860403],
  [-2.338653, -0.987985, -2.639782, -2.860702],
  [-2.345427, -0.981598, -2.640294, -2.873901],
  [-2.340318, -1.006489, -2.561132, -2.885397],  #BPS
  [-2.255134, -1.059781, -2.535165, -2.854148]
];




sub predict{
	my ($title,$seq) = @_;
	$don_line = substr($seq, 0, 9);
	$acc_line = substr($seq, $BPS_REG[0]-$SEQ_END, $BPS_REG[1]-$BPS_REG[0]+2);
	@d = split(//, $don_line);
	@a = split(//, $acc_line);
		
	($pd) = &mmScore($don_mat, 0, 1, @d);
		
	($pb_u12,$distba) = &mmScore($u12_acc_mat, 0, scalar(@a)-8+1,@a);
	($pb_u2) =  &mmScore($u2_acc_mat, 0, scalar(@a)-8+1,@a);
		
	$pb = $pb_u12 - $pb_u2;
	$distba = -$BPS_REG[0]-$distba-6;	
		
	return($pd,$pb,$distba);
}




sub mmScore{
	my ($mm, $beg,$len,@seq) = @_;
	
	my $max_s=-999999;
	my $s;
	my $max_p = -1;
	for(my $c=0;$c<$len;$c++){
		$s = 0;
		for(my $i=0;$i<scalar(@$mm);$i++){
			$s+=$mm->[$i][$naCode{$seq[$i+$beg+$c]}] if(defined($naCode{$seq[$i+$beg+$c]}));
		}
		if($s>$max_s){
			$max_s = $s;
			$max_p = $beg+$c;
		}
	}
	return ($max_s,$max_p);
}




require 'conf.pl';

my $sql = "select ClusterUID,annotationUID,org,intronNum,intronSeq from ClusterIntron2Gene where intronSeq != '' ";
my $ref = $NRIdbh->selectall_arrayref($sql);

my @dat = ();


###
#calc log odds scores
###

print "update ClusterIntron2Gene set intronType = 'U2' where intronType = 'U2';\n";
print "update ClusterIntron2Gene set intronType = 'GC' where substring(intronSeq,1,2) = 'GC' ;\n";
for (my $i=0;$i<scalar(@$ref);$i++){
	my @t = predict("",$ref->[$i]->[4]);

	$dat[$i] = \@t;
	if ($t[0] > $u12_min_donor & $t[1] > $u12_min_branch){
		print "update ClusterIntron2Gene set intronType = 'U12' where annotationUID = '$ref->[$i]->[1]' and org = '$ref->[$i]->[2]' and intronNum = '$ref->[$i]->[3]';\n";
	}

}
