#!/usr/bin/perl

use strict;
use Data::Dumper;
use JSON;
use CGI qw(:standard);

my $TABIX_PATH= '/usr/bin/tabix';
my $TABINDEXER = './tabindexer.pl';
my $DATA_DIR = '../data/';
my $INT_FILE = "$DATA_DIR/gene_targets.tab";
my $INT_FILE_IDX = "$INT_FILE.idx";
my $TI_PATH = "perl $TABINDEXER -t $INT_FILE -i $INT_FILE_IDX -s";
my $GENE_FILE = "$DATA_DIR/Homo_sapiens.GRCh37.75.genes.gtf.gz";
#my $ASSOC_FILE = "$DATA_DIR/hg19_gwas_ic_t1d_onengut_meta_4_18_1.gff.gz";
my $ASSOC_FILE = "$DATA_DIR/gb2_hg19_gwas_t1d_barrett_4_17_0.gff.gz";
my @GFF_FIELDS = ('chr','type','source','start','end','score','strand','phase','attributes');
#my $GENE2ENSG = "perl $TABINDEXER -t $DATA_DIR/ens2gene.tab -i $DATA_DIR/ens2gene.tab.idx -s";
my $GENE2ENSG = "perl $TABINDEXER -t $DATA_DIR/bid2gene.tab -i $DATA_DIR/bid2gene.tab.idx -s";
my $MAX_REGION_SIZE = 2e6;


my $cgi = CGI->new();

my %gene_colors = (
	'miscRNA' => 'ff7f0e',
	'protein_coding' => '1f77b4',
	'lincRNA' => '2ca02c',
	'snoRNA' => 'd62728',
	'antisense' => '9467bd',
	'miRNA' => 'e377c2',
	'snRNA' => '8c564b',
	'pseudogene' => '7f7f7f',
	'processed_transcript' => 'bcbd22'
);

if ($cgi->param('region')){
	my $region = $cgi->param('region');
	my ($c,$s,$e) = $region =~ /^(\d+):(\d+)-(\d+)$/;
	my $ret = ();
	$ret->{genes} = getGenes($c,$s,$e);
	$ret->{snps} = getSNPs($c,$s,$e);
	print header('application/json');
	print(encode_json($ret));
	
	exit(0);	
}


my $gene = $cgi->param('gene');
my $ensg;
## convert to ensg

my $cmd = "$GENE2ENSG $gene";
#print STDERR $cmd;

open(GENE,"$cmd |") || die "Cannot convert\n";
while(<GENE>){
	chomp;
	my @vals = split("\t",$_);
	$ensg = $vals[0];
	last;
}
print STDERR "ENS = ".$ensg."\n";
close(GENE);

my $tissue = $cgi->param('tissue');

sub makeRelative{
	my ($start,$end,$coord_fields,$ds) = @_;
	foreach my $k(@$coord_fields){
		#print STDERR join("--",$k,$ds->{$k},$start)."\n";
		my $coord = $ds->{$k} - $start + 1;
		#print STDERR join("--",$k,$ds->{$k},$start,$coord)."\n";
		$ds->{$k} = $coord > 0 ? ($coord > ($end-$start) ? $end-$start+1 : $coord) : 1;
	}
	#print "Next\n";
	return $ds;
}


sub getHiC{
	my ($gname,$tissue,$maxDist)=@_;
	my $cmd = "$TI_PATH $gname";
	#print STDERR $cmd;
	open(DATA,"$cmd |") or die "Error in tabindexer.pl request: $cmd\n";
	my @header=qw/cell bid chr start end gene type bait/;
	my %baits;
	my @targets;
	my @coords;
	my $ochr;
	my $bCoord;
	while(<DATA>){
		chomp;
		my ($t,$bid,$chr,$start,$stop,$gene,$type,$bait) = split("\t",$_);
		push @coords,($start,$stop);
		next unless $t eq $tissue;
		if($bait){
      #print "In here\n";
			$baits{'#bID'}=$bid;
			$baits{'bEnd'}=$stop;
			$baits{'bSt'}=$start;
			$bCoord=(($stop-$start)/2) + $start;
		}else{
			push @targets,{
				'oeEnd'=>$stop,
				'oeSt'=>$start
			};
		}
		$ochr=$chr;
	}
	my @results;
	my @so = grep{abs($_-$bCoord)<$maxDist}sort {$a<=>$b} @coords;
	foreach my $t(@targets){
		my %record = (%baits,%$t);
		my $tCoord = (($record{oeEnd}-$record{oeSt})/2) + $record{oeSt};
		
		if(abs($tCoord-$bCoord)<$maxDist){
			#print STDERR $so[-1]."\n";
			push @results, makeRelative($so[0],$so[-1],[qw/oeSt oeEnd bSt bEnd/],\%record);
		}
		#push @results, \%record;
	}
	
	
	
	my $region="$ochr:$so[0]-$so[$#so]";
	return({
			meta => { 
				ostart => $so[0],
				oend => $so[$#so],
				rstart => 1,
				rend => $so[$#so] - $so[0],
				rchr => $ochr
			},
			hic => \@results,
			region=>$region
	})
}


sub getGeneTabixData{
	my ($region,$datafile,$mand_list)=@_;
	my @field_list = @$mand_list;
	#print Dumper(\@field_list)."\n";
	my $cmd = join(" ",$TABIX_PATH,'-h',$datafile,$region);
	open(DATA,"$cmd |") or die "Error in tabix request $cmd\n";
	#print STDERR "$cmd\n";
	my @results;
	while(<DATA>){
		chomp;
		next if /^#/;
		my @tmp = split("\t",$_);
		my %res;
		for(my $i=0;$i<@field_list;$i++){
			## cope with attributes field
			if($field_list[$i] eq 'attributes'){
				($res{'gene_name'} = $tmp[$i]) =~s/.*gene_name\s"([^"]+)".*/\1/;
				($res{'gene_id'} = $tmp[$i]) =~s/.*gene_id\s"([^"]+)".*/\1/;
				$res{'bumpLevel'} = 0;
			}else{
				if($i == 1){
					my $t = ($tmp[$i] =~ /pseudo/) ? 'pseudogene' : $tmp[$i];
					my $c = (defined $gene_colors{$t}) ? $gene_colors{$t} : '17becf';
					$res{color} = $c;
				}
				$res{$field_list[$i]}=$tmp[$i];
			}
		}
		push @results,\%res;
  }
	return \@results;
}

sub getSNPTabixData{
	my ($region,$datafile,$mand_list)=@_;
	my @field_list = @$mand_list;
	#print Dumper(\@field_list)."\n";
	my $cmd = join(" ",$TABIX_PATH,'-h',$datafile,$region);
	open(DATA,"$cmd |") or die "Error in tabix request $cmd\n";
	#print STDERR "$cmd\n";
	my @results;
	while(<DATA>){
		chomp;
		next if /^#/;
		my @tmp = split("\t",$_);
		my %res;
		for(my $i=0;$i<@field_list;$i++){
			## cope with attributes field
			if($field_list[$i] eq 'attributes'){
				($res{'snp_name'} = $tmp[$i]) =~s/.*Name=([^;]+);.*/\1/;
				$res{'bumpLevel'} = 0;
			}else{
				$res{$field_list[$i]}=$tmp[$i];
			}
		}
		push @results,\%res;
  }
	return \@results;
}

sub getGenes{
	my ($chr,$start,$end) = @_;
	my $region = "$chr:$start-$end";
	my $genes = getGeneTabixData($region,$GENE_FILE,\@GFF_FIELDS);
	for(my $i=0;$i<@$genes;$i++){
		my $tmp = makeRelative($start,$end,['start','end'],$genes->[$i]);
		$genes->[$i]=$tmp;
	}

	return ($genes);
}

sub getSNPs{
	my ($chr,$start,$end) = @_;
	my $region = "$chr:$start-$end";
	my $snps = getSNPTabixData($region,$ASSOC_FILE,\@GFF_FIELDS);
	for(my $i=0;$i<@$snps;$i++){
		my $tmp = makeRelative($start,$end,['start','end'],$snps->[$i]);
		$snps->[$i]=$tmp;
	}

	return ($snps);
}


my $ret = getHiC($ensg,$tissue,2e6);
$ret->{genes} = getGenes($ret->{meta}->{rchr},$ret->{meta}->{ostart},$ret->{meta}->{oend});
$ret->{snps} = getSNPs($ret->{meta}->{rchr},$ret->{meta}->{ostart},$ret->{meta}->{oend});
print header('application/json');
print(encode_json($ret));

