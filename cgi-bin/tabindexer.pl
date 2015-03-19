#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

### FUNCTIONS

sub create_index{
	my ($tfile,$index,$column)=@_;
	my $head=0;
	my $tail=0;
  print "COlumn is $column\n";
	open(IN,$tfile);
	open(OUT,">$index");
	my %RES;
	while(<IN>){
    chomp;
		my @FIELDS = split("\t",$_);
		$tail = tell IN;
		push @{$RES{$FIELDS[$column]},},[$head, $tail-$head];
		$head = $tail;
	}
	close(IN);

	foreach my $gene(sort keys %RES){
		my @locs = map{join(",",@$_)}sort{$a->[0] <=> $b->[0]}@{$RES{$gene}};
		print OUT join("\t",$gene,@locs)."\n";
	}
	close(OUT);
}

sub get_lines{
	my ($tfile,$index,$search)=@_;
	open(DATA,"$tfile") || die "Cannot open $tfile\n";
	open(IDX,"look -b $search $index|");

	while(<IDX>){
		chomp;
		my @idf = split("\t",$_);
		shift @idf;
		foreach my $loc(@idf){
			my ($offset,$length)=split(",",$loc);
			seek DATA, $offset,0;
			my $line ='';
			read DATA, $line, $length;
			chomp($line);
			print $line . "\n" if ($line =~ /\b$search\b/);
			#print STDERR $line . "\n" if ($line =~ /\b$search\b/);
		}
	}
	close(IDX);
	close(DATA);
}


my $USAGE=<<EOL;
perl index.tab -[t]file -[c]olumn { -[s]earchterm -[f]orce -[i]ndex -[h]elp)

Simple tab delimited file indexer

	tfile - tab delimited file to index
	column - column to create index (zero based)
	index - index file to create or use
	force - force overwrite of given index
	searchterm - search string
	help - print this message

For help and suggestions please contact Olly Burren 
EOL

my ($tfile,$index,$column,$help,$force,$search,$ERROR_FLAG);

GetOptions(
	'tfile|t=s' => \$tfile,
	'help|h' => \$help,
	'column|c=i' => \$column,
	'searchterm|s=s' => \$search,
	'force|f'=> \$force,
	'index|i=s' => \$index);

if($help){
	print STDERR $USAGE."\n";
	exit(1);
}

if(! -e $tfile){
	print STDERR "Cannot locate $tfile\n";
	$ERROR_FLAG++;
}

if(! -e $index && !$column ){
	print STDERR "Index not found and cannot create as column is missing\n";
	$ERROR_FLAG++;
}

if(-e $index && $search){
	## this is a seek operation
	get_lines($tfile,$index,$search);
	exit(1);
}

if(-e $index && ! $force){
	print STDERR "Index file $index already exists use --force to recreate\n";
	exit(1);
}
print STDERR "CREATING INDEX\n";
unless($column>=0){
  print STDERR "Column not set aborting\n";
  exit(1);
}


create_index($tfile,$index,$column);
print STDERR "INDEX CREATED IN $index\n";


