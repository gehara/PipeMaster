#!/usr/bin/perl

$usage="
Reports the missing data attributes of a fasta file

missing_ms_report.pl

in=<file>
\t the input file with fasta sequences

";

if($#ARGV < 0){
    die $usage;
}


while($args = shift @ARGV){
    if($args =~ /in=(.*)/i){ $infile = $1; }
    else{ die "Argument $args is invalid. $usage\n"; }
}


open(IN, $infile) or die "Couldn't open $infile for reading\n";
@a = grep{!/^\s*$/}<IN>;

foreach (@a){
    $_ =~ s/\cM\n/\n/g;
}

chomp(@a);
close IN;

$gap='-';
$missing='N';

$readseq = 0;
$seqseg=();

@SEQ=();
@NAMES=();
for(my $j = 0; $j < @a; $j++){
    $line = $a[$j];
    #$line = uc($line);
    if($line !~ /^\s*>/ && $readseq == 0 && $j != $#a){
	next;
    }
    
    # beginning of the sequence
    if($line =~ /^\s*>(.*)/ && $readseq == 0 && $j != $#a){
	$seqname = $1;
	$readseq = 1;
	next;
    }
    if($line !~ /^\s*>/ && $readseq == 1 ){
	$line =~ s/\s+//;
	$seqseg .= $line;
    }
    if( ( ($line =~ /^\s*>/) && ($readseq == 1) && ($j != $#a) ) 
	|| ($j == $#a) 
	){
	push @SEQ, $seqseg;
	push @NAMES, $seqname;
	$readseq = 0;
	$seqseg = ();
	if( $j != $#a ){ $j--;}
    }
}

#check for the length of the sequences
$seq_length = 0;
$change = 0;
for my $seq (@SEQ){
    if( length($seq) > $seq_length ){ 
	$seq_length = length($seq);
	$change++;
	}
    if($change > 1){
	print STDERR "ERROR: The length of the sequence $keys is different".
	    " than the length of the rest ones. This is not allowed\n";
	for my $ind (0..$#SEQ){
	    print STDERR $NAMES[$ind], "\t", length($SEQ[$ind]), "\n";
	}
	exit;
    }
}


print STDERR "Length: $seq_length\n";


for my $ind (0..$#SEQ){
    my $seq = $SEQ[$ind];
    for(my $j=1; $j<=$seq_length; $j++){ # attention
	$c = substr($seq, $j-1, 1); # !! attention
	if( ($c eq $gap) || ($c eq $missing) ){
	    print $ind, "\t", ($j-1)/$seq_length, "\t", ($j+1)/$seq_length, "\n";
	    
	}
    }
}

