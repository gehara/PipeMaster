#!/usr/bin/perl -w

$usage = "

The program transforms the nex or fasta file format  to ms-like format.
It can use an outgroup if available. Then, derived mutations are marked with 1 and ancestral
state with 0.
If there is no outgroup, then the less frequent state is marked with 1 and the other with 0

=============================================================================================

usage:

fas2ms.pl fas=<file>  (this is suitable for fasta files)


--------------------------------------------------------------------------------------------

OPTIONAL

outgroup:\ta string that denotes the outgroup. It should be in the nex or fasta file file, and should be
uniq among the names of the sequences. 
For example if the outgroup is denoted as SIM_01, the argument could be outgroup=SIM, if noother sequence
contains the name SIM

==============================================================================================

";
$cnt = 0;
if($#ARGV < 0 ){die $usage;}

$nexusfile = "";
$infofile = "";
# $startposition = 1;
$warning = 0;

$folded = 0;

################### haplotypes options #############
$outformat = "";
$max_length = 0;
$invarout = 1;
###################################################
$lengthout = 0;


$outformat="hap";
$postype=0;
$fixed = 0;
$FASTA = 0;
$header = 0;
$NOREC = 0;

####################### filter outgroup ##############
# it filters out the sites that the outgroup shows some state and
# all the sequences in the populations show "-"

$filteroutgroup = 0;

######################################################
$outseq="outgroup";
$groups="";

$printmiss = 0;
while($args = shift @ARGV){
    if($args =~ /^nex=(.*)/i){$nexusfile = $1; $informat = "nexus";}
    elsif($args =~ /^fas=(.*)/i){$nexusfile = $1; $informat = "fasta";}
    elsif($args =~ /^info=(.*)/i){$infofile = $1;}
    elsif($args =~ /^outgroup=(.*)/i){$outgroup = $1; print STDERR "# $outgroup\n";}
    elsif($args =~ /^-FASTA/i){ $FASTA = 1; }
    elsif($args =~ /^-header/i){ $header = 1;}
    elsif($args =~ /^outformat=(.*)/i){$outformat = $1;}
    elsif($args =~ /^maxlength=(.*)/i){$max_length = $1;}
    elsif($args =~ /^-nomut/i){ $invarout = 0;}
    elsif($args =~ /^-fixed/i) { $fixed = 1 ; }
    elsif($args =~ /^TOREM=(.*)/i){ $TOR = $1;}
    #elsif($args =~ /^POS=(.*)/i){ $postype = $1; }
    #elsif($args =~ /^-LENGTH/i){ $lengthout = 1;}
    elsif($args =~ /^-NOREC/i){ $NOREC = 1;}
    elsif( $args =~ /^-filteroutgroup/i){ $filteroutgroup = 1; }
    elsif($args =~ /^groups=(.*)/i){ $groups = $1; }
    elsif($args =~ /^-printmiss/i){ $printmiss = 1; }
    else{die "Argument $args is invalid \n $usage";}
    
}

if($nexusfile eq ""){die "Nexus file is needed for the program to run\n";}


@toREMOVE = ();
@toREMOVE = split(/,/, uc($TOR) );
print STDERR "**** @toREMOVE ****\n";


# open the files
open(NEX, $nexusfile) or die "Coudln't open the nexus file $nexusfile for reading\n";
@A = grep{!/^\s*$/}<NEX>;

# get rid of the ^M characters
foreach (@A){
    $_ =~ s/\cM\n/\n/g;
}
chomp(@A);
close NEX;


if($nexusfile =~ /.*\/(.+)$/){$nexusfile = $1;}


# open the infofile and store the information in a hash
if($infofile ne ""){
    open(INFO, $infofile) or die "Couldn't open the info file for reading\n";
    @B = grep{!/^\s*$/}<INFO>;
    foreach(@B){
	$_ =~ s/\cM\n/\n/g;
    }
    chomp(@B);
    close INFO;
}


if($infofile ne ""){
    map{
	my @V = split /\s+/ ;
	$SPOSITION{$V[1]} = $V[0];
	if($outgroup eq "FILE"){
	  $OUTG_NAME{$V[1]} = $V[2];
	}
    }@B;
}

#parse the nexus file

# remove everything until 'matrix', included
if($informat =~ /nexus/i){
    $gap = 0;
    $missing = 0;
    
    while(my $line = shift @A){
	$cnt++;
	# the gap character
	if($line =~ /gap=(.)/i){$gap = $1;}
	if($line =~ /missing=\'(.)\'/i || $line =~ /missing=(.)/i ){$missing = $1;}
	if($line !~ /^\s*matrix\s*$/i){
	    next;
	}
	else{last;}
    }
    
    if(! $gap){
	print STDERR "warning:\t didn't find gap character information\n";
	$warning ++;
    }
    
    if(! $missing){
	print STDERR "warning:\t didn't find missing character information\n";
	$warning ++;
    }
    
    %SEQ = ();
    
    while(my $line = shift @A){
	$cnt++;
	if(($line =~ /^\s*\[.*\]/) || 
	   ($line =~ /\s*;/) || 
	   ($line =~ /end/) ){ next; }
	elsif($line =~ /([^\s]+)\s+(.*)/){ # ATTENTION
	    $name = $1; 
	    $name = uc( $name ); # transform to UPPERCASE, mostly for sorting reasons
	    $seqseg = $2;
	    $seqseg =~ s/\s+//g;
	}
	else{ 
	    print STDERR " warning: the line $line \n in  $cnt contains unacceptable characters\n";
	    $warning++;
	}
	
	# append the sequence into the hash
	
	# put the first segment of the sequence
	if(! $SEQ{$name}){ 
	    $SEQ{$name} = $seqseg;
	} 
	else{
	    $SEQ{$name} .= $seqseg;
	}
    }
}
elsif($informat =~ /fasta/i){
    $gap = '-';
    $missing = 'N';
    $readseq = 0;
    $seqseg = ();
    for(my $j = 0; $j < @A; $j++){
	$line = $A[$j];
	#$line = uc($line);
	if($line !~ /^\s*>/ && $readseq == 0 && $j != $#A){
	    next;
	}

	if($line =~ /^\s*>(.*)/ && $readseq == 0 && $j != $#A){
	    $seqname = $1;
	    $readseq = 1;
	    next;
	}
	if($line !~ /^\s*>/ && $readseq == 1 ){
	    $line =~ s/\s+//;
	    $seqseg .= $line;
	}
	if( ($line =~ /^\s*>/ && $readseq == 1 && $j != $#A) || $j == $#A){
	    $SEQ{$seqname} = $seqseg;
	    $readseq = 0;
	    $seqseg = ();
	    if( $j != $#A ){ $j--;}
	}
    }
}
	    

# remove sequences you wish
foreach my $rem (@toREMOVE){
    print STDERR $rem, "\n"; #P
    @REM = grep {/$rem/i} keys %SEQ;
    print STDERR "@REM\n"; #P
    foreach my $torem(@REM){
	print STDERR $torem, "\n"; #P
	delete $SEQ{$torem};
	print STDERR "Attention !!! Sequence: $torem is REMOVED !\n";
    }
}
		
print STDERR "Eventually, there are ".keys (%SEQ)." sequences\n";

#check for the length of the sequences
$seq_length = 0;
$change = 0;
for my $keys(keys %SEQ){
    if( length($SEQ{$keys}) > $seq_length ){ 
	$seq_length = length($SEQ{$keys});
	$change++;
	}
    if($change > 1){
	print STDERR "ERROR: The length of the sequence $keys is different".
	    " than the length of the rest ones. This is not allowed\n";
	for my $keys(keys %SEQ){
	    print STDERR $keys, "\t", length($SEQ{$keys}), "\n";
	}
	exit;
    }
}

if($lengthout){
    print $nexusfile, "\t", $seq_length, "\n";
    exit;
}
print STDERR "Length: $seq_length\n";



####################################   check for the outgroup
if($outgroup eq "FILE"){
    $outgroup = $OUTG_NAME{$nexusfile};
    # print STDERR $outgroup, "\n";
    # exit;
}

$nofoutgroup = 0;
%OUTG = ();
foreach my $keys(keys %SEQ){
    print STDERR $nofoutgroup, "\t", $outgroup, "\t", $keys, "\n";
    if($keys =~ /$outgroup/i){
	$nofoutgroup ++;
	print STDERR "OUTGROUP SEQUENCE $keys is found (match: $keys-$outgroup)\n";
	
	# put the outgroup in a new hash
	$OUTG{$keys} = $SEQ{$keys};
	# remove the outgroup from the hash
	delete $SEQ{$keys};
    }
}

if($nofoutgroup == 0){
    print STDERR "# No outgroup found\nThat means that the results are folded\n";
    $folded = 1;
}

if($nofoutgroup > 1){
    print STDERR "Warning: more than outgroup found -- TODO\n";

}


    
    
# transform the string into arrays
for my $key (keys %SEQ){
    @v = split(//,$SEQ{$key});
    $SEQ{$key} = [@v];
}

for my $key (keys %OUTG){
    @v = split(//, $OUTG{$key});
    $OUTG{$key} = [@v];
}

#decide about a single "average" outgroup and report the problematic sites
# the unified outgroup
@fout=();
# for every site
for(my $i=0; $i<$seq_length; $i++){
    %outinfo=();
    $noutinfo=0;
    foreach my $key (keys %OUTG){
	$site=$OUTG{$key}[$i];
	#print STDERR $site, "\n";
	if( ($site ne "-") && !(defined($outinfo{$site})) ){
	    $noutinfo++;
	    $outinfo{$site}=$site;
	}
    }
    # print STDERR "noutinfo: $noutinfo\t";
#     foreach my $key( keys %outinfo ){
# 	print $outinfo{$key}, "\n";
#     }

    # if more than one nucleotide is provided by the outgroup
    if($noutinfo > 1){
	%sinfo=();
	$nsinfo=0;
	foreach my $skey (keys %SEQ){
	    $ssite = $SEQ{$skey}[$i];
	    if( ($ssite ne "-") && (defined($outinfo{$ssite})) && !(defined($sinfo{$ssite})) ){
		$nsinfo++;
		$sinfo{$ssite} = $ssite;
	    }
	}
	# ambiguisity about the ancestral state
	if($nsinfo > 1){
	    print STDERR "Cannot decide for the site $i: outgroup: ";
	    foreach my $site (keys %outinfo){
		print STDERR "$outinfo{$site},";
	    }
	    print STDERR " sequences: ";
	    foreach my $site (keys %sinfo){
		print STDERR "$sinfo{$site},";
	    }
	    push @fout, "-";
	    print STDERR " Assume outgroup is -\n";
	}
	# only one site is common between the outgroup and the sequences
	elsif($nsinfo == 1){
	    my $l=0;
	    foreach my $site (keys %sinfo){
		push @fout, $sinfo{$site};
		$l++;
	    }
	    if($l > 1){
		print STDERR "ERROR l should be 1: l=$l\n";
	    }
	}
	#outgroup provides irrelevant information
	elsif($nsinfo < 1){
	    push @fout, "-";
	}
    }
    elsif($noutinfo == 1){
	my $l=0;
	$nsinfo=-2;
	foreach my $key (keys %outinfo){
	    push @fout, $outinfo{$key};
	    $l++;
	}
	if($l > 1){
	    print STDERR "ERROR l should be 1: l=$l\n";
	}
    }
    #outgroup provides irrelevant information
    elsif( $noutinfo < 1){
	$nsinfo=-3;
	push @fout, "-";
    }
    #print STDERR "$i\t$fout[$#fout]\t$#fout\t$nsinfo\t$noutinfo\n";
    #exit;
}


#reset the outgroup and build it again
%OUTG=();
$OUTG{$outseq} = [@fout];

if(scalar @fout != $seq_length){
    print STDERR "ERROR! the length of the ougroup is ", scalar @fout, " the length of the sequence is: $seq_length\n";
    exit;
}



# this is needed because later positions will be removed

# just take a key
$akey = "";
foreach my $key(keys %SEQ){
    $akey = $key;
    last;
}
    
#create the hash of positions
for(my $i = 0; $i < $seq_length; $i++){
    $position[$i] = $i;
}

# get rid of the sites that contain a gap but output them in a file

open(GAPS, ">GAPS.OUT") or die "Coudln't open GAPS\n";

# an array to store position of contiguous gaps
@contgaps=();

open(LOG1, ">>OUT.LOG1") or die "Coudln't open LOG1\n";
#print STDERR "gap: $gap\t missing: $missing\n";
$actlength = $seq_length;
#$initial_length = $seq_length;

$contflag = 0;

$recently_spliced_position = 0;
$prevar = 0;
for(my $i = 0; $i < @{$SEQ{$akey}}; $i++){
    $var = 0;
    
    foreach my $key (keys %SEQ){
	if( ($SEQ{$key}[$i] eq $gap) || ($SEQ{$key}[$i] eq $missing) ){
	    $var++;
	    
	    # splice for the sequences
	    # print $i, "\n";
	}
    }
    

    #print STDERR $position[$i], "\t", scalar (keys %SEQ), "\t$var\n";
    
    if(!$printmiss && $var > 0 ){
	$prevar = $var;
	if(!$contflag){
	    $contgaps[0] = $position[$i];
	    #$gapsflag = 1;
	    $contflag = 1;
	}
	foreach my $k (keys %SEQ){
	    splice(@{$SEQ{$k}}, $i, 1);
	}
	# splice for the outgroup
	foreach my $k (keys %OUTG){
	    splice(@{$OUTG{$k}}, $i, 1);
	}
	$actlength--;
	$recently_spliced_position = $position[$i];
	# the position of the gap is spliced out
	splice(@position, $i, 1);
	# exit the loop because you have checked everything
	# in the two small groups
	
	# make the index one less because the position 
	# has been removed
	$i--;
	$seq_length --;
	
    }
    elsif(!$printmiss && $var == 0){
	if($contflag){
	    $contgaps[1] = $recently_spliced_position;
	    if($prevar > 0 && $prevar < scalar(keys %SEQ)){
		print GAPS $contgaps[0], "\t", $contgaps[1], "\n";
	    }
	    $prevar = 0;
	    $contflag = 0;
	}
	
    }
    
    #print STDERR "At position $i we found $var gaps or missing\n";
    if( !$printmiss && ($var == keys (%SEQ)) && $filteroutgroup){
	print LOG1 "******************************  site $position[$i] from file $nexusfile is removed **********************\n";
	# remove this position from the position list
	for(my $j=$i+1; $j < @{ $SEQ{$akey} }; $j++){ # it is $i+1 because of the $i-- in the previous if-section
	    $position[$j] = $position[$j]-1;
	}
    }
    
    


}

print STDERR "-----------------------------------------\n";



# if information is given transform the positions
$abs_position = 1;
if($infofile ne ""){
    $abs_position = $SPOSITION{$nexusfile};
    print STDERR "ABSOLUTE POSITION: $abs_position\n";
}
for (my $i = 0; $i < @position; $i++){
    $position[$i] += $abs_position;
}

    
print STDERR " Final Position is $position[ $#position ]\n";
if($max_length == 0){ $max_length = $position[ $#position ];}

#####################################################################################################################
#########################################################  OUTPUT ###################################################
#####################################################################################################################


###################
## FASTA ##########
###################
if($FASTA){
    if($header){
	print keys(%SEQ)."\t".scalar @position."\t1\n";
    }
    
    foreach my $key (sort keys %SEQ){
	print ">$key\n";
	for(my $i = 0; $i < @{$SEQ{$akey}}; $i++){
	    print $SEQ{$key}[$i];
	}
	print "\n";
    }
    exit;
}
    

#####################
## BINARY STATES #### 
####################
%BinDERIVED=();
@REL_POSITION =();
$segsites=0;

    
    
    
@ANCESTOR = (); # the ancestor state is stored
%DERIVED = (); # the derived states are stored
    
$n = keys (%SEQ); @N = ();  # the size of a hash is given from the reference
@X = ();
$folded = 0;


for(my $i = 0; ($i < $seq_length) && !$printmiss; $i++){
    push(@N, $n);
    
    
    
    # check for missing data
    if( ! keys (%OUTG)  || ($OUTG{$outseq}[$i] eq $missing) ||($OUTG{$outseq}[$i] eq $gap)  ){
	
	# $N[$i]--;
	$folded = 1;
    }
    # check for missing data
    for my $key (keys %SEQ){
	if(!$printmiss && ($SEQ{$key}[$i] eq $missing) && ($SEQ{$key}[$i] eq $gap)){
	    print STDERR $i, "\t", $SEQ{$key}[$i], "\n";
	    $N[$i]--;
	}
    }
    #####################################
    ########### FOLDED  ################
	
	
    if($folded){
	%V = ();
	@SORT_OCC = ();
	
	for my $key(keys %SEQ){
	    # if missing then skip
	    if( ($SEQ{$key}[$i] eq $missing) || ($SEQ{$key}[$i] eq $gap) ){ print STDERR "MISSING: $i\t$SEQ{$key}[$i]\n"; next; }
	    $V{ $SEQ{$key}[$i] } ++;
	}
	@SORT_OCC = sort{ $V{$b} <=> $V{$a} } keys %V;
	push( @ANCESTOR, $SORT_OCC[0] );
	
	$derived_size = 0;
	for(my $j = 1; $j < @SORT_OCC; $j++){
	    push @{ $DERIVED{$i} }, $SORT_OCC[$j];
	    $derived_size += $V{ $SORT_OCC[$j] };
	}
	if(! $derived_size ){
	    push @{ $DERIVED{$i} }, "-";
	}
	
	
	
	push @X, $derived_size;
	
    }
    
    ### UNFOLDED ###
    ###############
    elsif($folded == 0){
	@SORT_ANC = ();
	%V = ();

	# find the ancestor state
	# this is the state
	# in the outgroup sequences
	# with the maximum appearance

	for my $key (keys %OUTG){
	    $V{ $OUTG{$key}[$i] } ++;
	}
	# based on the values sort the hash
	# the key with the maximum value is the ancestor one

	@SORT_ANC = sort{ $V{$b} <=> $V{$a} } keys %V; # descending sort
	push @ANCESTOR, $SORT_ANC[0];

	# we compare each state in the sequences
	# with the ancestor one
	# if it is different we have a polymorphism

	$derived_size = 0;
	%D = (); # store the states of the sample
	
	for my $key (keys %SEQ){
	    if($key ne $missing){
		$D{ $SEQ{$key}[$i] } ++;
	    }
	}
	for my $s ( keys %D){
	    if( $ANCESTOR[$i] !~ /$s/){
		$derived_size += $D{$s};
		push @{ $DERIVED{$i} }, $s;
	    }
	}
	if(! $derived_size){ push @{ $DERIVED{$i}}, "-" ;}

	# if all the nucleotides in the %SEQ are different than ancestor and there are more than one states
	# then the ancestor nucleotide gives no information, so the position is actually folded
	if( ($derived_size == keys (%SEQ) ) && ( keys(%D) > 1) ){
	   
	    $folded = 1;
	    splice(@ANCESTOR, $i, 1);
	    @{ $DERIVED{$i} } = ();
	    $derived_size = 0;
	    %D = ();
	     $i --;
	    next;
	}

	push @X, $derived_size;

    } # END OF UNFOLDED CASE

    # binary output
    if($outformat =~ /^hap/i || $NOREC ){

	my $suminvar = 0;
	for my $key(keys %SEQ){
	    if($ANCESTOR[$i] ne $SEQ{$key}[$i]) { $suminvar++;}
	}
	
	# print $suminvar, " ***** \n";
	if( (! $invarout && $fixed ) # no mutation has happened or mutation has happened
	    || ($invarout && $suminvar && $fixed) # mutation has happened (polymorphic or fixed) 
	    || ( $suminvar && ( $suminvar != keys(%SEQ) ) && !$fixed && $invarout ) # polymorphic
	    || ( ! $invarout && ($suminvar != keys(%SEQ) ) && !$fixed )
	    ) { # polymorphic

	    for my $key(sort keys %SEQ){
		if($ANCESTOR[$i] ne $SEQ{$key}[$i]) { 
		    push( @{ $BinDERIVED{ $key } }, 1); 
		#    print STDERR "1\n";
		}

		
#		if($ANCESTOR[$i] ne $SEQ{$key}[$i]) { push( @{ $BinDERIVED{ $SEQ{$key} } }, 1); }
		#else { push( @{ $BinDERIVED{ $SEQ{$key} } }, 0); }
		else { 
		    push( @{ $BinDERIVED{ $key } }, 0);
		 #   print STDERR "0\n";
		}
	    }
	    if($postype == 1){
		$rel_position = $position[$i];
	    }
	    elsif($postype == 0){
		$rel_position = $position[$i]/$max_length;
		
	    }
	    push @REL_POSITION, $rel_position;
	    $segsites++;
	    #print STDERR $position[$i], "\t", $rel_position, "\n";
	}

	
    }

    elsif($outformat =~ /^sf$/i){
	print $position[$i], "\t";
	print $X[$i], "\t";
	print $N[$i], "\t";
	print $folded, "\t";
	print $ANCESTOR[$i], "\t";
	foreach my $d ( @{ $DERIVED{$i} } ){
	    print $d;
	}
	print "\n";
    }
    
    $folded = 0;

}


# print with missing information
for(my $i = 0; ($i < $seq_length) && $printmiss; $i++){
    push(@N, $n);
    
    
    
    # check for missing data
    if( ! keys (%OUTG)  || ($OUTG{$outseq}[$i] eq $missing)||($OUTG{$outseq}[$i] eq $gap) ){
	
	# $N[$i]--;
	$folded = 1;
    }
    # check for missing data
    $var = 0;
    for my $key (keys %SEQ){
	if($printmiss && ( ($SEQ{$key}[$i] eq $missing) || ($SEQ{$key}[$i] eq $gap)) ){
	    #print STDERR $i, "\t", $SEQ{$key}[$i], "\n";
	    $N[$i]--;
	}
    }
    #####################################
    ########### FOLDED  ################
	   
    if($folded){
	%V = ();
	@SORT_OCC = ();
	
	for my $key(keys %SEQ){
	    # if missing then skip
	    #if( $SEQ{$key}[$i] eq $missing){ print STDERR "MISSING: $i\t$SEQ{$key}[$i]\n"; next; }
	    if( ($SEQ{$key}[$i] ne $missing) && ($SEQ{$key}[$i] ne $gap) ){
		$V{ $SEQ{$key}[$i] } ++;
	    }
	}
	
	@SORT_OCC = sort{ $V{$b} <=> $V{$a} } keys %V;
	if($V{$SORT_OCC[0]} > 0){
	    push( @ANCESTOR, $SORT_OCC[0] );
	}
	else{
	    push( @ANCESTOR, "-");
	}
	
	$derived_size = 0;
	for(my $j = 1; $j < @SORT_OCC; $j++){
	    push @{ $DERIVED{$i} }, $SORT_OCC[$j];
	    $derived_size += $V{ $SORT_OCC[$j] };
	}
	if(! $derived_size ){
	    push @{ $DERIVED{$i} }, "-";
	}
	
	
	
	push @X, $derived_size;
	
    }
	
    
    ### UNFOLDED ###
    ###############
    elsif($folded == 0){
	@SORT_ANC = ();
	%V = ();
	
	# find the ancestor state
	# this is the state
	# in the outgroup sequences
	# with the maximum appearance

	for my $key (keys %OUTG){
	    $V{ $OUTG{$key}[$i] } ++;
	}
	# based on the values sort the hash
	# the key with the maximum value is the ancestor one

	@SORT_ANC = sort{ $V{$b} <=> $V{$a} } keys %V; # descending sort
	push @ANCESTOR, $SORT_ANC[0];

	# we compare each state in the sequences
	# with the ancestor one
	# if it is different we have a polymorphism

	$derived_size = 0;
	%D = (); # store the states of the sample
	
	for my $key (keys %SEQ){
	    if(($key ne $missing) && ($key ne $gap)){
		$D{ $SEQ{$key}[$i] } ++;
	    }
	}
	for my $s ( keys %D){
	    if( $ANCESTOR[$i] !~ /$s/){
		$derived_size += $D{$s};
		push @{ $DERIVED{$i} }, $s;
	    }
	}
	if(! $derived_size){ push @{ $DERIVED{$i}}, "-" ;}

	# if all the nucleotides in the %SEQ are different than ancestor and there are more than one states
	# then the ancestor nucleotide gives no information, so the position is actually folded
	if( ($derived_size == keys (%SEQ) ) && ( keys(%D) > 1) ){
	   
	    $folded = 1;
	    splice(@ANCESTOR, $i, 1);
	    @{ $DERIVED{$i} } = ();
	    $derived_size = 0;
	    %D = ();
	     $i --;
	    next;
	}

	push @X, $derived_size;

    } # END OF UNFOLDED CASE
    $missingseqs=0;
    # binary output
    if($outformat =~ /^hap/i){

	my $suminvar = 0;
	for my $key(keys %SEQ){
	    if( ($ANCESTOR[$i] ne $SEQ{$key}[$i]) && ($SEQ{$key}[$i] ne $gap) && ($SEQ{$key}[$i] ne $missing)) { 
		$suminvar++;
	    }
	    elsif(($SEQ{$key}[$i] eq $gap) || ($SEQ{$key}[$i] eq $missing)){
		$missingseqs++;
	    }
	}
	
	# print $suminvar, " ***** \n";
	if( (! $invarout && $fixed ) # no mutation has happened or mutation has happened
	    || ($invarout && $suminvar && $fixed) # mutation has happened (polymorphic or fixed) 
	    || ( $suminvar && ( $suminvar != keys(%SEQ)-$missingseqs ) && !$fixed && $invarout ) # polymorphic
	    || ( ! $invarout && ($suminvar != keys(%SEQ)-$missingseqs ) && !$fixed )
	    ) { # polymorphic

	    for my $key(sort keys %SEQ){
		if( ($ANCESTOR[$i] ne $SEQ{$key}[$i]) && ($SEQ{$key}[$i] ne $gap) &&($SEQ{$key}[$i] ne $missing)){ 
		    push( @{ $BinDERIVED{ $key } }, 1); 
		#    print STDERR "1\n";
		}

		
#		if($ANCESTOR[$i] ne $SEQ{$key}[$i]) { push( @{ $BinDERIVED{ $SEQ{$key} } }, 1); }
		#else { push( @{ $BinDERIVED{ $SEQ{$key} } }, 0); }
		elsif(($ANCESTOR[$i] eq $SEQ{$key}[$i] )&& ($SEQ{$key}[$i] ne $gap) &&($SEQ{$key}[$i] ne $missing)) { 
		    push( @{ $BinDERIVED{ $key } }, 0);
		 #   print STDERR "0\n";
		}
		elsif($SEQ{$key}[$i] == '-'){
		    push( @{ $BinDERIVED{ $key } }, 2);
		}
	    }
	    if($postype == 1){
		$rel_position = $position[$i];
	    }
	    elsif($postype == 0){
		$rel_position = $position[$i]/$max_length;
		
	    }
	    push @REL_POSITION, $rel_position;
	    $segsites++;
	    #print STDERR $position[$i], "\t", $rel_position, "\n";
	}

	
    }

    elsif($outformat =~ /^sf$/i){
	print $position[$i], "\t";
	print $X[$i], "\t";
	print $N[$i], "\t";
	print $folded, "\t";
	print $ANCESTOR[$i], "\t";
	foreach my $d ( @{ $DERIVED{$i} } ){
	    print $d;
	}
	print "\n";
    }
    
    $folded = 0;

}



if($NOREC){ # find the blocks that no recombination is needed

    # Thu Mar 20 ... to be continued
    
    open(RECOMB, ">$nexusfile.RECOMB") or die "Coudln't open RECOMB file\n";
    
    # $s1 = "00";
    # $s2 = "01";
    # $s3 = "10";
    # $s4 = "11";
    %STATES = ();
    $end = 0;
    
    foreach my $key( keys %BinDERIVED){
	print STDERR $key, "\n";
    }
    print STDERR "$akey\t",scalar @{ $BinDERIVED{$akey} },"\n";
    for( my $ii = 0; $ii < $#{ $BinDERIVED{$akey} }; $ii++){
	for( my $iii = $ii + 1 ; $iii <  @{ $BinDERIVED{$akey} }; $iii++){
	    %STATES=();
	    foreach my $jj (keys %SEQ){
		print STDERR "$ii\t$iii\t$jj\n";
		$s = $BinDERIVED{$jj}[$ii].$BinDERIVED{$jj}[$iii];
		$STATES{$s} = 1;
	    }
	    
	    if(keys(%STATES) == 4 || $iii == $#{ $BinDERIVED{$akey} } ) {
		
		if( $iii == $#{ $BinDERIVED{$akey} } ){
		    print RECOMB $ii, "\t", $iii, "\n";
		    $end = 1;
		    last;
		}
		else{
		    print RECOMB $ii, "\t", $iii - 1, "\n";
		    last;
		}
	    }
	    
	}
	if( $end){ last;}
    }
    
    
    #exit;   
}

print STDERR "Final Position is : ", $position[ $#position ], "\n";

@GROUPS=();
if($groups ne ""){
    @GROUPS=split(/,/, $groups);
}

if($outformat =~ /^hap/i){
    #print "ms like $nexusfile\t".keys (%SEQ)."\t1 -s $segsites\t\#invar:$invarout\n1\t2\t3\n\n";
    if($groups ne ""){
	foreach my $g (@GROUPS){
	    open(G, ">$nexusfile.$g.ssw") or die "Could't open outfile group $g\n";
	    
	    print G "\n\n//\n";
	    print G "segsites: $segsites\n";
	    print  G "positions:\t";
	    my $p=0;
	    for my $ss (@REL_POSITION){
		#print STDERR $position[$p],", "; 
		$p++;
		print G "$ss\t";
	    }
	    print G "\n";
	    print STDERR "\n";
	    $ss = 0;
	    for my $key (sort keys %BinDERIVED ){
		if($key !~ /$g/){next; }
		print STDERR "**** ", $key, "\t", $g, "\n";
		$ss++;
		print STDERR $ss, "\t", $key, "\n";	
		$size = @{ $BinDERIVED{$key} };
		for (my $j = 0; $j < $size; $j++ ){
		    print G "$BinDERIVED{$key}[$j]";
		}
		print G "\n";
	    }
	}
    }
    else{
	print  "\n\n//\n";
	    print  "segsites: $segsites\n";
	    print   "positions:\t";
	    my $p=0;
	    for my $ss (@REL_POSITION){
		#print STDERR $position[$p],", "; 
		$p++;
		print  "$ss\t";
	    }
	    print  "\n";
	    print STDERR "\n";
	    $ss = 0;
	    for my $key (sort keys %BinDERIVED ){
		
		$ss++;
		print STDERR $ss, "\t", $key, "\n";	
		$size = @{ $BinDERIVED{$key} };
		for (my $j = 0; $j < $size; $j++ ){
		    print  "$BinDERIVED{$key}[$j]";
		}
		print  "\n";
	    }
    }
}



$theta_denom = 0;
for(my $i = 1; $i < keys (%SEQ); $i++){
    $theta_denom += 1/$i;
}

$theta_watterson = $segsites / $theta_denom / $actlength;

open(FINFO, ">$nexusfile.info") or die "Could't open for output\n";

print STDERR "theta_watterson_per_nucleotide: $theta_watterson\n"; 
print FINFO "$actlength\n";

close LOG1;
