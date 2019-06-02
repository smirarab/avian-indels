#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Initialize variables
############################################################################
# my($paupexec) = "paup4a147_osx"; # name of the paup executable

my($progname) = $0;
my($iter);
my($jter);
my($kter);
my($lter);
my($mter);
my($nter);

my($tempch);
my($tempvar);
my @temparray;

if ( @ARGV != 4 ) {
	print "Usage:\n  \$ $progname <dir> <locus #> <auth> <outdir>\n";
	print "  dir     = directory with alignments\n";
	print "  locus # = number of the locus\n";
	print "  auth    = authority file (list of taxa)\n";
	print "  outdir  = directory for output\n";
	print "exiting...\n";
	exit;
}

my($dir)=$ARGV[0];
my($locus)=$ARGV[1];
my($authfile)=$ARGV[2];
my($outputdir)=$ARGV[3];

if ( -d "$outputdir" ) {
	print "Data will be saved in directory $outputdir\n\n";
}
else {
	system("mkdir $outputdir");
	print "Directory $outputdir created.\n";
	print "Data will be saved in that directory\n\n";
}
# Root name of concatenated nexus file
my($catfile) = "$outputdir" . "/" ."uce_locus_" . "$locus" . "_indels";

############################################################################
# Read the input files
############################################################################
my($dataf) = "$dir" . "/" . "$locus" . "_s.fasta";

# Read the authority file listing all taxa
open (my $AUTHF, $authfile) or die "Could not open file $authfile for input.\n";
my @authinfo = <$AUTHF>; # Read the taxon list file
close($AUTHF) or die "Could not close file $authfile\n";
my($authlen) = $#authinfo + 1;

my($totaltaxa);
my($shortnamelen);
my($longnamelen);
my @shortnames;
my @longnames;
my @taxonpresent;

($totaltaxa,$shortnamelen,$longnamelen) = split(/\s+/, $authinfo[0]);

print "The complete file will have up to $totaltaxa taxa:\n";
for ($iter=0; $iter<$totaltaxa; $iter++) {
	($shortnames[$iter],$longnames[$iter]) = split(/\s+/, $authinfo[$iter+1]);
	$taxonpresent[$iter] = 0;
	print "  $shortnames[$iter] -- $longnames[$iter]\n";
}
print "(some taxa may not be included for this locus)\n\n";

# recode the amino acid data to simplify output from 2matrix
system("./aa2dummynt.pl $dataf temp.recoded.aa.txt");

# Read the recoded fasta file (temp.recoded.aa.txt)
open (my $FAF, "temp.recoded.aa.txt") or die "Could not open file temp.recoded.aa.txt for input.\n";
my @fastainfo = <$FAF>; # Read the input fasta file
close($FAF) or die "Could not close file temp.recoded.aa.txt\n";
my($fastalen) = $#fastainfo + 1;
$tempvar = $fastainfo[1];
my($nbp) = length($tempvar);
print "Alignment has $nbp sites\n\n";

############################################################################
# Perform the indel coding and create output files
############################################################################

# variables for the indel coding
my($ntax);
my($indelcharstart);
my($indelcharend);
my($nindel); # number of indel characters
my($taxname);
my($seq);
my($indelstring);
my($taxspaces);

my @indeldescription;	my($indellength);
# The indel description holds info about each indel. The order is:
#   [0] = character in 2x read file
#   [1] = dummy, simply says "indel"
#   [2] = first nucleotide position of indel
#   [3] = last nucleotide position of indel
#   [4]... additional rubbish

# open a log file for output
open (my $LOGF, ">$catfile.log.txt") or die "Could not open file $catfile.log.txt for output.\n";

print "locus $locus -- total length = $nbp\n";
print $LOGF "gap coding for locus $locus (total length = $nbp base pairs)\n";

# call 2matrix (new version of 2xread) to conduct the indel coding
$tempvar = "$outputdir" . "/" ."locus_" . "$locus" . "_uce";
system("./2matrix.pl -i temp.recoded.aa.txt -n $tempvar.2xread -o x");

# extract the indel data from the 2xread file --
#-----------------------------------------------------------------------------------------
# first, create a file listing the indels
system("grep \"_indel_\" $tempvar.2xread.ss > $tempvar.2xread.indel.list.txt");
system("sed \"s/{//g\" $tempvar.2xread.indel.list.txt > t.0");
system("mv t.0 $tempvar.2xread.indel.list.txt");
system("sed \"s/_indel_/ indel /g\" $tempvar.2xread.indel.list.txt > t.0");
system("mv t.0 $tempvar.2xread.indel.list.txt");
system("sed \"s/-/ /g\" $tempvar.2xread.indel.list.txt > t.0");
system("mv t.0 $tempvar.2xread.indel.list.txt");

#-----------------------------------------------------------------------------------------
# second, read the indel list
open (my $INDF, "$tempvar.2xread.indel.list.txt") or die "Could not open file $tempvar.2xread.indel.list.txt for input.\n";
my @indelinfo = <$INDF>; # Read the indel list file
close($INDF) or die "Could not close file $tempvar.2xread.indel.list.txt\n";
$nindel = $#indelinfo + 1;
print "locus $locus -- number of gap characters = $nindel\n\n";
print $LOGF "gap coding for locus $locus (total number of indels = $nindel)\n";
print $LOGF "data for locus\t$locus\t$nbp\t$nindel\t";
@temparray = split(/\s+/, $indelinfo[0]);
$indelcharstart = $temparray[0];	#$indelcharstart--;

#-----------------------------------------------------------------------------------------
# third, read the 2xread output file
open (my $SSF, "$tempvar.2xread.ss") or die "Could not open file $tempvar.2xread.ss for input.\n";
my @indeldata = <$SSF>;
close($SSF) or die "Could not close file $tempvar.2xread.ss\n";
($nter,$ntax) = split(/\s+/, $indeldata[3]);

#-----------------------------------------------------------------------------------------
# fouth, write the phylip and nexus format files
print "Writing data for $nindel indels for $ntax taxa\n";

if ( $nindel == 0 ) { 
	print $LOGF "0\t0\t0\t0\n";
	exit;
}

open (my $PHYF, ">$tempvar.indels.phy") or die "Could not open file $tempvar.indels.phy for output.\n";
open (my $LPHYF, ">$tempvar.indels.longnames.phy") or die "Could not open file $tempvar.indels.longnames.phy for output.\n";

print $PHYF "$ntax  $nindel\n";
print $LPHYF "$ntax  $nindel\n";

open (my $NEXF, ">$tempvar.indels.nex") or die "Could not open file $tempvar.indels.nex for output.\n";
open (my $LNEXF, ">$tempvar.indels.longnames.nex") or die "Could not open file $tempvar.indels.longnames.nex for output.\n";

# Output the nexus headers
print $NEXF "#NEXUS\n\n";
print $NEXF "[ Indel data (Simmons-Ochoterena coding) for uce of Jarvis et al. locus $locus ]\n";
print $NEXF "Begin data;\n";
print $NEXF "\tdimensions ntax=$ntax nchar=$nindel;\n";
print $NEXF "\tformat datatype=standard missing=? gap=- interleave;\n";
print $NEXF "Matrix\n\n";

print $LNEXF "#NEXUS\n\n";
print $LNEXF "[ Indel data (Simmons-Ochoterena coding) for uce of Jarvis et al. locus $locus ]\n";
print $LNEXF "Begin data;\n";
print $LNEXF "\tdimensions ntax=$ntax nchar=$nindel;\n";
print $LNEXF "\tformat datatype=standard missing=? gap=- interleave;\n";
print $LNEXF "Matrix\n\n";

# iterate through the taxa in the authority file
for ($iter=0; $iter<$totaltaxa; $iter++) {

	# iterate through the taxa in the 2xread file
	for ($jter=0; $jter<$ntax; $jter++) {
		$kter = $jter + 4;
		($taxname,$seq) = split(/\s+/, $indeldata[$kter]);
		chomp($seq);
		
		# output the sequence data
		if ( $taxname eq $shortnames[$iter] ) {
			# output the taxon names (pad longnames with spaces)
			print $PHYF "$taxname     ";	print $LPHYF "$longnames[$iter] ";
			print $NEXF "$taxname     ";	print $LNEXF "$longnames[$iter] ";
			$taxspaces = $longnamelen - length($longnames[$iter]) + 1;
			for ($lter=0; $lter<$taxspaces; $lter++) { print $LPHYF " "; print $LNEXF " "; }
			# now output the indel characters
			$indelstring = substr($seq,$indelcharstart,$nindel);
			print $PHYF "$indelstring\n";	print $LPHYF "$indelstring\n";
			print $NEXF "$indelstring\n";	print $LNEXF "$indelstring\n";
			# terminate the $jter loop
			$jter = $ntax;
		} # end if ( $taxname eq $shortnames[$iter] ...
		
	} # end for ($jter=0; $jter<$ntax; ...

} # end for ($iter=0; $iter<$totaltaxa; ...

close($PHYF) or die "Could not close file $tempvar.indels.phy\n";
close($LPHYF) or die "Could not close file $tempvar.indels.longnames.phy\n";

print $NEXF "\t;\nEnd;\n\n";	print $LNEXF "\t;\nEnd;\n\n";

# Now output a sets block, also count # of indels in each size class
#   note that I have kept the same variables as the cds coder (which was
#   based on codons, not nucleotides). The variables should be iterpreted
#   based on the charset info below.
my($gte2) = 0;	my($gte4) = 0;	my($gte10) = 0;	my($gte30) = 0;
print $NEXF "Begin sets;\n";
print $NEXF "\n\t[   charset len2plus   = >1 bp indels  ]\n";
print $NEXF "\t[   charset len10plus  = >9 bp indels  ]\n";
print $NEXF "\t[   charset len50plus  = >49 bp indels ]\n";
print $NEXF "\t[   charset len100plus = >99 bp indels ]\n";
print $NEXF "\tcharset locus_" . "$locus" . "_len2plus = ";
print $LNEXF "Begin sets;\n";
print $LNEXF "\n\t[   charset len2plus   = >1 bp indels  ]\n";
print $LNEXF "\t[   charset len10plus  = >9 bp indels  ]\n";
print $LNEXF "\t[   charset len50plus  = >49 bp indels ]\n";
print $LNEXF "\t[   charset len100plus = >99 bp indels ]\n";
print $LNEXF "\tcharset locus_" . "$locus" . "_len2plus = ";
for ($iter=0; $iter<$nindel; $iter++) {
	(@indeldescription) = split(/\s+/, $indelinfo[$iter]);
	$indellength = $indeldescription[3] - $indeldescription[2] + 1;
	if ( $indellength > 1 ) {
		$jter = $iter + 1;
		print $NEXF "$jter ";
		print $LNEXF "$jter ";
		$gte2++;
	}
}
print $NEXF ";\n";
print $LNEXF ";\n";
# indels of 10bp or longer
print $NEXF "\tcharset locus_" . "$locus" . "_len10plus = ";
print $LNEXF "\tcharset locus_" . "$locus" . "_len10plus = ";
for ($iter=0; $iter<$nindel; $iter++) {
	(@indeldescription) = split(/\s+/, $indelinfo[$iter]);
	$indellength = $indeldescription[3] - $indeldescription[2] + 1;
	if ( $indellength > 9 ) {
		$jter = $iter + 1;
		print $NEXF "$jter ";
		print $LNEXF "$jter ";
		$gte4++;
	}
}
print $NEXF ";\n";
print $LNEXF ";\n";
# indels of 50bp or longer
print $NEXF "\tcharset locus_" . "$locus" . "_len50plus = ";
print $LNEXF "\tcharset locus_" . "$locus" . "_len50plus = ";
for ($iter=0; $iter<$nindel; $iter++) {
	(@indeldescription) = split(/\s+/, $indelinfo[$iter]);
	$indellength = $indeldescription[3] - $indeldescription[2] + 1;
	if ( $indellength > 49 ) {
		$jter = $iter + 1;
		print $NEXF "$jter ";
		print $LNEXF "$jter ";
		$gte10++;
	}
}
print $NEXF ";\n";
print $LNEXF ";\n";
# indels of 100bp or longer
print $NEXF "\tcharset locus_" . "$locus" . "_len100plus = ";
print $LNEXF "\tcharset locus_" . "$locus" . "_len100plus = ";
for ($iter=0; $iter<$nindel; $iter++) {
	(@indeldescription) = split(/\s+/, $indelinfo[$iter]);
	$indellength = $indeldescription[3] - $indeldescription[2] + 1;
	if ( $indellength > 99 ) {
		$jter = $iter + 1;
		print $NEXF "$jter ";
		print $LNEXF "$jter ";
		$gte30++;
	}
}
print $NEXF ";\n";
print $LNEXF ";\n";
print $NEXF "\nEnd;\n\n";
print $LNEXF "\nEnd;\n\n";

close($NEXF) or die "Could not close file $tempvar.indels.nex\n";
close($LNEXF) or die "Could not close file $tempvar.indels.longnames.nex\n";

# see note above about the $gte2, $gte4, etc values
print $LOGF "$gte2\t$gte4\t$gte10\t$gte30\n";

# close the logfile and exit
close($LOGF) or die "Could not close file $catfile.log.txt\n";

exit;
