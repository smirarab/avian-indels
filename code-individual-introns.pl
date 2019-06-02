#!/usr/bin/env perl

use warnings;
use strict;
use Data::Dumper;

############################################################################
# Initialize variables
############################################################################
my($paupexec) = "paup4a147_osx"; # name of the paup executable

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
my($catfile) = "$outputdir" . "/" ."intronic_locus_" . "$locus" . "_indels";

############################################################################
# Read the input files
############################################################################
my($partf) = "$dir" . "/" . "$locus" . "/sate.removed.intron.noout-allgap.filtered.part";
my($dataf) = "$dir" . "/" . "$locus" . "/sate.removed.intron.noout.aligned-allgap.filtered";

system("sed \"s/-/ /g\" $partf > temp.parts.txt");

open (my $INF, "temp.parts.txt") or die "Could not open file temp.parts.txt for input.\n";
my @partinfo = <$INF>; # Read the taxon list file
close($INF) or die "Could not close file temp.parts.txt\n";
my($partnum) = $#partinfo + 1;

# Also read an authority file listing all taxa
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

############################################################################
# Calculate number of sites and the start and stop positions
############################################################################
(@temparray) = split(/\s+/, $partinfo[$partnum-1]);
my($length) = $temparray[4];

# variables for the indel coding
my @rawindeldata;
my($rawindellen);
my @indellist;
my($indellistlen);
my($ntax);
my($nchar);
my($nbp); # number of base pairs
my($nindel); # number of indel characters
my($taxname);
my($seq);
my($totalindelcharacters) = 0;

# open a log file for output
open (my $LOGF, ">$catfile.log.txt") or die "Could not open file $catfile.log.txt for output.\n";

print "locus $locus -- total length = $length\n\n";
print $LOGF "gap coding for locus $locus (total length = $length bp)\n";

# output each individual intron
for ($iter=0; $iter<$partnum; $iter++) {

  # prepare the fasta format alignments for indel coding --
	# extract each individual intron, output them as locus_X_intron_Y.aln.fasta
	(@temparray) = split(/\s+/, $partinfo[$iter]);
	print "intron $temparray[1] -- start = $temparray[3] -- end = $temparray[4]";
	$nbp = $temparray[4] - $temparray[3] + 1;
	print " (number of bp = $nbp";
	print ")\n";
	print $LOGF "-- locus $locus intron $temparray[1] -- length = $nbp bp ";
	$tempvar = "$outputdir" . "/" ."locus_" . "$locus" . "_intron_" . "$temparray[1]" . ".aln";
	system("./extract-intron $dataf $temparray[3] $temparray[4] $length > $tempvar.fasta");

	# call 2matrix (new version of 2xread) to conduct the indel coding
	system("./2matrix.pl -i $tempvar.fasta -n $tempvar.2xread -o x");
	system("./reformat-2xread-output.sh $tempvar.2xread.ss");
	
  # extract the indel data from the 2xread file --
	# first read the raw indel output from 2xread
	open (my $INDELF, "$tempvar.2xread.ss") or die "Could not open file $tempvar.2xread.ss for input.\n";
	@rawindeldata = <$INDELF>; # Read the taxon list file
	close($INDELF) or die "Could not close file $tempvar.2xread.ss\n";
	$rawindellen = $#rawindeldata + 1;
	($nchar,$ntax) = split(/\s+/, $rawindeldata[0]);
	print "   ntax = $ntax -- nchar = $nchar (includes nucleotides and indels)\n";
	$nindel = $nchar - $nbp;
	print "   number of indel characters = $nindel\n";
	print $LOGF " ngaps = $nindel\n";
	$totalindelcharacters = $totalindelcharacters + $nindel;
	print "   (total number of indel characters in concatenated file = $totalindelcharacters)\n\n";
	
	# then read the list of indels parsed from the 2xread output by reformat-2xread-output.sh
	open (my $LISTF, "$tempvar.2xread.ss.indel.list.txt") or die "Could not open file $tempvar.2xread.ss.indel.list.txt for input.\n";
	@indellist = <$LISTF>; # Read the indel list file
	close($LISTF) or die "Could not close file $tempvar.2xread.ss.indel.list.txt\n";
	$indellistlen = $#indellist + 1;
	
	# iterate through the list of indels, write file of indel characters
	open (my $GAPF, ">$tempvar.coded.gaps.phy") or die "Could not open file $tempvar.coded.gaps.phy for output.\n";
	print $GAPF "$ntax $nindel\n";
	for ($jter=0; $jter<$ntax; $jter++) {
		($taxname,$seq) = split(/\s+/, $rawindeldata[$jter+1]);
	# Remove comments to echo more rubbish (print that taxon names are found in the authority file) to the screen
	#	print " -- > $taxname";
		for ($lter=0; $lter<$totaltaxa; $lter++) {
			if ( $taxname eq $shortnames[$lter] ) { 
				$taxonpresent[$lter] = 1;
	#			print " -- found";
			}
		}
	#	print "\n";
		print $GAPF "$taxname     ";
		for ($kter=0; $kter<$indellistlen; $kter++) {
			# use data in indellist to extract the indels
			(@temparray) = split(/\s+/, $indellist[$kter]);
			$tempch = substr($seq,$temparray[0],1);
			print $GAPF "$tempch";
		} # end for ($kter=0; $kter<$indellistlen; ...
		print $GAPF "\n";
	} # end for ($jter=0; $jter<$ntax; ...
	close($GAPF) or die "Could not close file $tempvar.coded.gaps.phy\n";

	# iterate through the list of indels, write metadata for each indel
	# data in indellist --> position_in_seq "sequence_indel" position1 position2
	open (my $GAPINFOF, ">$tempvar.coded.gaps.info.txt") or die "Could not open file $tempvar.coded.gaps.info.txt for output.\n";
	print $GAPINFOF "indel_number\tstart\tend\tlength\n";
	for ($kter=0; $kter<$indellistlen; $kter++) {
		(@temparray) = split(/\s+/, $indellist[$kter]);
		$jter=$kter+1;
		print $GAPINFOF "$jter\t$temparray[2]\t$temparray[3]\t";
		$jter = $temparray[3] - $temparray[2] + 1;
		print $GAPINFOF "$jter\n";
	}
	close($GAPINFOF) or die "Could not close file $tempvar.coded.gaps.info.txt\n";
	
} # end for ($iter=0; $iter<$partnum; ...

print $LOGF "there are $totalindelcharacters indel characters coded for locus $locus\n\n";
close($LOGF) or die "Could not close file $catfile.log.txt\n";

############################################################################
# Now concatenate the introns from the same locus
############################################################################

my($totaltaxapresent) = 0;

# First echo the set of taxa that are present
print "\nTaxon presence for this locus is:\n";
for ($iter=0; $iter<$totaltaxa; $iter++) {
	if ( $taxonpresent[$iter] == 1) {
		print "  $shortnames[$iter] -- $longnames[$iter]\n";
		$totaltaxapresent++;
	}
	else {
		print "  $shortnames[$iter] -- $longnames[$iter] ***ABSENT***\n";
	}
}
print "\n";

print "There are a total of $totaltaxapresent taxa present ";
print "and a total of $totalindelcharacters indel characters\n\n";
print "Writing output files...";

# Name of concatenated nexus file
#my($catfile) = "$outputdir" . "/" ."intronic_locus_" . "$locus" . "_indels";
my @introngapinfo;
my($introngaplen);
my($introntaxa);
my($intronchars);
my($found);
my @intronnumbers;
my @intronstart;
my @intronend;

open (my $NEXF, ">$catfile.nex") or die "Could not open file $catfile.nex for output.\n";
open (my $NEXLONGF, ">$catfile.longname.nex") or die "Could not open file $catfile.longname.nex for output.\n";

print $NEXF "#NEXUS\n\n";
print $NEXF "[ Indel data (Simmons-Ochoterena coding) for Jarvis et al. intronic locus $locus ]\n";
print $NEXF "Begin data;\n";
print $NEXF "\tdimensions ntax=$totaltaxapresent nchar=$totalindelcharacters;\n";
print $NEXF "\tformat datatype=standard missing=? gap=- interleave;\n";
print $NEXF "Matrix\n\n";

print $NEXLONGF "#NEXUS\n\n";
print $NEXLONGF "[ Indel data (Simmons-Ochoterena coding) for Jarvis et al. intronic locus $locus ]\n";
print $NEXLONGF "Begin data;\n";
print $NEXLONGF "\tdimensions ntax=$totaltaxapresent nchar=$totalindelcharacters;\n";
print $NEXLONGF "\tformat datatype=standard missing=? gap=- interleave;\n";
print $NEXLONGF "Matrix\n\n";

# iterate through and output each individual intron
for ($iter=0; $iter<$partnum; $iter++) {

	# extract each individual intron
	(@temparray) = split(/\s+/, $partinfo[$iter]);
	$tempvar = "$outputdir" . "/" ."locus_" . "$locus" . "_intron_" . "$temparray[1]" . ".aln";
	
	# read the intron file
	open (my $INTRONGAPF, "$tempvar.coded.gaps.phy") or die "Could not open file $tempvar.coded.gaps.phy for input.\n";
	@introngapinfo = <$INTRONGAPF>; # Read the phylip format gap data file
	close($INTRONGAPF) or die "Could not close file $tempvar.coded.gaps.phy\n";
	($introntaxa,$intronchars) = split(/\s+/, $introngapinfo[0]);
	chomp($intronchars);
	# set the intron start and end
	$intronnumbers[$iter] = $temparray[1];
	if ( $iter == 0 ) {	$intronstart[$iter] = 1; }
	else { $intronstart[$iter] = $intronend[$iter-1] + 1; }
	if ( $iter == 0 ) {	$intronend[$iter] = $intronchars; }
	else { $intronend[$iter] = $intronend[$iter-1] + $intronchars; }
	
	if ( $intronchars > 0 ) {

		print $NEXF "[ Indel data for locus $locus intron $temparray[1] ]\n";
		print $NEXF "[ characters $intronstart[$iter] to $intronend[$iter] ";
		print $NEXF "($intronchars characters scored for $introntaxa taxa) ]\n";
	
		print $NEXLONGF "[ Indel data for locus $locus intron $temparray[1] (long taxon names) ]\n";
		print $NEXLONGF "[ characters $intronstart[$iter] to $intronend[$iter] ";
		print $NEXLONGF "($intronchars characters scored for $introntaxa taxa) ]\n";
	
		# iterate through the taxa and output in the order of the authority file
		for ($jter=0; $jter<$totaltaxa; $jter++) {
			if ( $taxonpresent[$jter] == 1) {
				$found=0;
				for ($kter=0; $kter<$introntaxa; $kter++) {
					($taxname,$seq) = split(/\s+/, $introngapinfo[$kter+1]);
					chomp($seq);
					if ( $taxname eq $shortnames[$jter] ) {
						print $NEXF "$taxname      $seq\n";
						$found=1;
						print $NEXLONGF "$longnames[$jter]";
						$nter = $longnamelen + 4 - length($longnames[$jter]) ;
						for ($mter=0; $mter<$nter; $mter++) { print $NEXLONGF " "; }
						print $NEXLONGF "$seq\n";
					}
				} # end for ($kter=0; $kter<$introntaxa; ...
				if ( $found == 0 ) {
					print $NEXF "$shortnames[$jter]      ";
					print $NEXLONGF "$longnames[$jter]";
					$nter = $longnamelen + 4 - length($longnames[$jter]) ;
					for ($mter=0; $mter<$nter; $mter++) { print $NEXLONGF " "; }
					for ($kter=0; $kter<$intronchars; $kter++) {
						print $NEXF "?";
						print $NEXLONGF "?";
					} #  end for ($kter=0; $kter<$intronchars; ...
					print $NEXF "\n";
					print $NEXLONGF "\n";
				}
			}
		}
		print $NEXF "\n";
		print $NEXLONGF "\n";
	
	} # if ( $intronchars > 0 ... (only print data if the partition has some indels
	else {
		print $NEXF "[ There are no indels in intron $temparray[1] of locus $locus ]\n\n";
		print $NEXLONGF "[ There are no indels in intron $temparray[1] of locus $locus ]\n\n";
	}

}

print $NEXF "\t;\nEnd;\n\n";
print $NEXLONGF "\t;\nEnd;\n\n";

# output a sets block
print $NEXF "Begin sets;\n";
print $NEXLONGF "Begin sets;\n";

# write raxml and iqtree partition files
open (my $RAXF, ">$catfile.raxml.parts") or die "Could not open file $catfile.raxml.parts for output.\n";
open (my $IQF, ">$catfile.iqtree.parts.nex") or die "Could not open file $catfile.iqtree.parts.nex for output.\n";

print $IQF "#NEXUS\n";
print $IQF "Begin sets;\n";

# write charsets for each intron
for ($iter=0; $iter<$partnum; $iter++) {
	if ( $intronend[$iter] > $intronstart[$iter] ) {
		print $NEXF "\tcharset locus_$locus";
		print $NEXF "_intron";
		print $NEXF "$intronnumbers[$iter] = $intronstart[$iter]-$intronend[$iter] ;\n";
		print $NEXLONGF "\tcharset locus_$locus";
		print $NEXLONGF "_intron";
		print $NEXLONGF "$intronnumbers[$iter] = $intronstart[$iter]-$intronend[$iter] ;\n";
		# raxml partition file
		print $RAXF "BIN, intron_$intronnumbers[$iter] = $intronstart[$iter]-$intronend[$iter]\n";
		# iq-tree partition file
		print $IQF "\tcharset locus_$locus";
		print $IQF "_intron";
		print $IQF "$intronnumbers[$iter] = $intronstart[$iter]-$intronend[$iter] ;\n";
	} # if ( $intronend[$iter] < $intronstart[$iter] ... 
	  # (only print the charsets for each intron if the intron has >0 indels)
}

print $IQF "End;\n";

close($RAXF) or die "Could not close file $catfile.raxml.parts\n";
close($IQF) or die "Could not close file $catfile.iqtree.parts.nex\n";

# array to hold gap length information
my @gaplen;
my @gaplendata;
my @gaplengthinfo;
my($gapinfofile);
my($gaplengthinfolen);

# iterate through the gap information for each intron
# fill the array @gaplen with the length of each gap
# also write the gap length information to a data file

# $GAPCOMBF is the list of gap info for all characters
# in the file, it has the following data:
#    1. character in matrix
#    2. intron number
#    3. indel number within intron
#    4. starting nucleotide within intron
#    5. ending nucleotide within intron
#    6. indel length
open (my $GAPCOMBF, ">$catfile.gaplengths.txt") or die "Could not open file $catfile.gaplengths.txt for output.\n";
print $GAPCOMBF "character\tintron\tindel_number\tstart\tend\tlength\n";

for ($iter=0; $iter<$partnum; $iter++) {
	
	# open the gap info file
	$gapinfofile = "$outputdir" . "/" . "locus_" . "$locus" . "_intron_" . "$intronnumbers[$iter]" . ".aln.coded.gaps.info.txt";	
	open (my $INFOGAPF, "$gapinfofile") or die "Could not open file $gapinfofile for input.\n";
	@gaplengthinfo = <$INFOGAPF>; # Read the gap length information file
	close($INFOGAPF) or die "Could not close file $gapinfofile\n";
	$gaplengthinfolen = $#gaplengthinfo + 1;
	
#	remove commented "print" statements to print more rubbish (information on gap lengths) to the screen
#	print "for intron $intronnumbers[$iter]\n";
	for ($jter=1; $jter<$gaplengthinfolen; $jter++) {
		$mter = $jter + $intronstart[$iter] - 2 ;
		(@gaplendata) = split(/\s+/, $gaplengthinfo[$jter]);
		$gaplen[$mter] = $gaplendata[3];
		chomp($gaplen[$mter]);
		print $GAPCOMBF "$mter\t$intronnumbers[$iter]\t$gaplendata[0]\t";
		print $GAPCOMBF "$gaplendata[1]\t$gaplendata[2]\t$gaplen[$mter]\n"
	}
#	print "\n";

}

close($GAPCOMBF) or die "Could not close file $catfile.gaplengths.txt\n";

print $NEXF "\t[ four charsets based on indel length: ]\n";
print $NEXF "\t[   charset len2plus   = >1 bp indels  ]\n";
print $NEXF "\t[   charset len10plus  = >9 bp indels  ]\n";
print $NEXF "\t[   charset len50plus  = >49 bp indels ]\n";
print $NEXF "\t[   charset len100plus = >99 bp indels ]\n";

print $NEXLONGF "\t[ four charsets based on indel length: ]\n";
print $NEXLONGF "\t[   charset len2plus   = >1 bp indels  ]\n";
print $NEXLONGF "\t[   charset len10plus  = >9 bp indels  ]\n";
print $NEXLONGF "\t[   charset len50plus  = >49 bp indels ]\n";
print $NEXLONGF "\t[   charset len100plus = >99 bp indels ]\n";

print $NEXF "\tcharset locus_$locus";
print $NEXF "_len2plus = ";
print $NEXLONGF "\tcharset locus_$locus";
print $NEXLONGF "_len2plus = ";
for ($iter=0; $iter<$totalindelcharacters; $iter++) {
	if ( $gaplen[$iter] > 1 ) {
		$mter = $iter + 1;
		print $NEXF "$mter ";
		print $NEXLONGF "$mter ";
	}
}
print $NEXF ";\n";
print $NEXF "\tcharset locus_$locus";
print $NEXF "_len10plus = ";
print $NEXLONGF ";\n";
print $NEXLONGF "\tcharset locus_$locus";
print $NEXLONGF "_len10plus = ";
for ($iter=0; $iter<$totalindelcharacters; $iter++) {
	if ( $gaplen[$iter] > 9 ) {
		$mter = $iter + 1;
		print $NEXF "$mter ";
		print $NEXLONGF "$mter ";
	}
}
print $NEXF ";\n";
print $NEXF "\tcharset locus_$locus";
print $NEXF "_len50plus = ";
print $NEXLONGF ";\n";
print $NEXLONGF "\tcharset locus_$locus";
print $NEXLONGF "_len50plus = ";
for ($iter=0; $iter<$totalindelcharacters; $iter++) {
	if ( $gaplen[$iter] > 49 ) {
		$mter = $iter + 1;
		print $NEXF "$mter ";
		print $NEXLONGF "$mter ";
	}
}
print $NEXF ";\n";
print $NEXF "\tcharset locus_$locus";
print $NEXF "_len100plus = ";
print $NEXLONGF ";\n";
print $NEXLONGF "\tcharset locus_$locus";
print $NEXLONGF "_len100plus = ";
for ($iter=0; $iter<$totalindelcharacters; $iter++) {
	if ( $gaplen[$iter] > 99 ) {
		$mter = $iter + 1;
		print $NEXF "$mter ";
		print $NEXLONGF "$mter ";
	}
}
print $NEXF ";\n";
print $NEXLONGF ";\n";

print $NEXF "End;\n\n";
print $NEXLONGF "End;\n\n";

close($NEXF) or die "Could not close file $catfile.nex\n";
close($NEXLONGF) or die "Could not close file $catfile.longname.nex\n";

system("cp $catfile.nex temp.nex");
open (my $APPF, ">>temp.nex") or die "Could not open file temp.nex for append.\n";
print $APPF "Begin paup;\n";
print $APPF "\texport file=$catfile.phy charsperline=all replace;\n";
print $APPF "\tquit;\n";
print $APPF "End;\n";
close($APPF) or die "Could not close file temp.nex\n";
system("$paupexec temp.nex > temp.paup.screen");

system("cp $catfile.longname.nex temp.longname.nex");
open (my $APPLONGF, ">>temp.longname.nex") or die "Could not open file temp.longname.nex for append.\n";
print $APPLONGF "Begin paup;\n";
print $APPLONGF "\texport file=$catfile.longname.phy charsperline=all replace;\n";
print $APPLONGF "\tquit;\n";
print $APPLONGF "End;\n";
close($APPLONGF) or die "Could not close file temp.longname.nex\n";
system("$paupexec temp.longname.nex > temp.paup.screen");

# write the list of output files
print "done\n\n";
print "Nexus file              -- intronic_locus_$locus" . "_indels.nex\n";
print "Phylip file             -- intronic_locus_$locus" . "_indels.phy\n";
print "Nexus file (longnames)  -- intronic_locus_$locus" . "_indels.longname.nex\n";
print "Phylip file (longnames) -- intronic_locus_$locus" . "_indels.longname.phy\n";
print "IQ-TREE partition file  -- intronic_locus_$locus" . "_indels.iqtree.parts.nex\n";
print "RAxML partition file    -- intronic_locus_$locus" . "_indels.raxml.parts\n";
print "indel data file         -- intronic_locus_$locus" . "_indels.gaplengths.txt\n\n";
print "Log file                -- intronic_locus_$locus" . "_indels.log.txt\n\n";

# clean up the temporary files
system("rm temp.parts.txt temp.nex temp.longname.nex temp.paup.screen");

