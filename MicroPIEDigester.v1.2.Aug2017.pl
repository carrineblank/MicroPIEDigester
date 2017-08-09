#!/usr/bin/perl -w

####################################################################
# This script creates a Mesquite readable nexus file from MicroPIE output
# file in comma-delimited text format .csv
# MicroPIEDigester version 1.2
# created by Carrine Blank on Sept-Oct 2016, Updated in Feb 2017, Jun 2017, and Aug 2017
#
# invoke by typing 'perl MicroPIEDigestor.pl infile'
#
# Assumes that the first column of the input filecontains the taxonomic names
# and that the names are identical to those in the destination Mesquite file.
# Ontologizes discrete text characters (like 'rod', 'rods', 'elongated cocci', 'bacilli', 'ovals').
####################################################################
use strict;
use warnings;
use utf8;

my $usage = "usage: perl MicroPIEDigester.pl infile\n";
my $infile = shift or die $usage;	# input Nexus MicroPIE file
my $unixfile = $infile."_unix";
my $unixfile2 = $infile."_tsv_unix";
my $outmessage = "The characters have been coded and saved to a Mesquite nexus file.\n";

#prepares the infile, converting the csv format to tab-delimited among other minor fixes
system ("LC_ALL=C tr '\r' '\n' <$infile >$unixfile");
system ("LC_ALL=C tr '\r' '\n' <$infile >$unixfile2");
open(IN, $unixfile) || die "Can't open $unixfile: $!\n";
open(OUT, '>', $unixfile2) || die "Can't open $unixfile2: $!\n";
while (my $row = <IN> ) {
	chomp $row;
	$row =~ s /,/\t/g; #converts all the commas to tabs
	$row =~ s/#/,/g; #converts all the # to commas
	$row =~ s/\|/,/g; #converts all the # to commas
	if ($row =~ /"Taxon"/) {
		$row =~ s/"//g;
	}
	else {
		$row =~ s /"/'/; #turns the first double to quote to a single quote around the taxa name
		$row =~ s /"/'/;#turns the second double to quote to a single quote around the taxa name
	}
	$row =~ s/\t"/	",/g; #puts a comma at the beginning of each cell
	$row =~ s/"\t/,"	/g; #puts a comma at the end of each cell
	print OUT "$row\n";
}
close IN;
close OUT;

#deposits character names into the array @charlist and asks user to select their desired character to process/code
open(IN, $unixfile2) || die "Can't open $unixfile2: $!\n";
my $i=1;
my $line = <IN>;
chomp $line;
$line =~ s/\t\t/\t/g;
my @charlist = split (/\t/,$line);
print "\n\nBelow is a list of character numbers and character types:\n";
#print OUT @charlist;
foreach (@charlist) {
	print $i++;
	print " $_\n";
	}
#asks user to choose a character to process and dumps that into $charnum
#cleans up the MicroPIE output a bit
print "\n\nEnter the CHARACTER NUMBER you wish to process:";
my $charnum = <STDIN>;
chomp $charnum;
my $desiredkey = $charnum - 1;
print "To confirm, your chosen character type is: $charlist[$desiredkey]\n";
print "Here is a list of your taxa and their raw character states:\n";

#saves the column with the header desiredkey into the outfile, which has the list of elements in the column
my $character = $charlist[$desiredkey];
my $rawmatrix = "raw.$character.txt";
local $, = "\t";
open (IN, '<', $unixfile2) or die $!;
open (OUT, '>', $rawmatrix) or die $!;
while (<IN>) {
	chomp;
	my $rawcharstates;
	my @columns = split /\t/;
	push (my @taxnames, $columns[0]);
	map {s/$/'/g; } @taxnames;  #add a single quote to the taxname
	map {s/^/'/g; } @taxnames;  #add a single quote to the taxname
	map {s/\s\'/'/g; } @taxnames;  #remove any spaces after the taxname
	push (my @rawcharstates, $columns[$desiredkey]);
	if (defined $columns[$desiredkey]) {
		map {s/, /,/g; } @rawcharstates;  #change # to commas
		map {s/#/,/g; } @rawcharstates;  #change # to commas
		map {s/\"//g; } @rawcharstates; #remove double quotes
		map {s/^\s//g; } @rawcharstates; # remove space at beginning of the string
		map {s/\s\s/ /g; } @rawcharstates; # remove double spaces
		map {s/^/,/g; } @rawcharstates; # puts a comma at the front of the string
		map {s/$/,/g; } @rawcharstates; # puts a comma at the end of the string
		map {s/,,/,/g; } @rawcharstates; # removes double commas
		map {s/no\s/not /gi; } @rawcharstates; # transforms no to not
		}
	else {
		@rawcharstates= ",";
		}
	print $columns[0], "@rawcharstates", "\n";
	print OUT "@taxnames", "@rawcharstates", "\n"; # outputs to $rawmatrix raw.charactername.txt
	}
sub median {
	my @vals = sort {$a <=> $b} @_;
	my $len = @vals;
	if ($len%2) { # odd
		return $vals[int($len/2)];
		}
	else { # even
		return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
		}
	}

###################
#process MOL G+C character
###################
if ($character eq "%g+c") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "molg+c.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @molgcdata = split /\t/, $line;
		push (@taxlabels, $molgcdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", "@taxlabels", "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'mol %gc Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 mol_%g+c\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @molgctable = split /\t/, $line;
		push (my @taxnames, $molgctable[0]);
		push (my @molgcdata, $molgctable[1]);
		map {s/\smol%//g; } @molgcdata; #gets rid of the text mol% in the output data
		map {s/^,//; } @molgcdata; #gets rid of first comma
		map {s/,$//; } @molgcdata; #gets rid of last comma
		local $, = "";
	my @mediangc;
	my $temp;
	foreach $temp (@molgcdata) {
		if ($temp =~ /±.*,/) { # 
			map {s/,.*//g; } @molgcdata; # get rid of the additional values
			}
		if ($temp =~ /-.*,/) { # 
			map {s/,.*//g; } @molgcdata; # get rid of the additional values
			}
		if ($temp =~ /,/) { # there are multiple values of mol gc data
			map {s/,.*//g; } @molgcdata; # get rid of the additional values
			@mediangc = @molgcdata;
			}
		elsif ($temp =~ /-/) { # mol gc data exists as a range
			my @medianvals = split /-/, $temp;
			my $mval1 = $medianvals[0];
			my $temp1 = "$mval1" + 0.0; # ensures that the value is a digit
			my $mval2 = $medianvals[1];
			my $temp2 = "$mval2" + 0.0; # ensures that the value is a digit
			my @temp3 = ($temp1, $temp2);
			my @mval3 = median(@temp3);
			my $mval4 = "@mval3" + 0.0; # ensures that the value is a digit
			@mediangc = "$mval4" + 0.0; 
				}
		elsif ($temp =~ /±/){ # mol gc data exists as a median with error
			my @molgcdata = split /±/, $temp;
			my $temp2 =$molgcdata[0];
			$temp2 =~ s/[^[:print:]]//g; #delete non-printable chars
			@mediangc = "$temp2" + 0.0; # ensures that mediangc is a digit
				}
		else { # molgc single value 
			@mediangc = @molgcdata;
			next;
			}
		}
		print OUT @taxnames, " ", @mediangc, "\n"; # prints to $nexusoutfile molg+c.nex
	}
	print OUT "\n\n\;\n\nEND\;\n";

	print $outmessage;
	unlink $rawmatrix;
}

###################
#process CELL SHAPE character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "cell shape") {
	my $homout = "hom.cellshape.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @cellshapes, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @cellshapes;
		map {s/,bacilli,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rod shaped,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rodshaped,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rod-shaped,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rod-like,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,sub-coccoidal,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,short bacilli,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,coccobacillary,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,coccobacilli,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,ellipsoidal,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,oval,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,ovals,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,ovoid,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,short rods,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,subcoccoidal,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,subspherical,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,barrel-shaped,/,cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,barreliform,/,cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,cylinders,/,cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,club-shaped,/,diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,corneform,/,diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,coryneform,/,diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,diphtheroidal,/,diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,coiled,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,coils,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,spirals,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,helical,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,helical spirals,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,may form coils,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,helix,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,spiral,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,spirillum,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,wavy,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,vibrio,/,vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,crescent,/,vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,crescent-shaped,/,vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,crescent shaped,/,vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,crooked,/,curved,/g; } @cellshapes; # synonyms of spirillum
		map {s/,curved,/,curved,/g; } @cellshapes; # synonyms of spirillum
		map {s/,curved rods,/,vibrioid,/g; } @cellshapes; # synonyms of spirillum
		map {s/,coccoidal,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,coccoid,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,cocci,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,coccoid bodies,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,globose,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,isodiametric,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,lobed,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,lobes,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,multi-lobed,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,polygonal-rounded,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,semi-globose,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,spherical,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,sphaerical,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,sphere-shaped,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,sphere shaped,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,spheres,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,spheroidal,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,sphaeroidal,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,spheroids,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,sphaeroids,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,strongly lobed,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,squares,/,cuboidal,/g; } @cellshapes; # synonyms of cuboidal
		map {s/,square,/,cuboidal,/g; } @cellshapes; # synonyms of cuboidal
		map {s/,discs,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disc,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disc-shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disc shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,discoid,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disks,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disk,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disk-shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,disk shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,plate-shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,plate shaped,/,discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,hemisphaerical,/,hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,hemispheres,/,hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,hemisphaeres,/,hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,pleiomorphic,/,pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,pleiomorphism,/,pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,pleomorphism,/,pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,polymophic,/,pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,polymorphous,/,pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,pear shaped,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,flask-shaped,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,flask shaped,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,flasks,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,flask,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,flask-like,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,pears,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,pear,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,pear-like,/,pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,polygonal shape,/,polygonal,/g; } @cellshapes; # synonyms of polygonal
		map {s/,appendaged,/,prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,appendages,/,prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,stalks,/,prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,stalked,/,prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,pyramids,/,pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,triangles,/,pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,triangular,/,pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,ring-like,/,ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,horseshoe-shaped,/,ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,circle-shaped,/,ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,fusiforms,/,spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,fusiform,/,spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,swollen,/,spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,curved ends,/,rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,tapering ends,/,rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,tapered ends,/,rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,tapered,/,rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,filaments,/,filamentous,/g; } @cellshapes; # synonyms of filamentous
		map {s/,flexous,/,flexible,/g; } @cellshapes; # synonyms of filamentous
		map {s/,flexuous,/,flexible,/g; } @cellshapes; # synonyms of filamentous
		map {s/,branches,/,branching,/g; } @cellshapes; # synonyms of filamentous
		map {s/,branched,/,branching,/g; } @cellshapes; # synonyms of filamentous
		map {s/,short chains,/,chains,/g; } @cellshapes; 
		map {s/,irregular,/,not regular,/g; } @cellshapes; 
		map {s/,pointed,/,pointed ends,/g; } @cellshapes; 
		#negation of ontology terms and synonyms
		map {s/,not irregular,/,regular,/g; } @cellshapes; 
		map {s/,unbranched,/,not branching,/g; } @cellshapes; # synonyms of filamentous
		map {s/,not short chains,/,not chains,/g; } @cellshapes; 
		map {s/,not bacilli,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rod shaped,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rodshaped,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rod-shaped,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rod-like,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not sub-coccoidal,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not short bacilli,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not coccobacillary,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not coccobacilli,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not ellipsoidal,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not oval,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not ovoid,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not short rods,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not subcoccoidal,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not subspherical,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not barrel-shaped,/,not cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,not barreliform,/,not cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,not cylinders,/,not cylindrical,/g; } @cellshapes; # synonyms of cylindrical
		map {s/,not club-shaped,/,not diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,not corneform,/,not diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,not coryneform,/,not diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,not diphtheroidal,/,not diphtheroid,/g; } @cellshapes; # synonyms of diphtheroid
		map {s/,not coiled,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not coils,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not helical,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not may form coils,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not helix,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not spiral,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not spirillum,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not wavy,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not vibrio,/,not vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,not crooked,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not curved,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not curved rods,/,not vibrioid,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not coccoidal,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not coccoid,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not coccoid bodies,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not globose,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not isodiametric,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not lobed,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not lobes,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not multi-lobed,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not polygonal-rounded,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not semi-globose,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not spherical,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not sphaerical,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not sphere-shaped,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not sphere shaped,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not spheres,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not spheroidal,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not sphaeroidal,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not spheroids,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not sphaeroids,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not strongly lobed,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not squares,/,not cuboidal,/g; } @cellshapes; # synonyms of cuboidal
		map {s/,not square,/,not cuboidal,/g; } @cellshapes; # synonyms of cuboidal
		map {s/,not discs,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disc,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disc-shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disc shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not discoid,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disks,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disk,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disk-shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not disk shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not plate-shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not plate shaped,/,not discoidal,/g; } @cellshapes; # synonyms of discoid
		map {s/,not hemisphaerical,/,not hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,not hemispheres,/,not hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,not hemisphaeres,/,not hemispherical,/g; } @cellshapes; # synonyms of hemispherical
		map {s/,not pleiomorphic,/,not pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,not polymophic,/,not pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,not polymorphous,/,not pleomorphic,/g; } @cellshapes; # synonyms of pleomorphic
		map {s/,not pear shaped,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not flask-shaped,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not flask shaped,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not flasks,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not flask,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not flask-like,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not pears,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not pear,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not pear-like,/,not pear-shaped,/g; } @cellshapes; # synonyms of pear-shaped
		map {s/,not polygonal shape,/,not polygonal,/g; } @cellshapes; # synonyms of polygonal
		map {s/,not appendaged,/,not prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,not appendages,/,not prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,not stalks,/,not prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,not stalked,/,not prosthecate,/g; } @cellshapes; # synonyms of prosthecate
		map {s/,not pyramids,/,not pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,not triangles,/,not pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,not triangular,/,not pyramidal,/g; } @cellshapes; # synonyms of pyramidal
		map {s/,not ring-like,/,not ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,not horseshoe-shaped,/,not ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,not circle-shaped,/,not ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,not fusiforms,/,not spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,not fusiform,/,not spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,not swollen,/,not spindle-shaped,/g; } @cellshapes; # synonyms of spindle-shaped
		map {s/,not curved ends,/,not rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,not tapering ends,/,not rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,not tapered ends,/,not rounded ends,/g; } @cellshapes; # synonyms of rounded ends
		map {s/,not filaments,/,not filamentous,/g; } @cellshapes; # synonyms of filamentous
		#further homologize more specialized language
		map {s/,irregular shapes,/,not regular,/g; } @cellshapes; # synonyms of bacillus
		map {s/,gram rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,few bacilli,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,oval bacillus,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,plump rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,slender rod,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,longer rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rather long bacilli,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,slender rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,small single bacillus,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,thin rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,rounded bacillus,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,long rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,wide rods,/,bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,single small oval,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,oval bacillus,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,small oval,/,coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,small coccoid forms,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,coccoid forms,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,small ovali,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,single small oval,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,small single oval,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,mostly circle-shaped,/,ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,slightly tapering ends/,rounded ends,/g; } @cellshapes; #synonyms of rounded ends
		map {s/,spiral forms,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,long filaments,/,filamentous,/g; } @cellshapes; # synonyms of filaments
		map {s/,often curved,/,curved,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,curved bacteria,/,curved,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,curved arrangements,/,curved,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,curved slightly,/,curved,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,very irregular,/,not regular,/g; } @cellshapes; # synonyms of irregular
		map {s/,regular shape,/,regular,/g; } @cellshapes; # synonyms of regular
		map {s/,nearly straight,/,straight,/g; } @cellshapes; 
		map {s/,nonflexible,/,not flexible,/g; } @cellshapes; #
		map {s/,rather long bacilli,/,bacillus,/g; } @cellshapes; #
		map {s/,long rather bacilli,/,bacillus,/g; } @cellshapes; #
		map {s/,uniseriate filaments,/,filamentous,/g; } @cellshapes; #
		map {s/,longer filaments,/,filamentous,/g; } @cellshapes; #
		map {s/,rod-shaped organism,/,bacillus,/g; } @cellshapes; 
		map {s/,solid staining bacillus,/,bacillus,/g; } @cellshapes; 
		map {s/,straight nearly,/,straight,/g; } @cellshapes; 
		map {s/,rod-shaped organism,/,bacillus,/g; } @cellshapes; 
		map {s/,long relatively daughter trichomes,/,filamentous,/g; } @cellshapes; 
		#negation of further homologize more specialized language
		map {s/,not oval bacillus,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not plump rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not slender rod,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not slender rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not small single bacillus,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not thin rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not rounded bacillus,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not long rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not wide rods,/,not bacillus,/g; } @cellshapes; # synonyms of bacillus
		map {s/,not single small oval,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not oval bacillus,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not small oval,/,not coccobacillus,/g; } @cellshapes; # synonyms of coccobacillus
		map {s/,not small coccoid forms,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not coccoid forms,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not small ovali,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not single small oval,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not small single oval,/,not coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,not mostly circle-shaped,/,not ring-shaped,/g; } @cellshapes; # synonyms of ring-shaped
		map {s/,not slightly tapering ends/,not rounded ends,/g; } @cellshapes; #synonyms of rounded ends
		map {s/,not spiral forms,/,not spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not long filaments,/,not filamentous,/g; } @cellshapes; # synonyms of filaments
		map {s/,not often curved,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not curved bacteria,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not curved arrangements,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not curved slightly,/,straight,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not very irregular,/,regular,/g; } @cellshapes; # synonyms of spirillum
		map {s/,not short,/,long,/g; } @cellshapes; #
		#split combined characters
		map {s/,slender flexible rods,/,flexible,bacillus,/g; } @cellshapes; 
		map {s/,straight rods,/,straight,bacillus,/g; } @cellshapes; 
		map {s/,bent slightly rods,/,not straight,bacillus,/g; } @cellshapes; 
		map {s/,curved slightly rods,/,vibrioid,/g; } @cellshapes; 
		map {s/,irregular slightly sides,/,not regular,/g; } @cellshapes; 
		map {s/,pleomorphic rods,/,pleomorphic,bacillus,/g; } @cellshapes; 
		map {s/,small single oval bacillus,/,coccus,bacillus,/g; } @cellshapes; 
		map {s/,flexible rods,/,flexible,bacillus,/g; } @cellshapes; 
		map {s/,flexible often rods,/,flexible,bacillus,/g; } @cellshapes; 
		map {s/,single flexible rods,/,flexible,bacillus,/g; } @cellshapes; 
		map {s/,coiled filaments,/,spirillum,filamentous,/g; } @cellshapes; 
		map {s/,never filaments,/,not filamentous,/g; } @cellshapes; 
		map {s/,irregular rods,/,not regular,bacillus,/g; } @cellshapes; 
		map {s/,multicellular filaments,/,multicellular,filamentous,/g; } @cellshapes; 
		map {s/,sinuous filaments,/,flexible,filamentous,/g; } @cellshapes; 
		map {s/,long sinuous filaments,/,long,flexible,filamentous,/g; } @cellshapes; 
		map {s/,lateral branches,/,branching,/g; } @cellshapes; 
		map {s/,slender flexible slender and flexible rods,/,flexible,bacillus,/g; } @cellshapes; 
		map {s/,long relatively rods,/,long,bacillus,/g; } @cellshapes; 
		map {s/,long rod-shaped bacteria,/,long,bacillus,/g; } @cellshapes; 
		map {s/,non-motile straight rods,/,straight,bacillus,/g; } @cellshapes; 
		map {s/,undulating filaments,/,spirillum,filamentous,/g; } @cellshapes; 
		map {s/,aerobic rods,/,bacillus,/g; } @cellshapes; 
		#negation of split combined characters
		map {s/,not slender flexible rods,/,not flexible,bacillus,/g; } @cellshapes; 
		map {s/,not straight rods,/,curved,bacillus,/g; } @cellshapes; 
		map {s/,not curved slightly rods,/,straight,bacillus,/g; } @cellshapes; 
		map {s/,not irregular slightly sides,/,regular,/g; } @cellshapes; 
		map {s/,not pleomorphic rods,/,not pleomorphic,bacillus,/g; } @cellshapes; 
		map {s/,not small single oval bacillus,/,not coccus,bacillus,/g; } @cellshapes; 
		map {s/,not flexible rods,/,not flexible,bacillus,/g; } @cellshapes; 
		map {s/,not flexible often rods,/,not flexible,bacillus,/g; } @cellshapes; 
		map {s/,not single flexible rods,/,not flexible,bacillus,/g; } @cellshapes; 
		#delete characters that aren't cell shapes or that aren't very useful
		map {s/,pairs,/,/g; } @cellshapes; 
		map {s/,elongate,/,/g; } @cellshapes; 
		map {s/,elongated,/,/g; } @cellshapes; 
		map {s/,chains,/,/g; } @cellshapes; 
		map {s/,flat,/,/g; } @cellshapes; 
		map {s/,zoogloea,/,/g; } @cellshapes; 
		map {s/,not zoogloea,/,/g; } @cellshapes; 
		map {s/,flat colonies,/,/g; } @cellshapes; 
		map {s/,filamentous tufts,/,/g; } @cellshapes; 
		map {s/,aerial conidia,/,/g; } @cellshapes; 
		map {s/,mostly very short,/,/g; } @cellshapes; 
		map {s/,not resting stages,/,/g; } @cellshapes; 
		map {s/,numerous filamentous tufts,/,/g; } @cellshapes; 
		map {s/,slender,/,/g; } @cellshapes; 
		map {s/,buds,/,/g; } @cellshapes; 
		map {s/,spherical slime layer,/,/g; } @cellshapes; 
		map {s/,rounded,/,/g; } @cellshapes; 
		map {s/,swellings,/,/g; } @cellshapes; 
		map {s/,singly,/,/g; } @cellshapes; 
		map {s/,tetrads,/,/g; } @cellshapes; 
		print OUT @taxlabels, @cellshapes, "\n"; # prints to $homout, hom.cellshape.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.cellshapes.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/cell shape//g; } @homcharlist; # gets rid of the label cell shape
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.cellshapes.txt
		}
#character discovery -sorts the characters and finds the unique characters
	sub uniq {
		my %seen;
		grep !$seen{$_}++, @_;
		}
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.cellshapes.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.cellshapes.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.cellshapes.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.cellshapes.txt
#temporarily rename charstates to label those that have been homologized with **		
	map {s/^bacillus/**bacillus/g; } @charlist2;
	map {s/^branching/**branching/g; } @charlist2;
	map {s/^coccobacillus/**coccobacillus/g; } @charlist2;
	map {s/^coccus/**coccus/g; } @charlist2;
	map {s/^curved/**curved/g; } @charlist2;
	map {s/^discoidal/**discoidal/g; } @charlist2;
	map {s/^cuboidal/**cuboidal/g; } @charlist2;
	map {s/^filamentous/**filamentous/g; } @charlist2;
	map {s/^flexible/**flexible/g; } @charlist2;
	map {s/^regular/**regular/g; } @charlist2;
	map {s/^long/**long/g; } @charlist2;
	map {s/^multicellular/**multicellular/g; } @charlist2;
	map {s/^pleomorphic/**pleomorphic/g; } @charlist2;
	map {s/^pointed ends/**pointed ends/g; } @charlist2;
	map {s/^prosthecate/**prosthecate/g; } @charlist2;
	map {s/^regular/**regular/g; } @charlist2;
	map {s/^ring-shaped/**ring-shaped/g; } @charlist2;
	map {s/^rounded ends/**rounded ends/g; } @charlist2;
	map {s/^short/**short/g; } @charlist2;
	map {s/^spindle-shaped/**spindle-shaped/g; } @charlist2;
	map {s/^spirillum/**spirillum/g; } @charlist2;
	map {s/^straight/**straight/g; } @charlist2;
	map {s/^trichomes/**trichomes/g; } @charlist2;
	map {s/^unicellular/**unicellular/g; } @charlist2;
	map {s/^vibrioid/**vibrioid/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.cellshapes.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.cellshapes.txt
		}
	my $temp6 = "temp6.cellshapes.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.cellshapes.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "cellshapes.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @cellshapedata = split (/\t/, $line);
		push (@taxnames, $cellshapedata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Shape Matrix'\;\n\tDIMENSIONS NCHAR=19\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'branching shape' \/  not_branching branching, ";
	print OUT "2 'shape flexibility' \/  'not flexible' flexible, ";
	print OUT "3 'shape regularity' \/  regular irregular, ";
	print OUT "4 'cell size' \/  short long, ";
	print OUT "5 multicellularity \/  unicellular multicellular, ";
	print OUT "6 'cell end shape' \/  'rounded ends' 'pointed ends', ";
	print OUT "7 'shape angularity' \/  'not curved' curved, ";
	print OUT "8 'bacillus shape' \/  'not bacillus' bacillus, ";
	print OUT "9 'vibrioid shape' \/  'not vibrioid' vibrioid, ";
	print OUT "10 'coccobacillus shape' \/  'not coccobacillus' coccobacillus, ";
	print OUT "11 'coccus shape' \/  'not coccus' coccus, ";
	print OUT "12 'spirillum shape' \/  'not spirillum' spirillum, ";
	print OUT "13 'discoidal shape' \/  'not discoidal' discoidal, ";
	print OUT "14 'spindle-shaped' \/  'not spindle-shaped' 'spindle-shaped', ";
	print OUT "15 'filamentous shape' \/  'not filamentous' filamentous, ";
	print OUT "16 'prosthecate shape' \/  'not prosthecate' prosthecate, ";
	print OUT "17 'ring-shaped' \/  'not ring-shaped' 'ring-shaped', ";
	print OUT "18 'cell shape' \/  bacillus vibrioid coccobacillus coccus spirillum discoidal spindle_shaped filamentous prosthecate ring_shaped";
	print OUT "19 cuboidal \/  'not cuboidal' cuboidal, ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @cellshapedata = split (/\t/, $line);
		push (my @taxnames, $cellshapedata[0]);

#code char 1 branching
		if ($line =~ /,branching,/) {
			print OUT @taxnames, "1";
			}
		else {
			print OUT @taxnames, "0";
		}
#code char 2 flexibility
		if ($line =~ /,flexible,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flexible,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 shape regularity
		if ($line =~ /,not regular,|,pleomorphic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,regular,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 cell size
		if ($line =~ /,short,|,not long,/) {
			print OUT "0";
			}
		elsif ($line =~ /,long,|,not short,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 5 multicellularity
		if ($line =~ /,multicellular,|,filamentous/) {
			print OUT "1";
			}
		elsif ($line =~ /,not multicellular,|,unicellular,/) {
			print OUT "1";
			}
		else {
			print OUT "0";
		}
#code char 6 end shape
		if ($line =~ /,rounded ends,|,not pointed ends,/) {
			print OUT "0";
			}
		elsif ($line =~ /,pointed ends,|,not rounded ends,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 7 shape angularity
		if ($line =~ /,straight,|,not curved,/) {
			print OUT "0";
			}
		elsif ($line =~ /,curved,|,vibrioid,|,spirillum,|,ring-shaped,|,not straight,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 8 bacillus shape
		if ($line =~ /,not bacillus,/) {
			print OUT "0";
			}
		elsif ($line =~ /,bacillus,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 9 vibrioid shape
		if ($line =~ /,not vibrioid,/) {
			print OUT "0";
			}
		elsif ($line =~ /,vibrioid,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 10 coccobacillus shape
		if ($line =~ /,not coccobacillus,/) {
			print OUT "0";
			}
		elsif ($line =~ /,coccobacillus,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 11 coccus shape
		if ($line =~ /,not coccus,/) {
			print OUT "0";
			}
		elsif ($line =~ /,coccus,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 12 spirillum shape
		if ($line =~ /,not spirillum,/) {
			print OUT "0";
			}
		elsif ($line =~ /,spirillum,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 13 discoidal shape
		if ($line =~ /,not discoidal,/) {
			print OUT "0";
			}
		elsif ($line =~ /,discoidal,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 14 spindle-shape shape
		if ($line =~ /,not spindle-shaped,/) {
			print OUT "0";
			}
		elsif ($line =~ /,spindle-shaped,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 15 filamentous shape
		if ($line =~ /,not filamentous,/) {
			print OUT "0";
			}
		elsif ($line =~ /,filamentous,|,trichomes,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 16 prosthecate shape
		if ($line =~ /,not prosthecate,/) {
			print OUT "0";
			}
		elsif ($line =~ /,prosthecate,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 17 ring-shaped shape
		if ($line =~ /,not ring-shaped,/) {
			print OUT "0";
			}
		elsif ($line =~ /,ring-shaped,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
#code char 18 morphology as a complex trait
		if ($line =~ /,spirillum.*filamentous.*vibrioid.*coccus,/) {
			print OUT "(1 3 4 7)";
			}
		elsif ($line =~ /,ring-shaped.*bacillus.*spirillum.*filamentous,/) {
			print OUT "(0 4 7 9)";
			}
		elsif ($line =~ /,filamentous.*spirillum.*coccus.*bacillus,/) {
			print OUT "(0 3 4 7)";
			}
		elsif ($line =~ /,bacillus.*spirillum.*filamentous,|,bacillus.*filamentous.*spirillum,|,spirillum.*bacillus.*filamentous,|,spirillum.*filamentous.*bacillus,|,filamentous.*spirillum.*bacillus,|,filamentous.*bacillus.*spirillum,/) {
			print OUT "(0 4 7)";
			}
		elsif ($line =~ /,bacillus.*vibrioid,|,vibrioid.*bacillus,/) {
			print OUT "(0 1)";
			}
		elsif ($line =~ /,bacillus.*coccobacillus,|,coccobacillus.*bacillus,/) {
			print OUT "(0 2)";
			}
		elsif ($line =~ /,bacillus.*coccus,|,coccus.*bacillus,/) {
			print OUT "(0 3)";
			}
		elsif ($line =~ /,coccobacillus.*coccus,|,coccus.*coccobacillus,/) {
			print OUT "(2 3)";
			}
		elsif ($line =~ /,bacillus.*filamentous,|,filamentous.*bacillus,/) {
			print OUT "(0 7)";
			}
		elsif ($line =~ /,bacillus.*spindle-shaped,|,spindle-shaped.*bacillus,/) {
			print OUT "(0 6)";
			}
		elsif ($line =~ /,filamentous.*spirillum,|,spirillum.*filamentous,/) {
			print OUT "(4 7)";
			}
		elsif ($line =~ /,bacillus.*coccobacillus.*filamentous,|,bacillus.*filamentous.*coccobacillus,|,coccobacillus.*bacillus.*filamentous,|,coccobacillus.*filamentous.*bacillus,|,filamentous.*coccobacillus.*bacillus,|,filamentous.*bacillus.*coccobacillus,/) {
			print OUT "(0 2 7)";
			}
		elsif ($line =~ /,bacillus.*discoidal,|,discoidal.*bacillus,/) {
			print OUT "(0 5)";
			}
		elsif ($line =~ /,bacillus.*ring-shaped,|,ring-shaped.*bacillus,/) {
			print OUT "(0 9)";
			}
		elsif ($line =~ /,bacillus,/) {
			print OUT "0";
			}
		elsif ($line =~ /,vibrioid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,coccobacillus,/) {
			print OUT "2";
			}
		elsif ($line =~ /,coccus,/) {
			print OUT "3";
			}
		elsif ($line =~ /,spirillum,/) {
			print OUT "4";
			}
		elsif ($line =~ /,discoidal,/) {
			print OUT "5";
			}
		elsif ($line =~ /,spindle-shaped,/) {
			print OUT "6";
			}
		elsif ($line =~ /,filamentous,|,trichomes,/) {
			print OUT "7";
			}
		elsif ($line =~ /,prosthecate,/) {
			print OUT "8";
			}
		elsif ($line =~ /,ring-shaped,/) {
			print OUT "9";
			}
		else {
			print OUT "?";
		}
#code char 19 cuboidal shape
		if ($line =~ /,not cuboidal,/) {
			print OUT "0";
			}
		elsif ($line =~ /,cuboidal,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}
		print OUT "\n";
	}
	print OUT "\n\;\n\nEND\;\n";
	print $outmessage;

	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}

###################
# cell length
###################
elsif ($character eq "cell length") {
#prepare nexus file
	my @celllengthdata;
	my @taxlabels;
	my $nexusoutfile = "celllength.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@celllengthdata = split (/\t/, $line);
		push (@taxlabels, $celllengthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Length Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'cell length'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @celllengthtable = split (/\t/, $line);
		push (my @taxnames, $celllengthtable[0]);
		push (my @celllengthdata, $celllengthtable[1]);
		map {s/\p{Greek}m//gi; } @celllengthdata; #gets rid of 'µm'
		map {s/nm//gi; } @celllengthdata; #gets rid of 'nm', which is a pdf-to-text conversion error
		map {s/pm//gi; } @celllengthdata; #gets rid of 'pm', which is a pdf-to-text conversion error
		map {s/\s//gi; } @celllengthdata; #gets rid of 'space'
		map {s/^,//; } @celllengthdata; #gets rid of first comma in the output data
		map {s/,$//; } @celllengthdata; #gets rid of last comma in the output data
	my @mediancelllength;
	my $temp;
	foreach $temp (@celllengthdata) {
		if ($temp =~ /±.*,/) { # 
			map {s/,.*//g; } @celllengthdata; # get rid of the additional values
			}
		if ($temp =~ /-.*,/) { # 
			map {s/,.*//g; } @celllengthdata; # get rid of the additional values
			}
		if ($temp =~ /,/) { # there are multiple values of mol gc data
			map {s/,.*//g; } @celllengthdata; # get rid of the additional values
			@mediancelllength = @celllengthdata;
			}
		elsif ($temp =~ /-/) { # mol gc data exists as a range
			$temp =~ s/[^[:print:]]//g; #delete non-printable chars
			my @medianvals = split /-/, $temp;
			my $mval1 = $medianvals[0];
			my $temp1 = "$mval1" + 0.0; # ensures that the value is a digit
			my $mval2 = $medianvals[1];
			my $temp2 = "$mval2" + 0.0; # ensures that the value is a digit
			my @temp3 = ($temp1, $temp2);
			my @mval3 = median(@temp3);
			my $mval4 = "@mval3" + 0.0; # ensures that the value is a digit
			@mediancelllength = "$mval4" + 0.0; 
				}
		elsif ($temp =~ /±/){ # mol gc data exists as a median with error
			my @celllengthdata = split /±/, $temp;
			my $temp2 =$celllengthdata[0];
			$temp2 =~ s/[^[:print:]]//g; #delete non-printable chars
			@mediancelllength = "$temp2" + 0.0; # ensures that mediangc is a digit
				}
		else { # molgc single value 
			$temp =~ s/[^[:print:]]//g; #delete non-printable chars
			@mediancelllength = $temp;
			next;
			}
		}
		local $, = "";
		print OUT @taxnames, " ", @mediancelllength, "\n"; # prints to $nexusoutfile meancelllength.nex
	}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}

###################
# cell width
###################
elsif ($character eq "cell width") {
#prepare nexus file
	my @cellwidthdata;
	my @taxlabels;
	my $nexusoutfile = "cellwidth.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@cellwidthdata = split (/\t/, $line);
		push (@taxlabels, $cellwidthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Median Cell Width Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'median cell width'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @cellwidthtable = split (/\t/, $line);
		push (my @taxnames, $cellwidthtable[0]);
		push (my @cellwidthdata, $cellwidthtable[1]);
		map {s/\p{Greek}m//gi; } @cellwidthdata; #gets rid of 'µm'
		map {s/nm//gi; } @cellwidthdata; #gets rid of 'nm', which is a pdf-to-text conversion error
		map {s/pm//gi; } @cellwidthdata; #gets rid of 'pm', which is a pdf-to-text conversion error
		map {s/\s//gi; } @cellwidthdata; #gets rid of 'space'
		map {s/^,//; } @cellwidthdata; #gets rid of first comma in the output data
		map {s/,$//; } @cellwidthdata; #gets rid of last comma in the output data
	my @mediancellwidth;
	my $temp;
	foreach $temp (@cellwidthdata) {
		if ($temp =~ /±.*,/) { # 
			map {s/,.*//g; } @cellwidthdata; # get rid of the additional values
			}
		if ($temp =~ /-.*,/) { # 
			map {s/,.*//g; } @cellwidthdata; # get rid of the additional values
			}
		if ($temp =~ /,/) { # there are multiple values of mol gc data
			map {s/,.*//g; } @cellwidthdata; # get rid of the additional values
			@mediancellwidth = @cellwidthdata;
			}
		elsif ($temp =~ /-/) { # mol gc data exists as a range
			$temp =~ s/[^[:print:]]//g; #delete non-printable chars
			my @medianvals = split /-/, $temp;
			my $mval1 = $medianvals[0];
			my $temp1 = "$mval1" + 0.0; # ensures that the value is a digit
			my $mval2 = $medianvals[1];
			my $temp2 = "$mval2" + 0.0; # ensures that the value is a digit
			my @temp3 = ($temp1, $temp2);
			my @mval3 = median(@temp3);
			my $mval4 = "@mval3" + 0.0; # ensures that the value is a digit
			@mediancellwidth = "$mval4" + 0.0; 
			}
		elsif ($temp =~ /±/) { # mol gc data exists as a median with error
			my @cellwidthdata = split /±/, $temp;
			my $temp2 =$cellwidthdata[0];
			$temp2 =~ s/[^[:print:]]//g; #delete non-printable chars
			@mediancellwidth = "$temp2" + 0.0; # ensures that mediangc is a digit
			}
		elsif ($temp =~ /[\d|\.]/) { # molgc single value 
			$temp =~ s/[^[:print:]]//g; #delete non-printable chars
			my $temp2 = "$temp" + 0.0;
			@mediancellwidth = $temp2;
			next;
			}
		else { # 
			next;
			}
		}
		local $, = "";
		print OUT @taxnames, " ", @mediancellwidth, "\n"; # prints to $nexusoutfile meancellwidth.nex
	}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}


###################
#mean cell length
###################
elsif ($character eq "mean cell length") {
#prepare nexus file
	my @meancelllengthdata;
	my @taxlabels;
	my $nexusoutfile = "meancelllength.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@meancelllengthdata = split (/\t/, $line);
		push (@taxlabels, $meancelllengthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Mean Cell Length Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'mean cell length'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @meancelllengthtable = split (/\t/, $line);
		push (my @taxnames, $meancelllengthtable[0]);
		push (my @meancelllengthdata, $meancelllengthtable[1]);
		map {s/,//g; } @meancelllengthdata;
		local $, = "";
		print OUT @taxnames, " ", @meancelllengthdata, "\n"; # prints to $nexusoutfile meancelllength.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#max cell length
###################
elsif ($character eq "max cell length") {
#prepare nexus file
	my @maxcelllengthdata;
	my @taxlabels;
	my $nexusoutfile = "maxcelllength.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@maxcelllengthdata = split (/\t/, $line);
		push (@taxlabels, $maxcelllengthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Max Cell Length Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'max cell length'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @maxcelllengthtable = split (/\t/, $line);
		push (my @taxnames, $maxcelllengthtable[0]);
		push (my @maxcelllengthdata, $maxcelllengthtable[1]);
		map {s/,//g; } @maxcelllengthdata;
		local $, = "";
		print OUT @taxnames, " ", @maxcelllengthdata, "\n"; # prints to $nexusoutfile maxcelllengthdata.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#mean cell width
###################
elsif ($character eq "mean cell width") {
#prepare nexus file
	my @meancellwidthdata;
	my @taxlabels;
	my $nexusoutfile = "meancellwidth.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@meancellwidthdata = split (/\t/, $line);
		push (@taxlabels, $meancellwidthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Mean Cell Width Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'mean cell width'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @meancellwidthtable = split (/\t/, $line);
		push (my @taxnames, $meancellwidthtable[0]);
		push (my @meancellwidthdata, $meancellwidthtable[1]);
		map {s/,//g; } @meancellwidthdata;
		local $, = "";
		print OUT @taxnames, " ", @meancellwidthdata, "\n"; # prints to $nexusoutfile meancellwidthdata.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#max cell width
###################
elsif ($character eq "max cell width") {
#prepare nexus file
	my @maxcellwidthdata;
	my @taxlabels;
	my $nexusoutfile = "maxcellwidth.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		@maxcellwidthdata = split (/\t/, $line);
		push (@taxlabels, $maxcellwidthdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Max Cell Width Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'max cell width'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @maxcellwidthtable = split (/\t/, $line);
		push (my @taxnames, $maxcellwidthtable[0]);
		push (my @maxcellwidthdata, $maxcellwidthtable[1]);
		map {s/,//g; } @maxcellwidthdata;
		local $, = "";
		print OUT @taxnames, " ", @maxcellwidthdata, "\n"; # prints to $nexusoutfile maxcellwidthdata.nex
		}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#cell relationships
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "cell relationships&aggregations") {
	my $homout = "hom.cellrelates.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @cellrelates, $columns[1]);
		#ontology terms and synonyms
		map {s/,clumps,/,clusters,/g; } @cellrelates; # synonyms of filamentous
		map {s/,filaments,/,filamentous,/g; } @cellrelates; # synonyms of filamentous
		map {s/,singly,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
		map {s/,single,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
		map {s/,unicells,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
		map {s/,streptobacilli,/,streptobacillus,/g; } @cellrelates; # synonyms of streptobacillus
		map {s/,streptococci,/,streptococcus,/g; } @cellrelates; # synonyms of streptococcus
		map {s/,chroococcoid cell clusters arise,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,chroococcoid cell clusters,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,chroococcoid stage,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,chroococcoidal,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,known chroococcal stages,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,chroococcal,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,chsroococcoid stage,/,chroococcoid,/g; } @cellrelates; # synonyms of chroococcoid
		map {s/,sarcinae,/,sarcina,/g; } @cellrelates; # synonyms of sarcina
		map {s/,sarcinoid,/,sarcina,/g; } @cellrelates; # synonyms of sarcina
		map {s/,staphylococci,/,staphylococcus,/g; } @cellrelates; # synonyms of staphylococcus
		map {s/,uniserial,/,uniseriate,/g; } @cellrelates; # synonyms of uniseriate
		map {s/,biserial,/,biseriate,/g; } @cellrelates; # synonyms of biseriate
		map {s/,multi-seriate,/,multiseriate,/g; } @cellrelates; # synonyms of multiseriate
		map {s/,multi-trichomous,/,multiseriate,/g; } @cellrelates; # synonyms of multiseriate
		map {s/,multiserial,/,multiseriate,/g; } @cellrelates; # synonyms of multiseriate
		map {s/,multitrichomous,/,multiseriate,/g; } @cellrelates; # synonyms of multiseriate
		map {s/,polyseriate,/,multiseriate,/g; } @cellrelates; # synonyms of multiseriate
		map {s/,diplobacilli,/,diplobacillus,/g; } @cellrelates; # synonyms of diplobacillus
		map {s/,diplococci,/,diplococcus,/g; } @cellrelates; # synonyms of diplococcus
		map {s/,dumbbell,/,diplococcus,/g; } @cellrelates; # synonyms of diplococcus
		map {s/,dumbbell-shaped,/,diplococcus,/g; } @cellrelates; # synonyms of diplococcus
		map {s/,palisades,/,palisade,/g; } @cellrelates; # synonyms of palisade
		map {s/,tetrads,/,tetrad,/g; } @cellrelates; # synonyms of tetrad
		print OUT @taxlabels, @cellrelates, "\n"; # prints to $homout, hom.cellrelates.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.cellrelates.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/cell relationships&aggregations//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.cellrelates.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.cellrelates.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.cellrelates.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.cellrelates.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.cellrelates.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/clusters/**clusters/g; } @charlist2;
	map {s/chains/**chains/g; } @charlist2;
	map {s/chroococcoid/**chroococcoid/g; } @charlist2;
	map {s/diplobacillus/**diplobacillus/g; } @charlist2;
	map {s/diplococcus/**diplococcus/g; } @charlist2;
	map {s/filamentous/**filamentous/g; } @charlist2;
	map {s/multicellular/**multicellular/g; } @charlist2;
	map {s/multiseriate/**multiseriate/g; } @charlist2;
	map {s/pairs/**pairs/g; } @charlist2;
	map {s/palisade/**palisade/g; } @charlist2;
	map {s/sarcina/**sarcina/g; } @charlist2;
	map {s/staphylobacillus/**staphylobacillus/g; } @charlist2;
	map {s/staphylococcus/**staphylococcus/g; } @charlist2;
	map {s/steptococcus/**steptococcus/g; } @charlist2;
	map {s/streptobacillus/**streptobacillus/g; } @charlist2;
	map {s/streptococcus/**streptococcus/g; } @charlist2;
	map {s/tetrad/**tetrad/g; } @charlist2;
	map {s/unicellular/**unicellular/g; } @charlist2;
	map {s/uniseriate/**uniseriate/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.cellrelates.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.cellrelates.txt
		}
	my $temp6 = "temp6.cellrelates.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.cellrelates.txt
		}
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "cellrelates.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @cellrelatesdata = split (/\t/, $line);
		push (@taxnames, $cellrelatesdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Relationships Matrix'\;\n\tDIMENSIONS NCHAR=18\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 unicellularity \/  'not unicellular' unicellular, ";
	print OUT "2 multicellularity \/  'not multicellular' multicellular, ";
	print OUT "3 filamentous \/  'not filamentous' filamentous, ";
	print OUT "4 streptobacillus \/  'not streptobacillus' streptobacillus, ";
	print OUT "5 streptococcus \/  'not streptococcus' streptococcus, ";
	print OUT "6 chroococcoid \/  'not chroococcoid' chroococcoid, ";
	print OUT "7 sarcina \/  'not sarcina' sarcina, ";
	print OUT "8 staphylococcus \/  'not staphylococcus' staphylococcus, ";
	print OUT "9 uniseriate \/  'not uniseriate' uniseriate, ";
	print OUT "10 multiseriate \/  'not multiseriate' multiseriate, ";
	print OUT "11 diplobacillus \/  'not diplobacillus' diplobacillus, ";
	print OUT "12 diplococcus \/  'not diplococcus' diplococcus, ";
	print OUT "13 palisades \/  'not palisades' palisades, ";
	print OUT "14 tetrads \/  'not tetrads' tetrads, ";
	print OUT "15 pairs \/  'not pairs' pairs, ";
	print OUT "16 cell_aggregation \/  unicellular pairs tetrads streptos multicellular, ";
	print OUT "17 chains \/  'not chains' chains, ";
	print OUT "18 clusters \/  'not clusters' clusters, ";

	print OUT " \;\n\tMATRIX\n";

	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @cellrelatesdata = split (/\t/, $line);
		push (my @taxnames, $cellrelatesdata[0]);
#code char 1 unicellularity
		if ($line =~ /,unicellular,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not unicellular,|,clusters,|,multicellular,|,filamentous,|,streptobacillus,|,streptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 multicellularity
		if ($line =~ /,multicellular,|,clusters,|,filamentous,|,streptobacillus,|,steptococcus,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplococcus,|,diplobacillus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,unicellular,|,not multicellular,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 filamentous
		if ($line =~ /,filamentous,|,uniseriate,|,multiseriate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not filamentous,|,unicellular,|,streptobacillus,|,streptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 streptobacillus
		if ($line =~ /,streptobacillus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not streptobacillus,|,filamentous,|,streptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 streptococcus
		if ($line =~ /,streptococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not streptococcus,|,filamentous,|,streptobacillus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 chroococcoid
		if ($line =~ /,chroococcoid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chroococcoid,|,pairs,|,filamentous,|,streptobacillus,|,streptococcus,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 sarcina
		if ($line =~ /,sarcina,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sarcina,|,unicellular,|,chroococcoid,|,filamentous,|,streptobacillus,|,steptococcus,|,chroococcoid,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 staphylococcus
		if ($line =~ /,staphylococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not staphylococcus,|,unicellular,|,chroococcoid,|,filamentous,|,streptobacillus,|,streptococcus,|,sarcina,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 uniseriate
		if ($line =~ /,uniseriate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not uniseriate,|,unicellular,|,streptobacillus,|,streptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 multiseriate
		if ($line =~ /,multiseriate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not multiseriate,|,unicellular,|,streptobacillus,|,steptococcus,|,chroococcoid,|,staphylococcus,|,uniseriate,|,sarcina,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 diplobacillus
		if ($line =~ /,diplobacillus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diplobacillus,|,filamentous,|,streptobacillus,|,steptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 diplococcus
		if ($line =~ /,diplococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diplococcus,|,filamentous,|,streptobacillus,|,steptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,palisade,|,tetrad,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 13 palisade
		if ($line =~ /,palisade,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not palisade,|,filamentous,|,streptobacillus,|,steptococcus,|,chroococcoid,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplococcus,|,diplobacillus,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 14 tetrad
		if ($line =~ /,tetrad,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetrad,|,unicellular,|,filamentous,|,streptobacillus,|,steptococcus,|,uniseriate,|,multiseriate,|,palisade,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 pairs
		if ($line =~ /,pairs,/) {
			print OUT "1";
			}
		elsif ($line =~ /,filamentous,|,streptobacillus,|,steptococcus,|,uniseriate,|,sarcina,|,tetrad,|,staphylococcus,|,staphylobacillus,|,multiseriate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 cell aggregation
		if ($line =~ /,unicellular,/) {
			print OUT "0";
			}
		elsif ($line =~ /,pairs,|,diplobacillus,|,diplococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,tetrad,/) {
			print OUT "2";
			}
		elsif ($line =~ /,streptobacillus,|,streptococcus,|,palisade,/) {
			print OUT "3";
			}
		elsif ($line =~ /,multicellular,|,clusters,|,filamentous,|,sarcina,|,multiseriate,|,uniseriate,|,chroococcoid,/) {
			print OUT "4";
			}
		else {
			print OUT "?";
		}
#code char 17 chains
		if ($line =~ /,chains,|,streptobacillus,|,streptococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chains,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 18 clusters
		if ($line =~ /,clusters,|,chroococcoid,|,sarcina,|,staphylococcus,/) {
			print OUT "1";
			}
		elsif ($line =~ /,filamentous,|,streptobacillus,|,uniseriate,|,tetrad,|,staphylobacillus,|,multiseriate,|,unicellular,|,pairs,|,diplobacillus,|,diplococcus,|,palisade,|,chains,|,streptococcus,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

		print OUT "\n";
		}
	print OUT "\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}	

###################
#process gram stain type character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "gram stain type") {
	my $homout = "hom.gramstain.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @gramstain, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @gramstain;
		map {s/,gram stain negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram staining negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram reaction negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gramnegative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram-staining-negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram-stain-negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram-negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative

		map {s/,gram stain positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram staining positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram reaction positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,grampositive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram-staining-positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram-stain-positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram-positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive

		map {s/,gram stain variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gram-variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gram staining variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gram reaction variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gramvariable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gram-staining-variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
		map {s/,gram indeterminate,/,gram variable,/gi; } @gramstain; # synonyms of gram variable

		map {s/,bi-polar gram stain,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,bipolar gram stain,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,bi-polar staining,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,bipolar staining,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,gram-staining-positive,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,bi polar staining,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain
		map {s/,not bi polar staining,/,bipolar gram stain,/gi; } @gramstain; # synonyms of bipolar gram stain

		print OUT @taxlabels, @gramstain, "\n"; # prints to $homout, hom.gramstain.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.gramstain.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/gram stain type//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.gramstain.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.gramstain.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.gramstain.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.gramstain.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.gramstain.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^bipolar gram stain/**bipolar gram stain/g; } @charlist2;
	map {s/^gram negative/**gram negative/g; } @charlist2;
	map {s/^gram positive/**gram positive/g; } @charlist2;
	map {s/^gram variable/**gram variable/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.gramstain.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.gramstain.txt
		}
	my $temp6 = "temp6.gramstain.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.gramstain.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "gramstain.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @gramstaindata = split (/\t/, $line);
		push (@taxnames, $gramstaindata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Gram Stain Matrix'\;\n\tDIMENSIONS NCHAR=5\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'gram staining' \/  'gram negative' 'gram positive' 'gram variable' bipolar_staining, ";
	print OUT "2 'gram negative' \/  'not gram negative' 'gram negative', ";
	print OUT "3 'gram positive' \/  'not gram positive' 'gram positive', ";
	print OUT "4 'gram variable' \/  'not gram variable' 'gram variable', ";
	print OUT "5 'bipolar gram stain' \/  'no bipolar staining' 'bipolar staining', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @gramstaindata = split (/\t/, $line);
		push (my @taxnames, $gramstaindata[0]);

#code char 1 gram stain
		if ($line =~ /,gram negative,/) {
			print OUT @taxnames, "0";
			}
		elsif ($line =~ /,gram positive,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,gram variable,/) {
			print OUT @taxnames, "2";
			}
		elsif ($line =~ /,bipolar gram stain,/) {
			print OUT @taxnames, "3";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 gram negative
		if ($line =~ /,gram negative,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gram negative,|,gram positive,|,gram variable,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 gram positive
		if ($line =~ /,gram positive,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gram positive,|,gram negative,|,gram variable,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 gram variable
		if ($line =~ /,gram variable,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gram variable,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 bipolar staining
		if ($line =~ /,bipolar gram stain,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not bipolar gram stain,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
		print OUT "\n";
		}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}

###################
#process internal features character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "internal features") {
	my $homout = "hom.intfeats.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @intfeats, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @intfeats;
		map {s/,central granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,course granules,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,fine granules,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,fine granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,globule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,granular inclusions,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,granular,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,granulated,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,intracellular globule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,large granules,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,large granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,numerous granules,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,single granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,singular granule,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,small electron-dense granules,/,granules,/gi; } @intfeats; # synonyms of 
		map {s/,solitary granules,/,granules,/gi; } @intfeats; # synonyms of 

		map {s/,metachromatic granules,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 
		map {s/,metachromatic granule,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 
		map {s/,polyphosphate granule,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 
		map {s/,polyphosphate,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 
		map {s/,volutin granules,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 
		map {s/,volutin granule,/,polyphosphate granules,/gi; } @intfeats; # synonyms of 

		map {s/,intracellular sulfur globule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,intracellular sulphur globule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulfur globule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulfur granule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulfur granule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulfur granules,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulfur granules,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulphur globule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulphur granule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulphur granule,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulphur granules,/,sulfur granules,/gi; } @intfeats; # synonyms of 
		map {s/,sulphur granules,/,sulfur granules,/gi; } @intfeats; # synonyms of 


		map {s/,phb,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate granules,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,polyhydroxybutyrate,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,polyhydroxybutyric acid,/,poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,starch,/,starch granules,/gi; } @intfeats; # synonyms of 
		map {s/,pha,/,poly-b-hydroxyalkanoate granules,/gi; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/gi; } @intfeats; # synonyms of 
		map {s/,poly-beta-hydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/gi; } @intfeats; # synonyms of 
		map {s/,polyhydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/gi; } @intfeats; # synonyms of 
		map {s/,glycogen granule,/,glycogen granules,/gi; } @intfeats; # synonyms of 
		map {s/,glycogen,/,glycogen granules,/gi; } @intfeats; # synonyms of 
		map {s/,glycogen particles,/,glycogen granules,/gi; } @intfeats; # synonyms of 
		map {s/,glycogen particle,/,glycogen granules,/gi; } @intfeats; # synonyms of 
		map {s/,cyanophycin,/,cyanophycin granules,/gi; } @intfeats; # synonyms of 

		map {s/,pairs gas vesicles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,gas vesicles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,vacuolated cells,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,vesicle-like microstructures,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,small vesicles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,small vacuole,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,small vacuoles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,small vesicles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,large gas vacuole,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,large vacuole,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,large vacuoles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,large vesciles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,gas vacuolation,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,vacuoles,/,gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,vesicular,/,gas vacuoles,/gi; } @intfeats; # synonyms of 

		map {s/,endospore,/,endospores,/gi; } @intfeats; # synonyms of 
		map {s/,resting cells,/,resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,fruiting bodies,/,resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,resting spores,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,spore-like structures,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,spore-former,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,spore-forming,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,sporulating,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,sporeforming,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,sporing,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,sporulating,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,sporeforming,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,heat-resistant endospores,/,endospores,/gi; } @intfeats; # synonyms of 
		map {s/,myxospore,/,myxospores,/gi; } @intfeats; # synonyms of 
		map {s/,spores,/,spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,endospore,/,endospores,/gi; } @intfeats; # synonyms of 

		map {s/,absent endospores,/,not endospores,/gi; } @intfeats; # synonyms of 
		map {s/,absent resting stages,/,not resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,absent spores,/,not endospores,/gi; } @intfeats; # synonyms of 
		map {s/,asporogenic,/,not resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,devoid of gas vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,devoid of vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,neither spores,/,not endospores,/gi; } @intfeats; # synonyms of 
		map {s/,neither spores,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non-spore-forming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non-sporeforming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non-sporing,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non-sporulating,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,nonspore-forming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,nonsporeforming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,nonsporing,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,nonsporulating,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,not fruiting bodies,/,not resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,not heat-resistant endospores,/,not endospores,/gi; } @intfeats; # synonyms of 
		map {s/,not resting cells,/,not resting stages,/gi; } @intfeats; # synonyms of 
		map {s/,not resting spores,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,not spore-like structures,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,not spores,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,not spore forming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,never spores,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non spore-forming,/,not spore-forming,/gi; } @intfeats; # synonyms of 
		map {s/,non-endospore-forming,/,not endospores,/gi; } @intfeats; # synonyms of 
		map {s/,not endospore,/,not endospores,/gi; } @intfeats; # synonyms of 

		map {s/,absent gas vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not gas vesicles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not pairs gas vesicles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not gas vesicles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not vacuolated cells,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not vesicle-like microstructures,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not small vesicles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not small vacuole,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not small vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not small vesicles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not large gas vacuole,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not large vacuole,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not large vacuoles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not large vesciles,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 
		map {s/,not poly-b-hydroxybutyrate,/,not poly-b-hydroxybutyrate granules,/gi; } @intfeats; # synonyms of 
		map {s/,not gas vacuolation,/,not gas vacuoles,/gi; } @intfeats; # synonyms of 

		map {s/,magnetite,/,magnetite,/gi; } @intfeats; # synonyms of 
		map {s/,pyrite,/,pyrite,/gi; } @intfeats; # synonyms of 

		map {s/,septa,/,septate,/gi; } @intfeats; # synonyms of 


		map {s/,xxx,/,xxx,/gi; } @intfeats; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @intfeats; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @intfeats; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @intfeats; # synonyms of 
#not internal features
		map {s/,carbon compound,/,/gi; } @intfeats; 
		map {s/,full,/,/gi; } @intfeats; 
		map {s/,solid,/,/gi; } @intfeats; 
		map {s/,no.*solid,/,/gi; } @intfeats; 
		map {s/,not agglutinate,/,/gi; } @intfeats; 
		map {s/,cellular,/,/gi; } @intfeats; 
		map {s/,complete,/,/gi; } @intfeats; 
		map {s/,complex,/,/gi; } @intfeats; 
		map {s/,no.*complex,/,/gi; } @intfeats; 
		map {s/,entire,/,/gi; } @intfeats; 
		map {s/,incomplete,/,/gi; } @intfeats; 
		map {s/,loose,/,/gi; } @intfeats; 
		map {s/,irregular,/,/gi; } @intfeats; 
		map {s/,open,/,/gi; } @intfeats; 
		map {s/,partial,/,/gi; } @intfeats; 
		map {s/,regular,/,/gi; } @intfeats; 
		map {s/,smooth,/,/gi; } @intfeats; 
		map {s/,.*smooth,/,/gi; } @intfeats; 
		map {s/,unbranched,/,/gi; } @intfeats; 
		map {s/,branching,/,/gi; } @intfeats; # synonyms of 
		map {s/,composite,/,/gi; } @intfeats; # synonyms of 
		map {s/,degenerate,/,/gi; } @intfeats; # synonyms of 
		map {s/,erose,/,/gi; } @intfeats; # synonyms of 
		map {s/,flagellate,/,/gi; } @intfeats; # synonyms of 
		map {s/,compound,/,/gi; } @intfeats; # synonyms of 
		map {s/,quinate,/,/gi; } @intfeats; # synonyms of 
		map {s/,sessile,/,/gi; } @intfeats; # synonyms of 
		map {s/,simple,/,/gi; } @intfeats; # synonyms of 
		map {s/,symmetry,/,/gi; } @intfeats; # synonyms of 

		map {s/,gram-positive,/,/gi; } @intfeats; # synonyms of 
		map {s/,C16,/,/gi; } @intfeats; # synonyms of 
		map {s/,singly,/,/gi; } @intfeats; # synonyms of 
		map {s/,major isoprenoid quinone,/,/gi; } @intfeats; # synonyms of 
		map {s/,xxx,/,/gi; } @intfeats; # synonyms of 
		map {s/,xxx,/,/gi; } @intfeats; # synonyms of 

		print OUT @taxlabels, @intfeats, "\n"; # prints to $homout, hom.intfeats.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.intfeats.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/internal features//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.intfeats.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.intfeats.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.intfeats.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.intfeats.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.intfeats.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^binary fission/**binary fission/g; } @charlist2;
	map {s/^cyanophycin granules/**cyanophycin granules/g; } @charlist2;
	map {s/^endospores/**endospores/g; } @charlist2;
	map {s/^gas vacuoles/**gas vacuoles/g; } @charlist2;
	map {s/^glycogen granules/**glycogen granules/g; } @charlist2;
	map {s/^granules/**granules/g; } @charlist2;
	map {s/^myxospores/**myxospores/g; } @charlist2;
	map {s/^poly-b-hydroxyalkanoate granules/**poly-b-hydroxyalkanoate granules/g; } @charlist2;
	map {s/^poly-b-hydroxybutyrate granules/**poly-b-hydroxybutyrate granules/g; } @charlist2;
	map {s/^polyphosphate granules/**polyphosphate granules/g; } @charlist2;
	map {s/^resting stages/**resting stages/g; } @charlist2;
	map {s/^spore-forming/**spore-forming/g; } @charlist2;
	map {s/^starch granules/**starch granules/g; } @charlist2;
	map {s/^sulfur granules/**sulfur granules/g; } @charlist2;

	map {s/^magnetite/**magnetite/g; } @charlist2;
	map {s/^pyrite/**pyrite/g; } @charlist2;
	map {s/^septate/**septate/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.intfeats.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.intfeats.txt
		}
	my $temp6 = "temp6.intfeats.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.intfeats.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "intfeats.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @intfeatsdata = split (/\t/, $line);
		push (@taxnames, $intfeatsdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Internal Features Matrix'\;\n\tDIMENSIONS NCHAR=14\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 resting_stages \/  'no resting stages' 'resting stages', ";
	print OUT "2 endospores \/  'no endospores' 'endospores present', ";
	print OUT "3 myxoxpores \/  'no myxoxpores' 'myxoxpores present', ";
	print OUT "4 'spore-forming' \/  'not spore-forming' 'spore-forming', ";
	print OUT "5 granules \/  'granules absent' 'granules present', ";
	print OUT "6 'starch granules' \/  'starch granules absent' 'starch granules present', ";
	print OUT "7 'poly-beta-hydroxybutyrate granules' \/  'PHB granules absent' 'PHB present', ";
	print OUT "8 'poly-beta-hydroxyalkanoate granules' \/  'PHA absent' 'PHA present', ";
	print OUT "9 'glycogen granules' \/  'glycogen absent' 'glycogen present', ";
	print OUT "10 'cyanophycin granules' \/  'cyanophycin absent' 'cyanophycin present', ";
	print OUT "11 'polyphosphate granules' \/  'polyphosphate absent' 'polyphosphate present', ";
	print OUT "12 'sulfur granules' \/  'sulfur granules absent' 'sulfur granules present', ";
	print OUT "13 'gas vacuoles' \/  'no gas vacuoles' 'gas vacuoles present', ";
	print OUT "14 'cell division' \/  'binary fission' present, ";
	print OUT "15 'iron minerals' \/  'no iron minerals' magnetite pyrite, ";
	print OUT "16 septation \/  'septation not noted' septate aseptate, ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @intfeatsdata = split (/\t/, $line);
		push (my @taxnames, $intfeatsdata[0]);

#code char 1 resting stages
		if ($line =~ /,resting stages,|,endospores,|,myxospores,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not resting stages,|,not spore-forming,|,not endospores,|,not myxospores,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
			}
#code char 2 endospores
		if ($line =~ /,endospores,|,spore-forming,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not endospores,|,not spore-forming,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 3 myxospores
		if ($line =~ /,myxospores,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not myxospores,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 4 spores
		if ($line =~ /,spore-forming,|,endospores,|,myxospores,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spore-forming,|not endospores,|,not myxospores,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 5 granules
		if ($line =~ /,granules,|,starch granules,|,poly-b-hydroxybutyrate granules,|,poly-b-hydroxyalkanoate granules,|,glycogen granules,|,cyanophycin granules,|,polyphosphate granules,|,sulfur granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 6 starch granules
		if ($line =~ /,starch granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not starch granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 7 poly-beta-hydroxybutyrate granules
		if ($line =~ /,poly-b-hydroxybutyrate granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not poly-b-hydroxybutyrate granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 8 poly-beta-hydroxyalkanoate granules
		if ($line =~ /,poly-b-hydroxyalkanoate granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not poly-b-hydroxyalkanoate granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 9 glycogen granules
		if ($line =~ /,glycogen granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycogen granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 10 cyanophycin granules
		if ($line =~ /,cyanophycin granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cyanophycin granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 11 polyphosphate granules
		if ($line =~ /,polyphosphate granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polyphosphate granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 12 sulfur granules
		if ($line =~ /,sulfur granules,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfur granules,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 13 gas vacuoles
		if ($line =~ /,gas vacuoles,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gas vacuoles,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 14 cell division
		if ($line =~ /,binary fission,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not binary fission,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 15 iron minerals
		if ($line =~ /,magnetite,/) {
			print OUT "1";
			}
		elsif ($line =~ /,pyrite,/) {
			print OUT "2";
			}
		else {
			print OUT "0";
			}
#code char 16 septation
		if ($line =~ /,septate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,aseptate,/) {
			print OUT "2";
			}
		else {
			print OUT "0";
			}

		print OUT "\n";
		}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}

###################
#process motility character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "motility") {
	my $homout = "hom.motility.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @motility, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @motility;
		map {s/,absent flagellar motility,/,not flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,absent flagellar,/,not flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,absent gliding motility,/,not gliding motility,/gi; } @motility; # synonyms of 
		map {s/,absent gliding,/,not gliding motility,/gi; } @motility; # synonyms of 
		map {s/,absent motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,actively motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,axial fibrils,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,axial filaments,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,bipolar polytrichous flagella,/,amphilophotrichous,/gi; } @motility; # synonyms of 
		map {s/,bundles of flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,devoid of flagellar motility,/,not flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,devoid of flagellar,/,not flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,endoflagella,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,endoflagellum,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,few flagella,/,flagella,not monotrichous,/gi; } @motility; # synonyms of 
		map {s/,flagellar flagellar,/,flagella,/gi; } @motility; # synonyms of 
		map {s/,generally non-motile,/,flagella,/gi; } @motility; # synonyms of 
		map {s/,gliding cells,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,gliding motilities,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,gliding types,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,gliding,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,gram-positive nonmotile nonmotile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,immotile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,lack flagella,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,lateral flagella,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,mainly non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,mature gliding stage,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,monopolar flagella,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,monopolar flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,monopolar polytrichous flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,monotrichous flagella,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,monotrichous monotrichous,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,motile arid motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,motile by gliding,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,motile cells,/,motile,/gi; } @motility; # synonyms of 
		map {s/,motile endospore-forming motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,motile motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,motile vegetative motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,motility,/,motile,/gi; } @motility; # synonyms of 
		map {s/,never motility,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,no.*motility,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,no.*non-motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,non motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-flagellated and non- motile,/,not flagella,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-flagellated and non-motile,/,not flagella,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-flagellated non-motile non-flagellated,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,non-flagellated non-motile,/,not flagella,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-flagellated,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,non-gliding,/,not gliding motility,/gi; } @motility; # synonyms of 
		map {s/,non-motile catalase-and non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-motile cells,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-motile non-gliding,/,not motile,not gliding motility,/gi; } @motility; # synonyms of 
		map {s/,non-motile non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-motile wide non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,non-spore-forming non-motile non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,nonflagellated,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,nonmotile nonmotile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,nonmotile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,nonmotile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,not flagellar flagellar,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,not flagellar,/,not flagella,/gi; } @motility; # synonyms of 
		map {s/,not gliding types,/,not gliding motility,/gi; } @motility; # synonyms of 
		map {s/,not motile motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,not motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,not motility,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,not non-motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,not swimming cells,/,not flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,numerous filamentous tufts,/,amphilophotrichous,/gi; } @motility; # synonyms of 
		map {s/,peculiar gliding motility,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,periplasmic fibrils,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,periplasmic flagella,/,periplasmic flagella,/gi; } @motility; # synonyms of 
		map {s/,peritrichous flagella,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,peritrichous flagella,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,peritrichously,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,polar bundle of flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar polar flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar polar flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar subpolar flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar tuft of flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,polar tufts of flagella,/,amphilophotrichous,/gi; } @motility; # synonyms of 
		map {s/,predominant motile motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,rare motile motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,retarded peritrichous flagella,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,semisolid motility,/,motile,/gi; } @motility; # synonyms of 
		map {s/,several flagella,/,flagella,not monotrichous,/gi; } @motility; # synonyms of 
		map {s/,several peritrichous flagella,/,peritrichous,/gi; } @motility; # synonyms of 
		map {s/,showing gliding motility,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,single flagella,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,single flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,single polar flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,single polar polar flagellum,/,monotrichous,/gi; } @motility; # synonyms of 
		map {s/,slightly motile,/,motile,/gi; } @motility; # synonyms of 
		map {s/,slow gliding motility,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,subpolar flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,swimming cells,/,flagellar motility,/gi; } @motility; # synonyms of 
		map {s/,terminal flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,tuft of polar flagella,/,lophotrichous,/gi; } @motility; # synonyms of 
		map {s/,tufts of polar flagella,/,amphilophotrichous,/gi; } @motility; # synonyms of 
		map {s/,tumbling,/,tumbling motility,/gi; } @motility; # synonyms of 
		map {s/,twitching motility,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,twitching,/,gliding motility,/gi; } @motility; # synonyms of 
		map {s/,usually non-motile,/,not motile,/gi; } @motility; # synonyms of 
		map {s/,usually nonmotile,/,not motile,/gi; } @motility; # synonyms of 

#not motility characters
		map {s/,filaments,/,/g; } @motility; # synonyms of 

		print OUT @taxlabels, @motility, "\n"; # prints to $homout, hom.motility.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.motility.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,motility,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.motility.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.motility.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.motility.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.motility.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.motility.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^amphilophotrichous/**amphilophotrichous/g; } @charlist2;
	map {s/^peritrichous/**peritrichous/g; } @charlist2;
	map {s/^amphitrichous/**amphitrichous/g; } @charlist2;
	map {s/^flagella/**flagella/g; } @charlist2;
	map {s/^flagellar motility/**flagellar motility/g; } @charlist2;
	map {s/^gliding motility/**gliding motility/g; } @charlist2;
	map {s/^lophotrichous/**lophotrichous/g; } @charlist2;
	map {s/^monotrichous/**monotrichous/g; } @charlist2;
	map {s/^motile/**motile/g; } @charlist2;
	map {s/^periplasmic flagella/**periplasmic flagella/g; } @charlist2;
	map {s/^polytrichous/**polytrichous/g; } @charlist2;
	map {s/^tumbling motility/**tumbling motility/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.motility.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.motility.txt
		}
	my $temp6 = "temp6.motility.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.motility.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "motility.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @motilitydata = split (/\t/, $line);
		push (@taxnames, $motilitydata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Motility Matrix'\;\n\tDIMENSIONS NCHAR=12\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 motility \/  'non-motile' motile, ";
	print OUT "2 flagellar_motility \/  'not flagellar motility' 'flagellar motility', ";
	print OUT "3 gliding_motility \/  'not gliding motility' 'gliding motility', ";
	print OUT "4 'tumbling motility' \/  'not tumbling motility' 'tumbling motility', ";
	print OUT "5 flagella \/  'no flagella' 'flagella', ";
	print OUT "6 monotrichous \/  'not monotrichous' monotrichous, ";
	print OUT "7 polytrichous \/  'not polytrichous' polytrichous, ";
	print OUT "8 lophotrichous \/  'not lophotrichous' lophotrichous, ";
	print OUT "9 amphitrichous \/  'not amphitrichous' amphitrichous, ";
	print OUT "10 amphilophotrichous \/  'not amphilophotrichous' amphilophotrichous, ";
	print OUT "11 'periplasmic flagella' \/  'periplasmic flagella absent' 'periplasmic flagella present', ";
	print OUT "12 peritrichous \/  'not peritrichous' peritrichous, ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @motilitydata = split (/\t/, $line);
		push (my @taxnames, $motilitydata[0]);

#code char 1 motility
		if ($line =~ /,motile,|,flagellar motility,|,gliding motility,|,tumbling motility,|,peritrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not motile,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
			}
#code char 2 flagellar motility
		if ($line =~ /,flagellar motility,|,peritrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flagellar motility,|,not flagella,|,not motile,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 3 gliding motility
		if ($line =~ /,gliding motility,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gliding motility,|,not motile,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 4 tumbling motility
		if ($line =~ /,tumbling motility,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tumbling motility,|,not motile,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 5 flagella
		if ($line =~ /,flagella,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,|,peritrichous/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flagella,|,not motile,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 6 monotrichous
		if ($line =~ /,monotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,|,peritrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 7 polytrichous
		if ($line =~ /,polytrichous,|,peritrichous,|,lophotrichous,|,amphilophotrichous,|,amphitrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polytrichous,|,monotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 8 lophotrichous
		if ($line =~ /,lophotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lophotrichous,|,monotrichous,|,polytrichous,|,amphitrichous,|,amphilophotrichous,|,peritrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 9 amphitrichous
		if ($line =~ /,amphitrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphitrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphilophotrichous,|,peritrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 10 amphilophotrichous
		if ($line =~ /,amphilophotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphilophotrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,not motile,|,peritrichous,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 11 periplasmic flagella
		if ($line =~ /,periplasmic flagella,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not periplasmic flagella,|,not motile,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 12 peritrichous
		if ($line =~ /,peritrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not peritrichous,|,monotrichous,|,lophotrichous,|,amphilophotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}


		print OUT "\n";
		}
	print OUT "\n\;\nEND\;\n";	
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}

###################
#process pigment compounds character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "pigment compounds") {
	my $homout = "hom.pigments.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @pigments, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @pigments;
		map {s/,diffusible pigment,/,pigment,diffusible,/gi; } @pigments; # synonyms of filamentous
		map {s/,no.*diffusible pigments,/,not pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,other pigmented enterococci,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,pigmentation,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,pigmented and unpigmented clones,/,pigment,not pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,yellow pigmentation,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,yellow pigmented,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,pigment compounds,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,multiple absorption maxima at.*,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous
		map {s/,xxx,/,pigment,/gi; } @pigments; # synonyms of filamentous

		print OUT @taxlabels, @pigments, "\n"; # prints to $homout, hom.pigments.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.pigments.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/pigment compounds//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.pigments.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.pigments.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.pigments.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.pigments.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.pigments.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/pigment/**pigment/g; } @charlist2;
	map {s/diffusible/**diffusible/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.pigments.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.pigments.txt
		}
	my $temp6 = "temp6.pigments.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.pigments.txt
		}
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "pigments.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @pigmentsdata = split (/\t/, $line);
		push (@taxnames, $pigmentsdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Pigment Matrix'\;\n\tDIMENSIONS NCHAR=2\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 pigments \/  'not pigments' pigments, ";
	print OUT "2 diffusibility \/  'not diffusible' diffusible, ";

	print OUT " \;\n\tMATRIX\n";

	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @pigmentsdata = split (/\t/, $line);
		push (my @taxnames, $pigmentsdata[0]);
#code char 1 pigments
		if ($line =~ /,pigments,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not pigments,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 diffusibility
		if ($line =~ /,diffusible,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not diffusible,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}

		print OUT "\n";
	}
	print OUT "\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}	

###################
#process salinity preference character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "salinity preference") {
	my $homout = "hom.salinityprefs.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @salinityprefs, $columns[1]);
		$_ = lc for @salinityprefs;
		#ontology terms and synonyms
		map {s/,Difco marine agar,/,marine,/gi; } @salinityprefs; # synonyms 
		map {s/,extremely halophilic,/,halophilic,/gi; } @salinityprefs; # synonyms
		map {s/,extremely halotolerant,/,halotolerant,/gi; } @salinityprefs; # synonyms
		map {s/,facultative halophile,/,halophilic,/gi; } @salinityprefs; # synonyms
		map {s/,marine environments,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,marine forms,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,marine lake,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,marine sand,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,marine solar saltern,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,marine sponge Plakortis simplex,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,moderately halophilic,/,halophilic,/gi; } @salinityprefs; # synonyms
		map {s/,moderately halotolerant,/,halotolerant,/gi; } @salinityprefs; # synonyms
		map {s/,no.* halotolerant,/,not halotolerant,/gi; } @salinityprefs; # synonyms
		map {s/,no.*marine bacteria,/,not marine,/gi; } @salinityprefs; # synonyms
		map {s/,no.*requires NaCl,/,not salt required,/gi; } @salinityprefs; # synonyms
		map {s/,no.*requires NaC1,/,not salt required,/gi; } @salinityprefs; # synonyms
		map {s/,requires NaCl,/,salt required,/gi; } @salinityprefs; # synonyms
		map {s/,requires NaC1,/,salt required,/gi; } @salinityprefs; # synonyms
		map {s/,no.*marine sediments,/,marine,/gi; } @salinityprefs; # synonyms
		map {s/,no.*moderately halophilic,/,not halophilic,/gi; } @salinityprefs; # synonyms
		map {s/,requires sodium,/,sodium required,/gi; } @salinityprefs; # synonyms
		map {s/,salinity preference,/,/gi; } @salinityprefs; # synonyms

		print OUT @taxlabels, @salinityprefs, "\n"; # prints to $homout, hom.salinityprefs.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.salinityprefs.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/salinity preference//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.salinityprefs.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.salinityprefs.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.salinityprefs.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.salinityprefs.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.salinityprefs.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/marine/**marine/g; } @charlist2;
	map {s/halophilic/**halophilic/g; } @charlist2;
	map {s/halotolerant/**halotolerant/g; } @charlist2;
	map {s/salt required/**salt required/g; } @charlist2;
	map {s/sodium required/**sodium required/g; } @charlist2;

	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.salinityprefs.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.salinityprefs.txt
		}
	my $temp6 = "temp6.salinityprefs.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.salinityprefs.txt
		}
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "salinityprefs.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @salinityprefsdata = split (/\t/, $line);
		push (@taxnames, $salinityprefsdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Pigment Matrix'\;\n\tDIMENSIONS NCHAR=7\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 marine \/  'not marine' marine, ";
	print OUT "2 halophilic \/  'not halophilic' halophilic, ";
	print OUT "3 halotolerant \/  'not halotolerant' halotolerant, ";
	print OUT "4 'salt required' \/  'no salt required' 'salt required', ";
	print OUT "5 'sodium required' \/  'no sodium required' 'sodium required', ";
	print OUT "6 'salt preference' \/  marine halotolerant halophilic, ";
	print OUT "7 'salt requirement' \/  'no salt required' 'sodium or salt required', ";

	print OUT " \;\n\tMATRIX\n";

	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @salinityprefsdata = split (/\t/, $line);
		push (my @taxnames, $salinityprefsdata[0]);
#code char 1 marine
		if ($line =~ /,marine,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not marine,|,halophilic,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 halophilic
		if ($line =~ /,halophilic,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not halophilic,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 3 halotolerant
		if ($line =~ /,halotolerant,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not halotolerant,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 4 salt required
		if ($line =~ /,salt required,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not salt required,|,not sodium required,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 5 sodium required
		if ($line =~ /,sodium required,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not sodium required,|,salt required,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 6 salt preference
		if ($line =~ /,halotolerant,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,halophilic,/) {
			print OUT @taxnames, "2";
			}
		elsif ($line =~ /,marine,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 7 salt requirement
		if ($line =~ /,salt required,|,sodium required,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not salt required,|,not sodium required,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}

		print OUT "\n";
	}
	print OUT "\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}	

###################
#process nacl minimum character
###################
elsif ($character eq "nacl minimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "naclmin.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @naclmindata = split (/\t/, $line);
		push (@taxlabels, $naclmindata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'NaCl Minimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'NaCl min'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @naclmintable = split (/\t/, $line);
		push (my @taxnames, $naclmintable[0]);
		push (my @naclmindata, $naclmintable[1]);
		local $, = "";
		map {s/,\d\smm,/,/g; } @naclmindata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\smmw\/v,/,/g; } @naclmindata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\.\d\smmw\/v,/,/g; } @naclmindata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\.\d\smm,/,/g; } @naclmindata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d.\d\d\smm,/,/g; } @naclmindata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\sdays,/,/gi; } @naclmindata; # not a salt min
		map {s/plus.*//gi; } @naclmindata; # get rid of any additional values
		map {s/^,//g; } @naclmindata; # get rid of first comma
		map {s/,$//g; } @naclmindata; # get rid of last comma
		map {s/,.*//g; } @naclmindata; # get rid of any additional values
		map {s/,nacl minimum,/,/gi; } @naclmindata; # not a salt min
		my $temp;
		#convert all concentration values to % NaCl and get rid of units
		foreach $temp (@naclmindata) {
			if ($temp =~ /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@naclmindata = "@val1" + 0.0;
				}	
			elsif ($temp =~ /0 nacl/i ) { # 
				@naclmindata = 0.0;
				}	
			elsif ($temp =~ /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@naclmindata = 0.05;
				}	
			elsif ($temp =~ /not required|absence of NaCl|without addition of NaCl/i ) {
				@naclmindata = 0.0;
				}	
			elsif ($temp =~ /sea.*salt/i) { # converts PERCENT SEA SALTS to % NaCl
				if ($temp =~ /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @naclmindata; # get rid of the first comma
					map {s/\,.*//; } @naclmindata; # get rid of any additional values
					map {s/%%.*//; } @naclmindata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmindata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif ($temp =~ /-/) { # % sea salt values with a range
					map {s/\,//; } @naclmindata; # get rid of the first comma
					map {s/\,.*//; } @naclmindata; # get rid of any additional values
					map {s/%.*//; } @naclmindata; # get rid of any additional values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/%.*//; } @mval2; # get rid of any additional values
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					my @mval3 = median(@temp3);
					my @mval4 = "@mval3" + 0.0;
					@naclmindata = "@mval4" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				else { # single value % sea salt value
					map {s/\,//; } @naclmindata; # get rid of the first comma
					map {s/\,.*//; } @naclmindata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmindata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif ($temp =~ /sea.*water/i) { #convert SEAWATER DILUTIONS to % Nacl
					if ($temp =~ /x|%/) { 
						if ($temp =~ /%/) { # when the seawater concentration reported as a per cent dilution
							my @naclvals = (split /%/, $temp);
							push (my @val1, $naclvals[0]);
							my @val2 = "@val1" + 0.0 ; # makes sure the string is a number
							my @val3 = @val2 / 100.0 ; # converts percent to a decimal
							@naclmindata = "@val3" * 3.5 ; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /x/i) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = (split /\s/, $temp);
							push (my @val1, $naclvals[0]);
							map {s/x.*//i; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@naclmindata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = (split /\//, $temp);
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0; #ensure that the numerator in the string is a nmber
							map {s/x.*//i; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0; #ensure that the denominator in the string is a number
							my @temp3 = @temp1 / @temp2;
							@naclmindata = "@temp3" * 3.5 ;
							}
						}
					else { # has elsif (@naclmindata = /^seawater$/) # assume it is pure seawater, undiluted seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@naclmindata = "@val1" + 0.0;
						}
				}
			elsif ($temp =~ /%.*NaCl|%.*w\/v|w\/v/i) { # PERCENT NaCl
				my @naclvals = split (/\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						$temp =~ s/%.*//; #removes the units
						$temp =~ s/w.*//i; #removes the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@naclmindata = median(@temp3);
						}
					else { 
						@naclmindata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif ($temp =~ / mm/i) { # convert MILLIMOLAR concentrations to % NaCl.  Assume mM Na+ means mM NaCl
					if ($temp =~ /-/) { # range of mm NaCl values
						$temp =~ s/\s.*mm//gi;
						$temp =~ s/\smm.*//gi;
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@naclmindata = median(@temp3);
						next;
						}
					else { # single % NaCl value
						$temp =~ s/\smm.*//gi;
						$temp = $temp + 0.0;
						@naclmindata = $temp; 
						}
				@naclmindata = "@naclmindata" * 5.8443; # converts M to g/100mL
				}
			elsif ($temp =~ / m/i) { # convert MOLAR concentrations to % NaCl.  Assume M Na+ means M NaCl
				if ($temp =~ /MgCl2|MgSO4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmindata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@naclmindata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								$temp =~ s/m.*//i; #removes the units
								my @medianvals = (split /-/, $temp);
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0; #ensures string is a number
								@naclmindata = "@val3" * 5.8443; # converts M to g/100mL
								}
							else {
								map {s/m\s.*//g; } @val1;
								my @val2 = "@val1" + 0.0; #ensures number is a string
								@naclmindata = "@val2" * 5.8443; # converts M to g/100mL
								}
							}
						}
					}
				}
			elsif ($temp =~ /\sg|g.*l|g.*salt.*l|g.*NaCl/i) { # convert G PER L to % NaCl
				my @naclvals = (split /\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/gl.*//gi; } @val1;
				map {s/g\sl.*//gi; } @val1;
				map {s/g\ssalt\sl.*//gi; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of values
						$temp =~ s/g.*//i; # remove the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@naclmindata = "@val3" * 0.10;
						}
					else { # single value
						my @val2 = "@val1" + 0.0;
						@naclmindata = "@val2" * 0.10;
						}
					}
				}
			elsif ($temp =~ /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if ($temp =~ /-/) { # RANGE of percent values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@naclmindata = median(@temp3);
					}
				else { 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmindata = "@val1" + 0.0; # single % NaCl value
					}
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @naclmindata, "\n"; # prints to $nexusoutfile naclmin.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#process nacl optimum character
###################
elsif ($character eq "nacl optimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "naclopt.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @nacloptdata = split (/\t/, $line);
		push (@taxlabels, $nacloptdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'NaCl Optimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'NaCl opt'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @naclopttable = split (/\t/, $line);
		push (my @taxnames, $naclopttable[0]);
		push (my @nacloptdata, $naclopttable[1]);
		map {s/,\d\smm,/,/g; } @nacloptdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\smmw\/v,/,/g; } @nacloptdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d.\d\smm,/,/g; } @nacloptdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d.\d\d\smm,/,/g; } @nacloptdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\sdays,/,/gi; } @nacloptdata; # not a salt min
		map {s/plus.*//gi; } @nacloptdata; # get rid of any additional values
		map {s/^,//g; } @nacloptdata; # get rid of first comma
		map {s/,$//g; } @nacloptdata; # get rid of last comma
		map {s/,.*//g; } @nacloptdata; # get rid of any additional values
		local $, = " ";
		my $temp;
		
		#convert all concentration values to % NaCl and get rid of units
		foreach $temp (@nacloptdata) {
			if ($temp =~ /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@nacloptdata = "@val1" + 0.0;
				}	
			elsif ($temp =~ /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@nacloptdata = 0.05;
				}	
			elsif ($temp =~ /not required|absence of NaCl|without addition of NaCl/i ) {
				@nacloptdata = 0.0;
				}	
			elsif ($temp =~ /sea.*salt/i) { # converts PERCENT SEA SALTS to % NaCl
				if ($temp =~ /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @nacloptdata; # get rid of the first comma
					map {s/\,.*//; } @nacloptdata; # get rid of any additional values
					map {s/%%.*//; } @nacloptdata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@nacloptdata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif ($temp =~ /-/) { # % sea salt values with a range
					map {s/\,//; } @nacloptdata; # get rid of the first comma
					map {s/\,.*//; } @nacloptdata; # get rid of any additional values
					map {s/%.*//; } @nacloptdata; # get rid of any additional values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/%.*//; } @mval2; # get rid of any additional values
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					my @mval3 = median(@temp3);
					my @mval4 = "@mval3" + 0.0;
					@nacloptdata = "@mval4" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				else { # single value % sea salt value
					map {s/\,//; } @nacloptdata; # get rid of the first comma
					map {s/\,.*//; } @nacloptdata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@nacloptdata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif ($temp =~ /sea.*water/i) { #convert SEAWATER DILUTIONS to % Nacl
					if ($temp =~ /x|%/) { 
						if ($temp =~ /%/) { # when the seawater concentration reported as a per cent dilution
							my @naclvals = (split /%/, $temp);
							push (my @val1, $naclvals[0]);
							my @val2 = "@val1" + 0.0 ; # makes sure the string is a number
							my @val3 = @val2 / 100.0 ; # converts percent to a decimal
							@nacloptdata = "@val3" * 3.5 ; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /x/i) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = (split /\s/, $temp);
							push (my @val1, $naclvals[0]);
							map {s/x.*//i; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@nacloptdata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = (split /\//, $temp);
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0; #ensure that the numerator in the string is a nmber
							map {s/x.*//i; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0; #ensure that the denominator in the string is a number
							my @temp3 = @temp1 / @temp2;
							@nacloptdata = "@temp3" * 3.5 ;
							}
						}
					else { # has elsif (@nacloptdata = /^seawater$/) # assume it is pure seawater, undiluted seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@nacloptdata = "@val1" + 0.0;
						}
				}
			elsif ($temp =~ /%.*NaCl|%.*w\/v|w\/v/i) { # PERCENT NaCl
				my @naclvals = split (/\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						$temp =~ s/%.*//; #removes the units
						$temp =~ s/w.*//i; #removes the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@nacloptdata = median(@temp3);
						}
					else { 
						@nacloptdata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif ($temp =~ /mm/i) { # convert MILLIMOLAR concentrations to % NaCl.  Assume mM Na+ means mM NaCl
					if ($temp =~ /-/) { # range of mm NaCl values
						$temp =~ s/\s.*mm//gi;
						$temp =~ s/\smm.*//gi;
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@nacloptdata = median(@temp3);
						next;
						}
					else { # single % NaCl value
						$temp =~ s/\smm.*//gi;
						$temp = $temp + 0.0;
						@nacloptdata = $temp; 
						}
				@nacloptdata = "@nacloptdata" * 5.8443; # converts M to g/100mL
				}
			elsif ($temp =~ / m/i) { # convert MOLAR concentrations to % NaCl.  Assume M Na+ means M NaCl
				if ($temp =~ /MgCl2|MgSO4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@nacloptdata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@nacloptdata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								$temp =~ s/m.*//i; #removes the units
								my @medianvals = (split /-/, $temp);
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0; #ensures string is a number
								@nacloptdata = "@val3" * 5.8443; # converts M to g/100mL
								}
							else {
								map {s/m\s.*//gi; } @val1;
								my @val2 = "@val1" + 0.0; #ensures number is a string
								@nacloptdata = "@val2" * 5.8443; # converts M to g/100mL
								}
							}
						}
					}
				}
			elsif ($temp =~ /\sg|g.*l|g.*salt.*l|g.*NaCl/i) { # convert G PER L to % NaCl
				my @naclvals = (split /\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/gl.*//gi; } @val1;
				map {s/g\sl.*//gi; } @val1;
				map {s/g\ssalt\sl.*//gi; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of values
						$temp =~ s/g.*//i; # remove the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@nacloptdata = "@val3" * 0.10;
						}
					else { # single value
						my @val2 = "@val1" + 0.0;
						@nacloptdata = "@val2" * 0.10;
						}
					}
				}
			elsif ($temp =~ /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if ($temp =~ /-/) { # RANGE of percent values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@nacloptdata = median(@temp3);
					}
				else { 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@nacloptdata = "@val1" + 0.0; # single % NaCl value
					}
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @nacloptdata, "\n"; # prints to $nexusoutfile naclopt.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#process nacl maximum character
###################
elsif ($character eq "nacl maximum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "naclmax.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @naclmaxdata = split (/\t/, $line);
		push (@taxlabels, $naclmaxdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'NaCl Maximum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'NaCl max'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @naclmaxtable = split (/\t/, $line);
		push (my @taxnames, $naclmaxtable[0]);
		push (my @naclmaxdata, $naclmaxtable[1]);
		map {s/,\d\smm,/,/g; } @naclmaxdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\smmw\/v,/,/g; } @naclmaxdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\.\d\smmw\/v,/,/g; } @naclmaxdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\.\d\smm,/,/g; } @naclmaxdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d.\d\d\smm,/,/g; } @naclmaxdata; # in the current version of microPIE colony sizes are incorrectly being picked up as salinity characters
		map {s/,\d\sdays,/,/gi; } @naclmaxdata; # not a salt min
		map {s/plus.*//gi; } @naclmaxdata; # get rid of any additional values
		map {s/^,//g; } @naclmaxdata; # get rid of first comma
		map {s/,$//g; } @naclmaxdata; # get rid of last comma
		map {s/,.*//g; } @naclmaxdata; # get rid of any additional values
		local $, = " ";
		my $temp;
		
		#convert all concentration values to % NaCl and get rid of units
		foreach $temp (@naclmaxdata) {
			if ($temp =~ /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@naclmaxdata = "@val1" + 0.0;
				}	
			elsif ($temp =~ /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@naclmaxdata = 0.05;
				}	
			elsif ($temp =~ /not required|absence of NaCl|without addition of NaCl/i ) {
				@naclmaxdata = 0.0;
				}	
			elsif ($temp =~ /sea.*salt/i) { # converts PERCENT SEA SALTS to % NaCl
				if ($temp =~ /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @naclmaxdata; # get rid of the first comma
					map {s/\,.*//; } @naclmaxdata; # get rid of any additional values
					map {s/%%.*//; } @naclmaxdata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmaxdata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif ($temp =~ /-/) { # % sea salt values with a range
					map {s/\,//; } @naclmaxdata; # get rid of the first comma
					map {s/\,.*//; } @naclmaxdata; # get rid of any additional values
					map {s/%.*//; } @naclmaxdata; # get rid of any additional values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/%.*//; } @mval2; # get rid of any additional values
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					my @mval3 = median(@temp3);
					my @mval4 = "@mval3" + 0.0;
					@naclmaxdata = "@mval4" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				else { # single value % sea salt value
					map {s/\,//; } @naclmaxdata; # get rid of the first comma
					map {s/\,.*//; } @naclmaxdata; # get rid of any additional values
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmaxdata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif ($temp =~ /sea.*water/i) { #convert SEAWATER DILUTIONS to % Nacl
					if ($temp =~ /x|%/) { 
						if ($temp =~ /%/) { # when the seawater concentration reported as a per cent dilution
							my @naclvals = (split /%/, $temp);
							push (my @val1, $naclvals[0]);
							my @val2 = "@val1" + 0.0 ; # makes sure the string is a number
							my @val3 = @val2 / 100.0 ; # converts percent to a decimal
							@naclmaxdata = "@val3" * 3.5 ; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /x/i) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = (split /\s/, $temp);
							push (my @val1, $naclvals[0]);
							map {s/x.*//i; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@naclmaxdata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						elsif ($temp =~ /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = (split /\//, $temp);
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0; #ensure that the numerator in the string is a nmber
							map {s/x.*//i; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0; #ensure that the denominator in the string is a number
							my @temp3 = @temp1 / @temp2;
							@naclmaxdata = "@temp3" * 3.5 ;
							}
						}
					else { # has elsif (@naclmaxdata = /^seawater$/) # assume it is pure seawater, undiluted seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@naclmaxdata = "@val1" + 0.0;
						}
				}
			elsif ($temp =~ /%.*NaCl|%.*w\/v|w\/v/i) { # PERCENT NaCl
				my @naclvals = split (/\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						$temp =~ s/%.*//; #removes the units
						$temp =~ s/w.*//i; #removes the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@naclmaxdata = median(@temp3);
						}
					else { 
						@naclmaxdata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif ($temp =~ /mm/i) { # convert MILLIMOLAR concentrations to % NaCl.  Assume mM Na+ means mM NaCl
					if ($temp =~ /-/) { # range of mm NaCl values
						$temp =~ s/\s.*mm//gi;
						$temp =~ s/\smm.*//gi;
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0; #ensure the first string in a range is a number
						my @temp2 = "@mval2" + 0.0; #ensure the second string in a range is a number
						my @temp3 = (@temp1, @temp2);
						@naclmaxdata = median(@temp3);
						next;
						}
					else { # single % NaCl value
						$temp =~ s/\smm.*//gi;
						$temp = $temp + 0.0;
						@naclmaxdata = $temp; 
						}
				@naclmaxdata = "@naclmaxdata" * 5.8443; # converts M to g/100mL
				}
			elsif ($temp =~ / m/i) { # convert MOLAR concentrations to % NaCl.  Assume M Na+ means M NaCl
				if ($temp =~ /MgCl2|MgSO4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmaxdata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@naclmaxdata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								$temp =~ s/m.*//i; #removes the units
								my @medianvals = (split /-/, $temp);
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0; #ensures string is a number
								@naclmaxdata = "@val3" * 5.8443; # converts M to g/100mL
								}
							else {
								map {s/m\s.*//gi; } @val1;
								my @val2 = "@val1" + 0.0; #ensures number is a string
								@naclmaxdata = "@val2" * 5.8443; # converts M to g/100mL
								}
							}
						}
					}
				}
			elsif ($temp =~ /\sg|g.*l|g.*salt.*l|g.*NaCl/i) { # convert G PER L to % NaCl
				my @naclvals = (split /\s/, $temp);
				push (my @val1, $naclvals[0]);
				map {s/gl.*//gi; } @val1;
				map {s/g\sl.*//gi; } @val1;
				map {s/g\ssalt\sl.*//gi; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of values
						$temp =~ s/g.*//i; # remove the units
						my @medianvals = (split /-/, $temp);
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@naclmaxdata = "@val3" * 0.10;
						}
					else { # single value
						my @val2 = "@val1" + 0.0;
						@naclmaxdata = "@val2" * 0.10;
						}
					}
				}
			elsif ($temp =~ /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if ($temp =~ /-/) { # RANGE of percent values
					my @medianvals = (split /-/, $temp);
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@naclmaxdata = median(@temp3);
					}
				else { 
					my @naclvals = (split /\s/, $temp);
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmaxdata = "@val1" + 0.0; # single % NaCl value
					}
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @naclmaxdata, "\n"; # prints to $nexusoutfile naclmax.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}



###################
#process pH minimum character
###################

elsif ($character eq "ph minimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "pHmin.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @pHmindata = split (/\t/, $line);
		push (@taxlabels, $pHmindata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'pH Minimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'pH minimum'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @pHmintable = split (/\t/, $line);
		push (my @taxnames, $pHmintable[0]);
		push (my @pHmindata, $pHmintable[1]);
		map {s/,not.*,/,/g; } @pHmindata; # delete any extraneous text
		map {s/\sNaCl//gi; } @pHmindata; # delete any extraneous text
		map {s/,.*days,/,/; } @pHmindata; #delete any extraneous text
		map {s/\sm,/,/; } @pHmindata; #delete any extraneous text
		map {s/^,//; } @pHmindata; #remove first comma
		map {s/,$//; } @pHmindata; #remove last comma
		map {s/around neutral/7.0/g; } @pHmindata; # synonyms of 
		foreach my $temp (@pHmindata) {
			if ($temp =~ /,/) { # if there is still two values separated by a comma
				my @testvals = split (/,/, $pHmindata[0]);
				my @firstval = splice @testvals, 0, 1;
				my @secondval = splice @testvals, 0, 1;
				my @firsttestval = @firstval;
				my @secondtestval = @secondval;
				map {s/-.*//g; } @firsttestval; # synonyms of 
				map {s/-.*//g; } @secondtestval; # synonyms of 
				if (grep {$_ > 12} @firsttestval) {#tests whether there is a credible pH value
					@pHmindata = @secondval;
					}
				elsif (grep {$_ < 12} @firsttestval) {#tests whether there is a credible pH value
					@pHmindata = @firstval;
					}
				elsif (grep {$_ > 12} @secondtestval) {#tests whether there is a credible pH value
					@pHmindata = @firstval;
					}
				elsif (grep {$_ < 12} @secondtestval) {#tests whether there is a credible pH value
					@pHmindata = @secondval;
					}
				else {
					next;
					}
				}
			}
		my $temp;
		foreach $temp (@pHmindata) {
			if ($temp =~ /-/) { #there is a range
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0; #ensures the first string is a number
				my @temp2 = "@val2" + 0.0; #ensures the first string is a number
				my @temp3 = (@temp1, @temp2); #takes the median of the two numbers
				@pHmindata = median(@temp3);
			}
			else { #there is a single value
				next;
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHmindata, "\n"; # prints to $nexusoutfile pHmin.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	}

###################
#process pH optimum character
###################
elsif ($character eq "ph optimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "pHopt.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @pHoptdata = split (/\t/, $line);
		push (@taxlabels, $pHoptdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'pH Optimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'pH optimum'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @pHopttable = split (/\t/, $line);
		push (my @taxnames, $pHopttable[0]);
		push (my @pHoptdata, $pHopttable[1]);
		map {s/,not.*,/,/g; } @pHoptdata; # delete any extraneous text
		map {s/\sNaCl//gi; } @pHoptdata; # delete any extraneous text
		map {s/,.*days,/,/; } @pHoptdata; #delete any extraneous text
		map {s/\sm,/,/; } @pHoptdata; #delete any extraneous text
		map {s/^,//; } @pHoptdata; #remove first comma
		map {s/,$//; } @pHoptdata; #remove last comma
		map {s/around neutral/7.0/g; } @pHoptdata; # synonyms of 
		foreach my $temp (@pHoptdata) {
			if ($temp =~ /,/) { # if there is still two values separated by a comma
				my @testvals = split (/,/, $pHoptdata[0]);
				my @firstval = splice @testvals, 0, 1;
				my @secondval = splice @testvals, 0, 1;
				my @firsttestval = @firstval;
				my @secondtestval = @secondval;
				map {s/-.*//g; } @firsttestval; # synonyms of 
				map {s/-.*//g; } @secondtestval; # synonyms of 
				if (grep {$_ > 12} @firsttestval) {#tests whether there is a credible pH value
					@pHoptdata = @secondval;
					}
				elsif (grep {$_ < 12} @firsttestval) {#tests whether there is a credible pH value
					@pHoptdata = @firstval;
					}
				elsif (grep {$_ > 12} @secondtestval) {#tests whether there is a credible pH value
					@pHoptdata = @firstval;
					}
				elsif (grep {$_ < 12} @secondtestval) {#tests whether there is a credible pH value
					@pHoptdata = @secondval;
					}
				else {
					next;
					}
				}
			}
		my $temp;
		foreach $temp (@pHoptdata) {
			if ($temp =~ /-/) { #there is a range
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0; #ensures the first string is a number
				my @temp2 = "@val2" + 0.0; #ensures the first string is a number
				my @temp3 = (@temp1, @temp2); #takes the median of the two numbers
				@pHoptdata = median(@temp3);
			}
			else { #there is a single value
				next;
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHoptdata, "\n"; # prints to $nexusoutfile pHopt.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}
###################
#process pH maximum character
###################
elsif ($character eq "ph maximum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "pHmax.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @pHmaxdata = split (/\t/, $line);
		push (@taxlabels, $pHmaxdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'pH Maximum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'pH maximum'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @pHmaxtable = split (/\t/, $line);
		push (my @taxnames, $pHmaxtable[0]);
		push (my @pHmaxdata, $pHmaxtable[1]);
		map {s/,not.*,/,/g; } @pHmaxdata; # delete any extraneous text
		map {s/\sNaCl//gi; } @pHmaxdata; # delete any extraneous text
		map {s/,.*days,/,/; } @pHmaxdata; #delete any extraneous text
		map {s/\sm,/,/; } @pHmaxdata; #delete any extraneous text
		map {s/^,//; } @pHmaxdata; #remove first comma
		map {s/,$//; } @pHmaxdata; #remove last comma
		map {s/around neutral/7.0/g; } @pHmaxdata; # synonyms of 
		foreach my $temp (@pHmaxdata) {
			if ($temp =~ /,/) { # if there is still two values separated by a comma
				my @testvals = split (/,/, $pHmaxdata[0]);
				my @firstval = splice @testvals, 0, 1;
				my @secondval = splice @testvals, 0, 1;
				my @firsttestval = @firstval;
				my @secondtestval = @secondval;
				map {s/-.*//g; } @firsttestval; # synonyms of 
				map {s/-.*//g; } @secondtestval; # synonyms of 
				if (grep {$_ > 12} @firsttestval) {#tests whether there is a credible pH value
					@pHmaxdata = @secondval;
					}
				elsif (grep {$_ < 12} @firsttestval) {#tests whether there is a credible pH value
					@pHmaxdata = @firstval;
					}
				elsif (grep {$_ > 12} @secondtestval) {#tests whether there is a credible pH value
					@pHmaxdata = @firstval;
					}
				elsif (grep {$_ < 12} @secondtestval) {#tests whether there is a credible pH value
					@pHmaxdata = @secondval;
					}
				else {
					next;
					}
				}
			}
		my $temp;
		foreach $temp (@pHmaxdata) {
			if ($temp =~ /-/) { #there is a range
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0; #ensures the first string is a number
				my @temp2 = "@val2" + 0.0; #ensures the first string is a number
				my @temp3 = (@temp1, @temp2); #takes the median of the two numbers
				@pHmaxdata = median(@temp3);
			}
			else { #there is a single value
				next;
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHmaxdata, "\n"; # prints to $nexusoutfile pHmax.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}
###################
#process temperature minimum character
###################
elsif ($character eq "temperature minimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "tempmin.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @tempmindata = split (/\t/, $line);
		push (@taxlabels, $tempmindata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Temperature Minimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'temperature min'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @tempmintable = split (/\t/, $line);
		push (my @taxnames, $tempmintable[0]);
		push (my @tempmindata, $tempmintable[1]);

		map {s/,not.*,/,/g; } @tempmindata; 
		map {s/^,//; } @tempmindata; #remove first comma
		map {s/,$//; } @tempmindata; # remove second comma
		map {s/,\S.*/,/g; } @tempmindata; # removes any second values 
		map {s/,\d.*/,/g; } @tempmindata; # removes any second values
		map {s/_c//gi; } @tempmindata; # removes any underscore celsius units
		map {s/\snacl//gi; } @tempmindata; # removes any underscore celsius units

		foreach my $temp (@tempmindata) {
			$temp =~ s/[^[:print:]]//g; #remove unprintable characters
			$temp =~ s/\sC//gi; #remove the celsius label
			$temp =~ s/,//g;
			if ($temp=~ /\d-\d/) { #this is a range of values
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@tempmindata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @tempmindata, "\n"; # prints to $nexusoutfile tempmin.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}
###################
#process temperature optimum character
###################
elsif ($character eq "temperature optimum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "tempopt.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @tempoptdata = split (/\t/, $line);
		push (@taxlabels, $tempoptdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Temperature Optimum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'temperature opt'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @tempopttable = split (/\t/, $line);
		push (my @taxnames, $tempopttable[0]);
		push (my @tempoptdata, $tempopttable[1]);

		map {s/,not.*,/,/g; } @tempoptdata; 
		map {s/^,//; } @tempoptdata; #remove first comma
		map {s/,$//; } @tempoptdata; # remove second comma
		map {s/,\S.*/,/g; } @tempoptdata; # removes any second values 
		map {s/,\d.*/,/g; } @tempoptdata; # removes any second values
		map {s/_c//gi; } @tempoptdata; # removes any underscore celsius units
		map {s/\snacl//gi; } @tempoptdata; # removes any underscore celsius units

		foreach my $temp (@tempoptdata) {
			$temp =~ s/[^[:print:]]//g; #remove unprintable characters
			$temp =~ s/\sC//gi; #remove the celsius label
			$temp =~ s/,//g;
			if ($temp=~ /\d-\d/) { #this is a range of values
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@tempoptdata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @tempoptdata, "\n"; # prints to $nexusoutfile tempopt.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}

###################
#process temperature maximum character
###################
elsif ($character eq "temperature maximum") {
#prepare nexus file
	my @taxlabels;
	my $nexusoutfile = "tempmax.nex";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @tempmaxdata = split (/\t/, $line);
		push (@taxlabels, $tempmaxdata[0]);
		}
	my $numtax = scalar(@taxlabels) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxlabels;
	local $, = " ";
	print OUT "\t\t", @taxlabels, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Temperature Maximum Matrix'\;\n\tDIMENSIONS NCHAR=1\;\n\tFORMAT DATATYPE \= CONTINUOUS INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'temperature max'\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @tempmaxtable = split (/\t/, $line);
		push (my @taxnames, $tempmaxtable[0]);
		push (my @tempmaxdata, $tempmaxtable[1]);

		map {s/,not.*,/,/g; } @tempmaxdata; 
		map {s/^,//; } @tempmaxdata; #remove first comma
		map {s/,$//; } @tempmaxdata; # remove second comma
		map {s/,\S.*/,/g; } @tempmaxdata; # removes any second values 
		map {s/,\d.*/,/g; } @tempmaxdata; # removes any second values
		map {s/_c//gi; } @tempmaxdata; # removes any underscore celsius units
		map {s/\snacl//gi; } @tempmaxdata; # removes any underscore celsius units

		foreach my $temp (@tempmaxdata) {
			$temp =~ s/[^[:print:]]//g; #remove unprintable characters
			$temp =~ s/\sC//gi; #remove the celsius label
			$temp =~ s/,//g;
			if ($temp=~ /\d-\d/) { #this is a range of values
				my @medianvals = (split /-/, $temp);
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@tempmaxdata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @tempmaxdata, "\n"; # prints to $nexusoutfile tempmax.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
}

###################
#process pigment compounds character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "pigment compounds") {
	my $homout = "hom.pigments.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @pigments, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @pigments;
		map {s/,a non-diffusible,/,not diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,brown diffusible pigment,/,brown,diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,bright yellow non-flexirubin-type pigments,/,yellow,not flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,carotenoid pigments,/,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,carotenoid zeaxanthin,/,carotenoids,zeaxanthin,pigment,/g; } @pigments; # synonyms of 
		map {s/,carotenoid-type,/,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,carotenoid-type pigments,/,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,diffusible pigment,/,diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,diffusible yellow pigment,/,diffusible,yellow,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin pigment,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin pigments,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin type,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin-like pigment,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin-type pigment,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubin-type pigments,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubins,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,intracellular carotenoid pigments,/,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,mainly zeaxanthin,/,zeaxanthin,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible and non-fluorescent yellow pigment,/,not diffusible,not fluorescent,yellow,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible carotenoid pigments,/,not diffusible,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible orange,/,not diffusible,orange,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible orange pigment,/,not diffusible,orange,pigment,/g; } @pigments; # synonyms of 
		map {s/,flexirubrin,/,flexirubin,/g; } @pigments; # synonyms of 

		map {s/,non-diffusible pigment,/,not diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible pigments,/,not diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible pink pigments,/,not diffusible,pink,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible yellow pigments,/,not diffusible,yellow,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible yellow carotenoid pigments,/,not diffusible,yellow,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,non-diffusible yellow pigments,/,not diffusible,yellow,pigment,/g; } @pigments; # synonyms of 
		map {s/,nondiffusible yellow carotenoid pigments,/,not diffusible,yellow,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,orange carotenoid-like pigments,/,orange,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,orange-red carotenoid,/,orange-red,carotenoids,pigment,/g; } @pigments; # synonyms of 
		map {s/,pigment compounds,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,pigmented,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,pigmented yellow-orange,/,yellow-orange,pigment,/g; } @pigments; # synonyms of 
		map {s/,reddish-brown diffusible pigment,/,reddish-brown,diffusible,pigment,/g; } @pigments; # synonyms of 
		map {s/,the flexirubin-type,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,typical flexirubin reaction,/,flexirubin,pigment,/g; } @pigments; # synonyms of 
		map {s/,yellow carotenoid pigments,/,yellow,carotenoids,pigment,/g; } @pigments; # synonyms of 

		map {s/,not flexirubrin,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not typical flexirubin reaction,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not the flexirubin-type,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not brown diffusible pigment,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,absent flexirubin pigments,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,absent flexirubin-type pigment,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,absent flexirubin-type pigments,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,flexirubin-negative,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not diffusible yellow pigment,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,not non-diffusible yellow,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,not non-diffusible orange pigment,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,not non-diffusible and non-fluorescent yellow pigment,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,not flexirubins,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not diffusible pigment,/,not pigment,/g; } @pigments; # synonyms of 
		map {s/,not flexirubin-type pigment,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not flexirubin-type pigments,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not flexirubin type,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not flexirubin pigment,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,not flexirubin pigments,/,not flexirubin,/g; } @pigments; # synonyms of 
		map {s/,absorption maxim.*,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,absorption spectrum.*,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,maximum absorption.*,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,452 and 478 nm,/,pigment,/g; } @pigments; # synonyms of 
		map {s/,470 and 400 nm,/,pigment,/g; } @pigments; # synonyms of 
#not pigment compounds
		map {s/,448,/,/g; } @pigments; # synonyms of 

		print OUT @taxlabels, @pigments, "\n"; # prints to $homout, hom.pigments.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.pigments.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,pigments,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.pigments.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.pigments.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.pigments.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.pigments.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.pigments.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^brown/**brown/g; } @charlist2;
	map {s/^carotenoids/**carotenoids/g; } @charlist2;
	map {s/^diffusible/**diffusible/g; } @charlist2;
	map {s/^flexirubin/**flexirubin/g; } @charlist2;
	map {s/^fluorescent/**fluorescent/g; } @charlist2;
	map {s/^orange-red/**orange-red/g; } @charlist2;
	map {s/^orange/**orange/g; } @charlist2;
	map {s/^pigment/**pigment/g; } @charlist2;
	map {s/^pink/**pink/g; } @charlist2;
	map {s/^red/**red/g; } @charlist2;
	map {s/^reddish-brown/**reddish-brown/g; } @charlist2;
	map {s/^yellow-orange/**yellow-orange/g; } @charlist2;
	map {s/^yellow/**yellow/g; } @charlist2;
	map {s/^zeaxanthin/**zeaxanthin/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.pigments.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.pigments.txt
		}
	my $temp6 = "temp6.pigments.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.pigments.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "pigments.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @pigmentsdata = split (/\t/, $line);
		push (@taxnames, $pigmentsdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Pigment Compounds Matrix'\;\n\tDIMENSIONS NCHAR=7\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 pigments \/  unpigmented pigmented, ";
	print OUT "2 'pigment diffusibility' \/  'pigments not diffusible' 'pigments diffusible', ";
	print OUT "3 'cell fluoresence' \/  'cells not fluorescent' 'cells fluorescent', ";
	print OUT "4 flexirubin \/  'flexirubin absent' 'flexirubin present', ";
	print OUT "5 carotenoids \/  'carotenoids absent' 'carotenoids present', ";
	print OUT "6 zeaxanthin \/  'zeaxanthin absent' 'zeaxanthin present', ";
	print OUT "7 'culture color' \/  'culture brown' 'culture reddish-brown' 'culture red' 'culture pink' 'culture orange-red' 'culture orange' 'culture yellow-orange' 'culture yellow', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @pigmentsdata = split (/\t/, $line);
		push (my @taxnames, $pigmentsdata[0]);

#code char 1 pigment
		if ($line =~ /,pigment,|,flexirubin,|,carotenoids,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not pigment,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 diffusibility
		if ($line =~ /,diffusible,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diffusible,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 fluorescence
		if ($line =~ /,fluorescent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fluorescent,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 flexirubin
		if ($line =~ /,flexirubin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flexirubin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 carotenoids
		if ($line =~ /,carotenoids,|,zeaxanthin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carotenoids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 zeaxanthin
		if ($line =~ /,zeaxanthin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not zeaxanthin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 colors
		if ($line =~ /,brown,/) {
			print OUT "0";
			}
		elsif ($line =~ /,reddish-brown,/) {
			print OUT "1";
			}
		elsif ($line =~ /,red,/) {
			print OUT "2";
			}
		elsif ($line =~ /,pink,/) {
			print OUT "3";
			}
		elsif ($line =~ /,orange-red,/) {
			print OUT "4";
			}
		elsif ($line =~ /,orange,/) {
			print OUT "5";
			}
		elsif ($line =~ /,yellow-orange,/) {
			print OUT "6";
			}
		elsif ($line =~ /,yellow,/) {
			print OUT "7";
			}
		else {
			print OUT "?";
		}

	print OUT "\n";

	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}
###################
#process aerophilicity character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "aerophilicity") {
	my $homout = "hom.aerophilicity.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @aerophilicity, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @aerophilicity;
		map {s/,aerobic,strict,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic,strict,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobe,strict,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,strict,aerobes,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,strict,anaerobes,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobe,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,microaerobically,/,microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @aerophilicity; # synonyms of 

		map {s/,aerobe,/,aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,aerobes,/,aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,aerobic growth,/,aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,aerobic metabolism,/,aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,aerobic strictly,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,an obligately anaerobic,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic 4 days incubation,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic conditions,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic facultatively organism,/,facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic growth,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic incubation,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic respiration,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic seal,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobic strictly conditions,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,atmospheric oxygen,/,aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,facultative aerobes,/,facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,facultative anaerobes,/,facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,facultatively aerobic,/,facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,facultatively anaerobic,/,facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,microaerophilic conditions,/,microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,microaerophilic microaerophilic,/,microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,obligate aerobe,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,obligate aerobes,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,obligately aerobic,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,obligately anaerobic,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,obligatory anaerobic,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,strictly aerobic,/,strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,strictly anaerobic,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 

		map {s/,not aerobe,/,not aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not aerobes,/,not aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not aerobic growth,/,not aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not aerobic metabolism,/,not aerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not aerobic strictly,/,not strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not an obligately anaerobic,/,not strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic 4 days incubation,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic conditions,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic facultatively organism,/,not facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic growth,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic incubation,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic respiration,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic seal,/,not anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,not anaerobic strictly conditions,/,not strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not atmospheric oxygen,/,not aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not facultative aerobes,/,not facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not facultative anaerobes,/,not facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not facultatively aerobic,/,not facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not facultatively anaerobic,/,not facultative anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not microaerophilic conditions,/,not microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,not microaerophilic microaerophilic,/,not microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,not obligate aerobe,/,not strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not obligate aerobes,/,not strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not obligately aerobic,/,not strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not obligately anaerobic,/,not strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not obligatory anaerobic,/,not strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not strictly aerobic,/,not strict aerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,not strictly anaerobic,/,not strict anaerobe,/g; } @aerophilicity; # synonyms of 

		map {s/,microaerobe,/,microaerophilic,/g; } @aerophilicity; # synonyms of 
		map {s/,anaerobes,/,anaerobic,/g; } @aerophilicity; # synonyms of 
		map {s/,fastidious anaerobe,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 
		map {s/,obligately anaerobes,/,strict anaerobe,/g; } @aerophilicity; # synonyms of 

		map {s/,strictly,/,/g; } @aerophilicity; # synonyms of 
		map {s/,strict,/,/g; } @aerophilicity; # synonyms of 
		print OUT @taxlabels, @aerophilicity, "\n"; # prints to $homout, hom.aerophilicity.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.aerophilicity.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,aerophilicity,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.aerophilicity.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.aerophilicity.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.aerophilicity.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.aerophilicity.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.aerophilicity.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^aerobic/**aerobic/g; } @charlist2;
	map {s/^aerotolerant/**aerotolerant/g; } @charlist2;
	map {s/^anaerobic/**anaerobic/g; } @charlist2;
	map {s/^facultative anaerobe/**facultative anaerobe/g; } @charlist2;
	map {s/^microaerophilic/**microaerophilic/g; } @charlist2;
	map {s/^strict aerobe/**strict aerobe/g; } @charlist2;
	map {s/^strict anaerobe/**strict anaerobe/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.aerophilicity.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.aerophilicity.txt
		}
	my $temp6 = "temp6.aerophilicity.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.aerophilicity.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "aerophilicity.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @aerophilicitydata = split (/\t/, $line);
		push (@taxnames, $aerophilicitydata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Aerophilicity Matrix'\;\n\tDIMENSIONS NCHAR=8\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 aerophilicity \/  aerobic microaerophilic anaerobic, ";
	print OUT "2 'strict aerobe' \/  'not strict aerobe' 'strict aerobe', ";
	print OUT "3 aerobic \/  'not aerobic' aerobic, ";
	print OUT "4 microaerophilic \/  'not microaerophilic' microaerophilic, ";
	print OUT "5 'facultative anaerobe' \/  'not facultative anaerobe' 'facultative anaerobe', ";
	print OUT "6 anaerobic \/  'not anaerobic' anaerobic, ";
	print OUT "7 'strict anaerobe' \/  'not strict anaerobe' 'strict anaerobe', ";
	print OUT "8 aerotolerant \/  'not aerotolerant' aerotolerant, ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @aerophilicitydata = split (/\t/, $line);
		push (my @taxnames, $aerophilicitydata[0]);

#code char 1 aerophilicity
		if ($line =~ /,facultative anaerobe,|,aerobic.*anaerobic,|,anaerobic.*aerobic,/) {
			print OUT @taxnames, "(0 2)";
			}
		elsif ($line =~ /,aerobic,|,strict aerobe,|,not anaerobic,/) {
			print OUT @taxnames, "0";
			}
		elsif ($line =~ /,microaerophilic,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,anaerobic,|,strict anaerobe,/) {
			print OUT @taxnames, "2";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 strict aerobe
		if ($line =~ /,strict aerobe,|,not anaerobic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not strict aerobe,|,aerobic.*anaerobic,|,anaerobic.*aerobic,|,strict anaerobe,|,facultative anaerobe,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 aerobic
		if ($line =~ /,aerobic,|,strict aerobe,|,microaerophilic,|,not anaerobic,|,facultative anaerobe,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aerobic,|,strict anaerobe,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 microaerophilic
		if ($line =~ /,microaerophilic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not microaerophilic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 facultative anaerobe
		if ($line =~ /,facultative anaerobe,|,aerobic.*anaerobic,|,anaerobic.*aerobic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not facultative anaerobe,|,strict aerobe,|,strict anaerobe,|,not anaerobic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 anaerobic
		if ($line =~ /,anaerobic,|,facultative anaerobe,|,strict anaerobe,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not anaerobic,|,strict aerobe,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 strict anaerobe
		if ($line =~ /,strict anaerobe,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not strict anaerobe,|,aerobic.*anaerobic,|,anaerobic.*aerobic,|,strict aerobe,|,not anaerobic,|,facultative anaerobe,|,aerobic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 aerotolerant
		if ($line =~ /,aerotolerant,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aerotolerant,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process vitamins and cofactors required for growth character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "vitamins and cofactors required for growth") {
	my $homout = "hom.vitcos.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @vitcos, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @vitcos;
		map {s/,[\d]+% chicken serum,/,serum,/gi; } @vitcos; # synonyms of 
		map {s/,.*alpha-toxin serum,/,serum,/gi; } @vitcos; # synonyms of 
		map {s/,.*serum bottles.*,/,serum,/gi; } @vitcos; # synonyms of 
		map {s/,1% ox-bile salts,/,ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,[\d]* bile,/,bile,/gi; } @vitcos; # synonyms of 
		map {s/,ox bile,/,ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,ox bile salts,/,ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,ox-bile salts,/,ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,20% bile,/,bile,/gi; } @vitcos; # synonyms of 
		map {s/,40% bile,/,bile,/gi; } @vitcos; # synonyms of 
		map {s/,bile salts,/,bile,/gi; } @vitcos; # synonyms of 
		map {s/,medium containing 20% bile,/,bile,/gi; } @vitcos; # synonyms of 
		map {s/,both vitamin b12,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,cobalamin,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,organic growth factors,/,growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,exogenous growth factor,/,growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,growth factor dependency,/,growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,growth factors,/,growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,hemin and coenzyme I,/,hemin,nad,/gi; } @vitcos; # synonyms of 
		map {s/,x.* + v factor,/,hemin,nad,/gi; } @vitcos; # synonyms of 
		map {s/,x+v factor disc.*,/,hemin,nad,/gi; } @vitcos; # synonyms of 
		map {s/,xv discs,/,hemin,nad,/gi; } @vitcos; # synonyms of 
		map {s/,xv disks,/,hemin,nad,/gi; } @vitcos; # synonyms of 
		map {s/,factor x,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,haem,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,haemin,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,x disc,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,x disk,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,x factor,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,x growth factor,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,x-factor,/,hemin,/gi; } @vitcos; # synonyms of 
		map {s/,factor v,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,v disc,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,v disk,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,v factor,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,v growth factor,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,v-factor,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,coenzyme I,/,nad,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin k,/,vitamin k1,/gi; } @vitcos; # synonyms of 
		map {s/,thiamine,/,thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b1,/,thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b-1,/,thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b-12,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,cyanocobalamin,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,cyanocobalamine,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,methylcobalamin,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,methylcobalamine,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,cobalamins,/,vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,.*pyridoxal hydrochloride,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin B1,/,thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b2,/,riboflavin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b3,/,niacin,/gi; } @vitcos; # synonyms of 
		map {s/,nicotinic acid,/,niacin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b5,/,pantothenic acid,/gi; } @vitcos; # synonyms of 
		map {s/,calcium pantothenate,/,pantothenic acid,/gi; } @vitcos; # synonyms of 
		map {s/,pyridoxine,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,pyridoxal,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,pyridoxal hydrochloride,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,.*% pyridoxal hydrochloride,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,pyridoxamine,/,vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b7,/,biotin,/gi; } @vitcos; # synonyms of 
		map {s/,vitamin b9,/,folic acid,/gi; } @vitcos; # synonyms of 
		map {s/,folate,/,folic acid,/gi; } @vitcos; # synonyms of 
		map {s/,water.*soluble vitamins,/,b vitamins,/gi; } @vitcos; # synonyms of 

		map {s/,no.*pyridoxal,/,not vitamin b6,/gi; } @vitcos; # synonyms of 
		map {s/,not.*thiamine,/,not thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*vitamin b1,/,not thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*vitamin b-1,/,not thiamin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*vitamin b-12,/,not vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,not.*vitamin k,/,not vitamin k1,/gi; } @vitcos; # synonyms of 
		map {s/,not.*2% chicken serum,/,not serum,/gi; } @vitcos; # synonyms of 
		map {s/,not.*1% ox-bile salts,/,not ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*ox bile,/,not ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*ox bile salts,/,not ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*ox-bile salts,/,not ox-bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*20% bile,/,not bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*40% bile,/,not bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*bile salts,/,not bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*medium containing 20% bile,/,not bile,/gi; } @vitcos; # synonyms of 
		map {s/,not.*both vitamin b12,/,not vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,not.*cobalamin,/,not vitamin b12,/gi; } @vitcos; # synonyms of 
		map {s/,not.*organic growth factors,/,not growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,not.*exogenous growth factor,/,not growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,not.*growth factor dependency,/,not growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,not.*growth factors,/,not growth factor,/gi; } @vitcos; # synonyms of 
		map {s/,not.*hemin and coenzyme I,/,not hemin,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x.* + v factor,/,not hemin,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x+v factor disc.*,/,not hemin,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*xv discs,/,not hemin,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*xv disks,/,not hemin,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*factor x,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*haem,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*haemin,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x disc,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x disk,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x factor,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x growth factor,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*x-factor,/,not hemin,/gi; } @vitcos; # synonyms of 
		map {s/,not.*factor v,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*v disc,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*v disk,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*v factor,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*v growth factor,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*v-factor,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,not.*coenzyme I,/,not nad,/gi; } @vitcos; # synonyms of 
		map {s/,neither pyridoxal hydrochloride nor satellitism,/,not vitamin b6,/gi; } @vitcos; # synonyms of 

#not vitamins or co-factors
		map {s/,clumping factor,/,/gi; } @vitcos; # synonyms of 
		map {s/,not.*clumping factor,/,/gi; } @vitcos; # synonyms of 

		print OUT @taxlabels, @vitcos, "\n"; # prints to $homout, hom.vitcos.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.vitcos.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,vitamins and cofactors required for growth,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.vitcos.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.vitcos.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.vitcos.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.vitcos.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.vitcos.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^bile/**bile/g; } @charlist2;
	map {s/^growth factor/**growth factor/g; } @charlist2;
	map {s/^hemin/**hemin/g; } @charlist2;
	map {s/^nad/**NAD/g; } @charlist2;
	map {s/^ox-bile/**ox-bile/g; } @charlist2;
	map {s/^serum/**serum/g; } @charlist2;
#vitamins
	map {s/^vitamins/**vitamins/g; } @charlist2;
	map {s/^vitamin k1/**vitamin k1/g; } @charlist2;
	map {s/^b vitamins/**b vitamins/g; } @charlist2;
	map {s/^thiamin/**thiamin/g; } @charlist2;#thiamin=vitamin B1
	map {s/^riboflavin/**riboflavin/g; } @charlist2; #riboflavin=vitamin B2
	map {s/^niacin/**niacin/g; } @charlist2; #niacin=vitamin B3
	map {s/^pantothenic acid/**pantothenic acid/g; } @charlist2; #pantothenic acid=vitamin B5
	map {s/^vitamin b6/**vitamin b6/g; } @charlist2; #pyridoxal hydrochloride=vitamin B6
	map {s/^biotin/**biotin/g; } @charlist2; #biotin=vitamin B7
	map {s/^folic acid/**folic acid/g; } @charlist2; #folic acid=vitamin B9
	map {s/^vitamin b12/**vitamin b12/g; } @charlist2; #cobalamin=vitamin B12

	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.vitcos.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.vitcos.txt
		}
	my $temp6 = "temp6.vitcos.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.vitcos.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "vitcos.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @vitcosdata = split (/\t/, $line);
		push (@taxnames, $vitcosdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Vitamins And Cofactors Matrix'\;\n\tDIMENSIONS NCHAR=17\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'vitamin requirement' \/  'vitamins not required' 'vitamins required', ";
	print OUT "2 'vitamin K1 requirement' \/  'vitamin K1 not required' 'vitamin K1 required', ";
	print OUT "3 'thiamin requirement' \/  'thiamin not required' 'thiamin required', ";
	print OUT "4 'vitamin B12 requirement' \/  'vitamin B21 not required' 'vitamin B12 required', ";
	print OUT "5 'growth factor requirement' \/  'growth factor not required' 'growth factor required', ";
	print OUT "6 'NAD as v-growth factor requirement' \/  'NAD not required' 'NAD required', ";
	print OUT "7 'hemin as x-growth factor requirement' \/  'hemin not required' 'hemin required', ";
	print OUT "8 'bile requirement' \/  'bile not required' 'bile required', ";
	print OUT "9 'ox-bile requirement' \/  'ox-bile not required' 'ox-bile required', ";
	print OUT "10 'serum requirement' \/  'serum not required' 'serum required', ";
	print OUT "11 'b vitamins' \/  'b vitamins not required' 'serum required', ";
	print OUT "12 'riboflavin requirement' \/  'riboflavin not required' 'riboflavin required', ";
	print OUT "13 'niacin requirement' \/  'niacin not required' 'niacin required', ";
	print OUT "14 'pantothenic acid requirement' \/  'pantothenic acid not required' 'pantothenic acid required', ";
	print OUT "15 'vitamin B6 requirement' \/  'vitamin B6 not required' 'vitamin B6 required', ";
	print OUT "16 'biotin requirement' \/  'biotin not required' 'biotin required', ";
	print OUT "17 'folic acid requirement' \/  'folic acid not required' 'folic acid required', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @vitcosdata = split (/\t/, $line);
		push (my @taxnames, $vitcosdata[0]);

#code char 1 vitamins
		if ($line =~ /,vitamins,|,vitamin k1,|,vitamin b12,|,thiamin,|,b vitamins,|,riboflavin,|,niacin,|,pantothenic acid,|,vitamin b6,|,biotin,|,folic acid,/) {
			print OUT @taxnames, "1";
			}
#		elsif ($line =~ /,not vitamins,/) {
#			print OUT @taxnames, "0";
#			}
		else {
			print OUT @taxnames, "0";
		}
#code char 2 vitamin k1
		if ($line =~ /,vitamin k1,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not vitamin k1,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 3 thiamin
		if ($line =~ /,thiamin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not thiamin,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 4 vitamin b12
		if ($line =~ /,vitamin b12,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not vitamin b12,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 5 growth factor
		if ($line =~ /,growth factor,|,nad,|,hemin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not growth factor,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 6 NAD
		if ($line =~ /,nad,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not nad,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 7 hemin
		if ($line =~ /,hemin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not hemin,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 8 bile
		if ($line =~ /,bile,|,ox-bile,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not bile,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 9 ox-bile
		if ($line =~ /,ox-bile,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not ox-bile,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 10 serum
		if ($line =~ /,serum,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not serum,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 11 B vitamins
		if ($line =~ /,b vitamins,|,thiamin,|,riboflavin,|,niacin,|,pantothenic acid,|,vitamin b6,|,biotin,|,folic acid,|,vitamin b12,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not b vitamins,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 12 riboflavin
		if ($line =~ /,riboflavin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not riboflavin,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 13 niacin
		if ($line =~ /,niacin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not niacin,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 14 pantothenic acid
		if ($line =~ /,pantothenic acid,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not pantothenic acid,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 15 vitamin B6
		if ($line =~ /,vitamin b6,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not vitamin b6,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 16 biotin
		if ($line =~ /,biotin,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not biotin,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}
#code char 17 folic acid
		if ($line =~ /,folic acid,/) {
			print OUT "1";
			}
#		elsif ($line =~ /,not folic acid,/) {
#			print OUT "0";
#			}
		else {
			print OUT "0";
		}



	print OUT "\n";

	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;

	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process antibiotic sensitivity character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "antibiotic sensitivity") {
	my $homout = "hom.antisens.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @antisens, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @antisens;

		map {s/[^[:print:]]//gi; } @antisens; # delete non-printable characters
		map {s/,\d g acriflavine ml1,/,acriflavine,/gi; } @antisens; # synonyms of 
		map {s/,acriflavine ml1,/,acriflavine,/gi; } @antisens; # synonyms of 
		map {s/,[\d|.]+% lysozyme,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,0.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,1000 g lysozyme ml−1,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,acetylspiramycin,/,spiramycin,/gi; } @antisens; # synonyms of 
		map {s/,actinomycin d,/,actinomycin,/gi; } @antisens; # synonyms of 
		map {s/,agent o.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,albamycin,/,novobiocin,/gi; } @antisens; # synonyms of 
		map {s/,aminoglycoside antibiotic tobramycin,/,tobramycin,/gi; } @antisens; # synonyms of 
		map {s/,amoxicillin\/clavulanic acid,/,amoxicillin,clavulanic acid,/gi; } @antisens; # synonyms of 
		map {s/,amoxycillin,/,amoxicillin,/gi; } @antisens; # synonyms of 
		map {s/,aniline green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,antibiotic tobramycin,/,tobramycin,/gi; } @antisens; # synonyms of 
		map {s/,aueromycin,/,chlortetracycline,/gi; } @antisens; # synonyms of 
		map {s/,augmentin,/,amoxicillin,clavulanic acid,/gi; } @antisens; # synonyms of 
		map {s/,aureomycin,/,chlortetracycline,/gi; } @antisens; # synonyms of 
		map {s/,benzaldehyde green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,benzylpenicilin,/,benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,benzylphenicillin,/,benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,carbenicillin 0,/,carbenicillin,/gi; } @antisens; # synonyms of 
		map {s/,cathomycin,/,novobiocin,/gi; } @antisens; # synonyms of 
		map {s/,cefaloridine,/,cephaloridine,/gi; } @antisens; # synonyms of 
		map {s/,cefalothin,/,cefalotin,/gi; } @antisens; # synonyms of 
		map {s/,cefazoline,/,cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,cefobid,/,cefoperazone,/gi; } @antisens; # synonyms of 
		map {s/,cephalexin,/,cefalexin,/gi; } @antisens; # synonyms of 
		map {s/,cephaloridin,/,cephaloridine,/gi; } @antisens; # synonyms of 
		map {s/,cephalosporins,/,cephalosporin,/gi; } @antisens; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/gi; } @antisens; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/gi; } @antisens; # synonyms of 
		map {s/,cephalotin,/,cefalotin,/gi; } @antisens; # synonyms of 
		map {s/,cephamezine vi,/,cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,cephamezine,/,cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,cephazolin,/,cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,Cetrimonium bromide,/,cetyltrimethylammonium bromide,/gi; } @antisens; # synonyms of 
		map {s/,chloramphenical,/,chloramphenicol,/gi; } @antisens; # synonyms of 
		map {s/,claventin,/,ticarcillin,clavulanic acid,/gi; } @antisens; # synonyms of 
		map {s/,co-trimoxazole,/,trimethoprim,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,compound o.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,cotrimoxazole,/,trimethoprim,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,CTAB,/,cetyltrimethylammonium bromide,/gi; } @antisens; # synonyms of 
		map {s/,dactinomycin,/,actinomycin,/gi; } @antisens; # synonyms of 
		map {s/,desferal,/,deferoxamine,/gi; } @antisens; # synonyms of 
		map {s/,diamond green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,emerald green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,fast green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,flagecidin,/,anisomycin,/gi; } @antisens; # synonyms of 
		map {s/,furadantin,/,nitrofurantoin,/gi; } @antisens; # ;synonyms of
		map {s/,gentam[iy]cin.a,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.b,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c1,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c1a,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c2,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.g,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.x,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentamycin,/,gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,gentian violet,/,crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,hexadecyltrimethylammonium bromide,/,cetyltrimethylammonium bromide,/gi; } @antisens; # synonyms of 
		map {s/,hexamethyl pararosaniline chloride,/,crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,imipinem,/,imipenem,/gi; } @antisens; # synonyms of 
		map {s/,isepamycin,/,isepamicin,/gi; } @antisens; # synonyms of 
		map {s/,kanamycin a,/,kanamycin,/gi; } @antisens; # synonyms of 
		map {s/,licomycin,/,lincomycin,/gi; } @antisens; # synonyms of 
		map {s/,lomefloxacin hydrochloride ,/,lomefloxacin,/gi; } @antisens; # synonyms of 
		map {s/,lysostaphin mic,/,lysostaphin,/gi; } @antisens; # synonyms of 
		map {s/,lysozome,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,lysozyme mic,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,malachite green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,mefoxin,/,cefoxitin,/gi; } @antisens; # synonyms of 
		map {s/,mepicycline penicillinate,/,penimepicycline,/gi; } @antisens; # synonyms of 
		map {s/,methyl violet 10b,/,crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,mynocycline,/,minocycline,/gi; } @antisens; # synonyms of 
		map {s/,nitrofurantoin 0,/,nitrofurantoin,/gi; } @antisens; # synonyms of 
		map {s/,novobiocin twenty strains,/,novobiocin,/gi; } @antisens; # synonyms of 
		map {s/,o.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,oleandomycin,/,oleandomycin,/gi; } @antisens; # synonyms of 
		map {s/,penecillin,/,penicillin,/gi; } @antisens; # synonyms of 
		map {s/,penicillin g mic,/,penicillin g,/gi; } @antisens; # synonyms of 
		map {s/,penicillin g u,/,penicillin g,/gi; } @antisens; # synonyms of 
		map {s/,penicillin g,/,benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,penicillin g.,/,penicillin g,/gi; } @antisens; # synonyms of 
		map {s/,penicillin v,/,phenoxymethylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,polyanetholsulfonic acid,/,sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,polym[iy]xin,/,polymyxins,/gi; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.m,/,polymyxin m,/gi; } @antisens; # synonyms of 
		map {s/,polymixin b.,/,polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,polymyxin b 0 u,/,polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,polymyxin b.,/,polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,polymyxin e,/,colistin,/gi; } @antisens; # synonyms of 
		map {s/,pristinamycin 11,/,pristinamycin,/gi; } @antisens; # synonyms of 
		map {s/,pristinamycine,/,pristinamycin,/gi; } @antisens; # synonyms of 
		map {s/,procaine penicillin,/,procaine benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,propamidine isothionate,/,propamidine,/gi; } @antisens; # synonyms of 
		map {s/,rifampin,/,rifampicin,/gi; } @antisens; # synonyms of 
		map {s/,smx,/,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,smz,/,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,sodium polyanetholsulfonate,/,sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,solid green,/,brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,spectinomycin 0,/,spectinomycin,/gi; } @antisens; # synonyms of 
		map {s/,sps,/,sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,streptomycin 2,/,streptomycin,/gi; } @antisens; # synonyms of 
		map {s/,sulfadiazin,/,sulfadiazine,/gi; } @antisens; # synonyms of 
		map {s/,sulfaisodimidine,/,sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,sulfamethin,/,sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,sulfamethoxazole-trimethoprim,/,sulfamethoxazole,trimethoprim,/gi; } @antisens; # synonyms of 
		map {s/,sulfamethoxazole\/trimethoprim,/,sulfamethoxazole,trimethoprim,/gi; } @antisens; # synonyms of 
		map {s/,sulfasomidine,/,sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,sulfuric diamide,/,sulfamides,/gi; } @antisens; # synonyms of 
		map {s/,tetracycline hydrochloride,/,tetracycline,/gi; } @antisens; # synonyms of 
		map {s/,tmp,/,trimethoprim,/gi; } @antisens; # synonyms of 
		map {s/,tmp\/smx,/,trimethoprim,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole 23.75,/,trimethoprim,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole,/,sulfamethoxazole,trimethoprim,/gi; } @antisens; # synonyms of 
		map {s/,vibriostatic agent 0.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,vibriostatic agent o.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,vibriostatic compound o.129,/,vibriostatic agent,/gi; } @antisens; # synonyms of 
		map {s/,bacitracine,/,bacitracin,/gi; } @antisens; # synonyms of 
		map {s/,[\d]+ g lysozyme ml1,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,lysozyme ml,/,lysozyme,/gi; } @antisens; # synonyms of 
		map {s/,63 phages,/,phages,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @antisens; # synonyms of 


		map {s/,not acetylspiramycin,/,not spiramycin,/gi; } @antisens; # synonyms of 
		map {s/,not actinomycin d,/,not actinomycin,/gi; } @antisens; # synonyms of 
		map {s/,not albamycin,/,not novobiocin,/gi; } @antisens; # synonyms of 
		map {s/,not amoxycillin,/,not amoxicillin,/gi; } @antisens; # synonyms of 
		map {s/,not aniline green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not aueromycin,/,not chlortetracycline,/gi; } @antisens; # synonyms of 
		map {s/,not aureomycin,/,not chlortetracycline,/gi; } @antisens; # synonyms of 
		map {s/,not benzaldehyde green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not benzylphenicillin,/,not benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,not cathomycin,/,not novobiocin,/gi; } @antisens; # synonyms of 
		map {s/,not cefaloridine,/,not cephaloridine,/gi; } @antisens; # synonyms of 
		map {s/,not cefazoline,/,not cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,not cefobid,/,not cefoperazone,/gi; } @antisens; # synonyms of 
		map {s/,not cephaloridin,/,not cephaloridine,/gi; } @antisens; # synonyms of 
		map {s/,not cephalosporins,/,not cephalosporin,/gi; } @antisens; # synonyms of 
		map {s/,not cephalothin,/,not cefalotin,/gi; } @antisens; # synonyms of 
		map {s/,not cephamezine vi,/,not cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,not cephamezine,/,not cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,not cephazolin,/,not cefazolin,/gi; } @antisens; # synonyms of 
		map {s/,not chloramphenical,/,not chloramphenicol,/gi; } @antisens; # synonyms of 
		map {s/,not co-trimoxazole,/,not trimethoprim,sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,not dactinomycin,/,not actinomycin,/gi; } @antisens; # synonyms of 
		map {s/,not diamond green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not emerald green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not fast green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not flagecidin,/,not anisomycin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.a,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.b,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c1,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c1a,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c2,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.g,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.x,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentamycin,/,not gentamicin,/gi; } @antisens; # synonyms of 
		map {s/,not gentian violet,/,not crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,not hexamethyl pararosaniline chloride,/,not crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,not kanamycin a,/,not kanamycin,/gi; } @antisens; # synonyms of 
		map {s/,not lomefloxacin hydrochloride ,/,not lomefloxacin,/gi; } @antisens; # synonyms of 
		map {s/,not malachite green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not mefoxin,/,not cefoxitin,/gi; } @antisens; # synonyms of 
		map {s/,not mepicycline penicillinate,/,not penimepicycline,/gi; } @antisens; # synonyms of 
		map {s/,not methyl violet 10b,/,not crystal violet,/gi; } @antisens; # synonyms of 
		map {s/,not oleandomycin,/,not oleandomycin,/gi; } @antisens; # synonyms of 
		map {s/,not penicillin g,/,not benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,not penicillin v,/,not phenoxymethylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,not polyanetholsulfonic acid,/,not sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin,/,not polymyxins,/gi; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/gi; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.m,/,not polymyxin m,/gi; } @antisens; # synonyms of 
		map {s/,not polymyxin e,/,not colistin,/gi; } @antisens; # synonyms of 
		map {s/,not procaine penicillin,/,not procaine benzylpenicillin,/gi; } @antisens; # synonyms of 
		map {s/,not rifampin,/,not rifampicin,/gi; } @antisens; # synonyms of 
		map {s/,not smx,/,not sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,not smz,/,not sulfamethoxazole,/gi; } @antisens; # synonyms of 
		map {s/,not sodium polyanetholsulfonate,/,not sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,not solid green,/,not brilliant green,/gi; } @antisens; # synonyms of 
		map {s/,not sps,/,not sodium polyanethol sulfonate,/gi; } @antisens; # synonyms of 
		map {s/,not sulfadiazin,/,not sulfadiazine,/gi; } @antisens; # synonyms of 
		map {s/,not sulfaisodimidine,/,not sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,not sulfamethin,/,not sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,not sulfasomidine,/,not sulfisomidine,/gi; } @antisens; # synonyms of 
		map {s/,not tmp,/,not trimethoprim,/gi; } @antisens; # synonyms of 

#not antibiotic terms
		map {s/,>[\d|]+ g,/,/gi; } @antisens; # synonyms of 
		map {s/,[\d|]+ g,/,/gi; } @antisens; # synonyms of 
		map {s/,[\d|]+ mg,/,/gi; } @antisens; # synonyms of 
		map {s/,.*%.*bile,/,/gi; } @antisens; # synonyms of 
		map {s/,.*%.*nacl,/,/gi; } @antisens; # synonyms of 
		map {s/,0.06,/,/gi; } @antisens; # synonyms of 
		map {s/,.*iu,/,/gi; } @antisens; # synonyms of 
		map {s/,.*oxgall,/,/gi; } @antisens; # synonyms of 
		map {s/,phages,/,/gi; } @antisens; # synonyms of 
		map {s/,23 s. aureus typing bacteriophages,/,/gi; } @antisens; # synonyms of 
		map {s/,[\d]+,/,/gi; } @antisens; # synonyms of 
		map {s/,[\d]+%,/,/gi; } @antisens; # synonyms of 
		map {s/,\Sg,/,/gi; } @antisens; # synonyms of 
		map {s/,0 iu,/,/gi; } @antisens; # synonyms of 
		map {s/,0.12.*,/,/gi; } @antisens; # synonyms of 
		map {s/,25.*,/,/gi; } @antisens; # synonyms of 
		map {s/,acid,/,/gi; } @antisens; # synonyms of 
		map {s/,addition,/,/gi; } @antisens; # synonyms of 
		map {s/,all,/,/gi; } @antisens; # synonyms of 
		map {s/,antibacterial compounds affecting protein synthesis,/,/gi; } @antisens; # synonyms of 
		map {s/,antibiotic resistant,/,/gi; } @antisens; # synonyms of 
		map {s/,antibiotic,/,/gi; } @antisens; # synonyms of 
		map {s/,antistaphylococcal antibiotics,/,/gi; } @antisens; # synonyms of 
		map {s/,discs,/,/gi; } @antisens; # synonyms of 
		map {s/,following antibiotics,/,/gi; } @antisens; # synonyms of 
		map {s/,following,/,/gi; } @antisens; # synonyms of 
		map {s/,g,/,/gi; } @antisens; # synonyms of 
		map {s/,heating,/,/gi; } @antisens; # synonyms of 
		map {s/,many antibiotics.*,/,/gi; } @antisens; # synonyms of 
		map {s/,mg,/,/gi; } @antisens; # synonyms of 
		map {s/,neo,/,/gi; } @antisens; # synonyms of 
		map {s/,other antibiotics,/,/gi; } @antisens; # synonyms of 
		map {s/,presence,/,/gi; } @antisens; # synonyms of 
		map {s/,same antibiotics,/,/gi; } @antisens; # synonyms of 
		map {s/,strains,/,/gi; } @antisens; # synonyms of 
		map {s/,this antibiotic,/,/gi; } @antisens; # synonyms of 
		map {s/,used antibiotics,/,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,/gi; } @antisens; # synonyms of 
		map {s/,xxx,/,/gi; } @antisens; # synonyms of 

		print OUT @taxlabels, @antisens, "\n"; # prints to $homout, hom.antisens.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.antisens.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,antibiotic sensitivity,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.antisens.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.antisens.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.antisens.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.antisens.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.antisens.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^actinomycin/**actinomycin/g; } @charlist2;
	map {s/^amikacin/**amikacin/g; } @charlist2;
	map {s/^aminoglycosides/**aminoglycosides/g; } @charlist2;
	map {s/^amoxicillin/**amoxicillin/g; } @charlist2;
	map {s/^amphotericin b/**amphotericin b/g; } @charlist2;
	map {s/^ampicillin/**ampicillin/g; } @charlist2;
	map {s/^anisomycin/**anisomycin/g; } @charlist2;
	map {s/^ansamycin/**ansamycin/g; } @charlist2;
	map {s/^antibiotic dyes/**antibiotic dyes/g; } @charlist2;
	map {s/^aphidicolin/**aphidicolin/g; } @charlist2;
	map {s/^bacitracin/**bacitracin/g; } @charlist2;
	map {s/^benzylpenicillin/**benzylpenicillin/g; } @charlist2;
	map {s/^beta-lactams/**beta-lactams/g; } @charlist2;
	map {s/^brilliant green/**brilliant green/g; } @charlist2;
	map {s/^carbenicillin/**carbenicillin/g; } @charlist2;
	map {s/^cefadroxil/**cefadroxil/g; } @charlist2;
	map {s/^cefalotin/**cefalotin/g; } @charlist2;
	map {s/^cefazolin/**cefazolin/g; } @charlist2;
	map {s/^cefoperazone/**cefoperazone/g; } @charlist2;
	map {s/^cefotaxime/**cefotaxime/g; } @charlist2;
	map {s/^cefoxitin/**cefoxitin/g; } @charlist2;
	map {s/^cefsulodin/**cefsulodin/g; } @charlist2;
	map {s/^ceftazidime/**ceftazidime/g; } @charlist2;
	map {s/^cefuroxime/**cefuroxime/g; } @charlist2;
	map {s/^cephaloridine/**cephaloridine/g; } @charlist2;
	map {s/^cephalosporin/**cephalosporin/g; } @charlist2;
	map {s/^cephalosporins/**cephalosporins/g; } @charlist2;
	map {s/^chloramphenicol/**chloramphenicol/g; } @charlist2;
	map {s/^chlortetracycline/**chlortetracycline/g; } @charlist2;
	map {s/^ciprofloxacin/**ciprofloxacin/g; } @charlist2;
	map {s/^clavulanic acid/**clavulanic acid/g; } @charlist2;
	map {s/^clindamycin/**clindamycin/g; } @charlist2;
	map {s/^clomocycline/**clomocycline/g; } @charlist2;
	map {s/^colistin/**colistin/g; } @charlist2;
	map {s/^crystal violet/**crystal violet/g; } @charlist2;
	map {s/^cyclic peptides/**cyclic peptides/g; } @charlist2;
	map {s/^cycloheximide/**cycloheximide/g; } @charlist2;
	map {s/^demeclocycline/**demeclocycline/g; } @charlist2;
	map {s/^dihydrostreptomycin/**dihydrostreptomycin/g; } @charlist2;
	map {s/^diterpenes/**diterpenes/g; } @charlist2;
	map {s/^doxycycline/**doxycycline/g; } @charlist2;
	map {s/^erythromycin/**erythromycin/g; } @charlist2;
	map {s/^fluoroquinolones/**fluoroquinolones/g; } @charlist2;
	map {s/^furazolidone/**furazolidone/g; } @charlist2;
	map {s/^fusidic acid/**fusidic acid/g; } @charlist2;
	map {s/^gentamicin/**gentamicin/g; } @charlist2;
	map {s/^glycopeptides/**glycopeptides/g; } @charlist2;
	map {s/^imipenem/**imipenem/g; } @charlist2;
	map {s/^kanamycin/**kanamycin/g; } @charlist2;
	map {s/^levofloxacin/**levofloxacin/g; } @charlist2;
	map {s/^lincomycin/**lincomycin/g; } @charlist2;
	map {s/^lomefloxacin/**lomefloxacin/g; } @charlist2;
	map {s/^lymecycline/**lymecycline/g; } @charlist2;
	map {s/^macrolides/**macrolides/g; } @charlist2;
	map {s/^meclocycline/**meclocycline/g; } @charlist2;
	map {s/^metacycline/**metacycline/g; } @charlist2;
	map {s/^metronidazole/**metronidazole/g; } @charlist2;
	map {s/^minocycline/**minocycline/g; } @charlist2;
	map {s/^nalidixic acid/**nalidixic acid/g; } @charlist2;
	map {s/^neomycin/**neomycin/g; } @charlist2;
	map {s/^nitrofurans/**nitrofurans/g; } @charlist2;
	map {s/^nitrofurantoin/**nitrofurantoin/g; } @charlist2;
	map {s/^nitroimidazole/**nitroimidazole/g; } @charlist2;
	map {s/^norfloxacin/**norfloxacin/g; } @charlist2;
	map {s/^novobiocin/**novobiocin/g; } @charlist2;
	map {s/^nucleosides/**nucleosides/g; } @charlist2;
	map {s/^oleandomycin/**oleandomycin/g; } @charlist2;
	map {s/^omadacycline/**omadacycline/g; } @charlist2;
	map {s/^organochlorine antibiotics/**organochlorine antibiotics/g; } @charlist2;
	map {s/^oxacillin/**oxacillin/g; } @charlist2;
	map {s/^oxytetracycline/**oxytetracycline/g; } @charlist2;
	map {s/^penicillin/**penicillin/g; } @charlist2;
	map {s/^penicillins/**penicillins/g; } @charlist2;
	map {s/^penimepicycline/**penimepicycline/g; } @charlist2;
	map {s/^phenoxymethylpenicillin/**phenoxymethylpenicillin/g; } @charlist2;
	map {s/^piperacillin/**piperacillin/g; } @charlist2;
	map {s/^piperazines/**piperazines/g; } @charlist2;
	map {s/^piperidines/**piperidines/g; } @charlist2;
	map {s/^polyketides/**polyketides/g; } @charlist2;
	map {s/^polymyxin b/**polymyxin b/g; } @charlist2;
	map {s/^polymyxin m/**polymyxin m/g; } @charlist2;
	map {s/^polymyxins/**polymyxins/g; } @charlist2;
	map {s/^polypeptides/**polypeptides/g; } @charlist2;
	map {s/^procaine benzylpenicillin/**procaine benzylpenicillin/g; } @charlist2;
	map {s/^puromycin/**puromycin/g; } @charlist2;
	map {s/^pyrimidiness/**pyrimidiness/g; } @charlist2;
	map {s/^pyrrolidines/**pyrrolidines/g; } @charlist2;
	map {s/^quinolones/**quinolones/g; } @charlist2;
	map {s/^rifampicin/**rifampicin/g; } @charlist2;
	map {s/^rolitetracycline/**rolitetracycline/g; } @charlist2;
	map {s/^roxithromycin/**roxithromycin/g; } @charlist2;
	map {s/^sodium polyanethol sulfonate/**sodium polyanethol sulfonate/g; } @charlist2;
	map {s/^spectinomycin/**spectinomycin/g; } @charlist2;
	map {s/^spiramycin/**spiramycin/g; } @charlist2;
	map {s/^steroid antibiotics/**steroid antibiotics/g; } @charlist2;
	map {s/^streptomycin/**streptomycin/g; } @charlist2;
	map {s/^sulfadiazine/**sulfadiazine/g; } @charlist2;
	map {s/^sulfamethoxazole/**sulfamethoxazole/g; } @charlist2;
	map {s/^sulfisomidine/**sulfisomidine/g; } @charlist2;
	map {s/^sulfonamides/**sulfonamides/g; } @charlist2;
	map {s/^tetracycline/**tetracycline/g; } @charlist2;
	map {s/^tetracyclines/**tetracyclines/g; } @charlist2;
	map {s/^ticarcillin/**ticarcillin/g; } @charlist2;
	map {s/^tobramycin/**tobramycin/g; } @charlist2;
	map {s/^trimethoprim/**trimethoprim/g; } @charlist2;
	map {s/^vancomycin/**vancomycin/g; } @charlist2;
#new additions
	map {s/^acriflavine/**acriflavine/g; } @charlist2;
	map {s/^arbekacin/**arbekacin/g; } @charlist2;
	map {s/^azlocillin/**azlocillin/g; } @charlist2;
	map {s/^aztreonam/**aztreonam/g; } @charlist2;
	map {s/^bacteriocin/**bacteriocin/g; } @charlist2;
	map {s/^boron/**boron/g; } @charlist2;
	map {s/^cadmium acetate/**cadmium acetate/g; } @charlist2;
	map {s/^cefalexin/**cefalexin/g; } @charlist2;
	map {s/^cefalotin/**cefalotin/g; } @charlist2;
	map {s/^cefamandole/**cefamandole/g; } @charlist2;
	map {s/^cefmetazole/**cefmetazole/g; } @charlist2;
	map {s/^cefotiam/**cefotiam/g; } @charlist2;
	map {s/^ceftiofur/**ceftiofur/g; } @charlist2;
	map {s/^ceftriaxone/**ceftriaxone/g; } @charlist2;
	map {s/^cetyltrimethylammonium bromide/**cetyltrimethylammonium bromide/g; } @charlist2;
	map {s/^chlorhexidine digluconate/**chlorhexidine digluconate/g; } @charlist2;
	map {s/^cloxacillin/**cloxacillin/g; } @charlist2;
	map {s/^deferoxamine/**deferoxamine/g; } @charlist2;
	map {s/^efrotomycin/**efrotomycin/g; } @charlist2;
	map {s/^enrofloxacin/**enrofloxacin/g; } @charlist2;
	map {s/^florfenicol/**florfenicol/g; } @charlist2;
	map {s/^flumequine/**flumequine/g; } @charlist2;
	map {s/^fosfomycin/**fosfomycin/g; } @charlist2;
	map {s/^isepamicin/**isepamicin/g; } @charlist2;
	map {s/^josamycin/**josamycin/g; } @charlist2;
	map {s/^linezolid/**linezolid/g; } @charlist2;
	map {s/^lysostaphin/**lysostaphin/g; } @charlist2;
	map {s/^lysozyme/**lysozyme/g; } @charlist2;
	map {s/^medicamycin/**medicamycin/g; } @charlist2;
	map {s/^mercuric nitrate/**mercuric nitrate/g; } @charlist2;
	map {s/^meropenem/**meropenem/g; } @charlist2;
	map {s/^methicillin/**methicillin/g; } @charlist2;
	map {s/^mupirocin/**mupirocin/g; } @charlist2;
	map {s/^netilmicin/**netilmicin/g; } @charlist2;
	map {s/^nitrofurazolidone/**nitrofurazolidone/g; } @charlist2;
	map {s/^nystatin/**nystatin/g; } @charlist2;
	map {s/^ofloxacin/**ofloxacin/g; } @charlist2;
	map {s/^optochin/**optochin/g; } @charlist2;
	map {s/^oxolinic acid/**oxolinic acid/g; } @charlist2;
	map {s/^pefloxacin/**pefloxacin/g; } @charlist2;
	map {s/^phages/**phages/g; } @charlist2;
	map {s/^phenylmethylsulfonyl fluoride/**phenylmethylsulfonyl fluoride/g; } @charlist2;
	map {s/^pristinamycin/**pristinamycin/g; } @charlist2;
	map {s/^propamidine/**propamidine/g; } @charlist2;
	map {s/^rifamycin sv/**rifamycin sv/g; } @charlist2;
	map {s/^rifamycin/**rifamycin/g; } @charlist2;
	map {s/^sulfamides/**sulfamides/g; } @charlist2; 
	map {s/^sulfisoxazole/**sulfisoxazole/g; } @charlist2;
	map {s/^teicoplanin/**teicoplanin/g; } @charlist2;
	map {s/^tigecycline/**tigecycline/g; } @charlist2;
	map {s/^troleandomycin/**troleandomycin/g; } @charlist2;
	map {s/^vibriostatic agent/**vibriostatic agent/g; } @charlist2;


	map {s/^xxx/**xxx/g; } @charlist2;
	map {s/^xxx/**xxx/g; } @charlist2;
	map {s/^xxx/**xxx/g; } @charlist2;
	map {s/^xxx/**xxx/g; } @charlist2;
	map {s/^xxx/**xxx/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.antisens.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.antisens.txt
		}
	my $temp6 = "temp6.antisens.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.antisens.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "antisens.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @antisensdata = split (/\t/, $line);
		push (@taxnames, $antisensdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Antibiotic Sensitivity Matrix'\;\n\tDIMENSIONS NCHAR=158\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 actinomycin \/  'not sensitive' sensitive, ";
	print OUT "2 amikacin \/  'not sensitive' sensitive, ";
	print OUT "3 amoxicillin \/  'not sensitive' sensitive, ";
	print OUT "4 'amphotericin B' \/  'not sensitive' sensitive, ";
	print OUT "5 ampicillin \/  'not sensitive' sensitive, ";
	print OUT "6 anisomycin \/  'not sensitive' sensitive, ";
	print OUT "7 ansamycin \/  'not sensitive' sensitive, ";
	print OUT "8 aphidicolin \/  'not sensitive' sensitive, ";
	print OUT "9 bacitracin \/  'not sensitive' sensitive, ";
	print OUT "10 benzylpenicillin \/  'not sensitive' sensitive, ";
	print OUT "11 'brilliant green' \/  'not sensitive' sensitive, ";
	print OUT "12 carbenicillin \/  'not sensitive' sensitive, ";
	print OUT "13 cefadroxil \/  'not sensitive' sensitive, ";
	print OUT "14 cefalotin \/  'not sensitive' sensitive, ";
	print OUT "15 cefazolin \/  'not sensitive' sensitive, ";
	print OUT "16 cefoperazone \/  'not sensitive' sensitive, ";
	print OUT "17 cefotaxime \/  'not sensitive' sensitive, ";
	print OUT "18 cefoxitin \/  'not sensitive' sensitive, ";
	print OUT "19 cefsulodin \/  'not sensitive' sensitive, ";
	print OUT "20 ceftazidime \/  'not sensitive' sensitive, ";
	print OUT "21 cefuroxime \/  'not sensitive' sensitive, ";
	print OUT "22 cephaloridine \/  'not sensitive' sensitive, ";
	print OUT "23 cephalosporin \/  'not sensitive' sensitive, ";
	print OUT "24 chloramphenicol \/  'not sensitive' sensitive, ";
	print OUT "25 chlortetracycline \/  'not sensitive' sensitive, ";
	print OUT "26 ciprofloxacin \/  'not sensitive' sensitive, ";
	print OUT "27 'clavulanic acid' \/  'not sensitive' sensitive, ";
	print OUT "28 clindamycin \/  'not sensitive' sensitive, ";
	print OUT "29 clomocycline \/  'not sensitive' sensitive, ";
	print OUT "30 colistin \/  'not sensitive' sensitive, ";
	print OUT "31 'crystal violet' \/  'not sensitive' sensitive, ";
	print OUT "32 cycloheximide \/  'not sensitive' sensitive, ";
	print OUT "33 demeclocycline \/  'not sensitive' sensitive, ";
	print OUT "34 dihydrostreptomycin \/  'not sensitive' sensitive, ";
	print OUT "35 doxycycline \/  'not sensitive' sensitive, ";
	print OUT "36 erythromycin \/  'not sensitive' sensitive, ";
	print OUT "37 furazolidone \/  'not sensitive' sensitive, ";
	print OUT "38 'fusidic acid' \/  'not sensitive' sensitive, ";
	print OUT "39 gentamicin \/  'not sensitive' sensitive, ";
	print OUT "40 imipenem \/  'not sensitive' sensitive, ";
	print OUT "41 kanamycin \/  'not sensitive' sensitive, ";
	print OUT "42 levofloxacin \/  'not sensitive' sensitive, ";
	print OUT "43 lincomycin \/  'not sensitive' sensitive, ";
	print OUT "44 lomefloxacin \/  'not sensitive' sensitive, ";
	print OUT "45 lymecycline \/  'not sensitive' sensitive, ";
	print OUT "46 meclocycline \/  'not sensitive' sensitive, ";
	print OUT "47 metacycline \/  'not sensitive' sensitive, ";
	print OUT "48 metronidazole \/  'not sensitive' sensitive, ";
	print OUT "49 minocycline \/  'not sensitive' sensitive, ";
	print OUT "50 'nalidixic acid' \/  'not sensitive' sensitive, ";
	print OUT "51 neomycin \/  'not sensitive' sensitive, ";
	print OUT "52 nitrofurantoin \/  'not sensitive' sensitive, ";
	print OUT "53 norfloxacin \/  'not sensitive' sensitive, ";
	print OUT "54 novobiocin \/  'not sensitive' sensitive, ";
	print OUT "55 oleandomycin \/  'not sensitive' sensitive, ";
	print OUT "56 omadacycline \/  'not sensitive' sensitive, ";
	print OUT "57 oxacillin \/  'not sensitive' sensitive, ";
	print OUT "58 oxytetracycline \/  'not sensitive' sensitive, ";
	print OUT "59 penicillin \/  'not sensitive' sensitive, ";
	print OUT "60 penicillins \/  'not sensitive' sensitive, ";
	print OUT "61 penimepicycline \/  'not sensitive' sensitive, ";
	print OUT "62 phenoxymethylpenicillin \/  'not sensitive' sensitive, ";
	print OUT "63 piperacillin \/  'not sensitive' sensitive, ";
	print OUT "64 'polymyxin B' \/  'not sensitive' sensitive, ";
	print OUT "65 'polymyxin M' \/  'not sensitive' sensitive, ";
	print OUT "66 polymyxins \/  'not sensitive' sensitive, ";
	print OUT "67 'procaine benzylpenicillin' \/  'not sensitive' sensitive, ";
	print OUT "68 puromycin \/  'not sensitive' sensitive, ";
	print OUT "69 rifampicin \/  'not sensitive' sensitive, ";
	print OUT "70 rolitetracycline \/  'not sensitive' sensitive, ";
	print OUT "71 roxithromycin \/  'not sensitive' sensitive, ";
	print OUT "72 'sodium polyanethol sulfonate' \/  'not sensitive' sensitive, ";
	print OUT "73 spectinomycin \/  'not sensitive' sensitive, ";
	print OUT "74 spiramycin \/  'not sensitive' sensitive, ";
	print OUT "75 streptomycin \/  'not sensitive' sensitive, ";
	print OUT "76 sulfadiazine \/  'not sensitive' sensitive, ";
	print OUT "77 sulfamethoxazole \/  'not sensitive' sensitive, ";
	print OUT "78 sulfisomidine \/  'not sensitive' sensitive, ";
	print OUT "79 tetracycline \/  'not sensitive' sensitive, ";
	print OUT "80 ticarcillin \/  'not sensitive' sensitive, ";
	print OUT "81 tobramycin \/  'not sensitive' sensitive, ";
	print OUT "82 trimethoprim \/  'not sensitive' sensitive, ";
	print OUT "83 vancomycin \/  'not sensitive' sensitive, ";
	print OUT "84 'aminoglycoside antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "85 'antibiotic dyes' \/  'not sensitive' sensitive, ";
	print OUT "86 'beta-lactam antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "87 'cephalosporin antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "88 'cyclic peptide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "89 'diterpene antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "90 'fluoroquinolone antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "91 'glycopeptide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "92 'macrolide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "93 'nitrofuran antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "94 'nitroimidazole antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "95 'nucleoside antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "96 'organochlorine antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "97 'piperazine antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "98 'piperidine antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "99 'polyketide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "100 'polypeptide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "101 'pyrimidine antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "102 'pyrrolidine antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "103 'quinolone antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "104 'steroid antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "105 'sulfonamide antibiotics' \/  'not sensitive' sensitive, ";
	print OUT "106 'tetracycline antibiotics' \/  'not sensitive' sensitive, ";

	print OUT "107 acriflavine \/  'not sensitive' sensitive, ";
	print OUT "108 arbekacin \/  'not sensitive' sensitive, ";
	print OUT "109 azlocillin \/  'not sensitive' sensitive, ";
	print OUT "110 aztreonam \/  'not sensitive' sensitive, ";
	print OUT "111 bacteriocin \/  'not sensitive' sensitive, ";
	print OUT "112 boron \/  'not sensitive' sensitive, ";
	print OUT "113 'cadmium acetate' \/  'not sensitive' sensitive, ";
	print OUT "114 cefalexin \/  'not sensitive' sensitive, ";
	print OUT "115 cefalotin \/  'not sensitive' sensitive, ";
	print OUT "116 cefamandole \/  'not sensitive' sensitive, ";
	print OUT "117 cefmetazole \/  'not sensitive' sensitive, ";
	print OUT "118 cefotiam \/  'not sensitive' sensitive, ";
	print OUT "129 ceftiofur \/  'not sensitive' sensitive, ";
	print OUT "120 ceftriaxone \/  'not sensitive' sensitive, ";
	print OUT "121 'cetyltrimethylammonium bromide' \/  'not sensitive' sensitive, ";
	print OUT "122 'chlorhexidine digluconate' \/  'not sensitive' sensitive, ";
	print OUT "123 cloxacillin \/  'not sensitive' sensitive, ";
	print OUT "124 deferoxamine \/  'not sensitive' sensitive, ";
	print OUT "125 efrotomycin \/  'not sensitive' sensitive, ";
	print OUT "126 enrofloxacin \/  'not sensitive' sensitive, ";
	print OUT "127 florfenicol \/  'not sensitive' sensitive, ";
	print OUT "128 flumequine \/  'not sensitive' sensitive, ";
	print OUT "139 fosfomycin \/  'not sensitive' sensitive, ";
	print OUT "130 isepamicin \/  'not sensitive' sensitive, ";
	print OUT "131 josamycin \/  'not sensitive' sensitive, ";
	print OUT "132 linezolid \/  'not sensitive' sensitive, ";
	print OUT "133 lysostaphin \/  'not sensitive' sensitive, ";
	print OUT "134 lysozyme \/  'not sensitive' sensitive, ";
	print OUT "135 medicamycin \/  'not sensitive' sensitive, ";
	print OUT "136 'mercuric nitrate' \/  'not sensitive' sensitive, ";
	print OUT "137 meropenem \/  'not sensitive' sensitive, ";
	print OUT "138 methicillin \/  'not sensitive' sensitive, ";
	print OUT "149 mupirocin \/  'not sensitive' sensitive, ";
	print OUT "140 netilmicin \/  'not sensitive' sensitive, ";
	print OUT "141 nitrofurazolidone \/  'not sensitive' sensitive, ";
	print OUT "142 nystatin \/  'not sensitive' sensitive, ";
	print OUT "143 ofloxacin \/  'not sensitive' sensitive, ";
	print OUT "144 optochin \/  'not sensitive' sensitive, ";
	print OUT "145 'oxolinic acid' \/  'not sensitive' sensitive, ";
	print OUT "146 pefloxacin \/  'not sensitive' sensitive, ";
	print OUT "147 phages \/  'not sensitive' sensitive, ";
	print OUT "148 'phenylmethylsulfonyl fluoride' \/  'not sensitive' sensitive, ";
	print OUT "159 pristinamycin \/  'not sensitive' sensitive, ";
	print OUT "150 propamidine \/  'not sensitive' sensitive, ";
	print OUT "151 'rifamycin sv' \/  'not sensitive' sensitive, ";
	print OUT "152 rifamycin \/  'not sensitive' sensitive, ";
	print OUT "153 sulfamides \/  'not sensitive' sensitive, ";
	print OUT "154 sulfisoxazole \/  'not sensitive' sensitive, ";
	print OUT "155 teicoplanin \/  'not sensitive' sensitive, ";
	print OUT "156 tigecycline \/  'not sensitive' sensitive, ";
	print OUT "157 troleandomycin \/  'not sensitive' sensitive, ";
	print OUT "158 'vibriostatic agent' \/  'not sensitive' sensitive, ";


	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @antisensdata = split (/\t/, $line);
		push (my @taxnames, $antisensdata[0]);

#code char 1 actinomycin
		if ($line =~ /,actinomycin,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not actinomycin,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 amikacin
		if ($line =~ /,amikacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amikacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 amoxicillin
		if ($line =~ /,amoxicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amoxicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 amphotericin b
		if ($line =~ /,amphotericin b,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphotericin b,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 ampicillin
		if ($line =~ /,ampicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ampicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 anisomycin
		if ($line =~ /,anisomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not anisomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 ansamycin
		if ($line =~ /,ansamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ansamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 aphidicolin
		if ($line =~ /,aphidicolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aphidicolin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 bacitracin
		if ($line =~ /,bacitracin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not bacitracin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 benzylpenicillin
		if ($line =~ /,benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not benzylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 brilliant green
		if ($line =~ /,brilliant green,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not brilliant green,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 carbenicillin
		if ($line =~ /,carbenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carbenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 13 cefadroxil
		if ($line =~ /,cefadroxil,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefadroxil,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 14 cefalotin
		if ($line =~ /,cefalotin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalotin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 cefazolin
		if ($line =~ /,cefazolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefazolin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 cefoperazone
		if ($line =~ /,cefoperazone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefoperazone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 17 cefotaxime
		if ($line =~ /,cefotaxime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefotaxime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 18 cefoxitin
		if ($line =~ /,cefoxitin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefoxitin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 19 cefsulodin
		if ($line =~ /,cefsulodin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefsulodin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 20 ceftazidime
		if ($line =~ /,ceftazidime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftazidime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 21 cefuroxime
		if ($line =~ /,cefuroxime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefuroxime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 22 cephaloridine
		if ($line =~ /,cephaloridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephaloridine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 23 cephalosporin
		if ($line =~ /,cephalosporin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephalosporin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 24 chloramphenicol
		if ($line =~ /,chloramphenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chloramphenicol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 25 chlortetracycline
		if ($line =~ /,chlortetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chlortetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 26 ciprofloxacin
		if ($line =~ /,ciprofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ciprofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 27 clavulanic acid
		if ($line =~ /,clavulanic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clavulanic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 28 clindamycin
		if ($line =~ /,clindamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clindamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 29 clomocycline
		if ($line =~ /,clomocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clomocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 30 colistin
		if ($line =~ /,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not colistin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 31 crystal violet
		if ($line =~ /,crystal violet,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not crystal violet,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 32 cycloheximide
		if ($line =~ /,cycloheximide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cycloheximide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 33 demeclocycline
		if ($line =~ /,demeclocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not demeclocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 34 dihydrostreptomycin
		if ($line =~ /,dihydrostreptomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dihydrostreptomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 35 doxycycline
		if ($line =~ /,doxycycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not doxycycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 36 erythromycin
		if ($line =~ /,erythromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erythromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 37 furazolidone
		if ($line =~ /,furazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not furazolidone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 38 fusidic acid
		if ($line =~ /,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fusidic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 39 gentamicin
		if ($line =~ /,gentamicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gentamicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 40 imipenem
		if ($line =~ /,imipenem,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not imipenem,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 41 kanamycin
		if ($line =~ /,kanamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not kanamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 42 levofloxacin
		if ($line =~ /,levofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not levofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 43 lincomycin
		if ($line =~ /,lincomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lincomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 44 lomefloxacin
		if ($line =~ /,lomefloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lomefloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 45 lymecycline
		if ($line =~ /,lymecycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lymecycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 46 meclocycline
		if ($line =~ /,meclocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meclocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 47 metacycline
		if ($line =~ /,metacycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not metacycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 48 metronidazole
		if ($line =~ /,metronidazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not metronidazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 49 xminocyclinexx
		if ($line =~ /,minocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not minocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 50 nalidixic acid
		if ($line =~ /,nalidixic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nalidixic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 51 neomycin
		if ($line =~ /,neomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not neomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 52 nitrofurantoin
		if ($line =~ /,nitrofurantoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurantoin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 53 norfloxacin
		if ($line =~ /,norfloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not norfloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 54 novobiocin
		if ($line =~ /,novobiocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not novobiocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 55 oleandomycin
		if ($line =~ /,oleandomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oleandomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 56 omadacycline
		if ($line =~ /,omadacycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not omadacycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 57 oxacillin
		if ($line =~ /,oxacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 58 oxytetracycline
		if ($line =~ /,oxytetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxytetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 59 penicillin
		if ($line =~ /,penicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 60 penicillins
		if ($line =~ /,penicillins,|,phenoxymethylpenicillin,|,benzylpenicillin,|,penicillin,|,procaine benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penicillins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 61 penimepicycline
		if ($line =~ /,penimepicycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penimepicycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 62 phenoxymethylpenicillin
		if ($line =~ /,phenoxymethylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenoxymethylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 63 piperacillin
		if ($line =~ /,piperacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 64 polymyxin b
		if ($line =~ /,polymyxin b,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxin b,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 65 polymyxin m
		if ($line =~ /,polymyxin m,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxin m,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 66 polymyxins
		if ($line =~ /,polymyxins,|,polymyxin m,|,polymyxin b,|,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 67 procaine benzylpenicillin
		if ($line =~ /,procaine benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not procaine benzylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 68 puromycin
		if ($line =~ /,puromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not puromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 69 rifampicin
		if ($line =~ /,rifampicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifampicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 70 rolitetracycline
		if ($line =~ /,rolitetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rolitetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 71 roxithromycin
		if ($line =~ /,roxithromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not roxithromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 72 sodium polyanethol sulfonate
		if ($line =~ /,sodium polyanethol sulfonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sodium polyanethol sulfonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 73 spectinomycin
		if ($line =~ /,spectinomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spectinomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 74 spiramycin
		if ($line =~ /,spiramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spiramycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 75 streptomycin
		if ($line =~ /,streptomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not streptomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 76 sulfadiazine
		if ($line =~ /,sulfadiazine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfadiazine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 77 sulfamethoxazole
		if ($line =~ /,sulfamethoxazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfamethoxazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 78 sulfisomidine
		if ($line =~ /,sulfisomidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfisomidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 79 tetracycline
		if ($line =~ /,tetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 80 ticarcillin
		if ($line =~ /,ticarcillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ticarcillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 81 tobramycin
		if ($line =~ /,tobramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tobramycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 82 trimethoprim
		if ($line =~ /,trimethoprim,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trimethoprim,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 83 vancomycin
		if ($line =~ /,vancomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not vancomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 84 aminoglycoside antibiotics
		if ($line =~ /,aminoglycosides,|,spectinomycin,|,dihydrostreptomycin,|,gentamicin,|,kanamycin,|,neomycin,|,streptomycin,|,tobramycin,|,amikacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aminoglycosides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 85 antibiotic dyes
		if ($line =~ /,antibiotic dyes,|,crystal violet,|,brilliant green,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not antibiotic dyes,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 86 beta-lactams
		if ($line =~ /,beta-lactams,|,amoxicillin,|,ampicillin,|,benzylpenicillin,|,carbenicillin,|,cephalosporin,|,clavulanic acid,|,imipenem,|,oxacillin,|,penicillin,|,phenoxymethylpenicillin,|,procaine benzylpenicillin,|,ticarcillin,|,cefoxitin,|,penicillins,|,piperacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not beta-lactams,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 87 cephalosporins
		if ($line =~ /,cephalosporins,|,cefadroxil,|,cefalotin,|,cefazolin,|,cefoperazone,|,cefotaxime,|,cefsulodin,|,ceftazidime,|,cefuroxime,|,cephaloridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephalosporins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 88 cyclic peptides
		if ($line =~ /,cyclic peptides,|,polymyxin b,|,polymyxin m,|,polymyxins,|,bacitracin,|,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cyclic peptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 89 diterpenes
		if ($line =~ /,diterpenes,|,aphidicolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diterpenes,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 90 fluoroquinolones
		if ($line =~ /,fluoroquinolones,|,levofloxacin,|,lomefloxacin,|,norfloxacin,|,ciprofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fluoroquinolones,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 91 glycopeptides
		if ($line =~ /,glycopeptides,|,vancomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycopeptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 92 macrolide antibiotics
		if ($line =~ /,macrolides,|,erythromycin,|,oleandomycin,|,roxithromycin,|,spiramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not macrolides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 93 nitrofuran antibiotics
		if ($line =~ /,nitrofurans,|,furazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurans,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 94 nitroimidazole antibiotics
		if ($line =~ /,nitroimidazole,|,metronidazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitroimidazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 95 nucleoside antibiotics
		if ($line =~ /,nucleosides,|,puromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nucleosides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 96 organochlorine antibiotics
		if ($line =~ /,organochlorine antibiotics,|,acriflavine,|,florfenicol,|,clindamycin,|,lincomycin,|,chloramphenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not organochlorine antibiotics,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 97 piperazines
		if ($line =~ /,piperazines,|,rifampicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperazines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 98 piperidines
		if ($line =~ /,piperidines,|,cycloheximide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperidines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 99 polyketides
		if ($line =~ /,polyketides,|,amphotericin b,|,ansamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polyketides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 100 polypeptide antibiotics
		if ($line =~ /,polypeptides,|,actinomycin,|,lysozyme,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polypeptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 101 pyrimidine antibiotics
		if ($line =~ /,pyrimidiness,|,trimethoprim,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrimidiness,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 102 pyrrolidine antibiotics
		if ($line =~ /,pyrrolidines,|,anisomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrrolidines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 103 quinolone antibiotics
		if ($line =~ /,quinolones,|,nalidixic acid,|,novobiocin,|,nitrofurantoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not quinolones,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 104 steroid antibiotics
		if ($line =~ /,steroid antibiotics,|,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not steroid antibiotics,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 105 sulfonamides
		if ($line =~ /,sulfonamides,|,sulfadiazine,|,sulfamethoxazole,|,sulfisomidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfonamides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 106 tetracyclines
		if ($line =~ /,tetracyclines,|,chlortetracycline,|,clomocycline,|,demeclocycline,|,doxycycline,|,lymecycline,|,meclocycline,|,metacycline,|,minocycline,|,omadacycline,|,oxytetracycline,|,penimepicycline,|,rolitetracycline,|,tetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetracyclines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

#code char 107 acriflavine
		if ($line =~ /,acriflavine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acriflavine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 108 arbekacin
		if ($line =~ /,arbekacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arbekacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 109 azlocillin
		if ($line =~ /,azlocillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not azlocillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 110 aztreonam
		if ($line =~ /,aztreonam,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aztreonam,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 111 bacteriocin
		if ($line =~ /,bacteriocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not bacteriocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 112 boron
		if ($line =~ /,boron,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not boron,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 113 cadmium acetate
		if ($line =~ /,cadmium acetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cadmium acetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 114 cefalexin
		if ($line =~ /,cefalexin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalexin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 115 cefalotin
		if ($line =~ /,cefalotin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalotin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 116 cefamandole
		if ($line =~ /,cefamandole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefamandole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 117 cefmetazole
		if ($line =~ /,cefmetazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefmetazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 118 cefotiam
		if ($line =~ /,cefotiam,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefotiam,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 119 ceftiofur
		if ($line =~ /,ceftiofur,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftiofur,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 120 ceftriaxone
		if ($line =~ /,ceftriaxone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftriaxone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 121 cetyltrimethylammonium bromide
		if ($line =~ /,cetyltrimethylammonium bromide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cetyltrimethylammonium bromide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 122 chlorhexidine digluconate
		if ($line =~ /,chlorhexidine digluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chlorhexidine digluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 123 cloxacillin
		if ($line =~ /,cloxacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cloxacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 124 deferoxamine
		if ($line =~ /,deferoxamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not deferoxamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 125 efrotomycin
		if ($line =~ /,efrotomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not efrotomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 126 enrofloxacin
		if ($line =~ /,enrofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not enrofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 127 florfenicol
		if ($line =~ /,florfenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not florfenicol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 128 flumequine
		if ($line =~ /,flumequine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flumequine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 129 fosfomycin
		if ($line =~ /,fosfomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fosfomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 130 isepamicin
		if ($line =~ /,isepamicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isepamicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 131 josamycin
		if ($line =~ /,josamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not josamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 132 linezolid
		if ($line =~ /,linezolid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not linezolid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 133 lysostaphin
		if ($line =~ /,lysostaphin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysostaphin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 134 lysozyme
		if ($line =~ /,lysozyme,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysozyme,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 135 medicamycin
		if ($line =~ /,medicamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not medicamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 136 mercuric nitrate
		if ($line =~ /,mercuric nitrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mercuric nitrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 137 meropenem
		if ($line =~ /,meropenem,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meropenem,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 138 methicillin
		if ($line =~ /,methicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 139 mupirocin
		if ($line =~ /,mupirocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mupirocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 140 netilmicin
		if ($line =~ /,netilmicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not netilmicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 141 nitrofurazolidone
		if ($line =~ /,nitrofurazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurazolidone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 142 nystatin
		if ($line =~ /,nystatin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nystatin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 143 ofloxacin
		if ($line =~ /,ofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 144 optochin
		if ($line =~ /,optochin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not optochin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 145 oxolinic acid
		if ($line =~ /,oxolinic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxolinic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 146 pefloxacin
		if ($line =~ /,pefloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pefloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 147 phages
		if ($line =~ /,phages,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phages,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 148 phenylmethylsulfonyl fluoride
		if ($line =~ /,phenylmethylsulfonyl fluoride,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylmethylsulfonyl fluoride,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 149 pristinamycin
		if ($line =~ /,pristinamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pristinamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 150 propamidine
		if ($line =~ /,propamidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propamidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 151 rifamycin sv
		if ($line =~ /,rifamycin sv,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifamycin sv,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 152 rifamycin
		if ($line =~ /,rifamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 153 sulfamides
		if ($line =~ /,sulfamides,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfamides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 154 sulfisoxazole
		if ($line =~ /,sulfisoxazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfisoxazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 155 teicoplanin
		if ($line =~ /,teicoplanin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not teicoplanin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 156 tigecycline
		if ($line =~ /,tigecycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tigecycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 157 troleandomycin
		if ($line =~ /,troleandomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not troleandomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 158 vibriostatic agent
		if ($line =~ /,vibriostatic agent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not vibriostatic agent,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}


	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process antibiotic resistant character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "antibiotic resistant") {
	my $homout = "hom.antiresis.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @antiresis, $columns[1]);
		$_ = lc for @antiresis;

		map {s/[^[:print:]]//gi; } @antiresis; # delete non-printable characters
		map {s/,\d g acriflavine ml1,/,acriflavine,/gi; } @antiresis; # synonyms of 
		map {s/,acriflavine ml1,/,acriflavine,/gi; } @antiresis; # synonyms of 
		map {s/,[\d|.]+% lysozyme,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,0.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,1000 g lysozyme ml−1,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,acetylspiramycin,/,spiramycin,/gi; } @antiresis; # synonyms of 
		map {s/,actinomycin d,/,actinomycin,/gi; } @antiresis; # synonyms of 
		map {s/,agent o.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,albamycin,/,novobiocin,/gi; } @antiresis; # synonyms of 
		map {s/,aminoglycoside antibiotic tobramycin,/,tobramycin,/gi; } @antiresis; # synonyms of 
		map {s/,amoxicillin\/clavulanic acid,/,amoxicillin,clavulanic acid,/gi; } @antiresis; # synonyms of 
		map {s/,amoxycillin,/,amoxicillin,/gi; } @antiresis; # synonyms of 
		map {s/,aniline green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,antibiotic tobramycin,/,tobramycin,/gi; } @antiresis; # synonyms of 
		map {s/,aueromycin,/,chlortetracycline,/gi; } @antiresis; # synonyms of 
		map {s/,augmentin,/,amoxicillin,clavulanic acid,/gi; } @antiresis; # synonyms of 
		map {s/,aureomycin,/,chlortetracycline,/gi; } @antiresis; # synonyms of 
		map {s/,benzaldehyde green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,benzylpenicilin,/,benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,benzylphenicillin,/,benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,carbenicillin 0,/,carbenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,cathomycin,/,novobiocin,/gi; } @antiresis; # synonyms of 
		map {s/,cefaloridine,/,cephaloridine,/gi; } @antiresis; # synonyms of 
		map {s/,cefalothin,/,cefalotin,/gi; } @antiresis; # synonyms of 
		map {s/,cefazoline,/,cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,cefobid,/,cefoperazone,/gi; } @antiresis; # synonyms of 
		map {s/,cephalexin,/,cefalexin,/gi; } @antiresis; # synonyms of 
		map {s/,cephaloridin,/,cephaloridine,/gi; } @antiresis; # synonyms of 
		map {s/,cephalosporins,/,cephalosporin,/gi; } @antiresis; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/gi; } @antiresis; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/gi; } @antiresis; # synonyms of 
		map {s/,cephalotin,/,cefalotin,/gi; } @antiresis; # synonyms of 
		map {s/,cephamezine vi,/,cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,cephamezine,/,cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,cephazolin,/,cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,Cetrimonium bromide,/,cetyltrimethylammonium bromide,/gi; } @antiresis; # synonyms of 
		map {s/,chloramphenical,/,chloramphenicol,/gi; } @antiresis; # synonyms of 
		map {s/,claventin,/,ticarcillin,clavulanic acid,/gi; } @antiresis; # synonyms of 
		map {s/,co-trimoxazole,/,trimethoprim,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,compound o.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,cotrimoxazole,/,trimethoprim,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,CTAB,/,cetyltrimethylammonium bromide,/gi; } @antiresis; # synonyms of 
		map {s/,dactinomycin,/,actinomycin,/gi; } @antiresis; # synonyms of 
		map {s/,desferal,/,deferoxamine,/gi; } @antiresis; # synonyms of 
		map {s/,diamond green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,emerald green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,fast green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,flagecidin,/,anisomycin,/gi; } @antiresis; # synonyms of 
		map {s/,furadantin,/,nitrofurantoin,/gi; } @antiresis; # ;synonyms of
		map {s/,gentam[iy]cin.a,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.b,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c1,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c1a,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c2,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.g,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.x,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentamycin,/,gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,gentian violet,/,crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,hexadecyltrimethylammonium bromide,/,cetyltrimethylammonium bromide,/gi; } @antiresis; # synonyms of 
		map {s/,hexamethyl pararosaniline chloride,/,crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,imipinem,/,imipenem,/gi; } @antiresis; # synonyms of 
		map {s/,isepamycin,/,isepamicin,/gi; } @antiresis; # synonyms of 
		map {s/,kanamycin a,/,kanamycin,/gi; } @antiresis; # synonyms of 
		map {s/,licomycin,/,lincomycin,/gi; } @antiresis; # synonyms of 
		map {s/,lomefloxacin hydrochloride ,/,lomefloxacin,/gi; } @antiresis; # synonyms of 
		map {s/,lysostaphin mic,/,lysostaphin,/gi; } @antiresis; # synonyms of 
		map {s/,lysozome,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,lysozyme mic,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,malachite green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,mefoxin,/,cefoxitin,/gi; } @antiresis; # synonyms of 
		map {s/,mepicycline penicillinate,/,penimepicycline,/gi; } @antiresis; # synonyms of 
		map {s/,methyl violet 10b,/,crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,mynocycline,/,minocycline,/gi; } @antiresis; # synonyms of 
		map {s/,nitrofurantoin 0,/,nitrofurantoin,/gi; } @antiresis; # synonyms of 
		map {s/,novobiocin twenty strains,/,novobiocin,/gi; } @antiresis; # synonyms of 
		map {s/,o.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,oleandomycin,/,oleandomycin,/gi; } @antiresis; # synonyms of 
		map {s/,penecillin,/,penicillin,/gi; } @antiresis; # synonyms of 
		map {s/,penicillin g mic,/,penicillin g,/gi; } @antiresis; # synonyms of 
		map {s/,penicillin g u,/,penicillin g,/gi; } @antiresis; # synonyms of 
		map {s/,penicillin g,/,benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,penicillin g.,/,penicillin g,/gi; } @antiresis; # synonyms of 
		map {s/,penicillin v,/,phenoxymethylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,polyanetholsulfonic acid,/,sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin,/,polymyxins,/gi; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.m,/,polymyxin m,/gi; } @antiresis; # synonyms of 
		map {s/,polymixin b.,/,polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,polymyxin b 0 u,/,polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,polymyxin b.,/,polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,polymyxin e,/,colistin,/gi; } @antiresis; # synonyms of 
		map {s/,pristinamycin 11,/,pristinamycin,/gi; } @antiresis; # synonyms of 
		map {s/,pristinamycine,/,pristinamycin,/gi; } @antiresis; # synonyms of 
		map {s/,procaine penicillin,/,procaine benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,propamidine isothionate,/,propamidine,/gi; } @antiresis; # synonyms of 
		map {s/,rifampin,/,rifampicin,/gi; } @antiresis; # synonyms of 
		map {s/,smx,/,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,smz,/,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,sodium polyanetholsulfonate,/,sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,solid green,/,brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,spectinomycin 0,/,spectinomycin,/gi; } @antiresis; # synonyms of 
		map {s/,sps,/,sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,streptomycin 2,/,streptomycin,/gi; } @antiresis; # synonyms of 
		map {s/,sulfadiazin,/,sulfadiazine,/gi; } @antiresis; # synonyms of 
		map {s/,sulfaisodimidine,/,sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,sulfamethin,/,sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,sulfamethoxazole-trimethoprim,/,sulfamethoxazole,trimethoprim,/gi; } @antiresis; # synonyms of 
		map {s/,sulfamethoxazole\/trimethoprim,/,sulfamethoxazole,trimethoprim,/gi; } @antiresis; # synonyms of 
		map {s/,sulfasomidine,/,sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,sulfuric diamide,/,sulfamides,/gi; } @antiresis; # synonyms of 
		map {s/,tetracycline hydrochloride,/,tetracycline,/gi; } @antiresis; # synonyms of 
		map {s/,tmp,/,trimethoprim,/gi; } @antiresis; # synonyms of 
		map {s/,tmp\/smx,/,trimethoprim,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole 23.75,/,trimethoprim,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole,/,sulfamethoxazole,trimethoprim,/gi; } @antiresis; # synonyms of 
		map {s/,vibriostatic agent 0.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,vibriostatic agent o.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,vibriostatic compound o.129,/,vibriostatic agent,/gi; } @antiresis; # synonyms of 
		map {s/,bacitracine,/,bacitracin,/gi; } @antiresis; # synonyms of 
		map {s/,[\d]+ g lysozyme ml1,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,lysozyme ml,/,lysozyme,/gi; } @antiresis; # synonyms of 
		map {s/,63 phages,/,phages,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @antiresis; # synonyms of 


		map {s/,not acetylspiramycin,/,not spiramycin,/gi; } @antiresis; # synonyms of 
		map {s/,not actinomycin d,/,not actinomycin,/gi; } @antiresis; # synonyms of 
		map {s/,not albamycin,/,not novobiocin,/gi; } @antiresis; # synonyms of 
		map {s/,not amoxycillin,/,not amoxicillin,/gi; } @antiresis; # synonyms of 
		map {s/,not aniline green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not aueromycin,/,not chlortetracycline,/gi; } @antiresis; # synonyms of 
		map {s/,not aureomycin,/,not chlortetracycline,/gi; } @antiresis; # synonyms of 
		map {s/,not benzaldehyde green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not benzylphenicillin,/,not benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,not cathomycin,/,not novobiocin,/gi; } @antiresis; # synonyms of 
		map {s/,not cefaloridine,/,not cephaloridine,/gi; } @antiresis; # synonyms of 
		map {s/,not cefazoline,/,not cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,not cefobid,/,not cefoperazone,/gi; } @antiresis; # synonyms of 
		map {s/,not cephaloridin,/,not cephaloridine,/gi; } @antiresis; # synonyms of 
		map {s/,not cephalosporins,/,not cephalosporin,/gi; } @antiresis; # synonyms of 
		map {s/,not cephalothin,/,not cefalotin,/gi; } @antiresis; # synonyms of 
		map {s/,not cephamezine vi,/,not cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,not cephamezine,/,not cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,not cephazolin,/,not cefazolin,/gi; } @antiresis; # synonyms of 
		map {s/,not chloramphenical,/,not chloramphenicol,/gi; } @antiresis; # synonyms of 
		map {s/,not co-trimoxazole,/,not trimethoprim,sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,not dactinomycin,/,not actinomycin,/gi; } @antiresis; # synonyms of 
		map {s/,not diamond green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not emerald green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not fast green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not flagecidin,/,not anisomycin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.a,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.b,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c1,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c1a,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c2,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.g,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.x,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentamycin,/,not gentamicin,/gi; } @antiresis; # synonyms of 
		map {s/,not gentian violet,/,not crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,not hexamethyl pararosaniline chloride,/,not crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,not kanamycin a,/,not kanamycin,/gi; } @antiresis; # synonyms of 
		map {s/,not lomefloxacin hydrochloride ,/,not lomefloxacin,/gi; } @antiresis; # synonyms of 
		map {s/,not malachite green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not mefoxin,/,not cefoxitin,/gi; } @antiresis; # synonyms of 
		map {s/,not mepicycline penicillinate,/,not penimepicycline,/gi; } @antiresis; # synonyms of 
		map {s/,not methyl violet 10b,/,not crystal violet,/gi; } @antiresis; # synonyms of 
		map {s/,not oleandomycin,/,not oleandomycin,/gi; } @antiresis; # synonyms of 
		map {s/,not penicillin g,/,not benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,not penicillin v,/,not phenoxymethylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,not polyanetholsulfonic acid,/,not sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin,/,not polymyxins,/gi; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/gi; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.m,/,not polymyxin m,/gi; } @antiresis; # synonyms of 
		map {s/,not polymyxin e,/,not colistin,/gi; } @antiresis; # synonyms of 
		map {s/,not procaine penicillin,/,not procaine benzylpenicillin,/gi; } @antiresis; # synonyms of 
		map {s/,not rifampin,/,not rifampicin,/gi; } @antiresis; # synonyms of 
		map {s/,not smx,/,not sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,not smz,/,not sulfamethoxazole,/gi; } @antiresis; # synonyms of 
		map {s/,not sodium polyanetholsulfonate,/,not sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,not solid green,/,not brilliant green,/gi; } @antiresis; # synonyms of 
		map {s/,not sps,/,not sodium polyanethol sulfonate,/gi; } @antiresis; # synonyms of 
		map {s/,not sulfadiazin,/,not sulfadiazine,/gi; } @antiresis; # synonyms of 
		map {s/,not sulfaisodimidine,/,not sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,not sulfamethin,/,not sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,not sulfasomidine,/,not sulfisomidine,/gi; } @antiresis; # synonyms of 
		map {s/,not tmp,/,not trimethoprim,/gi; } @antiresis; # synonyms of 

#not antibiotic terms
		map {s/,>[\d|]+ g,/,/gi; } @antiresis; # synonyms of 
		map {s/,[\d|]+ g,/,/gi; } @antiresis; # synonyms of 
		map {s/,[\d|]+ mg,/,/gi; } @antiresis; # synonyms of 
		map {s/,.*%.*bile,/,/gi; } @antiresis; # synonyms of 
		map {s/,.*%.*nacl,/,/gi; } @antiresis; # synonyms of 
		map {s/,0.06,/,/gi; } @antiresis; # synonyms of 
		map {s/,.*iu,/,/gi; } @antiresis; # synonyms of 
		map {s/,.*oxgall,/,/gi; } @antiresis; # synonyms of 
		map {s/,phages,/,/gi; } @antiresis; # synonyms of 
		map {s/,23 s. aureus typing bacteriophages,/,/gi; } @antiresis; # synonyms of 
		map {s/,[\d]+,/,/gi; } @antiresis; # synonyms of 
		map {s/,[\d]+%,/,/gi; } @antiresis; # synonyms of 
		map {s/,\Sg,/,/gi; } @antiresis; # synonyms of 
		map {s/,0 iu,/,/gi; } @antiresis; # synonyms of 
		map {s/,0.12.*,/,/gi; } @antiresis; # synonyms of 
		map {s/,25.*,/,/gi; } @antiresis; # synonyms of 
		map {s/,acid,/,/gi; } @antiresis; # synonyms of 
		map {s/,addition,/,/gi; } @antiresis; # synonyms of 
		map {s/,all,/,/gi; } @antiresis; # synonyms of 
		map {s/,antibacterial compounds affecting protein synthesis,/,/gi; } @antiresis; # synonyms of 
		map {s/,antibiotic resistant,/,/gi; } @antiresis; # synonyms of 
		map {s/,antibiotic,/,/gi; } @antiresis; # synonyms of 
		map {s/,antistaphylococcal antibiotics,/,/gi; } @antiresis; # synonyms of 
		map {s/,discs,/,/gi; } @antiresis; # synonyms of 
		map {s/,following antibiotics,/,/gi; } @antiresis; # synonyms of 
		map {s/,following,/,/gi; } @antiresis; # synonyms of 
		map {s/,g,/,/gi; } @antiresis; # synonyms of 
		map {s/,heating,/,/gi; } @antiresis; # synonyms of 
		map {s/,many antibiotics.*,/,/gi; } @antiresis; # synonyms of 
		map {s/,mg,/,/gi; } @antiresis; # synonyms of 
		map {s/,neo,/,/gi; } @antiresis; # synonyms of 
		map {s/,other antibiotics,/,/gi; } @antiresis; # synonyms of 
		map {s/,presence,/,/gi; } @antiresis; # synonyms of 
		map {s/,same antibiotics,/,/gi; } @antiresis; # synonyms of 
		map {s/,strains,/,/gi; } @antiresis; # synonyms of 
		map {s/,this antibiotic,/,/gi; } @antiresis; # synonyms of 
		map {s/,used antibiotics,/,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,/gi; } @antiresis; # synonyms of 
		map {s/,xxx,/,/gi; } @antiresis; # synonyms of 

		print OUT @taxlabels, @antiresis, "\n"; # prints to $homout, hom.antiresis.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.antiresis.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,antibiotic sensitivity,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.antiresis.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.antiresis.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.antiresis.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.antiresis.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.antiresis.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^actinomycin/**actinomycin/g; } @charlist2;
	map {s/^amikacin/**amikacin/g; } @charlist2;
	map {s/^aminoglycosides/**aminoglycosides/g; } @charlist2;
	map {s/^amoxicillin/**amoxicillin/g; } @charlist2;
	map {s/^amphotericin b/**amphotericin b/g; } @charlist2;
	map {s/^ampicillin/**ampicillin/g; } @charlist2;
	map {s/^anisomycin/**anisomycin/g; } @charlist2;
	map {s/^ansamycin/**ansamycin/g; } @charlist2;
	map {s/^antibiotic dyes/**antibiotic dyes/g; } @charlist2;
	map {s/^aphidicolin/**aphidicolin/g; } @charlist2;
	map {s/^bacitracin/**bacitracin/g; } @charlist2;
	map {s/^benzylpenicillin/**benzylpenicillin/g; } @charlist2;
	map {s/^beta-lactams/**beta-lactams/g; } @charlist2;
	map {s/^brilliant green/**brilliant green/g; } @charlist2;
	map {s/^carbenicillin/**carbenicillin/g; } @charlist2;
	map {s/^cefadroxil/**cefadroxil/g; } @charlist2;
	map {s/^cefalotin/**cefalotin/g; } @charlist2;
	map {s/^cefazolin/**cefazolin/g; } @charlist2;
	map {s/^cefoperazone/**cefoperazone/g; } @charlist2;
	map {s/^cefotaxime/**cefotaxime/g; } @charlist2;
	map {s/^cefoxitin/**cefoxitin/g; } @charlist2;
	map {s/^cefsulodin/**cefsulodin/g; } @charlist2;
	map {s/^ceftazidime/**ceftazidime/g; } @charlist2;
	map {s/^cefuroxime/**cefuroxime/g; } @charlist2;
	map {s/^cephaloridine/**cephaloridine/g; } @charlist2;
	map {s/^cephalosporin/**cephalosporin/g; } @charlist2;
	map {s/^cephalosporins/**cephalosporins/g; } @charlist2;
	map {s/^chloramphenicol/**chloramphenicol/g; } @charlist2;
	map {s/^chlortetracycline/**chlortetracycline/g; } @charlist2;
	map {s/^ciprofloxacin/**ciprofloxacin/g; } @charlist2;
	map {s/^clavulanic acid/**clavulanic acid/g; } @charlist2;
	map {s/^clindamycin/**clindamycin/g; } @charlist2;
	map {s/^clomocycline/**clomocycline/g; } @charlist2;
	map {s/^colistin/**colistin/g; } @charlist2;
	map {s/^crystal violet/**crystal violet/g; } @charlist2;
	map {s/^cyclic peptides/**cyclic peptides/g; } @charlist2;
	map {s/^cycloheximide/**cycloheximide/g; } @charlist2;
	map {s/^demeclocycline/**demeclocycline/g; } @charlist2;
	map {s/^dihydrostreptomycin/**dihydrostreptomycin/g; } @charlist2;
	map {s/^diterpenes/**diterpenes/g; } @charlist2;
	map {s/^doxycycline/**doxycycline/g; } @charlist2;
	map {s/^erythromycin/**erythromycin/g; } @charlist2;
	map {s/^fluoroquinolones/**fluoroquinolones/g; } @charlist2;
	map {s/^furazolidone/**furazolidone/g; } @charlist2;
	map {s/^fusidic acid/**fusidic acid/g; } @charlist2;
	map {s/^gentamicin/**gentamicin/g; } @charlist2;
	map {s/^glycopeptides/**glycopeptides/g; } @charlist2;
	map {s/^imipenem/**imipenem/g; } @charlist2;
	map {s/^kanamycin/**kanamycin/g; } @charlist2;
	map {s/^levofloxacin/**levofloxacin/g; } @charlist2;
	map {s/^lincomycin/**lincomycin/g; } @charlist2;
	map {s/^lomefloxacin/**lomefloxacin/g; } @charlist2;
	map {s/^lymecycline/**lymecycline/g; } @charlist2;
	map {s/^macrolides/**macrolides/g; } @charlist2;
	map {s/^meclocycline/**meclocycline/g; } @charlist2;
	map {s/^metacycline/**metacycline/g; } @charlist2;
	map {s/^metronidazole/**metronidazole/g; } @charlist2;
	map {s/^minocycline/**minocycline/g; } @charlist2;
	map {s/^nalidixic acid/**nalidixic acid/g; } @charlist2;
	map {s/^neomycin/**neomycin/g; } @charlist2;
	map {s/^nitrofurans/**nitrofurans/g; } @charlist2;
	map {s/^nitrofurantoin/**nitrofurantoin/g; } @charlist2;
	map {s/^nitroimidazole/**nitroimidazole/g; } @charlist2;
	map {s/^norfloxacin/**norfloxacin/g; } @charlist2;
	map {s/^novobiocin/**novobiocin/g; } @charlist2;
	map {s/^nucleosides/**nucleosides/g; } @charlist2;
	map {s/^oleandomycin/**oleandomycin/g; } @charlist2;
	map {s/^omadacycline/**omadacycline/g; } @charlist2;
	map {s/^organochlorine antibiotics/**organochlorine antibiotics/g; } @charlist2;
	map {s/^oxacillin/**oxacillin/g; } @charlist2;
	map {s/^oxytetracycline/**oxytetracycline/g; } @charlist2;
	map {s/^penicillin/**penicillin/g; } @charlist2;
	map {s/^penicillins/**penicillins/g; } @charlist2;
	map {s/^penimepicycline/**penimepicycline/g; } @charlist2;
	map {s/^phenoxymethylpenicillin/**phenoxymethylpenicillin/g; } @charlist2;
	map {s/^piperacillin/**piperacillin/g; } @charlist2;
	map {s/^piperazines/**piperazines/g; } @charlist2;
	map {s/^piperidines/**piperidines/g; } @charlist2;
	map {s/^polyketides/**polyketides/g; } @charlist2;
	map {s/^polymyxin b/**polymyxin b/g; } @charlist2;
	map {s/^polymyxin m/**polymyxin m/g; } @charlist2;
	map {s/^polymyxins/**polymyxins/g; } @charlist2;
	map {s/^polypeptides/**polypeptides/g; } @charlist2;
	map {s/^procaine benzylpenicillin/**procaine benzylpenicillin/g; } @charlist2;
	map {s/^puromycin/**puromycin/g; } @charlist2;
	map {s/^pyrimidiness/**pyrimidiness/g; } @charlist2;
	map {s/^pyrrolidines/**pyrrolidines/g; } @charlist2;
	map {s/^quinolones/**quinolones/g; } @charlist2;
	map {s/^rifampicin/**rifampicin/g; } @charlist2;
	map {s/^rolitetracycline/**rolitetracycline/g; } @charlist2;
	map {s/^roxithromycin/**roxithromycin/g; } @charlist2;
	map {s/^sodium polyanethol sulfonate/**sodium polyanethol sulfonate/g; } @charlist2;
	map {s/^spectinomycin/**spectinomycin/g; } @charlist2;
	map {s/^spiramycin/**spiramycin/g; } @charlist2;
	map {s/^steroid antibiotics/**steroid antibiotics/g; } @charlist2;
	map {s/^streptomycin/**streptomycin/g; } @charlist2;
	map {s/^sulfadiazine/**sulfadiazine/g; } @charlist2;
	map {s/^sulfamethoxazole/**sulfamethoxazole/g; } @charlist2;
	map {s/^sulfisomidine/**sulfisomidine/g; } @charlist2;
	map {s/^sulfonamides/**sulfonamides/g; } @charlist2;
	map {s/^tetracycline/**tetracycline/g; } @charlist2;
	map {s/^tetracyclines/**tetracyclines/g; } @charlist2;
	map {s/^ticarcillin/**ticarcillin/g; } @charlist2;
	map {s/^tobramycin/**tobramycin/g; } @charlist2;
	map {s/^trimethoprim/**trimethoprim/g; } @charlist2;
	map {s/^vancomycin/**vancomycin/g; } @charlist2;
#new additions
	map {s/^acriflavine/**acriflavine/g; } @charlist2;
	map {s/^arbekacin/**arbekacin/g; } @charlist2;
	map {s/^azlocillin/**azlocillin/g; } @charlist2;
	map {s/^aztreonam/**aztreonam/g; } @charlist2;
	map {s/^bacteriocin/**bacteriocin/g; } @charlist2;
	map {s/^boron/**boron/g; } @charlist2;
	map {s/^cadmium acetate/**cadmium acetate/g; } @charlist2;
	map {s/^cefalexin/**cefalexin/g; } @charlist2;
	map {s/^cefalotin/**cefalotin/g; } @charlist2;
	map {s/^cefamandole/**cefamandole/g; } @charlist2;
	map {s/^cefmetazole/**cefmetazole/g; } @charlist2;
	map {s/^cefotiam/**cefotiam/g; } @charlist2;
	map {s/^ceftiofur/**ceftiofur/g; } @charlist2;
	map {s/^ceftriaxone/**ceftriaxone/g; } @charlist2;
	map {s/^cetyltrimethylammonium bromide/**cetyltrimethylammonium bromide/g; } @charlist2;
	map {s/^chlorhexidine digluconate/**chlorhexidine digluconate/g; } @charlist2;
	map {s/^cloxacillin/**cloxacillin/g; } @charlist2;
	map {s/^deferoxamine/**deferoxamine/g; } @charlist2;
	map {s/^deferoxamine/**deferoxamine/g; } @charlist2;
	map {s/^efrotomycin/**efrotomycin/g; } @charlist2;
	map {s/^enrofloxacin/**enrofloxacin/g; } @charlist2;
	map {s/^enrofloxacin/**enrofloxacin/g; } @charlist2;
	map {s/^florfenicol/**florfenicol/g; } @charlist2;
	map {s/^flumequine/**flumequine/g; } @charlist2;
	map {s/^fosfomycin/**fosfomycin/g; } @charlist2;
	map {s/^fosfomycin/**fosfomycin/g; } @charlist2;
	map {s/^isepamicin/**isepamicin/g; } @charlist2;
	map {s/^josamycin/**josamycin/g; } @charlist2;
	map {s/^linezolid/**linezolid/g; } @charlist2;
	map {s/^lysostaphin/**lysostaphin/g; } @charlist2;
	map {s/^lysostaphin/**lysostaphin/g; } @charlist2;
	map {s/^lysozyme/**lysozyme/g; } @charlist2;
	map {s/^medicamycin/**medicamycin/g; } @charlist2;
	map {s/^mercuric nitrate/**mercuric nitrate/g; } @charlist2;
	map {s/^meropenem/**meropenem/g; } @charlist2;
	map {s/^methicillin/**methicillin/g; } @charlist2;
	map {s/^mupirocin/**mupirocin/g; } @charlist2;
	map {s/^netilmicin/**netilmicin/g; } @charlist2;
	map {s/^nitrofurazolidone/**nitrofurazolidone/g; } @charlist2;
	map {s/^nystatin/**nystatin/g; } @charlist2;
	map {s/^ofloxacin/**ofloxacin/g; } @charlist2;
	map {s/^optochin/**optochin/g; } @charlist2;
	map {s/^oxolinic acid/**oxolinic acid/g; } @charlist2;
	map {s/^pefloxacin/**pefloxacin/g; } @charlist2;
	map {s/^phages/**phages/g; } @charlist2;
	map {s/^phenylmethylsulfonyl fluoride/**phenylmethylsulfonyl fluoride/g; } @charlist2;
	map {s/^pristinamycin/**pristinamycin/g; } @charlist2;
	map {s/^propamidine/**propamidine/g; } @charlist2;
	map {s/^rifamycin sv/**rifamycin sv/g; } @charlist2;
	map {s/^rifamycin/**rifamycin/g; } @charlist2;
	map {s/^sulfamides/**sulfamides/g; } @charlist2; 
	map {s/^sulfisoxazole/**sulfisoxazole/g; } @charlist2;
	map {s/^teicoplanin/**teicoplanin/g; } @charlist2;
	map {s/^tigecycline/**tigecycline/g; } @charlist2;
	map {s/^troleandomycin/**troleandomycin/g; } @charlist2;
	map {s/^vibriostatic agent/**vibriostatic agent/g; } @charlist2;


	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.antiresis.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.antiresis.txt
		}
	my $temp6 = "temp6.antiresis.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.antiresis.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "antiresis.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @antiresisdata = split (/\t/, $line);
		push (@taxnames, $antiresisdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE' Antibiotic Resistance Matrix'\;\n\tDIMENSIONS NCHAR=158\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 actinomycin \/  sensitive resistant, ";
	print OUT "2 amikacin \/  sensitive resistant, ";
	print OUT "3 amoxicillin \/  sensitive resistant, ";
	print OUT "4 'amphotericin B' \/  sensitive resistant, ";
	print OUT "5 ampicillin \/  sensitive resistant, ";
	print OUT "6 anisomycin \/  sensitive resistant, ";
	print OUT "7 ansamycin \/  sensitive resistant, ";
	print OUT "8 aphidicolin \/  sensitive resistant, ";
	print OUT "9 bacitracin \/  sensitive resistant, ";
	print OUT "10 benzylpenicillin \/  sensitive resistant, ";
	print OUT "11 'brilliant green' \/  sensitive resistant, ";
	print OUT "12 carbenicillin \/  sensitive resistant, ";
	print OUT "13 cefadroxil \/  sensitive resistant, ";
	print OUT "14 cefalotin \/  sensitive resistant, ";
	print OUT "15 cefazolin \/  sensitive resistant, ";
	print OUT "16 cefoperazone \/  sensitive resistant, ";
	print OUT "17 cefotaxime \/  sensitive resistant, ";
	print OUT "18 cefoxitin \/  sensitive resistant, ";
	print OUT "19 cefsulodin \/  sensitive resistant, ";
	print OUT "20 ceftazidime \/  sensitive resistant, ";
	print OUT "21 cefuroxime \/  sensitive resistant, ";
	print OUT "22 cephaloridine \/  sensitive resistant, ";
	print OUT "23 cephalosporin \/  sensitive resistant, ";
	print OUT "24 chloramphenicol \/  sensitive resistant, ";
	print OUT "25 chlortetracycline \/  sensitive resistant, ";
	print OUT "26 ciprofloxacin \/  sensitive resistant, ";
	print OUT "27 'clavulanic acid' \/  sensitive resistant, ";
	print OUT "28 clindamycin \/  sensitive resistant, ";
	print OUT "29 clomocycline \/  sensitive resistant, ";
	print OUT "30 colistin \/  sensitive resistant, ";
	print OUT "31 'crystal violet' \/  sensitive resistant, ";
	print OUT "32 cycloheximide \/  sensitive resistant, ";
	print OUT "33 demeclocycline \/  sensitive resistant, ";
	print OUT "34 dihydrostreptomycin \/  sensitive resistant, ";
	print OUT "35 doxycycline \/  sensitive resistant, ";
	print OUT "36 erythromycin \/  sensitive resistant, ";
	print OUT "37 furazolidone \/  sensitive resistant, ";
	print OUT "38 'fusidic acid' \/  sensitive resistant, ";
	print OUT "39 gentamicin \/  sensitive resistant, ";
	print OUT "40 imipenem \/  sensitive resistant, ";
	print OUT "41 kanamycin \/  sensitive resistant, ";
	print OUT "42 levofloxacin \/  sensitive resistant, ";
	print OUT "43 lincomycin \/  sensitive resistant, ";
	print OUT "44 lomefloxacin \/  sensitive resistant, ";
	print OUT "45 lymecycline \/  sensitive resistant, ";
	print OUT "46 meclocycline \/  sensitive resistant, ";
	print OUT "47 metacycline \/  sensitive resistant, ";
	print OUT "48 metronidazole \/  sensitive resistant, ";
	print OUT "49 minocycline \/  sensitive resistant, ";
	print OUT "50 'nalidixic acid' \/  sensitive resistant, ";
	print OUT "51 neomycin \/  sensitive resistant, ";
	print OUT "52 nitrofurantoin \/  sensitive resistant, ";
	print OUT "53 norfloxacin \/  sensitive resistant, ";
	print OUT "54 novobiocin \/  sensitive resistant, ";
	print OUT "55 oleandomycin \/  sensitive resistant, ";
	print OUT "56 omadacycline \/  sensitive resistant, ";
	print OUT "57 oxacillin \/  sensitive resistant, ";
	print OUT "58 oxytetracycline \/  sensitive resistant, ";
	print OUT "59 penicillin \/  sensitive resistant, ";
	print OUT "60 penicillins \/  sensitive resistant, ";
	print OUT "61 penimepicycline \/  sensitive resistant, ";
	print OUT "62 phenoxymethylpenicillin \/  sensitive resistant, ";
	print OUT "63 piperacillin \/  sensitive resistant, ";
	print OUT "64 'polymyxin B' \/  sensitive resistant, ";
	print OUT "65 'polymyxin M' \/  sensitive resistant, ";
	print OUT "66 polymyxins \/  sensitive resistant, ";
	print OUT "67 'procaine benzylpenicillin' \/  sensitive resistant, ";
	print OUT "68 puromycin \/  sensitive resistant, ";
	print OUT "69 rifampicin \/  sensitive resistant, ";
	print OUT "70 rolitetracycline \/  sensitive resistant, ";
	print OUT "71 roxithromycin \/  sensitive resistant, ";
	print OUT "72 'sodium polyanethol sulfonate' \/  sensitive resistant, ";
	print OUT "73 spectinomycin \/  sensitive resistant, ";
	print OUT "74 spiramycin \/  sensitive resistant, ";
	print OUT "75 streptomycin \/  sensitive resistant, ";
	print OUT "76 sulfadiazine \/  sensitive resistant, ";
	print OUT "77 sulfamethoxazole \/  sensitive resistant, ";
	print OUT "78 sulfisomidine \/  sensitive resistant, ";
	print OUT "79 tetracycline \/  sensitive resistant, ";
	print OUT "80 ticarcillin \/  sensitive resistant, ";
	print OUT "81 tobramycin \/  sensitive resistant, ";
	print OUT "82 trimethoprim \/  sensitive resistant, ";
	print OUT "83 vancomycin \/  sensitive resistant, ";
	print OUT "84 'aminoglycoside antibiotics' \/  sensitive resistant, ";
	print OUT "85 'antibiotic dyes' \/  sensitive resistant, ";
	print OUT "86 'beta-lactam antibiotics' \/  sensitive resistant, ";
	print OUT "87 'cephalosporin antibiotics' \/  sensitive resistant, ";
	print OUT "88 'cyclic peptide antibiotics' \/  sensitive resistant, ";
	print OUT "89 'diterpene antibiotics' \/  sensitive resistant, ";
	print OUT "90 'fluoroquinolone antibiotics' \/  sensitive resistant, ";
	print OUT "91 'glycopeptide antibiotics' \/  sensitive resistant, ";
	print OUT "92 'macrolide antibiotics' \/  sensitive resistant, ";
	print OUT "93 'nitrofuran antibiotics' \/  sensitive resistant, ";
	print OUT "94 'nitroimidazole antibiotics' \/  sensitive resistant, ";
	print OUT "95 'nucleoside antibiotics' \/  sensitive resistant, ";
	print OUT "96 'organochlorine antibiotics' \/  sensitive resistant, ";
	print OUT "97 'piperazine antibiotics' \/  sensitive resistant, ";
	print OUT "98 'piperidine antibiotics' \/  sensitive resistant, ";
	print OUT "99 'polyketide antibiotics' \/  sensitive resistant, ";
	print OUT "100 'polypeptide antibiotics' \/  sensitive resistant, ";
	print OUT "101 'pyrimidine antibiotics' \/  sensitive resistant, ";
	print OUT "102 'pyrrolidine antibiotics' \/  sensitive resistant, ";
	print OUT "103 'quinolone antibiotics' \/  sensitive resistant, ";
	print OUT "104 'steroid antibiotics' \/  sensitive resistant, ";
	print OUT "105 'sulfonamide antibiotics' \/  sensitive resistant, ";
	print OUT "106 'tetracycline antibiotics' \/  sensitive resistant, ";

	print OUT "107 acriflavine \/  sensitive resistant, ";
	print OUT "108 arbekacin \/  sensitive resistant, ";
	print OUT "109 arbekacin \/  sensitive resistant, ";
	print OUT "110 azlocillin \/  sensitive resistant, ";
	print OUT "111 aztreonam \/  sensitive resistant, ";
	print OUT "112 bacteriocin \/  sensitive resistant, ";
	print OUT "113 boron \/  sensitive resistant, ";
	print OUT "114 'cadmium acetate' \/  sensitive resistant, ";
	print OUT "115 cefalexin \/  sensitive resistant, ";
	print OUT "116 cefalotin \/  sensitive resistant, ";
	print OUT "117 cefamandole \/  sensitive resistant, ";
	print OUT "118 cefmetazole \/  sensitive resistant, ";
	print OUT "119 cefotiam \/  sensitive resistant, ";
	print OUT "120 ceftiofur \/  sensitive resistant, ";
	print OUT "121 ceftriaxone \/  sensitive resistant, ";
	print OUT "122 'cetyltrimethylammonium bromide' \/  sensitive resistant, ";
	print OUT "123 'chlorhexidine digluconate' \/  sensitive resistant, ";
	print OUT "124 cloxacillin \/  sensitive resistant, ";
	print OUT "125 deferoxamine \/  sensitive resistant, ";
	print OUT "126 efrotomycin \/  sensitive resistant, ";
	print OUT "127 enrofloxacin \/  sensitive resistant, ";
	print OUT "128 florfenicol \/  sensitive resistant, ";
	print OUT "129 flumequine \/  sensitive resistant, ";
	print OUT "130 fosfomycin \/  sensitive resistant, ";
	print OUT "131 isepamicin \/  sensitive resistant, ";
	print OUT "132 josamycin \/  sensitive resistant, ";
	print OUT "133 linezolid \/  sensitive resistant, ";
	print OUT "134 lysostaphin \/  sensitive resistant, ";
	print OUT "135 lysozyme \/  sensitive resistant, ";
	print OUT "136 medicamycin \/  sensitive resistant, ";
	print OUT "137 'mercuric nitrate' \/  sensitive resistant, ";
	print OUT "138 meropenem \/  sensitive resistant, ";
	print OUT "139 methicillin \/  sensitive resistant, ";
	print OUT "140 mupirocin \/  sensitive resistant, ";
	print OUT "141 netilmicin \/  sensitive resistant, ";
	print OUT "142 nitrofurazolidone \/  sensitive resistant, ";
	print OUT "143 nystatin \/  sensitive resistant, ";
	print OUT "144 ofloxacin \/  sensitive resistant, ";
	print OUT "145 optochin \/  sensitive resistant, ";
	print OUT "146 'oxolinic acid' \/  sensitive resistant, ";
	print OUT "147 pefloxacin \/  sensitive resistant, ";
	print OUT "148 phages \/  sensitive resistant, ";
	print OUT "149 'phenylmethylsulfonyl fluoride' \/  sensitive resistant, ";
	print OUT "150 pristinamycin \/  sensitive resistant, ";
	print OUT "151 propamidine \/  sensitive resistant, ";
	print OUT "152 'rifamycin sv' \/  sensitive resistant, ";
	print OUT "153 rifamycin \/  sensitive resistant, ";
	print OUT "154 sulfamides \/  sensitive resistant, ";
	print OUT "155 sulfisoxazole \/  sensitive resistant, ";
	print OUT "156 teicoplanin \/  sensitive resistant, ";
	print OUT "157 tigecycline \/  sensitive resistant, ";
	print OUT "158 troleandomycin \/  sensitive resistant, ";
	print OUT "159 'vibriostatic agent' \/  sensitive resistant, ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @antiresisdata = split (/\t/, $line);
		push (my @taxnames, $antiresisdata[0]);

#code char 1 actinomycin
		if ($line =~ /,actinomycin,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not actinomycin,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 amikacin
		if ($line =~ /,amikacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amikacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 amoxicillin
		if ($line =~ /,amoxicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amoxicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 amphotericin b
		if ($line =~ /,amphotericin b,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphotericin b,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 ampicillin
		if ($line =~ /,ampicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ampicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 anisomycin
		if ($line =~ /,anisomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not anisomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 ansamycin
		if ($line =~ /,ansamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ansamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 aphidicolin
		if ($line =~ /,aphidicolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aphidicolin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 bacitracin
		if ($line =~ /,bacitracin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not bacitracin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 benzylpenicillin
		if ($line =~ /,benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not benzylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 brilliant green
		if ($line =~ /,brilliant green,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not brilliant green,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 carbenicillin
		if ($line =~ /,carbenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carbenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 13 cefadroxil
		if ($line =~ /,cefadroxil,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefadroxil,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 14 cefalotin
		if ($line =~ /,cefalotin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalotin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 cefazolin
		if ($line =~ /,cefazolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefazolin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 cefoperazone
		if ($line =~ /,cefoperazone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefoperazone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 17 cefotaxime
		if ($line =~ /,cefotaxime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefotaxime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 18 cefoxitin
		if ($line =~ /,cefoxitin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefoxitin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 19 cefsulodin
		if ($line =~ /,cefsulodin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefsulodin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 20 ceftazidime
		if ($line =~ /,ceftazidime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftazidime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 21 cefuroxime
		if ($line =~ /,cefuroxime,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefuroxime,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 22 cephaloridine
		if ($line =~ /,cephaloridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephaloridine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 23 cephalosporin
		if ($line =~ /,cephalosporin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephalosporin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 24 chloramphenicol
		if ($line =~ /,chloramphenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chloramphenicol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 25 chlortetracycline
		if ($line =~ /,chlortetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chlortetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 26 ciprofloxacin
		if ($line =~ /,ciprofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ciprofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 27 clavulanic acid
		if ($line =~ /,clavulanic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clavulanic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 28 clindamycin
		if ($line =~ /,clindamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clindamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 29 clomocycline
		if ($line =~ /,clomocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clomocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 30 colistin
		if ($line =~ /,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not colistin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 31 crystal violet
		if ($line =~ /,crystal violet,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not crystal violet,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 32 cycloheximide
		if ($line =~ /,cycloheximide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cycloheximide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 33 demeclocycline
		if ($line =~ /,demeclocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not demeclocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 34 dihydrostreptomycin
		if ($line =~ /,dihydrostreptomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dihydrostreptomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 35 doxycycline
		if ($line =~ /,doxycycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not doxycycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 36 erythromycin
		if ($line =~ /,erythromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erythromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 37 furazolidone
		if ($line =~ /,furazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not furazolidone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 38 fusidic acid
		if ($line =~ /,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fusidic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 39 gentamicin
		if ($line =~ /,gentamicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gentamicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 40 imipenem
		if ($line =~ /,imipenem,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not imipenem,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 41 kanamycin
		if ($line =~ /,kanamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not kanamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 42 levofloxacin
		if ($line =~ /,levofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not levofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 43 lincomycin
		if ($line =~ /,lincomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lincomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 44 lomefloxacin
		if ($line =~ /,lomefloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lomefloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 45 lymecycline
		if ($line =~ /,lymecycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lymecycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 46 meclocycline
		if ($line =~ /,meclocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meclocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 47 metacycline
		if ($line =~ /,metacycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not metacycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 48 metronidazole
		if ($line =~ /,metronidazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not metronidazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 49 xminocyclinexx
		if ($line =~ /,minocycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not minocycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 50 nalidixic acid
		if ($line =~ /,nalidixic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nalidixic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 51 neomycin
		if ($line =~ /,neomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not neomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 52 nitrofurantoin
		if ($line =~ /,nitrofurantoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurantoin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 53 norfloxacin
		if ($line =~ /,norfloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not norfloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 54 novobiocin
		if ($line =~ /,novobiocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not novobiocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 55 oleandomycin
		if ($line =~ /,oleandomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oleandomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 56 omadacycline
		if ($line =~ /,omadacycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not omadacycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 57 oxacillin
		if ($line =~ /,oxacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 58 oxytetracycline
		if ($line =~ /,oxytetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxytetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 59 penicillin
		if ($line =~ /,penicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 60 penicillins
		if ($line =~ /,penicillins,|,phenoxymethylpenicillin,|,benzylpenicillin,|,penicillin,|,procaine benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penicillins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 61 penimepicycline
		if ($line =~ /,penimepicycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not penimepicycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 62 phenoxymethylpenicillin
		if ($line =~ /,phenoxymethylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenoxymethylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 63 piperacillin
		if ($line =~ /,piperacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 64 polymyxin b
		if ($line =~ /,polymyxin b,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxin b,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 65 polymyxin m
		if ($line =~ /,polymyxin m,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxin m,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 66 polymyxins
		if ($line =~ /,polymyxins,|,polymyxin m,|,polymyxin b,|,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polymyxins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 67 procaine benzylpenicillin
		if ($line =~ /,procaine benzylpenicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not procaine benzylpenicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 68 puromycin
		if ($line =~ /,puromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not puromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 69 rifampicin
		if ($line =~ /,rifampicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifampicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 70 rolitetracycline
		if ($line =~ /,rolitetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rolitetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 71 roxithromycin
		if ($line =~ /,roxithromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not roxithromycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 72 sodium polyanethol sulfonate
		if ($line =~ /,sodium polyanethol sulfonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sodium polyanethol sulfonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 73 spectinomycin
		if ($line =~ /,spectinomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spectinomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 74 spiramycin
		if ($line =~ /,spiramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spiramycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 75 streptomycin
		if ($line =~ /,streptomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not streptomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 76 sulfadiazine
		if ($line =~ /,sulfadiazine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfadiazine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 77 sulfamethoxazole
		if ($line =~ /,sulfamethoxazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfamethoxazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 78 sulfisomidine
		if ($line =~ /,sulfisomidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfisomidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 79 tetracycline
		if ($line =~ /,tetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetracycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 80 ticarcillin
		if ($line =~ /,ticarcillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ticarcillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 81 tobramycin
		if ($line =~ /,tobramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tobramycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 82 trimethoprim
		if ($line =~ /,trimethoprim,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trimethoprim,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 83 vancomycin
		if ($line =~ /,vancomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not vancomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 84 aminoglycoside antibiotics
		if ($line =~ /,aminoglycosides,|,arbekacin,|,isepamicin,|,netilmicin,|,spectinomycin,|,dihydrostreptomycin,|,gentamicin,|,kanamycin,|,neomycin,|,streptomycin,|,tobramycin,|,amikacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aminoglycosides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 85 antibiotic dyes
		if ($line =~ /,antibiotic dyes,|,crystal violet,|,brilliant green,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not antibiotic dyes,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 86 beta-lactams
		if ($line =~ /,beta-lactams,|,cefalexin,|,cefalotin,|,cefamandole,|,cefmetazole,|,cefotiam,|,ceftriaxone,|,azlocillin,|,aztreonam,|,cloxacillin,|,methicillin,|,meropenem,|,amoxicillin,|,ampicillin,|,benzylpenicillin,|,carbenicillin,|,cephalosporin,|,clavulanic acid,|,imipenem,|,oxacillin,|,penicillin,|,phenoxymethylpenicillin,|,procaine benzylpenicillin,|,ticarcillin,|,cefoxitin,|,penicillins,|,piperacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not beta-lactams,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 87 cephalosporins
		if ($line =~ /,cephalosporins,|,cefalexin,|,cefalotin,|,cefamandole,|,cefmetazole,|,cefotiam,|,ceftriaxone,|,cefadroxil,|,cefalotin,|,cefazolin,|,cefoperazone,|,cefotaxime,|,cefsulodin,|,ceftazidime,|,cefuroxime,|,cephaloridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cephalosporins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 88 cyclic peptides
		if ($line =~ /,cyclic peptides,|,polymyxin b,|,polymyxin m,|,polymyxins,|,bacitracin,|,colistin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cyclic peptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 89 diterpenes
		if ($line =~ /,diterpenes,|,aphidicolin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diterpenes,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 90 fluoroquinolones
		if ($line =~ /,fluoroquinolones,|,enrofloxacin,|,levofloxacin,|,lomefloxacin,|,norfloxacin,|,ciprofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fluoroquinolones,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 91 glycopeptides
		if ($line =~ /,glycopeptides,|,vancomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycopeptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 92 macrolide antibiotics
		if ($line =~ /,macrolides,|,nystatin,|,josamycin,|,pristinamycin,|,troleandomycin,|,erythromycin,|,oleandomycin,|,roxithromycin,|,spiramycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not macrolides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 93 nitrofuran antibiotics
		if ($line =~ /,nitrofurans,|,furazolidone,|,nitrofurazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurans,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 94 nitroimidazole antibiotics
		if ($line =~ /,nitroimidazole,|,metronidazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitroimidazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 95 nucleoside antibiotics
		if ($line =~ /,nucleosides,|,puromycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nucleosides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 96 organochlorine antibiotics
		if ($line =~ /,organochlorine antibiotics,|,chlorhexidine digluconate,|,clindamycin,|,lincomycin,|,chloramphenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not organochlorine antibiotics,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 97 piperazines
		if ($line =~ /,piperazines,|,rifampicin,|,enrofloxacin,|,ofloxacin,|,pefloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperazines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 98 piperidines
		if ($line =~ /,piperidines,|,cycloheximide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not piperidines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 99 polyketides
		if ($line =~ /,polyketides,|,tigecycline,|,rifamycin,|,rifamycin sv,|,nystatin,|,amphotericin b,|,josamycin,|,pristinamycin,|,troleandomycin,|,ansamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polyketides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 100 polypeptide antibiotics
		if ($line =~ /,polypeptides,|,actinomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polypeptides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 101 pyrimidine antibiotics
		if ($line =~ /,pyrimidiness,|,trimethoprim,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrimidiness,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 102 pyrrolidine antibiotics
		if ($line =~ /,pyrrolidines,|,anisomycin,|,meropenem,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrrolidines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 103 quinolone antibiotics
		if ($line =~ /,quinolones,|,pefloxacin,|,oxolinic acid,|,pefloxacin,|,ofloxacin,|,enrofloxacin,|,nalidixic acid,|,novobiocin,|,nitrofurantoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not quinolones,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 104 steroid antibiotics
		if ($line =~ /,steroid antibiotics,|,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not steroid antibiotics,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 105 sulfonamides
		if ($line =~ /,sulfonamides,|,sulfisoxazole,|,sulfadiazine,|,sulfamethoxazole,|,sulfisomidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfonamides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 106 tetracyclines
		if ($line =~ /,tetracyclines,|,tigecycline,|,chlortetracycline,|,clomocycline,|,demeclocycline,|,doxycycline,|,lymecycline,|,meclocycline,|,metacycline,|,minocycline,|,omadacycline,|,oxytetracycline,|,penimepicycline,|,rolitetracycline,|,tetracycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetracyclines,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

#code char 107 acriflavine
		if ($line =~ /,acriflavine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acriflavine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 108 arbekacin
		if ($line =~ /,arbekacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arbekacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 109 azlocillin
		if ($line =~ /,azlocillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not azlocillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 110 aztreonam
		if ($line =~ /,aztreonam,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aztreonam,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 111 bacteriocin
		if ($line =~ /,bacteriocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not bacteriocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 112 boron
		if ($line =~ /,boron,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not boron,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 113 cadmium acetate
		if ($line =~ /,cadmium acetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cadmium acetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 114 cefalexin
		if ($line =~ /,cefalexin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalexin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 115 cefalotin
		if ($line =~ /,cefalotin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefalotin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 116 cefamandole
		if ($line =~ /,cefamandole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefamandole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 117 cefmetazole
		if ($line =~ /,cefmetazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefmetazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 118 cefotiam
		if ($line =~ /,cefotiam,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cefotiam,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 119 ceftiofur
		if ($line =~ /,ceftiofur,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftiofur,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 120 ceftriaxone
		if ($line =~ /,ceftriaxone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ceftriaxone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 121 cetyltrimethylammonium bromide
		if ($line =~ /,cetyltrimethylammonium bromide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cetyltrimethylammonium bromide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 122 chlorhexidine digluconate
		if ($line =~ /,chlorhexidine digluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chlorhexidine digluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 123 cloxacillin
		if ($line =~ /,cloxacillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cloxacillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 124 deferoxamine
		if ($line =~ /,deferoxamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not deferoxamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 125 efrotomycin
		if ($line =~ /,efrotomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not efrotomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 126 enrofloxacin
		if ($line =~ /,enrofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not enrofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 127 florfenicol
		if ($line =~ /,florfenicol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not florfenicol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 128 flumequine
		if ($line =~ /,flumequine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flumequine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 129 fosfomycin
		if ($line =~ /,fosfomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fosfomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 130 isepamicin
		if ($line =~ /,isepamicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isepamicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 131 josamycin
		if ($line =~ /,josamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not josamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 132 linezolid
		if ($line =~ /,linezolid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not linezolid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 133 lysostaphin
		if ($line =~ /,lysostaphin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysostaphin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 134 lysozyme
		if ($line =~ /,lysozyme,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysozyme,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 135 medicamycin
		if ($line =~ /,medicamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not medicamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 136 mercuric nitrate
		if ($line =~ /,mercuric nitrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mercuric nitrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 137 meropenem
		if ($line =~ /,meropenem,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meropenem,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 138 methicillin
		if ($line =~ /,methicillin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methicillin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 139 mupirocin
		if ($line =~ /,mupirocin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mupirocin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 140 netilmicin
		if ($line =~ /,netilmicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not netilmicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 141 nitrofurazolidone
		if ($line =~ /,nitrofurazolidone,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrofurazolidone,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 142 nystatin
		if ($line =~ /,nystatin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nystatin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 143 ofloxacin
		if ($line =~ /,ofloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ofloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 144 optochin
		if ($line =~ /,optochin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not optochin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 145 oxolinic acid
		if ($line =~ /,oxolinic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxolinic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 146 pefloxacin
		if ($line =~ /,pefloxacin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pefloxacin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 147 phages
		if ($line =~ /,phages,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phages,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 148 phenylmethylsulfonyl fluoride
		if ($line =~ /,phenylmethylsulfonyl fluoride,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylmethylsulfonyl fluoride,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 149 pristinamycin
		if ($line =~ /,pristinamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pristinamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 150 propamidine
		if ($line =~ /,propamidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propamidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 151 rifamycin sv
		if ($line =~ /,rifamycin sv,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifamycin sv,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 152 rifamycin
		if ($line =~ /,rifamycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rifamycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 153 sulfamides
		if ($line =~ /,sulfamides,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfamides,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 154 sulfisoxazole
		if ($line =~ /,sulfisoxazole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sulfisoxazole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 155 teicoplanin
		if ($line =~ /,teicoplanin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not teicoplanin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 156 tigecycline
		if ($line =~ /,tigecycline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tigecycline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 157 troleandomycin
		if ($line =~ /,troleandomycin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not troleandomycin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 158 vibriostatic agent
		if ($line =~ /,vibriostatic agent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not vibriostatic agent,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process colony shape character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "colony shape") {
	my $homout = "hom.colshape.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @colshape, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @colshape;

		map {s/[^[:print:]]//gi; } @colshape; # delete non-printable characters
		map {s/,an irregular,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,becoming rhizoid,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,circular and <1 mm,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular and ~0.5 mm,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular and 1-2 mm,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular and 2.0-3.0 mm,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular and creamy white,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular and low convex,/,circular,convex,/gi; } @colshape; # synonyms of 
		map {s/,circular and pale yellow,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular colony,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular droplets,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular p-agar colonies,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,circular.slightly irregular margins,/,circular,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,circular.slightly irregular,/,circular,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,collapsed center,/,crateriform,/gi; } @colshape; # synonyms of 
		map {s/,colony is round,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,convex and pale orange,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,convex and pale yellow,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,convex colonies,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,convex colony,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,convex elevation,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,crater-shaped,/,crateriform,/gi; } @colshape; # synonyms of 
		map {s/,crateriform colony,/,crateriform,/gi; } @colshape; # synonyms of 
		map {s/,cup-shaped,/,crateriform,/gi; } @colshape; # synonyms of 
		map {s/,curled colony,/,curled,/gi; } @colshape; # synonyms of 
		map {s/,dendritic colony,/,dendritic,/gi; } @colshape; # synonyms of 
		map {s/,dimpled,/,crateriform,/gi; } @colshape; # synonyms of 
		map {s/,disc.shaped,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,disk.shaped,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,domed,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,dry colony,/,dry,/gi; } @colshape; # synonyms of 
		map {s/,elevated and yellowish orange-coloured,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,elevated center,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,elevated,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,entire colony,/,entire,/gi; } @colshape; # synonyms of 
		map {s/,erose colony,/,erose,/gi; } @colshape; # synonyms of 
		map {s/,fine granulation,/,granular,/gi; } @colshape; # synonyms of 
		map {s/,flat and 1.1-3.8 mm,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,flat colonies,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,flat colony,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,flat edge,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,flat elevation,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,flat shape,/,flat,/gi; } @colshape; # synonyms of 
		map {s/,forming round,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,forms circular colonies,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,friable,/,dry,/gi; } @colshape; # synonyms of 
		map {s/,fried.egg,/,fried-egg,/gi; } @colshape; # synonyms of 
		map {s/,fried.eggs,/,fried-egg,/gi; } @colshape; # synonyms of 
		map {s/,fusiform,/,spindle,/gi; } @colshape; # synonyms of 
		map {s/,fuzzy,/,filamentous,/gi; } @colshape; # synonyms of 
		map {s/,granular center,/,granular,/gi; } @colshape; # synonyms of 
		map {s/,hemispherical knobs,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,high convex,/,raised,convex,/gi; } @colshape; # synonyms of 
		map {s/,irregular borders,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular centers,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular centres,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular edges,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular form,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular margins,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular shape,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,irregular,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,lobed margin,/,lobed,/gi; } @colshape; # synonyms of 
		map {s/,low convex,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,minute colonies,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,minute,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,more rhizoid,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,mucous,/,mucoid,/gi; } @colshape; # synonyms of 
		map {s/,notched,/,erose,/gi; } @colshape; # synonyms of 
		map {s/,peak,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,peaked,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,pin point,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,pulvinate colony,/,pulvinate,/gi; } @colshape; # synonyms of 
		map {s/,pulvinate elevation,/,pulvinate,/gi; } @colshape; # synonyms of 
		map {s/,punctate,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,punctiform colony,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,punctiforme,/,punctiform,/gi; } @colshape; # synonyms of 
		map {s/,raised at the center,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,raised at the centres,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,raised colony,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,raised in the center,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,raised in the centre,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,raised rhizoid growth,/,rhizoidal,raised,/gi; } @colshape; # synonyms of 
		map {s/,raised,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,rhizoid form,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,rhizoid form,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,rhizoid,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,round colonies,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,round,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,rounded,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,satellite behavior,/,satellite,/gi; } @colshape; # synonyms of 
		map {s/,satellite colonies,/,satellite,/gi; } @colshape; # synonyms of 
		map {s/,satellites,/,satellite,/gi; } @colshape; # synonyms of 
		map {s/,slightly convex and moderate yellow,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,slightly convex,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,slightly elevated center,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,slightly elevated centers,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,slightly elevated,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,slightly filamentous margins,/,filamentous,/gi; } @colshape; # synonyms of 
		map {s/,slightly irregular margins,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,slightly irregular,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,slightly peaked,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,slightly raised,/,raised,/gi; } @colshape; # synonyms of 
		map {s/,slightly sunken,/,sunken,/gi; } @colshape; # synonyms of 
		map {s/,slightly umbonate,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,slightly umobonate,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,sometimes feebly rhizoid,/,rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,spindle-shaped,/,spindle,/gi; } @colshape; # synonyms of 
		map {s/,spreading,/,effuse,/gi; } @colshape; # synonyms of 
		map {s/,thin,/,diffuse,/gi; } @colshape; # synonyms of 
		map {s/,toothed,/,erose,/gi; } @colshape; # synonyms of 
		map {s/,translucent and low convex,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,transparent and low convex,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,typically spindle-shaped,/,spindle,/gi; } @colshape; # synonyms of 
		map {s/,umbonate elevation,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,umobonate,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,undulate margin,/,undulate,/gi; } @colshape; # synonyms of 
		map {s/,uniformly round,/,circular,/gi; } @colshape; # synonyms of 
		map {s/,usually irregular,/,asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,veined,/,dendritic,/gi; } @colshape; # synonyms of 
		map {s/,very slightly umbonate,/,umbonate,/gi; } @colshape; # synonyms of 
		map {s/,wavy margin,/,curled,/gi; } @colshape; # synonyms of 
		map {s/,whitish and low convex,/,convex,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @colshape; # synonyms of 

		map {s/,not circular colony,/,not circular,/gi; } @colshape; # synonyms of 
		map {s/,not circular droplets,/,not circular,/gi; } @colshape; # synonyms of 
		map {s/,not collapsed center,/,not crateriform,/gi; } @colshape; # synonyms of 
		map {s/,not colony is round,/,not circular,/gi; } @colshape; # synonyms of 
		map {s/,not convex colony,/,not convex,/gi; } @colshape; # synonyms of 
		map {s/,not convex elevation,/,not convex,/gi; } @colshape; # synonyms of 
		map {s/,not crater-shaped,/,not crateriform,/gi; } @colshape; # synonyms of 
		map {s/,not crateriform colony,/,not crateriform,/gi; } @colshape; # synonyms of 
		map {s/,not cup-shaped,/,not crateriform,/gi; } @colshape; # synonyms of 
		map {s/,not curled colony,/,not curled,/gi; } @colshape; # synonyms of 
		map {s/,not dendritic colony,/,not dendritic,/gi; } @colshape; # synonyms of 
		map {s/,not dimpled,/,not crateriform,/gi; } @colshape; # synonyms of 
		map {s/,not disc.shaped,/,not flat,/gi; } @colshape; # synonyms of 
		map {s/,not disk.shaped,/,not flat,/gi; } @colshape; # synonyms of 
		map {s/,not domed,/,not convex,/gi; } @colshape; # synonyms of 
		map {s/,not dry colony,/,not dry,/gi; } @colshape; # synonyms of 
		map {s/,not elevated,/,not raised,/gi; } @colshape; # synonyms of 
		map {s/,not entire colony,/,not entire,/gi; } @colshape; # synonyms of 
		map {s/,not erose colony,/,not erose,/gi; } @colshape; # synonyms of 
		map {s/,not fine granulation,/,not granular,/gi; } @colshape; # synonyms of 
		map {s/,not flat colonies,/,not flat,/gi; } @colshape; # synonyms of 
		map {s/,not flat colony,/,not flat,/gi; } @colshape; # synonyms of 
		map {s/,not flat elevation,/,not flat,/gi; } @colshape; # synonyms of 
		map {s/,not friable,/,not dry,/gi; } @colshape; # synonyms of 
		map {s/,not fried.egg,/,not fried-egg,/gi; } @colshape; # synonyms of 
		map {s/,not fried.eggs,/,not fried-egg,/gi; } @colshape; # synonyms of 
		map {s/,not fusiform,/,not spindle,/gi; } @colshape; # synonyms of 
		map {s/,not fuzzy,/,not filamentous,/gi; } @colshape; # synonyms of 
		map {s/,not granular center,/,not granular,/gi; } @colshape; # synonyms of 
		map {s/,not hemispherical knobs,/,not convex,/gi; } @colshape; # synonyms of 
		map {s/,not irregular,/,not asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,not  irregular,/,not asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,not irregular borders,/,not asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,not  irregular borders,/,not asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,not lobed margin,/,not lobed,/gi; } @colshape; # synonyms of 
		map {s/,not low convex,/,not flat,convex,/gi; } @colshape; # synonyms of 
		map {s/,not more rhizoid,/,not rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,not mucous,/,not mucoid,/gi; } @colshape; # synonyms of 
		map {s/,not notched,/,not erose,/gi; } @colshape; # synonyms of 
		map {s/,not peak,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not peaked,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not pin point,/,not punctiform,/gi; } @colshape; # synonyms of 
		map {s/,not pulvinate colony,/,not pulvinate,/gi; } @colshape; # synonyms of 
		map {s/,not pulvinate elevation,/,not pulvinate,/gi; } @colshape; # synonyms of 
		map {s/,not punctate,/,not punctiform,/gi; } @colshape; # synonyms of 
		map {s/,not punctiform colony,/,not punctiform,/gi; } @colshape; # synonyms of 
		map {s/,not punctiforme,/,not punctiform,/gi; } @colshape; # synonyms of 
		map {s/,not raised at the center,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not raised at the centres,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not raised colony,/,not raised,/gi; } @colshape; # synonyms of 
		map {s/,not raised in the center,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not raised in the centre,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not raised,/,not raised,/gi; } @colshape; # synonyms of 
		map {s/,not rhizoid form,/,not rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,not rhizoid form,/,not rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,not rhizoid,/,not rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,not round,/,not circular,/gi; } @colshape; # synonyms of 
		map {s/,not rounded,/,not circular,/gi; } @colshape; # synonyms of 
		map {s/,not satellite behavior,/,not satellite,/gi; } @colshape; # synonyms of 
		map {s/,not satellite colonies,/,not satellite,/gi; } @colshape; # synonyms of 
		map {s/,not satellites,/,not satellite,/gi; } @colshape; # synonyms of 
		map {s/,not slightly convex,/,not convex,/gi; } @colshape; # synonyms of 
		map {s/,not slightly irregular,/,not asymmetrical,/gi; } @colshape; # synonyms of 
		map {s/,not slightly peaked,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not slightly raised,/,not raised,/gi; } @colshape; # synonyms of 
		map {s/,not slightly sunken,/,not sunken,/gi; } @colshape; # synonyms of 
		map {s/,not slightly umbonate,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not slightly umobonate,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not sometimes feebly rhizoid,/,not rhizoidal,/gi; } @colshape; # synonyms of 
		map {s/,not spindle-shaped,/,not spindle,/gi; } @colshape; # synonyms of 
		map {s/,not spreading,/,not effuse,/gi; } @colshape; # synonyms of 
		map {s/,not thin,/,not diffuse,/gi; } @colshape; # synonyms of 
		map {s/,not toothed,/,not erose,/gi; } @colshape; # synonyms of 
		map {s/,not umbonate elevation,/,not umbonate,/gi; } @colshape; # synonyms of 
		map {s/,not undulate margin,/,not undulate,/gi; } @colshape; # synonyms of 
		map {s/,not veined,/,not dendritic,/gi; } @colshape; # synonyms of 
		map {s/,not wavy margin,/,not curled,/gi; } @colshape; # synonyms of 

#not colony shape
		map {s/,numerous filamentous tufts,/,/gi; } @colshape; # synonyms of 
		map {s/,[\d]+ g,/,/gi; } @colshape; # synonyms of 
		map {s/,cefotaxime,/,/gi; } @colshape; # synonyms of 
		map {s/,ceftazidime,/,/gi; } @colshape; # synonyms of 
		map {s/,ceftriaxone,/,/gi; } @colshape; # synonyms of 
		map {s/,elevated temperature,/,/gi; } @colshape; # synonyms of 
		map {s/,fusidic acid,/,/gi; } @colshape; # synonyms of 
		map {s/,g,/,/gi; } @colshape; # synonyms of 
		map {s/,nitrofurantoin,/,/gi; } @colshape; # synonyms of 
		map {s/,pefloxacin,/,/gi; } @colshape; # synonyms of 
		map {s/,tetracycline,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 
		map {s/,xxx,/,/gi; } @colshape; # synonyms of 

		print OUT @taxlabels, @colshape, "\n"; # prints to $homout, hom.colshape.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.colshape.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,colony shape,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.colshape.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.colshape.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.colshape.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.colshape.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.colshape.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^asymmetrical/**asymmetrical/g; } @charlist2;
	map {s/^circular/**circular/g; } @charlist2;
	map {s/^convex/**convex/g; } @charlist2;
	map {s/^crateriform/**crateriform/g; } @charlist2;
	map {s/^curled/**curled/g; } @charlist2;
	map {s/^dendritic/**dendritic/g; } @charlist2;
	map {s/^dense center/**dense center/g; } @charlist2;
	map {s/^diffuse/**diffuse/g; } @charlist2;
	map {s/^dry/**dry/g; } @charlist2;
	map {s/^effuse/**effuse/g; } @charlist2;
	map {s/^entire/**entire/g; } @charlist2;
	map {s/^erose/**erose/g; } @charlist2;
	map {s/^filamentous/**filamentous/g; } @charlist2;
	map {s/^flat/**flat/g; } @charlist2;
	map {s/^fried-egg/**fried-egg/g; } @charlist2;
	map {s/^lobed/**lobed/g; } @charlist2;
	map {s/^pulvinate/**pulvinate/g; } @charlist2;
	map {s/^punctiform/**punctiform/g; } @charlist2;
	map {s/^pyramidal/**pyramidal/g; } @charlist2;
	map {s/^raised/**raised/g; } @charlist2;
	map {s/^rhizoidal/**rhizoidal/g; } @charlist2;
	map {s/^satellite/**satellite/g; } @charlist2;
	map {s/^spindle/**spindle/g; } @charlist2;
	map {s/^sunken/**sunken/g; } @charlist2;
	map {s/^umbonate/**umbonate/g; } @charlist2;
	map {s/^undulate/**undulate/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.colshape.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.colshape.txt
		}
	my $temp6 = "temp6.colshape.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.colshape.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "colshape.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @colshapedata = split (/\t/, $line);
		push (@taxnames, $colshapedata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Shape Matrix'\;\n\tDIMENSIONS NCHAR=26\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'circular colony' \/  'colonies not circular' 'colonies circular', ";
	print OUT "2 'entire colony' \/  'colonies not entire' 'colonies entire', ";
	print OUT "3 'convex colony' \/  'colonies not convex' 'colonies convex', ";
	print OUT "4 'crateriform colony' \/  'colonies not crateriform' 'colonies crateriform', ";
	print OUT "5 'curled colony' \/  'colonies not curled' 'colonies curled', ";
	print OUT "6 'dendritic colony' \/  'colonies not dendritic' 'colonies dendritic', ";
	print OUT "7 'erose colony' \/  'colonies not erose' 'colonies erose', ";
	print OUT "8 'filamentous colony' \/  'colonies not filamentous' 'colonies filamentous', ";
	print OUT "9 'flat colony' \/  'colonies not flat' 'colonies flat', ";
	print OUT "10 'lobed colony' \/  'colonies not lobed' 'colonies lobed', ";
	print OUT "11 'pulvinate colony' \/  'colonies not pulvinate' 'colonies pulvinate', ";
	print OUT "12 'punctiform colony' \/  'colonies not punctiform' 'colonies punctiform', ";
	print OUT "13 'raised colony' \/  'colonies not raised' 'colonies raised', ";
	print OUT "14 'spindle-shaped colony' \/  'colonies not spindle-shaped' 'colonies spindle-shaped', ";
	print OUT "15 'undulate colony' \/  'colonies not undulate' 'colonies undulate', ";
	print OUT "16 'satellite colonies' \/  'no satellite colonies' 'satellite colonies', ";
	print OUT "17 'sunken colony' \/  'colonies not sunken' 'colonies sunken', ";
	print OUT "18 'asymmetrical colony' \/  'colonies symmetrical' 'colonies asymmetrical', ";
	print OUT "19 'pyramidal colony' \/  'colonies not pyramidal' 'colonies pyramidal', ";
	print OUT "20 'umbonate colony' \/  'colonies not umbonate' 'colonies umbonate', ";
	print OUT "21 'rhizoidal colony' \/  'colonies not rhizoidal' 'colonies rhizoidal', ";
	print OUT "22 'fried-egg colony' \/  'colonies not fried-egg' 'colonies fried-egg', ";
	print OUT "23 'dry colony' \/  'colonies not dry' 'colonies dry', ";
	print OUT "24 'diffuse colony' \/  'colonies not diffuse' 'colonies diffuse', ";
	print OUT "25 'effuse colony' \/  'colonies not effuse' 'colonies effuse', ";
	print OUT "26 'dense center colony' \/  'colonies lack dense center ' 'colonies have a dense center', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @colshapedata = split (/\t/, $line);
		push (my @taxnames, $colshapedata[0]);

#code char 1 circular
		if ($line =~ /,circular,|,entire,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,circular,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 entire
		if ($line =~ /,entire,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not entire,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 convex
		if ($line =~ /,convex,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not convex,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 crateriform
		if ($line =~ /,crateriform,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not crateriform,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 curled
		if ($line =~ /,curled,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not curled,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 dendritic
		if ($line =~ /,dendritic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dendritic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 erose
		if ($line =~ /,erose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 filamentous
		if ($line =~ /,filamentous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not filamentous,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 flat
		if ($line =~ /,flat,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flat,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 lobed
		if ($line =~ /,lobed,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lobed,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 pulvinate
		if ($line =~ /,pulvinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pulvinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 punctiform
		if ($line =~ /,punctiform,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not punctiform,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 13 raised
		if ($line =~ /,raised,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not raised,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 14 spindle
		if ($line =~ /,spindle,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not spindle,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 undulate
		if ($line =~ /,undulate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not undulate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 satellite
		if ($line =~ /,satellite,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not satellite,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 17 sunken
		if ($line =~ /,sunken,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sunken,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 18 asymmetrical
		if ($line =~ /,asymmetrical,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not asymmetrical,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 19 pyramidal
		if ($line =~ /,pyramidal,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyramidal,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 20 umbonate
		if ($line =~ /,umbonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not umbonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 21 rhizoidal
		if ($line =~ /,rhizoidal,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rhizoidal,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 22 fried-egg
		if ($line =~ /,fried-egg,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fried-egg,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 23 dry
		if ($line =~ /,dry,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dry,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 24 diffuse
		if ($line =~ /,diffuse,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diffuse,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 25 effuse
		if ($line =~ /,effuse,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not effuse,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 26 dense center
		if ($line =~ /,dense center,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dense center,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

	print OUT "\n";

	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process colony margin character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "colony margin") {
	my $homout = "hom.colmargin.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @colmargin, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @colmargin;
		map {s/,a-haemolysis zone,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-haemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-haemolytic zone,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-haemolytic,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-hemolysis zone,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-hemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,a-hemolytic zone,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,almost entire margins,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,alpha-haemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-haemolytic,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-hemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-hemolytic,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-like haemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-like haemolytic,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-like hemolysis,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,alpha-like hemolytic,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-haemolysis zone,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-haemolysis zone,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-haemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-haemolystic zone,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-haemolytic,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-hemolysis zone,/,a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-hemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,b-hemolystic zone,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-haemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-haemolytic,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-hemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-hemolytic,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-like haemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-like haemolytic,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,beta-like hemolysis,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,haemolysis,/,hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,haemolytic,/,hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,hemolysis,/,hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,diffuse,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,effuse,/,effuse margins,/g; } @colmargin; # synonyms of 
		map {s/,entire,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,entire edge,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,entire edges,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,entire margin,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,entire translucent margins,/,entire margins,clear margins,/g; } @colmargin; # synonyms of 
		map {s/,entirely margined,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,erose,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,flame-like edge,/,rhizoidal margins,/g; } @colmargin; # synonyms of 
		map {s/,irregular,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,irregular margins,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,lobate edge,/,lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,margins irregular,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,non-haemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,non-hemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,nonhemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,nonhaemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,regular edges,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,round margins,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,scalloped margins,/,lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,slightly erose,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,slightly spreading edges,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,spreading,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,spreading edge,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,spreading edge possessing lighter pigmentation,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,spreading margin,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,undulate,/,undulate margins,/g; } @colmargin; # synonyms of 
		map {s/,uneven edges,/,undulate margins,/g; } @colmargin; # synonyms of 
		map {s/,weakly b-haemolytic,/,b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,an entire,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,beveled edge,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,circular,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,circular.slightly irregular margins,/,entire margins,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,continuous margins,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,crenate edges,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,diffuse edges,/,diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,entire and 1-2 mm,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,erose margin,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,filamentous,/,rhizoidal margins,/g; } @colmargin; # synonyms of 
		map {s/,flat edge,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,fuzzy edges,/,rhizoidal margins,/g; } @colmargin; # synonyms of 
		map {s/,indented margins,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,irregular edges,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,lobate,/,lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,regular margins,/,symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,scalloped edges,/,lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,slightly filamentous margins,/,rhizoidal margins,/g; } @colmargin; # synonyms of 
		map {s/,slightly irregular edge,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,slightly irregular margins,/,not symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,small regular margins,/,symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,smooth edges,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,smooth margins,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,undulating margins,/,undulate margins,/g; } @colmargin; # synonyms of 
		map {s/,usually entire,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,very pronounced crenate edges,/,erose margins,/g; } @colmargin; # synonyms of 
		map {s/,whole margins,/,entire margins,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colmargin; # synonyms of 

		map {s/,not a-haemolysis zone,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not a-haemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not a-haemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not a-hemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not a-hemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not almost entire margins,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-haemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-hemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-haemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-hemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-like haemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-like haemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-like hemolytic,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not alpha-like hemolysis,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-haemolysis zone,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-haemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-haemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-hemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-hemolysis zone,/,not a-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-hemolysis zone,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not b-hemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-haemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-haemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-hemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-hemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-like haemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-like haemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-like hemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not beta-like hemolysis,/,not b-hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not haemolysis,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not haemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not hemolytic,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not hemolysis,/,not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not diffuse,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not effuse,/,not effuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not entire,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not entire edge,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not entire edges,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not entire margin,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not entire translucent margins,/,not entire margins,clear margins,/g; } @colmargin; # synonyms of 
		map {s/,not entirely margined,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not erose,/,not erose margins,/g; } @colmargin; # synonyms of 
		map {s/,not flame-like edge,/,not rhizoidal margins,/g; } @colmargin; # synonyms of 
		map {s/,not irregular,/,symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,not irregular margins,/,symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,not lobate edge,/,not lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,not margins irregular,/,symmetric margins,/g; } @colmargin; # synonyms of 
		map {s/,not non-haemolytic,/,not not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not non-hemolytic,/,not not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not nonhemolytic,/,not not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not nonhaemolytic,/,not not hemolytic,/g; } @colmargin; # synonyms of 
		map {s/,not regular edges,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not round margins,/,not entire margins,/g; } @colmargin; # synonyms of 
		map {s/,not scalloped margins,/,not lobed margins,/g; } @colmargin; # synonyms of 
		map {s/,not slightly erose,/,not erose margins,/g; } @colmargin; # synonyms of 
		map {s/,not slightly spreading edges,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not spreading,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not spreading edge,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not spreading edge possessing lighter pigmentation,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not spreading margin,/,not diffuse margins,/g; } @colmargin; # synonyms of 
		map {s/,not undulate,/,not undulate margins,/g; } @colmargin; # synonyms of 
		map {s/,not uneven edges,/,not undulate margins,/g; } @colmargin; # synonyms of 
		map {s/,not weakly b-haemolytic,/,not b-hemolytic,/g; } @colmargin; # synonyms of 

#not colony margin
		map {s/,margins,/,/g; } @colmargin; # synonyms of 
		map {s/,numerous filamentous tufts,/,/g; } @colmargin; # synonyms of 
		map {s/,colony edges,/,/g; } @colmargin; # synonyms of 
		map {s/,colony margins,/,/g; } @colmargin; # synonyms of 
		map {s/,cream edge,/,/g; } @colmargin; # synonyms of 
		map {s/,edge,/,/g; } @colmargin; # synonyms of 
		map {s/,edges,/,/g; } @colmargin; # synonyms of 
		map {s/,entire cell,/,/g; } @colmargin; # synonyms of 
		map {s/,entire internal structure,/,/g; } @colmargin; # synonyms of 
		map {s/,entire surface,/,/g; } @colmargin; # synonyms of 
		map {s/,entire t,/,/g; } @colmargin; # synonyms of 
		map {s/,low-orange edge,/,/g; } @colmargin; # synonyms of 
		map {s/,margin,/,/g; } @colmargin; # synonyms of 
		map {s/,translucent edges,/,/g; } @colmargin; # synonyms of 
		map {s/,viscous,/,/g; } @colmargin; # synonyms of 
		map {s/,viscous texture,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 
		map {s/,xxx,/,/g; } @colmargin; # synonyms of 

		print OUT @taxlabels, @colmargin, "\n"; # prints to $homout, hom.colmargin.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.colmargin.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,colony margin,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.colmargin.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.colmargin.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.colmargin.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.colmargin.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.colmargin.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^a-hemolytic/**a-hemolytic/g; } @charlist2;
	map {s/^b-hemolytic/**b-hemolytic/g; } @charlist2;
	map {s/^clear margins/**clear margins/g; } @charlist2;
	map {s/^diffuse margins/**diffuse margins/g; } @charlist2;
	map {s/^effuse margins/**effuse margins/g; } @charlist2;
	map {s/^entire margins/**entire margins/g; } @charlist2;
	map {s/^erose margins/**erose margins/g; } @charlist2;
	map {s/^hemolytic/**hemolytic/g; } @charlist2;
	map {s/^lobed margins/**lobed margins/g; } @charlist2;
	map {s/^rhizoidal margins/**rhizoidal margins/g; } @charlist2;
	map {s/^symmetric margin/**symmetric margin/g; } @charlist2;
	map {s/^undulate margins/**undulate margins/g; } @charlist2;
	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.colmargin.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.colmargin.txt
		}
	my $temp6 = "temp6.colmargin.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.colmargin.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "colmargin.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @colmargindata = split (/\t/, $line);
		push (@taxnames, $colmargindata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Margin Matrix'\;\n\tDIMENSIONS NCHAR=12\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'hemolytic colony' \/  'not hemolytic' hemolytic, ";
	print OUT "2 'alpha hemolytic colony' \/  'not alpha hemoltyic' 'alpha hemolytic', ";
	print OUT "3 'beta hemolytic colony' \/  'not beta hemolytic' 'beta hemolytic', ";
	print OUT "4 'symmetric colony margins' \/  'asymmetric margins' 'symmetric margins', ";
	print OUT "5 'clear margins' \/  'not clear margins' 'clear margins', ";
	print OUT "6 'entire margins' \/  'not entire margins' 'entire margins', ";
	print OUT "7 'diffuse margins' \/  'not diffuse margins' 'diffuse margins', ";
	print OUT "8 'effuse margins' \/  'not effuse margins' 'effuse margins', ";
	print OUT "9 'erose margins' \/  'not erose margins' 'erose margins', ";
	print OUT "10 'lobed margins' \/  'not lobed margins' 'lobed margins', ";
	print OUT "11 'rhizoidal margins' \/  'not rhizoidal margins' 'rhizoidal margins', ";
	print OUT "12 'undulate margins' \/  'not undulate margins' 'undulate margins', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @colmargindata = split (/\t/, $line);
		push (my @taxnames, $colmargindata[0]);

#code char 1 hemolytic
		if ($line =~ /,hemolytic,|,a-hemolytic,|,b-hemolytic,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not hemolytic,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 a-hemolytic
		if ($line =~ /,a-hemolytic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-hemolytic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 b-hemolytic
		if ($line =~ /,b-hemolytic,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-hemolytic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 symmetric colony margin
		if ($line =~ /,symmetric margin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not symmetric margin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 clear margins
		if ($line =~ /,clear margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not clear margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 entire margins
		if ($line =~ /,entire margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not entire margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 diffuse margins
		if ($line =~ /,diffuse margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not diffuse margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 effuse margins
		if ($line =~ /,effuse margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not effuse margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 erose margins
		if ($line =~ /,erose margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erose margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 lobed margins
		if ($line =~ /,lobed margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lobed margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 rhizoidal margins
		if ($line =~ /,rhizoidal margin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rhizoidal margin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 undulate margins
		if ($line =~ /,undulate margins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not undulate margins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process colony texture character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "colony texture") {
	my $homout = "hom.coltext.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @coltext, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @coltext;
		map {s/,a smooth,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,a stringy,/,ropy texture,/g; } @coltext; # synonyms of 
		map {s/,aromatic,/,aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,butryaceous,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,buttery texture,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,butyrous consistency,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,butyrous,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,concentric ridges,/,ridged texture,/g; } @coltext; # synonyms of 
		map {s/,concentric,/,ridged texture,/g; } @coltext; # synonyms of 
		map {s/,do not adhere to the agar,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,dry surfaces,/,dry texture,/g; } @coltext; # synonyms of 
		map {s/,elastic-gummy consistency,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,finely granular surfaces,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,flaky,/,flaky texture,/g; } @coltext; # synonyms of 
		map {s/,glistening,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,granular appearance,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,granular textures,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,granular,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,grayish concentric ring pattern,/,ridged texture,/g; } @coltext; # synonyms of 
		map {s/,growth tends to be mainly below the surface,/,sunken texture,/g; } @coltext; # synonyms of 
		map {s/,hard to emulsify,/,not easily emulsified,/g; } @coltext; # synonyms of 
		map {s/,moist,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,mucoid,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,non-adherent,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,occasionally granular,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,occasionally smooth,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,opalescence,/,opalescent texture,/g; } @coltext; # synonyms of 
		map {s/,opaque and granular surface texture,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,polished,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,relatively smooth,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,rough matt surface,/,dull texture,/g; } @coltext; # synonyms of 
		map {s/,rough surface,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,rough surfaces,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,rough,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,shiny,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,slightly aromatic odour,/,aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,slightly rough and matt surfaces,/,rough texture,dull texture,/g; } @coltext; # synonyms of 
		map {s/,slightly rough,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,slime,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,slimy,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,small shallow pits,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,smooth and 1-2 mm diameter,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth and glistening surfaces,/,smooth texture,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,smooth and pale yellow,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth colonies,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth dull surface,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth dull,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth glossy surfaces,/,smooth texture,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,smooth surface,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth surfaces,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth texture,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth-surfaced,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,soft,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,sometimes stringy,/,ropy texture,/g; } @coltext; # synonyms of 
		map {s/,sticky consistency,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,sticky,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,stringy,/,ropy texture,/g; } @coltext; # synonyms of 
		map {s/,strongly adherent,/,adherent texture,/g; } @coltext; # synonyms of 
		map {s/,sunken,/,sunken texture,/g; } @coltext; # synonyms of 
		map {s/,tough,/,adherent texture,/g; } @coltext; # synonyms of 
		map {s/,translucent colonies,/,translucent colony,/g; } @coltext; # synonyms of 
		map {s/,translucent,/,translucent colony,/g; } @coltext; # synonyms of 
		map {s/,usually sticky,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,viscid consistency,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,viscid,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,watery,/,watery colony,/g; } @coltext; # synonyms of 
		map {s/,white granular deposit,/,granular texture,/g; } @coltext; # synonyms of 
		map {s/,wrinkled,/,wrinkled texture,/g; } @coltext; # synonyms of 
		map {s/,smooth glistening surfaces,/,smooth texture,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @coltext; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @coltext; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @coltext; # synonyms of 

		map {s/,easy to emulsify,/,easily emulsified,/g; } @coltext; # synonyms of 
		map {s/,no  a smooth,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not a smooth,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not a stringy,/,not ropy texture,/g; } @coltext; # synonyms of 
		map {s/,not aromatic,/,not aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,not butryaceous,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not buttery texture,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not butyrous consistency,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not butyrous,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not do not adhere to the agar,/,not not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not elastic-gummy consistency,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not flaky,/,not flaky texture,/g; } @coltext; # synonyms of 
		map {s/,not glistening,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not growth tends to be mainly below the surface,/,not sunken texture,/g; } @coltext; # synonyms of 
		map {s/,not mucoid,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not non-adherent,/,not not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not occasionally smooth,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not opalescence,/,not opalescent texture,/g; } @coltext; # synonyms of 
		map {s/,not polished,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not rough matt surface,/,not dull texture,/g; } @coltext; # synonyms of 
		map {s/,not rough surface,/,not rough texture,/g; } @coltext; # synonyms of 
		map {s/,not shiny,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not slightly aromatic odour,/,not aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,not slime,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not slimy,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not small shallow pits,/,not rough texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth surface,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth surfaces,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth-surfaced,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not soft,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not sometimes stringy,/,ropy texture,/g; } @coltext; # synonyms of 
		map {s/,not stringy,/,not ropy texture,/g; } @coltext; # synonyms of 
		map {s/,not sunken,/,not sunken texture,/g; } @coltext; # synonyms of 
		map {s/,not tough,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not translucent colonies,/,not translucent colony,/g; } @coltext; # synonyms of 
		map {s/,not translucent,/,not translucent colony,/g; } @coltext; # synonyms of 
		map {s/,not viscid consistency,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not viscid,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not watery,/,not watery colony,/g; } @coltext; # synonyms of 

		map {s/,filamentous,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,fuzzy edges,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,slightly filamentous margins,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,slightly filamentous,/,not smooth texture,/g; } @coltext; # synonyms of 

#not colony texture
		map {s/,bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,compact centre,/,/g; } @coltext; # synonyms of 
		map {s/,contoured,/,/g; } @coltext; # synonyms of 
		map {s/,dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,entire,/,/g; } @coltext; # synonyms of 
		map {s/,hemolytic activity,/,/g; } @coltext; # synonyms of 
		map {s/,hemolytic effect,/,/g; } @coltext; # synonyms of 
		map {s/,hemolytic effects,/,/g; } @coltext; # synonyms of 
		map {s/,hemolytic,/,/g; } @coltext; # synonyms of 
		map {s/,incompletely haemolytic,/,/g; } @coltext; # synonyms of 
		map {s/,irregular,/,/g; } @coltext; # synonyms of 
		map {s/,non-haemolytic and 1-3 mm,/,/g; } @coltext; # synonyms of 
		map {s/,non-haemolytic,/,/g; } @coltext; # synonyms of 
		map {s/,non-hemolytic,/,/g; } @coltext; # synonyms of 
		map {s/,not bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not yellow colony,/,/g; } @coltext; # synonyms of 
		map {s/,not yellow,/,/g; } @coltext; # synonyms of 
		map {s/,numerous filamentous tufts,/,/g; } @coltext; # synonyms of 
		map {s/,satellite colonies,/,/g; } @coltext; # synonyms of 
		map {s/,slightly swollen sporangia,/,/g; } @coltext; # synonyms of 
		map {s/,slightly swollen,/,/g; } @coltext; # synonyms of 
		map {s/,smooth edges,/,/g; } @coltext; # synonyms of 
		map {s/,smooth sediment,/,/g; } @coltext; # synonyms of 
		map {s/,stringy sediment,/,/g; } @coltext; # synonyms of 
		map {s/,subterminal swollen sporangia,/,/g; } @coltext; # synonyms of 
		map {s/,swollen sporangia,/,/g; } @coltext; # synonyms of 
		map {s/,swollen sporangium,/,/g; } @coltext; # synonyms of 
		map {s/,swollen,/,/g; } @coltext; # synonyms of 
		map {s/,yellow colony,/,/g; } @coltext; # synonyms of 
		map {s/,yellow,/,/g; } @coltext; # synonyms of 
		map {s/,yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not swollen,/,/g; } @coltext; # synonyms of 
		map {s/,not  swollen,/,/g; } @coltext; # synonyms of 
		map {s/,xxx,/,/g; } @coltext; # synonyms of 

		print OUT @taxlabels, @coltext, "\n"; # prints to $homout, hom.coltext.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.coltext.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,colony texture,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.coltext.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.coltext.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.coltext.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.coltext.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.coltext.txt

#temporarily rename charstates to label those that have been homologized		
	map {s/^adherent texture/**adherent texture/g; } @charlist2;
	map {s/^aromatic colony/**aromatic colony/g; } @charlist2;
	map {s/^butyrous texture/**butyrous texture/g; } @charlist2;
	map {s/^dull texture/**dull texture/g; } @charlist2;
	map {s/^easily emulsified/**easily emulsified/g; } @charlist2;
	map {s/^flaky texture/**flaky texture/g; } @charlist2;
	map {s/glistening texture/**glistening texture/g; } @charlist2;
	map {s/^opalescent texture/**opalescent texture/g; } @charlist2;
	map {s/^opaque colony/**opaque colony/g; } @charlist2;
	map {s/^ropy texture/**ropy texture/g; } @charlist2;
	map {s/^rough texture/**rough texture/g; } @charlist2;
	map {s/^slimy texture/**slimy texture/g; } @charlist2;
	map {s/^smooth texture/**smooth texture/g; } @charlist2;
	map {s/^sunken texture/**sunken texture/g; } @charlist2;
	map {s/^translucent colony/**translucent colony/g; } @charlist2;
	map {s/^transparent colony/**transparent colony/g; } @charlist2;
	map {s/^viscid colony/**viscid colony/g; } @charlist2;
	map {s/^watery colony/**watery colony/g; } @charlist2;
#new
	map {s/^dry texture/**dry texture/g; } @charlist2;
	map {s/^granular texture/**granular texture/g; } @charlist2;
	map {s/^ridged texture/**ridged texture/g; } @charlist2;
	map {s/^wrinkled texture/**wrinkled texture/g; } @charlist2;

	print "\n\nBelow is your list of homologized characters:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.coltext.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.coltext.txt
		}
	my $temp6 = "temp6.coltext.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.coltext.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "coltext.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @coltextdata = split (/\t/, $line);
		push (@taxnames, $coltextdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Texture Matrix'\;\n\tDIMENSIONS NCHAR=17\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'butyrous texture' \/  'not butyrous texture' 'butyrous texture', ";
	print OUT "2 'adherent texture' \/  'not adherent texture' 'adherent texture', ";
	print OUT "3 'emulsification texture' \/  'not easily emulsified' 'easily emulsified', ";
	print OUT "4 'flaky texture' \/  'not flaky texture' 'flaky texture', ";
	print OUT "5 'glistening texture' \/  'dull texture' 'glistening texture', ";
	print OUT "6 'slimy texture' \/  'not slimy texture' 'slimy texture', ";
	print OUT "7 'opalescent texture' \/  'not opalescent texture' 'opalescent texture', ";
	print OUT "8 'aromatic colony' \/  'not aromatic colony' 'aromatic colony', ";
	print OUT "9 'smooth texture' \/  'rough texture' 'smooth texture', ";
	print OUT "10 'ropy texture' \/  'not ropy texture' 'ropy texture', ";
	print OUT "11 'translucent colony' \/  'opaque colony' 'translucent colony', ";
	print OUT "12 'viscid texture' \/  'watery texture' 'viscid texture', ";
	print OUT "13 'sunken texture' \/  'not sunken texture' 'sunken texture', ";
	print OUT "14 'dry texture' \/  'not dry texture' 'dry texture', ";
	print OUT "15 'granular texture' \/  'not granular texture' 'granular texture', ";
	print OUT "16 'ridged texture' \/  'not ridged texture' 'ridged texture', ";
	print OUT "17 'wrinkled texture' \/  'not wrinkled texture' 'wrinkled texture', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @coltextdata = split (/\t/, $line);
		push (my @taxnames, $coltextdata[0]);

#code char 1 butyrous texture
		if ($line =~ /,butyrous texture,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not butyrous texture,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 adherent texture
		if ($line =~ /,adherent texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adherent texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 easily emulsified
		if ($line =~ /,easily emulsified,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not easily emulsified,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 flaky texture
		if ($line =~ /,flaky texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flaky texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 glistening texture
		if ($line =~ /,glistening texture,|,not dull texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glistening texture,|,dull texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 slimy texture
		if ($line =~ /,slimy texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not slimy texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 opalescent texture
		if ($line =~ /,opalescent texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not opalescent texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 aromatic colony
		if ($line =~ /,aromatic colony,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aromatic colony,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 smooth texture
		if ($line =~ /,smooth texture,|,not rough texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not smooth texture,|,rough texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 ropy texture
		if ($line =~ /,ropy texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ropy texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 translucent colony
		if ($line =~ /,translucent colony,|,transparent colony,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not translucent colony|,not transparent colony,|,opaque colony,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 viscid colony
		if ($line =~ /,viscid colony,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not viscid colony,|,watery colony,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 sunken texture
		if ($line =~ /,sunken texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sunken texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}

#code char 13 dry texture
		if ($line =~ /,dry texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dry texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 14 granular texture
		if ($line =~ /,granular texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not granular texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 ridged texture
		if ($line =~ /,ridged texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ridged texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 wrinkled texture
		if ($line =~ /,wrinkled texture,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not wrinkled texture,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}




	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}
###################
#process colony color character
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "colony color") {
	my $homout = "hom.colcolor.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @colcolor, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @colcolor;
		map {s/,a creamy,/,white,/g; } @colcolor; # synonyms of 
		map {s/,beige,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,black pigment,/,black,/g; } @colcolor; # synonyms of 
		map {s/,bright golden yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,bright orange,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,bright pink,/,dark,red,/g; } @colcolor; # synonyms of 
		map {s/,bright yellow pigment,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,bright yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,bright,/,dark,/g; } @colcolor; # synonyms of 
		map {s/,brown deposit,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brown,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,buff,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,clear margins,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,clear,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,colorless,/,white,/g; } @colcolor; # synonyms of 
		map {s/,cream,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy,/,white,/g; } @colcolor; # synonyms of 
		map {s/,dark brown,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,dark orange,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,dark pink,/,dark,red,/g; } @colcolor; # synonyms of 
		map {s/,dark reddish orange,/,dark,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,dark yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,dark-orange-coloured,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,deep orange-brown,/,dark,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,deep,/,dark,/g; } @colcolor; # synonyms of 
		map {s/,dirty yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,dorothy white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,entire translucent margins,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,glitstening,/,glistening,/g; } @colcolor; # synonyms of 
		map {s/,golden-yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]y-white,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]y,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]yish-white,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]yish,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,greenish gr[ae]y,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,iridescence,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,iridescent lustre,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,light brown,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,light orange,/,pale,orange,/g; } @colcolor; # synonyms of 
		map {s/,light pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,light purple,/,pale,violet,/g; } @colcolor; # synonyms of 
		map {s/,light red,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,light salmon pink,/,pale,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,light yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,light-pink-coloured,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,light,/,pale,/g; } @colcolor; # synonyms of 
		map {s/,matt,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,medium yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,metallic sheen,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,moderately reddish-orange,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,much pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,non-watersoluble pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,off-white centre,/,white,/g; } @colcolor; # synonyms of 
		map {s/,off-white-grey,/,white,/g; } @colcolor; # synonyms of 
		map {s/,off-white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,orange coloured,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange fluorescence,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange pigment,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange-coloured,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange\/yellow,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,orangish red,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,pale orange-yellow,/,pale,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,pale orange,/,pale,orange,/g; } @colcolor; # synonyms of 
		map {s/,pale pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pale pinkish gray,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pale speckled-pink,/,pale,mottled,red,/g; } @colcolor; # synonyms of 
		map {s/,pale.yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,pigmentation,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,pigmented creamy white colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,pink-coloured,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pinkish,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pronounced metallic tinge,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,purple,/,violet,/g; } @colcolor; # synonyms of 
		map {s/,reddish brown,/,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,reddish-brown diffusible pigment,/,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,reddish-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,reddish,/,red,/g; } @colcolor; # synonyms of 
		map {s/,reddish.orange,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,rose-colored,/,red,/g; } @colcolor; # synonyms of 
		map {s/,rose-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,rusty-orange-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,saffron-colored,/,red,/g; } @colcolor; # synonyms of 
		map {s/,salmon-pink,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,salmon,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,semi-opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,semi-transparent,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,semiopaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,semitranslucent,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,shiny gr[ae]y periphery,/,gray,glistening,/g; } @colcolor; # synonyms of 
		map {s/,shiny gr[ae]y periphery,/,shiny,gray,/g; } @colcolor; # synonyms of 
		map {s/,shiny,/,glistening,/g; } @colcolor; # synonyms of 
		map {s/,slightly brown,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,a slight yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slightly motlted,/,mottled,/g; } @colcolor; # synonyms of 
		map {s/,slightly mottled,/,pale,mottled,/g; } @colcolor; # synonyms of 
		map {s/,slightly opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,slightly pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,slightly,/,pale,/g; } @colcolor; # synonyms of 
		map {s/,smooth white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,tan,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,translucent colonies,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,translucent yellow,/,translucent,yellow,/g; } @colcolor; # synonyms of 
		map {s/,translucent yellowish,/,translucent,yellow,/g; } @colcolor; # synonyms of 
		map {s/,translucent.whitish,/,translucent,white,/g; } @colcolor; # synonyms of 
		map {s/,trasnparent,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,water-soluble brown pigment,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,white fluorescent tubes,/,fluorescent,white,/g; } @colcolor; # synonyms of 
		map {s/,white-gr[ae]yish,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,whitegr[ae]yish,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,whitish,/,white,/g; } @colcolor; # synonyms of 
		map {s/,yellow colonies,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow flexibacteria,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow-pigmented,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow.orange,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellowish orange,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellowish,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,a yellowish,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,almost opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,beige aerial and substrate mycelium,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,black colonies,/,black,/g; } @colcolor; # synonyms of 
		map {s/,blue iridescence,/,blue,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,blueish green,/,blue,green,/g; } @colcolor; # synonyms of 
		map {s/,bluish appearance,/,blue,/g; } @colcolor; # synonyms of 
		map {s/,bluish grey,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,bright yellow color,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,brown-red pigment,/,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish aesculin-degrading colonies,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish colonies,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish pigment,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish soluble pigment,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brownish white,/,light,brown,/g; } @colcolor; # synonyms of 
		map {s/,chalk white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,circular and creamy white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,circular and pale yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,clear colonies,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,cloudy,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,colored,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,coloured,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,convex and pale orange,/,pale,orange,/g; } @colcolor; # synonyms of 
		map {s/,convex and pale yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,cream colour,/,white,/g; } @colcolor; # synonyms of 
		map {s/,cream pigmentation,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy grey-coloured,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,creamy pigment,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy white and 0.6 mm,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy white-pigmented,/,white,/g; } @colcolor; # synonyms of 
		map {s/,different pigment intensities,/,mottled,/g; } @colcolor; # synonyms of 
		map {s/,distinct yellow ring pattern,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,dull surface,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,dull yellow-orange,/,dull,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,eggshell,/,white,/g; } @colcolor; # synonyms of 
		map {s/,elevated and yellowish orange-coloured,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,faintly glistening,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glistening and 1-3 mm,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glistening surface,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glistening surfaces,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy appearance,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy textured surfaces,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy centres,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy surface,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy surfaces,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,glossy white,/,glossy,white,/g; } @colcolor; # synonyms of 
		map {s/,gray and gray-white rings,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,gray rings,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,grayish concentric ring pattern,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,grayish white,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,grayish yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,greyish white,/,light,gray,/g; } @colcolor; # synonyms of 
		map {s/,ivory,/,white,/g; } @colcolor; # synonyms of 
		map {s/,light gray,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,light grey,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,light orange-yellow,/,pale,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,light orange-yellow coloured,/,pale,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,light yellow color,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,light yellow-gray pigment,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,light yellowish-brown,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,light-brown pigment,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,lustreless and translucent white,/,not glossy,translucent,white,/g; } @colcolor; # synonyms of 
		map {s/,milky,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,milky white,/,translucent,white,/g; } @colcolor; # synonyms of 
		map {s/,minute white colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,olive green,/,dark,green,/g; } @colcolor; # synonyms of 
		map {s/,opalescent whitish,/,white,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,opaque and creamy and 2-5 mm,/,opaque,white,/g; } @colcolor; # synonyms of 
		map {s/,opaque and granular surface texture,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,opaque centers,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,opaque centres,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,opaque colonies,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,opaque white,/,opaque,white,/g; } @colcolor; # synonyms of 
		map {s/,opaque white and 2-6 mm,/,opaque,white,/g; } @colcolor; # synonyms of 
		map {s/,orangish brown,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,originally opaque gray-white,/,opaque,gray,/g; } @colcolor; # synonyms of 
		map {s/,other pigmented enterococci,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,pale pink,/,light,red,/g; } @colcolor; # synonyms of 
		map {s/,pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,pigmented and unpigmented clones,/,pigmented,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,pink tint,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,producing brilliantly colored colonies,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,producing deep pink colonies,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,slight darkening,/,dark,/g; } @colcolor; # synonyms of 
		map {s/,slight glossy surface,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,slight yellow tint,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slight yellowish tint,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slightly convex and moderate yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slightly dull,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,slightly translucent,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,slightly transparent,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,slightly yellowish,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slightly yellowish pigment,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,slightly yellowish tint,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,small white colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,smooth and pale yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,smooth dull,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,smooth dull surface,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,smooth glossy surfaces,/,glossy,/g; } @colcolor; # synonyms of 
		map {s/,soluble pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,sparse white aerial mycelia,/,white,/g; } @colcolor; # synonyms of 
		map {s/,tinted,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,translucent and creamy white,/,translucent,white,/g; } @colcolor; # synonyms of 
		map {s/,translucent and low convex,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,translucent edges,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,translucent fringe,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,transparent and low convex,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,usually bright yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,usually opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,usually pigmented,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,very slight tint,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,water insoluble pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,white aerial mycelium,/,white,/g; } @colcolor; # synonyms of 
		map {s/,white and <0.5 mm,/,white,/g; } @colcolor; # synonyms of 
		map {s/,white and substrate mycelium,/,white,/g; } @colcolor; # synonyms of 
		map {s/,white centre,/,white,/g; } @colcolor; # synonyms of 
		map {s/,white colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,white colour,/,white,/g; } @colcolor; # synonyms of 
		map {s/,clear colonies,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,white mycelia,/,white,/g; } @colcolor; # synonyms of 
		map {s/,whitish and low convex,/,white,/g; } @colcolor; # synonyms of 
		map {s/,whitish colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,yellow and 1-3 mm,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow and 1.0-1.5 mm,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow pigment,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow pigmentation,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow pigmented,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow ring,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow-gray pigment,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,yellow-orange,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellow-orange pigment,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellowish center,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellowish gray,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,yellowish orange pigment,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellowish pigment,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellowish tint,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellowish white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,soluble pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colcolor; # synonyms of 

		map {s/,non-pigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,non-translucent,/,not translucent,/g; } @colcolor; # synonyms of 
		map {s/,nonpigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,not black pigment,/,not black,/g; } @colcolor; # synonyms of 
		map {s/,not bright yellow pigment,/,not yellow,/g; } @colcolor; # synonyms of 
		map {s/,not clear colonies,/,not translucent,/g; } @colcolor; # synonyms of 
		map {s/,not pigment,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,not soluble pigment,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,not white fluorescent tubes,/,not fluorescent,not white,/g; } @colcolor; # synonyms of 
		map {s/,unpigmented colonies,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,unpigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,usually unpigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,without diffusible pigment,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,iridescent centre,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,not black colonies,/,not black,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,xxx,/g; } @colcolor; # synonyms of 

#not a colony color
		map {s/,0.1%methylene blue milk,/,/g; } @colcolor; # synonyms of 
		map {s/,5% horse blood,/,/g; } @colcolor; # synonyms of 
		map {s/,5% sheep blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,5% sheep blood,/,/g; } @colcolor; # synonyms of 
		map {s/,bhi blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,bhia streak tube,/,/g; } @colcolor; # synonyms of 
		map {s/,bhia streak tubes,/,/g; } @colcolor; # synonyms of 
		map {s/,black halo,/,/g; } @colcolor; # synonyms of 
		map {s/,blood agar plates,/,/g; } @colcolor; # synonyms of 
		map {s/,blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,not blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,blood,/,/g; } @colcolor; # synonyms of 
		map {s/,blue halo,/,/g; } @colcolor; # synonyms of 
		map {s/,bovine blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,bovine blood,/,/g; } @colcolor; # synonyms of 
		map {s/,brain heart infusion agar streak tubes,/,/g; } @colcolor; # synonyms of 
		map {s/,butyrous,/,/g; } @colcolor; # synonyms of 
		map {s/,clear \w,/,/g; } @colcolor; # synonyms of 
		map {s/,clear egg yolk medium,/,/g; } @colcolor; # synonyms of 
		map {s/,clear hemolysin activity,/,/g; } @colcolor; # synonyms of 
		map {s/,clear medium,/,/g; } @colcolor; # synonyms of 
		map {s/,clear supernatants,/,/g; } @colcolor; # synonyms of 
		map {s/,clear zone,/,/g; } @colcolor; # synonyms of 
		map {s/,clear zones,/,/g; } @colcolor; # synonyms of 
		map {s/,clearer,/,/g; } @colcolor; # synonyms of 
		map {s/,colombia blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,columbia blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,columbia blood,/,/g; } @colcolor; # synonyms of 
		map {s/,columbia horse blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,columbia sheep blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,crystal violet agar,/,/g; } @colcolor; # synonyms of 
		map {s/,glistening,/,/g; } @colcolor; # synonyms of 
		map {s/,horse blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,horse blood,/,/g; } @colcolor; # synonyms of 
		map {s/,human blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,human blood,/,/g; } @colcolor; # synonyms of 
		map {s/,inoculation streak,/,/g; } @colcolor; # synonyms of 
		map {s/,light growth,/,/g; } @colcolor; # synonyms of 
		map {s/,light,/,/g; } @colcolor; # synonyms of 
		map {s/,lighter and darker colors,/,/g; } @colcolor; # synonyms of 
		map {s/,litmus milk,/,/g; } @colcolor; # synonyms of 
		map {s/,not litmus milk,/,/g; } @colcolor; # synonyms of 
		map {s/,mastitic milk,/,/g; } @colcolor; # synonyms of 
		map {s/,methyl red indicator,/,/g; } @colcolor; # synonyms of 
		map {s/,methyl red,/,/g; } @colcolor; # synonyms of 
		map {s/,milk,/,/g; } @colcolor; # synonyms of 
		map {s/,neutral red,/,/g; } @colcolor; # synonyms of 
		map {s/,not 0.1%methylene blue milk,/,/g; } @colcolor; # synonyms of 
		map {s/,not blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,not glistening,/,/g; } @colcolor; # synonyms of 
		map {s/,not rabbit blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,not sheep blood,/,/g; } @colcolor; # synonyms of 
		map {s/,oblique light,/,/g; } @colcolor; # synonyms of 
		map {s/,oblique transmitted light,/,/g; } @colcolor; # synonyms of 
		map {s/,obliquely transmitted light,/,/g; } @colcolor; # synonyms of 
		map {s/,orange juice,/,/g; } @colcolor; # synonyms of 
		map {s/,punctate colonies,/,/g; } @colcolor; # synonyms of 
		map {s/,purple agar base,/,/g; } @colcolor; # synonyms of 
		map {s/,purple broth medium,/,/g; } @colcolor; # synonyms of 
		map {s/,rabbit blood agar plates,/,/g; } @colcolor; # synonyms of 
		map {s/,semi-transparent a-haemolysis,/,/g; } @colcolor; # synonyms of 
		map {s/,sheep blood agar,/,/g; } @colcolor; # synonyms of 
		map {s/,sheep blood,/,/g; } @colcolor; # synonyms of 
		map {s/,streak tubes,/,/g; } @colcolor; # synonyms of 
		map {s/,surrounding bovine blood agar medium,/,/g; } @colcolor; # synonyms of 
		map {s/,the deeper,/,/g; } @colcolor; # synonyms of 
		map {s/,white mice,/,/g; } @colcolor; # synonyms of 
		map {s/,white soapflakes,/,/g; } @colcolor; # synonyms of 
		map {s/,not blue halo,/,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,/g; } @colcolor; # synonyms of 
		map {s/,xxx,/,/g; } @colcolor; # synonyms of 
print @colcolor, "\n\n";

		print OUT @taxlabels, @colcolor, "\n"; # prints to $homout, hom.colcolor.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.colcolor.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/,colony color,//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.colcolor.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.colcolor.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.colcolor.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.colcolor.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.colcolor.txt
#temporarily rename charstates to label those that have been homologized		
	map {s/^pigmented/**pigmented/g; } @charlist2;

	map {s/^dark/**dark/g; } @charlist2;
	map {s/^pale/**pale/g; } @charlist2;

	map {s/^dull/**dull/g; } @charlist2;
	map {s/^mottled/**mottled/g; } @charlist2;

	map {s/^fluorescent/**fluorescent/g; } @charlist2;
	map {s/^iridescent/**iridescent/g; } @charlist2;

	map {s/^opaque/**opaque/g; } @charlist2;
	map {s/^translucent/**translucent/g; } @charlist2;
	map {s/^transparent/**transparent/g; } @charlist2;

	map {s/^black/**black/g; } @charlist2;
	map {s/^brown/**brown/g; } @charlist2;
	map {s/^gray/**gray/g; } @charlist2;
	map {s/^orange-red/**orange-red/g; } @charlist2;
	map {s/^orange/**orange/g; } @charlist2;
	map {s/^reddish-brown/**reddish-brown/g; } @charlist2;
	map {s/^red/**red/g; } @charlist2;
	map {s/^violet/**violet/g; } @charlist2;
	map {s/^white/**white/g; } @charlist2;
	map {s/^yellow-orange/**yellow-orange/g; } @charlist2;
	map {s/^yellow/**yellow/g; } @charlist2;
#new
	map {s/^blue/**blue/g; } @charlist2;
	map {s/^green/**green/g; } @charlist2;
	map {s/^glossy/**glossy/g; } @charlist2;

	print "\n\nBelow is your list of homologized characters:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.colcolor.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.colcolor.txt
		}
	my $temp6 = "temp6.colcolor.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.colcolor.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "colcolor.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @colcolordata = split (/\t/, $line);
		push (@taxnames, $colcolordata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Color Matrix'\;\n\tDIMENSIONS NCHAR=21\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 'brown colonies' \/  'colonies not brown' 'brown colonies', ";
	print OUT "2 'reddish-brown colonies' \/  'colonies not reddish-brown' 'reddish-brown colonies', ";
	print OUT "3 'red colonies' \/  'colonies not red' 'red colonies', ";
	print OUT "4 'orange-red colonies' \/  'colonies not orange-red' 'orange-red colonies', ";
	print OUT "5 'orange colonies' \/  'colonies not orange' 'orange colonies', ";
	print OUT "6 'yellow-orange colonies' \/  'colonies not yellow-orange' 'yellow-orange colonies', ";
	print OUT "7 'yellow colonies' \/  'colonies not yellow' 'yellow colonies', ";
	print OUT "8 'violet colonies' \/  'colonies not violet' 'violet colonies', ";
	print OUT "9 'black colonies' \/  'colonies not black' 'black colonies', ";
	print OUT "10 'gray colonies' \/  'colonies not gray' 'gray colonies', ";
	print OUT "11 'white colonies' \/  'colonies not white' 'white colonies', ";
	print OUT "12 'colony optical quality' \/  glistening 'dull colonies', ";
	print OUT "13 'colony opacity' \/  transparent translucent 'opaque colonies', ";
	print OUT "14 'colony fluorescence' \/  'colonies not fluorescent' 'fluorescent colonies', ";
	print OUT "15 'mottled colonies' \/  'colonies not mottled' 'mottled colonies', ";
	print OUT "16 'iridescent colonies' \/  'colonies not iridescent' 'iridescent colonies', ";
	print OUT "17 'colony color saturation' \/  'pale color saturation' 'medium color saturation' 'dark color saturation', ";
	print OUT "18 'colony pigmentation' \/  'colonies not pigmented' 'pigmented colonies', ";
#new
	print OUT "19 'blue colonies' \/  'colonies not blue' 'blue colonies', ";
	print OUT "20 'green colonies' \/  'colonies not green' 'green colonies', ";
	print OUT "21 'glossy colonies' \/  'colonies not glossy' 'glossy colonies', ";


	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @colcolordata = split (/\t/, $line);
		push (my @taxnames, $colcolordata[0]);

#code char 1 brown
		if ($line =~ /,brown,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not brown,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 reddish-brown
		if ($line =~ /,reddish-brown,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not reddish-brown,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 3 red
		if ($line =~ /,red,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not red,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 4 orange-red
		if ($line =~ /,orange-red,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not orange-red,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 5 orange
		if ($line =~ /,orange,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not orange,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 6 yellow-orange
		if ($line =~ /,yellow-orange,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not yellow-orange,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 7 yellow
		if ($line =~ /,yellow,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not yellow,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 8 violet
		if ($line =~ /,violet,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not violet,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 9 black
		if ($line =~ /,black,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not black,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 10 gray
		if ($line =~ /,gray,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gray,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 11 white
		if ($line =~ /,white,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not white,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 12 optical quality
		if ($line =~ /,dull,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dull,|,glistening,|,glossy,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 13 opacity
		if ($line =~ /,transparent.*translucent,|,translucent.*transparent,/) {
			print OUT "(0 1)";
			}
		elsif ($line =~ /,transparent.*opaque,|,opaque.*transparent,/) {
			print OUT "(0 2)";
			}
		elsif ($line =~ /,translucent.*opaque,|,opaque.*translucent,/) {
			print OUT "(1 2)";
			}
		elsif ($line =~ /,transparent,/) {
			print OUT "0";
			}
		elsif ($line =~ /,translucent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,opaque,/) {
			print OUT "2";
			}
		else {
			print OUT "?";
		}
#code char 14 fluorescent
		if ($line =~ /,fluorescent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fluorescent,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 15 mottled
		if ($line =~ /,mottled,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mottled,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 16 iridescent
		if ($line =~ /,iridescent,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not iridescent,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 17 color saturation
		if ($line =~ /,pale,/) {
			print OUT "0";
			}
		elsif ($line =~ /,dark,/) {
			print OUT "2";
			}
		elsif ($line =~ /,pigmented,|,blue,|,green,|,brown,|,reddish-brown,|,red,|,orange-red,|,orange,|,yellow-orange,|,yellow,|,violet,|,black,|,gray,/) {
			print OUT "1";
		}
		else {
			print OUT "?";
		}
#code char 18 pigmentation
		if ($line =~ /,not pigmented,/) {
			print OUT "0";
			}
		elsif ($line =~ /,pigmented,|,blue,|,green,|,brown,|,reddish-brown,|,red,|,orange-red,|,orange,|,yellow-orange,|,yellow,|,violet,|,black,|,gray,/) {
			print OUT "1";
			}
		else {
			print OUT "?";
		}

#code char 19 blue
		if ($line =~ /,blue,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not blue,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 20 green
		if ($line =~ /,green,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not green,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}
#code char 21 glossy
		if ($line =~ /,glossy,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glossy,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
		}


	print OUT "\n";
	}
	print OUT "\n\;\nEND\;\n";
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
}

###################
#process tests positive
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "tests positive") {
	my $homout = "tests.positive.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		$line =~ s/\xce\xb2/b/g; #changes greek betas β to b
		$line =~ s/\xc3\x9f/b/g; #changes greek eszett β to b
		$line =~ s/\xce\xb1/a/g;#changes greek alphas α to a
		$line =~ s/\xce\x91/a/g;#changes greek capital alphas to a
		$line =~ s/\xce\x92/b/g;#changes greek capital betas to b
		$line =~ s/\xc2\xb1/ /g;#changes plus or minus to space
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @testspos, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @testspos;
		map {s/,\d\.\d\%\s,/,,/gi; } @testspos; # synonyms of 
		map {s/,\d\.\d\d\%\s,/,,/gi; } @testspos; # synonyms of 


	#substrates used and/or fermented
		map {s/,.*% fructose,/,fructose,/gi; } @testspos; # synonyms of 
		map {s/,.*% glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,.*% l-lactate,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,.*% l-lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,.*% potassium tellurite,/,tellurite,/gi; } @testspos; # synonyms of 
		map {s/,.*% sucrose,/,sucrose,/gi; } @testspos; # synonyms of 
		map {s/,.*% tween 80,/,tween 80,/gi; } @testspos; # synonyms of 
		map {s/,.*mm succinate,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d\d% formate,/,formate,/gi; } @testspos; # synonyms of 
		map {s/,1 2 propanediol,/,1-2-propanediol,/gi; } @testspos; # synonyms of 
		map {s/,1 2-propanediol,/,1-2-propanediol,/gi; } @testspos; # synonyms of 
		map {s/,1 3-butanediol,/,1-3-butanediol,/gi; } @testspos; # synonyms of 
		map {s/,1 3-propanediol,/,1-3-propanediol,/gi; } @testspos; # synonyms of 
		map {s/,1 g glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,2 3 butanediol,/,2-3-butanediol,/gi; } @testspos; # synonyms of 
		map {s/,2 3-butanediol,/,2-3-butanediol,/gi; } @testspos; # synonyms of 
		map {s/,2 moles lactate,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,2-keto-gluconate,/,2-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,2-ketoglutarate,/,2-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,20 mm pyruvate,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,25 mm leucine,/,leucine,/gi; } @testspos; # synonyms of 
		map {s/,25 mm phenylalanine,/,phenylalanine,/gi; } @testspos; # synonyms of 
		map {s/,25 mm threonine,/,threonine,/gi; } @testspos; # synonyms of 
		map {s/,25 mm valine,/,valine,/gi; } @testspos; # synonyms of 
		map {s/,3-hydroxybenzoic acid,/,3-hydroxybenzoate,/gi; } @testspos; # synonyms of 
		map {s/,3-hydroxybutyrate,/,3-hydroxybutyrate,/gi; } @testspos; # synonyms of 
		map {s/,4-hydroxybenzoic acid,/,4-hydroxybenzoate,/gi; } @testspos; # synonyms of 
		map {s/,5-keto-d gluconate,/,5-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,a-d-glucose 1-phosphate,/,glucose 1-phosphate,/gi; } @testspos; # synonyms of 
		map {s/,a-d-glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,a-ketoglutaric acid,/,a-ketoglutarate,/gi; } @testspos; # synonyms of 
		map {s/,a-ketovaleric acid,/,a-ketovalerate,/gi; } @testspos; # synonyms of 
		map {s/,acetate andethanol,/,acetate,ethanol,/gi; } @testspos; # synonyms of 
		map {s/,acetate strip test,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,acetate synthesis,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,acetate utilization,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,acetate,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,acetic acid,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,acetic and lactic acid,/,acetate,lactate,/gi; } @testspos; # synonyms of 
		map {s/,acetoacetic acid,/,acetoacetate,/gi; } @testspos; # synonyms of 
		map {s/,additional glutamic acid,/,glutamate,/gi; } @testspos; # synonyms of 
		map {s/,adenosine 5,/,adenosine 5\'-monophosphate,/gi; } @testspos; # synonyms of 
		map {s/,adipate,/,adipate,/gi; } @testspos; # synonyms of 
		map {s/,adipic acid,/,adipate,/gi; } @testspos; # synonyms of 
		map {s/,alanine,/,alanine,/gi; } @testspos; # synonyms of 
		map {s/,amygdaline,/,amygdalin,/gi; } @testspos; # synonyms of 
		map {s/,anaerobically d-ribose,/,ribose,/gi; } @testspos; # synonyms of 
		map {s/,anaerobically glucose,/,ribose,/gi; } @testspos; # synonyms of 
		map {s/,approximately 1 mol l-alanine,/,alanine,/gi; } @testspos; # synonyms of 
		map {s/,arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,aspartic acid,/,aspartate,/gi; } @testspos; # synonyms of 
		map {s/,b-methyl-d-glucoside,/,methyl b-d-glucoside,/gi; } @testspos; # synonyms of 
		map {s/,b-phenylethylamine,/,phenylethylamine,/gi; } @testspos; # synonyms of 
		map {s/,beijerinckii raffinose,/,raffinose,/gi; } @testspos; # synonyms of 
		map {s/,biomass using malonate,/,malonate,/gi; } @testspos; # synonyms of 
		map {s/,both d-and l-lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,butyrate plus co2,/,butyrate,co2,/gi; } @testspos; # synonyms of 
		map {s/,butyric acid,/,butyrate,/gi; } @testspos; # synonyms of 
		map {s/,caprate,/,caprate,/gi; } @testspos; # synonyms of 
		map {s/,capric acid assimilation,/,caprate,/gi; } @testspos; # synonyms of 
		map {s/,capric acid,/,caprate,/gi; } @testspos; # synonyms of 
		map {s/,carbohydrate derivatives including glucose,/,carbohydrates,glucose,/gi; } @testspos; # synonyms of 
		map {s/,carbohydrates including arabinose,/,carbohydrates,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,carbohydrates ribose,/,carbohydrates,ribose,/gi; } @testspos; # synonyms of 
		map {s/,cell suspensions ferment glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose fermentation,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose help,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,chemoorganoheterotroph utilizing cellobiose,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,citrate fermentation,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,citrate hydrolysis,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,citrate utilization,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,citrate,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,citric acid,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,corresponding fatty acid plus acetic acid,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,d cellobiose,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,d fructose,/,fructose,/gi; } @testspos; # synonyms of 
		map {s/,d glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,d lactose,/,lactose,/gi; } @testspos; # synonyms of 
		map {s/,d mannitol,/,mannitol,/gi; } @testspos; # synonyms of 
		map {s/,d mannose,/,mannose,/gi; } @testspos; # synonyms of 
		map {s/,d trehalose,/,trehalose,/gi; } @testspos; # synonyms of 
		map {s/,d xylose,/,xylose,/gi; } @testspos; # synonyms of 
		map {s/,d-adonitol,/,ribitol,/gi; } @testspos; # synonyms of 
		map {s/,d-alanine,/,alanine,/gi; } @testspos; # synonyms of 
		map {s/,d-amygdalin,/,amygdalin,/gi; } @testspos; # synonyms of 
		map {s/,d-arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,d-arabitol,/,arabitol,/gi; } @testspos; # synonyms of 
		map {s/,d-aspartic acid,/,aspartate,/gi; } @testspos; # synonyms of 
		map {s/,d-cellobiose,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,d-fructose 6-phosphate,/,fructose 6-phosphate,/gi; } @testspos; # synonyms of 
		map {s/,d-fructose-6-po4,/,fructose 6-phosphate,/gi; } @testspos; # synonyms of 
		map {s/,d-fructose,/,fructose,/gi; } @testspos; # synonyms of 
		map {s/,d-fucose,/,fucose,/gi; } @testspos; # synonyms of 
		map {s/,d-galactose,/,galactose,/gi; } @testspos; # synonyms of 
		map {s/,d-galacturonate,/,galacturonate,/gi; } @testspos; # synonyms of 
		map {s/,d-galacturonic acid,/,galacturonate,/gi; } @testspos; # synonyms of 
		map {s/,d-gluconate,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,d-gluconic acid,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,d-glucose 1-phosphate,/,glucose 1-phosphate,/gi; } @testspos; # synonyms of 
		map {s/,d-glucose fermentation,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,d-glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,d-glucuronate,/,glucuronate,/gi; } @testspos; # synonyms of 
		map {s/,d-glucuronic acid,/,glucuronate,/gi; } @testspos; # synonyms of 
		map {s/,d-lactic acid methyl ester,/,lactic acid methyl ester,/gi; } @testspos; # synonyms of 
		map {s/,d-lactose,/,lactose,/gi; } @testspos; # synonyms of 
		map {s/,d-lyxose,/,lyxose,/gi; } @testspos; # synonyms of 
		map {s/,d-malate,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,d-malic acid,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,d-maltose,/,maltose,/gi; } @testspos; # synonyms of 
		map {s/,d-mannitol,/,mannitol,/gi; } @testspos; # synonyms of 
		map {s/,d-mannose,/,mannose,/gi; } @testspos; # synonyms of 
		map {s/,d-melezitose,/,melezitose,/gi; } @testspos; # synonyms of 
		map {s/,d-melibiose,/,melibiose,/gi; } @testspos; # synonyms of 
		map {s/,d-psicose,/,psicose,/gi; } @testspos; # synonyms of 
		map {s/,d-raffinose,/,raffinose,/gi; } @testspos; # synonyms of 
		map {s/,d-ribitol,/,ribitol,/gi; } @testspos; # synonyms of 
		map {s/,d-ribose,/,ribose,/gi; } @testspos; # synonyms of 
		map {s/,d-saccharose,/,sucrose,/gi; } @testspos; # synonyms of 
		map {s/,d-salicin,/,salicin,/gi; } @testspos; # synonyms of 
		map {s/,d-sorbitol,/,sorbitol,/gi; } @testspos; # synonyms of 
		map {s/,d-tagatose,/,tagatose,/gi; } @testspos; # synonyms of 
		map {s/,d-trehalose,/,trehalose,/gi; } @testspos; # synonyms of 
		map {s/,d-turanose,/,turanose,/gi; } @testspos; # synonyms of 
		map {s/,d-xylose,/,xylose,/gi; } @testspos; # synonyms of 
		map {s/,dextrin crystals,/,dextrin,/gi; } @testspos; # synonyms of 
		map {s/,dextrin utilization,/,dextrin,/gi; } @testspos; # synonyms of 
		map {s/,disodium succinate,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,dl-carnitine,/,carnitine,/gi; } @testspos; # synonyms of 
		map {s/,dl-lactate,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,either l-arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,esculin nor starch,/,aesculin,starch,/gi; } @testspos; # synonyms of 
		map {s/,exogenous citrate,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,fewer strains ferment cellobiose,/,cellobiose,/gi; } @testspos; # synonyms of 
		map {s/,formate being major end products,/,formate,/gi; } @testspos; # synonyms of 
		map {s/,formate yield succinate,/,formate,succinate,/gi; } @testspos; # synonyms of 
		map {s/,formic acid,/,formate,/gi; } @testspos; # synonyms of 
		map {s/,formic and lactic acid,/,formate,lactate,/gi; } @testspos; # synonyms of 
		map {s/,galactitol,/,dulcitol,/gi; } @testspos; # synonyms of 
		map {s/,gluconic acid sodium salt,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,gluconic acid,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,glucose acidification,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,glucose fermentation products,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,glucose fermentation,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,glucose metabolism,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,glucose supplementation,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,glucuronic acid,/,glucuronate,/gi; } @testspos; # synonyms of 
		map {s/,glutamic acid,/,glutamate,/gi; } @testspos; # synonyms of 
		map {s/,glycerol 1 3-propanediol,/,glycerol,1-3-propanediol,/gi; } @testspos; # synonyms of 
		map {s/,glycogen utilization,/,glycogen,/gi; } @testspos; # synonyms of 
		map {s/,glycogene,/,glycogen,/gi; } @testspos; # synonyms of 
		map {s/,h2s production,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,high xylan content,/,xylan,/gi; } @testspos; # synonyms of 
		map {s/,hydrogen sulfide,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,hydrogen sulphide production,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses gelatin,/,gelatin,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses starch,/,starch,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses tween 40,/,tween 40,/gi; } @testspos; # synonyms of 
		map {s/,including acetate,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,including arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,indole production,/,indole,/gi; } @testspos; # synonyms of 
		map {s/,itaconic acid,/,itaconate,/gi; } @testspos; # synonyms of 
		map {s/,keto butyric acid,/,2-ketobutyrate,/gi; } @testspos; # synonyms of 
		map {s/,l-alaninamide,/,alaninamide,/gi; } @testspos; # synonyms of 
		map {s/,l-alanine,/,alanine,/gi; } @testspos; # synonyms of 
		map {s/,L-alanine,/,alanine,/gi; } @testspos; # synonyms of 
		map {s/,l-alanyl-glycine,/,alanyl glycine,/gi; } @testspos; # synonyms of 
		map {s/,l-arabinose d-maltose,/,arabinose,maltose,/gi; } @testspos; # synonyms of 
		map {s/,l-arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,l-arabitol,/,arabitol,/gi; } @testspos; # synonyms of 
		map {s/,l-arginine,/,arginine,/gi; } @testspos; # synonyms of 
		map {s/,l-asparagine,/,asparagine,/gi; } @testspos; # synonyms of 
		map {s/,l-aspartate,/,aspartate,/gi; } @testspos; # synonyms of 
		map {s/,l-aspartic acid,/,aspartate,/gi; } @testspos; # synonyms of 
		map {s/,l-fucose,/,fucose,/gi; } @testspos; # synonyms of 
		map {s/,l-glutamate,/,glutamate,/gi; } @testspos; # synonyms of 
		map {s/,l-glutamic acid,/,glutamate,/gi; } @testspos; # synonyms of 
		map {s/,l-glutamine,/,glutamine,/gi; } @testspos; # synonyms of 
		map {s/,l-histidine,/,histidine,/gi; } @testspos; # synonyms of 
		map {s/,l-hydroxyproline,/,hydroxyproline,/gi; } @testspos; # synonyms of 
		map {s/,l-lactate isomers,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,l-lactate,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,l-lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,l-leucine,/,leucine,/gi; } @testspos; # synonyms of 
		map {s/,l-malate,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,l-malic acid,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,l-methionine,/,methionine,/gi; } @testspos; # synonyms of 
		map {s/,l-ornithine,/,ornithine,/gi; } @testspos; # synonyms of 
		map {s/,l-phenylalanine,/,phenylalanine,/gi; } @testspos; # synonyms of 
		map {s/,l-proline,/,proline,/gi; } @testspos; # synonyms of 
		map {s/,l-pyroglutamic acid,/,pyroglutamate,/gi; } @testspos; # synonyms of 
		map {s/,l-rhamnose fermentation,/,rhamnose,/gi; } @testspos; # synonyms of 
		map {s/,l-rhamnose,/,rhamnose,/gi; } @testspos; # synonyms of 
		map {s/,l-serine,/,serine,/gi; } @testspos; # synonyms of 
		map {s/,l-sorbose,/,sorbose,/gi; } @testspos; # synonyms of 
		map {s/,l-threonine,/,threonine,/gi; } @testspos; # synonyms of 
		map {s/,l-valine,/,valine,/gi; } @testspos; # synonyms of 
		map {s/,l-xylose,/,xylose,/gi; } @testspos; # synonyms of 
		map {s/,lactate either,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,lactic acid group,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,lactic acid isomers,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,lactic acid production,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,lactose help,/,lactose,/gi; } @testspos; # synonyms of 
		map {s/,little glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,m-inositol,/,myo-inositol,/gi; } @testspos; # synonyms of 
		map {s/,mainly acetate,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,mainly butyrate,/,butyrate,/gi; } @testspos; # synonyms of 
		map {s/,mainly lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,major lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,malic acid assimilation,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,malic acid,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,malonate utilization,/,malonate,/gi; } @testspos; # synonyms of 
		map {s/,maltose assimilation,/,maltose,/gi; } @testspos; # synonyms of 
		map {s/,mannose fermentation,/,mannose,/gi; } @testspos; # synonyms of 
		map {s/,mannose support weak growth,/,mannose,/gi; } @testspos; # synonyms of 
		map {s/,methyl a-d-glucopyranoside,/,methyl a-d-glucoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl a-d-glucose,/,3-methyl d-glucose,/gi; } @testspos; # synonyms of 
		map {s/,methyl a-d-mannopyranoside,/,methyl a-d-mannoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl b-d-galactopyranoside,/,methyl b-d-galactoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl b-d-xylopyranoside,/,methyl b-d-xyloside,/gi; } @testspos; # synonyms of 
		map {s/,methyl pyruvate,/,methyl pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,mol glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,mole glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,most strains ferment glucose,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,most strains ferment mannitol,/,mannitol,/gi; } @testspos; # synonyms of 
		map {s/,mucic acid,/,mucate,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-b-glucosaminidase activities,/,n-acetyl-b-d-glucosamine,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-d-glucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-glucosamine assimilation,/,xxx,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-glucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-l-glutamic acid,/,n-acetyl-l-glutamate,/gi; } @testspos; # synonyms of 
		map {s/,n-acetylglucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testspos; # synonyms of 
		map {s/,negative and gelatin hydrolysis,/,gelatin,/gi; } @testspos; # synonyms of 
		map {s/,negative trehalose fermentation,/,trehalose,/gi; } @testspos; # synonyms of 
		map {s/,nor fumarate,/,fumarate,/gi; } @testspos; # synonyms of 
		map {s/,oat spelt xylan,/,xylan,/gi; } @testspos; # synonyms of 
		map {s/,only a-ketovaleric acid,/,a-ketovalerate,/gi; } @testspos; # synonyms of 
		map {s/,only acetic acid,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,only d-ribose,/,ribose,/gi; } @testspos; # synonyms of 
		map {s/,only gelatin,/,gelatin,/gi; } @testspos; # synonyms of 
		map {s/,only glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,only l-lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,only lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,organic acids pyruvate,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,ornithine using moeller,/,ornithine,/gi; } @testspos; # synonyms of 
		map {s/,p-hydroxyphenylacetic acid,/,4-hydroxyphenylacetate,/gi; } @testspos; # synonyms of 
		map {s/,phenyl acetic acid,/,phenylacetate,/gi; } @testspos; # synonyms of 
		map {s/,phenylacetic acid,/,phenylacetate,/gi; } @testspos; # synonyms of 
		map {s/,poorly glycerol,/,glycerol,/gi; } @testspos; # synonyms of 
		map {s/,potassium 2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,potassium 5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testspos; # synonyms of 
		map {s/,potassium gluconate assimi-lation,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,potassium gluconate assimilation,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,potassium gluconate,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,potassium tellurite,/,tellurite,/gi; } @testspos; # synonyms of 
		map {s/,prefers pyruvate,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,propionic acid,/,propionate,/gi; } @testspos; # synonyms of 
		map {s/,putrescine fermentation,/,putrescine,/gi; } @testspos; # synonyms of 
		map {s/,pyroglutamic acid,/,pyroglutamate,/gi; } @testspos; # synonyms of 
		map {s/,pyruvate fermentation,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,pyruvate utilization,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,pyruvic acid methyl ester,/,methyl pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,pyruvic acid,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,quinic acid,/,quinate,/gi; } @testspos; # synonyms of 
		map {s/,raffinose cultures,/,raffinose,/gi; } @testspos; # synonyms of 
		map {s/,raffinose d-fructose,/,raffinose,fructose,/gi; } @testspos; # synonyms of 
		map {s/,reduces fumarate,/,fumarate,/gi; } @testspos; # synonyms of 
		map {s/,rhamnose utilization,/,rhamnose,/gi; } @testspos; # synonyms of 
		map {s/,s226 fermented glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,saccharose,/,sucrose,/gi; } @testspos; # synonyms of 
		map {s/,salicine d-sorbitol,/,salicin,sorbitol,/gi; } @testspos; # synonyms of 
		map {s/,sodium acetate,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,sodium citrate,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,sodium gluconate,/,gluconate,/gi; } @testspos; # synonyms of 
		map {s/,sodium glutamate,/,glutamate,/gi; } @testspos; # synonyms of 
		map {s/,sodium hippurate rind esculin,/,hippurate,aesculin,/gi; } @testspos; # synonyms of 
		map {s/,sodium hippurate,/,hippurate,/gi; } @testspos; # synonyms of 
		map {s/,sodium l-lactate,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,sodium malate,/,malate,/gi; } @testspos; # synonyms of 
		map {s/,sodium pyruvate,/,pyruvate,/gi; } @testspos; # synonyms of 
		map {s/,sodium succinate,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,solely l-lactic acid,/,lactate,/gi; } @testspos; # synonyms of 
		map {s/,soluble potato starch,/,starch,/gi; } @testspos; # synonyms of 
		map {s/,soluble starch,/,starch,/gi; } @testspos; # synonyms of 
		map {s/,sometimes formate,/,formate,/gi; } @testspos; # synonyms of 
		map {s/,sometimes succinate,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,starch hydrolysis,/,starch,/gi; } @testspos; # synonyms of 
		map {s/,starch nor esculin,/,starch,aesculin,/gi; } @testspos; # synonyms of 
		map {s/,starch utilization,/,starch,/gi; } @testspos; # synonyms of 
		map {s/,strains ferment arabinose,/,arabinose,/gi; } @testspos; # synonyms of 
		map {s/,strains ferment fructose,/,fructose,/gi; } @testspos; # synonyms of 
		map {s/,strains ferment lactose,/,lactose,/gi; } @testspos; # synonyms of 
		map {s/,strains weakly ferment glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,suberic acid,/,suberate,/gi; } @testspos; # synonyms of 
		map {s/,succinic acid monomethyl ester,/,monomethyl succinate,/gi; } @testspos; # synonyms of 
		map {s/,succinic acid,/,succinate,/gi; } @testspos; # synonyms of 
		map {s/,sucrose fermentation,/,sucrose,/gi; } @testspos; # synonyms of 
		map {s/,thymidine 5,/,thymidine 5\'-monophosphate,/gi; } @testspos; # synonyms of 
		map {s/,trisodium citrate assimilation,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,trisodium citrate tests,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,trisodium citrate,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,tween 80 hydrolysis variable,/,tween 80,/gi; } @testspos; # synonyms of 
		map {s/,tween 80 hydrolysis,/,tween 80,/gi; } @testspos; # synonyms of 
		map {s/,tween 80 type,/,tween 80,/gi; } @testspos; # synonyms of 
		map {s/,uridine 5,/,uridine 5\'-monophosphate,/gi; } @testspos; # synonyms of 
		map {s/,using glycerol,/,glycerol,/gi; } @testspos; # synonyms of 
		map {s/,usually acetate,/,acetate,/gi; } @testspos; # synonyms of 
		map {s/,usually formate,/,formate,/gi; } @testspos; # synonyms of 
		map {s/,usually glucose,/,glucose,/gi; } @testspos; # synonyms of 
		map {s/,utilizes citrate,/,citrate,/gi; } @testspos; # synonyms of 
		map {s/,valeric acid,/,valerate,/gi; } @testspos; # synonyms of 
		map {s/,valine,/,valine,/gi; } @testspos; # synonyms of 
		map {s/,weakly saccharolytic and ferment fructose,/,saccharolytic,fructose,/gi; } @testspos; # synonyms of 
		map {s/,xylan fermentation,/,xylan,/gi; } @testspos; # synonyms of 
		map {s/,xylan ribose,/,xylan,ribose,/gi; } @testspos; # synonyms of 
		map {s/,xylan support growth,/,xylan,/gi; } @testspos; # synonyms of 
		map {s/,xylose fermentation,/,xylose,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testspos; # synonyms of 

	#inorganics
		map {s/,aerobic nitrate reduction,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,hydrolysed and nitrate reduction,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate and indole production,/,nitrate reduction,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,nitrate and nitrite reduction,/,nitrate reduction,nitrite reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reductase a.,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reductase activities,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reductase,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction activity,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction reaction,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction reactions,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction test,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction tests,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction variable,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction-positive,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate respiration,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate supplementation,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrite and nitrate reduction,/,nitrate reduction,nitrite reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrite reductase activity,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nitrite reductase,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,nor nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,reduced nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,reducing nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,sodium nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,using nitrate,/,nitrate reduction,/gi; } @testspos; # synonyms of 

		map {s/,h2s formation,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,h2s gas,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,hydrogen sulfide production,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,hydrogen sulphide,/,h2s,/gi; } @testspos; # synonyms of 
		map {s/,moderate h2s,/,h2s,/gi; } @testspos; # synonyms of 

	#enzymes present
		map {s/,aesculin ferric citrate reaction,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,aesculin ferric citrate,/,aesculin hydrolysi,/gi; } @testspos; # synonyms of 
		map {s/,aesculin ferric,/,aesculin hydrolysi,/gi; } @testspos; # synonyms of 
		map {s/,aesculin hydrolysi,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,aesculin hydrolysis,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,aesculin,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,bile . aesculin test,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,esculin ferric citrate,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,esculin hydrolysis,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,esculin,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,ferric citrate,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses aesculin,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses esculin,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,positive aesculin hydrolysis,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,positive aesculin reaction,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,type strain hydrolyses aesculin,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 
		map {s/,aesculin dihydrolase,/,aesculin hydrolysis,/gi; } @testspos; # synonyms of 

		map {s/,cytochrome oxidase activities,/,cytochrome oxidase,/gi; } @testspos; # synonyms of 

		map {s/,2-naphthyl phosphate,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,2-naphthyl phosphate,/,phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,acid phosphatase activities,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,acid phosphatase no3,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,acid phosphatase,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,acidic phosphatase,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,bis-p-nitrophenyl-phosphate,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,generally weak alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,inconstant alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,moderate alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,moderate alkaline phosphatase,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,only acid phosphatase,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,p-nitrophenyl phosphate,/,phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,phosphatase acid activity,/,phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,phosphatase activity,/,phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,phosphatase reaction,/,phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,phosphatase-positive and ferment lactose,/,phosphatase,lactose,/gi; } @testspos; # synonyms of 
		map {s/,present and acid phosphatase,/,acid phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,weak alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 
		map {s/,weak alkaline phosphatase,/,alkaline phosphatase,/gi; } @testspos; # synonyms of 

		map {s/,c14 lipase,/,lipase c14,/gi; } @testspos; # synonyms of 
		map {s/,lipase c-14,/,lipase c14,/gi; } @testspos; # synonyms of 
		map {s/,lipase c 14,/,lipase c14,/gi; } @testspos; # synonyms of 
		map {s/,2-naphthyl myristate,/,lipase c14,/gi; } @testspos; # synonyms of 

		map {s/,2-naphthyl butyrate,/,esterase c4,/gi; } @testspos; # synonyms of 
		map {s/,c-4 esterase,/,esterase c4,/gi; } @testspos; # synonyms of 
		map {s/,c4 esterase,/,esterase c4,/gi; } @testspos; # synonyms of 
		map {s/,esterase c-4,/,esterase c4,/gi; } @testspos; # synonyms of 
		map {s/,esterase polymorphism,/,esterase,/gi; } @testspos; # synonyms of 
		map {s/,esterase production,/,esterase,/gi; } @testspos; # synonyms of 
		map {s/,esterase,/,esterase c4,/gi; } @testspos; # synonyms of 
		map {s/,negative including esterase,/,esterase,/gi; } @testspos; # synonyms of 

		map {s/,2-naphthyl caprylate,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,c-8 esterase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,c8 esterase lipase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,c8 esterase lipase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,c8 esterase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,c8 lipase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,ester lipase c-8,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,ester lipase c8,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,ester lipase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,esterase lipase c8,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,esterase lipase,/,lipase c8,/gi; } @testspos; # synonyms of 
		map {s/,present and esterase lipase,/,lipase c8,/gi; } @testspos; # synonyms of 

		map {s/,egg yolk lecithinase,/,lecithinase,/gi; } @testspos; # synonyms of 
		map {s/,lecithinase activity,/,lecithinase,/gi; } @testspos; # synonyms of 
		map {s/,egg-yolk lecithinase,/,lecithinase,/gi; } @testspos; # synonyms of 

		map {s/,gelatin hydrolysis,/,gelatinase,/gi; } @testspos; # synonyms of 
		map {s/,gelatin liquefaction,/,gelatinase,/gi; } @testspos; # synonyms of 
		map {s/,gelatinase,/,gelatinase,/gi; } @testspos; # synonyms of 
		map {s/,gelatinase activity,/,gelatinase,/gi; } @testspos; # synonyms of 
		map {s/,gelatinase activities,/,gelatinase,/gi; } @testspos; # synonyms of 
		map {s/,gelatin decomposition,/,gelatinase,/gi; } @testspos; # synonyms of 

		map {s/,n-acetyl-beta-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-b-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testspos; # synonyms of 

	#amino acid decarboxylases
		map {s/,arginine decarboxylase activity,/,arginine decarboxylase,/gi; } @testspos; # synonyms of 
		map {s/,l-ornithine decarboxylases,/,ornithine decarboxylase,/gi; } @testspos; # synonyms of 
		map {s/,lysine & ornithine decarboxylases,/,ornithine decarboxylase,lysine decarboxylase,/gi; } @testspos; # synonyms of 
		map {s/,lysine decarboxylase tests,/,lysine decarboxylase,/gi; } @testspos; # synonyms of 
		map {s/,ornithine decarboxylase activities,/,ornithine decarboxylase,/gi; } @testspos; # synonyms of 
		map {s/,ornithine decarboxylases,/,ornithine decarboxylase,/gi; } @testspos; # synonyms of 

	#amino acid deaminases
		map {s/,alanine dihydrolase activity,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,alanine dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dehydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase 2,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase activities,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase activity,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase hydrolysis,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase oxidase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase production,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase tests,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,arginine-dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,glutamate dehydrogenase,/,glutamate deaminase,/gi; } @testspos; # synonyms of 
		map {s/,indole activities,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole and nitrate reduction activity,/,tryptophanase,nitrate reduction,/gi; } @testspos; # synonyms of 
		map {s/,indole hydrolysis,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole negative,/,not tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole negative,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole nor h2s,/,tryptophanase,h2s,/gi; } @testspos; # synonyms of 
		map {s/,indole positive,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole reactions,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole test,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,indole tests,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,inducible tryptophanase,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,l-arginine dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,l-arginine hydrolysis,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,l-glutamate dehydrogenase,/,glutamate deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine deaminase activities,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine deaminase activity,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine deaminase tests,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine deamination,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine dehydrogenase,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,phenylalanine hydrolysis,/,phenylalanine deaminase,/gi; } @testspos; # synonyms of 
		map {s/,positive indole reaction,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,produces indole,/,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,tryptophan arylamidase,/,tryptophan deaminase,/gi; } @testspos; # synonyms of 
		map {s/,tryptophan deaminase activities,/,tryptophan deaminase,/gi; } @testspos; # synonyms of 
		map {s/,tryptophan deaminase activity,/,tryptophan deaminase,/gi; } @testspos; # synonyms of 
		map {s/,tryptophane deaminase,/,tryptophan deaminase,/gi; } @testspos; # synonyms of 
		map {s/,tyrosine arylamidase activities,/,tyrosine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testspos; # synonyms of 

	#amino acid arylamidases
		map {s/,alanine amd,/,alanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,arginine arylamidase activity,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,arginine arylamidase tests,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,arginine arylamidase,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,arginine dihydrolase hydrolysis,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,cystine aminopeptidase,/,cystine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,cystine arylamidase activities,/,cystine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,cystine arylamidase activity,/,cystine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,enzymes arginine arylamidase,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,glycine 7-amido-4-methylcoumarin,/,glycine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,glycine amd,/,glycine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,histidine 7-amido-4-methylcoumarin,/,histidine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,histidine amd,/,histidine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,histidine arylamidase activities,/,histidine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,histidine arylamidase activity,/,histidine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,inconstant alanine arylamidase activities,/,alanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-alanine amino-peptidase,/,alanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-arginine 7-amido-4-methylcoumarin,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-arginine amc,/,arginine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-aspartate arylamidase,/,aspartate arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-isoleucine arylamidase,/,isoleucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-phenylalanine 7-amido-4-methylcoumarin,/,phenylalanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-phenylalanine 7-amido-4,/,phenylalanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-phenylalanine amc,/,phenylalanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-proline arylamidase,/,proline arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-pyroglutamic acid 7-amido-4-methyl coumarin,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-pyroglutamic acid 7-amido-4-methylcoumarin,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-pyroglutamic acid amc,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-valine 7-amido-4-methylcoumarin,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,l-valine amc,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine aminopeptidase activity,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine aminopeptidase,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylamidase activities,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylamidase activity,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylamidase reaction,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylamidase tests,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylamidases,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylaminidase,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine arylaminidases,/,leucine arylaminidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine b-naphthylamide,/,leucine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine b-naphthylamide,/,leucine arylaminidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine p-nitroanilide,/,leucine arylaminidase,/gi; } @testspos; # synonyms of 
		map {s/,leucine pnitroanilide,/,leucine arylaminidase,/gi; } @testspos; # synonyms of 
		map {s/,phenyl alanine arylamidase,/,phenylalanine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,positive pyrrolidonyl arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,proline aminopeptidase,/,proline arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyrazinamidase activities,/,pyrazinamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyrazinamidase activity,/,pyrazinamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyroglutamate arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyroglutamic acid arylamidase activities,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyroglutamic acid arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyrrolidonyl arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,serine arylamidase activities,/,serine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,using commercial api systems arginine dihydrolase,/,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,valine aminopeptidase,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,valine arylamidase activities,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,valine arylamidase pyrrolidonyl arylamidase,/,valine arylamidase,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,valine arylamidases,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,valine arylamidasetrypsin activities,/,valine arylamidase,trypsin,/gi; } @testspos; # synonyms of 
		map {s/,valine arylaminidase,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,valine arylarnidase,/,valine arylamidase,/gi; } @testspos; # synonyms of 
		map {s/,weak pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testspos; # synonyms of 

		map {s/,coagulase activity,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,fixed coagulase,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,free coagulase,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,rapid coagulase test,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,test tube coagulase test reaction,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,test tube coagulase test,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,tube coagulase test,/,coagulase,/gi; } @testspos; # synonyms of 
		map {s/,tube coagulase,/,coagulase,/gi; } @testspos; # synonyms of 

		map {s/,benzidine,/,iron-phorphyrins,/gi; } @testspos; # synonyms of 
		map {s/,benzidine test,/,iron-phorphyrins,/gi; } @testspos; # synonyms of 
		map {s/,benzidine tests,/,iron-phorphyrins,/gi; } @testspos; # synonyms of 
		map {s/,benzidine test positive,/,iron-phorphyrins,/gi; } @testspos; # synonyms of 
		map {s/,benzidine test positive.thirty-one percent,/,iron-phorphyrins,/gi; } @testspos; # synonyms of 

		map {s/,naphthol phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphthol-as-bi-phosphate,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphthol-as-bi-phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphthol-as-bi-phosphohydrolase activities,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphthol-as-bi-phosphohydrolase activity,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphtol-as-bi-phosphohydrolase activities,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,naphtol-as-bi-phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,napthol-as-bi-phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,weak naphthol-as-bi-phosphohydrolase activity,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testspos; # synonyms of 

		map {s/,moderate phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 
		map {s/,phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testspos; # synonyms of 

		map {s/,indophenol oxidase,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,kovacs oxidase test,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,modified oxidase activity,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase activities,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase activity tests,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase activity,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase production,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase reaction,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase reactions,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase test,/,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,oxidase tests,/,oxidase,/gi; } @testspos; # synonyms of 

		map {s/,collagenase activity,/,collagenase,/gi; } @testspos; # synonyms of 

		map {s/,arylsulfatase production,/,arylsulfatase,/gi; } @testspos; # synonyms of 

		map {s/,weak catalase activity,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase activities,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase activity,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase production,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase reaction,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase reactions,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase test,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase tests,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase variable,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,catalase-positive and weakly oxidase,/,catalase,oxidase,/gi; } @testspos; # synonyms of 
		map {s/,feeble catalase activity,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,h2s catalase,/,h2s,catalase,/gi; } @testspos; # synonyms of 
		map {s/,moderate catalase activity,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,occasionally very weak catalase activity,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,positive catalase reaction,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,positive catalase,/,catalase,/gi; } @testspos; # synonyms of 
		map {s/,produce catalase,/,catalase,/gi; } @testspos; # synonyms of 

		map {s/,produces urease,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,strongest urease reaction,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urea hydrolysis,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease activities,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease activity,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease production,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease reaction,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease results,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease test,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,urease tests,/,urease,/gi; } @testspos; # synonyms of 

		map {s/,a-chymotrypsin activities,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,a-chymotrypsin,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,alpha-chymotrypsin,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,chymotrypsin activities,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,chymotrypsin activity,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,present and weak a-chymotrypsin activity,/,chymotrypsin,/gi; } @testspos; # synonyms of 
		map {s/,trypsin activities,/,trypsin,/gi; } @testspos; # synonyms of 

		map {s/,2-naphthyl-b-d-galactopyranoside,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,2-nitrophenyl-b-d-galactopyranoside,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-fucosidase activities,/,a-fucosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-fucosidase activity,/,a-fucosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-fucosidase production isnegative,/,not a-fucosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-glucosidase activities,/,a-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-glucosidase activity,/,a-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-mannosidase activity,/,a-mannosidase,/gi; } @testspos; # synonyms of 
		map {s/,a-mannosidase reaction,/,a-mannosidase,/gi; } @testspos; # synonyms of 
		map {s/,addition b-galactosidase,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,alpha-galactosidase,/,a-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,alpha-glucosidase,/,a-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,alpha-mannosidase,/,a-mannosidase,/gi; } @testspos; # synonyms of 
		map {s/,amylase plate test,/,amylase,/gi; } @testspos; # synonyms of 
		map {s/,amylase test,/,amylase,/gi; } @testspos; # synonyms of 
		map {s/,b-galactosidase 6-phosphatase,/,b-galactosidase 6-phosphate,/gi; } @testspos; # synonyms of 
		map {s/,b-galactosidase activities,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,b-galactosidase activity,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,b-galactosidase arginine,/,b-galactosidase,arginine,/gi; } @testspos; # synonyms of 
		map {s/,b-glucoronidase,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucosidase activities,/,b-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucosidase activity,/,b-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucosidase reactions,/,b-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucuronidase activities,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucuronidase activity,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,b-glucuronidase test,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,beta-galactosidase,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,beta-glucosidase,/,b-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,beta-glucuronidase,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,beta-mannosidase,/,b-mannosidase,/gi; } @testspos; # synonyms of 
		map {s/,biphosphate aldolase,/,fructose bisphosphate aldolase,/gi; } @testspos; # synonyms of 
		map {s/,diphosphate aldolase,/,fructose bisphosphate aldolase,/gi; } @testspos; # synonyms of 
		map {s/,enzymes b-glucosidase,/,b-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,exhibits a-glucosidase activity,/,a-glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,fucosidase activities,/,fucosidase,/gi; } @testspos; # synonyms of 
		map {s/,glucosidase activities,/,glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,glucosidase activity,/,glucosidase,/gi; } @testspos; # synonyms of 
		map {s/,hippurate hydrolysis,/,hippurate,/gi; } @testspos; # synonyms of 
		map {s/,hydrolyses xylan,/,xylanase,/gi; } @testspos; # synonyms of 
		map {s/,including xylanase,/,xylanase,/gi; } @testspos; # synonyms of 
		map {s/,liquefying amylase,/,amylase,/gi; } @testspos; # synonyms of 
		map {s/,methyl alpha-d-mannoside,/,methyl a-d-mannoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl b-d-glucopyranoside,/,methyl b-d-glucoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl beta-d-xyloside,/,methyl b-d-xyloside,/gi; } @testspos; # synonyms of 
		map {s/,methyl-a-d-glucopyranoside,/,methyl a-d-glucoside,/gi; } @testspos; # synonyms of 
		map {s/,methyl-a-d-glucoside,/,methyl a-d-glucoside,/gi; } @testspos; # synonyms of 
		map {s/,moderate b-galactosidase activity,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,n-acetyl-beta-D-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testspos; # synonyms of 
		map {s/,naphthol-as-bi-b-d-glucuronide,/,b-glucuronidase,/gi; } @testspos; # synonyms of 
		map {s/,o-nitrophenyl-b-d-galactoside,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,onpg hydrolysis,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,onpg reaction,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,onpg test,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,onpg,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,para-nitro-d-methyl galactose,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,producing b-galactosidase,/,b-galactosidase,/gi; } @testspos; # synonyms of 
		map {s/,saccharifying amylase,/,amylase,/gi; } @testspos; # synonyms of 
		map {s/,weak amylase production,/,amylase,/gi; } @testspos; # synonyms of 
		map {s/,xylan-inducible xylanase,/,xylanase,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testspos; # synonyms of 

		map {s/,deoxyribonuclease activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,deoxyribonuclease,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,dnase activities,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,dnase activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,moderate dnase activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,positive dnase activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,produce urease,/,urease,/gi; } @testspos; # synonyms of 
		map {s/,strong dnase,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,urease arginine dihydrolase,/,urease,arginine dihydrolase,/gi; } @testspos; # synonyms of 
		map {s/,weak deoxyribonuclease activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,weak dnase,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,weak extracellular deoxyribonuclease activity,/,dnase,/gi; } @testspos; # synonyms of 
		map {s/,weakly positive dnase,/,dnase,/gi; } @testspos; # synonyms of 

		map {s/,negative hyaluronidase activity,/,hyaluronidase,/gi; } @testspos; # synonyms of 

		map {s/,weak caseinase activity,/,caseinase,/gi; } @testspos; # synonyms of 
		map {s/,caseinase activities,/,caseinase,/gi; } @testspos; # synonyms of 

		map {s/,tetrazolium,/,tetrazolium reduction,/gi; } @testspos; # synonyms of 
		map {s/,tetrazolium red,/,tetrazolium reduction,/gi; } @testspos; # synonyms of 
		map {s/,tetrazolium chloride,/,tetrazolium reduction,/gi; } @testspos; # synonyms of 

		map {s/,lipase activities,/,lipase,/gi; } @testspos; # synonyms of 
		map {s/,lipase activity,/,lipase,/gi; } @testspos; # synonyms of 
		map {s/,lipase nor lecithinase,/,lipase,lecithinase,/gi; } @testspos; # synonyms of 
		map {s/,lipase reaction,/,lipase,/gi; } @testspos; # synonyms of 
		map {s/,lipolytic activity,/,lipase,/gi; } @testspos; # synonyms of 
		map {s/,strong and true lipolytic activity,/,lipase,/gi; } @testspos; # synonyms of 

		map {s/,nadp-dependent 6-phosphogluconate dehydrogenase,/,6-phosphogluconate dehydrogenase,/gi; } @testspos; # synonyms of 

	#metabolic products
		map {s/,methyl red,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red test,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red tests,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red indicator,/,methyl red acidification,/gi; } @testspos; # synonyms of 

		map {s/,acetoin production,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,acetoin,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,acetylmethylcarbinol production,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,acetylmethylcarbinol,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer reaction,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer reaction,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer reactions,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer reactions,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer test,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red and voges-proskauer test,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red reaction,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red reactions,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl red testing,/,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,methyl-red and voges-proskauer reactions,/,methyl red acidification,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,moderate acetylmethylcarbinol reactions,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,produces acetoin,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer and indole negative,/,acetoin,tryptophanase,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer and methylred reactions,/,acetoin,methyl red acidification,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer broth range,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer broth ranges,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer reaction,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer reactions,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer reagents,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer test range,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer test,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer tests,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,vp reaction,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,vp test,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,vp tests,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,vp,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,weak and voges-proskauer reaction,/,acetoin,/gi; } @testspos; # synonyms of 
		map {s/,weak vp reaction,/,acetoin,/gi; } @testspos; # synonyms of 

		map {s/,a-haemolysis,/,a-hemolysis,/gi; } @testspos; # synonyms of 
		map {s/,alpha hemolysis,/,a-hemolysis,/gi; } @testspos; # synonyms of 
		map {s/,alpha-haemolytic,/,a-hemolysis,/gi; } @testspos; # synonyms of 
		map {s/,alpha-hemolysis,/,a-hemolysis,/gi; } @testspos; # synonyms of 
		map {s/,b-hemolysis,/,b-haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,b-hemolysis,/,b-haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,beta-hemolytic,/,b-haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,blood haemolysis,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,double hemolysis zone,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,hemolysis activity,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,hemolysis,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,partial haemolysis,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,usually nonhemolytic,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,weak hemolysis activity,/,haemolysis,/gi; } @testspos; # synonyms of 
		map {s/,weak hemolysis,/,haemolysis,/gi; } @testspos; # synonyms of 

		map {s/,methylene blue,/,methylene blue reduction,/gi; } @testspos; # synonyms of 

		map {s/,koh test,/,flexirubin,/gi; } @testspos; # synonyms of 

		map {s/,indoxyl acetate reaction,/,indoxyl acetate hydrolysis,/gi; } @testspos; # synonyms of 

		map {s/,tetrathionate,/,tetrathionate reduction,/gi; } @testspos; # synonyms of 
		map {s/,resazurin,/,resazurin reduction,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testspos; # synonyms of 
#not a test
		map {s/,xxx,/,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,/gi; } @testspos; # synonyms of 
		map {s/,xxx,/,/gi; } @testspos; # synonyms of 
		map {s/,Gen,/,/gi; } @testspos; # synonyms of 
		map {s/,\d g yeast extract l-1or plus 10 mol glucose l-1,/,/gi; } @testspos; # synonyms of 
		map {s/,\d percent bile agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d-day buffered glucose peptone broth culture,/,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d ala,/,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d\d% agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d% agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d% nacl broth,/,/gi; } @testspos; # synonyms of 
		map {s/,\d.\d% nacl,/,/gi; } @testspos; # synonyms of 
		map {s/,\d\d% bile agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d% agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d% blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d% nacl,/,/gi; } @testspos; # synonyms of 
		map {s/,\d% sucrose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,\d% sucrose medium,/,/gi; } @testspos; # synonyms of 
		map {s/,0.36 mb agar,/,/gi; } @testspos; # synonyms of 
		map {s/,0.36 r2a agar,/,/gi; } @testspos; # synonyms of 
		map {s/,1 glu,/,/gi; } @testspos; # synonyms of 
		map {s/,1% peptone-1% yeast extract-1% glucose broth cultures,/,/gi; } @testspos; # synonyms of 
		map {s/,2,/,/gi; } @testspos; # synonyms of 
		map {s/,2216e agar,/,/gi; } @testspos; # synonyms of 
		map {s/,5% horse-blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,acetate containing media,/,/gi; } @testspos; # synonyms of 
		map {s/,actinomyces agar,/,/gi; } @testspos; # synonyms of 
		map {s/,aesculin bile agar,/,/gi; } @testspos; # synonyms of 
		map {s/,agar containing nacl,/,/gi; } @testspos; # synonyms of 
		map {s/,agar disk diffusion tests,/,/gi; } @testspos; # synonyms of 
		map {s/,agar media,/,/gi; } @testspos; # synonyms of 
		map {s/,agar medium containing bmc,/,/gi; } @testspos; # synonyms of 
		map {s/,agar medium,/,/gi; } @testspos; # synonyms of 
		map {s/,agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,agar roll tubes containing high numbers,/,/gi; } @testspos; # synonyms of 
		map {s/,agar roll tubes containing large numbers,/,/gi; } @testspos; # synonyms of 
		map {s/,agar shake tubes,/,/gi; } @testspos; # synonyms of 
		map {s/,agar slopes containing either cellobiose,/,/gi; } @testspos; # synonyms of 
		map {s/,ala,/,/gi; } @testspos; # synonyms of 
		map {s/,alkaline nutrient agar,/,/gi; } @testspos; # synonyms of 
		map {s/,alkaline soybean agar,/,/gi; } @testspos; # synonyms of 
		map {s/,amd,/,/gi; } @testspos; # synonyms of 
		map {s/,aminot acids glutamic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,and\/or eg agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,arthrobacter agar,/,/gi; } @testspos; # synonyms of 
		map {s/,azide-containing enterococcal selective agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bacitracin,/,/gi; } @testspos; # synonyms of 
		map {s/,basal sml0 containing glucose,/,/gi; } @testspos; # synonyms of 
		map {s/,betaine agar grown cells,/,/gi; } @testspos; # synonyms of 
		map {s/,bhi agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bhi-blood agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,bile . esculin agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bile aesculin agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bile agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bile.aesculin agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bile.esculin agar,/,/gi; } @testspos; # synonyms of 
		map {s/,blood agar base no. 2,/,/gi; } @testspos; # synonyms of 
		map {s/,blood agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,blood-containing agar,/,/gi; } @testspos; # synonyms of 
		map {s/,bovine blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,brain heart infusion agar,/,/gi; } @testspos; # synonyms of 
		map {s/,broth containing \d.\d% nacl,/,/gi; } @testspos; # synonyms of 
		map {s/,broth containing d-fructose,/,/gi; } @testspos; # synonyms of 
		map {s/,broth cultures containing glucose,/,/gi; } @testspos; # synonyms of 
		map {s/,c agar medium colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,c agar medium,/,/gi; } @testspos; # synonyms of 
		map {s/,casein agar,/,/gi; } @testspos; # synonyms of 
		map {s/,caso agar,/,/gi; } @testspos; # synonyms of 
		map {s/,catalase negative,/,/gi; } @testspos; # synonyms of 
		map {s/,cdc anaerobe 5% sheep blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,cell wall teichoic acid ontained glycerol,/,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose broth,/,/gi; } @testspos; # synonyms of 
		map {s/,cellobiose broths,/,/gi; } @testspos; # synonyms of 
		map {s/,cellulose agar medium colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,cellulose agar medium,/,/gi; } @testspos; # synonyms of 
		map {s/,cellulose agar surface colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,cellulose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,cetrimide agar,/,/gi; } @testspos; # synonyms of 
		map {s/,cgy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,chapman agar,/,/gi; } @testspos; # synonyms of 
		map {s/,citrate agar,/,/gi; } @testspos; # synonyms of 
		map {s/,ckonkey agar,/,/gi; } @testspos; # synonyms of 
		map {s/,columbia agar,/,/gi; } @testspos; # synonyms of 
		map {s/,columbia blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,columbia sheep-blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,complex standard media acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,containing 0.006% phenol red,/,/gi; } @testspos; # synonyms of 
		map {s/,containing 5% xylose,/,/gi; } @testspos; # synonyms of 
		map {s/,cyc agar,/,/gi; } @testspos; # synonyms of 
		map {s/,e.g. blood trypticase soy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,eg agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,eg agar,/,/gi; } @testspos; # synonyms of 
		map {s/,egg.yolk agar,/,/gi; } @testspos; # synonyms of 
		map {s/,erythromycin,/,/gi; } @testspos; # synonyms of 
		map {s/,facultatively aerobic and heterofermentative lactic acid bacterium,/,/gi; } @testspos; # synonyms of 
		map {s/,ferment .*,/,/gi; } @testspos; # synonyms of 
		map {s/,fermentative growth using pyruvate,/,/gi; } @testspos; # synonyms of 
		map {s/,fermented.*,/,/gi; } @testspos; # synonyms of 
		map {s/,fermenting.*,/,/gi; } @testspos; # synonyms of 
		map {s/,ferments.*,/,/gi; } @testspos; # synonyms of 
		map {s/,gel electrophoretic patterns,/,/gi; } @testspos; # synonyms of 
		map {s/,glu,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose 0.85 butyrate.*,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose agar deep cultures,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose broth,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose broths,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose media,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose medium,/,/gi; } @testspos; # synonyms of 
		map {s/,glucose peptone-yeast extract broth,/,/gi; } @testspos; # synonyms of 
		map {s/,glutamate media,/,/gi; } @testspos; # synonyms of 
		map {s/,glycerol media,/,/gi; } @testspos; # synonyms of 
		map {s/,glycerol yields acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,gyp agar colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,gyp agar slant,/,/gi; } @testspos; # synonyms of 
		map {s/,gyp agar,/,/gi; } @testspos; # synonyms of 
		map {s/,h2.co2 plus acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,h2.co2.1 mm acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,hemin,/,/gi; } @testspos; # synonyms of 
		map {s/,horse-blood trypticase soy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,horse.blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,horseblood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,iron agar,/,/gi; } @testspos; # synonyms of 
		map {s/,isp 2 agar,/,/gi; } @testspos; # synonyms of 
		map {s/,kaa agar,/,/gi; } @testspos; # synonyms of 
		map {s/,kanamycin aesculin azide agar,/,/gi; } @testspos; # synonyms of 
		map {s/,kanamycin,/,/gi; } @testspos; # synonyms of 
		map {s/,kanamycin.aesculin.azide agar,/,/gi; } @testspos; # synonyms of 
		map {s/,kf streptococcus agar,/,/gi; } @testspos; # synonyms of 
		map {s/,l-rhamnose broth,/,/gi; } @testspos; # synonyms of 
		map {s/,l-tyrosine-containing marine agar,/,/gi; } @testspos; # synonyms of 
		map {s/,lab lemco agar,/,/gi; } @testspos; # synonyms of 
		map {s/,lacking hemolysis,/,/gi; } @testspos; # synonyms of 
		map {s/,lb agar,/,/gi; } @testspos; # synonyms of 
		map {s/,le ruyet et al.,/,/gi; } @testspos; # synonyms of 
		map {s/,lead acetate agar,/,/gi; } @testspos; # synonyms of 
		map {s/,lead acetate paper,/,/gi; } @testspos; # synonyms of 
		map {s/,least 1.5% agar,/,/gi; } @testspos; # synonyms of 
		map {s/,lincomycin,/,/gi; } @testspos; # synonyms of 
		map {s/,lipid films,/,/gi; } @testspos; # synonyms of 
		map {s/,litmus milk,/,/gi; } @testspos; # synonyms of 
		map {s/,m.l. gen n. bromii,/,/gi; } @testspos; # synonyms of 
		map {s/,M.L. gen n. bromii,/,/gi; } @testspos; # synonyms of 
		map {s/,m.r.s. agar,/,/gi; } @testspos; # synonyms of 
		map {s/,m17 glucose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,macconkey agar,/,/gi; } @testspos; # synonyms of 
		map {s/,maltose amygdalin,/,/gi; } @testspos; # synonyms of 
		map {s/,maltose broth,/,/gi; } @testspos; # synonyms of 
		map {s/,mannose containing lipid,/,/gi; } @testspos; # synonyms of 
		map {s/,marine agar 2216 media,/,/gi; } @testspos; # synonyms of 
		map {s/,marine agar 2216,/,/gi; } @testspos; # synonyms of 
		map {s/,marine agar 2218,/,/gi; } @testspos; # synonyms of 
		map {s/,meat agar,/,/gi; } @testspos; # synonyms of 
		map {s/,meat broth agar,/,/gi; } @testspos; # synonyms of 
		map {s/,media containing \d% glucose,/,/gi; } @testspos; # synonyms of 
		map {s/,media containing nitrate,/,/gi; } @testspos; # synonyms of 
		map {s/,medium containing glucose plus yeast extract,/,/gi; } @testspos; # synonyms of 
		map {s/,mineral media containing pyruvate,/,/gi; } @testspos; # synonyms of 
		map {s/,mineral medium containing 10 mm acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,minocycline,/,/gi; } @testspos; # synonyms of 
		map {s/,mm10 blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,mm10 sheep blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,mrs agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,mrs agar,/,/gi; } @testspos; # synonyms of 
		map {s/,msg agar,/,/gi; } @testspos; # synonyms of 
		map {s/,mty agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,na containing 4% nacl,/,/gi; } @testspos; # synonyms of 
		map {s/,nacl agar,/,/gi; } @testspos; # synonyms of 
		map {s/,nalidixic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,negative aesculin hydrolysis,/,/gi; } @testspos; # synonyms of 
		map {s/,negative and aesculin hydrolysis,/,/gi; } @testspos; # synonyms of 
		map {s/,negative and methyl red test,/,/gi; } @testspos; # synonyms of 
		map {s/,negative nitrate reduction,/,/gi; } @testspos; # synonyms of 
		map {s/,negative voges-proskauer reaction,/,/gi; } @testspos; # synonyms of 
		map {s/,nitrate reduction-negative,/,/gi; } @testspos; # synonyms of 
		map {s/,non-haemolytic,/,/gi; } @testspos; # synonyms of 
		map {s/,non-hemolytic,/,/gi; } @testspos; # synonyms of 
		map {s/,nonhemolytic,/,/gi; } @testspos; # synonyms of 
		map {s/,novobiocin,/,/gi; } @testspos; # synonyms of 
		map {s/,nutrient agar colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,nutrient agar cultures,/,/gi; } @testspos; # synonyms of 
		map {s/,nutrient agar growth,/,/gi; } @testspos; # synonyms of 
		map {s/,nutrient agar,/,/gi; } @testspos; # synonyms of 
		map {s/,nutrient agar,/,/gi; } @testspos; # synonyms of 
		map {s/,old culture,/,/gi; } @testspos; # synonyms of 
		map {s/,old cultures,/,/gi; } @testspos; # synonyms of 
		map {s/,old,/,/gi; } @testspos; # synonyms of 
		map {s/,ox gall,/,/gi; } @testspos; # synonyms of 
		map {s/,oxidase negative,/,/gi; } @testspos; # synonyms of 
		map {s/,p agar,/,/gi; } @testspos; # synonyms of 
		map {s/,parasites-of man,/,/gi; } @testspos; # synonyms of 
		map {s/,pc skim milk agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,pea azide agar,/,/gi; } @testspos; # synonyms of 
		map {s/,peptone yeast glucose broth,/,/gi; } @testspos; # synonyms of 
		map {s/,peptone-yeast broth medium containing glucose,/,/gi; } @testspos; # synonyms of 
		map {s/,peptone-yeast extract-glucose deep agar cultures,/,/gi; } @testspos; # synonyms of 
		map {s/,ph in voges-proskauer broth,/,/gi; } @testspos; # synonyms of 
		map {s/,phosphatidyl ethanolamine,/,/gi; } @testspos; # synonyms of 
		map {s/,phosphatidyl glycerol,/,/gi; } @testspos; # synonyms of 
		map {s/,positive and .*,/,/gi; } @testspos; # synonyms of 
		map {s/,potato agar,/,/gi; } @testspos; # synonyms of 
		map {s/,potato dextrose agar,/,/gi; } @testspos; # synonyms of 
		map {s/,potato-glucose agar media,/,/gi; } @testspos; # synonyms of 
		map {s/,ppyg agar,/,/gi; } @testspos; # synonyms of 
		map {s/,pras media+ 0.5% tween 80,/,/gi; } @testspos; # synonyms of 
		map {s/,predominantly l-lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,predominately l-lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,primarily lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,producing acetate,/,/gi; } @testspos; # synonyms of 
		map {s/,producing exclusively lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,producing l-lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,producing lactic acid,/,/gi; } @testspos; # synonyms of 
		map {s/,py agar,/,/gi; } @testspos; # synonyms of 
		map {s/,py-1% glucose medium ph 6.8,/,/gi; } @testspos; # synonyms of 
		map {s/,pye agar,/,/gi; } @testspos; # synonyms of 
		map {s/,pyi agar,/,/gi; } @testspos; # synonyms of 
		map {s/,pyt broth medium containing esculin,/,/gi; } @testspos; # synonyms of 
		map {s/,r2a agar,/,/gi; } @testspos; # synonyms of 
		map {s/,raffinose liquid media,/,/gi; } @testspos; # synonyms of 
		map {s/,rgm-rumen fluid plus glucose,/,/gi; } @testspos; # synonyms of 
		map {s/,salmonella-shigella agar,/,/gi; } @testspos; # synonyms of 
		map {s/,schaedler sheep blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,semi-solid agar medium,/,/gi; } @testspos; # synonyms of 
		map {s/,sheep blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,sheep-blood agar,/,/gi; } @testspos; # synonyms of 
		map {s/,skimmed-milk agar,/,/gi; } @testspos; # synonyms of 
		map {s/,slanetz & bartley agar,/,/gi; } @testspos; # synonyms of 
		map {s/,slanetz-bartley agar,/,/gi; } @testspos; # synonyms of 
		map {s/,soft 0.5% mrs agar,/,/gi; } @testspos; # synonyms of 
		map {s/,solid agar containing cellulose,/,/gi; } @testspos; # synonyms of 
		map {s/,spot,/,/gi; } @testspos; # synonyms of 
		map {s/,standard mrs agar,/,/gi; } @testspos; # synonyms of 
		map {s/,swcm agar,/,/gi; } @testspos; # synonyms of 
		map {s/,tgy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,thallous acetate-tetrazolium agar,/,/gi; } @testspos; # synonyms of 
		map {s/,thermophilic and hydrolyses agar,/,/gi; } @testspos; # synonyms of 
		map {s/,thin agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,todd-hewitt agar,/,/gi; } @testspos; # synonyms of 
		map {s/,troleandomycin,/,/gi; } @testspos; # synonyms of 
		map {s/,trypticase soy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,trypticase soy broth agar,/,/gi; } @testspos; # synonyms of 
		map {s/,trypticase soy sheep blood agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,trypticase-soy-sheep blood agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,tryptone soy agar,/,/gi; } @testspos; # synonyms of 
		map {s/,tryptone soya agar,/,/gi; } @testspos; # synonyms of 
		map {s/,tryptose soya agar colonies,/,/gi; } @testspos; # synonyms of 
		map {s/,tsa agar,/,/gi; } @testspos; # synonyms of 
		map {s/,tsb agar,/,/gi; } @testspos; # synonyms of 
		map {s/,tween 80 agar,/,/gi; } @testspos; # synonyms of 
		map {s/,vibriostatic agent o.129,/,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer broth,/,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer medium,/,/gi; } @testspos; # synonyms of 
		map {s/,voges-proskauer negative,/,/gi; } @testspos; # synonyms of 
		map {s/,vy\/2 agar,/,/gi; } @testspos; # synonyms of 
		map {s/,wet agar plates,/,/gi; } @testspos; # synonyms of 
		map {s/,wilkins agar,/,/gi; } @testspos; # synonyms of 
		map {s/,xr2a agarxx,/,/gi; } @testspos; # synonyms of 
		map {s/,xylan agar medium,/,/gi; } @testspos; # synonyms of 
		map {s/,xylan agar,/,/gi; } @testspos; # synonyms of 
		map {s/,ycfa medium containing 0.2% glucose,/,/gi; } @testspos; # synonyms of 

		print OUT @taxlabels, @testspos, "\n"; # prints to $homout, hom.testspos.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.testspos.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/tests positive//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.testspos.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.testspos.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.testspos.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.testspos.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.testspos.txt
#temporarily rename charstates to label those that have been homologized		
	#substrates used and/or fermented
	map {s/^2-aminoethanol/**2-aminoethanol/gi } @charlist2;#1
	map {s/^2-keto-d-gluconate/**2-keto-d-gluconate/gi } @charlist2;#2
	map {s/^2-ketobutyrate/**2-ketobutyrate/gi } @charlist2;#3
	map {s/^2-3-butanediol/**2-3-butanediol/gi } @charlist2;#4
	map {s/^1-3-butanediol/**1-3-butanediol/gi } @charlist2;#5
	map {s/^3-hydroxybenzoate/**3-hydroxybenzoate/gi } @charlist2;#6
	map {s/^3-hydroxybutyrate/**3-hydroxybutyrate/gi } @charlist2;#7
	map {s/^3-methyl d-glucose/**3-methyl d-glucose/gi } @charlist2;#8
	map {s/^4-hydroxybenzoate/**4-hydroxybenzoate/gi } @charlist2;#9
	map {s/^4-hydroxyphenylacetate/**4-hydroxyphenylacetate/gi } @charlist2;#10
	map {s/^5-keto-d-gluconate/**5-keto-d-gluconate/gi } @charlist2;#12
	map {s/^a-ketoglutarate/**a-ketoglutarate/gi } @charlist2;#13
	map {s/^a-ketovalerate/**a-ketovalerate/gi } @charlist2;#14
	map {s/^acetate/**acetate/gi } @charlist2;#15
	map {s/^acetoacetate/**acetoacetate/gi } @charlist2;#16
	map {s/^adenosine 5\'-monophosphate/**adenosine 5\'-monophosphate/gi } @charlist2;#17
	map {s/^adenosine/**adenosine/gi } @charlist2;#18
	map {s/^adipate/**adipate/gi } @charlist2;#19
	map {s/^adonitol/**adonitol/gi } @charlist2;#20
	map {s/^agar/**agar/gi } @charlist2;#21
	map {s/^alaninamide/**alaninamide/gi } @charlist2;#22
	map {s/^alanine/**alanine/gi } @charlist2;#23
	map {s/^alanyl glycine/**alanyl glycine/gi } @charlist2;#24
	map {s/^amygdalin/**amygdalin/gi } @charlist2;#25
	map {s/^arabinose/**arabinose/gi } @charlist2;#26
	map {s/^arabitol/**arabitol/gi } @charlist2;#27
	map {s/^arbutin/**arbutin/gi } @charlist2;#28
	map {s/^arginine/**arginine/gi } @charlist2;#29
	map {s/^asparagine/**asparagine/gi } @charlist2;#30
	map {s/^aspartate/**aspartate/gi } @charlist2;#31
	map {s/^b-cyclodextrin/**b-cyclodextrin/gi } @charlist2;#32
	map {s/^b-xylosidase/**b-xylosidase/gi } @charlist2;#34
	map {s/^butyrate/**butyrate/gi } @charlist2;#35
	map {s/^caprate/**caprate/gi } @charlist2;#36
	map {s/^carbohydrates/**carbohydrates/gi } @charlist2;#37
	map {s/^carboxymethyl cellulose/**carboxymethyl cellulose/gi } @charlist2;#38
	map {s/^carnitine/**carnitine/gi } @charlist2;#39
	map {s/^cellobiose/**cellobiose/gi } @charlist2;#40
	map {s/^cis-aconitate/**cis-aconitate/gi } @charlist2;#41
	map {s/^citrate/**citrate/gi } @charlist2;#42
	map {s/^co2/**co2/gi } @charlist2;#43
	map {s/^cyclodextrin/**cyclodextrin/gi } @charlist2;#44
	map {s/^dextrin/**dextrin/gi } @charlist2;#45
	map {s/^dulcitol/**dulcitol/gi } @charlist2;#46
	map {s/^erythritol/**erythritol/gi } @charlist2;#47
	map {s/^ethanol/**ethanol/gi } @charlist2;#48
	map {s/^formate/**formate/gi } @charlist2;#49
	map {s/^fructose 6-phosphate/**fructose 6-phosphate/gi } @charlist2;#50
	map {s/^fructose/**fructose/gi } @charlist2;#51
	map {s/^fucose/**fucose/gi } @charlist2;#52
	map {s/^fumarate/**fumarate/gi } @charlist2;#53
	map {s/^fusidic acid/**fusidic acid/gi } @charlist2;#54
	map {s/^galactose/**galactose/gi } @charlist2;#55
	map {s/^galacturonate/**galacturonate/gi } @charlist2;#56
	map {s/^gelatin/**gelatin/gi } @charlist2;#57
	map {s/^gentiobiose/**gentiobiose/gi } @charlist2;#58
	map {s/^gluconate/**gluconate/gi } @charlist2;#59
	map {s/^glucose 1-phosphate/**glucose 1-phosphate/gi } @charlist2;#60
	map {s/^glucose/**glucose/gi } @charlist2;#61
	map {s/^glucuronate/**glucuronate/gi } @charlist2;#62
	map {s/^glutamate/**glutamate/gi } @charlist2;#63
	map {s/^glutamic acid polymer/**glutamic acid polymer/gi } @charlist2;#64
	map {s/^glutamine/**glutamine/gi } @charlist2;#65
	map {s/^glutamyl glutamic acid/**glutamyl glutamic acid/gi } @charlist2;#66
	map {s/^glycerol phosphate/**glycerol phosphate/gi } @charlist2;#68
	map {s/^glycerol/**glycerol/gi } @charlist2;#69
	map {s/^glycogen/**glycogen/gi } @charlist2;#70
	map {s/^glycyl l-glutamic acid/**glycyl l-glutamic acid/gi } @charlist2;#71
	map {s/^glycyl-l-proline/**glycyl-l-proline/gi } @charlist2;#72
	map {s/^glyoxylate/**glyoxylate/gi } @charlist2;#73
	map {s/^guanidine/**guanidine/gi } @charlist2;#74
	map {s/^h2s/**h2s/gi } @charlist2;#75
	map {s/^hippurate/**hippurate/gi } @charlist2;#76
	map {s/^histidine/**histidine/gi } @charlist2;#77
	map {s/^indole/**indole/gi } @charlist2;#78
	map {s/^inosine/**inosine/gi } @charlist2;#79
	map {s/^inositol/**inositol/gi } @charlist2;#80
	map {s/^inulin/**inulin/gi } @charlist2;#81
	map {s/^itaconate/**itaconate/gi } @charlist2;#82
	map {s/^lactate/**lactate/gi } @charlist2;#83
	map {s/^lactic acid methyl ester/**lactic acid methyl ester/gi } @charlist2;#84
	map {s/^lactose/**lactose/gi } @charlist2;#85
	map {s/^lactulose/**lactulose/gi } @charlist2;#86
	map {s/^leucine/**leucine/gi } @charlist2;#87
	map {s/^lyxose/**lyxose/gi } @charlist2;#88
	map {s/^malate/**malate/gi } @charlist2;#89
	map {s/^malonate/**malonate/gi } @charlist2;#90
	map {s/^maltose/**maltose/gi } @charlist2;#92
	map {s/^maltotriose/**maltotiose/gi } @charlist2;#93
	map {s/^mannan/**mannan/gi } @charlist2;#94
	map {s/^mannitol/**mannitol/gi } @charlist2;#95
	map {s/^mannose/**mannose/gi } @charlist2;#96
	map {s/^melezitose/**melezitose/gi } @charlist2;#97
	map {s/^melibiose/**melibiose/gi } @charlist2;#98
	map {s/^meso-tartrate/**meso-tartrate/gi } @charlist2;#99
	map {s/^methionine/**methionine/gi } @charlist2;#100
	map {s/^methyl a-d-galactoside/**methyl a-d-galactoside/gi } @charlist2;#101
	map {s/^methyl a-d-glucoside/**methyl a-d-glucoside/gi } @charlist2;#102
	map {s/^methyl b-d-galactoside/**methyl b-d-galactoside/gi } @charlist2;#103
	map {s/^methyl b-d-glucoside/**methyl b-d-glucoside/gi } @charlist2;#104
	map {s/^methyl pyruvate/**methyl pyruvate/gi } @charlist2;#105
	map {s/^monomethyl succinate/**monomethyl succinate/gi } @charlist2;#106
	map {s/^mucate/**mucate/gi } @charlist2;#107
	map {s/^myo-inositol/**myo-inositol/gi } @charlist2;#108
	map {s/^n-acetyl-b-d-glucosamine/**n-acetyl-b-d-glucosamine/gi } @charlist2;#109
	map {s/^n-acetyl-b-d-mannosamine/**n-acetyl-b-d-mannosamine/gi } @charlist2;#110
	map {s/^n-acetyl-l-glutamate/**n-acetyl-l-glutamate/gi } @charlist2;#111
	map {s/^ornithine/**ornithine/gi } @charlist2;#112
	map {s/^palatinose/**palatinose/gi } @charlist2;#113
	map {s/^pectin/**pectin/gi } @charlist2;#114
	map {s/^phenylacetate/**phenylacetate/gi } @charlist2;#115
	map {s/^phenylalanine/**phenylalanine/gi } @charlist2;#116
	map {s/^phenylethylamine/**phenylethylamine/gi } @charlist2;#117
	map {s/^proline/**proline/gi } @charlist2;#118
	map {s/^propionate/**propionate/gi } @charlist2;#119
	map {s/^psicose/**psicose/gi } @charlist2;#120
	map {s/^pullulan/**pullulan/gi } @charlist2;#121
	map {s/^putrescine/**putrescine/gi } @charlist2;#122
	map {s/^pyroglutamate/**pyroglutamate/gi } @charlist2;#123
	map {s/^pyruvate/**pyruvate/gi } @charlist2;#124
	map {s/^quinate/**quinate/gi } @charlist2;#126
	map {s/^raffinose/**raffinose/gi } @charlist2;#127
	map {s/^rhamnose/**rhamnose/gi } @charlist2;#128
	map {s/^ribitol/**ribitol/gi } @charlist2;#129
	map {s/^ribose/**ribose/gi } @charlist2;#130
	map {s/^saccharolytic/**saccharolytic/gi } @charlist2;#131
	map {s/^salicin/**salicin/gi } @charlist2;#133
	map {s/^sedoheptulosan/**sedoheptulosan/gi } @charlist2;#134
	map {s/^serine/**serine/gi } @charlist2;#135
	map {s/^sorbitol/**sorbitol/gi } @charlist2;#136
	map {s/^sorbose/**sorbose/gi } @charlist2;#137
	map {s/^stachyose/**stachyose/gi } @charlist2;#138
	map {s/^starch/**starch/gi } @charlist2;#139
	map {s/^suberate/**suberate/gi } @charlist2;#140
	map {s/^succinamic acid/**succinamic acid/gi } @charlist2;#141
	map {s/^succinate/**succinate/gi } @charlist2;#142
	map {s/^cysteine/**cysteine/gi } @charlist2;#143
	map {s/^sucrose/**sucrose/gi } @charlist2;#144
	map {s/^tagatose/**tagatose/gi } @charlist2;#145
	map {s/^tartrate/**tartrate/gi } @charlist2;#146
	map {s/^tellurite/**tellurite/gi } @charlist2;#147
	map {s/^threonine/**threonine/gi } @charlist2;#148
	map {s/^thymidine 5\'-monophosphate/**thymidine 5\'-monophosphate/gi } @charlist2;#149
	map {s/^thymidine/**thymidine/gi } @charlist2;#150
	map {s/^trehalose/**trehalose/gi } @charlist2;#151
	map {s/^glycine/**glycine/gi } @charlist2;#152
	map {s/^turanose/**turanose/gi } @charlist2;#153
	map {s/^tween 40/**tween 40/gi } @charlist2;#154
	map {s/^tween 80/**tween 80/gi } @charlist2;#155
	map {s/^uridine/**uridine/gi } @charlist2;#156
	map {s/^valerate/**valerate/gi } @charlist2;#157
	map {s/^valine/**valine/gi } @charlist2;#158
	map {s/^xylan/**xylan/gi } @charlist2;#159
	map {s/^xylitol/**xylitol/gi } @charlist2;#160
	map {s/^xylose/**xylose/gi } @charlist2;#161

	map {s/^isoleucine/**isoleucine/gi } @charlist2;#162
	map {s/^lysine/**lysine/gi } @charlist2;#163
	map {s/^tryptophan/**tryptophan/gi } @charlist2;#164
	map {s/^selenocysteine/**selenocysteine/gi } @charlist2;#165
	map {s/^pyrrolysine/**pyrrolysine/gi } @charlist2;#166

	#inorganics
	map {s/^nitrate reduction/**nitrate reduction/gi } @charlist2;#167
	map {s/^nitrogen fixation/**nitrogen fixation/gi } @charlist2;#168
	map {s/^nitrite reduction/**nitrite reduction/gi } @charlist2;#169

	#enzymes present or enzymatic substrates
	map {s/^6-phospho-b-galactosidase/**6-phospho-b-galactosidase/gi } @charlist2;#170
	map {s/^6-phosphogluconate dehydrogenase/**6-phosphogluconate dehydrogenase/gi } @charlist2;#171
	map {s/^a-arabinosidase/**a-arabinosidase/gi } @charlist2;#172
	map {s/^a-fucosidase/**a-fucosidase/gi } @charlist2;#173
	map {s/^a-galactosidase/**a-galactosidase/gi } @charlist2;#174
	map {s/^a-glucosidase/**a-glucosidase/gi } @charlist2;#175
	map {s/^a-mannosidase/**a-mannosidase/gi } @charlist2;#176
	map {s/^a-methyl glucosidase/**a-methyl glucosidase/gi } @charlist2;#177
	map {s/^a-methylmannoside/**a-methylmannoside/gi } @charlist2;#178
	map {s/^acid phosphatase/**acid phosphatase/gi } @charlist2;#179
	map {s/^aesculin hydrolysis/**aesculin hydrolysis/gi } @charlist2;#180
	map {s/^alanine arylamidase/**alanine arylamidase/gi } @charlist2;#181
	map {s/^alkaline phosphatase/**alkaline phosphatase/gi } @charlist2;#182
	map {s/^alpha-glutamate arylamidase/**alpha-glutamate arylamidase/gi } @charlist2;#183
	map {s/^amylase/**amylase/gi } @charlist2;#184
	map {s/^arginine arylamidase/**arginine arylamidase/gi } @charlist2;#185
	map {s/^arginine decarboxylase/**arginine decarboxylase/gi } @charlist2;#186
	map {s/^arginine dihydrolase/**arginine dihydrolase/gi } @charlist2;#187
	map {s/^arylsulfatase/**arylsulfatase/gi } @charlist2;#188
	map {s/^asparagine arylamidase/**asparagine arylamidase/gi } @charlist2;#189
	map {s/^aspartate arylamidase/**aspartate arylamidase/gi } @charlist2;#190
	map {s/^b-galactosidase/**b-galactosidase/gi } @charlist2;#191
	map {s/^b-glucosidase/**b-glucosidase/gi } @charlist2;#192
	map {s/^b-glucuronidase/**b-glucuronidase/gi } @charlist2;#193
	map {s/^b-mannosidase/**b-mannosidase/gi } @charlist2;#194
	map {s/^casein hydrolysis/**casein hydrolysis/gi } @charlist2;#195
	map {s/^caseinase/**caseinase/gi } @charlist2;#196
	map {s/^catalase/**catalase/gi } @charlist2;#197
	map {s/^cellobiose/**cellobiose/gi } @charlist2;#198
	map {s/^chymotrypsin/**chymotrypsin/gi } @charlist2;#199
	map {s/^coagulase/**coagulase/gi } @charlist2;#200
	map {s/^collagenase/**collagenase/gi } @charlist2;#201
	map {s/^cysteine arylamidase/**cysteine arylamidase/gi } @charlist2;#202
	map {s/^cystine arylamidase/**cystine arylamidase/gi } @charlist2;#203
	map {s/^cytochrome oxidase/**cytochrome oxidase/gi } @charlist2;#204
	map {s/^dnase/**dnase/gi } @charlist2;#205
	map {s/^esterase c4/**esterase c4/gi } @charlist2;#206
	map {s/^fructose bisphosphate aldolase/**fructose bisphosphate aldolase/gi } @charlist2;#207
	map {s/^fucosidase/**fucosidase/gi } @charlist2;#208
	map {s/^gelatinase/**gelatinase/gi } @charlist2;#209
	map {s/^glucosidase/**glucosidase/gi } @charlist2;#210
	map {s/^glutamate arylamidase/**glutamate arylamidase/gi } @charlist2;#211
	map {s/^glutamate deaminase/**glutamate deaminase/gi } @charlist2;#212
	map {s/^glutamic acid decarboxylase/**glutamic acid decarboxylase/gi } @charlist2;#213
	map {s/^glutamine arylamidase/**glutamine arylamidase/gi } @charlist2;#214
	map {s/^glutamyl glutamic acid arylamidase/**glutamyl glutamic acid arylamidase/gi } @charlist2;#215
	map {s/^glycine arylamidase/**glycine arylamidase/gi } @charlist2;#216
	map {s/^hippurate hydrolysis/**hippurate hydrolysis/gi } @charlist2;#217
	map {s/^histidine arylamidase/**histidine arylamidase/gi } @charlist2;#218
	map {s/^hyaluronidase/**hyaluronidase/gi } @charlist2;#219
	map {s/^hydroxy proline arylamidase/**hydroxy proline arylamidase/gi } @charlist2;#220
	map {s/^iron-phorphyrins/**iron-phorphyrins/gi } @charlist2;#221
	map {s/^isoleucine arylamidase/**isoleucine arylamidase/gi } @charlist2;#222
	map {s/^l-lactic acid racemase/**l-lactic acid racemase/gi } @charlist2;#223
	map {s/^lecithinase/**lecithinase/gi } @charlist2;#224
	map {s/^leucine arylamidase/**leucine arylamidase/gi } @charlist2;#225
	map {s/^leucyl glycine arylamidase/**leucyl glycine arylamidase/gi } @charlist2;#226
	map {s/^lipase c14/**lipase c14/gi } @charlist2;#227
	map {s/^lipase c8/**lipase c8/gi } @charlist2;#228
	map {s/^lipase/**lipase/gi } @charlist2;#229
	map {s/^lysine arylamidase/**lysine arylamidase/gi } @charlist2;#230
	map {s/^lysine deaminase/**lysine deaminase/gi } @charlist2;#231
	map {s/^lysine decarboxylase/**lysine decarboxylase/gi } @charlist2;#232
	map {s/^methionine arylamidase/**methionine arylamidase/gi } @charlist2;#233
	map {s/^methyl a-d-mannoside/**methyl a-d-mannoside/gi } @charlist2;#234
	map {s/^methyl b-d-xyloside/**methyl b-d-xyloside/gi } @charlist2;#235
	map {s/^n-acetyl-b-d-glucosaminidase/**n-acetyl-b-d-glucosaminidase/gi } @charlist2;#236
	map {s/^nadh oxidase/**nadh oxidase/gi } @charlist2;#237
	map {s/^naphthol-as-bi-phosphohydrolase/**naphthol-as-bi-phosphohydrolase/gi } @charlist2;#238
	map {s/^ornithine arylamidase/**ornithine arylamidase/gi } @charlist2;#239
	map {s/^ornithine decarboxylase/**ornithine decarboxylase/gi } @charlist2;#240
	map {s/^oxidase/**oxidase/gi } @charlist2;#241
	map {s/^para-phenylalanine deaminase/**para-phenylalanine deaminase/gi } @charlist2;#242
	map {s/^phenylalanine arylamidase/**phenylalanine arylamidase/gi } @charlist2;#243
	map {s/^phenylalanine deaminase/**phenylalanine deaminase/gi } @charlist2;#244
	map {s/^phosphatase/**phosphatase/gi } @charlist2;#245
	map {s/^proline arylamidase/**proline arylamidase/gi } @charlist2;#246
	map {s/^protease/**protease/gi } @charlist2;#247
	map {s/^pyrazinamidase/**pyrazinamidase/gi } @charlist2;#248
	map {s/^pyroglutamic acid arylamidase/**pyroglutamic acid arylamidase/gi } @charlist2;#249
	map {s/^serine arylamidase/**serine arylamidase/gi } @charlist2;#251
	map {s/^threonine arylamidase/**threonine arylamidase/gi } @charlist2;#252
	map {s/^trypsin/**trypsin/gi } @charlist2;#253
	map {s/^tryptophan arylamidase/**tryptophan arylamidase/gi } @charlist2;#254
	map {s/^tryptophan deaminase/**tryptophan deaminase/gi } @charlist2;#255
	map {s/^tryptophanase/**tryptophanase/gi } @charlist2;#256
	map {s/^tyrosine arylamidase/**tyrosine arylamidase/gi } @charlist2;#257
	map {s/^urease/**urease/gi } @charlist2;#258
	map {s/^valine arylamidase/**valine arylamidase/gi } @charlist2;#259
	map {s/^xylanase/**xylanase/gi } @charlist2;#260

	#metabolic products
	map {s/^a-hemolysis/**a-hemolysis/gi } @charlist2;#261
	map {s/^acetoin/**acetoin/gi } @charlist2;#262
	map {s/^b-haemolysis/**b-haemolysis/gi } @charlist2;#263
	map {s/^ellman/**ellman/gi } @charlist2;#264
	map {s/^flexirubin/**flexirubin/gi } @charlist2;#265
	map {s/^haemolysis/**haemolysis/gi } @charlist2;#266
	map {s/^indoxyl acetate hydrolysis/**indoxyl acetate hydrolysis/gi } @charlist2;#267
	map {s/^methyl red acidification/**methyl red acidification/gi } @charlist2;#268
	map {s/^methylene blue reduction/**methylene blue reduction/gi } @charlist2;#269
	map {s/^resazurin reduction/**resazurin reduction/gi } @charlist2;#270
	map {s/^tetrathionate reduction/**tetrathionate reduction/gi } @charlist2;#271
	map {s/^tetrazolium reduction/**tetrazolium reduction/gi } @charlist2;#272
	map {s/^carboxymethyl cellulase/**carboxymethyl cellulase/gi } @charlist2;#273

	map {s/^1-3-propanediol/**1-3-propanediol/gi } @charlist2; #274
	map {s/^1-2-propanediol/**1-2-propanediol/gi } @charlist2; #275
	map {s/^propanediol/**propanediol/gi } @charlist2; #276

	map {s/^phosphatidylinositol phospholipase c/**phosphatidylinositol phospholipase c/gi } @charlist2;#277
	map {s/^lecithinase/**lecithinase/gi } @charlist2;#278
	map {s/^lactamide/**lactamide/gi } @charlist2;#279
	map {s/^hydroxyproline/**hydroxyproline/gi } @charlist2;#280
	map {s/^glycyl tryptophan arylamidase/**glycyl tryptophan arylamidase/gi } @charlist2;#281
	
	map {s/^esterase/**esterase/gi } @charlist2;#283


	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.testspos.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.testspos.txt
		}
	my $temp6 = "temp6.testspos.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.testspos.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "testspos.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @testsposdata = split (/\t/, $line);
		push (@taxnames, $testsposdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Tests Positive Matrix'\;\n\tDIMENSIONS NCHAR=283\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 '2-aminoethanol' \/  'no 2-aminoethanol' '2-aminoethanol', ";
	print OUT "2 '2-keto-d-gluconate' \/  'no 2-keto-d-gluconate' '2-keto-d-gluconate', ";
	print OUT "3 '2-ketobutyrate' \/  'no 2-ketobutyrate' '2-ketobutyrate', ";
	print OUT "4 '2,3-butanediol' \/  'no 2,3-butanediol' '2,3-butanediol', ";
	print OUT "5 '1,3-butanediol' \/  'no 1,3-butanediol' '1,3-butanediol', ";
	print OUT "6 '3-hydroxybenzoate' \/  'no 3-hydroxybenzoate' '3-hydroxybenzoate', ";
	print OUT "7 '3-hydroxybutyrate' \/  'no 3-hydroxybutyrate' '3-hydroxybutyrate', ";
	print OUT "8 '3-methyl d-glucose' \/  'no 3-methyl d-glucose' '3-methyl d-glucose', ";
	print OUT "9 '4-hydroxybenzoate' \/  'no 4-hydroxybenzoate' '4-hydroxybenzoate', ";
	print OUT "10 '4-hydroxyphenylacetate' \/  'no 4-hydroxyphenylacetate' '4-hydroxyphenylacetate', ";

	print OUT "11 'alcohols' \/  'no alcohols' 'alcohols', ";

	print OUT "12 '5-keto-d-gluconate' \/  'no 5-keto-d-gluconate' '5-keto-d-gluconate', ";
	print OUT "13 'a-ketoglutarate' \/  'no a-ketoglutarate' 'a-ketoglutarate', ";
	print OUT "14 'a-ketovalerate' \/  'no a-ketovalerate' 'a-ketovalerate', ";
	print OUT "15 'acetate' \/  'no acetate' 'acetate', ";
	print OUT "16 'acetoacetate' \/  'no acetoacetate' 'acetoacetate', ";
	print OUT "17 'adenosine 5prime-monophosphate' \/  'no adenosine 5prime-monophosphate' 'adenosine 5prime-monophosphate', ";
	print OUT "18 'adenosine' \/  'no adenosine' 'adenosine', ";
	print OUT "19 'adipate' \/  'no adipate' 'adipate', ";

	print OUT "20 'organonitrogen compounds' \/  'no organonitrogen compounds' 'organonitrogen compounds', ";

	print OUT "21 'agar' \/  'no agar' 'agar', ";
	print OUT "22 'alaninamide' \/  'no alaninamide' 'alaninamide', ";
	print OUT "23 'alanine' \/  'no alanine' 'alanine', ";
	print OUT "24 'alanyl glycine' \/  'no alanyl glycine' 'alanyl glycine', ";
	print OUT "25 'amygdalin' \/  'no amygdalin' 'amygdalin', ";
	print OUT "26 'arabinose' \/  'no arabinose' 'arabinose', ";
	print OUT "27 'arabitol' \/  'no arabitol' 'arabitol', ";
	print OUT "28 'arbutin' \/  'no arbutin' 'arbutin', ";
	print OUT "29 'arginine' \/  'no arginine' 'arginine', ";
	print OUT "30 'asparagine' \/  'no asparagine' 'asparagine', ";
	print OUT "31 'aspartate' \/  'no aspartate' 'aspartate', ";
	print OUT "32 'b-cyclodextrin' \/  'no b-cyclodextrin' 'b-cyclodextrin', ";

	print OUT "33 'fatty acids' \/  'no fatty acids' 'fatty acids', ";

	print OUT "34 'b-xylosidase' \/  'no b-xylosidase' 'b-xylosidase', ";
	print OUT "35 'butyrate' \/  'no butyrate' 'butyrate', ";
	print OUT "36 'caprate' \/  'no caprate' 'caprate', ";
	print OUT "37 'carbohydrates' \/  'no carbohydrates' 'carbohydrates', ";
	print OUT "38 'carboxymethyl cellulose' \/  'no carboxymethyl cellulose' 'carboxymethyl cellulose', ";
	print OUT "39 'carnitine' \/  'no carnitine' 'carnitine', ";
	print OUT "40 'cellobiose' \/  'no cellobiose' 'cellobiose', ";
	print OUT "41 'cis-aconitate' \/  'no cis-aconitate' 'cis-aconitate', ";
	print OUT "42 'citrate' \/  'no citrate' 'citrate', ";
	print OUT "43 'co2' \/  'no co2' 'co2', ";
	print OUT "44 'cyclodextrin' \/  'no cyclodextrin' 'cyclodextrin', ";
	print OUT "45 'dextrin' \/  'no dextrin' 'dextrin', ";
	print OUT "46 'dulcitol' \/  'no dulcitol' 'dulcitol', ";
	print OUT "47 'erythritol' \/  'no erythritol' 'erythritol', ";
	print OUT "48 'ethanol' \/  'no ethanol' 'ethanol', ";
	print OUT "49 'formate' \/  'no formate' 'formate', ";
	print OUT "50 'fructose 6-phosphate' \/  'no fructose 6-phosphate' 'fructose 6-phosphate', ";
	print OUT "51 'fructose' \/  'no fructose' 'fructose', ";
	print OUT "52 'fucose' \/  'no fucose' 'fucose', ";
	print OUT "53 'fumarate' \/  'no fumarate' 'fumarate', ";
	print OUT "54 'fusidic acid' \/  'no fusidic acid' 'fusidic acid', ";
	print OUT "55 'galactose' \/  'no galactose' 'galactose', ";
	print OUT "56 'galacturonate' \/  'no galacturonate' 'galacturonate', ";
	print OUT "57 'gelatin' \/  'no gelatin' 'gelatin', ";
	print OUT "58 'gentiobiose' \/  'no gentiobiose' 'gentiobiose', ";
	print OUT "59 'gluconate' \/  'no gluconate' 'gluconate', ";
	print OUT "60 'glucose 1-phosphate' \/  'no glucose 1-phosphate' 'glucose 1-phosphate', ";
	print OUT "61 'glucose' \/  'no glucose' 'glucose', ";
	print OUT "62 'glucuronate' \/  'no glucuronate' 'glucuronate', ";
	print OUT "63 'glutamate' \/  'no glutamate' 'glutamate', ";
	print OUT "64 'glutamic acid polymer' \/  'no glutamic acid polymer' 'glutamic acid polymer', ";
	print OUT "65 'glutamine' \/  'no glutamine' 'glutamine', ";
	print OUT "66 'glutamyl glutamic acid' \/  'no glutamyl glutamic acid' 'glutamyl glutamic acid', ";

	print OUT "67 'amino acids' \/  'no amino acids' 'amino acids', ";

	print OUT "68 'glycerol phosphate' \/  'no glycerol phosphate' 'glycerol phosphate', ";
	print OUT "69 'glycerol' \/  'no glycerol' 'glycerol', ";
	print OUT "70 'glycogen' \/  'no glycogen' 'glycogen', ";
	print OUT "71 'glycyl l-glutamic acid' \/  'no glycyl l-glutamic acid' 'glycyl l-glutamic acid', ";
	print OUT "72 'glycyl-l-proline' \/  'no glycyl-l-proline' 'glycyl-l-proline', ";
	print OUT "73 'glyoxylate' \/  'no glyoxylate' 'glyoxylate', ";
	print OUT "74 'guanidine' \/  'no guanidine' 'guanidine', ";
	print OUT "75 'h2s' \/  'no h2s' 'h2s', ";
	print OUT "76 'hippurate' \/  'no hippurate' 'hippurate', ";
	print OUT "77 'histidine' \/  'no histidine' 'histidine', ";
	print OUT "78 'indole' \/  'no indole' 'indole', ";
	print OUT "79 'inosine' \/  'no inosine' 'inosine', ";
	print OUT "80 'inositol' \/  'no inositol' 'inositol', ";
	print OUT "81 'inulin' \/  'no inulin' 'inulin', ";
	print OUT "82 'itaconate' \/  'no itaconate' 'itaconate', ";
	print OUT "83 'lactate' \/  'no lactate' 'lactate', ";
	print OUT "84 'lactic acid methyl ester' \/  'no lactic acid methyl ester' 'lactic acid methyl ester', ";
	print OUT "85 'lactose' \/  'no lactose' 'lactose', ";
	print OUT "86 'lactulose' \/  'no lactulose' 'lactulose', ";
	print OUT "87 'leucine' \/  'no leucine' 'leucine', ";
	print OUT "88 'lyxose' \/  'no lyxose' 'lyxose', ";
	print OUT "89 'malate' \/  'no malate' 'malate', ";
	print OUT "90 'malonate' \/  'no malonate' 'malonate', ";

	print OUT "91 'carboxylic acids' \/  'no carboxylic acids' 'carboxylic acids', ";

	print OUT "92 'maltose' \/  'no maltose' 'maltose', ";
	print OUT "93 'maltotiose' \/  'no maltotiose' 'maltotiose', ";
	print OUT "94 'mannan' \/  'no mannan' 'mannan', ";
	print OUT "95 'mannitol' \/  'no mannitol' 'mannitol', ";
	print OUT "96 'mannose' \/  'no mannose' 'mannose', ";
	print OUT "97 'melezitose' \/  'no melezitose' 'melezitose', ";
	print OUT "98 'melibiose' \/  'no melibiose' 'melibiose', ";
	print OUT "99 'meso-tartrate' \/  'no meso-tartrate' 'meso-tartrate', ";
	print OUT "100 'methionine' \/  'no methionine' 'methionine', ";
	print OUT "101 'methyl a-d-galactoside' \/  'no methyl a-d-galactoside' 'methyl a-d-galactoside', ";
	print OUT "102 'methyl a-d-glucoside' \/  'no methyl a-d-glucoside' 'methyl a-d-glucoside', ";
	print OUT "103 'methyl b-d-galactoside' \/  'no methyl b-d-galactoside' 'methyl b-d-galactoside', ";
	print OUT "104 'methyl b-d-glucoside' \/  'no methyl b-d-glucoside' 'methyl b-d-glucoside', ";
	print OUT "105 'methyl pyruvate' \/  'no methyl pyruvate' 'methyl pyruvate', ";
	print OUT "106 'monomethyl succinate' \/  'no monomethyl succinate' 'monomethyl succinate', ";
	print OUT "107 'mucate' \/  'no mucate' 'mucate', ";
	print OUT "108 'myo-inositol' \/  'no myo-inositol' 'myo-inositol', ";
	print OUT "109 'n-acetyl-b-d-glucosamine' \/  'no n-acetyl-b-d-glucosamine' 'n-acetyl-b-d-glucosamine', ";
	print OUT "110 'n-acetyl-b-d-mannosamine' \/  'no n-acetyl-b-d-mannosamine' 'n-acetyl-b-d-mannosamine', ";
	print OUT "111 'n-acetyl-l-glutamate' \/  'no n-acetyl-l-glutamate' 'n-acetyl-l-glutamate', ";
	print OUT "112 'ornithine' \/  'no ornithine' 'ornithine', ";
	print OUT "113 'palatinose' \/  'no palatinose' 'palatinose', ";
	print OUT "114 'pectin' \/  'no pectin' 'pectin', ";
	print OUT "115 'phenylacetate' \/  'no phenylacetate' 'phenylacetate', ";
	print OUT "116 'phenylalanine' \/  'no phenylalanine' 'phenylalanine', ";
	print OUT "117 'phenylethylamine' \/  'no phenylethylamine' 'phenylethylamine', ";
	print OUT "118 'proline' \/  'no proline' 'proline', ";
	print OUT "119 'propionate' \/  'no propionate' 'propionate', ";
	print OUT "120 'psicose' \/  'no psicose' 'psicose', ";
	print OUT "121 'pullulan' \/  'no pullulan' 'pullulan', ";
	print OUT "122 'putrescine' \/  'no putrescine' 'putrescine', ";
	print OUT "123 'pyroglutamate' \/  'no pyroglutamate' 'pyroglutamate', ";
	print OUT "124 'pyruvate' \/  'no pyruvate' 'pyruvate', ";

	print OUT "125 'nucleobases' \/  'no nucleobases' 'nucleobases', ";

	print OUT "126 'quinate' \/  'no quinate' 'quinate', ";
	print OUT "127 'raffinose' \/  'no raffinose' 'raffinose', ";
	print OUT "128 'rhamnose' \/  'no rhamnose' 'rhamnose', ";
	print OUT "129 'ribitol' \/  'no ribitol' 'ribitol', ";
	print OUT "130 'ribose' \/  'no ribose' 'ribose', ";
	print OUT "131 'saccharolytic' \/  'no saccharolytic' 'saccharolytic', ";

	print OUT "132 'saccharolytic enzymes' \/  'no saccharolytic enzymes' 'saccharolytic enzymes', ";

	print OUT "133 'salicin' \/  'no salicin' 'salicin', ";
	print OUT "134 'sedoheptulosan' \/  'no sedoheptulosan' 'sedoheptulosan', ";
	print OUT "135 'serine' \/  'no serine' 'serine', ";
	print OUT "136 'sorbitol' \/  'no sorbitol' 'sorbitol', ";
	print OUT "137 'sorbose' \/  'no sorbose' 'sorbose', ";
	print OUT "138 'stachyose' \/  'no stachyose' 'stachyose', ";
	print OUT "139 'starch' \/  'no starch' 'starch', ";
	print OUT "140 'suberate' \/  'no suberate' 'suberate', ";
	print OUT "141 'succinamic acid' \/  'no succinamic acid' 'succinamic acid', ";
	print OUT "142 'succinate' \/  'no succinate' 'succinate', ";

	print OUT "143 'cysteine' \/  'no cysteine' 'cysteine', ";

	print OUT "144 'sucrose' \/  'no sucrose' 'sucrose', ";
	print OUT "145 'tagatose' \/  'no tagatose' 'tagatose', ";
	print OUT "146 'tartrate' \/  'no tartrate' 'tartrate', ";
	print OUT "147 'tellurite' \/  'no tellurite' 'tellurite', ";
	print OUT "148 'threonine' \/  'no threonine' 'threonine', ";
	print OUT "149 'thymidine 5prime-monophosphate' \/  'no thymidine 5prime-monophosphate' 'thymidine 5prime-monophosphate', ";
	print OUT "150 'thymidine' \/  'no thymidine' 'thymidine', ";
	print OUT "151 'trehalose' \/  'no trehalose' 'trehalose', ";

	print OUT "152 'glycine' \/  'no glycine' 'glycine', ";

	print OUT "153 'turanose' \/  'no turanose' 'turanose', ";
	print OUT "154 'tween 40' \/  'no tween 40' 'tween 40', ";
	print OUT "155 'tween 80' \/  'no tween 80' 'tween 80', ";
	print OUT "156 'uridine' \/  'no uridine' 'uridine', ";
	print OUT "157 'valerate' \/  'no valerate' 'valerate', ";
	print OUT "158 'valine' \/  'no valine' 'valine', ";
	print OUT "159 'xylan' \/  'no xylan' 'xylan', ";
	print OUT "160 'xylitol' \/  'no xylitol' 'xylitol', ";
	print OUT "161 'xylose' \/  'no xylose' 'xylose', ";

	print OUT "162 'isoleucine' \/  'no isoleucine' 'isoleucine', ";

	print OUT "163 'lysine' \/  'no lysine' 'lysine', ";
	print OUT "164 'tryptophan' \/  'no tryptophan' 'tryptophan', ";
	print OUT "165 'selenocysteine' \/  'no selenocysteine' 'selenocysteine', ";
	print OUT "166 'pyrrolysine' \/  'no pyrrolysine' 'pyrrolysine', ";

	print OUT "167 'nitrate reduction' \/  'no nitrate reduction' 'nitrate reduction', ";
	print OUT "168 'nitrogen fixation' \/  'no nitrogen fixation' 'nitrogen fixation', ";
	print OUT "169 'nitrite reduction' \/  'no nitrite reduction' 'nitrite reduction', ";

	print OUT "170 '6-phospho-b-galactosidase' \/  'no 6-phospho-b-galactosidase' '6-phospho-b-galactosidase', ";
	print OUT "171 '6-phosphogluconate dehydrogenase' \/  'no 6-phosphogluconate dehydrogenase' '6-phosphogluconate dehydrogenase', ";
	print OUT "172 'a-arabinosidase' \/  'no a-arabinosidase' 'a-arabinosidase', ";
	print OUT "173 'a-fucosidase' \/  'no a-fucosidase' 'a-fucosidase', ";
	print OUT "174 'a-galactosidase' \/  'no a-galactosidase' 'a-galactosidase', ";
	print OUT "175 'a-glucosidase' \/  'no a-glucosidase' 'a-glucosidase', ";
	print OUT "176 'a-mannosidase' \/  'no a-mannosidase' 'a-mannosidase', ";
	print OUT "177 'a-methyl glucosidase' \/  'no a-methyl glucosidase' 'a-methyl glucosidase', ";
	print OUT "178 'a-methylmannoside' \/  'no a-methylmannoside' 'a-methylmannoside', ";
	print OUT "179 'acid phosphatase' \/  'no acid phosphatase' 'acid phosphatase', ";
	print OUT "180 'aesculin hydrolysis' \/  'no aesculin hydrolysis' 'aesculin hydrolysis', ";
	print OUT "181 'alanine arylamidase' \/  'no alanine arylamidase' 'alanine arylamidase', ";
	print OUT "182 'alkaline phosphatase' \/  'no alkaline phosphatase' 'alkaline phosphatase', ";
	print OUT "183 'alpha-glutamate arylamidase' \/  'no alpha-glutamate arylamidase' 'alpha-glutamate arylamidase', ";
	print OUT "184 'amylase' \/  'no amylase' 'amylase', ";
	print OUT "185 'arginine arylamidase' \/  'no arginine arylamidase' 'arginine arylamidase', ";
	print OUT "186 'arginine decarboxylase' \/  'no arginine decarboxylase' 'arginine decarboxylase', ";
	print OUT "187 'arginine dihydrolase' \/  'no arginine dihydrolase' 'arginine dihydrolase', ";
	print OUT "188 'arylsulfatase' \/  'no arylsulfatase' 'arylsulfatase', ";
	print OUT "189 'asparagine arylamidase' \/  'no asparagine arylamidase' 'asparagine arylamidase', ";
	print OUT "190 'aspartate arylamidase' \/  'no aspartate arylamidase' 'aspartate arylamidase', ";
	print OUT "191 'b-galactosidase' \/  'no b-galactosidase' 'b-galactosidase', ";
	print OUT "192 'b-glucosidase' \/  'no b-glucosidase' 'b-glucosidase', ";
	print OUT "193 'b-glucuronidase' \/  'no b-glucuronidase' 'b-glucuronidase', ";
	print OUT "194 'b-mannosidase' \/  'no b-mannosidase' 'b-mannosidase', ";
	print OUT "195 'casein hydrolysis' \/  'no casein hydrolysis' 'casein hydrolysis', ";
	print OUT "196 'caseinase' \/  'no caseinase' 'caseinase', ";
	print OUT "197 'catalase' \/  'no catalase' 'catalase', ";
	print OUT "198 'cellobiose' \/  'no cellobiose' 'cellobiose', ";
	print OUT "199 'chymotrypsin' \/  'no chymotrypsin' 'chymotrypsin', ";
	print OUT "200 'coagulase' \/  'no coagulase' 'coagulase', ";
	print OUT "201 'collagenase' \/  'no collagenase' 'collagenase', ";
	print OUT "202 'cysteine arylamidase' \/  'no cysteine arylamidase' 'cysteine arylamidase', ";
	print OUT "203 'cystine arylamidase' \/  'no cystine arylamidase' 'cystine arylamidase', ";
	print OUT "204 'cytochrome oxidase' \/  'no cytochrome oxidase' 'cytochrome oxidase', ";
	print OUT "205 'dnase' \/  'no dnase' 'dnase', ";
	print OUT "206 'esterase c4' \/  'no esterase c4' 'esterase c4', ";
	print OUT "207 'fructose bisphosphate aldolase' \/  'no fructose bisphosphate aldolase' 'fructose bisphosphate aldolase', ";
	print OUT "208 'fucosidase' \/  'no fucosidase' 'fucosidase', ";
	print OUT "209 'gelatinase' \/  'no gelatinase' 'gelatinase', ";
	print OUT "210 'glucosidase' \/  'no glucosidase' 'glucosidase', ";
	print OUT "211 'glutamate arylamidase' \/  'no glutamate arylamidase' 'glutamate arylamidase', ";
	print OUT "212 'glutamate deaminase' \/  'no glutamate deaminase' 'glutamate deaminase', ";
	print OUT "213 'glutamic acid decarboxylase' \/  'no glutamic acid decarboxylase' 'glutamic acid decarboxylase', ";
	print OUT "214 'glutamine arylamidase' \/  'no glutamine arylamidase' 'glutamine arylamidase', ";
	print OUT "215 'glutamyl glutamic acid arylamidase' \/  'no glutamyl glutamic acid arylamidase' 'glutamyl glutamic acid arylamidase', ";
	print OUT "216 'glycine arylamidase' \/  'no glycine arylamidase' 'glycine arylamidase', ";
	print OUT "217 'hippurate hydrolysis' \/  'no hippurate hydrolysis' 'hippurate hydrolysis', ";
	print OUT "218 'histidine arylamidase' \/  'no histidine arylamidase' 'histidine arylamidase', ";
	print OUT "219 'hyaluronidase' \/  'no hyaluronidase' 'hyaluronidase', ";
	print OUT "220 'hydroxy proline arylamidase' \/  'no hydroxy proline arylamidase' 'hydroxy proline arylamidase', ";
	print OUT "221 'iron-phorphyrins' \/  'no iron-phorphyrins' 'iron-phorphyrins', ";
	print OUT "222 'isoleucine arylamidase' \/  'no isoleucine arylamidase' 'isoleucine arylamidase', ";
	print OUT "223 'l-lactic acid racemase' \/  'no l-lactic acid racemase' 'l-lactic acid racemase', ";
	print OUT "224 'lecithinase' \/  'no lecithinase' 'lecithinase', ";
	print OUT "225 'leucine arylamidase' \/  'no leucine arylamidase' 'leucine arylamidase', ";
	print OUT "226 'leucyl glycine arylamidase' \/  'no leucyl glycine arylamidase' 'leucyl glycine arylamidase', ";
	print OUT "227 'lipase c14' \/  'no lipase c14' 'lipase c14', ";
	print OUT "228 'lipase c8' \/  'no lipase c8' 'lipase c8', ";
	print OUT "229 'lipase' \/  'no lipase' 'lipase', ";
	print OUT "230 'lysine arylamidase' \/  'no lysine arylamidase' 'lysine arylamidase', ";
	print OUT "231 'lysine deaminase' \/  'no lysine deaminase' 'lysine deaminase', ";
	print OUT "232 'lysine decarboxylase' \/  'no lysine decarboxylase' 'lysine decarboxylase', ";
	print OUT "233 'methionine arylamidase' \/  'no methionine arylamidase' 'methionine arylamidase', ";
	print OUT "234 'methyl a-d-mannoside' \/  'no methyl a-d-mannoside' 'methyl a-d-mannoside', ";
	print OUT "235 'methyl b-d-xyloside' \/  'no methyl b-d-xyloside' 'methyl b-d-xyloside', ";
	print OUT "236 'n-acetyl-b-d-glucosaminidase' \/  'no n-acetyl-b-d-glucosaminidase' 'n-acetyl-b-d-glucosaminidase', ";
	print OUT "237 'nadh oxidase' \/  'no nadh oxidase' 'nadh oxidase', ";
	print OUT "238 'naphthol-as-bi-phosphohydrolase' \/  'no naphthol-as-bi-phosphohydrolase' 'naphthol-as-bi-phosphohydrolase', ";
	print OUT "239 'ornithine arylamidase' \/  'no ornithine arylamidase' 'ornithine arylamidase', ";
	print OUT "240 'ornithine decarboxylase' \/  'no ornithine decarboxylase' 'ornithine decarboxylase', ";
	print OUT "241 'oxidase' \/  'no oxidase' 'oxidase', ";
	print OUT "242 'para-phenylalanine deaminase' \/  'no para-phenylalanine deaminase' 'para-phenylalanine deaminase', ";
	print OUT "243 'phenylalanine arylamidase' \/  'no phenylalanine arylamidase' 'phenylalanine arylamidase', ";
	print OUT "244 'phenylalanine deaminase' \/  'no phenylalanine deaminase' 'phenylalanine deaminase', ";
	print OUT "245 'phosphatase' \/  'no phosphatase' 'phosphatase', ";
	print OUT "246 'proline arylamidase' \/  'no proline arylamidase' 'proline arylamidase', ";
	print OUT "247 'protease' \/  'no protease' 'protease', ";
	print OUT "248 'pyrazinamidase' \/  'no pyrazinamidase' 'pyrazinamidase', ";
	print OUT "249 'pyroglutamic acid arylamidase' \/  'no pyroglutamic acid arylamidase' 'pyroglutamic acid arylamidase', ";
	print OUT "250 'proteolytic enzymes' \/  'no proteolytic enzymes' 'proteolytic enzymes', ";
	print OUT "251 'serine arylamidase' \/  'no serine arylamidase' 'serine arylamidase', ";
	print OUT "252 'threonine arylamidase' \/  'no threonine arylamidase' 'threonine arylamidase', ";
	print OUT "253 'trypsin' \/  'no trypsin' 'trypsin', ";
	print OUT "254 'tryptophan arylamidase' \/  'no tryptophan arylamidase' 'tryptophan arylamidase', ";
	print OUT "255 'tryptophan deaminase' \/  'no tryptophan deaminase' 'tryptophan deaminase', ";
	print OUT "256 'tryptophanase' \/  'no tryptophanase' 'tryptophanase', ";
	print OUT "257 'tyrosine arylamidase' \/  'no tyrosine arylamidase' 'tyrosine arylamidase', ";
	print OUT "258 'urease' \/  'no urease' 'urease', ";
	print OUT "259 'valine arylamidase' \/  'no valine arylamidase' 'valine arylamidase', ";
	print OUT "260 'xylanase' \/  'no xylanase' 'xylanase', ";

	print OUT "261 'a-hemolysis' \/  'no a-hemolysis' 'a-hemolysis', ";
	print OUT "262 'acetoin' \/  'no acetoin' 'acetoin', ";
	print OUT "263 'b-haemolysis' \/  'no b-haemolysis' 'b-haemolysis', ";
	print OUT "264 'ellman' \/  'no ellman' 'ellman', ";
	print OUT "265 'flexirubin' \/  'no flexirubin' 'flexirubin', ";
	print OUT "266 'haemolysis' \/  'no haemolysis' 'haemolysis', ";
	print OUT "267 'indoxyl acetate hydrolysis' \/  'no indoxyl acetate hydrolysis' 'indoxyl acetate hydrolysis', ";
	print OUT "268 'methyl red acidification' \/  'no methyl red acidification' 'methyl red acidification', ";
	print OUT "269 'methylene blue reduction' \/  'no methylene blue reduction' 'methylene blue reduction', ";
	print OUT "270 'resazurin reduction' \/  'no resazurin reduction' 'resazurin reduction', ";
	print OUT "271 'tetrathionate reduction' \/  'no tetrathionate reduction' 'tetrathionate reduction', ";
	print OUT "272 'tetrazolium reduction' \/  'no tetrazolium reduction' 'tetrazolium reduction', ";
	print OUT "273 'carboxymethyl cellulase' \/  'no carboxymethyl cellulase' 'carboxymethyl cellulase', ";

	print OUT "274 '1,3-propanediol' \/  'no 1,3-propanediol' '1,3-propanediol', ";
	print OUT "275 '1,2-propanediol' \/  'no 1,2-propanediol' '1,2-propanediol', ";
	print OUT "276 'propanediol' \/  'no propanediol' 'propanediol', ";
	print OUT "277 'phosphatidylinositol phospholipase c' \/  'no phosphatidylinositol phospholipase c' 'phosphatidylinositol phospholipase c', ";
	print OUT "278 'lecithinase' \/  'no lecithinase' 'lecithinase', ";
	print OUT "279 'lactamide' \/  'no lactamide' 'lactamide', ";
	print OUT "280 'hydroxyproline' \/  'no hydroxyproline' 'hydroxyproline', ";
	print OUT "281 'glycyl tryptophan arylamidase' \/  'no glycyl tryptophan arylamidase' 'glycyl tryptophan arylamidase', ";
	print OUT "282 'lipolytic enzymes' \/  'no lipolytic enzymes' 'lipolytic enzymes', ";
	print OUT "283 'esterase' \/  'no esterase' 'esterase', ";
#	print OUT "284 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "285 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "286 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "287 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "288 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "289 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "290 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "291 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "292 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "293 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "294 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "295 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "296 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "297 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "298 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "299 'xxx' \/  'no xxx' 'xxx', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @testsposdata = split (/\t/, $line);
		push (my @taxnames, $testsposdata[0]);

#code char 1 2-aminoethanol
		if ($line =~ /,2-aminoethanol,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not 2-aminoethanol,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
			}
#code char 2 2-keto-d-gluconate
		if ($line =~ /,2-keto-d-gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-keto-d-gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 3 2-ketobutyrate
		if ($line =~ /,2-ketobutyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-ketobutyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 4 2,3-butanediol
		if ($line =~ /,2-3-butanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-3-butanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 5 1,3-butanediol
		if ($line =~ /,1-3-butanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-3-butanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 6 3-hydroxybenzoate
		if ($line =~ /,3-hydroxybenzoate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-hydroxybenzoate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 7 3-hydroxybutyrate
		if ($line =~ /,3-hydroxybutyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-hydroxybutyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 8 3-methyl d-glucose
		if ($line =~ /,3-methyl d-glucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-methyl d-glucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 9 4-hydroxybenzoate
		if ($line =~ /,4-hydroxybenzoate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 4-hydroxybenzoate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 10 4-hydroxyphenylacetate
		if ($line =~ /,4-hydroxyphenylacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 4-hydroxyphenylacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 11 alcohols
		if ($line =~ /,1-3-propanediol,|,1-2-propanediol,|,propanediol,|,lactamide,|,alcohols,|,2-aminoethanol,|,2-3-butanediol,|,1-3-butanediol,|,arabitol,|,dulcitol,|,ethythritol,|,ethanol,|,glycerol phosphate,|,glycerol,|,inositol,|,mannitol,|,myo-inositol,|,ribitol,|,salicin,|,sorbitol,|,xylitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alcohols,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 12 5-keto-d-gluconate
		if ($line =~ /,5-keto-d-gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 5-keto-d-gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 13 a-ketoglutarate
		if ($line =~ /,a-ketoglutarate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-ketoglutarate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 14 a-ketovalerate
		if ($line =~ /,a-ketovalerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-ketovalerate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 15 acetate
		if ($line =~ /,acetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 16 acetoacetate
		if ($line =~ /,acetoacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetoacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 17 adenosine 5'-monophosphate
		if ($line =~ /,adenosine 5\'-monophosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adenosine 5\'-monophosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 18 adenosine
		if ($line =~ /,adenosine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adenosine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 19 adipate
		if ($line =~ /,adipate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adipate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 20 organonitrogen compounds
		if ($line =~ /,organonitrogen compounds,|,2-aminoethanol,|,alaninamide,|,alanine,|,alanyl glycine,|,amygdalin,|,asparagine,|,aspartate,|,carnitine,|,glutamic acid polymer,|,glutamine,|,glutamyl glutamic acid,|,glycyl l-glutamic acid,|,glycyl-l-proline,|,guanidine,|,hippurate,|,histidine,|,leucine,|,methionine,|,n-acetyl-b-d-glucosamine,|,n-acetyl-b-d-mannosamine,|,n-acetyl-l-glutamate,|,ornithine,|,phenylalanine,|,phenylethylamine,|,proline,|,putrescine,|,pyroglutamate,|,serine,|,threonine,|,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not organonitrogen compounds,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 21 agar
		if ($line =~ /,agar,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not agar,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 22 alaninamide
		if ($line =~ /,alaninamide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alaninamide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 23 alanine
		if ($line =~ /,alanine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 24 alanyl glycine
		if ($line =~ /,alanyl glycine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanyl glycine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 25 amygdalin
		if ($line =~ /,amygdalin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amygdalin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 26 arabinose
		if ($line =~ /,arabinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arabinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 27 arabitol
		if ($line =~ /,arabitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arabitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 28 arbutin
		if ($line =~ /,arbutin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arbutin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 29 arginine
		if ($line =~ /,arginine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 30 asparagine
		if ($line =~ /,asparagine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not asparagine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 31 aspartate
		if ($line =~ /,aspartate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aspartate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 32 b-cyclodextrin
		if ($line =~ /,b-cyclodextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-cyclodextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 33 fatty acids
		if ($line =~ /,fatty acids,|,2-ketyobutyrate,|,3-hydroxybutyrate,|,a-ketovalerate,|,butyrate,|,propionate,|,valerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fatty acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 34 b-xylosidase
		if ($line =~ /,b-xylosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-xylosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 35 butyrate
		if ($line =~ /,butyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not butyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 36 caprate
		if ($line =~ /,caprate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not caprate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 37 carbohydrates
		if ($line =~ /,carbohydrates,|,3-methyl d-glucose,|,agar,|,amygdalin,|,arabinose,|,arabitol,|,arbutin,|,b-cyclodextrin,|,carboxymethyl cellulose,|,cellobiose,|,cyclodextrin,|,dextrin,|,dulcitol,|,erythritol,|,fructose 6-phosphate,|,fructose,|,fucose,|,galactose,|,galacturonate,|,gelatin,|,gentiobiose,|,gluconate,|,glucose 1-phosphate,|,glucose,|,glucuronate,|,glycerol phosphate,|,glycerol,|,glycogen,|,inulin,|,lactose,|,lactulose,|,lyxose,|,maltose,|,maltotriose,|,mannan,|,mannitol,|,mannose,|,melezitose,|,melibiose,|,meso-tartrate,|,methyl a-d-galactoside,|,methyl a-d-glucoside,|,methyl b-d-glactoside,|,methyl b-d-glucoside,|,mucate,|,n-acetyl-b-d-glucosamine,|,n-acetyl-b-d-mannosamine,|,palatinose,|,pectin,|,psicose,|,pullulan,|,raffinose,|,rhamnose,|,ribitol,|,ribose,|,salicin,|,sedoheptulosan,|,sorbitol,|,sorbose,|,stachyose,|,starch,|,sucrose,|,tagatose,|,tartrate,|,trehalose,|,turanose,|,xylan,|,xylitol,|,xylose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carbohydrates,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 38 carboxymethyl cellulose
		if ($line =~ /,carboxymethyl cellulose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxymethyl cellulose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 39 carnitine
		if ($line =~ /,carnitine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carnitine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 40 cellobiose
		if ($line =~ /,cellobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cellobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 41 cis-aconitate
		if ($line =~ /,cis-aconitate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cis-aconitate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 42 citrate
		if ($line =~ /,citrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not citrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 43 co2
		if ($line =~ /,co2,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not co2,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 44 cyclodextrin
		if ($line =~ /,cyclodextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cyclodextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 45 dextrin
		if ($line =~ /,dextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 46 dulcitol
		if ($line =~ /,dulcitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dulcitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 47 erythritol
		if ($line =~ /,erythritol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erythritol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 48 ethanol
		if ($line =~ /,ethanol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ethanol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 49 formate
		if ($line =~ /,formate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not formate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 50 fructose 6-phosphate
		if ($line =~ /,fructose 6-phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose 6-phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 51 fructose
		if ($line =~ /,fructose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 52 fucose
		if ($line =~ /,fucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 53 fumarate
		if ($line =~ /,fumarate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fumarate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 54 fusidic acid
		if ($line =~ /,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fusidic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 55 galactose
		if ($line =~ /,galactose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not galactose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 56 galacturonate
		if ($line =~ /,galacturonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not galacturonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 57 gelatin
		if ($line =~ /,gelatin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gelatin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 58 gentiobiose
		if ($line =~ /,gentiobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gentiobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 59 gluconate
		if ($line =~ /,gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 60 glucose 1-phosphate
		if ($line =~ /,glucose 1-phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucose 1-phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 61 glucose
		if ($line =~ /,glucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 62 glucuronate
		if ($line =~ /,glucuronate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucuronate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 63 glutamate
		if ($line =~ /,glutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 64 glutamic acid polymer
		if ($line =~ /,glutamic acid polymer,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamic acid polymer,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 65 glutamine
		if ($line =~ /,glutamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 66 glutamyl glutamic acid
		if ($line =~ /,glutamyl glutamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamyl glutamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 67 amino acids
		if ($line =~ /,hydroxyproline,|,amino acids,|,alaninamide,|,alanine,|,alanyl glycine,|,arginine,|,asparagine,|,aspartate,|,glutamate,|,glutamic acid polymer,|,glutamine,|,glutamyl glutamic acid,|,glycine,|,glycyl l-glutamic acid,|,glycyl-l-proline,|,histidine,|,isoleucine,|,leucine,|,lysine,|,methionine,|,ornithine,|,phenylalanine,|,proline,|,pyroglutamate,|,pyrrolysine,|,selenocysteine,|,serine,|,threonine,|,tryptophan,|,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amino acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 68 glycerol phosphate
		if ($line =~ /,glycerol phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycerol phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 69 glycerol
		if ($line =~ /,glycerol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycerol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 70 glycogen
		if ($line =~ /,glycogen,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycogen,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 71 glycyl l-glutamic acid
		if ($line =~ /,glycyl l-glutamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl l-glutamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 72 glycyl-l-proline
		if ($line =~ /,glycyl-l-proline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl-l-proline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 73 glyoxylate
		if ($line =~ /,glyoxylate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glyoxylate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 74 guanidine
		if ($line =~ /,guanidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not guanidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 75 h2s
		if ($line =~ /,h2s,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not h2s,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 76 hippurate
		if ($line =~ /,hippurate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hippurate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 77 histidine
		if ($line =~ /,histidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not histidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 78 indole
		if ($line =~ /,indole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not indole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 79 inosine
		if ($line =~ /,inosine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inosine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 80 inositol
		if ($line =~ /,inositol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inositol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 81 inulin
		if ($line =~ /,inulin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inulin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 82 itaconate
		if ($line =~ /,itaconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not itaconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 83 lactate
		if ($line =~ /,lactate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 84 lactic acid methyl ester
		if ($line =~ /,lactic acid methyl ester,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactic acid methyl ester,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 85 lactose
		if ($line =~ /,lactose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 86 lactulose
		if ($line =~ /,lactulose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactulose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 87 leucine
		if ($line =~ /,leucine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 88 lyxose
		if ($line =~ /,lyxose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lyxose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 89 malate
		if ($line =~ /,malate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not malate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 90 malonate
		if ($line =~ /,malonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not malonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 91 carboxylic acids
		if ($line =~ /,hydroxyproline,|,carboxylic acids,|,2-keto-d-gluconate,|,2-ketobutyrate,|,3-hydroxybenzoate,|,3-hydroxybutyrate,|,4-hydroxybenzoate,|,4-hydroxyphenylacetate,|,5-keto-d-gluconate,|,a-ketoglutarate,|,a-ketovalerate,|,acetate,|,acetoacetate,|,adipate,|,butyrate,|,cis-aconitate,|,citrate,|,formate,|,fumarate,|,fusidic acid,|,galacturonate,|,gluconate,|,glucuronate,|,glutamate,|,glyoxylate,|,itaconate,|,lactate,|,malate,|,malonate,|,meso-tartrate,|,methyl pyruvate,|,monomethyl succinate,|,mucate,|,phenylacetate,|,pyruvate,|,quinate,|,suberate,|,succinamic acid,|,tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxylic acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 92 maltose
		if ($line =~ /,maltose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not maltose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 93 maltotiose
		if ($line =~ /,maltotiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not maltotiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 94 mannan
		if ($line =~ /,mannan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 95 mannitol
		if ($line =~ /,mannitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 96 mannose
		if ($line =~ /,mannose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 97 melezitose
		if ($line =~ /,melezitose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not melezitose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 98 melibiose
		if ($line =~ /,melibiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not melibiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 99 meso-tartrate
		if ($line =~ /,meso-tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meso-tartrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 100 methionine
		if ($line =~ /,methionine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methionine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 101 methyl a-d-galactoside
		if ($line =~ /,methyl a-d-galactoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-galactoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 102 methyl a-d-glucoside
		if ($line =~ /,methyl a-d-glucoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-glucoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 103 methyl b-d-galactoside
		if ($line =~ /,methyl b-d-galactoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-galactoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 104 methyl b-d-glucoside
		if ($line =~ /,methyl b-d-glucoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-glucoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 105 methyl pyruvate
		if ($line =~ /,methyl pyruvate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl pyruvate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 106 monomethyl succinate
		if ($line =~ /,monomethyl succinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not monomethyl succinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 107 mucate
		if ($line =~ /,mucate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mucate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 108 myo-inositol
		if ($line =~ /,myo-inositol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not myo-inositol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 109 n-acetyl-b-d-glucosamine
		if ($line =~ /,n-acetyl-b-d-glucosamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-glucosamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 110 n-acetyl-b-d-mannosamine
		if ($line =~ /,n-acetyl-b-d-mannosamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-mannosamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 111 n-acetyl-l-glutamate
		if ($line =~ /,n-acetyl-l-glutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-l-glutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 112 ornithine
		if ($line =~ /,ornithine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 113 palatinose
		if ($line =~ /,palatinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not palatinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 114 pectin
		if ($line =~ /,pectin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pectin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 115 phenylacetate
		if ($line =~ /,phenylacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 116 phenylalanine
		if ($line =~ /,phenylalanine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 117 phenylethylamine
		if ($line =~ /,phenylethylamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylethylamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 118 proline
		if ($line =~ /,proline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not proline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 119 propionate
		if ($line =~ /,propionate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propionate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 120 psicose
		if ($line =~ /,psicose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not psicose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 121 pullulan
		if ($line =~ /,pullulan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pullulan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 122 putrescine
		if ($line =~ /,putrescine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not putrescine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 123 pyroglutamate
		if ($line =~ /,pyroglutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyroglutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 124 pyruvate
		if ($line =~ /,pyruvate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyruvate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 125 nucleobases
		if ($line =~ /,nucleobases,|,adenosine 5\'-monophosphate,|,adenosine,|,thymidine 5\'-monophosphate,|,thymidine,|,uridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nucleobases,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 126 quinate
		if ($line =~ /,quinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not quinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 127 raffinose
		if ($line =~ /,raffinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not raffinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 128 rhamnose
		if ($line =~ /,rhamnose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rhamnose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 129 ribitol
		if ($line =~ /,ribitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ribitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 130 ribose
		if ($line =~ /,ribose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ribose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 131 saccharolytic
		if ($line =~ /,saccharolytic,|,agar,|,amygdaline,|,b-cyclodextrin,|,carboxymethyl cellulose,|,cellobiose,|,cyclodextrin,|,dextrin,|,gelatin,|,gentiobiose,|,glycogen,|,inulin,|,lactose,|,lactulose,|,maltose,|,maltotriose,|,mannan,|,melezitose,|,melibiose,|,palatinose,|,pectin,|,pullulan,|,raffinose,|,stachyose,|,starch,|,trehalose,|,turanose,|,xylan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not saccharolytic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 132 saccharolytic enzymes
		if ($line =~ /,6-phospho-b-galactosidase,|,a-arabinosidase,|,a-fucosidase,|,a-galactosidase,|,a-glucosidase,|,a-mannosidase,|,a-methyl glucosidase,|,amylase,|,b-galactosidase,|,b-glucosidase,|,b-glucuronidase,|,b-mannosidase,|,fucosidase,|,glucosidase,|,hyaluronidase,|,n-acetyl-b-d-glucosaminidase,|,xylanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 133 salicin
		if ($line =~ /,salicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not salicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 134 sedoheptulosan
		if ($line =~ /,sedoheptulosan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sedoheptulosan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 135 serine
		if ($line =~ /,serine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not serine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 136 sorbitol
		if ($line =~ /,sorbitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sorbitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 137 sorbose
		if ($line =~ /,sorbose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sorbose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 138 stachyose
		if ($line =~ /,stachyose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not stachyose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 139 starch
		if ($line =~ /,starch,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not starch,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 140 suberate
		if ($line =~ /,suberate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not suberate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 141 succinamic acid
		if ($line =~ /,succinamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not succinamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 142 succinate
		if ($line =~ /,succinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not succinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 143 cysteine
		if ($line =~ /,cysteine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cysteine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 144 sucrose
		if ($line =~ /,sucrose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sucrose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 145 tagatose
		if ($line =~ /,tagatose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tagatose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 146 tartrate
		if ($line =~ /,tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tartrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 147 tellurite
		if ($line =~ /,tellurite,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tellurite,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 148 threonine
		if ($line =~ /,threonine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not threonine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 149 thymidine 5'-monophosphate
		if ($line =~ /,thymidine 5\'-monophosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not thymidine 5\'-monophosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 150 thymidine
		if ($line =~ /,thymidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not thymidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 151 trehalose
		if ($line =~ /,trehalose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trehalose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 152 glycine
		if ($line =~ /,glycine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 153 turanose
		if ($line =~ /,turanose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not turanose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 154 tween 40
		if ($line =~ /,tween 40,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tween 40,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 155 tween 80
		if ($line =~ /,tween 80,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tween 80,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 156 uridine
		if ($line =~ /,uridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not uridine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 157 valerate
		if ($line =~ /,valerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valerate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 158 valine
		if ($line =~ /,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 159 xylan
		if ($line =~ /,xylan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 160 xylitol
		if ($line =~ /,xylitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 161 xylose
		if ($line =~ /,xylose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 162 isoleucine
		if ($line =~ /,isoleucine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isoleucine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 163 lysine
		if ($line =~ /,lysine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 164 tryptophan
		if ($line =~ /,tryptophan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 165 selenocysteine
		if ($line =~ /,selenocysteine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not selenocysteine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 166 pyrrolysine
		if ($line =~ /,pyrrolysine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrrolysine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 167 nitrate reduction
		if ($line =~ /,nitrate reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrate reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 168 nitrogen fixation
		if ($line =~ /,nitrogen fixation,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrogen fixation,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 169 nitrite reduction
		if ($line =~ /,nitrite reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrite reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 170 6-phospho-b-galactosidase
		if ($line =~ /,6-phospho-b-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 6-phospho-b-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 171 6-phosphogluconate dehydrogenase
		if ($line =~ /,6-phosphogluconate dehydrogenase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 6-phosphogluconate dehydrogenase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 172 a-arabinosidase
		if ($line =~ /,a-arabinosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-arabinosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 173 a-fucosidase
		if ($line =~ /,a-fucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-fucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 174 a-galactosidase
		if ($line =~ /,a-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 175 a-glucosidase
		if ($line =~ /,a-glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 176 a-mannosidase
		if ($line =~ /,a-mannosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-mannosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 177 a-methyl glucosidase
		if ($line =~ /,a-methyl glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-methyl glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 178 a-methylmannoside
		if ($line =~ /,a-methylmannoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-methylmannoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 179 acid phosphatase
		if ($line =~ /,acid phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acid phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 180 aesculin hydrolysis
		if ($line =~ /,aesculin hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aesculin hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 181 alanine arylamidase
		if ($line =~ /,alanine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 182 alkaline phosphatase
		if ($line =~ /,alkaline phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alkaline phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 183 alpha-glutamate arylamidase
		if ($line =~ /,alpha-glutamate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alpha-glutamate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 184 amylase
		if ($line =~ /,amylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 185 arginine arylamidase
		if ($line =~ /,arginine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 186 arginine decarboxylase
		if ($line =~ /,arginine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 187 arginine dihydrolase
		if ($line =~ /,arginine dihydrolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine dihydrolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 188 arylsulfatase
		if ($line =~ /,arylsulfatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arylsulfatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 189 asparagine arylamidase
		if ($line =~ /,asparagine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not asparagine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 190 aspartate arylamidase
		if ($line =~ /,aspartate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aspartate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 191 b-galactosidase
		if ($line =~ /,b-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 192 b-glucosidase
		if ($line =~ /,b-glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 193 b-glucuronidase
		if ($line =~ /,b-glucuronidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-glucuronidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 194 b-mannosidase
		if ($line =~ /,b-mannosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-mannosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 195 casein hydrolysis
		if ($line =~ /,casein hydrolysis,|,caseinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not casein hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 196 caseinase
		if ($line =~ /,caseinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not caseinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 197 catalase
		if ($line =~ /,catalase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not catalase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 198 cellobiose
		if ($line =~ /,cellobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cellobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 199 chymotrypsin
		if ($line =~ /,chymotrypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chymotrypsin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 200 coagulase
		if ($line =~ /,coagulase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not coagulase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 201 collagenase
		if ($line =~ /,collagenase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not collagenase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 202 cysteine arylamidase
		if ($line =~ /,cysteine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cysteine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 203 cystine arylamidase
		if ($line =~ /,cystine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cystine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 204 cytochrome oxidase
		if ($line =~ /,cytochrome oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cytochrome oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 205 dnase
		if ($line =~ /,dnase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dnase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 206 esterase c4
		if ($line =~ /,esterase c4,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not esterase c4,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 207 fructose bisphosphate aldolase
		if ($line =~ /,fructose bisphosphate aldolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose bisphosphate aldolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 208 fucosidase
		if ($line =~ /,fucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 209 gelatinase
		if ($line =~ /,gelatinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gelatinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 210 glucosidase
		if ($line =~ /,glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 211 glutamate arylamidase
		if ($line =~ /,glutamate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 212 glutamate deaminase
		if ($line =~ /,glutamate deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 213 glutamic acid decarboxylase
		if ($line =~ /,glutamic acid decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamic acid decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 214 glutamine arylamidase
		if ($line =~ /,glutamine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 215 glutamyl glutamic acid arylamidase
		if ($line =~ /,glutamyl glutamic acid arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamyl glutamic acid arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 216 glycine arylamidase
		if ($line =~ /,glycine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 217 hippurate hydrolysis
		if ($line =~ /,hippurate hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hippurate hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 218 histidine arylamidase
		if ($line =~ /,histidine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not histidine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 219 hyaluronidase
		if ($line =~ /,hyaluronidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hyaluronidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 220 hydroxy proline arylamidase
		if ($line =~ /,hydroxy proline arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hydroxy proline arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 221 iron-phorphyrins
		if ($line =~ /,iron-phorphyrins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not iron-phorphyrins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 222 isoleucine arylamidase
		if ($line =~ /,isoleucine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isoleucine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 223 l-lactic acid racemase
		if ($line =~ /,l-lactic acid racemase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not l-lactic acid racemase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 224 lecithinase
		if ($line =~ /,lecithinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lecithinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 225 leucine arylamidase
		if ($line =~ /,leucine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 226 leucyl glycine arylamidase
		if ($line =~ /,leucyl glycine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucyl glycine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 227 lipase c14
		if ($line =~ /,lipase c14,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase c14,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 228 lipase c8
		if ($line =~ /,lipase c8,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase c8,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 229 lipase
		if ($line =~ /,lipase,|,lipase c8,|,lipase c14,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 230 lysine arylamidase
		if ($line =~ /,lysine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 231 lysine deaminase
		if ($line =~ /,lysine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 232 lysine decarboxylase
		if ($line =~ /,lysine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 233 methionine arylamidase
		if ($line =~ /,methionine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methionine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 234 methyl a-d-mannoside
		if ($line =~ /,methyl a-d-mannoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-mannoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 235 methyl b-d-xyloside
		if ($line =~ /,methyl b-d-xyloside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-xyloside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 236 n-acetyl-b-d-glucosaminidase
		if ($line =~ /,n-acetyl-b-d-glucosaminidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-glucosaminidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 237 nadh oxidase
		if ($line =~ /,nadh oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nadh oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 238 naphthol-as-bi-phosphohydrolase
		if ($line =~ /,naphthol-as-bi-phosphohydrolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not naphthol-as-bi-phosphohydrolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 239 ornithine arylamidase
		if ($line =~ /,ornithine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 240 ornithine decarboxylase
		if ($line =~ /,ornithine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 241 oxidase
		if ($line =~ /,oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 242 para-phenylalanine deaminase
		if ($line =~ /,para-phenylalanine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not para-phenylalanine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 243 phenylalanine arylamidase
		if ($line =~ /,phenylalanine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 244 phenylalanine deaminase
		if ($line =~ /,phenylalanine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 245 phosphatase
		if ($line =~ /,phosphatase,,|,alkaline phosphatase,|,acid phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 246 proline arylamidase
		if ($line =~ /,proline arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not proline arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 247 protease
		if ($line =~ /,protease,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not protease,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 248 pyrazinamidase
		if ($line =~ /,pyrazinamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrazinamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 249 pyroglutamic acid arylamidase
		if ($line =~ /,pyroglutamic acid arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyroglutamic acid arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 250 proteolytic enzymes
		if ($line =~ /,caseinase,|,casein hydrolysis,|,chymotrypsin,|,coagulase,|,collagenase,|,gelatinase,|,protease,|,trypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 251 serine arylamidase
		if ($line =~ /,serine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not serine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 252 threonine arylamidase
		if ($line =~ /,threonine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not threonine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 253 trypsin
		if ($line =~ /,trypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trypsin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 254 tryptophan arylamidase
		if ($line =~ /,tryptophan arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 255 tryptophan deaminase
		if ($line =~ /,tryptophan deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 256 tryptophanase
		if ($line =~ /,tryptophanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophanase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 257 tyrosine arylamidase
		if ($line =~ /,tyrosine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tyrosine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 258 urease
		if ($line =~ /,urease,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not urease,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 259 valine arylamidase
		if ($line =~ /,valine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 260 xylanase
		if ($line =~ /,xylanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylanase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 261 a-hemolysis
		if ($line =~ /,a-hemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-hemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 262 acetoin
		if ($line =~ /,acetoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetoin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 263 b-haemolysis
		if ($line =~ /,b-haemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-haemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 264 ellman
		if ($line =~ /,ellman,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ellman,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 265 flexirubin
		if ($line =~ /,flexirubin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flexirubin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 266 haemolysis
		if ($line =~ /,haemolysis,|,a-haemolysis,|,b-haemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not haemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 267 indoxyl acetate hydrolysis
		if ($line =~ /,indoxyl acetate hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not indoxyl acetate hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 268 methyl red acidification
		if ($line =~ /,methyl red acidification,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl red acidification,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 269 methylene blue reduction
		if ($line =~ /,methylene blue reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methylene blue reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 270 resazurin reduction
		if ($line =~ /,resazurin reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not resazurin reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 271 tetrathionate reduction
		if ($line =~ /,tetrathionate reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetrathionate reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 272 tetrazolium reduction
		if ($line =~ /,tetrazolium reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetrazolium reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 273 carboxymethyl cellulase
		if ($line =~ /,carboxymethyl cellulase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxymethyl cellulase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 274 1,3-propanediol
		if ($line =~ /,1-3-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-3-propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 275 1,2-propanediol
		if ($line =~ /,1-2-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-2-propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 276 propanediol
		if ($line =~ /,propanediol,|,1-2-propanediol,|,1-3-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 277 phosphatidylinositol phospholipase c
		if ($line =~ /,phosphatidylinositol phospholipase c,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phosphatidylinositol phospholipase c,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 278 lecithinase
		if ($line =~ /,lecithinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lecithinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 279 lactamide
		if ($line =~ /,lactamide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactamide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 280 hydroxyproline
		if ($line =~ /,hydroxyproline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hydroxyproline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 281 glycyl tryptophan arylamidase
		if ($line =~ /,glycyl tryptophan arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl tryptophan arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 282 lipolytic enzymes
		if ($line =~ /,esterase c4,|,lecithinase,|,lipase c14,|,lipase c8,|,lipase,|,esterase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 283 esterase
		if ($line =~ /,esterase c4,|,esterase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}

		print OUT "\n";
		}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;

	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}

###################
#process tests negative
###################
#first discover the character states by homologizing them to the MicrO ontology
elsif ($character eq "tests negative") {
	my $homout = "tests.negative.txt";
	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>', $homout) or die $!;
	local $, = "\t";
	while (my $line = <IN> ) { # pushes the elements into an array and homologizes the terms in the array
		chomp $line;
		$line =~ s/\xce\xb2/b/g; #changes greek betas β to b
		$line =~ s/\xc3\x9f/b/g; #changes greek eszett β to b
		$line =~ s/\xce\xb1/a/g;#changes greek alphas α to a
		$line =~ s/\xce\x91/a/g;#changes greek capital alphas to a
		$line =~ s/\xce\x92/b/g;#changes greek capital betas to b
		$line =~ s/\xc2\xb1/ /g;#changes plus or minus to space
		my @columns = split /\t/, $line;
		push (my @taxlabels, $columns[0]);
		push (my @testsneg, $columns[1]);
		#ontology terms and synonyms
		$_ = lc for @testsneg;
		map {s/,\d\.\d\%\s,/,,/gi; } @testsneg; # synonyms of 
		map {s/,\d\.\d\d\%\s,/,,/gi; } @testsneg; # synonyms of 


	#substrates used and/or fermented
		map {s/,.*% fructose,/,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,.*% glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,.*% l-lactate,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,.*% l-lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,.*% potassium tellurite,/,tellurite,/gi; } @testsneg; # synonyms of 
		map {s/,.*% sucrose,/,sucrose,/gi; } @testsneg; # synonyms of 
		map {s/,.*% tween 80,/,tween 80,/gi; } @testsneg; # synonyms of 
		map {s/,.*mm succinate,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d\d% formate,/,formate,/gi; } @testsneg; # synonyms of 
		map {s/,1 2 propanediol,/,1-2-propanediol,/gi; } @testsneg; # synonyms of 
		map {s/,1 2-propanediol,/,1-2-propanediol,/gi; } @testsneg; # synonyms of 
		map {s/,1 3-butanediol,/,1-3-butanediol,/gi; } @testsneg; # synonyms of 
		map {s/,1 3-propanediol,/,1-3-propanediol,/gi; } @testsneg; # synonyms of 
		map {s/,1 g glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,2 3 butanediol,/,2-3-butanediol,/gi; } @testsneg; # synonyms of 
		map {s/,2 3-butanediol,/,2-3-butanediol,/gi; } @testsneg; # synonyms of 
		map {s/,2 moles lactate,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,2-keto-gluconate,/,2-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,2-ketoglutarate,/,2-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,20 mm pyruvate,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,25 mm leucine,/,leucine,/gi; } @testsneg; # synonyms of 
		map {s/,25 mm phenylalanine,/,phenylalanine,/gi; } @testsneg; # synonyms of 
		map {s/,25 mm threonine,/,threonine,/gi; } @testsneg; # synonyms of 
		map {s/,25 mm valine,/,valine,/gi; } @testsneg; # synonyms of 
		map {s/,3-hydroxybenzoic acid,/,3-hydroxybenzoate,/gi; } @testsneg; # synonyms of 
		map {s/,3-hydroxybutyrate,/,3-hydroxybutyrate,/gi; } @testsneg; # synonyms of 
		map {s/,4-hydroxybenzoic acid,/,4-hydroxybenzoate,/gi; } @testsneg; # synonyms of 
		map {s/,5-keto-d gluconate,/,5-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,a-d-glucose 1-phosphate,/,glucose 1-phosphate,/gi; } @testsneg; # synonyms of 
		map {s/,a-d-glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,a-ketoglutaric acid,/,a-ketoglutarate,/gi; } @testsneg; # synonyms of 
		map {s/,a-ketovaleric acid,/,a-ketovalerate,/gi; } @testsneg; # synonyms of 
		map {s/,acetate andethanol,/,acetate,ethanol,/gi; } @testsneg; # synonyms of 
		map {s/,acetate strip test,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,acetate synthesis,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,acetate utilization,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,acetate,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,acetic acid,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,acetic and lactic acid,/,acetate,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,acetoacetic acid,/,acetoacetate,/gi; } @testsneg; # synonyms of 
		map {s/,additional glutamic acid,/,glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,adenosine 5,/,adenosine 5\'-monophosphate,/gi; } @testsneg; # synonyms of 
		map {s/,adipate,/,adipate,/gi; } @testsneg; # synonyms of 
		map {s/,adipic acid,/,adipate,/gi; } @testsneg; # synonyms of 
		map {s/,alanine,/,alanine,/gi; } @testsneg; # synonyms of 
		map {s/,amygdaline,/,amygdalin,/gi; } @testsneg; # synonyms of 
		map {s/,anaerobically d-ribose,/,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,anaerobically glucose,/,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,approximately 1 mol l-alanine,/,alanine,/gi; } @testsneg; # synonyms of 
		map {s/,arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,aspartic acid,/,aspartate,/gi; } @testsneg; # synonyms of 
		map {s/,b-methyl-d-glucoside,/,methyl b-d-glucoside,/gi; } @testsneg; # synonyms of 
		map {s/,b-phenylethylamine,/,phenylethylamine,/gi; } @testsneg; # synonyms of 
		map {s/,beijerinckii raffinose,/,raffinose,/gi; } @testsneg; # synonyms of 
		map {s/,biomass using malonate,/,malonate,/gi; } @testsneg; # synonyms of 
		map {s/,both d-and l-lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,butyrate plus co2,/,butyrate,co2,/gi; } @testsneg; # synonyms of 
		map {s/,butyric acid,/,butyrate,/gi; } @testsneg; # synonyms of 
		map {s/,caprate,/,caprate,/gi; } @testsneg; # synonyms of 
		map {s/,capric acid assimilation,/,caprate,/gi; } @testsneg; # synonyms of 
		map {s/,capric acid,/,caprate,/gi; } @testsneg; # synonyms of 
		map {s/,carbohydrate derivatives including glucose,/,carbohydrates,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,carbohydrates including arabinose,/,carbohydrates,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,carbohydrates ribose,/,carbohydrates,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,cell suspensions ferment glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose fermentation,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose help,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,chemoorganoheterotroph utilizing cellobiose,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,citrate fermentation,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,citrate hydrolysis,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,citrate utilization,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,citrate,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,citric acid,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,corresponding fatty acid plus acetic acid,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,d cellobiose,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,d fructose,/,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,d glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,d lactose,/,lactose,/gi; } @testsneg; # synonyms of 
		map {s/,d mannitol,/,mannitol,/gi; } @testsneg; # synonyms of 
		map {s/,d mannose,/,mannose,/gi; } @testsneg; # synonyms of 
		map {s/,d trehalose,/,trehalose,/gi; } @testsneg; # synonyms of 
		map {s/,d xylose,/,xylose,/gi; } @testsneg; # synonyms of 
		map {s/,d-adonitol,/,ribitol,/gi; } @testsneg; # synonyms of 
		map {s/,d-alanine,/,alanine,/gi; } @testsneg; # synonyms of 
		map {s/,d-amygdalin,/,amygdalin,/gi; } @testsneg; # synonyms of 
		map {s/,d-arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,d-arabitol,/,arabitol,/gi; } @testsneg; # synonyms of 
		map {s/,d-aspartic acid,/,aspartate,/gi; } @testsneg; # synonyms of 
		map {s/,d-cellobiose,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,d-fructose 6-phosphate,/,fructose 6-phosphate,/gi; } @testsneg; # synonyms of 
		map {s/,d-fructose-6-po4,/,fructose 6-phosphate,/gi; } @testsneg; # synonyms of 
		map {s/,d-fructose,/,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,d-fucose,/,fucose,/gi; } @testsneg; # synonyms of 
		map {s/,d-galactose,/,galactose,/gi; } @testsneg; # synonyms of 
		map {s/,d-galacturonate,/,galacturonate,/gi; } @testsneg; # synonyms of 
		map {s/,d-galacturonic acid,/,galacturonate,/gi; } @testsneg; # synonyms of 
		map {s/,d-gluconate,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,d-gluconic acid,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,d-glucose 1-phosphate,/,glucose 1-phosphate,/gi; } @testsneg; # synonyms of 
		map {s/,d-glucose fermentation,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,d-glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,d-glucuronate,/,glucuronate,/gi; } @testsneg; # synonyms of 
		map {s/,d-glucuronic acid,/,glucuronate,/gi; } @testsneg; # synonyms of 
		map {s/,d-lactic acid methyl ester,/,lactic acid methyl ester,/gi; } @testsneg; # synonyms of 
		map {s/,d-lactose,/,lactose,/gi; } @testsneg; # synonyms of 
		map {s/,d-lyxose,/,lyxose,/gi; } @testsneg; # synonyms of 
		map {s/,d-malate,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,d-malic acid,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,d-maltose,/,maltose,/gi; } @testsneg; # synonyms of 
		map {s/,d-mannitol,/,mannitol,/gi; } @testsneg; # synonyms of 
		map {s/,d-mannose,/,mannose,/gi; } @testsneg; # synonyms of 
		map {s/,d-melezitose,/,melezitose,/gi; } @testsneg; # synonyms of 
		map {s/,d-melibiose,/,melibiose,/gi; } @testsneg; # synonyms of 
		map {s/,d-psicose,/,psicose,/gi; } @testsneg; # synonyms of 
		map {s/,d-raffinose,/,raffinose,/gi; } @testsneg; # synonyms of 
		map {s/,d-ribitol,/,ribitol,/gi; } @testsneg; # synonyms of 
		map {s/,d-ribose,/,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,d-saccharose,/,sucrose,/gi; } @testsneg; # synonyms of 
		map {s/,d-salicin,/,salicin,/gi; } @testsneg; # synonyms of 
		map {s/,d-sorbitol,/,sorbitol,/gi; } @testsneg; # synonyms of 
		map {s/,d-tagatose,/,tagatose,/gi; } @testsneg; # synonyms of 
		map {s/,d-trehalose,/,trehalose,/gi; } @testsneg; # synonyms of 
		map {s/,d-turanose,/,turanose,/gi; } @testsneg; # synonyms of 
		map {s/,d-xylose,/,xylose,/gi; } @testsneg; # synonyms of 
		map {s/,dextrin crystals,/,dextrin,/gi; } @testsneg; # synonyms of 
		map {s/,dextrin utilization,/,dextrin,/gi; } @testsneg; # synonyms of 
		map {s/,disodium succinate,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,dl-carnitine,/,carnitine,/gi; } @testsneg; # synonyms of 
		map {s/,dl-lactate,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,either l-arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,esculin nor starch,/,aesculin,starch,/gi; } @testsneg; # synonyms of 
		map {s/,exogenous citrate,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,fewer strains ferment cellobiose,/,cellobiose,/gi; } @testsneg; # synonyms of 
		map {s/,formate being major end products,/,formate,/gi; } @testsneg; # synonyms of 
		map {s/,formate yield succinate,/,formate,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,formic acid,/,formate,/gi; } @testsneg; # synonyms of 
		map {s/,formic and lactic acid,/,formate,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,galactitol,/,dulcitol,/gi; } @testsneg; # synonyms of 
		map {s/,gluconic acid sodium salt,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,gluconic acid,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,glucose acidification,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,glucose fermentation products,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,glucose fermentation,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,glucose metabolism,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,glucose supplementation,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,glucuronic acid,/,glucuronate,/gi; } @testsneg; # synonyms of 
		map {s/,glutamic acid,/,glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,glycerol 1 3-propanediol,/,glycerol,1-3-propanediol,/gi; } @testsneg; # synonyms of 
		map {s/,glycogen utilization,/,glycogen,/gi; } @testsneg; # synonyms of 
		map {s/,glycogene,/,glycogen,/gi; } @testsneg; # synonyms of 
		map {s/,h2s production,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,high xylan content,/,xylan,/gi; } @testsneg; # synonyms of 
		map {s/,hydrogen sulfide,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,hydrogen sulphide production,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses gelatin,/,gelatin,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses starch,/,starch,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses tween 40,/,tween 40,/gi; } @testsneg; # synonyms of 
		map {s/,including acetate,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,including arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,indole production,/,indole,/gi; } @testsneg; # synonyms of 
		map {s/,itaconic acid,/,itaconate,/gi; } @testsneg; # synonyms of 
		map {s/,keto butyric acid,/,2-ketobutyrate,/gi; } @testsneg; # synonyms of 
		map {s/,l-alaninamide,/,alaninamide,/gi; } @testsneg; # synonyms of 
		map {s/,l-alanine,/,alanine,/gi; } @testsneg; # synonyms of 
		map {s/,L-alanine,/,alanine,/gi; } @testsneg; # synonyms of 
		map {s/,l-alanyl-glycine,/,alanyl glycine,/gi; } @testsneg; # synonyms of 
		map {s/,l-arabinose d-maltose,/,arabinose,maltose,/gi; } @testsneg; # synonyms of 
		map {s/,l-arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,l-arabitol,/,arabitol,/gi; } @testsneg; # synonyms of 
		map {s/,l-arginine,/,arginine,/gi; } @testsneg; # synonyms of 
		map {s/,l-asparagine,/,asparagine,/gi; } @testsneg; # synonyms of 
		map {s/,l-aspartate,/,aspartate,/gi; } @testsneg; # synonyms of 
		map {s/,l-aspartic acid,/,aspartate,/gi; } @testsneg; # synonyms of 
		map {s/,l-fucose,/,fucose,/gi; } @testsneg; # synonyms of 
		map {s/,l-glutamate,/,glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,l-glutamic acid,/,glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,l-glutamine,/,glutamine,/gi; } @testsneg; # synonyms of 
		map {s/,l-histidine,/,histidine,/gi; } @testsneg; # synonyms of 
		map {s/,l-hydroxyproline,/,hydroxyproline,/gi; } @testsneg; # synonyms of 
		map {s/,l-lactate isomers,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,l-lactate,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,l-lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,l-leucine,/,leucine,/gi; } @testsneg; # synonyms of 
		map {s/,l-malate,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,l-malic acid,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,l-methionine,/,methionine,/gi; } @testsneg; # synonyms of 
		map {s/,l-ornithine,/,ornithine,/gi; } @testsneg; # synonyms of 
		map {s/,l-phenylalanine,/,phenylalanine,/gi; } @testsneg; # synonyms of 
		map {s/,l-proline,/,proline,/gi; } @testsneg; # synonyms of 
		map {s/,l-pyroglutamic acid,/,pyroglutamate,/gi; } @testsneg; # synonyms of 
		map {s/,l-rhamnose fermentation,/,rhamnose,/gi; } @testsneg; # synonyms of 
		map {s/,l-rhamnose,/,rhamnose,/gi; } @testsneg; # synonyms of 
		map {s/,l-serine,/,serine,/gi; } @testsneg; # synonyms of 
		map {s/,l-sorbose,/,sorbose,/gi; } @testsneg; # synonyms of 
		map {s/,l-threonine,/,threonine,/gi; } @testsneg; # synonyms of 
		map {s/,l-valine,/,valine,/gi; } @testsneg; # synonyms of 
		map {s/,l-xylose,/,xylose,/gi; } @testsneg; # synonyms of 
		map {s/,lactate either,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,lactic acid group,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,lactic acid isomers,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,lactic acid production,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,lactose help,/,lactose,/gi; } @testsneg; # synonyms of 
		map {s/,little glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,m-inositol,/,myo-inositol,/gi; } @testsneg; # synonyms of 
		map {s/,mainly acetate,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,mainly butyrate,/,butyrate,/gi; } @testsneg; # synonyms of 
		map {s/,mainly lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,major lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,malic acid assimilation,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,malic acid,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,malonate utilization,/,malonate,/gi; } @testsneg; # synonyms of 
		map {s/,maltose assimilation,/,maltose,/gi; } @testsneg; # synonyms of 
		map {s/,mannose fermentation,/,mannose,/gi; } @testsneg; # synonyms of 
		map {s/,mannose support weak growth,/,mannose,/gi; } @testsneg; # synonyms of 
		map {s/,methyl a-d-glucopyranoside,/,methyl a-d-glucoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl a-d-glucose,/,3-methyl d-glucose,/gi; } @testsneg; # synonyms of 
		map {s/,methyl a-d-mannopyranoside,/,methyl a-d-mannoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl b-d-galactopyranoside,/,methyl b-d-galactoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl b-d-xylopyranoside,/,methyl b-d-xyloside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl pyruvate,/,methyl pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,mol glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,mole glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,most strains ferment glucose,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,most strains ferment mannitol,/,mannitol,/gi; } @testsneg; # synonyms of 
		map {s/,mucic acid,/,mucate,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-b-glucosaminidase activities,/,n-acetyl-b-d-glucosamine,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-d-glucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-glucosamine assimilation,/,xxx,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-glucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-l-glutamic acid,/,n-acetyl-l-glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetylglucosamine,/,n-acetyl-b-d-glucosamine,/gi; } @testsneg; # synonyms of 
		map {s/,negative and gelatin hydrolysis,/,gelatin,/gi; } @testsneg; # synonyms of 
		map {s/,negative trehalose fermentation,/,trehalose,/gi; } @testsneg; # synonyms of 
		map {s/,nor fumarate,/,fumarate,/gi; } @testsneg; # synonyms of 
		map {s/,oat spelt xylan,/,xylan,/gi; } @testsneg; # synonyms of 
		map {s/,only a-ketovaleric acid,/,a-ketovalerate,/gi; } @testsneg; # synonyms of 
		map {s/,only acetic acid,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,only d-ribose,/,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,only gelatin,/,gelatin,/gi; } @testsneg; # synonyms of 
		map {s/,only glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,only l-lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,only lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,organic acids pyruvate,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,ornithine using moeller,/,ornithine,/gi; } @testsneg; # synonyms of 
		map {s/,p-hydroxyphenylacetic acid,/,4-hydroxyphenylacetate,/gi; } @testsneg; # synonyms of 
		map {s/,phenyl acetic acid,/,phenylacetate,/gi; } @testsneg; # synonyms of 
		map {s/,phenylacetic acid,/,phenylacetate,/gi; } @testsneg; # synonyms of 
		map {s/,poorly glycerol,/,glycerol,/gi; } @testsneg; # synonyms of 
		map {s/,potassium 2-ketogluconate,/,2-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,potassium 5-ketogluconate,/,5-keto-d-gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,potassium gluconate assimi-lation,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,potassium gluconate assimilation,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,potassium gluconate,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,potassium tellurite,/,tellurite,/gi; } @testsneg; # synonyms of 
		map {s/,prefers pyruvate,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,propionic acid,/,propionate,/gi; } @testsneg; # synonyms of 
		map {s/,putrescine fermentation,/,putrescine,/gi; } @testsneg; # synonyms of 
		map {s/,pyroglutamic acid,/,pyroglutamate,/gi; } @testsneg; # synonyms of 
		map {s/,pyruvate fermentation,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,pyruvate utilization,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,pyruvic acid methyl ester,/,methyl pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,pyruvic acid,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,quinic acid,/,quinate,/gi; } @testsneg; # synonyms of 
		map {s/,raffinose cultures,/,raffinose,/gi; } @testsneg; # synonyms of 
		map {s/,raffinose d-fructose,/,raffinose,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,reduces fumarate,/,fumarate,/gi; } @testsneg; # synonyms of 
		map {s/,rhamnose utilization,/,rhamnose,/gi; } @testsneg; # synonyms of 
		map {s/,s226 fermented glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,saccharose,/,sucrose,/gi; } @testsneg; # synonyms of 
		map {s/,salicine d-sorbitol,/,salicin,sorbitol,/gi; } @testsneg; # synonyms of 
		map {s/,sodium acetate,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium citrate,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium gluconate,/,gluconate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium glutamate,/,glutamate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium hippurate rind esculin,/,hippurate,aesculin,/gi; } @testsneg; # synonyms of 
		map {s/,sodium hippurate,/,hippurate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium l-lactate,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium malate,/,malate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium pyruvate,/,pyruvate,/gi; } @testsneg; # synonyms of 
		map {s/,sodium succinate,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,solely l-lactic acid,/,lactate,/gi; } @testsneg; # synonyms of 
		map {s/,soluble potato starch,/,starch,/gi; } @testsneg; # synonyms of 
		map {s/,soluble starch,/,starch,/gi; } @testsneg; # synonyms of 
		map {s/,sometimes formate,/,formate,/gi; } @testsneg; # synonyms of 
		map {s/,sometimes succinate,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,starch hydrolysis,/,starch,/gi; } @testsneg; # synonyms of 
		map {s/,starch nor esculin,/,starch,aesculin,/gi; } @testsneg; # synonyms of 
		map {s/,starch utilization,/,starch,/gi; } @testsneg; # synonyms of 
		map {s/,strains ferment arabinose,/,arabinose,/gi; } @testsneg; # synonyms of 
		map {s/,strains ferment fructose,/,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,strains ferment lactose,/,lactose,/gi; } @testsneg; # synonyms of 
		map {s/,strains weakly ferment glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,suberic acid,/,suberate,/gi; } @testsneg; # synonyms of 
		map {s/,succinic acid monomethyl ester,/,monomethyl succinate,/gi; } @testsneg; # synonyms of 
		map {s/,succinic acid,/,succinate,/gi; } @testsneg; # synonyms of 
		map {s/,sucrose fermentation,/,sucrose,/gi; } @testsneg; # synonyms of 
		map {s/,thymidine 5,/,thymidine 5\'-monophosphate,/gi; } @testsneg; # synonyms of 
		map {s/,trisodium citrate assimilation,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,trisodium citrate tests,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,trisodium citrate,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,tween 80 hydrolysis variable,/,tween 80,/gi; } @testsneg; # synonyms of 
		map {s/,tween 80 hydrolysis,/,tween 80,/gi; } @testsneg; # synonyms of 
		map {s/,tween 80 type,/,tween 80,/gi; } @testsneg; # synonyms of 
		map {s/,uridine 5,/,uridine 5\'-monophosphate,/gi; } @testsneg; # synonyms of 
		map {s/,using glycerol,/,glycerol,/gi; } @testsneg; # synonyms of 
		map {s/,usually acetate,/,acetate,/gi; } @testsneg; # synonyms of 
		map {s/,usually formate,/,formate,/gi; } @testsneg; # synonyms of 
		map {s/,usually glucose,/,glucose,/gi; } @testsneg; # synonyms of 
		map {s/,utilizes citrate,/,citrate,/gi; } @testsneg; # synonyms of 
		map {s/,valeric acid,/,valerate,/gi; } @testsneg; # synonyms of 
		map {s/,valine,/,valine,/gi; } @testsneg; # synonyms of 
		map {s/,weakly saccharolytic and ferment fructose,/,saccharolytic,fructose,/gi; } @testsneg; # synonyms of 
		map {s/,xylan fermentation,/,xylan,/gi; } @testsneg; # synonyms of 
		map {s/,xylan ribose,/,xylan,ribose,/gi; } @testsneg; # synonyms of 
		map {s/,xylan support growth,/,xylan,/gi; } @testsneg; # synonyms of 
		map {s/,xylose fermentation,/,xylose,/gi; } @testsneg; # synonyms of 

	#inorganics
		map {s/,aerobic nitrate reduction,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolysed and nitrate reduction,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate and indole production,/,nitrate reduction,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate and nitrite reduction,/,nitrate reduction,nitrite reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reductase a.,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reductase activities,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reductase,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction activity,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction reaction,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction reactions,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction test,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction tests,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction variable,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction-positive,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate respiration,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate supplementation,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrite and nitrate reduction,/,nitrate reduction,nitrite reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrite reductase activity,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nitrite reductase,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,nor nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,reduced nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,reducing nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,sodium nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,using nitrate,/,nitrate reduction,/gi; } @testsneg; # synonyms of 

		map {s/,h2s formation,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,h2s gas,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,hydrogen sulfide production,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,hydrogen sulphide,/,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,moderate h2s,/,h2s,/gi; } @testsneg; # synonyms of 

	#enzymes present
		map {s/,aesculin ferric citrate reaction,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin ferric citrate,/,aesculin hydrolysi,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin ferric,/,aesculin hydrolysi,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin hydrolysi,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin hydrolysis,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,bile . aesculin test,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,esculin ferric citrate,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,esculin hydrolysis,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,esculin,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,ferric citrate,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses aesculin,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses esculin,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,positive aesculin hydrolysis,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,positive aesculin reaction,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,type strain hydrolyses aesculin,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin dihydrolase,/,aesculin hydrolysis,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 

		map {s/,cytochrome oxidase activities,/,cytochrome oxidase,/gi; } @testsneg; # synonyms of 

		map {s/,2-naphthyl phosphate,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,2-naphthyl phosphate,/,phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,acid phosphatase activities,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,acid phosphatase no3,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,acid phosphatase,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,acidic phosphatase,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,bis-p-nitrophenyl-phosphate,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,generally weak alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,inconstant alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,moderate alkaline phosphatase activity,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,moderate alkaline phosphatase,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,only acid phosphatase,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,p-nitrophenyl phosphate,/,phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatase acid activity,/,phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatase activity,/,phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatase reaction,/,phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatase-positive and ferment lactose,/,phosphatase,lactose,/gi; } @testsneg; # synonyms of 
		map {s/,present and acid phosphatase,/,acid phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,weak alkaline phosphatase activities,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 
		map {s/,weak alkaline phosphatase,/,alkaline phosphatase,/gi; } @testsneg; # synonyms of 

		map {s/,2-naphthyl myristate,/,lipase c14,/gi; } @testsneg; # synonyms of 
		map {s/,c14 lipase,/,lipase c14,/gi; } @testsneg; # synonyms of 
		map {s/,lipase c 14,/,lipase c14,/gi; } @testsneg; # synonyms of 
		map {s/,lipase c-14,/,lipase c14,/gi; } @testsneg; # synonyms of 

		map {s/,2-naphthyl butyrate,/,esterase c4,/gi; } @testsneg; # synonyms of 
		map {s/,c-4 esterase,/,esterase c4,/gi; } @testsneg; # synonyms of 
		map {s/,c4 esterase,/,esterase c4,/gi; } @testsneg; # synonyms of 
		map {s/,esterase c-4,/,esterase c4,/gi; } @testsneg; # synonyms of 
		map {s/,esterase polymorphism,/,esterase,/gi; } @testsneg; # synonyms of 
		map {s/,esterase production,/,esterase,/gi; } @testsneg; # synonyms of 
		map {s/,esterase,/,esterase c4,/gi; } @testsneg; # synonyms of 
		map {s/,negative including esterase,/,esterase,/gi; } @testsneg; # synonyms of 

		map {s/,2-naphthyl caprylate,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,c-8 esterase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,c8 esterase lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,C8 esterase lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,c8 esterase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,c8 lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,ester lipase c-8,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,ester lipase c8,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,ester lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,esterase lipase c8,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,esterase lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 
		map {s/,present and esterase lipase,/,lipase c8,/gi; } @testsneg; # synonyms of 

		map {s/,egg yolk lecithinase,/,lecithinase,/gi; } @testsneg; # synonyms of 
		map {s/,egg-yolk lecithinase,/,lecithinase,/gi; } @testsneg; # synonyms of 
		map {s/,lecithinase activity,/,lecithinase,/gi; } @testsneg; # synonyms of 

		map {s/,gelatin decomposition,/,gelatinase,/gi; } @testsneg; # synonyms of 
		map {s/,gelatin hydrolysis,/,gelatinase,/gi; } @testsneg; # synonyms of 
		map {s/,gelatin liquefaction,/,gelatinase,/gi; } @testsneg; # synonyms of 
		map {s/,gelatinase activities,/,gelatinase,/gi; } @testsneg; # synonyms of 
		map {s/,gelatinase activity,/,gelatinase,/gi; } @testsneg; # synonyms of 
		map {s/,gelatinase,/,gelatinase,/gi; } @testsneg; # synonyms of 

		map {s/,n-acetyl-beta-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-b-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testsneg; # synonyms of 

	#amino acid decarboxylases
		map {s/,arginine decarboxylase activity,/,arginine decarboxylase,/gi; } @testsneg; # synonyms of 
		map {s/,l-ornithine decarboxylases,/,ornithine decarboxylase,/gi; } @testsneg; # synonyms of 
		map {s/,lysine & ornithine decarboxylases,/,ornithine decarboxylase,lysine decarboxylase,/gi; } @testsneg; # synonyms of 
		map {s/,lysine decarboxylase tests,/,lysine decarboxylase,/gi; } @testsneg; # synonyms of 
		map {s/,ornithine decarboxylase activities,/,ornithine decarboxylase,/gi; } @testsneg; # synonyms of 
		map {s/,ornithine decarboxylases,/,ornithine decarboxylase,/gi; } @testsneg; # synonyms of 

	#amino acid deaminases
		map {s/,alanine dihydrolase activity,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,alanine dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dehydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase 2,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase activities,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase activity,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase hydrolysis,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase oxidase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase production,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase tests,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine-dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,glutamate dehydrogenase,/,glutamate deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,indole activities,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole and nitrate reduction activity,/,tryptophanase,nitrate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,indole hydrolysis,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole negative,/,not tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole negative,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole nor h2s,/,tryptophanase,h2s,/gi; } @testsneg; # synonyms of 
		map {s/,indole positive,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole reactions,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole test,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,indole tests,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,inducible tryptophanase,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,l-arginine dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,l-arginine hydrolysis,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,l-glutamate dehydrogenase,/,glutamate deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine deaminase activities,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine deaminase activity,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine deaminase tests,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine deamination,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine dehydrogenase,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,phenylalanine hydrolysis,/,phenylalanine deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,positive indole reaction,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,produces indole,/,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,tryptophan arylamidase,/,tryptophan deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,tryptophan deaminase activities,/,tryptophan deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,tryptophan deaminase activity,/,tryptophan deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,tryptophane deaminase,/,tryptophan deaminase,/gi; } @testsneg; # synonyms of 
		map {s/,tyrosine arylamidase activities,/,tyrosine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 

	#amino acid arylamidases
		map {s/,alanine amd,/,alanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine arylamidase activity,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine arylamidase tests,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine arylamidase,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,arginine dihydrolase hydrolysis,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,cystine aminopeptidase,/,cystine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,cystine arylamidase activities,/,cystine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,cystine arylamidase activity,/,cystine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,enzymes arginine arylamidase,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,glycine 7-amido-4-methylcoumarin,/,glycine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,glycine amd,/,glycine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,histidine 7-amido-4-methylcoumarin,/,histidine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,histidine amd,/,histidine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,histidine arylamidase activities,/,histidine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,histidine arylamidase activity,/,histidine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,inconstant alanine arylamidase activities,/,alanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-alanine amino-peptidase,/,alanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-arginine 7-amido-4-methylcoumarin,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-arginine amc,/,arginine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-aspartate arylamidase,/,aspartate arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-isoleucine arylamidase,/,isoleucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-phenylalanine 7-amido-4-methylcoumarin,/,phenylalanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-phenylalanine 7-amido-4,/,phenylalanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-phenylalanine amc,/,phenylalanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-proline arylamidase,/,proline arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-pyroglutamic acid 7-amido-4-methyl coumarin,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-pyroglutamic acid 7-amido-4-methylcoumarin,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-pyroglutamic acid amc,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-valine 7-amido-4-methylcoumarin,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,l-valine amc,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine aminopeptidase activity,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine aminopeptidase,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylamidase activities,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylamidase activity,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylamidase reaction,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylamidase tests,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylamidases,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylaminidase,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine arylaminidases,/,leucine arylaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine b-naphthylamide,/,leucine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine b-naphthylamide,/,leucine arylaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine p-nitroanilide,/,leucine arylaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,leucine pnitroanilide,/,leucine arylaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,phenyl alanine arylamidase,/,phenylalanine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,positive pyrrolidonyl arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,proline aminopeptidase,/,proline arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyrazinamidase activities,/,pyrazinamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyrazinamidase activity,/,pyrazinamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyroglutamate arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyroglutamic acid arylamidase activities,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyroglutamic acid arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyrrolidonyl arylamidase activity,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,serine arylamidase activities,/,serine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,using commercial api systems arginine dihydrolase,/,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,valine aminopeptidase,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylamidase activities,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylamidase pyrrolidonyl arylamidase,/,valine arylamidase,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylamidases,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylamidasetrypsin activities,/,valine arylamidase,trypsin,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylaminidase,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,valine arylarnidase,/,valine arylamidase,/gi; } @testsneg; # synonyms of 
		map {s/,weak pyrrolidonyl arylamidase,/,pyroglutamic acid arylamidase,/gi; } @testsneg; # synonyms of 

		map {s/,coagulase activity,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,fixed coagulase,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,free coagulase,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,rapid coagulase test,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,test tube coagulase test reaction,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,test tube coagulase test,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,tube coagulase test,/,coagulase,/gi; } @testsneg; # synonyms of 
		map {s/,tube coagulase,/,coagulase,/gi; } @testsneg; # synonyms of 

		map {s/,benzidine,/,iron-phorphyrins,/gi; } @testsneg; # synonyms of 
		map {s/,benzidine test,/,iron-phorphyrins,/gi; } @testsneg; # synonyms of 
		map {s/,benzidine tests,/,iron-phorphyrins,/gi; } @testsneg; # synonyms of 
		map {s/,benzidine test positive,/,iron-phorphyrins,/gi; } @testsneg; # synonyms of 
		map {s/,benzidine test positive.thirty-one percent,/,iron-phorphyrins,/gi; } @testsneg; # synonyms of 

		map {s/,naphthol phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphthol-as-bi-phosphate,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphthol-as-bi-phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphthol-as-bi-phosphohydrolase activities,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphthol-as-bi-phosphohydrolase activity,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphtol-as-bi-phosphohydrolase activities,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,naphtol-as-bi-phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,napthol-as-bi-phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphohydrolase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,weak naphthol-as-bi-phosphohydrolase activity,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 

		map {s/,moderate phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,phosphoamidase,/,naphthol-as-bi-phosphohydrolase,/gi; } @testsneg; # synonyms of 

		map {s/,indophenol oxidase,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,kovacs oxidase test,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,modified oxidase activity,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase activities,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase activity tests,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase activity,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase production,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase reaction,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase reactions,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase test,/,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase tests,/,oxidase,/gi; } @testsneg; # synonyms of 

		map {s/,collagenase activity,/,collagenase,/gi; } @testsneg; # synonyms of 

		map {s/,arylsulfatase production,/,arylsulfatase,/gi; } @testsneg; # synonyms of 

		map {s/,catalase activities,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase activity,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase production,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase reaction,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase reactions,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase test,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase tests,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase variable,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,catalase-positive and weakly oxidase,/,catalase,oxidase,/gi; } @testsneg; # synonyms of 
		map {s/,feeble catalase activity,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,h2s catalase,/,h2s,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,moderate catalase activity,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,occasionally very weak catalase activity,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,positive catalase reaction,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,positive catalase,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,produce catalase,/,catalase,/gi; } @testsneg; # synonyms of 
		map {s/,weak catalase activity,/,catalase,/gi; } @testsneg; # synonyms of 

		map {s/,produces urease,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,strongest urease reaction,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urea hydrolysis,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease activities,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease activity,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease production,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease reaction,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease results,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease test,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,urease tests,/,urease,/gi; } @testsneg; # synonyms of 

		map {s/,a-chymotrypsin activities,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,a-chymotrypsin,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-chymotrypsin,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,chymotrypsin activities,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,chymotrypsin activity,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,present and weak a-chymotrypsin activity,/,chymotrypsin,/gi; } @testsneg; # synonyms of 
		map {s/,trypsin activities,/,trypsin,/gi; } @testsneg; # synonyms of 

		map {s/,2-naphthyl-b-d-galactopyranoside,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,2-nitrophenyl-b-d-galactopyranoside,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-fucosidase activities,/,a-fucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-fucosidase activity,/,a-fucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-fucosidase production isnegative,/,not a-fucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-glucosidase activities,/,a-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-glucosidase activity,/,a-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-mannosidase activity,/,a-mannosidase,/gi; } @testsneg; # synonyms of 
		map {s/,a-mannosidase reaction,/,a-mannosidase,/gi; } @testsneg; # synonyms of 
		map {s/,addition b-galactosidase,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-galactosidase,/,a-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-glucosidase,/,a-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-mannosidase,/,a-mannosidase,/gi; } @testsneg; # synonyms of 
		map {s/,amylase plate test,/,amylase,/gi; } @testsneg; # synonyms of 
		map {s/,amylase test,/,amylase,/gi; } @testsneg; # synonyms of 
		map {s/,b-galactosidase 6-phosphatase,/,b-galactosidase 6-phosphate,/gi; } @testsneg; # synonyms of 
		map {s/,b-galactosidase activities,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-galactosidase activity,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-galactosidase arginine,/,b-galactosidase,arginine,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucoronidase,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucosidase activities,/,b-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucosidase activity,/,b-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucosidase reactions,/,b-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucuronidase activities,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucuronidase activity,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,b-glucuronidase test,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,beta-galactosidase,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,beta-glucosidase,/,b-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,beta-glucuronidase,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,beta-mannosidase,/,b-mannosidase,/gi; } @testsneg; # synonyms of 
		map {s/,biphosphate aldolase,/,fructose bisphosphate aldolase,/gi; } @testsneg; # synonyms of 
		map {s/,diphosphate aldolase,/,fructose bisphosphate aldolase,/gi; } @testsneg; # synonyms of 
		map {s/,enzymes b-glucosidase,/,b-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,exhibits a-glucosidase activity,/,a-glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,fucosidase activities,/,fucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,glucosidase activities,/,glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,glucosidase activity,/,glucosidase,/gi; } @testsneg; # synonyms of 
		map {s/,hippurate hydrolysis,/,hippurate,/gi; } @testsneg; # synonyms of 
		map {s/,hydrolyses xylan,/,xylanase,/gi; } @testsneg; # synonyms of 
		map {s/,including xylanase,/,xylanase,/gi; } @testsneg; # synonyms of 
		map {s/,liquefying amylase,/,amylase,/gi; } @testsneg; # synonyms of 
		map {s/,methyl alpha-d-mannoside,/,methyl a-d-mannoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl b-d-glucopyranoside,/,methyl b-d-glucoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl beta-d-xyloside,/,methyl b-d-xyloside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl-a-d-glucopyranoside,/,methyl a-d-glucoside,/gi; } @testsneg; # synonyms of 
		map {s/,methyl-a-d-glucoside,/,methyl a-d-glucoside,/gi; } @testsneg; # synonyms of 
		map {s/,moderate b-galactosidase activity,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,n-acetyl-beta-D-glucosaminidase,/,n-acetyl-b-d-glucosaminidase,/gi; } @testsneg; # synonyms of 
		map {s/,naphthol-as-bi-b-d-glucuronide,/,b-glucuronidase,/gi; } @testsneg; # synonyms of 
		map {s/,o-nitrophenyl-b-d-galactoside,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,onpg hydrolysis,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,onpg reaction,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,onpg test,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,onpg,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,para-nitro-d-methyl galactose,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,producing b-galactosidase,/,b-galactosidase,/gi; } @testsneg; # synonyms of 
		map {s/,saccharifying amylase,/,amylase,/gi; } @testsneg; # synonyms of 
		map {s/,weak amylase production,/,amylase,/gi; } @testsneg; # synonyms of 
		map {s/,xylan-inducible xylanase,/,xylanase,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 

		map {s/,deoxyribonuclease activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,deoxyribonuclease,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,dnase activities,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,dnase activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,moderate dnase activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,positive dnase activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,produce urease,/,urease,/gi; } @testsneg; # synonyms of 
		map {s/,strong dnase,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,urease arginine dihydrolase,/,urease,arginine dihydrolase,/gi; } @testsneg; # synonyms of 
		map {s/,weak deoxyribonuclease activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,weak dnase,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,weak extracellular deoxyribonuclease activity,/,dnase,/gi; } @testsneg; # synonyms of 
		map {s/,weakly positive dnase,/,dnase,/gi; } @testsneg; # synonyms of 

		map {s/,negative hyaluronidase activity,/,hyaluronidase,/gi; } @testsneg; # synonyms of 

		map {s/,weak caseinase activity,/,caseinase,/gi; } @testsneg; # synonyms of 
		map {s/,caseinase activities,/,caseinase,/gi; } @testsneg; # synonyms of 

		map {s/,tetrazolium,/,tetrazolium reduction,/gi; } @testsneg; # synonyms of 
		map {s/,tetrazolium red,/,tetrazolium reduction,/gi; } @testsneg; # synonyms of 
		map {s/,tetrazolium chloride,/,tetrazolium reduction,/gi; } @testsneg; # synonyms of 

		map {s/,lipase activities,/,lipase,/gi; } @testsneg; # synonyms of 
		map {s/,lipase activity,/,lipase,/gi; } @testsneg; # synonyms of 
		map {s/,lipase nor lecithinase,/,lipase,lecithinase,/gi; } @testsneg; # synonyms of 
		map {s/,lipase reaction,/,lipase,/gi; } @testsneg; # synonyms of 
		map {s/,lipolytic activity,/,lipase,/gi; } @testsneg; # synonyms of 
		map {s/,strong and true lipolytic activity,/,lipase,/gi; } @testsneg; # synonyms of 

		map {s/,nadp-dependent 6-phosphogluconate dehydrogenase,/,6-phosphogluconate dehydrogenase,/gi; } @testsneg; # synonyms of 

	#metabolic products
		map {s/,methyl red,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red test,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red tests,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red indicator,/,methyl red acidification,/gi; } @testsneg; # synonyms of 

		map {s/,acetoin production,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,acetoin,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,acetylmethylcarbinol production,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,acetylmethylcarbinol,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer reaction,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer reaction,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer reactions,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer reactions,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer test,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red and voges-proskauer test,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red reaction,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red reactions,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl red testing,/,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,methyl-red and voges-proskauer reactions,/,methyl red acidification,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,moderate acetylmethylcarbinol reactions,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,produces acetoin,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer and indole negative,/,acetoin,tryptophanase,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer and methylred reactions,/,acetoin,methyl red acidification,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer broth range,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer broth ranges,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer reaction,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer reactions,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer reagents,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer test range,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer test,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer tests,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,vp reaction,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,vp test,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,vp tests,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,vp,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,weak and voges-proskauer reaction,/,acetoin,/gi; } @testsneg; # synonyms of 
		map {s/,weak vp reaction,/,acetoin,/gi; } @testsneg; # synonyms of 

		map {s/,a-haemolysis,/,a-hemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,alpha hemolysis,/,a-hemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-haemolytic,/,a-hemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,alpha-hemolysis,/,a-hemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,b-hemolysis,/,b-haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,b-hemolysis,/,b-haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,beta-hemolytic,/,b-haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,blood haemolysis,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,double hemolysis zone,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,hemolysis activity,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,hemolysis,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,partial haemolysis,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,usually nonhemolytic,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,weak hemolysis activity,/,haemolysis,/gi; } @testsneg; # synonyms of 
		map {s/,weak hemolysis,/,haemolysis,/gi; } @testsneg; # synonyms of 

		map {s/,methylene blue,/,methylene blue reduction,/gi; } @testsneg; # synonyms of 

		map {s/,koh test,/,flexirubin,/gi; } @testsneg; # synonyms of 

		map {s/,indoxyl acetate reaction,/,indoxyl acetate hydrolysis,/gi; } @testsneg; # synonyms of 

		map {s/,tetrathionate,/,tetrathionate reduction,/gi; } @testsneg; # synonyms of 
		map {s/,resazurin,/,resazurin reduction,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,xxx,/gi; } @testsneg; # synonyms of 
#not a test
		map {s/,xxx,/,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,/gi; } @testsneg; # synonyms of 
		map {s/,xxx,/,/gi; } @testsneg; # synonyms of 
		map {s/,Gen,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d g yeast extract l-1or plus 10 mol glucose l-1,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d percent bile agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d-day buffered glucose peptone broth culture,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d ala,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d\d% agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d% agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d% nacl broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d.\d% nacl,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d\d% bile agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d% agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d% blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d% nacl,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d% sucrose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,\d% sucrose medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,0.36 mb agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,0.36 r2a agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,1 glu,/,/gi; } @testsneg; # synonyms of 
		map {s/,1% peptone-1% yeast extract-1% glucose broth cultures,/,/gi; } @testsneg; # synonyms of 
		map {s/,2,/,/gi; } @testsneg; # synonyms of 
		map {s/,2216e agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,5% horse-blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,acetate containing media,/,/gi; } @testsneg; # synonyms of 
		map {s/,actinomyces agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,aesculin bile agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar containing nacl,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar disk diffusion tests,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar media,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar medium containing bmc,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar roll tubes containing high numbers,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar roll tubes containing large numbers,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar shake tubes,/,/gi; } @testsneg; # synonyms of 
		map {s/,agar slopes containing either cellobiose,/,/gi; } @testsneg; # synonyms of 
		map {s/,ala,/,/gi; } @testsneg; # synonyms of 
		map {s/,alkaline nutrient agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,alkaline soybean agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,amd,/,/gi; } @testsneg; # synonyms of 
		map {s/,aminot acids glutamic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,and\/or eg agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,arthrobacter agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,azide-containing enterococcal selective agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bacitracin,/,/gi; } @testsneg; # synonyms of 
		map {s/,basal sml0 containing glucose,/,/gi; } @testsneg; # synonyms of 
		map {s/,betaine agar grown cells,/,/gi; } @testsneg; # synonyms of 
		map {s/,bhi agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bhi-blood agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,bile . esculin agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bile aesculin agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bile agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bile.aesculin agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bile.esculin agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,blood agar base no. 2,/,/gi; } @testsneg; # synonyms of 
		map {s/,blood agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,blood-containing agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,bovine blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,brain heart infusion agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,broth containing \d.\d% nacl,/,/gi; } @testsneg; # synonyms of 
		map {s/,broth containing d-fructose,/,/gi; } @testsneg; # synonyms of 
		map {s/,broth cultures containing glucose,/,/gi; } @testsneg; # synonyms of 
		map {s/,c agar medium colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,c agar medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,casein agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,caso agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,catalase negative,/,/gi; } @testsneg; # synonyms of 
		map {s/,cdc anaerobe 5% sheep blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,cell wall teichoic acid ontained glycerol,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellobiose broths,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellulose agar medium colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellulose agar medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellulose agar surface colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,cellulose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,cetrimide agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,cgy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,chapman agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,citrate agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,ckonkey agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,columbia agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,columbia blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,columbia sheep-blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,complex standard media acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,containing 0.006% phenol red,/,/gi; } @testsneg; # synonyms of 
		map {s/,containing 5% xylose,/,/gi; } @testsneg; # synonyms of 
		map {s/,cyc agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,e.g. blood trypticase soy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,eg agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,eg agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,egg.yolk agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,erythromycin,/,/gi; } @testsneg; # synonyms of 
		map {s/,facultatively aerobic and heterofermentative lactic acid bacterium,/,/gi; } @testsneg; # synonyms of 
		map {s/,ferment .*,/,/gi; } @testsneg; # synonyms of 
		map {s/,fermentative growth using pyruvate,/,/gi; } @testsneg; # synonyms of 
		map {s/,fermented.*,/,/gi; } @testsneg; # synonyms of 
		map {s/,fermenting.*,/,/gi; } @testsneg; # synonyms of 
		map {s/,ferments.*,/,/gi; } @testsneg; # synonyms of 
		map {s/,gel electrophoretic patterns,/,/gi; } @testsneg; # synonyms of 
		map {s/,glu,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose 0.85 butyrate.*,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose agar deep cultures,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose broths,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose media,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,glucose peptone-yeast extract broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,glutamate media,/,/gi; } @testsneg; # synonyms of 
		map {s/,glycerol media,/,/gi; } @testsneg; # synonyms of 
		map {s/,glycerol yields acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,gyp agar colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,gyp agar slant,/,/gi; } @testsneg; # synonyms of 
		map {s/,gyp agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,h2.co2 plus acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,h2.co2.1 mm acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,hemin,/,/gi; } @testsneg; # synonyms of 
		map {s/,horse-blood trypticase soy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,horse.blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,horseblood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,iron agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,isp 2 agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,kaa agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,kanamycin aesculin azide agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,kanamycin,/,/gi; } @testsneg; # synonyms of 
		map {s/,kanamycin.aesculin.azide agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,kf streptococcus agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,l-rhamnose broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,l-tyrosine-containing marine agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,lab lemco agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,lacking hemolysis,/,/gi; } @testsneg; # synonyms of 
		map {s/,lb agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,le ruyet et al.,/,/gi; } @testsneg; # synonyms of 
		map {s/,lead acetate agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,lead acetate paper,/,/gi; } @testsneg; # synonyms of 
		map {s/,least 1.5% agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,lincomycin,/,/gi; } @testsneg; # synonyms of 
		map {s/,lipid films,/,/gi; } @testsneg; # synonyms of 
		map {s/,litmus milk,/,/gi; } @testsneg; # synonyms of 
		map {s/,m.l. gen n. bromii,/,/gi; } @testsneg; # synonyms of 
		map {s/,M.L. gen n. bromii,/,/gi; } @testsneg; # synonyms of 
		map {s/,m.r.s. agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,m17 glucose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,macconkey agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,maltose amygdalin,/,/gi; } @testsneg; # synonyms of 
		map {s/,maltose broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,mannose containing lipid,/,/gi; } @testsneg; # synonyms of 
		map {s/,marine agar 2216 media,/,/gi; } @testsneg; # synonyms of 
		map {s/,marine agar 2216,/,/gi; } @testsneg; # synonyms of 
		map {s/,marine agar 2218,/,/gi; } @testsneg; # synonyms of 
		map {s/,meat agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,meat broth agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,media containing \d% glucose,/,/gi; } @testsneg; # synonyms of 
		map {s/,media containing nitrate,/,/gi; } @testsneg; # synonyms of 
		map {s/,medium containing glucose plus yeast extract,/,/gi; } @testsneg; # synonyms of 
		map {s/,mineral media containing pyruvate,/,/gi; } @testsneg; # synonyms of 
		map {s/,mineral medium containing 10 mm acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,minocycline,/,/gi; } @testsneg; # synonyms of 
		map {s/,mm10 blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,mm10 sheep blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,mrs agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,mrs agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,msg agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,mty agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,na containing 4% nacl,/,/gi; } @testsneg; # synonyms of 
		map {s/,nacl agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,nalidixic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,negative aesculin hydrolysis,/,/gi; } @testsneg; # synonyms of 
		map {s/,negative and aesculin hydrolysis,/,/gi; } @testsneg; # synonyms of 
		map {s/,negative and methyl red test,/,/gi; } @testsneg; # synonyms of 
		map {s/,negative nitrate reduction,/,/gi; } @testsneg; # synonyms of 
		map {s/,negative voges-proskauer reaction,/,/gi; } @testsneg; # synonyms of 
		map {s/,nitrate reduction-negative,/,/gi; } @testsneg; # synonyms of 
		map {s/,non-haemolytic,/,/gi; } @testsneg; # synonyms of 
		map {s/,non-hemolytic,/,/gi; } @testsneg; # synonyms of 
		map {s/,nonhemolytic,/,/gi; } @testsneg; # synonyms of 
		map {s/,novobiocin,/,/gi; } @testsneg; # synonyms of 
		map {s/,nutrient agar colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,nutrient agar cultures,/,/gi; } @testsneg; # synonyms of 
		map {s/,nutrient agar growth,/,/gi; } @testsneg; # synonyms of 
		map {s/,nutrient agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,nutrient agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,old culture,/,/gi; } @testsneg; # synonyms of 
		map {s/,old cultures,/,/gi; } @testsneg; # synonyms of 
		map {s/,old,/,/gi; } @testsneg; # synonyms of 
		map {s/,ox gall,/,/gi; } @testsneg; # synonyms of 
		map {s/,oxidase negative,/,/gi; } @testsneg; # synonyms of 
		map {s/,p agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,parasites-of man,/,/gi; } @testsneg; # synonyms of 
		map {s/,pc skim milk agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,pea azide agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,peptone yeast glucose broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,peptone-yeast broth medium containing glucose,/,/gi; } @testsneg; # synonyms of 
		map {s/,peptone-yeast extract-glucose deep agar cultures,/,/gi; } @testsneg; # synonyms of 
		map {s/,ph in voges-proskauer broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatidyl ethanolamine,/,/gi; } @testsneg; # synonyms of 
		map {s/,phosphatidyl glycerol,/,/gi; } @testsneg; # synonyms of 
		map {s/,positive and .*,/,/gi; } @testsneg; # synonyms of 
		map {s/,potato agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,potato dextrose agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,potato-glucose agar media,/,/gi; } @testsneg; # synonyms of 
		map {s/,ppyg agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,pras media+ 0.5% tween 80,/,/gi; } @testsneg; # synonyms of 
		map {s/,predominantly l-lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,predominately l-lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,primarily lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,producing acetate,/,/gi; } @testsneg; # synonyms of 
		map {s/,producing exclusively lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,producing l-lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,producing lactic acid,/,/gi; } @testsneg; # synonyms of 
		map {s/,py agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,py-1% glucose medium ph 6.8,/,/gi; } @testsneg; # synonyms of 
		map {s/,pye agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,pyi agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,pyt broth medium containing esculin,/,/gi; } @testsneg; # synonyms of 
		map {s/,r2a agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,raffinose liquid media,/,/gi; } @testsneg; # synonyms of 
		map {s/,rgm-rumen fluid plus glucose,/,/gi; } @testsneg; # synonyms of 
		map {s/,salmonella-shigella agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,schaedler sheep blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,semi-solid agar medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,sheep blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,sheep-blood agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,skimmed-milk agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,slanetz & bartley agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,slanetz-bartley agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,soft 0.5% mrs agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,solid agar containing cellulose,/,/gi; } @testsneg; # synonyms of 
		map {s/,spot,/,/gi; } @testsneg; # synonyms of 
		map {s/,standard mrs agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,swcm agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,tgy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,thallous acetate-tetrazolium agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,thermophilic and hydrolyses agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,thin agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,todd-hewitt agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,troleandomycin,/,/gi; } @testsneg; # synonyms of 
		map {s/,trypticase soy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,trypticase soy broth agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,trypticase soy sheep blood agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,trypticase-soy-sheep blood agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,tryptone soy agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,tryptone soya agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,tryptose soya agar colonies,/,/gi; } @testsneg; # synonyms of 
		map {s/,tsa agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,tsb agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,tween 80 agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,vibriostatic agent o.129,/,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer broth,/,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,voges-proskauer negative,/,/gi; } @testsneg; # synonyms of 
		map {s/,vy\/2 agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,wet agar plates,/,/gi; } @testsneg; # synonyms of 
		map {s/,wilkins agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,xr2a agarxx,/,/gi; } @testsneg; # synonyms of 
		map {s/,xylan agar medium,/,/gi; } @testsneg; # synonyms of 
		map {s/,xylan agar,/,/gi; } @testsneg; # synonyms of 
		map {s/,ycfa medium containing 0.2% glucose,/,/gi; } @testsneg; # synonyms of 

		print OUT @taxlabels, @testsneg, "\n"; # prints to $homout, hom.testsneg.txt
		}	

#character discovery - puts all the characters in a single line, gets rid of "not"s in characters
	my $temp2 = "temp2.testsneg.txt";
	open (IN, '<', $homout) or die $!; # opens up the list of homologized terms
	open (OUT, '>', $temp2) or die $!; 
	local $, = "\t";	
	while (my $line = <IN> ) { # pushes the elements into an array, sorts them, retains only the unique ones
		chomp $line;
		my @unsortedlist = split /\t/, $line;
		push (my @homcharlist, $unsortedlist[1]);
		map {s/,not /,/g; } @homcharlist; # gets rid of the word "not" in the beginning of characters
		map {s/tests negative//g; } @homcharlist; # gets rid of the character name label
		map {s/^,//g; } @homcharlist; # gets rid of the comma at the beginning
		map {s/,$//g; } @homcharlist; # gets rid of the comma at the end
		map {s/,/\t/g; } @homcharlist; # converts commas to tabs
		print OUT @homcharlist, "\t"; # prints to $temp2, temp2.testsneg.txt
		}
#character discovery -sorts the characters and finds the unique characters
	my $p = 1;
	my $m = 1;
	my $temp3 = "temp3.testsneg.txt";	
	open (IN, '<', $temp2) or die $!;
	open (OUT, '>', $temp3) or die $!;
	my $line = <IN>;
	chomp $line;
	$line =~ s/\t\t/\t/g;
	$line =~ s/$/\t/;
	my @values = split /\t/, $line;
	my @filtered = uniq(@values);
	@filtered = sort(@filtered);
	print OUT @filtered; # prints to $temp3, temp3.testsneg.txt
#character discovery -prints out the homologized characters
	my $r = 1;
	my $temp4 = "temp4.testsneg.txt";
	open (IN, '<', $temp3) or die $!;
	open (OUT, '>', $temp4) or die $!;
	$line = <IN>;
	chomp $line;
	$line =~ s/^\t//;
	my @charlist2 = split (/\t/,$line);
	print OUT "@charlist2", "\n"; # prints to $temp4 temp4.testsneg.txt
#temporarily rename charstates to label those that have been homologized		
	#substrates used and/or fermented
	map {s/^2-aminoethanol/**2-aminoethanol/gi } @charlist2;#1
	map {s/^2-keto-d-gluconate/**2-keto-d-gluconate/gi } @charlist2;#2
	map {s/^2-ketobutyrate/**2-ketobutyrate/gi } @charlist2;#3
	map {s/^2-3-butanediol/**2-3-butanediol/gi } @charlist2;#4
	map {s/^1-3-butanediol/**1-3-butanediol/gi } @charlist2;#5
	map {s/^3-hydroxybenzoate/**3-hydroxybenzoate/gi } @charlist2;#6
	map {s/^3-hydroxybutyrate/**3-hydroxybutyrate/gi } @charlist2;#7
	map {s/^3-methyl d-glucose/**3-methyl d-glucose/gi } @charlist2;#8
	map {s/^4-hydroxybenzoate/**4-hydroxybenzoate/gi } @charlist2;#9
	map {s/^4-hydroxyphenylacetate/**4-hydroxyphenylacetate/gi } @charlist2;#10
	map {s/^5-keto-d-gluconate/**5-keto-d-gluconate/gi } @charlist2;#12
	map {s/^a-ketoglutarate/**a-ketoglutarate/gi } @charlist2;#13
	map {s/^a-ketovalerate/**a-ketovalerate/gi } @charlist2;#14
	map {s/^acetate/**acetate/gi } @charlist2;#15
	map {s/^acetoacetate/**acetoacetate/gi } @charlist2;#16
	map {s/^adenosine 5\'-monophosphate/**adenosine 5\'-monophosphate/gi } @charlist2;#17
	map {s/^adenosine/**adenosine/gi } @charlist2;#18
	map {s/^adipate/**adipate/gi } @charlist2;#19
	map {s/^adonitol/**adonitol/gi } @charlist2;#20
	map {s/^agar/**agar/gi } @charlist2;#21
	map {s/^alaninamide/**alaninamide/gi } @charlist2;#22
	map {s/^alanine/**alanine/gi } @charlist2;#23
	map {s/^alanyl glycine/**alanyl glycine/gi } @charlist2;#24
	map {s/^amygdalin/**amygdalin/gi } @charlist2;#25
	map {s/^arabinose/**arabinose/gi } @charlist2;#26
	map {s/^arabitol/**arabitol/gi } @charlist2;#27
	map {s/^arbutin/**arbutin/gi } @charlist2;#28
	map {s/^arginine/**arginine/gi } @charlist2;#29
	map {s/^asparagine/**asparagine/gi } @charlist2;#30
	map {s/^aspartate/**aspartate/gi } @charlist2;#31
	map {s/^b-cyclodextrin/**b-cyclodextrin/gi } @charlist2;#32
	map {s/^b-xylosidase/**b-xylosidase/gi } @charlist2;#34
	map {s/^butyrate/**butyrate/gi } @charlist2;#35
	map {s/^caprate/**caprate/gi } @charlist2;#36
	map {s/^carbohydrates/**carbohydrates/gi } @charlist2;#37
	map {s/^carboxymethyl cellulose/**carboxymethyl cellulose/gi } @charlist2;#38
	map {s/^carnitine/**carnitine/gi } @charlist2;#39
	map {s/^cellobiose/**cellobiose/gi } @charlist2;#40
	map {s/^cis-aconitate/**cis-aconitate/gi } @charlist2;#41
	map {s/^citrate/**citrate/gi } @charlist2;#42
	map {s/^co2/**co2/gi } @charlist2;#43
	map {s/^cyclodextrin/**cyclodextrin/gi } @charlist2;#44
	map {s/^dextrin/**dextrin/gi } @charlist2;#45
	map {s/^dulcitol/**dulcitol/gi } @charlist2;#46
	map {s/^erythritol/**erythritol/gi } @charlist2;#47
	map {s/^ethanol/**ethanol/gi } @charlist2;#48
	map {s/^formate/**formate/gi } @charlist2;#49
	map {s/^fructose 6-phosphate/**fructose 6-phosphate/gi } @charlist2;#50
	map {s/^fructose/**fructose/gi } @charlist2;#51
	map {s/^fucose/**fucose/gi } @charlist2;#52
	map {s/^fumarate/**fumarate/gi } @charlist2;#53
	map {s/^fusidic acid/**fusidic acid/gi } @charlist2;#54
	map {s/^galactose/**galactose/gi } @charlist2;#55
	map {s/^galacturonate/**galacturonate/gi } @charlist2;#56
	map {s/^gelatin/**gelatin/gi } @charlist2;#57
	map {s/^gentiobiose/**gentiobiose/gi } @charlist2;#58
	map {s/^gluconate/**gluconate/gi } @charlist2;#59
	map {s/^glucose 1-phosphate/**glucose 1-phosphate/gi } @charlist2;#60
	map {s/^glucose/**glucose/gi } @charlist2;#61
	map {s/^glucuronate/**glucuronate/gi } @charlist2;#62
	map {s/^glutamate/**glutamate/gi } @charlist2;#63
	map {s/^glutamic acid polymer/**glutamic acid polymer/gi } @charlist2;#64
	map {s/^glutamine/**glutamine/gi } @charlist2;#65
	map {s/^glutamyl glutamic acid/**glutamyl glutamic acid/gi } @charlist2;#66
	map {s/^glycerol phosphate/**glycerol phosphate/gi } @charlist2;#68
	map {s/^glycerol/**glycerol/gi } @charlist2;#69
	map {s/^glycogen/**glycogen/gi } @charlist2;#70
	map {s/^glycyl l-glutamic acid/**glycyl l-glutamic acid/gi } @charlist2;#71
	map {s/^glycyl-l-proline/**glycyl-l-proline/gi } @charlist2;#72
	map {s/^glyoxylate/**glyoxylate/gi } @charlist2;#73
	map {s/^guanidine/**guanidine/gi } @charlist2;#74
	map {s/^h2s/**h2s/gi } @charlist2;#75
	map {s/^hippurate/**hippurate/gi } @charlist2;#76
	map {s/^histidine/**histidine/gi } @charlist2;#77
	map {s/^indole/**indole/gi } @charlist2;#78
	map {s/^inosine/**inosine/gi } @charlist2;#79
	map {s/^inositol/**inositol/gi } @charlist2;#80
	map {s/^inulin/**inulin/gi } @charlist2;#81
	map {s/^itaconate/**itaconate/gi } @charlist2;#82
	map {s/^lactate/**lactate/gi } @charlist2;#83
	map {s/^lactic acid methyl ester/**lactic acid methyl ester/gi } @charlist2;#84
	map {s/^lactose/**lactose/gi } @charlist2;#85
	map {s/^lactulose/**lactulose/gi } @charlist2;#86
	map {s/^leucine/**leucine/gi } @charlist2;#87
	map {s/^lyxose/**lyxose/gi } @charlist2;#88
	map {s/^malate/**malate/gi } @charlist2;#89
	map {s/^malonate/**malonate/gi } @charlist2;#90
	map {s/^maltose/**maltose/gi } @charlist2;#92
	map {s/^maltotriose/**maltotiose/gi } @charlist2;#93
	map {s/^mannan/**mannan/gi } @charlist2;#94
	map {s/^mannitol/**mannitol/gi } @charlist2;#95
	map {s/^mannose/**mannose/gi } @charlist2;#96
	map {s/^melezitose/**melezitose/gi } @charlist2;#97
	map {s/^melibiose/**melibiose/gi } @charlist2;#98
	map {s/^meso-tartrate/**meso-tartrate/gi } @charlist2;#99
	map {s/^methionine/**methionine/gi } @charlist2;#100
	map {s/^methyl a-d-galactoside/**methyl a-d-galactoside/gi } @charlist2;#101
	map {s/^methyl a-d-glucoside/**methyl a-d-glucoside/gi } @charlist2;#102
	map {s/^methyl b-d-galactoside/**methyl b-d-galactoside/gi } @charlist2;#103
	map {s/^methyl b-d-glucoside/**methyl b-d-glucoside/gi } @charlist2;#104
	map {s/^methyl pyruvate/**methyl pyruvate/gi } @charlist2;#105
	map {s/^monomethyl succinate/**monomethyl succinate/gi } @charlist2;#106
	map {s/^mucate/**mucate/gi } @charlist2;#107
	map {s/^myo-inositol/**myo-inositol/gi } @charlist2;#108
	map {s/^n-acetyl-b-d-glucosamine/**n-acetyl-b-d-glucosamine/gi } @charlist2;#109
	map {s/^n-acetyl-b-d-mannosamine/**n-acetyl-b-d-mannosamine/gi } @charlist2;#110
	map {s/^n-acetyl-l-glutamate/**n-acetyl-l-glutamate/gi } @charlist2;#111
	map {s/^ornithine/**ornithine/gi } @charlist2;#112
	map {s/^palatinose/**palatinose/gi } @charlist2;#113
	map {s/^pectin/**pectin/gi } @charlist2;#114
	map {s/^phenylacetate/**phenylacetate/gi } @charlist2;#115
	map {s/^phenylalanine/**phenylalanine/gi } @charlist2;#116
	map {s/^phenylethylamine/**phenylethylamine/gi } @charlist2;#117
	map {s/^proline/**proline/gi } @charlist2;#118
	map {s/^propionate/**propionate/gi } @charlist2;#119
	map {s/^psicose/**psicose/gi } @charlist2;#120
	map {s/^pullulan/**pullulan/gi } @charlist2;#121
	map {s/^putrescine/**putrescine/gi } @charlist2;#122
	map {s/^pyroglutamate/**pyroglutamate/gi } @charlist2;#123
	map {s/^pyruvate/**pyruvate/gi } @charlist2;#124
	map {s/^quinate/**quinate/gi } @charlist2;#126
	map {s/^raffinose/**raffinose/gi } @charlist2;#127
	map {s/^rhamnose/**rhamnose/gi } @charlist2;#128
	map {s/^ribitol/**ribitol/gi } @charlist2;#129
	map {s/^ribose/**ribose/gi } @charlist2;#130
	map {s/^saccharolytic/**saccharolytic/gi } @charlist2;#131
	map {s/^salicin/**salicin/gi } @charlist2;#133
	map {s/^sedoheptulosan/**sedoheptulosan/gi } @charlist2;#134
	map {s/^serine/**serine/gi } @charlist2;#135
	map {s/^sorbitol/**sorbitol/gi } @charlist2;#136
	map {s/^sorbose/**sorbose/gi } @charlist2;#137
	map {s/^stachyose/**stachyose/gi } @charlist2;#138
	map {s/^starch/**starch/gi } @charlist2;#139
	map {s/^suberate/**suberate/gi } @charlist2;#140
	map {s/^succinamic acid/**succinamic acid/gi } @charlist2;#141
	map {s/^succinate/**succinate/gi } @charlist2;#142
	map {s/^cysteine/**cysteine/gi } @charlist2;#143
	map {s/^sucrose/**sucrose/gi } @charlist2;#144
	map {s/^tagatose/**tagatose/gi } @charlist2;#145
	map {s/^tartrate/**tartrate/gi } @charlist2;#146
	map {s/^tellurite/**tellurite/gi } @charlist2;#147
	map {s/^threonine/**threonine/gi } @charlist2;#148
	map {s/^thymidine 5\'-monophosphate/**thymidine 5\'-monophosphate/gi } @charlist2;#149
	map {s/^thymidine/**thymidine/gi } @charlist2;#150
	map {s/^trehalose/**trehalose/gi } @charlist2;#151
	map {s/^glycine/**glycine/gi } @charlist2;#152
	map {s/^turanose/**turanose/gi } @charlist2;#153
	map {s/^tween 40/**tween 40/gi } @charlist2;#154
	map {s/^tween 80/**tween 80/gi } @charlist2;#155
	map {s/^uridine/**uridine/gi } @charlist2;#156
	map {s/^valerate/**valerate/gi } @charlist2;#157
	map {s/^valine/**valine/gi } @charlist2;#158
	map {s/^xylan/**xylan/gi } @charlist2;#159
	map {s/^xylitol/**xylitol/gi } @charlist2;#160
	map {s/^xylose/**xylose/gi } @charlist2;#161

	map {s/^isoleucine/**isoleucine/gi } @charlist2;#162
	map {s/^lysine/**lysine/gi } @charlist2;#163
	map {s/^tryptophan/**tryptophan/gi } @charlist2;#164
	map {s/^selenocysteine/**selenocysteine/gi } @charlist2;#165
	map {s/^pyrrolysine/**pyrrolysine/gi } @charlist2;#166

	#inorganics
	map {s/^nitrate reduction/**nitrate reduction/gi } @charlist2;#167
	map {s/^nitrogen fixation/**nitrogen fixation/gi } @charlist2;#168
	map {s/^nitrite reduction/**nitrite reduction/gi } @charlist2;#169

	#enzymes present or enzymatic substrates
	map {s/^6-phospho-b-galactosidase/**6-phospho-b-galactosidase/gi } @charlist2;#170
	map {s/^6-phosphogluconate dehydrogenase/**6-phosphogluconate dehydrogenase/gi } @charlist2;#171
	map {s/^a-arabinosidase/**a-arabinosidase/gi } @charlist2;#172
	map {s/^a-fucosidase/**a-fucosidase/gi } @charlist2;#173
	map {s/^a-galactosidase/**a-galactosidase/gi } @charlist2;#174
	map {s/^a-glucosidase/**a-glucosidase/gi } @charlist2;#175
	map {s/^a-mannosidase/**a-mannosidase/gi } @charlist2;#176
	map {s/^a-methyl glucosidase/**a-methyl glucosidase/gi } @charlist2;#177
	map {s/^a-methylmannoside/**a-methylmannoside/gi } @charlist2;#178
	map {s/^acid phosphatase/**acid phosphatase/gi } @charlist2;#179
	map {s/^aesculin hydrolysis/**aesculin hydrolysis/gi } @charlist2;#180
	map {s/^alanine arylamidase/**alanine arylamidase/gi } @charlist2;#181
	map {s/^alkaline phosphatase/**alkaline phosphatase/gi } @charlist2;#182
	map {s/^alpha-glutamate arylamidase/**alpha-glutamate arylamidase/gi } @charlist2;#183
	map {s/^amylase/**amylase/gi } @charlist2;#184
	map {s/^arginine arylamidase/**arginine arylamidase/gi } @charlist2;#185
	map {s/^arginine decarboxylase/**arginine decarboxylase/gi } @charlist2;#186
	map {s/^arginine dihydrolase/**arginine dihydrolase/gi } @charlist2;#187
	map {s/^arylsulfatase/**arylsulfatase/gi } @charlist2;#188
	map {s/^asparagine arylamidase/**asparagine arylamidase/gi } @charlist2;#189
	map {s/^aspartate arylamidase/**aspartate arylamidase/gi } @charlist2;#190
	map {s/^b-galactosidase/**b-galactosidase/gi } @charlist2;#191
	map {s/^b-glucosidase/**b-glucosidase/gi } @charlist2;#192
	map {s/^b-glucuronidase/**b-glucuronidase/gi } @charlist2;#193
	map {s/^b-mannosidase/**b-mannosidase/gi } @charlist2;#194
	map {s/^casein hydrolysis/**casein hydrolysis/gi } @charlist2;#195
	map {s/^caseinase/**caseinase/gi } @charlist2;#196
	map {s/^catalase/**catalase/gi } @charlist2;#197
	map {s/^cellobiose/**cellobiose/gi } @charlist2;#198
	map {s/^chymotrypsin/**chymotrypsin/gi } @charlist2;#199
	map {s/^coagulase/**coagulase/gi } @charlist2;#200
	map {s/^collagenase/**collagenase/gi } @charlist2;#201
	map {s/^cysteine arylamidase/**cysteine arylamidase/gi } @charlist2;#202
	map {s/^cystine arylamidase/**cystine arylamidase/gi } @charlist2;#203
	map {s/^cytochrome oxidase/**cytochrome oxidase/gi } @charlist2;#204
	map {s/^dnase/**dnase/gi } @charlist2;#205
	map {s/^esterase c4/**esterase c4/gi } @charlist2;#206
	map {s/^fructose bisphosphate aldolase/**fructose bisphosphate aldolase/gi } @charlist2;#207
	map {s/^fucosidase/**fucosidase/gi } @charlist2;#208
	map {s/^gelatinase/**gelatinase/gi } @charlist2;#209
	map {s/^glucosidase/**glucosidase/gi } @charlist2;#210
	map {s/^glutamate arylamidase/**glutamate arylamidase/gi } @charlist2;#211
	map {s/^glutamate deaminase/**glutamate deaminase/gi } @charlist2;#212
	map {s/^glutamic acid decarboxylase/**glutamic acid decarboxylase/gi } @charlist2;#213
	map {s/^glutamine arylamidase/**glutamine arylamidase/gi } @charlist2;#214
	map {s/^glutamyl glutamic acid arylamidase/**glutamyl glutamic acid arylamidase/gi } @charlist2;#215
	map {s/^glycine arylamidase/**glycine arylamidase/gi } @charlist2;#216
	map {s/^hippurate hydrolysis/**hippurate hydrolysis/gi } @charlist2;#217
	map {s/^histidine arylamidase/**histidine arylamidase/gi } @charlist2;#218
	map {s/^hyaluronidase/**hyaluronidase/gi } @charlist2;#219
	map {s/^hydroxy proline arylamidase/**hydroxy proline arylamidase/gi } @charlist2;#220
	map {s/^iron-phorphyrins/**iron-phorphyrins/gi } @charlist2;#221
	map {s/^isoleucine arylamidase/**isoleucine arylamidase/gi } @charlist2;#222
	map {s/^l-lactic acid racemase/**l-lactic acid racemase/gi } @charlist2;#223
	map {s/^lecithinase/**lecithinase/gi } @charlist2;#224
	map {s/^leucine arylamidase/**leucine arylamidase/gi } @charlist2;#225
	map {s/^leucyl glycine arylamidase/**leucyl glycine arylamidase/gi } @charlist2;#226
	map {s/^lipase c14/**lipase c14/gi } @charlist2;#227
	map {s/^lipase c8/**lipase c8/gi } @charlist2;#228
	map {s/^lipase/**lipase/gi } @charlist2;#229
	map {s/^lysine arylamidase/**lysine arylamidase/gi } @charlist2;#230
	map {s/^lysine deaminase/**lysine deaminase/gi } @charlist2;#231
	map {s/^lysine decarboxylase/**lysine decarboxylase/gi } @charlist2;#232
	map {s/^methionine arylamidase/**methionine arylamidase/gi } @charlist2;#233
	map {s/^methyl a-d-mannoside/**methyl a-d-mannoside/gi } @charlist2;#234
	map {s/^methyl b-d-xyloside/**methyl b-d-xyloside/gi } @charlist2;#235
	map {s/^n-acetyl-b-d-glucosaminidase/**n-acetyl-b-d-glucosaminidase/gi } @charlist2;#236
	map {s/^nadh oxidase/**nadh oxidase/gi } @charlist2;#237
	map {s/^naphthol-as-bi-phosphohydrolase/**naphthol-as-bi-phosphohydrolase/gi } @charlist2;#238
	map {s/^ornithine arylamidase/**ornithine arylamidase/gi } @charlist2;#239
	map {s/^ornithine decarboxylase/**ornithine decarboxylase/gi } @charlist2;#240
	map {s/^oxidase/**oxidase/gi } @charlist2;#241
	map {s/^para-phenylalanine deaminase/**para-phenylalanine deaminase/gi } @charlist2;#242
	map {s/^phenylalanine arylamidase/**phenylalanine arylamidase/gi } @charlist2;#243
	map {s/^phenylalanine deaminase/**phenylalanine deaminase/gi } @charlist2;#244
	map {s/^phosphatase/**phosphatase/gi } @charlist2;#245
	map {s/^proline arylamidase/**proline arylamidase/gi } @charlist2;#246
	map {s/^protease/**protease/gi } @charlist2;#247
	map {s/^pyrazinamidase/**pyrazinamidase/gi } @charlist2;#248
	map {s/^pyroglutamic acid arylamidase/**pyroglutamic acid arylamidase/gi } @charlist2;#249

	map {s/^serine arylamidase/**serine arylamidase/gi } @charlist2;#251
	map {s/^threonine arylamidase/**threonine arylamidase/gi } @charlist2;#252
	map {s/^trypsin/**trypsin/gi } @charlist2;#253
	map {s/^tryptophan arylamidase/**tryptophan arylamidase/gi } @charlist2;#254
	map {s/^tryptophan deaminase/**tryptophan deaminase/gi } @charlist2;#255
	map {s/^tryptophanase/**tryptophanase/gi } @charlist2;#256
	map {s/^tyrosine arylamidase/**tyrosine arylamidase/gi } @charlist2;#257
	map {s/^urease/**urease/gi } @charlist2;#258
	map {s/^valine arylamidase/**valine arylamidase/gi } @charlist2;#259
	map {s/^xylanase/**xylanase/gi } @charlist2;#260

	#metabolic products
	map {s/^a-hemolysis/**a-hemolysis/gi } @charlist2;#261
	map {s/^acetoin/**acetoin/gi } @charlist2;#262
	map {s/^b-haemolysis/**b-haemolysis/gi } @charlist2;#263
	map {s/^ellman/**ellman/gi } @charlist2;#264
	map {s/^flexirubin/**flexirubin/gi } @charlist2;#265
	map {s/^haemolysis/**haemolysis/gi } @charlist2;#266
	map {s/^indoxyl acetate hydrolysis/**indoxyl acetate hydrolysis/gi } @charlist2;#267
	map {s/^methyl red acidification/**methyl red acidification/gi } @charlist2;#268
	map {s/^methylene blue reduction/**methylene blue reduction/gi } @charlist2;#269
	map {s/^resazurin reduction/**resazurin reduction/gi } @charlist2;#270
	map {s/^tetrathionate reduction/**tetrathionate reduction/gi } @charlist2;#271
	map {s/^tetrazolium reduction/**tetrazolium reduction/gi } @charlist2;#272
	map {s/^carboxymethyl cellulase/**carboxymethyl cellulase/gi } @charlist2;#273

	map {s/^1-3-propanediol/**1-3-propanediol/gi } @charlist2; #274
	map {s/^1-2-propanediol/**1-2-propanediol/gi } @charlist2; #275
	map {s/^propanediol/**propanediol/gi } @charlist2; #276

	map {s/^phosphatidylinositol phospholipase c/**phosphatidylinositol phospholipase c/gi } @charlist2;#277
	map {s/^lecithinase/**lecithinase/gi } @charlist2;#278
	map {s/^lactamide/**lactamide/gi } @charlist2;#279
	map {s/^hydroxyproline/**hydroxyproline/gi } @charlist2;#280
	map {s/^glycyl tryptophan arylamidase/**glycyl tryptophan arylamidase/gi } @charlist2;#281

	map {s/^esterase/**esterase/gi } @charlist2;#283


	print "\n\nBelow is your list of homologized characters states for the character $character:\n";
	print "Characters with ** have been homologized.  Those without ** have to be added to the perl script using map statements.\n\n";
	foreach (@charlist2) {
		print $r++;
		print " $_\n";
		}
	print "\n";		

#prepare for coding characters by removing duplicate homologized characters
	my $temp5 = "temp5.testsneg.txt";
	open (IN, '<', $homout) or die $!;
	open (OUT, '>', $temp5) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/\t,/\t/g;
		$line =~ s/,\t//g;
		$line =~ s/\t/,/g;
		my @charstates = split (/,/, $line);
		my @filteredstates = uniq(@charstates);
#		map {s/\t/,/g; } @filteredstates; 
		local $, = ",";
		print OUT @filteredstates, "\n";# prints to $temp5 temp5.testsneg.txt
		}
	my $temp6 = "temp6.testsneg.txt";
	open (IN, '<', $temp5) or die $!;
	open (OUT, '>', $temp6) or die $!;
	while ($line = <IN>) {
		chomp $line;
		$line =~ s/,/\t,/;
		print OUT $line, "\n"; # prints to $temp6 temp6.testsneg.txt
		}
	
#prepare nexus file
	my @taxnames;
	my $nexusoutfile = "testsneg.nex";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		my @testsnegdata = split (/\t/, $line);
		push (@taxnames, $testsnegdata[0]);
		}
	my $numtax = scalar(@taxnames) - 1;
	print OUT "#NEXUS\n\nBEGIN TAXA\;\n\tTITLE Taxa\;\n\tDIMENSIONS NTAX=$numtax\;\n\tTAXLABELS\n";
	shift @taxnames;
	local $, = " ";
	print OUT "\t\t", @taxnames, "\n" ;
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Tests Negative Matrix'\;\n\tDIMENSIONS NCHAR=283\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
	print OUT "CHARSTATELABELS\n\t\t";
	print OUT "1 '2-aminoethanol' \/  'no 2-aminoethanol' '2-aminoethanol', ";
	print OUT "2 '2-keto-d-gluconate' \/  'no 2-keto-d-gluconate' '2-keto-d-gluconate', ";
	print OUT "3 '2-ketobutyrate' \/  'no 2-ketobutyrate' '2-ketobutyrate', ";
	print OUT "4 '2,3-butanediol' \/  'no 2,3-butanediol' '2,3-butanediol', ";
	print OUT "5 '1,3-butanediol' \/  'no 1,3-butanediol' '1,3-butanediol', ";
	print OUT "6 '3-hydroxybenzoate' \/  'no 3-hydroxybenzoate' '3-hydroxybenzoate', ";
	print OUT "7 '3-hydroxybutyrate' \/  'no 3-hydroxybutyrate' '3-hydroxybutyrate', ";
	print OUT "8 '3-methyl d-glucose' \/  'no 3-methyl d-glucose' '3-methyl d-glucose', ";
	print OUT "9 '4-hydroxybenzoate' \/  'no 4-hydroxybenzoate' '4-hydroxybenzoate', ";
	print OUT "10 '4-hydroxyphenylacetate' \/  'no 4-hydroxyphenylacetate' '4-hydroxyphenylacetate', ";

	print OUT "11 'alcohols' \/  'no alcohols' 'alcohols', ";

	print OUT "12 '5-keto-d-gluconate' \/  'no 5-keto-d-gluconate' '5-keto-d-gluconate', ";
	print OUT "13 'a-ketoglutarate' \/  'no a-ketoglutarate' 'a-ketoglutarate', ";
	print OUT "14 'a-ketovalerate' \/  'no a-ketovalerate' 'a-ketovalerate', ";
	print OUT "15 'acetate' \/  'no acetate' 'acetate', ";
	print OUT "16 'acetoacetate' \/  'no acetoacetate' 'acetoacetate', ";
	print OUT "17 'adenosine 5prime-monophosphate' \/  'no adenosine 5prime-monophosphate' 'adenosine 5prime-monophosphate', ";
	print OUT "18 'adenosine' \/  'no adenosine' 'adenosine', ";
	print OUT "19 'adipate' \/  'no adipate' 'adipate', ";

	print OUT "20 'organonitrogen compounds' \/  'no organonitrogen compounds' 'organonitrogen compounds', ";

	print OUT "21 'agar' \/  'no agar' 'agar', ";
	print OUT "22 'alaninamide' \/  'no alaninamide' 'alaninamide', ";
	print OUT "23 'alanine' \/  'no alanine' 'alanine', ";
	print OUT "24 'alanyl glycine' \/  'no alanyl glycine' 'alanyl glycine', ";
	print OUT "25 'amygdalin' \/  'no amygdalin' 'amygdalin', ";
	print OUT "26 'arabinose' \/  'no arabinose' 'arabinose', ";
	print OUT "27 'arabitol' \/  'no arabitol' 'arabitol', ";
	print OUT "28 'arbutin' \/  'no arbutin' 'arbutin', ";
	print OUT "29 'arginine' \/  'no arginine' 'arginine', ";
	print OUT "30 'asparagine' \/  'no asparagine' 'asparagine', ";
	print OUT "31 'aspartate' \/  'no aspartate' 'aspartate', ";
	print OUT "32 'b-cyclodextrin' \/  'no b-cyclodextrin' 'b-cyclodextrin', ";

	print OUT "33 'fatty acids' \/  'no fatty acids' 'fatty acids', ";

	print OUT "34 'b-xylosidase' \/  'no b-xylosidase' 'b-xylosidase', ";
	print OUT "35 'butyrate' \/  'no butyrate' 'butyrate', ";
	print OUT "36 'caprate' \/  'no caprate' 'caprate', ";
	print OUT "37 'carbohydrates' \/  'no carbohydrates' 'carbohydrates', ";
	print OUT "38 'carboxymethyl cellulose' \/  'no carboxymethyl cellulose' 'carboxymethyl cellulose', ";
	print OUT "39 'carnitine' \/  'no carnitine' 'carnitine', ";
	print OUT "40 'cellobiose' \/  'no cellobiose' 'cellobiose', ";
	print OUT "41 'cis-aconitate' \/  'no cis-aconitate' 'cis-aconitate', ";
	print OUT "42 'citrate' \/  'no citrate' 'citrate', ";
	print OUT "43 'co2' \/  'no co2' 'co2', ";
	print OUT "44 'cyclodextrin' \/  'no cyclodextrin' 'cyclodextrin', ";
	print OUT "45 'dextrin' \/  'no dextrin' 'dextrin', ";
	print OUT "46 'dulcitol' \/  'no dulcitol' 'dulcitol', ";
	print OUT "47 'erythritol' \/  'no erythritol' 'erythritol', ";
	print OUT "48 'ethanol' \/  'no ethanol' 'ethanol', ";
	print OUT "49 'formate' \/  'no formate' 'formate', ";
	print OUT "50 'fructose 6-phosphate' \/  'no fructose 6-phosphate' 'fructose 6-phosphate', ";
	print OUT "51 'fructose' \/  'no fructose' 'fructose', ";
	print OUT "52 'fucose' \/  'no fucose' 'fucose', ";
	print OUT "53 'fumarate' \/  'no fumarate' 'fumarate', ";
	print OUT "54 'fusidic acid' \/  'no fusidic acid' 'fusidic acid', ";
	print OUT "55 'galactose' \/  'no galactose' 'galactose', ";
	print OUT "56 'galacturonate' \/  'no galacturonate' 'galacturonate', ";
	print OUT "57 'gelatin' \/  'no gelatin' 'gelatin', ";
	print OUT "58 'gentiobiose' \/  'no gentiobiose' 'gentiobiose', ";
	print OUT "59 'gluconate' \/  'no gluconate' 'gluconate', ";
	print OUT "60 'glucose 1-phosphate' \/  'no glucose 1-phosphate' 'glucose 1-phosphate', ";
	print OUT "61 'glucose' \/  'no glucose' 'glucose', ";
	print OUT "62 'glucuronate' \/  'no glucuronate' 'glucuronate', ";
	print OUT "63 'glutamate' \/  'no glutamate' 'glutamate', ";
	print OUT "64 'glutamic acid polymer' \/  'no glutamic acid polymer' 'glutamic acid polymer', ";
	print OUT "65 'glutamine' \/  'no glutamine' 'glutamine', ";
	print OUT "66 'glutamyl glutamic acid' \/  'no glutamyl glutamic acid' 'glutamyl glutamic acid', ";

	print OUT "67 'amino acids' \/  'no amino acids' 'amino acids', ";

	print OUT "68 'glycerol phosphate' \/  'no glycerol phosphate' 'glycerol phosphate', ";
	print OUT "69 'glycerol' \/  'no glycerol' 'glycerol', ";
	print OUT "70 'glycogen' \/  'no glycogen' 'glycogen', ";
	print OUT "71 'glycyl l-glutamic acid' \/  'no glycyl l-glutamic acid' 'glycyl l-glutamic acid', ";
	print OUT "72 'glycyl-l-proline' \/  'no glycyl-l-proline' 'glycyl-l-proline', ";
	print OUT "73 'glyoxylate' \/  'no glyoxylate' 'glyoxylate', ";
	print OUT "74 'guanidine' \/  'no guanidine' 'guanidine', ";
	print OUT "75 'h2s' \/  'no h2s' 'h2s', ";
	print OUT "76 'hippurate' \/  'no hippurate' 'hippurate', ";
	print OUT "77 'histidine' \/  'no histidine' 'histidine', ";
	print OUT "78 'indole' \/  'no indole' 'indole', ";
	print OUT "79 'inosine' \/  'no inosine' 'inosine', ";
	print OUT "80 'inositol' \/  'no inositol' 'inositol', ";
	print OUT "81 'inulin' \/  'no inulin' 'inulin', ";
	print OUT "82 'itaconate' \/  'no itaconate' 'itaconate', ";
	print OUT "83 'lactate' \/  'no lactate' 'lactate', ";
	print OUT "84 'lactic acid methyl ester' \/  'no lactic acid methyl ester' 'lactic acid methyl ester', ";
	print OUT "85 'lactose' \/  'no lactose' 'lactose', ";
	print OUT "86 'lactulose' \/  'no lactulose' 'lactulose', ";
	print OUT "87 'leucine' \/  'no leucine' 'leucine', ";
	print OUT "88 'lyxose' \/  'no lyxose' 'lyxose', ";
	print OUT "89 'malate' \/  'no malate' 'malate', ";
	print OUT "90 'malonate' \/  'no malonate' 'malonate', ";

	print OUT "91 'carboxylic acids' \/  'no carboxylic acids' 'carboxylic acids', ";

	print OUT "92 'maltose' \/  'no maltose' 'maltose', ";
	print OUT "93 'maltotiose' \/  'no maltotiose' 'maltotiose', ";
	print OUT "94 'mannan' \/  'no mannan' 'mannan', ";
	print OUT "95 'mannitol' \/  'no mannitol' 'mannitol', ";
	print OUT "96 'mannose' \/  'no mannose' 'mannose', ";
	print OUT "97 'melezitose' \/  'no melezitose' 'melezitose', ";
	print OUT "98 'melibiose' \/  'no melibiose' 'melibiose', ";
	print OUT "99 'meso-tartrate' \/  'no meso-tartrate' 'meso-tartrate', ";
	print OUT "100 'methionine' \/  'no methionine' 'methionine', ";
	print OUT "101 'methyl a-d-galactoside' \/  'no methyl a-d-galactoside' 'methyl a-d-galactoside', ";
	print OUT "102 'methyl a-d-glucoside' \/  'no methyl a-d-glucoside' 'methyl a-d-glucoside', ";
	print OUT "103 'methyl b-d-galactoside' \/  'no methyl b-d-galactoside' 'methyl b-d-galactoside', ";
	print OUT "104 'methyl b-d-glucoside' \/  'no methyl b-d-glucoside' 'methyl b-d-glucoside', ";
	print OUT "105 'methyl pyruvate' \/  'no methyl pyruvate' 'methyl pyruvate', ";
	print OUT "106 'monomethyl succinate' \/  'no monomethyl succinate' 'monomethyl succinate', ";
	print OUT "107 'mucate' \/  'no mucate' 'mucate', ";
	print OUT "108 'myo-inositol' \/  'no myo-inositol' 'myo-inositol', ";
	print OUT "109 'n-acetyl-b-d-glucosamine' \/  'no n-acetyl-b-d-glucosamine' 'n-acetyl-b-d-glucosamine', ";
	print OUT "110 'n-acetyl-b-d-mannosamine' \/  'no n-acetyl-b-d-mannosamine' 'n-acetyl-b-d-mannosamine', ";
	print OUT "111 'n-acetyl-l-glutamate' \/  'no n-acetyl-l-glutamate' 'n-acetyl-l-glutamate', ";
	print OUT "112 'ornithine' \/  'no ornithine' 'ornithine', ";
	print OUT "113 'palatinose' \/  'no palatinose' 'palatinose', ";
	print OUT "114 'pectin' \/  'no pectin' 'pectin', ";
	print OUT "115 'phenylacetate' \/  'no phenylacetate' 'phenylacetate', ";
	print OUT "116 'phenylalanine' \/  'no phenylalanine' 'phenylalanine', ";
	print OUT "117 'phenylethylamine' \/  'no phenylethylamine' 'phenylethylamine', ";
	print OUT "118 'proline' \/  'no proline' 'proline', ";
	print OUT "119 'propionate' \/  'no propionate' 'propionate', ";
	print OUT "120 'psicose' \/  'no psicose' 'psicose', ";
	print OUT "121 'pullulan' \/  'no pullulan' 'pullulan', ";
	print OUT "122 'putrescine' \/  'no putrescine' 'putrescine', ";
	print OUT "123 'pyroglutamate' \/  'no pyroglutamate' 'pyroglutamate', ";
	print OUT "124 'pyruvate' \/  'no pyruvate' 'pyruvate', ";

	print OUT "125 'nucleobases' \/  'no nucleobases' 'nucleobases', ";

	print OUT "126 'quinate' \/  'no quinate' 'quinate', ";
	print OUT "127 'raffinose' \/  'no raffinose' 'raffinose', ";
	print OUT "128 'rhamnose' \/  'no rhamnose' 'rhamnose', ";
	print OUT "129 'ribitol' \/  'no ribitol' 'ribitol', ";
	print OUT "130 'ribose' \/  'no ribose' 'ribose', ";
	print OUT "131 'saccharolytic' \/  'no saccharolytic' 'saccharolytic', ";

	print OUT "132 'saccharolytic enzymes' \/  'no saccharolytic enzymes' 'saccharolytic enzymes', ";

	print OUT "133 'salicin' \/  'no salicin' 'salicin', ";
	print OUT "134 'sedoheptulosan' \/  'no sedoheptulosan' 'sedoheptulosan', ";
	print OUT "135 'serine' \/  'no serine' 'serine', ";
	print OUT "136 'sorbitol' \/  'no sorbitol' 'sorbitol', ";
	print OUT "137 'sorbose' \/  'no sorbose' 'sorbose', ";
	print OUT "138 'stachyose' \/  'no stachyose' 'stachyose', ";
	print OUT "139 'starch' \/  'no starch' 'starch', ";
	print OUT "140 'suberate' \/  'no suberate' 'suberate', ";
	print OUT "141 'succinamic acid' \/  'no succinamic acid' 'succinamic acid', ";
	print OUT "142 'succinate' \/  'no succinate' 'succinate', ";

	print OUT "143 'cysteine' \/  'no cysteine' 'cysteine', ";

	print OUT "144 'sucrose' \/  'no sucrose' 'sucrose', ";
	print OUT "145 'tagatose' \/  'no tagatose' 'tagatose', ";
	print OUT "146 'tartrate' \/  'no tartrate' 'tartrate', ";
	print OUT "147 'tellurite' \/  'no tellurite' 'tellurite', ";
	print OUT "148 'threonine' \/  'no threonine' 'threonine', ";
	print OUT "149 'thymidine 5prime-monophosphate' \/  'no thymidine 5prime-monophosphate' 'thymidine 5prime-monophosphate', ";
	print OUT "150 'thymidine' \/  'no thymidine' 'thymidine', ";
	print OUT "151 'trehalose' \/  'no trehalose' 'trehalose', ";

	print OUT "152 'glycine' \/  'no glycine' 'glycine', ";

	print OUT "153 'turanose' \/  'no turanose' 'turanose', ";
	print OUT "154 'tween 40' \/  'no tween 40' 'tween 40', ";
	print OUT "155 'tween 80' \/  'no tween 80' 'tween 80', ";
	print OUT "156 'uridine' \/  'no uridine' 'uridine', ";
	print OUT "157 'valerate' \/  'no valerate' 'valerate', ";
	print OUT "158 'valine' \/  'no valine' 'valine', ";
	print OUT "159 'xylan' \/  'no xylan' 'xylan', ";
	print OUT "160 'xylitol' \/  'no xylitol' 'xylitol', ";
	print OUT "161 'xylose' \/  'no xylose' 'xylose', ";

	print OUT "162 'isoleucine' \/  'no isoleucine' 'isoleucine', ";

	print OUT "163 'lysine' \/  'no lysine' 'lysine', ";
	print OUT "164 'tryptophan' \/  'no tryptophan' 'tryptophan', ";
	print OUT "165 'selenocysteine' \/  'no selenocysteine' 'selenocysteine', ";
	print OUT "166 'pyrrolysine' \/  'no pyrrolysine' 'pyrrolysine', ";

	print OUT "167 'nitrate reduction' \/  'no nitrate reduction' 'nitrate reduction', ";
	print OUT "168 'nitrogen fixation' \/  'no nitrogen fixation' 'nitrogen fixation', ";
	print OUT "169 'nitrite reduction' \/  'no nitrite reduction' 'nitrite reduction', ";

	print OUT "170 '6-phospho-b-galactosidase' \/  'no 6-phospho-b-galactosidase' '6-phospho-b-galactosidase', ";
	print OUT "171 '6-phosphogluconate dehydrogenase' \/  'no 6-phosphogluconate dehydrogenase' '6-phosphogluconate dehydrogenase', ";
	print OUT "172 'a-arabinosidase' \/  'no a-arabinosidase' 'a-arabinosidase', ";
	print OUT "173 'a-fucosidase' \/  'no a-fucosidase' 'a-fucosidase', ";
	print OUT "174 'a-galactosidase' \/  'no a-galactosidase' 'a-galactosidase', ";
	print OUT "175 'a-glucosidase' \/  'no a-glucosidase' 'a-glucosidase', ";
	print OUT "176 'a-mannosidase' \/  'no a-mannosidase' 'a-mannosidase', ";
	print OUT "177 'a-methyl glucosidase' \/  'no a-methyl glucosidase' 'a-methyl glucosidase', ";
	print OUT "178 'a-methylmannoside' \/  'no a-methylmannoside' 'a-methylmannoside', ";
	print OUT "179 'acid phosphatase' \/  'no acid phosphatase' 'acid phosphatase', ";
	print OUT "180 'aesculin hydrolysis' \/  'no aesculin hydrolysis' 'aesculin hydrolysis', ";
	print OUT "181 'alanine arylamidase' \/  'no alanine arylamidase' 'alanine arylamidase', ";
	print OUT "182 'alkaline phosphatase' \/  'no alkaline phosphatase' 'alkaline phosphatase', ";
	print OUT "183 'alpha-glutamate arylamidase' \/  'no alpha-glutamate arylamidase' 'alpha-glutamate arylamidase', ";
	print OUT "184 'amylase' \/  'no amylase' 'amylase', ";
	print OUT "185 'arginine arylamidase' \/  'no arginine arylamidase' 'arginine arylamidase', ";
	print OUT "186 'arginine decarboxylase' \/  'no arginine decarboxylase' 'arginine decarboxylase', ";
	print OUT "187 'arginine dihydrolase' \/  'no arginine dihydrolase' 'arginine dihydrolase', ";
	print OUT "188 'arylsulfatase' \/  'no arylsulfatase' 'arylsulfatase', ";
	print OUT "189 'asparagine arylamidase' \/  'no asparagine arylamidase' 'asparagine arylamidase', ";
	print OUT "190 'aspartate arylamidase' \/  'no aspartate arylamidase' 'aspartate arylamidase', ";
	print OUT "191 'b-galactosidase' \/  'no b-galactosidase' 'b-galactosidase', ";
	print OUT "192 'b-glucosidase' \/  'no b-glucosidase' 'b-glucosidase', ";
	print OUT "193 'b-glucuronidase' \/  'no b-glucuronidase' 'b-glucuronidase', ";
	print OUT "194 'b-mannosidase' \/  'no b-mannosidase' 'b-mannosidase', ";
	print OUT "195 'casein hydrolysis' \/  'no casein hydrolysis' 'casein hydrolysis', ";
	print OUT "196 'caseinase' \/  'no caseinase' 'caseinase', ";
	print OUT "197 'catalase' \/  'no catalase' 'catalase', ";
	print OUT "198 'cellobiose' \/  'no cellobiose' 'cellobiose', ";
	print OUT "199 'chymotrypsin' \/  'no chymotrypsin' 'chymotrypsin', ";
	print OUT "200 'coagulase' \/  'no coagulase' 'coagulase', ";
	print OUT "201 'collagenase' \/  'no collagenase' 'collagenase', ";
	print OUT "202 'cysteine arylamidase' \/  'no cysteine arylamidase' 'cysteine arylamidase', ";
	print OUT "203 'cystine arylamidase' \/  'no cystine arylamidase' 'cystine arylamidase', ";
	print OUT "204 'cytochrome oxidase' \/  'no cytochrome oxidase' 'cytochrome oxidase', ";
	print OUT "205 'dnase' \/  'no dnase' 'dnase', ";
	print OUT "206 'esterase c4' \/  'no esterase c4' 'esterase c4', ";
	print OUT "207 'fructose bisphosphate aldolase' \/  'no fructose bisphosphate aldolase' 'fructose bisphosphate aldolase', ";
	print OUT "208 'fucosidase' \/  'no fucosidase' 'fucosidase', ";
	print OUT "209 'gelatinase' \/  'no gelatinase' 'gelatinase', ";
	print OUT "210 'glucosidase' \/  'no glucosidase' 'glucosidase', ";
	print OUT "211 'glutamate arylamidase' \/  'no glutamate arylamidase' 'glutamate arylamidase', ";
	print OUT "212 'glutamate deaminase' \/  'no glutamate deaminase' 'glutamate deaminase', ";
	print OUT "213 'glutamic acid decarboxylase' \/  'no glutamic acid decarboxylase' 'glutamic acid decarboxylase', ";
	print OUT "214 'glutamine arylamidase' \/  'no glutamine arylamidase' 'glutamine arylamidase', ";
	print OUT "215 'glutamyl glutamic acid arylamidase' \/  'no glutamyl glutamic acid arylamidase' 'glutamyl glutamic acid arylamidase', ";
	print OUT "216 'glycine arylamidase' \/  'no glycine arylamidase' 'glycine arylamidase', ";
	print OUT "217 'hippurate hydrolysis' \/  'no hippurate hydrolysis' 'hippurate hydrolysis', ";
	print OUT "218 'histidine arylamidase' \/  'no histidine arylamidase' 'histidine arylamidase', ";
	print OUT "219 'hyaluronidase' \/  'no hyaluronidase' 'hyaluronidase', ";
	print OUT "220 'hydroxy proline arylamidase' \/  'no hydroxy proline arylamidase' 'hydroxy proline arylamidase', ";
	print OUT "221 'iron-phorphyrins' \/  'no iron-phorphyrins' 'iron-phorphyrins', ";
	print OUT "222 'isoleucine arylamidase' \/  'no isoleucine arylamidase' 'isoleucine arylamidase', ";
	print OUT "223 'l-lactic acid racemase' \/  'no l-lactic acid racemase' 'l-lactic acid racemase', ";
	print OUT "224 'lecithinase' \/  'no lecithinase' 'lecithinase', ";
	print OUT "225 'leucine arylamidase' \/  'no leucine arylamidase' 'leucine arylamidase', ";
	print OUT "226 'leucyl glycine arylamidase' \/  'no leucyl glycine arylamidase' 'leucyl glycine arylamidase', ";
	print OUT "227 'lipase c14' \/  'no lipase c14' 'lipase c14', ";
	print OUT "228 'lipase c8' \/  'no lipase c8' 'lipase c8', ";
	print OUT "229 'lipase' \/  'no lipase' 'lipase', ";
	print OUT "230 'lysine arylamidase' \/  'no lysine arylamidase' 'lysine arylamidase', ";
	print OUT "231 'lysine deaminase' \/  'no lysine deaminase' 'lysine deaminase', ";
	print OUT "232 'lysine decarboxylase' \/  'no lysine decarboxylase' 'lysine decarboxylase', ";
	print OUT "233 'methionine arylamidase' \/  'no methionine arylamidase' 'methionine arylamidase', ";
	print OUT "234 'methyl a-d-mannoside' \/  'no methyl a-d-mannoside' 'methyl a-d-mannoside', ";
	print OUT "235 'methyl b-d-xyloside' \/  'no methyl b-d-xyloside' 'methyl b-d-xyloside', ";
	print OUT "236 'n-acetyl-b-d-glucosaminidase' \/  'no n-acetyl-b-d-glucosaminidase' 'n-acetyl-b-d-glucosaminidase', ";
	print OUT "237 'nadh oxidase' \/  'no nadh oxidase' 'nadh oxidase', ";
	print OUT "238 'naphthol-as-bi-phosphohydrolase' \/  'no naphthol-as-bi-phosphohydrolase' 'naphthol-as-bi-phosphohydrolase', ";
	print OUT "239 'ornithine arylamidase' \/  'no ornithine arylamidase' 'ornithine arylamidase', ";
	print OUT "240 'ornithine decarboxylase' \/  'no ornithine decarboxylase' 'ornithine decarboxylase', ";
	print OUT "241 'oxidase' \/  'no oxidase' 'oxidase', ";
	print OUT "242 'para-phenylalanine deaminase' \/  'no para-phenylalanine deaminase' 'para-phenylalanine deaminase', ";
	print OUT "243 'phenylalanine arylamidase' \/  'no phenylalanine arylamidase' 'phenylalanine arylamidase', ";
	print OUT "244 'phenylalanine deaminase' \/  'no phenylalanine deaminase' 'phenylalanine deaminase', ";
	print OUT "245 'phosphatase' \/  'no phosphatase' 'phosphatase', ";
	print OUT "246 'proline arylamidase' \/  'no proline arylamidase' 'proline arylamidase', ";
	print OUT "247 'protease' \/  'no protease' 'protease', ";
	print OUT "248 'pyrazinamidase' \/  'no pyrazinamidase' 'pyrazinamidase', ";
	print OUT "249 'pyroglutamic acid arylamidase' \/  'no pyroglutamic acid arylamidase' 'pyroglutamic acid arylamidase', ";
	print OUT "250 'proteolytic enzymes' \/  'no proteolytic enzymes' 'proteolytic enzymes', ";
	print OUT "251 'serine arylamidase' \/  'no serine arylamidase' 'serine arylamidase', ";
	print OUT "252 'threonine arylamidase' \/  'no threonine arylamidase' 'threonine arylamidase', ";
	print OUT "253 'trypsin' \/  'no trypsin' 'trypsin', ";
	print OUT "254 'tryptophan arylamidase' \/  'no tryptophan arylamidase' 'tryptophan arylamidase', ";
	print OUT "255 'tryptophan deaminase' \/  'no tryptophan deaminase' 'tryptophan deaminase', ";
	print OUT "256 'tryptophanase' \/  'no tryptophanase' 'tryptophanase', ";
	print OUT "257 'tyrosine arylamidase' \/  'no tyrosine arylamidase' 'tyrosine arylamidase', ";
	print OUT "258 'urease' \/  'no urease' 'urease', ";
	print OUT "259 'valine arylamidase' \/  'no valine arylamidase' 'valine arylamidase', ";
	print OUT "260 'xylanase' \/  'no xylanase' 'xylanase', ";

	print OUT "261 'a-hemolysis' \/  'no a-hemolysis' 'a-hemolysis', ";
	print OUT "262 'acetoin' \/  'no acetoin' 'acetoin', ";
	print OUT "263 'b-haemolysis' \/  'no b-haemolysis' 'b-haemolysis', ";
	print OUT "264 'ellman' \/  'no ellman' 'ellman', ";
	print OUT "265 'flexirubin' \/  'no flexirubin' 'flexirubin', ";
	print OUT "266 'haemolysis' \/  'no haemolysis' 'haemolysis', ";
	print OUT "267 'indoxyl acetate hydrolysis' \/  'no indoxyl acetate hydrolysis' 'indoxyl acetate hydrolysis', ";
	print OUT "268 'methyl red acidification' \/  'no methyl red acidification' 'methyl red acidification', ";
	print OUT "269 'methylene blue reduction' \/  'no methylene blue reduction' 'methylene blue reduction', ";
	print OUT "270 'resazurin reduction' \/  'no resazurin reduction' 'resazurin reduction', ";
	print OUT "271 'tetrathionate reduction' \/  'no tetrathionate reduction' 'tetrathionate reduction', ";
	print OUT "272 'tetrazolium reduction' \/  'no tetrazolium reduction' 'tetrazolium reduction', ";
	print OUT "273 'carboxymethyl cellulase' \/  'no carboxymethyl cellulase' 'carboxymethyl cellulase', ";

	print OUT "274 '1,3-propanediol' \/  'no 1,3-propanediol' '1,3-propanediol', ";
	print OUT "275 '1,2-propanediol' \/  'no 1,2-propanediol' '1,2-propanediol', ";
	print OUT "276 'propanediol' \/  'no propanediol' 'propanediol', ";
	print OUT "277 'phosphatidylinositol phospholipase c' \/  'no phosphatidylinositol phospholipase c' 'phosphatidylinositol phospholipase c', ";
	print OUT "278 'lecithinase' \/  'no lecithinase' 'lecithinase', ";
	print OUT "279 'lactamide' \/  'no lactamide' 'lactamide', ";
	print OUT "280 'hydroxyproline' \/  'no hydroxyproline' 'hydroxyproline', ";
	print OUT "281 'glycyl tryptophan arylamidase' \/  'no glycyl tryptophan arylamidase' 'glycyl tryptophan arylamidase', ";
	print OUT "282 'lipolytic enzymes' \/  'no lipolytic enzymes' 'lipolytic enzymes', ";
	print OUT "283 'esterase' \/  'no esterase' 'esterase', ";
#	print OUT "284 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "285 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "286 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "287 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "288 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "289 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "290 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "291 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "292 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "293 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "294 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "295 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "296 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "297 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "298 'xxx' \/  'no xxx' 'xxx', ";
#	print OUT "299 'xxx' \/  'no xxx' 'xxx', ";

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Taxon.*/) {
			next;
			}
		my @testsnegdata = split (/\t/, $line);
		push (my @taxnames, $testsnegdata[0]);

#code char 1 2-aminoethanol
		if ($line =~ /,2-aminoethanol,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not 2-aminoethanol,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
			}
#code char 2 2-keto-d-gluconate
		if ($line =~ /,2-keto-d-gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-keto-d-gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 3 2-ketobutyrate
		if ($line =~ /,2-ketobutyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-ketobutyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 4 2,3-butanediol
		if ($line =~ /,2-3-butanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 2-3-butanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 5 1,3-butanediol
		if ($line =~ /,1-3-butanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-3-butanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 6 3-hydroxybenzoate
		if ($line =~ /,3-hydroxybenzoate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-hydroxybenzoate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 7 3-hydroxybutyrate
		if ($line =~ /,3-hydroxybutyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-hydroxybutyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 8 3-methyl d-glucose
		if ($line =~ /,3-methyl d-glucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 3-methyl d-glucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 9 4-hydroxybenzoate
		if ($line =~ /,4-hydroxybenzoate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 4-hydroxybenzoate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 10 4-hydroxyphenylacetate
		if ($line =~ /,4-hydroxyphenylacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 4-hydroxyphenylacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 11 alcohols
		if ($line =~ /,1-3-propanediol,|,1-2-propanediol,|,propanediol,|,lactamide,|,alcohols,|,2-aminoethanol,|,2-3-butanediol,|,1-3-butanediol,|,arabitol,|,dulcitol,|,ethythritol,|,ethanol,|,glycerol phosphate,|,glycerol,|,inositol,|,mannitol,|,myo-inositol,|,ribitol,|,salicin,|,sorbitol,|,xylitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alcohols,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 12 5-keto-d-gluconate
		if ($line =~ /,5-keto-d-gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 5-keto-d-gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 13 a-ketoglutarate
		if ($line =~ /,a-ketoglutarate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-ketoglutarate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 14 a-ketovalerate
		if ($line =~ /,a-ketovalerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-ketovalerate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 15 acetate
		if ($line =~ /,acetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 16 acetoacetate
		if ($line =~ /,acetoacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetoacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 17 adenosine 5'-monophosphate
		if ($line =~ /,adenosine 5\'-monophosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adenosine 5\'-monophosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 18 adenosine
		if ($line =~ /,adenosine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adenosine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 19 adipate
		if ($line =~ /,adipate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not adipate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 20 organonitrogen compounds
		if ($line =~ /,organonitrogen compounds,|,2-aminoethanol,|,alaninamide,|,alanine,|,alanyl glycine,|,amygdalin,|,asparagine,|,aspartate,|,carnitine,|,glutamic acid polymer,|,glutamine,|,glutamyl glutamic acid,|,glycyl l-glutamic acid,|,glycyl-l-proline,|,guanidine,|,hippurate,|,histidine,|,leucine,|,methionine,|,n-acetyl-b-d-glucosamine,|,n-acetyl-b-d-mannosamine,|,n-acetyl-l-glutamate,|,ornithine,|,phenylalanine,|,phenylethylamine,|,proline,|,putrescine,|,pyroglutamate,|,serine,|,threonine,|,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not organonitrogen compounds,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 21 agar
		if ($line =~ /,agar,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not agar,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 22 alaninamide
		if ($line =~ /,alaninamide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alaninamide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 23 alanine
		if ($line =~ /,alanine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 24 alanyl glycine
		if ($line =~ /,alanyl glycine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanyl glycine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 25 amygdalin
		if ($line =~ /,amygdalin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amygdalin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 26 arabinose
		if ($line =~ /,arabinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arabinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 27 arabitol
		if ($line =~ /,arabitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arabitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 28 arbutin
		if ($line =~ /,arbutin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arbutin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 29 arginine
		if ($line =~ /,arginine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 30 asparagine
		if ($line =~ /,asparagine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not asparagine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 31 aspartate
		if ($line =~ /,aspartate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aspartate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 32 b-cyclodextrin
		if ($line =~ /,b-cyclodextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-cyclodextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 33 fatty acids
		if ($line =~ /,fatty acids,|,2-ketyobutyrate,|,3-hydroxybutyrate,|,a-ketovalerate,|,butyrate,|,propionate,|,valerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fatty acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 34 b-xylosidase
		if ($line =~ /,b-xylosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-xylosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 35 butyrate
		if ($line =~ /,butyrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not butyrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 36 caprate
		if ($line =~ /,caprate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not caprate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 37 carbohydrates
		if ($line =~ /,carbohydrates,|,3-methyl d-glucose,|,agar,|,amygdalin,|,arabinose,|,arabitol,|,arbutin,|,b-cyclodextrin,|,carboxymethyl cellulose,|,cellobiose,|,cyclodextrin,|,dextrin,|,dulcitol,|,erythritol,|,fructose 6-phosphate,|,fructose,|,fucose,|,galactose,|,galacturonate,|,gelatin,|,gentiobiose,|,gluconate,|,glucose 1-phosphate,|,glucose,|,glucuronate,|,glycerol phosphate,|,glycerol,|,glycogen,|,inulin,|,lactose,|,lactulose,|,lyxose,|,maltose,|,maltotriose,|,mannan,|,mannitol,|,mannose,|,melezitose,|,melibiose,|,meso-tartrate,|,methyl a-d-galactoside,|,methyl a-d-glucoside,|,methyl b-d-glactoside,|,methyl b-d-glucoside,|,mucate,|,n-acetyl-b-d-glucosamine,|,n-acetyl-b-d-mannosamine,|,palatinose,|,pectin,|,psicose,|,pullulan,|,raffinose,|,rhamnose,|,ribitol,|,ribose,|,salicin,|,sedoheptulosan,|,sorbitol,|,sorbose,|,stachyose,|,starch,|,sucrose,|,tagatose,|,tartrate,|,trehalose,|,turanose,|,xylan,|,xylitol,|,xylose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carbohydrates,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 38 carboxymethyl cellulose
		if ($line =~ /,carboxymethyl cellulose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxymethyl cellulose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 39 carnitine
		if ($line =~ /,carnitine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carnitine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 40 cellobiose
		if ($line =~ /,cellobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cellobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 41 cis-aconitate
		if ($line =~ /,cis-aconitate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cis-aconitate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 42 citrate
		if ($line =~ /,citrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not citrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 43 co2
		if ($line =~ /,co2,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not co2,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 44 cyclodextrin
		if ($line =~ /,cyclodextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cyclodextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 45 dextrin
		if ($line =~ /,dextrin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dextrin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 46 dulcitol
		if ($line =~ /,dulcitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dulcitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 47 erythritol
		if ($line =~ /,erythritol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not erythritol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 48 ethanol
		if ($line =~ /,ethanol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ethanol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 49 formate
		if ($line =~ /,formate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not formate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 50 fructose 6-phosphate
		if ($line =~ /,fructose 6-phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose 6-phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 51 fructose
		if ($line =~ /,fructose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 52 fucose
		if ($line =~ /,fucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 53 fumarate
		if ($line =~ /,fumarate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fumarate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 54 fusidic acid
		if ($line =~ /,fusidic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fusidic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 55 galactose
		if ($line =~ /,galactose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not galactose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 56 galacturonate
		if ($line =~ /,galacturonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not galacturonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 57 gelatin
		if ($line =~ /,gelatin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gelatin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 58 gentiobiose
		if ($line =~ /,gentiobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gentiobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 59 gluconate
		if ($line =~ /,gluconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gluconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 60 glucose 1-phosphate
		if ($line =~ /,glucose 1-phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucose 1-phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 61 glucose
		if ($line =~ /,glucose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 62 glucuronate
		if ($line =~ /,glucuronate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucuronate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 63 glutamate
		if ($line =~ /,glutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 64 glutamic acid polymer
		if ($line =~ /,glutamic acid polymer,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamic acid polymer,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 65 glutamine
		if ($line =~ /,glutamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 66 glutamyl glutamic acid
		if ($line =~ /,glutamyl glutamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamyl glutamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 67 amino acids
		if ($line =~ /,hydroxyproline,|,amino acids,|,alaninamide,|,alanine,|,alanyl glycine,|,arginine,|,asparagine,|,aspartate,|,glutamate,|,glutamic acid polymer,|,glutamine,|,glutamyl glutamic acid,|,glycine,|,glycyl l-glutamic acid,|,glycyl-l-proline,|,histidine,|,isoleucine,|,leucine,|,lysine,|,methionine,|,ornithine,|,phenylalanine,|,proline,|,pyroglutamate,|,pyrrolysine,|,selenocysteine,|,serine,|,threonine,|,tryptophan,|,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amino acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 68 glycerol phosphate
		if ($line =~ /,glycerol phosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycerol phosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 69 glycerol
		if ($line =~ /,glycerol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycerol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 70 glycogen
		if ($line =~ /,glycogen,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycogen,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 71 glycyl l-glutamic acid
		if ($line =~ /,glycyl l-glutamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl l-glutamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 72 glycyl-l-proline
		if ($line =~ /,glycyl-l-proline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl-l-proline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 73 glyoxylate
		if ($line =~ /,glyoxylate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glyoxylate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 74 guanidine
		if ($line =~ /,guanidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not guanidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 75 h2s
		if ($line =~ /,h2s,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not h2s,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 76 hippurate
		if ($line =~ /,hippurate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hippurate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 77 histidine
		if ($line =~ /,histidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not histidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 78 indole
		if ($line =~ /,indole,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not indole,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 79 inosine
		if ($line =~ /,inosine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inosine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 80 inositol
		if ($line =~ /,inositol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inositol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 81 inulin
		if ($line =~ /,inulin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not inulin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 82 itaconate
		if ($line =~ /,itaconate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not itaconate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 83 lactate
		if ($line =~ /,lactate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 84 lactic acid methyl ester
		if ($line =~ /,lactic acid methyl ester,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactic acid methyl ester,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 85 lactose
		if ($line =~ /,lactose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 86 lactulose
		if ($line =~ /,lactulose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactulose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 87 leucine
		if ($line =~ /,leucine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 88 lyxose
		if ($line =~ /,lyxose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lyxose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 89 malate
		if ($line =~ /,malate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not malate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 90 malonate
		if ($line =~ /,malonate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not malonate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 91 carboxylic acids
		if ($line =~ /,hydroxyproline,|,carboxylic acids,|,2-keto-d-gluconate,|,2-ketobutyrate,|,3-hydroxybenzoate,|,3-hydroxybutyrate,|,4-hydroxybenzoate,|,4-hydroxyphenylacetate,|,5-keto-d-gluconate,|,a-ketoglutarate,|,a-ketovalerate,|,acetate,|,acetoacetate,|,adipate,|,butyrate,|,cis-aconitate,|,citrate,|,formate,|,fumarate,|,fusidic acid,|,galacturonate,|,gluconate,|,glucuronate,|,glutamate,|,glyoxylate,|,itaconate,|,lactate,|,malate,|,malonate,|,meso-tartrate,|,methyl pyruvate,|,monomethyl succinate,|,mucate,|,phenylacetate,|,pyruvate,|,quinate,|,suberate,|,succinamic acid,|,tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxylic acids,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 92 maltose
		if ($line =~ /,maltose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not maltose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 93 maltotiose
		if ($line =~ /,maltotiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not maltotiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 94 mannan
		if ($line =~ /,mannan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 95 mannitol
		if ($line =~ /,mannitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 96 mannose
		if ($line =~ /,mannose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mannose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 97 melezitose
		if ($line =~ /,melezitose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not melezitose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 98 melibiose
		if ($line =~ /,melibiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not melibiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 99 meso-tartrate
		if ($line =~ /,meso-tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not meso-tartrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 100 methionine
		if ($line =~ /,methionine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methionine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 101 methyl a-d-galactoside
		if ($line =~ /,methyl a-d-galactoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-galactoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 102 methyl a-d-glucoside
		if ($line =~ /,methyl a-d-glucoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-glucoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 103 methyl b-d-galactoside
		if ($line =~ /,methyl b-d-galactoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-galactoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 104 methyl b-d-glucoside
		if ($line =~ /,methyl b-d-glucoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-glucoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 105 methyl pyruvate
		if ($line =~ /,methyl pyruvate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl pyruvate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 106 monomethyl succinate
		if ($line =~ /,monomethyl succinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not monomethyl succinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 107 mucate
		if ($line =~ /,mucate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not mucate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 108 myo-inositol
		if ($line =~ /,myo-inositol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not myo-inositol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 109 n-acetyl-b-d-glucosamine
		if ($line =~ /,n-acetyl-b-d-glucosamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-glucosamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 110 n-acetyl-b-d-mannosamine
		if ($line =~ /,n-acetyl-b-d-mannosamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-mannosamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 111 n-acetyl-l-glutamate
		if ($line =~ /,n-acetyl-l-glutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-l-glutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 112 ornithine
		if ($line =~ /,ornithine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 113 palatinose
		if ($line =~ /,palatinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not palatinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 114 pectin
		if ($line =~ /,pectin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pectin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 115 phenylacetate
		if ($line =~ /,phenylacetate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylacetate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 116 phenylalanine
		if ($line =~ /,phenylalanine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 117 phenylethylamine
		if ($line =~ /,phenylethylamine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylethylamine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 118 proline
		if ($line =~ /,proline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not proline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 119 propionate
		if ($line =~ /,propionate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propionate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 120 psicose
		if ($line =~ /,psicose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not psicose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 121 pullulan
		if ($line =~ /,pullulan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pullulan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 122 putrescine
		if ($line =~ /,putrescine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not putrescine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 123 pyroglutamate
		if ($line =~ /,pyroglutamate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyroglutamate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 124 pyruvate
		if ($line =~ /,pyruvate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyruvate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 125 nucleobases
		if ($line =~ /,nucleobases,|,adenosine 5\'-monophosphate,|,adenosine,|,thymidine 5\'-monophosphate,|,thymidine,|,uridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nucleobases,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 126 quinate
		if ($line =~ /,quinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not quinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 127 raffinose
		if ($line =~ /,raffinose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not raffinose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 128 rhamnose
		if ($line =~ /,rhamnose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not rhamnose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 129 ribitol
		if ($line =~ /,ribitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ribitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 130 ribose
		if ($line =~ /,ribose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ribose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 131 saccharolytic
		if ($line =~ /,saccharolytic,|,agar,|,amygdaline,|,b-cyclodextrin,|,carboxymethyl cellulose,|,cellobiose,|,cyclodextrin,|,dextrin,|,gelatin,|,gentiobiose,|,glycogen,|,inulin,|,lactose,|,lactulose,|,maltose,|,maltotriose,|,mannan,|,melezitose,|,melibiose,|,palatinose,|,pectin,|,pullulan,|,raffinose,|,stachyose,|,starch,|,trehalose,|,turanose,|,xylan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not saccharolytic,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 132 saccharolytic enzymes
		if ($line =~ /,6-phospho-b-galactosidase,|,a-arabinosidase,|,a-fucosidase,|,a-galactosidase,|,a-glucosidase,|,a-mannosidase,|,a-methyl glucosidase,|,amylase,|,b-galactosidase,|,b-glucosidase,|,b-glucuronidase,|,b-mannosidase,|,fucosidase,|,glucosidase,|,hyaluronidase,|,n-acetyl-b-d-glucosaminidase,|,xylanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 133 salicin
		if ($line =~ /,salicin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not salicin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 134 sedoheptulosan
		if ($line =~ /,sedoheptulosan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sedoheptulosan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 135 serine
		if ($line =~ /,serine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not serine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 136 sorbitol
		if ($line =~ /,sorbitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sorbitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 137 sorbose
		if ($line =~ /,sorbose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sorbose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 138 stachyose
		if ($line =~ /,stachyose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not stachyose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 139 starch
		if ($line =~ /,starch,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not starch,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 140 suberate
		if ($line =~ /,suberate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not suberate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 141 succinamic acid
		if ($line =~ /,succinamic acid,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not succinamic acid,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 142 succinate
		if ($line =~ /,succinate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not succinate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 143 cysteine
		if ($line =~ /,cysteine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cysteine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 144 sucrose
		if ($line =~ /,sucrose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not sucrose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 145 tagatose
		if ($line =~ /,tagatose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tagatose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 146 tartrate
		if ($line =~ /,tartrate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tartrate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 147 tellurite
		if ($line =~ /,tellurite,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tellurite,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 148 threonine
		if ($line =~ /,threonine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not threonine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 149 thymidine 5'-monophosphate
		if ($line =~ /,thymidine 5\'-monophosphate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not thymidine 5\'-monophosphate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 150 thymidine
		if ($line =~ /,thymidine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not thymidine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 151 trehalose
		if ($line =~ /,trehalose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trehalose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 152 glycine
		if ($line =~ /,glycine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 153 turanose
		if ($line =~ /,turanose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not turanose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 154 tween 40
		if ($line =~ /,tween 40,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tween 40,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 155 tween 80
		if ($line =~ /,tween 80,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tween 80,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 156 uridine
		if ($line =~ /,uridine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not uridine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 157 valerate
		if ($line =~ /,valerate,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valerate,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 158 valine
		if ($line =~ /,valine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 159 xylan
		if ($line =~ /,xylan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 160 xylitol
		if ($line =~ /,xylitol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylitol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 161 xylose
		if ($line =~ /,xylose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 162 isoleucine
		if ($line =~ /,isoleucine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isoleucine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 163 lysine
		if ($line =~ /,lysine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 164 tryptophan
		if ($line =~ /,tryptophan,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 165 selenocysteine
		if ($line =~ /,selenocysteine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not selenocysteine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 166 pyrrolysine
		if ($line =~ /,pyrrolysine,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrrolysine,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 167 nitrate reduction
		if ($line =~ /,nitrate reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrate reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 168 nitrogen fixation
		if ($line =~ /,nitrogen fixation,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrogen fixation,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 169 nitrite reduction
		if ($line =~ /,nitrite reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nitrite reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 170 6-phospho-b-galactosidase
		if ($line =~ /,6-phospho-b-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 6-phospho-b-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 171 6-phosphogluconate dehydrogenase
		if ($line =~ /,6-phosphogluconate dehydrogenase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 6-phosphogluconate dehydrogenase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 172 a-arabinosidase
		if ($line =~ /,a-arabinosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-arabinosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 173 a-fucosidase
		if ($line =~ /,a-fucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-fucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 174 a-galactosidase
		if ($line =~ /,a-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 175 a-glucosidase
		if ($line =~ /,a-glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 176 a-mannosidase
		if ($line =~ /,a-mannosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-mannosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 177 a-methyl glucosidase
		if ($line =~ /,a-methyl glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-methyl glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 178 a-methylmannoside
		if ($line =~ /,a-methylmannoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-methylmannoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 179 acid phosphatase
		if ($line =~ /,acid phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acid phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 180 aesculin hydrolysis
		if ($line =~ /,aesculin hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aesculin hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 181 alanine arylamidase
		if ($line =~ /,alanine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alanine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 182 alkaline phosphatase
		if ($line =~ /,alkaline phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alkaline phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 183 alpha-glutamate arylamidase
		if ($line =~ /,alpha-glutamate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not alpha-glutamate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 184 amylase
		if ($line =~ /,amylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 185 arginine arylamidase
		if ($line =~ /,arginine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 186 arginine decarboxylase
		if ($line =~ /,arginine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 187 arginine dihydrolase
		if ($line =~ /,arginine dihydrolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arginine dihydrolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 188 arylsulfatase
		if ($line =~ /,arylsulfatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not arylsulfatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 189 asparagine arylamidase
		if ($line =~ /,asparagine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not asparagine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 190 aspartate arylamidase
		if ($line =~ /,aspartate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not aspartate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 191 b-galactosidase
		if ($line =~ /,b-galactosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-galactosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 192 b-glucosidase
		if ($line =~ /,b-glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 193 b-glucuronidase
		if ($line =~ /,b-glucuronidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-glucuronidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 194 b-mannosidase
		if ($line =~ /,b-mannosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-mannosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 195 casein hydrolysis
		if ($line =~ /,casein hydrolysis,|,caseinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not casein hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 196 caseinase
		if ($line =~ /,caseinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not caseinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 197 catalase
		if ($line =~ /,catalase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not catalase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 198 cellobiose
		if ($line =~ /,cellobiose,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cellobiose,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 199 chymotrypsin
		if ($line =~ /,chymotrypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not chymotrypsin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 200 coagulase
		if ($line =~ /,coagulase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not coagulase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 201 collagenase
		if ($line =~ /,collagenase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not collagenase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 202 cysteine arylamidase
		if ($line =~ /,cysteine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cysteine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 203 cystine arylamidase
		if ($line =~ /,cystine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cystine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 204 cytochrome oxidase
		if ($line =~ /,cytochrome oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not cytochrome oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 205 dnase
		if ($line =~ /,dnase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not dnase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 206 esterase c4
		if ($line =~ /,esterase c4,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not esterase c4,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 207 fructose bisphosphate aldolase
		if ($line =~ /,fructose bisphosphate aldolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fructose bisphosphate aldolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 208 fucosidase
		if ($line =~ /,fucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not fucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 209 gelatinase
		if ($line =~ /,gelatinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not gelatinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 210 glucosidase
		if ($line =~ /,glucosidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glucosidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 211 glutamate arylamidase
		if ($line =~ /,glutamate arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 212 glutamate deaminase
		if ($line =~ /,glutamate deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamate deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 213 glutamic acid decarboxylase
		if ($line =~ /,glutamic acid decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamic acid decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 214 glutamine arylamidase
		if ($line =~ /,glutamine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 215 glutamyl glutamic acid arylamidase
		if ($line =~ /,glutamyl glutamic acid arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glutamyl glutamic acid arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 216 glycine arylamidase
		if ($line =~ /,glycine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 217 hippurate hydrolysis
		if ($line =~ /,hippurate hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hippurate hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 218 histidine arylamidase
		if ($line =~ /,histidine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not histidine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 219 hyaluronidase
		if ($line =~ /,hyaluronidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hyaluronidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 220 hydroxy proline arylamidase
		if ($line =~ /,hydroxy proline arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hydroxy proline arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 221 iron-phorphyrins
		if ($line =~ /,iron-phorphyrins,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not iron-phorphyrins,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 222 isoleucine arylamidase
		if ($line =~ /,isoleucine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not isoleucine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 223 l-lactic acid racemase
		if ($line =~ /,l-lactic acid racemase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not l-lactic acid racemase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 224 lecithinase
		if ($line =~ /,lecithinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lecithinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 225 leucine arylamidase
		if ($line =~ /,leucine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 226 leucyl glycine arylamidase
		if ($line =~ /,leucyl glycine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not leucyl glycine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 227 lipase c14
		if ($line =~ /,lipase c14,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase c14,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 228 lipase c8
		if ($line =~ /,lipase c8,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase c8,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 229 lipase
		if ($line =~ /,lipase,|,lipase c8,|,lipase c14,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lipase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 230 lysine arylamidase
		if ($line =~ /,lysine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 231 lysine deaminase
		if ($line =~ /,lysine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 232 lysine decarboxylase
		if ($line =~ /,lysine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lysine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 233 methionine arylamidase
		if ($line =~ /,methionine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methionine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 234 methyl a-d-mannoside
		if ($line =~ /,methyl a-d-mannoside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl a-d-mannoside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 235 methyl b-d-xyloside
		if ($line =~ /,methyl b-d-xyloside,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl b-d-xyloside,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 236 n-acetyl-b-d-glucosaminidase
		if ($line =~ /,n-acetyl-b-d-glucosaminidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not n-acetyl-b-d-glucosaminidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 237 nadh oxidase
		if ($line =~ /,nadh oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not nadh oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 238 naphthol-as-bi-phosphohydrolase
		if ($line =~ /,naphthol-as-bi-phosphohydrolase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not naphthol-as-bi-phosphohydrolase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 239 ornithine arylamidase
		if ($line =~ /,ornithine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 240 ornithine decarboxylase
		if ($line =~ /,ornithine decarboxylase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ornithine decarboxylase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 241 oxidase
		if ($line =~ /,oxidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not oxidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 242 para-phenylalanine deaminase
		if ($line =~ /,para-phenylalanine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not para-phenylalanine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 243 phenylalanine arylamidase
		if ($line =~ /,phenylalanine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 244 phenylalanine deaminase
		if ($line =~ /,phenylalanine deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phenylalanine deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 245 phosphatase
		if ($line =~ /,phosphatase,,|,alkaline phosphatase,|,acid phosphatase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phosphatase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 246 proline arylamidase
		if ($line =~ /,proline arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not proline arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 247 protease
		if ($line =~ /,protease,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not protease,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 248 pyrazinamidase
		if ($line =~ /,pyrazinamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyrazinamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 249 pyroglutamic acid arylamidase
		if ($line =~ /,pyroglutamic acid arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not pyroglutamic acid arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 250 proteolytic enzymes
		if ($line =~ /,caseinase,|,casein hydrolysis,|,chymotrypsin,|,coagulase,|,collagenase,|,gelatinase,|,protease,|,trypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 251 serine arylamidase
		if ($line =~ /,serine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not serine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 252 threonine arylamidase
		if ($line =~ /,threonine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not threonine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 253 trypsin
		if ($line =~ /,trypsin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not trypsin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 254 tryptophan arylamidase
		if ($line =~ /,tryptophan arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 255 tryptophan deaminase
		if ($line =~ /,tryptophan deaminase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophan deaminase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 256 tryptophanase
		if ($line =~ /,tryptophanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tryptophanase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 257 tyrosine arylamidase
		if ($line =~ /,tyrosine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tyrosine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 258 urease
		if ($line =~ /,urease,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not urease,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 259 valine arylamidase
		if ($line =~ /,valine arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not valine arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 260 xylanase
		if ($line =~ /,xylanase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xylanase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 261 a-hemolysis
		if ($line =~ /,a-hemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not a-hemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 262 acetoin
		if ($line =~ /,acetoin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not acetoin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 263 b-haemolysis
		if ($line =~ /,b-haemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not b-haemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 264 ellman
		if ($line =~ /,ellman,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not ellman,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 265 flexirubin
		if ($line =~ /,flexirubin,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not flexirubin,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 266 haemolysis
		if ($line =~ /,haemolysis,|,a-haemolysis,|,b-haemolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not haemolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 267 indoxyl acetate hydrolysis
		if ($line =~ /,indoxyl acetate hydrolysis,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not indoxyl acetate hydrolysis,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 268 methyl red acidification
		if ($line =~ /,methyl red acidification,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methyl red acidification,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 269 methylene blue reduction
		if ($line =~ /,methylene blue reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not methylene blue reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 270 resazurin reduction
		if ($line =~ /,resazurin reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not resazurin reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 271 tetrathionate reduction
		if ($line =~ /,tetrathionate reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetrathionate reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 272 tetrazolium reduction
		if ($line =~ /,tetrazolium reduction,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not tetrazolium reduction,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 273 carboxymethyl cellulase
		if ($line =~ /,carboxymethyl cellulase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not carboxymethyl cellulase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 274 1,3-propanediol
		if ($line =~ /,1-3-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-3-propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 275 1,2-propanediol
		if ($line =~ /,1-2-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not 1-2-propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 276 propanediol
		if ($line =~ /,propanediol,|,1-2-propanediol,|,1-3-propanediol,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not propanediol,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 277 phosphatidylinositol phospholipase c
		if ($line =~ /,phosphatidylinositol phospholipase c,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not phosphatidylinositol phospholipase c,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 278 lecithinase
		if ($line =~ /,lecithinase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lecithinase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 279 lactamide
		if ($line =~ /,lactamide,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lactamide,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 280 hydroxyproline
		if ($line =~ /,hydroxyproline,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not hydroxyproline,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 281 glycyl tryptophan arylamidase
		if ($line =~ /,glycyl tryptophan arylamidase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not glycyl tryptophan arylamidase,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 282 lipolytic enzymes
		if ($line =~ /,esterase c4,|,lecithinase,|,lipase c14,|,lipase c8,|,lipase,|,esterase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 283 esterase
		if ($line =~ /,esterase c4,|,esterase,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not xxx,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}

		print OUT "\n";
		}
	print OUT "\n\;\nEND\;\n";
	print $outmessage;

	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
	}




else {
	print "NOTE*********************************************************************************************************\n";
	print "NOTE********************The Digester is not able to process this character at present\.***********************\n";
	print "NOTE*********************************************************************************************************\n";
	}
