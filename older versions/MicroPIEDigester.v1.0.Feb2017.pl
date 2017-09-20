#!/usr/bin/perl -w

####################################################################
# This script creates a Mesquite readable nexus file from MicroPIE output
# file in tab-delimited text format
# Carrine Blank September through October 2016
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

my $usage = "usage: perl MC2MesquiteDiscrete.pl infile\n";
my $infile = shift or die $usage;	# input Nexus MicroPIE file
my $unixfile = $infile."_unix";

#deposits character names into the array @charlist and asks user to select their desired character to process/code
system ("LC_ALL=C tr '\r' '\n' <$infile >$unixfile");
open(IN, $unixfile) || die "Can't open $unixfile: $!\n";
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
open (IN, '<', $unixfile) or die $!;
open (OUT, '>', $rawmatrix) or die $!;
while (my $line = <IN> ) {
	chomp $line;
	my $rawcharstates;
	my @columns = split /\t/, $line;
	push (my @taxnames, $columns[0]);
	push (my @rawcharstates, $columns[$desiredkey]);
	$_ = lc for @rawcharstates;
	map {s/, /,/g; } @rawcharstates;  #change # to commas
	map {s/#/,/g; } @rawcharstates;  #change # to commas
	map {s/\"//g; } @rawcharstates; #remove double quotes
	map {s/^\s//g; } @rawcharstates; # remove space at beginning of the string
	map {s/\s\s/ /g; } @rawcharstates; # remove double spaces
	map {s/^/,/g; } @rawcharstates; # puts a comma at the front of the string
	map {s/$/,/g; } @rawcharstates; # puts a comma at the end of the string
	map {s/,,/,/g; } @rawcharstates; # removes double commas
	map {s/no\s/not /g; } @rawcharstates; # transforms no to not
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
if ($character eq "median mol %g+c") {
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
	print OUT "1 median_mol_%g+c\;\n\tMATRIX\n";

	open (IN, '<', $rawmatrix) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @molgctable = split /\t/, $line;
		push (my @taxnames, $molgctable[0]);
		push (my @molgcdata, $molgctable[1]);
		map {s/,//g; } @molgcdata;
		local $, = "";
		print OUT @taxnames, " ", @molgcdata, "\n"; # prints to $nexusoutfile molg+c.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
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
		map {s/,helical,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,helical spirals,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,may form coils,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,helix,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,spiral,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,spirillum,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,wavy,/,spirillum,/g; } @cellshapes; # synonyms of spirillum
		map {s/,vibrio,/,vibrioid,/g; } @cellshapes; # synonyms of vibrioid
		map {s/,crooked,/,curved,/g; } @cellshapes; # synonyms of spirillum
		map {s/,curved,/,curved,/g; } @cellshapes; # synonyms of spirillum
		map {s/,curved rods,/,vibrioid,/g; } @cellshapes; # synonyms of spirillum
		map {s/,coccoidal,/,coccus,/g; } @cellshapes; # synonyms of coccus
		map {s/,coccoid,/,coccus,/g; } @cellshapes; # synonyms of coccus
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
		map {s/,filaments,/,filamentous,/g; } @cellshapes; # synonyms of filamentous
		map {s/,flexous,/,flexible,/g; } @cellshapes; # synonyms of filamentous
		map {s/,flexuous,/,flexible,/g; } @cellshapes; # synonyms of filamentous
		map {s/,branches,/,branching,/g; } @cellshapes; # synonyms of filamentous
		map {s/,short chains,/,chains,/g; } @cellshapes; 
		map {s/,irregular,/,not regular,/g; } @cellshapes; 
		#negation of ontology terms and synonyms
		map {s/,not irregular,/,regular,/g; } @cellshapes; 
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Shape Matrix'\;\n\tDIMENSIONS NCHAR=18\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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
	print OUT " \;\n\tMATRIX\n";
#	print OUT "18 charname \/  absence presence, "; # template
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		print OUT "\n";
	}
	print OUT "\n\;\n\nEND\;\n";
	
	unlink $rawmatrix;
	unlink $homout;
	unlink $temp2;
	unlink $temp3;
	unlink $temp4;
	unlink $temp5;
	unlink $temp6;
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
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /Tree Name.*/) {
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
	unlink $rawmatrix;
	}

###################
#process cell relationships character
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
		map {s/,filaments,/,filamentous,/g; } @cellrelates; # synonyms of filamentous
		map {s/,singly,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
		map {s/,single,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
		map {s/,unicellular,/,unicellular,/g; } @cellrelates; # synonyms of unicellular
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Relationships Matrix'\;\n\tDIMENSIONS NCHAR=17\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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
	print OUT " \;\n\tMATRIX\n";

	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @cellrelatesdata = split (/\t/, $line);
		push (my @taxnames, $cellrelatesdata[0]);
#code char 1 unicellularity
		if ($line =~ /,unicellular,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not unicellular,|,multicellular,|,filamentous,|,streptobacillus,|,streptococcus,|,chroococcoid,|,sarcina,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplobacillus,|,diplococcus,|,palisade,|,tetrad,/) {
			print OUT @taxnames, "0";
		}
		else {
			print OUT @taxnames, "?";
		}
#code char 2 multicellularity
		if ($line =~ /,multicellular,|,filamentous,|,streptobacillus,|,steptococcus,|,staphylococcus,|,uniseriate,|,multiseriate,|,diplococcus,|,diplobacillus,/) {
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
		elsif ($line =~ /,multicellular,|,filamentous,|,sarcina,|,multiseriate,|,uniseriate,|,chroococcoid,/) {
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
		print OUT "\n";
		}
	print OUT "\;\nEND\;\n";
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
		map {s/,gram stain negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram staining negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram reaction negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gramnegative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative
		map {s/,gram-staining-negative,/,gram negative,/gi; } @gramstain; # synonyms of gram negative

		map {s/,gram stain positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram staining positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram reaction positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,grampositive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive
		map {s/,gram-staining-positive,/,gram positive,/gi; } @gramstain; # synonyms of gram positive

		map {s/,gram stain variable,/,gram variable,/gi; } @gramstain; # synonyms of gram variable
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
		if ($line =~ /Tree Name.*/) {
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
		map {s/,central granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,course granules,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,fine granules,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,fine granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,globule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,granular inclusions,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,granular,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,granulated,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,intracellular globule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,large granules,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,large granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,numerous granules,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,single granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,singular granule,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,small electron-dense granules,/,granules,/g; } @intfeats; # synonyms of 
		map {s/,solitary granules,/,granules,/g; } @intfeats; # synonyms of 

		map {s/,metachromatic granules,/,polyphosphate granules,/g; } @intfeats; # synonyms of 
		map {s/,metachromatic granule,/,polyphosphate granules,/g; } @intfeats; # synonyms of 
		map {s/,polyphosphate granule,/,polyphosphate granules,/g; } @intfeats; # synonyms of 
		map {s/,polyphosphate,/,polyphosphate granules,/g; } @intfeats; # synonyms of 
		map {s/,volutin granules,/,polyphosphate granules,/g; } @intfeats; # synonyms of 
		map {s/,volutin granule,/,polyphosphate granules,/g; } @intfeats; # synonyms of 

		map {s/,intracellular sulfur globule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,intracellular sulphur globule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulfur globule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulfur granule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulfur granule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulfur granules,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulfur granules,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulphur globule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulphur granule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulphur granule,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulphur granules,/,sulfur granules,/g; } @intfeats; # synonyms of 
		map {s/,sulphur granules,/,sulfur granules,/g; } @intfeats; # synonyms of 


		map {s/,PHB,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxybutyrate granules,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,polyhydroxybutyrate,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,polyhydroxybutyric acid,/,poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,starch,/,starch granules,/g; } @intfeats; # synonyms of 
		map {s/,PHA,/,poly-b-hydroxyalkanoate granules,/g; } @intfeats; # synonyms of 
		map {s/,poly-b-hydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/g; } @intfeats; # synonyms of 
		map {s/,poly-beta-hydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/g; } @intfeats; # synonyms of 
		map {s/,polyhydroxyalkanoate,/,poly-b-hydroxyalkanoate granules,/g; } @intfeats; # synonyms of 
		map {s/,glycogen granule,/,glycogen granules,/g; } @intfeats; # synonyms of 
		map {s/,glycogen,/,glycogen granules,/g; } @intfeats; # synonyms of 
		map {s/,glycogen particles,/,glycogen granules,/g; } @intfeats; # synonyms of 
		map {s/,glycogen particle,/,glycogen granules,/g; } @intfeats; # synonyms of 
		map {s/,cyanophycin,/,cyanophycin granules,/g; } @intfeats; # synonyms of 

		map {s/,pairs gas vesicles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,gas vesicles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,vacuolated cells,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,vesicle-like microstructures,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,small vesicles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,small vacuole,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,small vacuoles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,small vesicles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,large gas vacuole,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,large vacuole,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,large vacuoles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,large vesciles,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,gas vacuolation,/,gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,vacuoles,/,gas vacuoles,/g; } @intfeats; # synonyms of 

		map {s/,resting cells,/,resting stages,/g; } @intfeats; # synonyms of 
		map {s/,fruiting bodies,/,resting stages,/g; } @intfeats; # synonyms of 
		map {s/,resting spores,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,spore-like structures,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,spore-former,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,spore-forming,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,sporulating,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,sporeforming,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,sporing,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,sporulating,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,sporeforming,/,spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,heat-resistant endospores,/,endospores,/g; } @intfeats; # synonyms of 
		map {s/,myxospore,/,myxospores,/g; } @intfeats; # synonyms of 
		map {s/,spores,/,spore-forming,/g; } @intfeats; # synonyms of 

		map {s/,not vacuoles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not spore-like structures,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,not resting cells,/,not resting stages,/g; } @intfeats; # synonyms of 
		map {s/,absent resting stages,/,not resting stages,/g; } @intfeats; # synonyms of 
		map {s/,asporogenic,/,not resting stages,/g; } @intfeats; # synonyms of 
		map {s/,devoid of gas vacuoles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,devoid of vacuoles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,non-spore-forming,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,non-sporeforming,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,non-sporing,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,non-sporulating,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,nonspore-forming,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,nonsporeforming,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,nonsporing,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,nonsporulating,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,not resting spores,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,not spores,/,not spore-forming,/g; } @intfeats; # synonyms of 
		map {s/,not heat-resistant endospores,/,not endospores,/g; } @intfeats; # synonyms of 
		map {s/,not fruiting bodies,/,not resting stages,/g; } @intfeats; # synonyms of 

		map {s/,not gas vesicles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not pairs gas vesicles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not gas vesicles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not vacuolated cells,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not vesicle-like microstructures,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not small vesicles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not small vacuole,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not small vacuoles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not small vesicles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not large gas vacuole,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not large vacuole,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not large vacuoles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not large vesciles,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
		map {s/,not poly-b-hydroxybutyrate,/,not poly-b-hydroxybutyrate granules,/g; } @intfeats; # synonyms of 
		map {s/,not gas vacuolation,/,not gas vacuoles,/g; } @intfeats; # synonyms of 
#not internal features
		map {s/,carbon compound,/,/g; } @intfeats; 

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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /,endospores,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not endospores,/) {
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
		map {s/,gliding,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,gliding cells,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,gliding motilities,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,gliding types,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,mature gliding stage,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,motile by gliding,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,motile cells,/,motile,/g; } @motility; # synonyms of 
		map {s/,peculiar gliding motility,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,showing gliding motility,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,slow gliding motility,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,swimming cells,/,flagellar motility,/g; } @motility; # synonyms of 
		map {s/,axial fibrils,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,axial filaments,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,endoflagella,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,endoflagellum,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,periplasmic fibrils,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,periplasmic flagella,/,periplasmic flagella,/g; } @motility; # synonyms of 
		map {s/,tumbling,/,tumbling motility,/g; } @motility; # synonyms of 
		map {s/,twitching,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,twitching motility,/,gliding motility,/g; } @motility; # synonyms of 
		map {s/,monopolar flagella,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,monopolar flagellum,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,monotrichous flagella,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,polar flagellum,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,single flagella,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,single flagellum,/,monotrichous,/g; } @motility; # synonyms of 
		map {s/,peritrichously,/,peritrichous,/g; } @motility; # synonyms of 
		map {s/,polar flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,monopolar polytrichous flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,polar bundle of flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,polar tuft of flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,tuft of polar flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,bundles of flagella,/,lophotrichous,/g; } @motility; # synonyms of 
		map {s/,tufts of polar flagella,/,amphilophotrichous,/g; } @motility; # synonyms of 
		map {s/,numerous filamentous tufts,/,amphilophotrichous,/g; } @motility; # synonyms of 
		map {s/,polar tufts of flagella,/,amphilophotrichous,/g; } @motility; # synonyms of 
		map {s/,bipolar polytrichous flagella,/,amphilophotrichous,/g; } @motility; # synonyms of 

		map {s/,non-gliding,/,not gliding motility,/g; } @motility; # synonyms of 
		map {s/,absent flagellar,/,not flagellar motility,/g; } @motility; # synonyms of 
		map {s/,absent flagellar motility,/,not flagellar motility,/g; } @motility; # synonyms of 
		map {s/,absent gliding,/,not gliding motility,/g; } @motility; # synonyms of 
		map {s/,devoid of flagellar,/,not flagellar motility,/g; } @motility; # synonyms of 
		map {s/,devoid of flagellar motility,/,not flagellar motility,/g; } @motility; # synonyms of 
		map {s/,immotile,/,not motile,/g; } @motility; # synonyms of 
		map {s/,never motility,/,not motile,/g; } @motility; # synonyms of 
		map {s/,non-motile,/,not motile,/g; } @motility; # synonyms of 
		map {s/,non-motile cells,/,not motile,/g; } @motility; # synonyms of 
		map {s/,non motile,/,not motile,/g; } @motility; # synonyms of 
		map {s/,non-motile non-gliding,/,not motile,not gliding motility,/g; } @motility; # synonyms of 
		map {s/,nonmotile,/,not motile,/g; } @motility; # synonyms of 
		map {s/,absent gliding motility,/,not gliding motility,/g; } @motility; # synonyms of 
		map {s/,not swimming cells,/,not flagellar motility,/g; } @motility; # synonyms of 
		map {s/,not motility,/,not motile,/g; } @motility; # synonyms of 

		map {s/,lack flagella,/,not flagella,/g; } @motility; # synonyms of 
		map {s/,non-flagellated,/,not flagella,/g; } @motility; # synonyms of 
		map {s/,non-flagellated and non- motile,/,not flagella,not motile,/g; } @motility; # synonyms of 
		map {s/,non-flagellated and non-motile,/,not flagella,not motile,/g; } @motility; # synonyms of 
		map {s/,non-flagellated non-motile,/,not flagella,not motile,/g; } @motility; # synonyms of 
		map {s/,nonflagellated,/,not flagella,/g; } @motility; # synonyms of 
		map {s/,not gliding types,/,not gliding motility,/g; } @motility; # synonyms of 
		map {s/,not non-motile,/,motile,/g; } @motility; # synonyms of 
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Cell Motility Matrix'\;\n\tDIMENSIONS NCHAR=11\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @motilitydata = split (/\t/, $line);
		push (my @taxnames, $motilitydata[0]);

#code char 1 motility
		if ($line =~ /,motile,|,flagellar motility,|,gliding motility,|,tumbling motility,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,/) {
			print OUT @taxnames, "1";
			}
		elsif ($line =~ /,not motile,/) {
			print OUT @taxnames, "0";
			}
		else {
			print OUT @taxnames, "?";
			}
#code char 2 flagellar motility
		if ($line =~ /,flagellar motility,/) {
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
		if ($line =~ /,flagella,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,/) {
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
		elsif ($line =~ /,not monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 7 polytrichous
		if ($line =~ /,polytrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not polytrichous,|,monotrichous,|,lophotrichous,|,amphitrichous,|,amphilophotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 8 lophotrichous
		if ($line =~ /,lophotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not lophotrichous,|,monotrichous,|,polytrichous,|,amphitrichous,|,amphilophotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 9 amphitrichous
		if ($line =~ /,amphitrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphitrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphilophotrichous,|,not motile,|,not flagella,/) {
			print OUT "0";
			}
		else {
			print OUT "?";
			}
#code char 10 amphilophotrichous
		if ($line =~ /,amphilophotrichous,/) {
			print OUT "1";
			}
		elsif ($line =~ /,not amphilophotrichous,|,monotrichous,|,polytrichous,|,lophotrichous,|,amphitrichous,|,not motile,|,not flagella,/) {
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @naclmintable = split (/\t/, $line);
		push (my @taxnames, $naclmintable[0]);
		push (my @naclmindata, $naclmintable[1]);
		local $, = "";
		map {s/plus.*//g; } @naclmindata; # get rid of any additional values
		map {s/,//g; } @naclmindata; # get rid of any additional values
		#convert all concentration values to % NaCl and get rid of units
		foreach (@naclmindata) {
			if (@naclmindata = /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@naclmindata = "@val1" + 0.0;
				}	
			elsif (@naclmindata = /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@naclmindata = 0.05;
				}	
			elsif (@naclmindata = /not required|absence of NaCl|without addition of NaCl/i ) {
				@naclmindata = 0.0;
				}	
			elsif (@naclmindata = /sea\ssalt/) { # converts PERCENT SEA SALTS to % NaCl
				if (@naclmindata = /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @naclmindata; # get rid of the first comma
					map {s/\,.*//; } @naclmindata; # get rid of any additional values
					map {s/%%.*//; } @naclmindata; # get rid of any additional values
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmindata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif (@naclmindata = /-/) { # % sea salt values with a range
					map {s/\,//; } @naclmindata; # get rid of the first comma
					map {s/\,.*//; } @naclmindata; # get rid of any additional values
					map {s/%.*//; } @naclmindata; # get rid of any additional values
					my @medianvals = split /-/, $_;
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
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmindata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif (@naclmindata = /\% NaCl|\%NaCl|\%w\/v NaCl|\%w\/vNaCl|\%w\/v|\% w\/v nacl/i) { # PERCENT NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/\%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						@naclmindata = median(@temp3);
						}
					else { 
						@naclmindata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif (@naclmindata = / mm/i) { # convert MILLIMOLAR concentrations to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/mm.*//g; } @val1;
				my @val2 = "@val1" + 0.0;
				my @val3 = "@val2" * .001;
				my @val4 = "@val3" * 6.6;
				@naclmindata = "@val4" - 0.7546;
				}
			elsif (@naclmindata = / m/i) { # convert MOLAR concentrations to % NaCl
				if (@naclmindata = /mgcl2|mgso4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmindata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@naclmindata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								my @medianvals = split /-/, $_;
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0;
								my @val4 = "@val3" * 6.6;
								@naclmindata = "@val4" - 0.7546;
								}
							else {
								map {s/m\s.*//g; } @val1;
								my @val2 = "@val1" + 0.0;
								my @val3 = "@val2" * 6.6;
								@naclmindata = "@val3" - 0.7546;
								}
							}
						}
					}
				}
			elsif (@naclmindata = /gl|g l|g salt l/) { # convert G PER L to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/gl.*//g; } @val1;
				map {s/g\sl.*//g; } @val1;
				map {s/g\ssalt\sl.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) {
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@naclmindata = "@val3" * 0.10;
						}
					else {
						my @val2 = "@val1" + 0.0;
						@naclmindata = "@val2" * 0.10;
						}
					}
				}
			elsif (@naclmindata = /seawater/) { #convert SEAWATER DILUTIONS to % Nacl
					if (@naclmindata = /x|%/) {
						if (@naclmindata = /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = split /\//, $_;
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0;
							map {s/x.*//; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0;
							my @temp3 = "@temp1"/"@temp2";
							@naclmindata = "@temp3" * 3.5 ;
							}
						elsif (@naclmindata = /x/) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = split /\s/, $_;
							push (my @val1, $naclvals[0]);
							map {s/x.*//; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@naclmindata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						}
					else { # has elsif (@naclmindata = /^seawater$/) # pure seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@naclmindata = "@val1" + 0.0;
						}
				}
			elsif (@naclmindata = /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if (@naclmindata = /-/) { # RANGE of percent values
					my @medianvals = split /-/, $_;
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@naclmindata = median(@temp3);
					}
				else { 
					my @naclvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @naclopttable = split (/\t/, $line);
		push (my @taxnames, $naclopttable[0]);
		push (my @nacloptdata, $naclopttable[1]);
		map {s/,//; } @nacloptdata; # get rid of the first comma
		map {s/,.*//; } @nacloptdata; # get rid of any additional values
		map {s/plus.*//g; } @nacloptdata; # get rid of any additional values
		local $, = " ";
		#convert all concentration values to % NaCl and get rid of units
		foreach (@nacloptdata) {
			if (@nacloptdata = /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@nacloptdata = "@val1" + 0.0;
				}	
			elsif (@nacloptdata = /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@nacloptdata = 0.05;
				}	
			elsif (@nacloptdata = /not required|absence of NaCl|without addition of NaCl/i ) {
				@nacloptdata = 0.0;
				}	
			elsif (@nacloptdata = /sea\ssalt/) { # converts PERCENT SEA SALTS to % NaCl
				if (@nacloptdata = /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @nacloptdata; # get rid of the first comma
					map {s/\,.*//; } @nacloptdata; # get rid of any additional values
					map {s/%%.*//; } @nacloptdata; # get rid of any additional values
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@nacloptdata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif (@nacloptdata = /-/) { # % sea salt values with a range
					map {s/\,//; } @nacloptdata; # get rid of the first comma
					map {s/\,.*//; } @nacloptdata; # get rid of any additional values
					map {s/%.*//; } @nacloptdata; # get rid of any additional values
					my @medianvals = split /-/, $_;
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
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@nacloptdata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif (@nacloptdata = /\% NaCl|\%NaCl|\%w\/v NaCl|\%w\/vNaCl|\%w\/v|\% w\/v nacl/i) { # PERCENT NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/\%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						@nacloptdata = median(@temp3);
						}
					else { 
						@nacloptdata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif (@nacloptdata = / mm/i) { # convert MILLIMOLAR concentrations to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/mm.*//g; } @val1;
				my @val2 = "@val1" + 0.0;
				my @val3 = "@val2" * .001;
				my @val4 = "@val3" * 6.6;
				@nacloptdata = "@val4" - 0.7546;
				}
			elsif (@nacloptdata = / m/i) { # convert MOLAR concentrations to % NaCl
				if (@nacloptdata = /mgcl2|mgso4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@nacloptdata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@nacloptdata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								my @medianvals = split /-/, $_;
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0;
								my @val4 = "@val3" * 6.6;
								@nacloptdata = "@val4" - 0.7546;
								}
							else {
								map {s/m\s.*//g; } @val1;
								my @val2 = "@val1" + 0.0;
								my @val3 = "@val2" * 6.6;
								@nacloptdata = "@val3" - 0.7546;
								}
							}
						}
					}
				}
			elsif (@nacloptdata = /gl|g l|g salt l/) { # convert G PER L to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/gl.*//g; } @val1;
				map {s/g\sl.*//g; } @val1;
				map {s/g\ssalt\sl.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) {
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@nacloptdata = "@val3" * 0.10;
						}
					else {
						my @val2 = "@val1" + 0.0;
						@nacloptdata = "@val2" * 0.10;
						}
					}
				}
			elsif (@nacloptdata = /seawater/) { #convert SEAWATER DILUTIONS to % Nacl
					if (@nacloptdata = /x|%/) {
						if (@nacloptdata = /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = split /\//, $_;
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0;
							map {s/x.*//; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0;
							my @temp3 = "@temp1"/"@temp2";
							@nacloptdata = "@temp3" * 3.5 ;
							}
						elsif (@nacloptdata = /x/) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = split /\s/, $_;
							push (my @val1, $naclvals[0]);
							map {s/x.*//; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@nacloptdata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						}
					else { # has elsif (@nacloptdata = /^seawater$/) # pure seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@nacloptdata = "@val1" + 0.0;
						}
				}
			elsif (@nacloptdata = /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if (@nacloptdata = /-/) { # RANGE of percent values
					my @medianvals = split /-/, $_;
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@nacloptdata = median(@temp3);
					}
				else { 
					my @naclvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @naclmaxtable = split (/\t/, $line);
		push (my @taxnames, $naclmaxtable[0]);
		push (my @naclmaxdata, $naclmaxtable[1]);
		map {s/,//; } @naclmaxdata; # get rid of the first comma
		map {s/,.*//; } @naclmaxdata; # get rid of any additional values
		map {s/plus.*//g; } @naclmaxdata; # get rid of any additional values
		local $, = " ";
		#convert all concentration values to % NaCl and get rid of units
		foreach (@naclmaxdata) {
			if (@naclmaxdata = /not in absence|not 0|not in the absence of nacl/i ) { # this is an arbitrary setting
				my @val1 = 0.01;
				@naclmaxdata = "@val1" + 0.0;
				}	
			elsif (@naclmaxdata = /not freshwater/i ) { # sets salinity to be at the highest level of freshwater salinity
				@naclmaxdata = 0.05;
				}	
			elsif (@naclmaxdata = /not required|absence of NaCl|without addition of NaCl/i ) {
				@naclmaxdata = 0.0;
				}	
			elsif (@naclmaxdata = /sea\ssalt/) { # converts PERCENT SEA SALTS to % NaCl
				if (@naclmaxdata = /\%\%/) { # per mil sea salts concentration
					map {s/\,//; } @naclmaxdata; # get rid of the first comma
					map {s/\,.*//; } @naclmaxdata; # get rid of any additional values
					map {s/%%.*//; } @naclmaxdata; # get rid of any additional values
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/%%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmaxdata = "@val2" * 0.084; # sea salt contains ~8.4%% NaCl
					}
				elsif (@naclmaxdata = /-/) { # % sea salt values with a range
					map {s/\,//; } @naclmaxdata; # get rid of the first comma
					map {s/\,.*//; } @naclmaxdata; # get rid of any additional values
					map {s/%.*//; } @naclmaxdata; # get rid of any additional values
					my @medianvals = split /-/, $_;
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
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%//g; } @val1;
					my @val2 = "@val1" + 0.0;
					@naclmaxdata = "@val2" * 0.84; # sea salt contains ~84% w/w NaCl
					}
				}
			elsif (@naclmaxdata = /\% NaCl|\%NaCl|\%w\/v NaCl|\%w\/vNaCl|\%w\/v|\% w\/v nacl/i) { # PERCENT NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/\%.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) { # range of % NaCl values
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						@naclmaxdata = median(@temp3);
						}
					else { 
						@naclmaxdata = "@val1" + 0.0; # single % NaCl value
						}
					}
				}
			elsif (@naclmaxdata = / mm/i) { # convert MILLIMOLAR concentrations to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/mm.*//g; } @val1;
				my @val2 = "@val1" + 0.0;
				my @val3 = "@val2" * .001;
				my @val4 = "@val3" * 6.6;
				@naclmaxdata = "@val4" - 0.7546;
				}
			elsif (@naclmaxdata = / m/i) { # convert MOLAR concentrations to % NaCl
				if (@naclmaxdata = /mgcl2|mgso4/i ) { # assumption: treat % MgCl or % MgSO4 as % NaCl 
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/\%.*//g; } @val1; # 
					@naclmaxdata = "@val1" + 0.0; # single % "NaCl" value
					}
				else {
					my @naclvals = split /\s/, $_;
					push (my @val1, $naclvals[0]);
					map {s/m\s.*//g; } @val1;
					if ("@val1" eq 0) { # if value is 0 M
						@naclmaxdata = 0.0;
						}
					else { # if value is not 0 M
						foreach (@naclvals) {
							if (@naclvals = /-/) {
								my @medianvals = split /-/, $_;
								push (my @mval1, $medianvals[0]);
								push (my @mval2, $medianvals[1]);
								my @temp1 = "@mval1" + 0.0;
								my @temp2 = "@mval2" + 0.0;
								my @temp3 = (@temp1, @temp2);
								my @val2 = median(@temp3);
								my @val3 = "@val2" + 0.0;
								my @val4 = "@val3" * 6.6;
								@naclmaxdata = "@val4" - 0.7546;
								}
							else {
								map {s/m\s.*//g; } @val1;
								my @val2 = "@val1" + 0.0;
								my @val3 = "@val2" * 6.6;
								@naclmaxdata = "@val3" - 0.7546;
								}
							}
						}
					}
				}
			elsif (@naclmaxdata = /gl|g l|g salt l/) { # convert G PER L to % NaCl
				my @naclvals = split /\s/, $_;
				push (my @val1, $naclvals[0]);
				map {s/gl.*//g; } @val1;
				map {s/g\sl.*//g; } @val1;
				map {s/g\ssalt\sl.*//g; } @val1;
				foreach (@naclvals) {
					if (@naclvals = /-/) {
						my @medianvals = split /-/, $_;
						push (my @mval1, $medianvals[0]);
						push (my @mval2, $medianvals[1]);
						my @temp1 = "@mval1" + 0.0;
						my @temp2 = "@mval2" + 0.0;
						my @temp3 = (@temp1, @temp2);
						my @val2 = median(@temp3);
						my @val3 = "@val2" + 0.0;
						@naclmaxdata = "@val3" * 0.10;
						}
					else {
						my @val2 = "@val1" + 0.0;
						@naclmaxdata = "@val2" * 0.10;
						}
					}
				}
			elsif (@naclmaxdata = /seawater/) { #convert SEAWATER DILUTIONS to % Nacl
					if (@naclmaxdata = /x|%/) {
						if (@naclmaxdata = /\//) { # when the seawater concentration is a fraction like 2/3x
							my @medianvals = split /\//, $_;
							push (my @mval1, $medianvals[0]);
							push (my @mval2, $medianvals[1]);
							my @temp1 = "@mval1" + 0.0;
							map {s/x.*//; } @mval2; # remove the x
							my @temp2 = "@mval2" + 0.0;
							my @temp3 = "@temp1"/"@temp2";
							@naclmaxdata = "@temp3" * 3.5 ;
							}
						elsif (@naclmaxdata = /x/) { # when the seawater concentration is a whole number like 2x or 0.5x
							my @naclvals = split /\s/, $_;
							push (my @val1, $naclvals[0]);
							map {s/x.*//; } @val1; # remove the x
							my @val2 = "@val1" + 0.0;
							@naclmaxdata = "@val2" * 3.5; # seawater is 3.5% salinity
							}
						}
					else { # has elsif (@naclmaxdata = /^seawater$/) # pure seawater
						my @val1 = 3.5; # seawater is 3.5% salinity
						@naclmaxdata = "@val1" + 0.0;
						}
				}
			elsif (@naclmaxdata = /%/) { # reports a PERCENT but lacking units - assume this is % NaCl
				if (@naclmaxdata = /-/) { # RANGE of percent values
					my @medianvals = split /-/, $_;
					push (my @mval1, $medianvals[0]);
					push (my @mval2, $medianvals[1]);
					map {s/\%.*//g; } @mval2; # 
					my @temp1 = "@mval1" + 0.0;
					my @temp2 = "@mval2" + 0.0;
					my @temp3 = (@temp1, @temp2);
					@naclmaxdata = median(@temp3);
					}
				else { 
					my @naclvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @pHmintable = split (/\t/, $line);
		push (my @taxnames, $pHmintable[0]);
		push (my @pHmindata, $pHmintable[1]);
		map {s/,not.*,/,/g; } @pHmindata; # synonyms of 
		map {s/,//; } @pHmindata;
		map {s/,\d.\d-\d.\d,/,/g; } @pHmindata; # synonyms of 
		map {s/,\d.\d//g; } @pHmindata; # synonyms of 
		map {s/around neutral/7.0/g; } @pHmindata; # synonyms of 
		map {s/,//g; } @pHmindata;
		map {s/-/ /g; } @pHmindata;
		foreach (@pHmindata) {
			if ($_ =~ /\d.\d\s\d.\d/) {
				my @medianvals = split /\s/, $_;
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@pHmindata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHmindata, "\n"; # prints to $nexusoutfile pHmin.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @pHopttable = split (/\t/, $line);
		push (my @taxnames, $pHopttable[0]);
		push (my @pHoptdata, $pHopttable[1]);
		#ontology terms and synonyms
		map {s/,not.*,/,/g; } @pHoptdata; # synonyms of 
		map {s/,//; } @pHoptdata;
		map {s/,\d.\d-\d.\d,/,/g; } @pHoptdata; # synonyms of 
		map {s/,\d.\d//g; } @pHoptdata; # synonyms of 
		map {s/around neutral/7.0/g; } @pHoptdata; # synonyms of 
		map {s/,//g; } @pHoptdata;
		map {s/-/ /g; } @pHoptdata;
		foreach (@pHoptdata) {
			if ($_ =~ /\d.\d\s\d.\d/) {
				my @medianvals = split /\s/, $_;
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@pHoptdata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHoptdata, "\n"; # prints to $nexusoutfile pHopt.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @pHmaxtable = split (/\t/, $line);
		push (my @taxnames, $pHmaxtable[0]);
		push (my @pHmaxdata, $pHmaxtable[1]);
		#ontology terms and synonyms
		map {s/,not.*,/,/g; } @pHmaxdata; # synonyms of 
		map {s/,//; } @pHmaxdata;
		map {s/,\d.\d-\d.\d,/,/g; } @pHmaxdata; # synonyms of 
		map {s/,\d.\d//g; } @pHmaxdata; # synonyms of 
		map {s/around neutral/7.0/g; } @pHmaxdata; # synonyms of 
		map {s/,//g; } @pHmaxdata;
		map {s/-/ /g; } @pHmaxdata;
		foreach (@pHmaxdata) {
			if ($_ =~ /\d.\d\s\d.\d/) {
				my @medianvals = split /\s/, $_;
				push (my @val1, $medianvals[0]);
				push (my @val2, $medianvals[1]);
				my @temp1 = "@val1" + 0.0;
				my @temp2 = "@val2" + 0.0;
				my @temp3 = (@temp1, @temp2);
				@pHmaxdata = median(@temp3);
				}
			}
		local $, = "";
		print OUT @taxnames, " ", @pHmaxdata, "\n"; # prints to $nexusoutfile pHmax.nex
		}
	print OUT "\n\n\;\n\nEND\;\n";
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @tempmintable = split (/\t/, $line);
		push (my @taxnames, $tempmintable[0]);
		push (my @tempmindata, $tempmintable[1]);
		#ontology terms and synonyms
		map {s/,not.*,/,/g; } @tempmindata; 
		map {s/,//; } @tempmindata;
		map {s/,\S.*/,/g; } @tempmindata; # synonyms of 
		map {s/,\d.*/,/g; } @tempmindata; # synonyms of 
		map {s/_c//g; } @tempmindata; # synonyms of 
		map {s/,//g; } @tempmindata;
		foreach (@tempmindata) {
			map {s/-/ /g; } @tempmindata;
			if ($_ =~ /\d+\s\d+/) {
				my @medianvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @tempopttable = split (/\t/, $line);
		push (my @taxnames, $tempopttable[0]);
		push (my @tempoptdata, $tempopttable[1]);
		#ontology terms and synonyms
		map {s/,not.*,/,/g; } @tempoptdata; 
		map {s/,//; } @tempoptdata;
		map {s/,\S.*/,/g; } @tempoptdata; # synonyms of 
		map {s/,\d.*/,/g; } @tempoptdata; # synonyms of 
		map {s/_c//g; } @tempoptdata; # synonyms of 
		map {s/,//g; } @tempoptdata;
		foreach (@tempoptdata) {
			map {s/-/ /g; } @tempoptdata;
			if ($_ =~ /\d+\s\d+/) {
				my @medianvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @tempmaxtable = split (/\t/, $line);
		push (my @taxnames, $tempmaxtable[0]);
		push (my @tempmaxdata, $tempmaxtable[1]);
		#ontology terms and synonyms
		map {s/,not.*,/,/g; } @tempmaxdata; 
		map {s/,//; } @tempmaxdata;
		map {s/,\S.*/,/g; } @tempmaxdata; # synonyms of 
		map {s/,\d.*/,/g; } @tempmaxdata; # synonyms of 
		map {s/_c//g; } @tempmaxdata; # synonyms of 
		map {s/,//g; } @tempmaxdata;
		foreach (@tempmaxdata) {
			map {s/-/ /g; } @tempmaxdata;
			if ($_ =~ /\d+\s\d+/) {
				my @medianvals = split /\s/, $_;
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
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /Tree Name.*/) {
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
		map {s/,2% chicken serum,/,serum,/g; } @vitcos; # synonyms of 
		map {s/,1% ox-bile salts,/,ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,ox bile,/,ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,ox bile salts,/,ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,ox-bile salts,/,ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,20% bile,/,bile,/g; } @vitcos; # synonyms of 
		map {s/,40% bile,/,bile,/g; } @vitcos; # synonyms of 
		map {s/,bile salts,/,bile,/g; } @vitcos; # synonyms of 
		map {s/,medium containing 20% bile,/,bile,/g; } @vitcos; # synonyms of 
		map {s/,both vitamin b12,/,vitamin b12,/g; } @vitcos; # synonyms of 
		map {s/,cobalamin,/,vitamin b12,/g; } @vitcos; # synonyms of 
		map {s/,organic growth factors,/,growth factor,/g; } @vitcos; # synonyms of 
		map {s/,exogenous growth factor,/,growth factor,/g; } @vitcos; # synonyms of 
		map {s/,growth factor dependency,/,growth factor,/g; } @vitcos; # synonyms of 
		map {s/,growth factors,/,growth factor,/g; } @vitcos; # synonyms of 
		map {s/,hemin and coenzyme I,/,hemin,nad,/g; } @vitcos; # synonyms of 
		map {s/,x.* + v factor,/,hemin,nad,/g; } @vitcos; # synonyms of 
		map {s/,x+v factor disc.*,/,hemin,nad,/g; } @vitcos; # synonyms of 
		map {s/,xv discs,/,hemin,nad,/g; } @vitcos; # synonyms of 
		map {s/,xv disks,/,hemin,nad,/g; } @vitcos; # synonyms of 
		map {s/,factor x,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,haem,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,haemin,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,x disc,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,x disk,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,x factor,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,x growth factor,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,x-factor,/,hemin,/g; } @vitcos; # synonyms of 
		map {s/,factor v,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,v disc,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,v disk,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,v factor,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,v growth factor,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,v-factor,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,coenzyme I,/,nad,/g; } @vitcos; # synonyms of 
		map {s/,vitamin k,/,vitamin k1,/g; } @vitcos; # synonyms of 
		map {s/,thiamine,/,thiamin,/g; } @vitcos; # synonyms of 
		map {s/,vitamin b1,/,thiamin,/g; } @vitcos; # synonyms of 
		map {s/,vitamin b-1,/,thiamin,/g; } @vitcos; # synonyms of 
		map {s/,vitamin b-12,/,vitamin b12,/g; } @vitcos; # synonyms of 
 
		map {s/,not thiamine,/,not thiamin,/g; } @vitcos; # synonyms of 
		map {s/,not vitamin b1,/,not thiamin,/g; } @vitcos; # synonyms of 
		map {s/,not vitamin b-1,/,not thiamin,/g; } @vitcos; # synonyms of 
		map {s/,not vitamin b-12,/,not vitamin b12,/g; } @vitcos; # synonyms of 
		map {s/,not vitamin k,/,not vitamin k1,/g; } @vitcos; # synonyms of 
		map {s/,not 2% chicken serum,/,not serum,/g; } @vitcos; # synonyms of 
		map {s/,not 1% ox-bile salts,/,not ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,not ox bile,/,not ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,not ox bile salts,/,not ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,not ox-bile salts,/,not ox-bile,/g; } @vitcos; # synonyms of 
		map {s/,not 20% bile,/,not bile,/g; } @vitcos; # synonyms of 
		map {s/,not 40% bile,/,not bile,/g; } @vitcos; # synonyms of 
		map {s/,not bile salts,/,not bile,/g; } @vitcos; # synonyms of 
		map {s/,not medium containing 20% bile,/,not bile,/g; } @vitcos; # synonyms of 
		map {s/,not both vitamin b12,/,not vitamin b12,/g; } @vitcos; # synonyms of 
		map {s/,not cobalamin,/,not vitamin b12,/g; } @vitcos; # synonyms of 
		map {s/,not organic growth factors,/,not growth factor,/g; } @vitcos; # synonyms of 
		map {s/,not exogenous growth factor,/,not growth factor,/g; } @vitcos; # synonyms of 
		map {s/,not growth factor dependency,/,not growth factor,/g; } @vitcos; # synonyms of 
		map {s/,not growth factors,/,not growth factor,/g; } @vitcos; # synonyms of 
		map {s/,not hemin and coenzyme I,/,not hemin,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not x.* + v factor,/,not hemin,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not x+v factor disc.*,/,not hemin,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not xv discs,/,not hemin,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not xv disks,/,not hemin,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not factor x,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not haem,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not haemin,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not x disc,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not x disk,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not x factor,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not x growth factor,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not x-factor,/,not hemin,/g; } @vitcos; # synonyms of 
		map {s/,not factor v,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not v disc,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not v disk,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not v factor,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not v growth factor,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not v-factor,/,not nad,/g; } @vitcos; # synonyms of 
		map {s/,not coenzyme I,/,not nad,/g; } @vitcos; # synonyms of 

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
	map {s/^thiamin/**thiamin/g; } @charlist2;
	map {s/^vitamin b12/**vitamin b12/g; } @charlist2;
	map {s/^vitamin k1/**vitamin k1/g; } @charlist2;
	map {s/^vitamins/**vitamins/g; } @charlist2;
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Vitamins And Cofactors Matrix'\;\n\tDIMENSIONS NCHAR=10\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
			next;
			}
		my @vitcosdata = split (/\t/, $line);
		push (my @taxnames, $vitcosdata[0]);

#code char 1 vitamins
		if ($line =~ /,vitamins,|,vitamin k1,|,vitamin b12,|,thiamin,/) {
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
		map {s/,acetylspiramycin,/,spiramycin,/g; } @antisens; # synonyms of 
		map {s/,dactinomycin,/,actinomycin,/g; } @antisens; # synonyms of 
		map {s/,actinomycin d,/,actinomycin,/g; } @antisens; # synonyms of 
		map {s/,amoxycillin,/,amoxicillin,/g; } @antisens; # synonyms of 
		map {s/,amoxicillin\/clavulanic acid,/,amoxicillin,clavulanic acid,/g; } @antisens; # synonyms of 
		map {s/,benzylphenicillin,/,benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,penicillin g,/,benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,malachite green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,emerald green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,solid green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,diamond green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,aniline green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,benzaldehyde green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,fast green,/,brilliant green,/g; } @antisens; # synonyms of 
		map {s/,cefazoline,/,cefazolin,/g; } @antisens; # synonyms of 
		map {s/,cephazolin,/,cefazolin,/g; } @antisens; # synonyms of 
		map {s/,cefobid,/,cefoperazone,/g; } @antisens; # synonyms of 
		map {s/,mefoxin,/,cefoxitin,/g; } @antisens; # synonyms of 
		map {s/,cephaloridin,/,cephaloridine,/g; } @antisens; # synonyms of 
		map {s/,cefaloridine,/,cephaloridine,/g; } @antisens; # synonyms of 
		map {s/,cephalosporins,/,cephalosporin,/g; } @antisens; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/g; } @antisens; # synonyms of 
		map {s/,chloramphenical,/,chloramphenicol,/g; } @antisens; # synonyms of 
		map {s/,aureomycin,/,chlortetracycline,/g; } @antisens; # synonyms of 
		map {s/,co-trimoxazole,/,trimethoprim,sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,tmp\/smx,/,trimethoprim,sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,polymyxin e,/,colistin,/g; } @antisens; # synonyms of 
		map {s/,gentian violet,/,crystal violet,/g; } @antisens; # synonyms of 
		map {s/,methyl violet 10b,/,crystal violet,/g; } @antisens; # synonyms of 
		map {s/,hexamethyl pararosaniline chloride,/,crystal violet,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.a,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.b,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.x,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.g,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c1,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c1a,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,gentam[iy]cin.c2,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,kanamycin a,/,kanamycin,/g; } @antisens; # synonyms of 
		map {s/,lomefloxacin hydrochloride ,/,lomefloxacin,/g; } @antisens; # synonyms of 
		map {s/,albamycin,/,novobiocin,/g; } @antisens; # synonyms of 
		map {s/,cathomycin,/,novobiocin,/g; } @antisens; # synonyms of 
		map {s/,oleandomycin,/,oleandomycin,/g; } @antisens; # synonyms of 
		map {s/,penicillin v,/,phenoxymethylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,procaine penicillin,/,procaine benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/g; } @antisens; # synonyms of 
		map {s/,polym[iy]xin,/,polymyxins,/g; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/g; } @antisens; # synonyms of 
		map {s/,polym[iy]xin.m,/,polymyxin m,/g; } @antisens; # synonyms of 
		map {s/,rifampin,/,rifampicin,/g; } @antisens; # synonyms of 
		map {s/,sps,/,sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,polyanetholsulfonic acid,/,sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,sodium polyanetholsulfonate,/,sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,sulfadiazin,/,sulfadiazine,/g; } @antisens; # synonyms of 
		map {s/,sulfamethoxazole-trimethoprim,/,sulfamethoxazole,trimethoprim,/g; } @antisens; # synonyms of 
		map {s/,sulfamethoxazole\/trimethoprim,/,sulfamethoxazole,trimethoprim,/g; } @antisens; # synonyms of 
		map {s/,smz,/,sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,smx,/,sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,sulfasomidine,/,sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,sulfamethin,/,sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,sulfaisodimidine,/,sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,tmp,/,trimethoprim,/g; } @antisens; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole,/,sulfamethoxazole,trimethoprim,/g; } @antisens; # synonyms of 
		map {s/,mepicycline penicillinate,/,penimepicycline,/g; } @antisens; # synonyms of 
		map {s/,flagecidin,/,anisomycin,/g; } @antisens; # synonyms of 
		map {s/,aueromycin,/,chlortetracycline,/g; } @antisens; # synonyms of 
		map {s/,cephamezine,/,cefazolin,/g; } @antisens; # synonyms of 
		map {s/,cephamezine vi,/,cefazolin,/g; } @antisens; # synonyms of 
		map {s/,gentamycin,/,gentamicin,/g; } @antisens; # synonyms of 
		map {s/,penecillin,/,penicillin,/g; } @antisens; # synonyms of 

		map {s/,not acetylspiramycin,/,not spiramycin,/g; } @antisens; # synonyms of 
		map {s/,not dactinomycin,/,not actinomycin,/g; } @antisens; # synonyms of 
		map {s/,not actinomycin d,/,not actinomycin,/g; } @antisens; # synonyms of 
		map {s/,not amoxycillin,/,not amoxicillin,/g; } @antisens; # synonyms of 
		map {s/,not benzylphenicillin,/,not benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,not penicillin g,/,not benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,not malachite green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not emerald green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not solid green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not diamond green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not aniline green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not benzaldehyde green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not fast green,/,not brilliant green,/g; } @antisens; # synonyms of 
		map {s/,not cefazoline,/,not cefazolin,/g; } @antisens; # synonyms of 
		map {s/,not cephazolin,/,not cefazolin,/g; } @antisens; # synonyms of 
		map {s/,not cefobid,/,not cefoperazone,/g; } @antisens; # synonyms of 
		map {s/,not mefoxin,/,not cefoxitin,/g; } @antisens; # synonyms of 
		map {s/,not cephaloridin,/,not cephaloridine,/g; } @antisens; # synonyms of 
		map {s/,not cefaloridine,/,not cephaloridine,/g; } @antisens; # synonyms of 
		map {s/,not cephalosporins,/,not cephalosporin,/g; } @antisens; # synonyms of 
		map {s/,not cephalothin,/,not cefalotin,/g; } @antisens; # synonyms of 
		map {s/,not chloramphenical,/,not chloramphenicol,/g; } @antisens; # synonyms of 
		map {s/,not aureomycin,/,not chlortetracycline,/g; } @antisens; # synonyms of 
		map {s/,not co-trimoxazole,/,not trimethoprim,sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,not polymyxin e,/,not colistin,/g; } @antisens; # synonyms of 
		map {s/,not gentian violet,/,not crystal violet,/g; } @antisens; # synonyms of 
		map {s/,not methyl violet 10b,/,not crystal violet,/g; } @antisens; # synonyms of 
		map {s/,not hexamethyl pararosaniline chloride,/,not crystal violet,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.a,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.b,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.x,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.g,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c1,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c1a,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not gentam[iy]cin.c2,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not kanamycin a,/,not kanamycin,/g; } @antisens; # synonyms of 
		map {s/,not lomefloxacin hydrochloride ,/,not lomefloxacin,/g; } @antisens; # synonyms of 
		map {s/,not albamycin,/,not novobiocin,/g; } @antisens; # synonyms of 
		map {s/,not cathomycin,/,not novobiocin,/g; } @antisens; # synonyms of 
		map {s/,not oleandomycin,/,not oleandomycin,/g; } @antisens; # synonyms of 
		map {s/,not penicillin v,/,not phenoxymethylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,not procaine penicillin,/,not procaine benzylpenicillin,/g; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/g; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin,/,not polymyxins,/g; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/g; } @antisens; # synonyms of 
		map {s/,not polym[iy]xin.m,/,not polymyxin m,/g; } @antisens; # synonyms of 
		map {s/,not rifampin,/,not rifampicin,/g; } @antisens; # synonyms of 
		map {s/,not sps,/,not sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,not polyanetholsulfonic acid,/,not sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,not sodium polyanetholsulfonate,/,not sodium polyanethol sulfonate,/g; } @antisens; # synonyms of 
		map {s/,not sulfadiazin,/,not sulfadiazine,/g; } @antisens; # synonyms of 
		map {s/,not smz,/,not sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,not smx,/,not sulfamethoxazole,/g; } @antisens; # synonyms of 
		map {s/,not sulfasomidine,/,not sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,not sulfamethin,/,not sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,not sulfaisodimidine,/,not sulfisomidine,/g; } @antisens; # synonyms of 
		map {s/,not tmp,/,not trimethoprim,/g; } @antisens; # synonyms of 
		map {s/,not mepicycline penicillinate,/,not penimepicycline,/g; } @antisens; # synonyms of 
		map {s/,not flagecidin,/,not anisomycin,/g; } @antisens; # synonyms of 
		map {s/,not aueromycin,/,not chlortetracycline,/g; } @antisens; # synonyms of 
		map {s/,not cephamezine,/,not cefazolin,/g; } @antisens; # synonyms of 
		map {s/,not cephamezine vi,/,not cefazolin,/g; } @antisens; # synonyms of 
		map {s/,not gentamycin,/,not gentamicin,/g; } @antisens; # synonyms of 
		map {s/,not penecillin,/,not penicillin,/g; } @antisens; # synonyms of 

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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Antibiotic Sensitivity Matrix'\;\n\tDIMENSIONS NCHAR=106\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /,organochlorine antibiotics,|,clindamycin,|,lincomycin,|,chloramphenicol,/) {
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
		#ontology terms and synonyms
		map {s/,acetylspiramycin,/,spiramycin,/g; } @antiresis; # synonyms of 
		map {s/,dactinomycin,/,actinomycin,/g; } @antiresis; # synonyms of 
		map {s/,actinomycin d,/,actinomycin,/g; } @antiresis; # synonyms of 
		map {s/,amoxycillin,/,amoxicillin,/g; } @antiresis; # synonyms of 
		map {s/,amoxicillin\/clavulanic acid,/,amoxicillin,clavulanic acid,/g; } @antiresis; # synonyms of 
		map {s/,benzylphenicillin,/,benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,penicillin g,/,benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,malachite green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,emerald green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,solid green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,diamond green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,aniline green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,benzaldehyde green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,fast green,/,brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,cefazoline,/,cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,cephazolin,/,cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,cefobid,/,cefoperazone,/g; } @antiresis; # synonyms of 
		map {s/,mefoxin,/,cefoxitin,/g; } @antiresis; # synonyms of 
		map {s/,cephaloridin,/,cephaloridine,/g; } @antiresis; # synonyms of 
		map {s/,cefaloridine,/,cephaloridine,/g; } @antiresis; # synonyms of 
		map {s/,cephalosporins,/,cephalosporin,/g; } @antiresis; # synonyms of 
		map {s/,cephalothin,/,cefalotin,/g; } @antiresis; # synonyms of 
		map {s/,chloramphenical,/,chloramphenicol,/g; } @antiresis; # synonyms of 
		map {s/,aureomycin,/,chlortetracycline,/g; } @antiresis; # synonyms of 
		map {s/,co-trimoxazole,/,trimethoprim,sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,tmp\/smx,/,trimethoprim,sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,polymyxin e,/,colistin,/g; } @antiresis; # synonyms of 
		map {s/,gentian violet,/,crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,methyl violet 10b,/,crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,hexamethyl pararosaniline chloride,/,crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.a,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.b,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.x,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.g,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c1,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c2,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,gentam[iy]cin.c1a,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,kanamycin a,/,kanamycin,/g; } @antiresis; # synonyms of 
		map {s/,lomefloxacin hydrochloride ,/,lomefloxacin,/g; } @antiresis; # synonyms of 
		map {s/,albamycin,/,novobiocin,/g; } @antiresis; # synonyms of 
		map {s/,cathomycin,/,novobiocin,/g; } @antiresis; # synonyms of 
		map {s/,oleandomycin,/,oleandomycin,/g; } @antiresis; # synonyms of 
		map {s/,penicillin v,/,phenoxymethylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,procaine penicillin,/,procaine benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/g; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin,/,polymyxins,/g; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.b,/,polymyxin b,/g; } @antiresis; # synonyms of 
		map {s/,polym[iy]xin.m,/,polymyxin m,/g; } @antiresis; # synonyms of 
		map {s/,rifampin,/,rifampicin,/g; } @antiresis; # synonyms of 
		map {s/,sps,/,sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,polyanetholsulfonic acid,/,sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,sodium polyanetholsulfonate,/,sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,sulfadiazin,/,sulfadiazine,/g; } @antiresis; # synonyms of 
		map {s/,sulfamethoxazole-trimethoprim,/,sulfamethoxazole,trimethoprim,/g; } @antiresis; # synonyms of 
		map {s/,sulfamethoxazole\/trimethoprim,/,sulfamethoxazole,trimethoprim,/g; } @antiresis; # synonyms of 
		map {s/,smz,/,sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,smx,/,sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,sulfasomidine,/,sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,sulfamethin,/,sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,sulfaisodimidine,/,sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,tmp,/,trimethoprim,/g; } @antiresis; # synonyms of 
		map {s/,trimethoprim.sulfamethoxazole,/,sulfamethoxazole,trimethoprim,/g; } @antiresis; # synonyms of 
		map {s/,mepicycline penicillinate,/,penimepicycline,/g; } @antiresis; # synonyms of 
		map {s/,flagecidin,/,anisomycin,/g; } @antiresis; # synonyms of 
		map {s/,aueromycin,/,chlortetracycline,/g; } @antiresis; # synonyms of 
		map {s/,cephamezine,/,cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,cephamezine vi,/,cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,gentamycin,/,gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,penecillin,/,penicillin,/g; } @antiresis; # synonyms of 

		map {s/,not acetylspiramycin,/,not spiramycin,/g; } @antiresis; # synonyms of 
		map {s/,not dactinomycin,/,not actinomycin,/g; } @antiresis; # synonyms of 
		map {s/,not actinomycin d,/,not actinomycin,/g; } @antiresis; # synonyms of 
		map {s/,not amoxycillin,/,not amoxicillin,/g; } @antiresis; # synonyms of 
		map {s/,not benzylphenicillin,/,not benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,not penicillin g,/,not benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,not malachite green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not emerald green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not solid green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not diamond green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not aniline green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not benzaldehyde green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not fast green,/,not brilliant green,/g; } @antiresis; # synonyms of 
		map {s/,not cefazoline,/,not cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,not cephazolin,/,not cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,not cefobid,/,not cefoperazone,/g; } @antiresis; # synonyms of 
		map {s/,not mefoxin,/,not cefoxitin,/g; } @antiresis; # synonyms of 
		map {s/,not cephaloridin,/,not cephaloridine,/g; } @antiresis; # synonyms of 
		map {s/,not cefaloridine,/,not cephaloridine,/g; } @antiresis; # synonyms of 
		map {s/,not cephalosporins,/,not cephalosporin,/g; } @antiresis; # synonyms of 
		map {s/,not cephalothin,/,not cefalotin,/g; } @antiresis; # synonyms of 
		map {s/,not chloramphenical,/,not chloramphenicol,/g; } @antiresis; # synonyms of 
		map {s/,not aureomycin,/,not chlortetracycline,/g; } @antiresis; # synonyms of 
		map {s/,not co-trimoxazole,/,not trimethoprim,sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,not polymyxin e,/,not colistin,/g; } @antiresis; # synonyms of 
		map {s/,not gentian violet,/,not crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,not methyl violet 10b,/,not crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,not hexamethyl pararosaniline chloride,/,not crystal violet,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.a,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.b,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.x,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.g,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c1,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c1a,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not gentam[iy]cin.c2,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not kanamycin a,/,not kanamycin,/g; } @antiresis; # synonyms of 
		map {s/,not lomefloxacin hydrochloride ,/,not lomefloxacin,/g; } @antiresis; # synonyms of 
		map {s/,not albamycin,/,not novobiocin,/g; } @antiresis; # synonyms of 
		map {s/,not cathomycin,/,not novobiocin,/g; } @antiresis; # synonyms of 
		map {s/,not oleandomycin,/,not oleandomycin,/g; } @antiresis; # synonyms of 
		map {s/,not penicillin v,/,not phenoxymethylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,not procaine penicillin,/,not procaine benzylpenicillin,/g; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/g; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin,/,not polymyxins,/g; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.b,/,not polymyxin b,/g; } @antiresis; # synonyms of 
		map {s/,not polym[iy]xin.m,/,not polymyxin m,/g; } @antiresis; # synonyms of 
		map {s/,not rifampin,/,not rifampicin,/g; } @antiresis; # synonyms of 
		map {s/,not sps,/,not sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,not polyanetholsulfonic acid,/,not sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,not sodium polyanetholsulfonate,/,not sodium polyanethol sulfonate,/g; } @antiresis; # synonyms of 
		map {s/,not sulfadiazin,/,not sulfadiazine,/g; } @antiresis; # synonyms of 
		map {s/,not smz,/,not sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,not smx,/,not sulfamethoxazole,/g; } @antiresis; # synonyms of 
		map {s/,not sulfasomidine,/,not sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,not sulfamethin,/,not sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,not sulfaisodimidine,/,not sulfisomidine,/g; } @antiresis; # synonyms of 
		map {s/,not tmp,/,not trimethoprim,/g; } @antiresis; # synonyms of 
		map {s/,not mepicycline penicillinate,/,not penimepicycline,/g; } @antiresis; # synonyms of 
		map {s/,not flagecidin,/,not anisomycin,/g; } @antiresis; # synonyms of 
		map {s/,not aueromycin,/,not chlortetracycline,/g; } @antiresis; # synonyms of 
		map {s/,not cephamezine,/,not cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,not cephamezine vi,/,not cefazolin,/g; } @antiresis; # synonyms of 
		map {s/,not gentamycin,/,not gentamicin,/g; } @antiresis; # synonyms of 
		map {s/,not penecillin,/,not penicillin,/g; } @antiresis; # synonyms of 

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
		map {s/,antibiotic resistant,//g; } @homcharlist; # gets rid of the character name label
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE' Antibiotic Resistance Matrix'\;\n\tDIMENSIONS NCHAR=106\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /,organochlorine antibiotics,|,clindamycin,|,lincomycin,|,chloramphenicol,/) {
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
		map {s/,colony is round,/,circular,/g; } @colshape; # synonyms of 
		map {s/,round,/,circular,/g; } @colshape; # synonyms of 
		map {s/,circular colony,/,circular,/g; } @colshape; # synonyms of 
		map {s/,circular droplets,/,circular,/g; } @colshape; # synonyms of 
		map {s/,entire colony,/,entire,/g; } @colshape; # synonyms of 
		map {s/,convex elevation,/,convex,/g; } @colshape; # synonyms of 
		map {s/,domed,/,convex,/g; } @colshape; # synonyms of 
		map {s/,convex colony,/,convex,/g; } @colshape; # synonyms of 
		map {s/,collapsed center,/,crateriform,/g; } @colshape; # synonyms of 
		map {s/,crater-shaped,/,crateriform,/g; } @colshape; # synonyms of 
		map {s/,cup-shaped,/,crateriform,/g; } @colshape; # synonyms of 
		map {s/,dimpled,/,crateriform,/g; } @colshape; # synonyms of 
		map {s/,crateriform colony,/,crateriform,/g; } @colshape; # synonyms of 
		map {s/,wavy margin,/,curled,/g; } @colshape; # synonyms of 
		map {s/,curled colony,/,curled,/g; } @colshape; # synonyms of 
		map {s/,veined,/,dendritic,/g; } @colshape; # synonyms of 
		map {s/,dendritic colony,/,dendritic,/g; } @colshape; # synonyms of 
		map {s/,erose colony,/,erose,/g; } @colshape; # synonyms of 
		map {s/,notched,/,erose,/g; } @colshape; # synonyms of 
		map {s/,toothed,/,erose,/g; } @colshape; # synonyms of 
		map {s/,fuzzy,/,filamentous,/g; } @colshape; # synonyms of 
		map {s/,flat elevation,/,flat,/g; } @colshape; # synonyms of 
		map {s/,flat colony,/,flat,/g; } @colshape; # synonyms of 
		map {s/,lobed margin,/,lobed,/g; } @colshape; # synonyms of 
		map {s/,pulvinate elevation,/,pulvinate,/g; } @colshape; # synonyms of 
		map {s/,pulvinate colony,/,pulvinate,/g; } @colshape; # synonyms of 
		map {s/,punctiform colony,/,punctiform,/g; } @colshape; # synonyms of 
		map {s/,punctiforme,/,punctiform,/g; } @colshape; # synonyms of 
		map {s/,raised colony,/,raised,/g; } @colshape; # synonyms of 
		map {s/,fusiform,/,spindle,/g; } @colshape; # synonyms of 
		map {s/,spindle-shaped,/,spindle,/g; } @colshape; # synonyms of 
		map {s/,undulate margin,/,undulate,/g; } @colshape; # synonyms of 
		map {s/,satellite behavior,/,satellite,/g; } @colshape; # synonyms of 
		map {s/,satellites,/,satellite,/g; } @colshape; # synonyms of 
		map {s/,satellite colonies,/,satellite,/g; } @colshape; # synonyms of 
		map {s/,fine granulation,/,granular,/g; } @colshape; # synonyms of 
		map {s/,granular center,/,granular,/g; } @colshape; # synonyms of 
		map {s/,mucous,/,mucoid,/g; } @colshape; # synonyms of 
		map {s/,disc.shaped,/,flat,/g; } @colshape; # synonyms of 
		map {s/,disk.shaped,/,flat,/g; } @colshape; # synonyms of 
		map {s/,elevated,/,raised,/g; } @colshape; # synonyms of 
		map {s/,flat colonies,/,flat,/g; } @colshape; # synonyms of 
		map {s/,hemispherical knobs,/,convex,/g; } @colshape; # synonyms of 
		map {s/,high convex,/,raised,convex,/g; } @colshape; # synonyms of 
		map {s/,irregular,/,asymmetrical,/g; } @colshape; # synonyms of 
		map {s/,low convex,/,flat,convex,/g; } @colshape; # synonyms of 
		map {s/,raised at the center,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,raised at the centres,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,raised in the center,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,umbonate elevation,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,rhizoid form,/,rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,raised rhizoid growth,/,rhizoidal,raised,/g; } @colshape; # synonyms of 
		map {s/,rhizoid form,/,rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,rhizoid,/,rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,more rhizoid,/,rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,peak,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,peaked,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,pin point,/,punctiform,/g; } @colshape; # synonyms of 
		map {s/,punctate,/,punctiform,/g; } @colshape; # synonyms of 
		map {s/,raised,/,raised,/g; } @colshape; # synonyms of 
		map {s/,raised in the centre,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,rounded,/,circular,/g; } @colshape; # synonyms of 
		map {s/,slightly convex,/,convex,/g; } @colshape; # synonyms of 
		map {s/,slightly irregular,/,asymmetrical,/g; } @colshape; # synonyms of 
		map {s/,slightly peaked,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,slightly raised,/,raised,/g; } @colshape; # synonyms of 
		map {s/,slightly sunken,/,sunken,/g; } @colshape; # synonyms of 
		map {s/,slightly umbonate,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,slightly umobonate,/,umbonate,/g; } @colshape; # synonyms of 
		map {s/,sometimes feebly rhizoid,/,rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,spreading,/,effuse,/g; } @colshape; # synonyms of 
		map {s/,thin,/,diffuse,/g; } @colshape; # synonyms of 
		map {s/,friable,/,dry,/g; } @colshape; # synonyms of 
		map {s/,dry colony,/,dry,/g; } @colshape; # synonyms of 
		map {s/,fried.egg,/,fried-egg,/g; } @colshape; # synonyms of 
		map {s/,fried.eggs,/,fried-egg,/g; } @colshape; # synonyms of 

		map {s/,not friable,/,not dry,/g; } @colshape; # synonyms of 
		map {s/,not dry colony,/,not dry,/g; } @colshape; # synonyms of 
		map {s/,not fried.egg,/,not fried-egg,/g; } @colshape; # synonyms of 
		map {s/,not fried.eggs,/,not fried-egg,/g; } @colshape; # synonyms of 
		map {s/,not colony is round,/,not circular,/g; } @colshape; # synonyms of 
		map {s/,not round,/,not circular,/g; } @colshape; # synonyms of 
		map {s/,not circular colony,/,not circular,/g; } @colshape; # synonyms of 
		map {s/,not circular droplets,/,not circular,/g; } @colshape; # synonyms of 
		map {s/,not entire colony,/,not entire,/g; } @colshape; # synonyms of 
		map {s/,not convex elevation,/,not convex,/g; } @colshape; # synonyms of 
		map {s/,not domed,/,not convex,/g; } @colshape; # synonyms of 
		map {s/,not convex colony,/,not convex,/g; } @colshape; # synonyms of 
		map {s/,not collapsed center,/,not crateriform,/g; } @colshape; # synonyms of 
		map {s/,not crater-shaped,/,not crateriform,/g; } @colshape; # synonyms of 
		map {s/,not cup-shaped,/,not crateriform,/g; } @colshape; # synonyms of 
		map {s/,not dimpled,/,not crateriform,/g; } @colshape; # synonyms of 
		map {s/,not crateriform colony,/,not crateriform,/g; } @colshape; # synonyms of 
		map {s/,not wavy margin,/,not curled,/g; } @colshape; # synonyms of 
		map {s/,not curled colony,/,not curled,/g; } @colshape; # synonyms of 
		map {s/,not veined,/,not dendritic,/g; } @colshape; # synonyms of 
		map {s/,not dendritic colony,/,not dendritic,/g; } @colshape; # synonyms of 
		map {s/,not erose colony,/,not erose,/g; } @colshape; # synonyms of 
		map {s/,not notched,/,not erose,/g; } @colshape; # synonyms of 
		map {s/,not toothed,/,not erose,/g; } @colshape; # synonyms of 
		map {s/,not fuzzy,/,not filamentous,/g; } @colshape; # synonyms of 
		map {s/,not flat elevation,/,not flat,/g; } @colshape; # synonyms of 
		map {s/,not flat colony,/,not flat,/g; } @colshape; # synonyms of 
		map {s/,not lobed margin,/,not lobed,/g; } @colshape; # synonyms of 
		map {s/,not pulvinate elevation,/,not pulvinate,/g; } @colshape; # synonyms of 
		map {s/,not pulvinate colony,/,not pulvinate,/g; } @colshape; # synonyms of 
		map {s/,not punctiform colony,/,not punctiform,/g; } @colshape; # synonyms of 
		map {s/,not punctiforme,/,not punctiform,/g; } @colshape; # synonyms of 
		map {s/,not raised colony,/,not raised,/g; } @colshape; # synonyms of 
		map {s/,not fusiform,/,not spindle,/g; } @colshape; # synonyms of 
		map {s/,not spindle-shaped,/,not spindle,/g; } @colshape; # synonyms of 
		map {s/,not undulate margin,/,not undulate,/g; } @colshape; # synonyms of 
		map {s/,not satellite behavior,/,not satellite,/g; } @colshape; # synonyms of 
		map {s/,not satellites,/,not satellite,/g; } @colshape; # synonyms of 
		map {s/,not satellite colonies,/,not satellite,/g; } @colshape; # synonyms of 
		map {s/,not fine granulation,/,not granular,/g; } @colshape; # synonyms of 
		map {s/,not granular center,/,not granular,/g; } @colshape; # synonyms of 
		map {s/,not mucous,/,not mucoid,/g; } @colshape; # synonyms of 
		map {s/,not disc.shaped,/,not flat,/g; } @colshape; # synonyms of 
		map {s/,not disk.shaped,/,not flat,/g; } @colshape; # synonyms of 
		map {s/,not elevated,/,not raised,/g; } @colshape; # synonyms of 
		map {s/,not flat colonies,/,not flat,/g; } @colshape; # synonyms of 
		map {s/,not hemispherical knobs,/,not convex,/g; } @colshape; # synonyms of 
		map {s/,not irregular,/,not asymmetrical,/g; } @colshape; # synonyms of 
		map {s/,not low convex,/,not flat,convex,/g; } @colshape; # synonyms of 
		map {s/,not raised at the center,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not raised at the centres,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not raised in the center,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not umbonate elevation,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not rhizoid form,/,not rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,not rhizoid form,/,not rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,not rhizoid,/,not rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,not more rhizoid,/,not rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,not peak,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not peaked,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not pin point,/,not punctiform,/g; } @colshape; # synonyms of 
		map {s/,not punctate,/,not punctiform,/g; } @colshape; # synonyms of 
		map {s/,not raised,/,not raised,/g; } @colshape; # synonyms of 
		map {s/,not raised in the centre,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not rounded,/,not circular,/g; } @colshape; # synonyms of 
		map {s/,not slightly convex,/,not convex,/g; } @colshape; # synonyms of 
		map {s/,not slightly irregular,/,not asymmetrical,/g; } @colshape; # synonyms of 
		map {s/,not slightly peaked,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not slightly raised,/,not raised,/g; } @colshape; # synonyms of 
		map {s/,not slightly sunken,/,not sunken,/g; } @colshape; # synonyms of 
		map {s/,not slightly umbonate,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not slightly umobonate,/,not umbonate,/g; } @colshape; # synonyms of 
		map {s/,not sometimes feebly rhizoid,/,not rhizoidal,/g; } @colshape; # synonyms of 
		map {s/,not spreading,/,not effuse,/g; } @colshape; # synonyms of 
		map {s/,not thin,/,not diffuse,/g; } @colshape; # synonyms of 

#not colony shape
		map {s/,numerous filamentous tufts,/,/g; } @colshape; # synonyms of 

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
		if ($line =~ /Tree Name.*/) {
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
		if ($line =~ /Tree Name.*/) {
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
		map {s/,butryaceous,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,buttery texture,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,butyrous,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,butyrous consistency,/,butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,do not adhere to the agar,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,non-adherent,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,tough,/,adherent texture,/g; } @coltext; # synonyms of 
		map {s/,flaky,/,flaky texture,/g; } @coltext; # synonyms of 
		map {s/,hard to emulsify,/,not easily emulsified,/g; } @coltext; # synonyms of 
		map {s/,glistening,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,shiny,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,slimy,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,slime,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,mucoid,/,slimy texture,/g; } @coltext; # synonyms of 
		map {s/,polished,/,glistening texture,/g; } @coltext; # synonyms of 
		map {s/,opalescence,/,opalescent texture,/g; } @coltext; # synonyms of 
		map {s/,rough matt surface,/,dull texture,/g; } @coltext; # synonyms of 
		map {s/,slightly aromatic odour,/,aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,aromatic,/,aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,small shallow pits,/,rough texture,/g; } @coltext; # synonyms of 
		map {s/,smooth,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth surface,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth surfaces,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,smooth-surfaced,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,soft,/,smooth texture,/g; } @coltext; # synonyms of 
		map {s/,stringy,/,ropy texture,/g; } @coltext; # synonyms of 
		map {s/,sunken,/,sunken texture,/g; } @coltext; # synonyms of 
		map {s/,growth tends to be mainly below the surface,/,sunken texture,/g; } @coltext; # synonyms of 
		map {s/,translucent,/,translucent colony,/g; } @coltext; # synonyms of 
		map {s/,translucent colonies,/,translucent colony,/g; } @coltext; # synonyms of 
		map {s/,elastic-gummy consistency,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,viscid,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,viscid consistency,/,viscid colony,/g; } @coltext; # synonyms of 
		map {s/,watery,/,watery colony,/g; } @coltext; # synonyms of 
		map {s/,rough surface,/,rough texture,/g; } @coltext; # synonyms of 

		map {s/,not rough surface,/,not rough texture,/g; } @coltext; # synonyms of 
		map {s/,not butryaceous,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not buttery texture,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not butyrous,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not butyrous consistency,/,not butyrous texture,/g; } @coltext; # synonyms of 
		map {s/,not do not adhere to the agar,/,not not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not non-adherent,/,not not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not tough,/,not adherent texture,/g; } @coltext; # synonyms of 
		map {s/,not flaky,/,not flaky texture,/g; } @coltext; # synonyms of 
		map {s/,easy to emulsify,/,easily emulsified,/g; } @coltext; # synonyms of 
		map {s/,not glistening,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not shiny,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not slimy,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not slime,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not mucoid,/,not slimy texture,/g; } @coltext; # synonyms of 
		map {s/,not polished,/,not glistening texture,/g; } @coltext; # synonyms of 
		map {s/,not opalescence,/,not opalescent texture,/g; } @coltext; # synonyms of 
		map {s/,not rough matt surface,/,not dull texture,/g; } @coltext; # synonyms of 
		map {s/,not slightly aromatic odour,/,not aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,not aromatic,/,not aromatic colony,/g; } @coltext; # synonyms of 
		map {s/,not small shallow pits,/,not rough texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth surface,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth surfaces,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not smooth-surfaced,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not soft,/,not smooth texture,/g; } @coltext; # synonyms of 
		map {s/,not stringy,/,not ropy texture,/g; } @coltext; # synonyms of 
		map {s/,not sunken,/,not sunken texture,/g; } @coltext; # synonyms of 
		map {s/,not growth tends to be mainly below the surface,/,not sunken texture,/g; } @coltext; # synonyms of 
		map {s/,not translucent,/,not translucent colony,/g; } @coltext; # synonyms of 
		map {s/,not translucent colonies,/,not translucent colony,/g; } @coltext; # synonyms of 
		map {s/,not elastic-gummy consistency,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not viscid,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not viscid consistency,/,not viscid colony,/g; } @coltext; # synonyms of 
		map {s/,not watery,/,not watery colony,/g; } @coltext; # synonyms of 

#not colony texture
		map {s/,compact centre,/,/g; } @coltext; # synonyms of 
		map {s/,satellite colonies,/,/g; } @coltext; # synonyms of 
		map {s/,yellow colony,/,/g; } @coltext; # synonyms of 
		map {s/,not yellow colony,/,/g; } @coltext; # synonyms of 
		map {s/,slightly filamentous,/,/g; } @coltext; # synonyms of 
		map {s/,irregular,/,/g; } @coltext; # synonyms of 
		map {s/,numerous filamentous tufts,/,/g; } @coltext; # synonyms of 
		map {s/,bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,yellow,/,/g; } @coltext; # synonyms of 
		map {s/,bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not bright yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not dirty yellow,/,/g; } @coltext; # synonyms of 
		map {s/,not yellow,/,/g; } @coltext; # synonyms of 

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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Texture Matrix'\;\n\tDIMENSIONS NCHAR=13\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		map {s/,beige,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,slightly brown,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,water-soluble brown pigment,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,tan,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,brown,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,brown deposit,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,buff,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,dark brown,/,brown,/g; } @colcolor; # synonyms of 
		map {s/,reddish-brown diffusible pigment,/,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,reddish brown,/,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,deep orange-brown,/,dark,reddish-brown,/g; } @colcolor; # synonyms of 
		map {s/,bright pink,/,dark,red,/g; } @colcolor; # synonyms of 
		map {s/,dark pink,/,dark,red,/g; } @colcolor; # synonyms of 
		map {s/,light pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,light red,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,light-pink-coloured,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pale pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pale pinkish gray,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pale speckled-pink,/,pale,mottled,red,/g; } @colcolor; # synonyms of 
		map {s/,pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pink-coloured,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,pinkish,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,reddish,/,red,/g; } @colcolor; # synonyms of 
		map {s/,reddish-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,rose-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,rose-colored,/,red,/g; } @colcolor; # synonyms of 
		map {s/,rusty-orange-coloured,/,red,/g; } @colcolor; # synonyms of 
		map {s/,saffron-colored,/,red,/g; } @colcolor; # synonyms of 
		map {s/,slightly pink,/,pale,red,/g; } @colcolor; # synonyms of 
		map {s/,dark reddish orange,/,dark,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,light salmon pink,/,pale,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,moderately reddish-orange,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,orangish red,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,reddish.orange,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,salmon,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,salmon-pink,/,orange-red,/g; } @colcolor; # synonyms of 
		map {s/,bright orange,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,dark orange,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,dark-orange-coloured,/,dark,orange,/g; } @colcolor; # synonyms of 
		map {s/,light orange,/,pale,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange coloured,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange fluorescence,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange pigment,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,orange-coloured,/,orange,/g; } @colcolor; # synonyms of 
		map {s/,pale orange,/,pale,orange,/g; } @colcolor; # synonyms of 
		map {s/,yellowish orange,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,yellow.orange,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,pale orange-yellow,/,pale,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,orange\/yellow,/,yellow-orange,/g; } @colcolor; # synonyms of 
		map {s/,bright golden yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,bright yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,bright yellow pigment,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,dark yellow,/,dark,yellow,/g; } @colcolor; # synonyms of 
		map {s/,dirty yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,golden-yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,light yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,medium yellow,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,pale.yellow,/,pale,yellow,/g; } @colcolor; # synonyms of 
		map {s/,translucent yellow,/,translucent,yellow,/g; } @colcolor; # synonyms of 
		map {s/,translucent yellowish,/,translucent,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow colonies,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow flexibacteria,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellow-pigmented,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,yellowish,/,yellow,/g; } @colcolor; # synonyms of 
		map {s/,purple,/,violet,/g; } @colcolor; # synonyms of 
		map {s/,light purple,/,pale,violet,/g; } @colcolor; # synonyms of 
		map {s/,black pigment,/,black,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]y,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]yish,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,greenish gr[ae]y,/,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]y-white,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,gr[ae]yish-white,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,shiny gr[ae]y periphery,/,shiny,gray,/g; } @colcolor; # synonyms of 
		map {s/,white-gr[ae]yish,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,whitegr[ae]yish,/,pale,gray,/g; } @colcolor; # synonyms of 
		map {s/,cream,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy,/,white,/g; } @colcolor; # synonyms of 
		map {s/,creamy white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,dorothy white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,off-white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,off-white centre,/,white,/g; } @colcolor; # synonyms of 
		map {s/,off-white-grey,/,white,/g; } @colcolor; # synonyms of 
		map {s/,pigmented creamy white colonies,/,white,/g; } @colcolor; # synonyms of 
		map {s/,smooth white,/,white,/g; } @colcolor; # synonyms of 
		map {s/,translucent.whitish,/,translucent,white,/g; } @colcolor; # synonyms of 
		map {s/,white fluorescent tubes,/,fluorescent,white,/g; } @colcolor; # synonyms of 
		map {s/,whitish,/,white,/g; } @colcolor; # synonyms of 
		map {s/,matt,/,dull,/g; } @colcolor; # synonyms of 
		map {s/,shiny,/,glistening,/g; } @colcolor; # synonyms of 
		map {s/,shiny gr[ae]y periphery,/,gray,glistening,/g; } @colcolor; # synonyms of 
		map {s/,glitstening,/,glistening,/g; } @colcolor; # synonyms of 
		map {s/,trasnparent,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,semi-transparent,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,clear,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,clear margins,/,transparent,/g; } @colcolor; # synonyms of 
		map {s/,colorless,/,white,/g; } @colcolor; # synonyms of 
		map {s/,entire translucent margins,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,semitranslucent,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,translucent colonies,/,translucent,/g; } @colcolor; # synonyms of 
		map {s/,semiopaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,semi-opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,slightly opaque,/,opaque,/g; } @colcolor; # synonyms of 
		map {s/,slightly motlted,/,mottled,/g; } @colcolor; # synonyms of 
		map {s/,pronounced metallic tinge,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,metallic sheen,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,iridescence,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,iridescent lustre,/,iridescent,/g; } @colcolor; # synonyms of 
		map {s/,light,/,pale,/g; } @colcolor; # synonyms of 
		map {s/,slightly,/,pale,/g; } @colcolor; # synonyms of 
		map {s/,deep,/,dark,/g; } @colcolor; # synonyms of 
		map {s/,bright,/,dark,/g; } @colcolor; # synonyms of 
		map {s/,much pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,pigmentation,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,non-watersoluble pigment,/,pigmented,/g; } @colcolor; # synonyms of 
		map {s/,light brown,/,pale,brown,/g; } @colcolor; # synonyms of 
		map {s/,slightly mottled,/,pale,mottled,/g; } @colcolor; # synonyms of 

		map {s/,not white fluorescent tubes,/,not fluorescent,not white,/g; } @colcolor; # synonyms of 
		map {s/,not pigment,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,not bright yellow pigment,/,not yellow,/g; } @colcolor; # synonyms of 
		map {s/,not black pigment,/,not black,/g; } @colcolor; # synonyms of 
		map {s/,without diffusible pigment,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,non-pigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,nonpigmented,/,not pigmented,/g; } @colcolor; # synonyms of 
		map {s/,non-translucent,/,not translucent,/g; } @colcolor; # synonyms of 

#not colony color
		map {s/,glistening,/,/g; } @colcolor; # synonyms of 
		map {s/,not glistening,/,/g; } @colcolor; # synonyms of 
		map {s/,blood agar,/,/g; } @colcolor; # synonyms of 
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
	map {s/^black/**black/g; } @charlist2;
	map {s/^brown/**brown/g; } @charlist2;
	map {s/^dark/**dark/g; } @charlist2;
	map {s/^dull/**dull/g; } @charlist2;
	map {s/^fluorescent/**fluorescent/g; } @charlist2;
	map {s/^gray/**gray/g; } @charlist2;
	map {s/^iridescent/**iridescent/g; } @charlist2;
	map {s/^mottled/**mottled/g; } @charlist2;
	map {s/^opaque/**opaque/g; } @charlist2;
	map {s/^orange-red/**orange-red/g; } @charlist2;
	map {s/^orange/**orange/g; } @charlist2;
	map {s/^pale/**pale/g; } @charlist2;
	map {s/^pigmented/**pigmented/g; } @charlist2;
	map {s/^reddish-brown/**reddish-brown/g; } @charlist2;
	map {s/^red/**red/g; } @charlist2;
	map {s/^translucent/**translucent/g; } @charlist2;
	map {s/^transparent/**transparent/g; } @charlist2;
	map {s/^violet/**violet/g; } @charlist2;
	map {s/^white/**white/g; } @charlist2;
	map {s/^yellow-orange/**yellow-orange/g; } @charlist2;
	map {s/^yellow/**yellow/g; } @charlist2;
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
	print OUT "\;\nEND\;\n\n\nBEGIN CHARACTERS\;\n\tTITLE 'Colony Color Matrix'\;\n\tDIMENSIONS NCHAR=18\;\n\tFORMAT DATATYPE \= STANDARD INTERLEAVE GAP \= \- MISSING \= \?\;\n\t";
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

	print OUT " \;\n\tMATRIX\n";
	open (IN, '<', $temp6) or die $!;
	open (OUT, '>>', $nexusoutfile) or die $!;
	while ($line = <IN>) {
		chomp $line;
		if ($line =~ /Tree Name.*/) {
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
		elsif ($line =~ /,not dull,|,glistening,/) {
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
		elsif ($line =~ /,pigmented,|,brown,|,reddish-brown,|,red,|,orange-red,|,orange,|,yellow-orange,|,yellow,|,violet,|,black,|,gray,/) {
			print OUT "1";
		}
		else {
			print OUT "?";
		}
#code char 18 pigmentation
		if ($line =~ /,not pigmented,/) {
			print OUT "0";
			}
		elsif ($line =~ /,pigmented,|,brown,|,reddish-brown,|,red,|,orange-red,|,orange,|,yellow-orange,|,yellow,|,violet,|,black,|,gray,/) {
			print OUT "1";
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






else {
	print "NOTE********************you picked a character that has been finished yet**********************\n";
	}
