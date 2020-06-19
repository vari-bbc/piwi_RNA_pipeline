#!/usr/bin/perl -w
use strict;
use List::Util qw(shuffle);

my @gtf = split(',',shift @ARGV);#
my $ngtf=scalar @gtf;
my @samples = split(',',shift @ARGV);#
#~ my @_samples = split(',',shift @ARGV);#
#~ my @samples = shuffle(@_samples);
my $nsamples=scalar @samples;
my $strand = shift @ARGV;#
my $indir = shift @ARGV;
my $fileext = shift @ARGV;
my $outfile = $indir . "assigned_counts." . $strand . ".txt";
die print "Error! strand: $strand\n" if $strand ne 'sense' && $strand ne 'antisense';



my $testing=0;

## Report some information
print STDOUT <<EOF;
Information:
	Number of samples: $nsamples
	Samples: @samples
	Number of GTFs: $ngtf
	GTF assignment: @gtf 
	Strand: $strand
	In directory: $indir
	Outfile: $outfile
	File extension: $fileext

EOF

#~ exit

my %reads=(); # if a read is initialized, then it is ignored
my %data=(); # sample -> feature -> count
my %summary=();
my %all_reads=(); # this is initialized for all reads so that we have a complete count even if not ever mapped
# The sum of all assigments for 1 read is 1

foreach my $s (@samples){
	print "\t$s\n";
	for my $g (@gtf){
		my $infile = $indir . $s . '.' . $strand . '.' . $g . '.' . $fileext;
		print "\t\t$g: $infile\n";
		my $n=0;

		# now open the infile
		#~ open(my $in,"awk \'\$2==\"Assigned\"\' $infile|");
		open(my $in,"<",$infile);
		while(<$in>){
			chomp;
			# print "\t\t\t$_\n";
			my($read,$assigned,$n_assignments,$feature_string)=split('\t',$_);
			if($assigned eq 'Assigned'){
				my @features = split(',',$feature_string);
				my $value = sprintf("%.2f",(1/$n_assignments));
				# print "\t\t\t$read -> @features\n";
				
				# If the read does not exist, then it is OK to be counted.
				# Distribute counts to the read's assigned feature(s )
				if (! exists $reads{$read}){
					foreach my $f (@features){
						#~ print "Assigning $read to $f with $value\n";
						$data{"$g\__$f"}{$s}+=$value; # put the gtf source in the file as well
						$n++; 
					}
					$reads{$read}=(); # initialize this reads as being assigned already
				}
				# even if it does exists, need to keep track
				$summary{$g}++ # this counts all assignments but fails to track where multiply assigned reads are being assigned
			}else{ # read NOT assigned
				$all_reads{$read}=();
			}
			
			if($testing){
				if($n>100 or $g eq 'mitochondrion_genome'){
					last;
				}
			} # testing
		}
		print "\t\t\tAssigned $n reads\n";
	}
}


open(my $out,'>',$outfile) || die print "$outfile not opened!!\n";
foreach my $s (@samples){ print $out "\t$s";}
print $out "\n";
# Now print the counts matrix
foreach my $feature (keys %data){ # random 
	print $out "$feature";
	foreach my $s (@samples){
		#~ if($testing){print STDOUT "$s\n"}
		if(exists $data{$feature}{$s}){
			print $out "\t" . sprintf("%.0f",$data{$feature}{$s});
		}else{
			#~ if($testing){print STDOUT "\t$feature NOT FOUND\n"}
			print $out "\t0";
		}
	}
	print $out "\n";
}

## print summary
print "="x100,"\n";
print "SUMMARY:\n";
foreach my $gtf (keys %summary){
	print "\t$gtf: ",$summary{$gtf},"\n";
}
