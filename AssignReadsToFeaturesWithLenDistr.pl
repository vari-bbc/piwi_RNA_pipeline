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
my $trimmeddir = shift @ARGV;
my $trimmed_suffix = shift @ARGV;
my $outfile = $indir . "length_distribution." . $strand . ".txt";
my $outfile_w_feature = $indir . "length_distribution." . $strand . ".withFeatures.txt";
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
	Outfile with features: $outfile
	File extension: $fileext

EOF

#~ exit;
my %data=(); # sample -> read -> gtf
my %data_wFeature=(); # sample -> read -> feature
my %fd=(); # final data
my %fd_wFeature=(); # final data
my %already_assigned=();

foreach my $s (@samples){
	print "\t$s\n";
	my $trimmed_fq = $trimmeddir . $s . $trimmed_suffix;
	#~ print "Trimmed: $trimmed_fq\n"; next;
	for my $g (@gtf){
		# next;
		my $infile = $indir . $s . '.' . $strand . '.' . $g . '.' . $fileext;
		print "\t\t$g: $infile\n";
		my $n=0;

		# now open the file that has the read assignments
		open(my $in,"<",$infile);
		while(<$in>){
			chomp;
			# print "\t\t\t$_\n";
			my($read,$assigned,$n_assignments,$feature_string)=split('\t',$_);
			
			
			# edit 5/4/2020 to exlude LSU-rRNA_Dme and SSU-rRNA_Dme -- DO NOT COUNT THESE
			if($feature_string eq 'LSU-rRNA_Dme'){
				$already_assigned{$read}=(); # initialize this read as being assigned already, so it does not get assigned again
				next;
			}
			if($feature_string eq 'SSU-rRNA_Dme'){
				$already_assigned{$read}=(); # initialize this read as being assigned already, so it does not get assigned again
				next;
			}
			
			
			# edit 4/28/2020
			# if $g is repeatmasker, divide this into simple and complex repeats depending on the feature_string
			my $gNew;
			if($g eq 'repeatmasker'){
				# get simple repeats
				if($feature_string =~/^\(/ || $feature_string =~/\-rich/){
					$gNew = 'simple_repeat';
				}else{
					$gNew = 'complex_repeat'
				}
			}else{
				$gNew = $g;
			}
			
			if($assigned eq 'Assigned'){
				# If the read does not exist, then it is OK to be counted for this gtf...
				if (! exists $already_assigned{$read}){
					$already_assigned{$read}=(); # initialize this read as being assigned already, so it does not get assigned again
					#~ print "Assigned $g feature as a $gNew\t\t$feature_string\n";
					$data{$read}=$gNew;
					$data_wFeature{$read}=$feature_string;
					$n++;
				}
			}else{ # read NOT assigned ... do nothing

			}

			if($testing){
				if($n>100 or $gNew eq 'mitochondrion_genome'){
					last;
				}
			} # testing
		}
		close($in);
		print "\t\t\tAssigned $n reads\n";
	}
	# Now you have all of the reads assigned, so it makes sense to go ahead and quantify
	# 1) first  nucleotide and 2) length of the reads assigned to each annotation from
	# the trimmed data

	# open the trimmed data
	my $in;
	if($trimmed_fq=~m/.gz$/){open($in,"zcat $trimmed_fq|");}
	if($trimmed_fq!~m/.gz$/){open($in,"cat $trimmed_fq|");}
	my $l=0;
	my $read;
	my $first_nt;
	my $len;
	while(<$in>){
		chomp;
		if($l==0){
			my @header= split / /,$_;
			$read=$header[0];
			$read=~s/^\@//;
			$l=1;
			next;
		}
		if($l==1){
			$first_nt = substr($_,0,1);
			$len = length($_);

			# print STDOUT "READ: $read\n\t$_\n\t$first_nt\n\t$len\n";

			## if the read is assigned ...
			if(exists ${data}{$read}){
				#~ print "$read assigned to $data{$read} & $data_wFeature{$read}\n";
				#~ print "\t"
				#~ print "Found $read -> $data{$read}\n";
				$fd{$s}{	$data{$read}	}{$first_nt}{$len}++; # count the number of reads assigned to este gtf by first nt and length
				$fd_wFeature{$s}{	$data{$read}	}{	$data_wFeature{$read}	}{$first_nt}{$len}++; # count the number of reads assigned to este gtf by first nt and length
			}
			$l=2;
			next;
		}
		if($l==2){
			$l=3;
			next;
		}
		if($l==3){
			$l=0;
			next;
		}
	}
}



### print the data
open(my $out,'>',$outfile)|| die print "Cannot open $outfile!\n";

print $out "sample\tannotation\tfirstNt\tlength\tcount\n";

foreach my $sample (keys %fd){
	foreach my $anno (sort keys %{$fd{$sample}}){
		foreach my $firstNt (sort keys %{$fd{$sample}{$anno}}){
			foreach my $len (sort keys %{$fd{$sample}{$anno}{$firstNt}}){
				print $out "$sample\t$anno\t$firstNt\t$len\t$fd{$sample}{$anno}{$firstNt}{$len}\n";
			}
		}
	}
}

### print the data
open(my $out2,'>',$outfile_w_feature)|| die print "Cannot open $outfile_w_feature!\n";

print $out2 "sample\tannotation\tfeature\tfirstNt\tlength\tcount\n";

foreach my $sample (keys %fd_wFeature){
	foreach my $anno (sort keys %{$fd_wFeature{$sample}}){
		foreach my $feat (sort keys %{$fd_wFeature{$sample}{$anno}}){
			foreach my $firstNt (sort keys %{$fd_wFeature{$sample}{$anno}{$feat}}){
				foreach my $len (sort keys %{$fd_wFeature{$sample}{$anno}{$feat}{$firstNt}}){
					#~ print STDOUT "$sample\t$anno\t$feat\t$firstNt\t$len\t$fd_wFeature{$sample}{$anno}{$feat}{$firstNt}{$len}\n";
					print $out2 "$sample\t$anno\t$feat\t$firstNt\t$len\t$fd_wFeature{$sample}{$anno}{$feat}{$firstNt}{$len}\n";
				}
			}
		}
	}
}

print "Done $0!\n";
