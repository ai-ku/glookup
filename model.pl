#!/usr/bin/perl -w
use strict;
warn '$Id: model.pl,v 1.12 2007/03/04 13:14:41 dyuret Exp dyuret $' ."\n";

require 'gngram.pl';
use Getopt::Long;
my $verbose = 0;
my $cachefile;
my $ngram = 1;
my @A;
my @B;
GetOptions('cache=s' => \$cachefile,
	   'verbose' => \$verbose,
	   'ngram=i' => \$ngram,
# We took gbits to gngram.pl, the parameters can no longer be set.
# 	   'a2=f' => \$A[2],
# 	   'a3=f' => \$A[3],
# 	   'a4=f' => \$A[4],
# 	   'a5=f' => \$A[5],
# 	   'b2=f' => \$B[2],
# 	   'b3=f' => \$B[3],
# 	   'b4=f' => \$B[4],
# 	   'b5=f' => \$B[5]
);

my $nword;
my $nbits;
my $nline;

while(<>) {
    $nline++;
    my @s = ('<S>', gtokenize($_), '</S>');
    my $skip = 0;
    for (my $i = 1; $i < $#s; $i++) {
	if (gngram($s[$i]) == 0) { 
	    warn "Warning[$nline]: [$s[$i]] unknown, skipping sentence.\n";
	    $skip = 1; last; 
	}
    }
    next if $skip;
    for (my $i = 1; $i < $#s; $i++) {
	$nword++;
	my $b = gbits(\@s, $i, $ngram);
	$nbits += $b;
	if ($verbose) {
	    print $s[$i];
	    for (my $n = 1; $n < $ngram; $n++) {
		printf "\t%.4f", gbits(\@s, $i, $n);
	    }
	    printf "\t%.4f\n", $b;
	}
    }
    print "\n" if $verbose;
}

print $nword . ' ' . ($nbits/$nword) . "\n";

