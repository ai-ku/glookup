#!/usr/bin/perl -w
use strict;
require 'gngram.pl';
ginit();
my $nword;
my $nbits;
my $nline;
my $logtotal = log2($::GTotal);
my $SM = 0.96;

while(<>) {
    $nline++;
    s/([a-z]) n't/$1n't/g;
    s/(\w)([-\/])(\w)/$1 $2 $3/gi;
    s/(\w)([-.;,\/\'\%\)]) /$1 $2 /gi;
    s/(\S)(\'[s]) /$1 $2 /gi;
    s/ ([\$\%\(\`])(\w)/ $1 $2/g;
    my @s = ('<S>', split(), '</S>');
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
	$nbits += bits(\@s, $i);
    }
}

printf "%d %f\n", $nword, $nbits/$nword;

sub log2 { log($_[0])/log(2); }

sub bits {
    my ($s, $i) = @_;
    my $n1 = gngram($s->[$i-1]);
    my $n2 = gngram($s->[$i]);
    my $n12 = gngram("$s->[$i-1] $s->[$i]");
    my $p12 = $n12/$n1;
    my $p2 = $n2/$::GTotal;
    return -log2($SM*$p12 + (1-$SM)*$p2);
}
