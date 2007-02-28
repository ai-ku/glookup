#!/usr/bin/perl -w
use strict;
require 'gngram.pl';
ginit();
my $nword;
my $nbits;
my $logtotal = log2($::GTotal);

while(<>) {
    s/([a-z]) n't/$1n't/g;
    s/(\w)([-\/])(\w)/$1 $2 $3/gi;
    s/(\w)([-.;,\/\'\%\)]) /$1 $2 /gi;
    s/(\S)(\'[s]) /$1 $2 /gi;
    s/ ([\$\%\(\`])(\w)/ $1 $2/g;
    my @s = ('<S>', split(), '</S>');
    my $skip = 0;
    for (my $i = 1; $i < $#s; $i++) {
	if (gngram($s[$i]) == 0) { 
	    warn "Warning[$.]: [$s[$i]] unknown, skipping sentence.\n";
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
    return $logtotal - log2(gngram($s->[$i]));
}
