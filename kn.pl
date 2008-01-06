#!/usr/bin/perl -w
use strict;
warn q{$Id: kn.pl,v 1.4 2007/12/02 10:37:48 dyuret Exp dyuret $ }."\n";
my $DEBUG = 0;

# TODO:
# - test using backoff vs regular counts when context has zero count
# - do something better about oov words
# - test the kn discount coefficient estimates given in paper
# - test using n0x vs n0x_

# kn(s,i,n): returns the probability of the i'th word in sentence s
# according to an n-gram model.  n defaults to 5.  The first token in
# the sentence array should be <S>.  The model needs to be initialized
# (kn_init) with a file of counts before use.  If this is not
# performed, then the program simply writes out the patterns that
# would be needed for its computation into the hash %main::kn_patterns

# BUG: For unknown words we assume a count of 100, which is half of
# the minimum count for google tokens.
# Change this: unknown words should return 0 probability, the caller
# should handle them, replacing with <UNK> if necessary.

# The program uses three types of counts:
# n0: is the total count of the ngrams matching a given pattern
# n1: is the number of distinct ngrams matching a given pattern
# n2: is the number of distinct words that appear on the right
# side of a pattern in the form of "x _".  X itself may have
# wildcards; this is needed for kn smoothing.

my %kn0;
my %kn1;
my %kn2;

# The format of the file that kn_init uses to initalize should have
# the following tab separated fields:

# PATTERN <tab> N0 [<tab> N1 [<tab> N2]]

# The pattern is an ngram string, tokens separated by spaces, possibly
# some wildcards.  N0 always exists, N1 only exists if the pattern
# contains a wildcard, and N2 only exists if the pattern has more than
# one wildcard with one at the end.  N1 or N2 may be missing if N0==0.

sub kn_init {
    my $path = shift;
    $DEBUG = shift;
    undef $main::kn_patterns;
    return if not defined $path;
    open(FP, $path) or die $!;
    while(<FP>) {
	chomp;
	my ($pat, $n0, $n1, $n2) = split(/\t/);
	die if not defined $pat;
	die if not defined $n0;
	$kn0{$pat} = $n0;
	$kn1{$pat} = $n1 if defined $n1;
	$kn2{$pat} = $n2 if defined $n2;
	print STDERR '.' unless $. % 100000;
    }
    close(FP);
}

# The following parameters have been optimized on the Brown corpus:

my @D40 = ( undef, undef, 28.11096928644257128068, 33.98733562666299928484, 39.16948427744310745367, 39.81172633310000762000 );
my @D1 = ( undef, undef, .85948550385118495566, .96034527604501211555, .99565930386648919162, undef );

my $ANY = '_';
my $NGRAM = 5;

sub kn {
    my ($s, $i, $n) = @_;
    $n = $NGRAM if not defined $n;

    # Handle exceptions

    if (($i < 0) || ($i > $#{$s})) {
	die "kn: Position out of range";

    } elsif ($n <= 0) {
	die "kn: Bad n [$n]";

    } elsif ($s->[$i] eq $ANY) {
	# The probability of seeing any word in position i is 1
	return 1;

    } elsif ($n > $NGRAM) {
	return kn($s, $i, $NGRAM);

    } elsif ($n > $i + 1) {
	return kn($s, $i, $i+1);

    } 

    # Calculate probability: in the following x stands for the context
    # and y stands for the target word.

    my $p;

    if ($n == 1) {
	my $ny = kn0($s->[$i]);
	if ($ny == 0) {
	    # BUG: Any word below 200 count is excluded from the google
	    # data.  We will give such words a count of 100.  Note that
	    # this may happen because of tokenization errors.  It is
	    # not probabilistically sound.
	    $ny = 100;
	}
	my $n_ = kn0($ANY);
	$p = $ny / $n_;

    } else {

	# See if the n-1 word context x is any good

	my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
	my $n0x_ = kn0("$x _");
	if (($n0x_ == 0)
	    # Backoff when context count is zero
	    # BUG: instead of backing off here, we should try using the
	    # <UNK> token for unknown words in context.
	    or ((" $x " =~ / [;!?.] /) and
		not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
	    # bad context: google data does not have any regular tokens
	    # following these four characters.
	    or (not %kn0)
	    # If kn_init has not been called we do not know if we need
	    # to backoff or not, so we might as well do every order.
	    ) {

	    # This is my original implementation, use the unmodified
	    # counts when backing off because of a zero context.
	    $p = kn($s, $i, $n-1); 

	    # This is the srilm implementation, always use the
	    # modified context when backing off.
	    # my $p = kn_backoff($s, $i, $n-1); 

	    return $p if %kn0;
	}

	# First approximation

	my $n0xy = kn0("$x $s->[$i]");
	warn "kn[$n][0]=$n0xy/$n0x_=" .($n0xy/$n0x_)."\n" if $DEBUG;
	my $D = $D40[$n];   # parameter to subtract from the n0xy count
	$n0xy -= $D if $n0xy > 0;	# we only subtract if it is positive
	$p = $n0xy / $n0x_;	# first approximation
	warn "kn[$n][1]=($n0xy-$D)/$n0x_=$p\n" if $DEBUG;

	# Implement Kneser-Ney smoothing

	# Missing count: The source of the missing count is the $D we
	# subtract for each y with n0xy > 0.  The following gives us
	# the missing count from the D subtractions:

	my $n1x_ = kn1("$x $ANY"); # number of distinct words following x
	my $mc1 = $D * $n1x_;	# missing count from D subtractions

	# DEPRECATED:
	# By using n0x_ instead of n0x we do not need mc2 any more.
	#
	# Second source of the missing count is that the google data has
	# filtered out all ngrams below the count of 40.  Thus there is a
	# missing count between all ngrams that start with x and the count
	# of the n-1 gram x itself:
	# my $nx_ = kn0("$x $ANY");
	# my $mc2 = $nx - $nx_;

	# Add the kn smoothing term:

	$p += ($mc1 / $n0x_) * kn_backoff($s, $i, $n-1);
    }

    warn "kn[$n]=$p\n" if $DEBUG;
    return $p;
}

sub kn_backoff {
    my ($s, $i, $n) = @_;
    my $p;

    if ($n == 1) {
	my $n1_y = kn1("$ANY $s->[$i]");
	if ($n1_y == 0) {
	    # Making this 1 is exactly what would happen if we smoothed with a
	    # 0-order 1/V model.  It is not probabilistically sound.
	    $n1_y = 1;
	}
	my $n1__ = kn1("$ANY $ANY");
        $p = $n1_y / $n1__;

    } else {

	my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
	my $n1_x_ = kn1("$ANY $x $ANY");

	# See if the context has non-zero count:
	if ($n1_x_ == 0) {
	    # BUG: instead of backing off here, we should try using the
	    # <UNK> token for unknown words in context.
	    return kn_backoff($s, $i, $n-1);
	}    

	# First approximation

	my $n1_xy = kn1("$ANY $x $s->[$i]");
	my $D = $D1[$n];	# parameter to subtract
	$n1_xy -= $D if $n1_xy > 0; # only subtract if positive count
	$p = $n1_xy / $n1_x_;	# first approximation
	warn "kn_backoff[$n][0]=$n1_xy/$n1_x_=$p\n" if $DEBUG;

	# Implement Kneser-Ney smoothing

	# Missing count: if we add up n1_xy for all y, the total will
	# be less than n1_x_, and we will use the difference for
	# smoothing.  The source of the missing count is the $D we
	# subtract for each y with n1_xy > 0.  How many such y are
	# there?  We can't use n1x_, because of the filtering of
	# ngrams with counts less than 40, n1x_ is not equal to n2_x_.
	# We can't use n1_x_ because that counts each axb for distinct
	# a and b.  We need a special counter for the number of y with
	# n1_xy > 0.  That is the n2_x_ count.

	my $n2_x_ = kn2("$ANY $x $ANY");
	my $mc = $D * $n2_x_;

	# Add the kn smoothing term:

	$p += ($mc / $n1_x_) * kn_backoff($s, $i, $n-1);
    }
    warn "kn_backoff[$n]=$p\n" if $DEBUG;
    return $p;
}

sub kn0 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $main::kn_patterns{$pat}++;
	return 1;
    } elsif (defined $kn0{$pat}) {
	return $kn0{$pat};
    } else {
	die "Error: n0 for pattern not found [$pat]";
    }
}

sub kn1 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $main::kn_patterns{$pat}++;
	return 1;
    } elsif (defined $kn1{$pat}) {
	return $kn1{$pat};
    } elsif (defined $kn0{$pat}) {
	if ($kn0{$pat} == 0) {
	    return 0;
	} elsif ($pat !~ /_/) {
	    return 1;
	} else {
	    die "Error: n1 for pattern not found [$pat]";
	}
    } else {
	die "Error: n0 for pattern not found [$pat]";
    }
}

sub kn2 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $main::kn_patterns{$pat}++;
	return 1;
    } elsif (defined $kn2{$pat}) {
	return $kn2{$pat};
    } elsif ($pat !~ /^.*_$/) {
	die "Error: n2 pattern does not end with wildcard [$pat]";
    } elsif (defined $kn1{$pat}) {
	if ($pat =~ /^[^_]*_$/) {
	    return $kn1{$pat};
	} else {
	    die "Error: n2 for pattern not found [$pat]";
	}
    } elsif (defined $kn0{$pat}) {
	if ($kn0{$pat} == 0) {
	    return 0;
	} else {
	    die "Error: n1 for pattern not found [$pat]";
	}
    } else {
	die "Error: n0 for pattern not found [$pat]";
    }
}

use Data::Dumper;
require "gtokenize.pl";

sub kn_test {
    my ($path, $ngram) = @_;
    $ngram = $NGRAM if not defined $ngram;
    my %cnt;

    # We assume a plain text file with one sentence per line
    open(FP, $path) or die $!;
    while(<FP>) {
	next unless /\S/;
	$cnt{sentences}++;
	my @s = ('<S>', gtokenize($_), '</S>');
	for (my $i = 1; $i <= $#s; $i++) {
	    $cnt{words}++;
	    print $s[$i] if %kn0;
	    print '?' if %kn0 and kn0($s[$i]) == 0;
 	    for (my $n = 1; $n < $ngram; $n++) {
		my $p =  kn(\@s, $i, $n);
		my $b = -log($p)/log(2);
 		printf "\t%.4f", $b if %kn0;
 	    }
	    my $p = kn(\@s, $i, $ngram);
	    my $b = -log($p)/log(2);
	    $cnt{bits} += $b;
	    printf "\t%.4f\n", $b if %kn0;
	}
    }
    close(FP);
    if (not %kn0) {
	print "$_\n" for keys %main::kn_patterns;
    }
    warn Dumper(\%cnt);
}

1;
