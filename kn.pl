#!/usr/bin/perl -w
use strict;
warn q{$Id: kn.pl,v 1.6 2008/01/29 11:50:53 dyuret Exp dyuret $ }."\n";
my $DEBUG = 0;
my $ANY = '_';
my $NGRAM = 5;
my $MINNGRAM = 40;
my %myngram;
my %MISSING;

# The following parameters have been optimized on the Brown corpus:

my @KNMOD = 			# knmod: 1k=7.80248865053247 full=7.85466149603376
    (
     1.4351794, 		# coef=(A+1)/(A+40) for kn2 (7.8536 if 0)
     1.3009039, 2.0350656, 2.4501021, # C2..C4 (x n1modx_) for kn2; D=1 for kn2 (8.0358 if 0)
     3.0263603, 2.9146123, 3.1023827, 3.5312749, # C2..C5 (x mc) for kn (8.8599 if 0)
     0.17869642, 3.3607215e-05, 0.23590772, 0.22686745 # D2..D5 (x40) for kn (7.8035 if 0)
     );

# TODO:
# - test using backoff vs regular counts when context has zero count
# - do something better about oov words
# - test the kn discount coefficient estimates given in paper
# - test using n0x vs n0x_
# - test using always modified counts.

# kn(s,i,n): returns the probability of the i'th word in sentence s
# according to an n-gram model.  n defaults to 5.  The first token in
# the sentence array should be <S>.  The model needs to be initialized
# (kn_init) with a file of counts before use.  If this is not
# performed, then the program simply writes out the patterns that
# would be needed for its computation into the hash %myngram

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
    my $y = $s->[$i];

    if ($n == 1) {
	my $ny = kn0($y);
	if ($ny == 0) {
	    # BUG: Any word below 200 count is excluded from the google
	    # data.  We will give such words a count of 100.  Note that
	    # this may happen because of tokenization errors.  It is
	    # not probabilistically sound.
	    $ny = 100;
	}
	my $n_ = kn0($ANY);
	$p = $ny / $n_;
	warn "kn($y,$i,$n): $ny/$n_=$p\n" if $DEBUG;

    } else {

	# See if the n-1 word context x is any good

	my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
	my $n0x_ = kn0("$x _");
	my $n0x  = kn0($x);
	my $mc = $n0x - $n0x_;

	if (($n0x == 0)
	    # Backoff when context count is zero
	    # BUG: instead of backing off here, we should try using the
	    # <UNK> token for unknown words in context.
	    or ((" $x " =~ / [;!?.] /) and
		not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
	    # bad context: google data does not have any regular tokens
	    # following these four characters.
	    or (not %kn0)
	    # If kn_init has not been called we do not know if we need
	    # to backoff or not, so we might as well do every order.
	    ) {

	    # This is my original implementation, use the unmodified
	    # counts when backing off because of a zero context.
	    $p = kn($s, $i, $n-1); 

	    warn("kn($y,$i,$n): nx=$n0x nx_=$n0x_ skip [$x]\n") if $DEBUG;

	    # This is the srilm implementation, always use the
	    # modified context when backing off.
	    # my $p = kn_backoff($s, $i, $n-1); 

	    return $p if %kn0;
	}

	my $p0 = kn_backoff($s, $i, $n-1);
	warn("kn($y,$i,$n): kn0 = $p0\n") if $DEBUG;

	# First approximation

	my $n0xy = kn0("$x $y");
	warn("kn($y,$i,$n): ML: nxy = $n0xy / nx = $n0x => ".($n0xy/$n0x)."\n") if $DEBUG;
	warn("kn($y,$i,$n): nx=$n0x nx_=$n0x_ mc=$mc\n") if $DEBUG;

	my $D = $MINNGRAM * $KNMOD[$n+6];   # parameter to subtract from the n0xy count
	warn("kn($y,$i,$n): D = $D\n") if $DEBUG;
	$n0xy -= $D if $n0xy > 0;	# we only subtract if it is positive
	# $p = $n0xy / $n0x_;	# first approximation
	#warn "kn[$n][1]=($n0xy-$D)/$n0x_=$p\n" if $DEBUG;

	# Implement Kneser-Ney smoothing

	# Missing count: The source of the missing count is the $D we
	# subtract for each y with n0xy > 0.  The following gives us
	# the missing count from the D subtractions:

	my $n1x_ = kn1("$x $ANY"); # number of distinct words following x
	# my $mc1 = $D * $n1x_;	# missing count from D subtractions
	warn("kn($y,$i,$n): n1x_ = $n1x_\n") if $DEBUG;
	warn("kn($y,$i,$n): x = [$x] nx = $n0x\n") if $DEBUG;
	warn("kn($y,$i,$n): mc = $mc\n") if $DEBUG;

	# Do some more missing count discounting
	my $extra = $KNMOD[$n+2] * $mc;
	$extra = 1 if 0 == $extra;

	$p = ($n0xy + $p0 * ($extra + $n1x_ * $D)) / ($n0x_ + $extra);
    }

    warn("kn($y,$i,$n): kn = $p\n") if $DEBUG;
    return $p;
}

sub kn_backoff {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = ($KNMOD[0] + 1)/($KNMOD[0] + 40);
    my $p;

    if ($n == 1) {
	my $n1_y = kn1("$ANY $y");
	if ($n1_y == 0) {
	    warn("kn2($y,$i,$n): n1_y($y) == 0\n") if $DEBUG;
	    # Making this 1 is exactly what would happen if we smoothed with a
	    # 0-order 1/V model.  It is not probabilistically sound.
	    $n1_y = 1;
	}
	my $n1__ = kn1("$ANY $ANY");
        $p = $n1_y / $n1__;
	warn("kn2($y,$i,$n): n1_y=$n1_y / n1__=$n1__ => $p\n") if $DEBUG;

    } else {

	my $p0 = kn_backoff($s, $i, $n-1);

	my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
	my $n1_x_ = kn1mod("$ANY $x $ANY", $coef);

	# See if the context has non-zero count:
	if (($n1_x_ == 0)
	    # Backoff when context count is zero
	    # BUG: instead of backing off here, we should try using the
	    # <UNK> token for unknown words in context.
	    or ((" $x " =~ / [;!?.] /) and
		not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
	    # bad context: google data does not have any regular tokens
	    # following these four characters.
	    ) {
	    warn("kn2($y,$i,$n): n1_x_ = $n1_x_ skip [$x]\n") if $DEBUG;
	    return $p0;
	}    

	# First approximation

	my $n1_xy = kn1mod("$ANY $x $y", $coef);
	die "n1_xy=$n1_xy" if ($n1_xy > 0 and $n1_xy < 1);
	warn("kn2($y,$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $DEBUG;

	my $D = 1;	# parameter to subtract
	warn("kn2($y,$i,$n): D = $D\n") if $DEBUG;
	$n1_xy -= $D if $n1_xy > 0; # only subtract if positive count
	# $p = $n1_xy / $n1_x_;	# first approximation
	# warn "kn_backoff[$n][0]=$n1_xy/$n1_x_=$p\n" if $DEBUG;

	my $n1x_ = kn1("$x _");
	warn("kn2($y,$i,$n): n1x_ = $n1x_\n") if $DEBUG;
	my $C = $KNMOD[$n-1];
	my $n1modx_ = kn1mod("$x _");

	$p = ($n1_xy + $p0 * ($C * $n1modx_ + $D * $n1x_)) / ($C * $n1modx_ + $n1_x_);

    }
    warn("kn2($y,$i,$n): kn2 = $p\n") if $DEBUG;
    return $p;
}

sub kn0 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $myngram{$pat}++;
    } 
    if (defined $kn0{$pat}) {
	return $kn0{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless %MISSING;
	$MISSING{$pat}++;
	return 0;
    }
}

sub kn1 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $myngram{$pat}++;
    } 
    if (defined $kn1{$pat}) {
	return $kn1{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless %MISSING;
	$MISSING{$pat}++;
	return 0;
    }
}

sub kn2 {
    my $pat = shift;
    if (not %kn0) {
	print "$pat\n" unless $myngram{$pat}++;
    } 
    if (defined $kn2{$pat}) {
	return $kn2{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless %MISSING;
	$MISSING{$pat}++;
	return 0;
    }
}

sub kn1mod {
    my ($pat, $coef) = @_;
    if (defined $coef) {
	die "Bad n1mod coefficient $coef" if $coef < 1/40;
    } else {
	$coef = 1/40;
    }
    my $n1 = kn1($pat);
    my $n0 = kn0($pat);
    $pat =~ s/\s*_\s*//;
    my $n0x = kn0($pat);
    unless ($n0x >= $n0 and $n0 >= $n1) {
	warn "kn1mod: Inconsistent count [$pat] n1=$n1 n0=$n0 n0x=$n0x\n";
	return $n1;
    }
    return $n1 + $coef * ($n0x - $n0);
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
	for (my $i = 1; $i < $#s; $i++) {
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
    if (%kn0 and %MISSING) {
	print STDERR "There were missing patterns:\n";
	for my $pat (keys %MISSING) {
	    print STDERR "$pat\n";
	}
    }
    $cnt{entropy} = $cnt{bits} / $cnt{words};
    warn Dumper(\%cnt);
}

1;
