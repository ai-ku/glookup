#!/usr/bin/perl -w
warn '$Id: model.pl,v 1.14 2007/04/05 23:18:02 dyuret Exp dyuret $' ."\n";

use strict;
use Getopt::Long;
use PDL;
use PDL::Opt::Simplex;
require 'gngram.pl';
require 'fileio.pl';

# Command line options:

my $verbose = 0;
my $optimize = 0;
my $cachefile;
my $ngram = 5;
my $random = 0;
my $zeroes = 0;
# smoothing = { linear, product, baseline }
my $smoothing = 'linear';
my @A = (undef, undef, 8.00, 6.00, 5.40, 4.89);
my @B = (undef, undef, 0.00, 0.00, 0.00, 0.00);
my @C = (undef, undef, 0.01, 0.01, 0.01, 0.01);

GetOptions('cache=s' => \$cachefile,
           'verbose' => \$verbose,
           'optimize' => \$optimize,
           'random' => \$random,
           'zeroes' => \$zeroes,
	   'ngram=i' => \$ngram,
           'smoothing=s' => \$smoothing,
 	   'a2=f' => \$A[2],
 	   'a3=f' => \$A[3],
 	   'a4=f' => \$A[4],
	   'a5=f' => \$A[5],
 	   'b2=f' => \$B[2],
 	   'b3=f' => \$B[3],
 	   'b4=f' => \$B[4],
 	   'b5=f' => \$B[5],
 	   'c2=f' => \$C[2],
 	   'c3=f' => \$C[3],
 	   'c4=f' => \$C[4],
 	   'c5=f' => \$C[5],
);

my %init_fn = ('linear' => \&init_linear, 'product' => \&init_product, 'baseline' => \&init_baseline);
my %score_fn = ('linear' => \&score_linear, 'product' => \&score_product, 'baseline' => \&score_baseline);

# Main:

my $GTotal = ginit($cachefile);
warn "gtotal=$GTotal\n";
my $corpus = read_corpus();
warn "corpus=" . scalar(@$corpus) . " sentences\n";
my $nline = 0;
my $nscore = 0;
my $infinity = 1E9;
my $epsilon = 1E-6;
$optimize ? optimize() : ngram();

# Subroutines:

sub optimize {
    die "Choose ngram >= 2 for optimization"
	unless $ngram >= 2;
    my $init = &{$init_fn{$smoothing}}();
    warn "initial score:\n";
    my $initscore = &{$score_fn{$smoothing}}($init);
    my $initsize = 0.1;
    my $minsize = 1E-4;
    my $maxiter = 1E6;
    my ($optimum, $ssize) = simplex($init, $initsize, $minsize,
				    $maxiter, \&score, \&display);
}

sub ngram {
    my ($nword, $nbits) = (0, 0);
    for my $s (@$corpus) {
	my ($b, $w) = process_sentence($s);
	$nbits += $b; $nword += $w;
    }
    my $avgbits = $nbits / $nword;
    print "$nword $avgbits\n"
	if not $optimize;
    return $avgbits;
}

sub init_linear {
    my $init;
    my $ndims = 2 * $ngram - 2;
    if ($zeroes) {
	$init = zeroes($ndims);
    } elsif ($random) {
	$init = 10 * random($ndims);
    } else {
	my @x;
	for (my $i = 2; $i <= $ngram; $i++) {
	    $x[2 * $i - 4] = $A[$i];
	    $x[2 * $i - 3] = $B[$i];
	}
	$init = pdl(@x);
    }
    return $init;
}

sub init_product {
    my $init;
    my $ndims = $ngram - 1;
    $init = $zeroes ? zeroes($ndims) 
	: $random ? (10 * random($ndims))
	: pdl(@A[2 .. $ngram]);
    return $init;
}

sub init_baseline {
    my $init;
    my $ndims = $ngram - 1;
    $init = $zeroes ? zeroes($ndims) 
	: $random ? random($ndims)
	: pdl(@C[2 .. $ngram]);
    return $init;
}

sub score {
    my $x = shift;
    my ($ndims, $npoints) = $x->dims;
    if (defined $npoints) {
	my @score;
	for (my $i = 0; $i < $npoints; $i++) {
	    my $xi = $x->slice(":,($i)");
	    $score[$i] = &{$score_fn{$smoothing}}($xi);
	}
	return pdl(@score);
    } else {
	return pdl(&{$score_fn{$smoothing}}($x));
    }
}

sub score_linear {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $ngram; $i++) {
	$A[$i] = nonnegative($x->at(2 * $i - 4));
	$B[$i] = nonnegative($x->at(2 * $i - 3));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_product {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $ngram; $i++) {
	$A[$i] = nonnegative($x->at($i - 2));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_baseline {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $ngram; $i++) {
	$C[$i] = zero_one($x->at($i - 2));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub nonnegative {
    my $x = shift;
    $x >= 0 ? $x : 0;
}

sub zero_one {
    my $x = shift;
    ($x < $epsilon) ? $epsilon 
	: (($x > (1 - $epsilon)) ? (1 - $epsilon) : $x);
}

sub display {
    my ($simp, $vals, $ssize) = @_;
    my $best = min($vals);
    warn "best=$best ssize=$ssize\n";
}

sub process_sentence {
    my ($s) = @_;
#    warn "Processing [" . join(' ', @$s) . "]\n";
    $nline++;
    my $nword;
    my $nbits;
    for (my $i = 1; $i < $#{$s}; $i++) {
	$nword++;
	if (gngram($s->[$i]) == 0) { 
	    warn "Warning[$nline]: [$s->[$i]] unknown, skipping sentence.\n"
		if not $optimize;
	    return (0, 0);
	}
	my $b = bits($s, $i, $ngram);
	$nbits += $b;
	if ($verbose) {
	    print $s->[$i];
	    for (my $n = 1; $n < $ngram; $n++) {
		printf "\t%.4f", bits($s, $i, $n);
	    }
	    printf "\t%.4f\n", $b;
	}
    }
    print "\n" if $verbose;
    return ($nbits, $nword);
}

sub read_corpus {
    my @corpus;
    readfile('', sub {
	return unless /\S/;
	my @s = ('<S>', gtokenize($_), '</S>');
	push @corpus, \@s;
    });
    return \@corpus;
}

sub bits {
    my ($s, $i, $n) = @_;
    $n = 1 if not defined $n;
    if ($n == 1) {
	my $g = gngram($s->[$i]);

	# BUG: Any word below 200 count is excluded from the google
	# data.  We will give such words a count of 100.  Note that
	# this may happen because of tokenization errors, and it is
	# not probabilistically sound.

	# die "Unknown word [$s->[$i]]" if $g == 0;
	$g = 100 if $g == 0;

	return log2($GTotal / $g);
    } elsif ($n > $i + 1) {
	return bits($s, $i, $i+1);
    }

    my $a = join(' ', @{$s}[($i-$n+1) .. ($i-1)]); # a = n-1 word prefix
    my $ga = gngram($a);	# ga = count of a
    my $x = bits($s, $i, $n-1);	# x = lower order model bits
    if ($ga == 0) { 
#	    warn "Warning: Zero a-count[$a]\n" 
#		if $n <= 3 and $a =~ /[^\w ]/;		# check what is going on with punctuation
	return $x;		# return lower order result
    }
    my $px = exp2(-$x);		# px = lower order model probability
    my $b = $a . " $s->[$i]";	# b = all n words
    my $gb = gngram($b);	# gb = count of b
    my $c = $a . " :?:";	# c = all ngrams that start with a
    my $gc = gngram($c);	# gc = count of c
    my $missing_count = $ga - $gc; # missing_count = occurances of (a) 
    #   that are missing from ngram data
    my $pb;

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.  In fact the only three cases are [!], [?], [;].  I am
# going to add [.] to the list for the same reason even though there
# seem to be bigrams starting with [.].

    my $bad_context = 0;
    for my $tok (@{$s}[($i-$n+1) .. ($i-1)]) {
	if ($tok =~ /^[;!?.]$/) {
	    $bad_context = 1; last;
	}
    }

    if ($bad_context) {
	$pb = $px;

    } elsif ($missing_count > 0) { # apply our smoothing formula
	if ($smoothing eq 'linear') {
	    my $extra = $A[$n] * $missing_count + exp10($B[$n] - 1.0); # want 0 when b=0
	    $pb = ($gb + $px * ($missing_count + $extra)) 
		/ ($ga + $extra);
	} elsif ($smoothing eq 'product') {
	    my $extra = $A[$n] * $missing_count;
	    $pb = ($gb + $px * ($missing_count + $extra)) 
		/ ($ga + $extra);
	} elsif ($smoothing eq 'baseline') {
	    $pb = $C[$n] * $px;
	    $pb += (1 - $C[$n]) * ($gb / $gc)
		if $gc > 0;
	} else {
	    die "Unknown smoothing $smoothing";
	}

    } elsif ($missing_count == 0) {
#	    warn "Warning: Zero missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	$pb = ($gb + $px) / ($ga + 1);

    } elsif ($missing_count < 0) {

# This is a bug in gngram.  This means there are more instances of
# words following "A" than there are instances of "A" by itself.

	warn "Warning: Negative missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	$pb = ($gb + $px) / ($gc + 1);
    }
    return ($pb > 0 ? -log2($pb) : $infinity);
}


=pod

We will use *.left files to do one-count smoothing.  MacKay and
Peto argue that a reasonable form for a smoothed distribution is:

P_one(w[i] | w[i-n+1..i-1]) = c(w[i-n+1..i]) + a P_one(w[i] | w[i-n+2..i-1])
                              ----------------------------------------------
                                  c(w[i-n+1..i-1]) + a

Where the parameter "a" can be interpreted as the extra counts
added to the given distribution and these extra counts are
distributed as the lower order distribution.  Good-Turing suggests
that the ideal value for "a" should be proportional to the number
of words with exactly one count in the given distribution.

a = g [ n1(w[i-n+1..i-1]) + b ]    ; g, b constants

n1(w[i-n+1..i-1]) = |wi : c(w[i-n+1..i]) = 1|

For us there are no n-grams with count = 1.  Google throws away
all ngrams with count < 40.  However we have a reasonable
alternative: we can count the ngrams that were thrown away by
looking at:

c(w[i-n+1..i-1]) - Sum_w[i] c(w[i-n+1..i])

It is the second term of this expression that the .left files will
give us.  We can use the difference as the "a" multiplier in the
nominator.  We do not need an "a" multiplier in the denominator
because the "a" count is really missing from the n-gram data.

What if the difference is 0?  We then use the regular one count
formula with a = 1.

In the implementation above we use the A and B arrays for parameters.
We first find out the missing count:

missing_count = c(w[i-n+1..i-1]) - Sum_w[i] c(w[i-n+1..i])

This many counts we can use for free because they exist in the
denominator but not the nominator.  We add some extra for more
smoothing:

extra = An * missing_count + 10^Bn

Then the smoothed probability becomes simply:

P_one(w[i] | w[i-n+1..i-1]) = c(w[i-n+1..i]) + (missing_count + extra) P_one(w[i] | w[i-n+2..i-1])
                              --------------------------------------------------------------------
                                  c(w[i-n+1..i-1]) + extra

=cut

