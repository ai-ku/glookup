#!/usr/bin/perl -w
warn q{$Id: model.pl,v 2.4 2007/04/16 11:31:50 dyuret Exp dyuret $ } ."\n";

use strict;
use Getopt::Long;
use PDL;
use PDL::Opt::Simplex;
require 'gngram.pl';
require 'fileio.pl';
require 'dhc.pl';

# Command line options:

my $verbose = 0;
my $optimize = 0;
my $dhc = 0;
my $patterns;			# output patterns
my $cachefile;
my $ngram = 5;
my $random = 0;
my $zeroes = 0;
# smoothing = { linear, product, baseline, kn }
my $smoothing = 'linear';
my @A = (undef, undef, 6.3712181,  6.2403967, 6.2855943,  5.375136); # 8.06058787993131
my @B = (undef, undef, 0.00,       0.00,      2.4973338,  2.457501); 
my @C = (undef, undef, 0.12244049, 0.4886369, 0.74636033, 0.83561995); # 8.22092294839358
my @D = (undef, undef, 6.7131229,  5.9414447, 6.5528203,  5.7060572); # 8.06083590891475
my @KN = (undef, undef, 0, 0, 0, 0);

GetOptions('cache=s' => \$cachefile,
           'verbose' => \$verbose,
           'optimize' => \$optimize,
	   'dhc' => \$dhc,
	   'patterns' => \$patterns,
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
 	   'd2=f' => \$D[2],
 	   'd3=f' => \$D[3],
 	   'd4=f' => \$D[4],
 	   'd5=f' => \$D[5],
	   'kn2=f' => \$KN[2],
	   'kn3=f' => \$KN[3],
	   'kn4=f' => \$KN[4],
	   'kn5=f' => \$KN[5],
);

my %init_fn = ('linear' => \&init_linear, 'product' => \&init_product, 'baseline' => \&init_baseline, 'kn' => \&init_kn);
my %score_fn = ('linear' => \&score_linear, 'product' => \&score_product, 'baseline' => \&score_baseline, 'kn' => \&score_kn);

# Main:

my $GTotal = $patterns ? 1024908267229 : ginit($cachefile);
warn "gtotal=$GTotal\n";
my $corpus = read_corpus();
warn "corpus=" . scalar(@$corpus) . " sentences\n";
my $nline = 0;
my $nscore = 0;
my $debug = 0;
my $infinity = 1E9;
my $epsilon = 1E-6;
my $MINNGRAM = 40;
my %myngram;
$optimize ? optimize() : ngram();

# Subroutines:

sub optimize {
    die "Choose ngram >= 2 for optimization"
	unless $ngram >= 2;
    my $init = &{$init_fn{$smoothing}}();
#    warn "init=$init\n";
    warn "initial score:\n";
    my $initscore = &{$score_fn{$smoothing}}($init);
    my $initsize = 1.0;
    my $minsize = 1E-4;
    my $maxiter = 1E6;
    if ($dhc) {
	my ($optimum, $final) = dhc(\&score, $init);
	warn "$optimum <= $final\n";
    } else {
	my ($optimum, $ssize) = simplex($init, $initsize, $minsize,
					$maxiter, \&score, \&display);
    }
}

sub ngram {
    my ($nword, $nbits) = (0, 0);
    for my $s (@$corpus) {
	my ($b, $w) = process_sentence($s);
	$nbits += $b; $nword += $w;
    }
    my $avgbits = $nbits / $nword;
    warn "$nword $avgbits\n"
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
	: pdl(@D[2 .. $ngram]);
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

sub init_kn {
    my $init;
    my $ndims = $ngram - 1;
    $init = $zeroes ? zeroes($ndims) 
	: $random ? log(1/random($ndims)-1)
	: pdl(@KN[2 .. $ngram]);
    return $init;
}

sub score {
    my $x = shift;
#    warn "score($x)\n";
    my ($ndims, $npoints) = $x->dims;
#    warn "ndims=$ndims npoints=$npoints\n";
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
	$D[$i] = nonnegative($x->at($i - 2));
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

sub score_kn {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $ngram; $i++) {
	$KN[$i] = $x->at($i - 2);
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
    if (!$patterns) {
	for (my $i = 1; $i < $#{$s}; $i++) {
	    if (myngram($s->[$i]) == 0) { 
		warn "Warning[$nline]: [$s->[$i]] unknown, skipping sentence.\n"
		    if not $optimize;
		return (0, 0);
	    }
	}
    }
    for (my $i = 1; $i < $#{$s}; $i++) {
	$nword++;
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
    $n = 5 if not defined $n;
    if ($smoothing =~ /^kn/) {
	my $pb = kn($s, $i, $n);
	die "Kneser-Ney returned $pb" if $pb <= 0;
	return -log2($pb);
    }
    if ($n == 1) {
	my $g = myngram($s->[$i]);

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
    my $ga = myngram($a);	# ga = count of a
    my $x = bits($s, $i, $n-1);	# x = lower order model bits
    if ($ga == 0) { 
#	    warn "Warning: Zero a-count[$a]\n" 
#		if $n <= 3 and $a =~ /[^\w ]/;		# check what is going on with punctuation
	return $x;		# return lower order result
    }
    my $px = exp2(-$x);		# px = lower order model probability
    my $b = $a . " $s->[$i]";	# b = all n words
    my $gb = myngram($b);	# gb = count of b
    my $c = $a . " :?:";	# c = all ngrams that start with a
    my $gc = myngram($c);	# gc = count of c
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
	    my $extra = $D[$n] * $missing_count;
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


sub kn {
    my ($s, $i, $n) = @_;
    my ($nxy, $x, $nx_, $n1x_, $kn, $D);
    $n = 5 if not defined $n;
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {

 	# BUG: no need for 0 order model, if the word is known, this
 	# is exactly equivalent to k/N.  If the word is not known it
 	# is cheating.
	
 	# return 1/myngram('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $g = myngram($s->[$i]);

	# BUG: Any word below 200 count is excluded from the google
	# data.  We will give such words a count of 100.  Note that
	# this may happen because of tokenization errors, and it is
	# not probabilistically sound.

	# die "Unknown word [$s->[$i]]" if $g == 0;
	$g = 100 if $g == 0;

	return $g / $GTotal;
	
    } elsif ($n > 5) {
	return kn($s, $i, 5);
    } elsif ($n > $i + 1) {
	return kn($s, $i, $i+1);
    }    
    $x = ($n > 1) ? join(' ', @{$s}[($i-$n+1) .. ($i-1)]).' ' : '';
    $nx_ = myngram($x . ':?:');
    warn("kn($s->[$i],$i,$n): nx_ = $nx_\n") if $debug;
    if ($nx_ == 0) {
	return kn($s, $i, $n-1);
    } elsif (" $x" =~ / [;!?.] /) { # bad context
	return kn($s, $i, $n-1);
    }	
    $n1x_ = myngram($x . '_');
    warn("kn($s->[$i],$i,$n): n1x_ = $n1x_\n") if $debug;
    $nxy = myngram($x . $s->[$i]);
    warn("kn($s->[$i],$i,$n): nxy = $nxy\n") if $debug;
    # $D = 1/(1+exp(-$KN[$n]));
    $D = $MINNGRAM/(1+exp(-$KN[$n]));
    warn("kn($s->[$i],$i,$n): D = $D\n") if $debug;
    $nxy -= $D if $nxy > 0;
    #$kn = $nxy / $nx_ + ($n1x_ * $D / $nx_) * kn0($s, $i, $n-1);

    # new formula:
    my ($nx, $mc);
    if ($n == 1) {
	$nx = $nx_;
	$mc = 0;
    } else {
	$nx = myngram($x);
	$mc = $nx - $nx_;
    }
    warn("kn($s->[$i],$i,$n): x = [$x] nx = $nx\n") if $debug;
    warn("kn($s->[$i],$i,$n): mc = $mc\n") if $debug;
    $kn = $nxy / $nx + (($mc + $n1x_ * $D) / $nx) * kn0($s, $i, $n-1);

    warn("kn($s->[$i],$i,$n): kn = $kn\n") if $debug;
    if ($patterns) { 
	# if outputting patterns do the lower order models as well
	kn($s, $i, $n-1); 
    }
    return $kn;
}

# BUG: The n1 counts are grossly underestimated whereas the n counts
# are exact.  This especially effects the second term in kn which is a
# n1/n term.  Must derive corrected n1 counts from the missing count.

sub kn0 {
    my ($s, $i, $n) = @_;
    my ($x, $n1_xy, $n1_x_, $n1x_, $D, $kn0);
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {
	# return 1/myngram('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = myngram('_ ' . $s->[$i]);
	if ($n1_y == 0) {
	    # warn "n1_y($s->[$i]) == 0";
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = myngram('_ _');
	return $n1_y / $n1__;
    }
    $x = ($n > 1) ? join(' ', @{$s}[($i-$n+1) .. ($i-1)]).' ' : '';
    $n1_x_ = myngram('_ ' . $x . '_');
    warn("kn0($s->[$i],$i,$n): n1_x_ = $n1_x_\n") if $debug;
    die "KN0 zero denominator [$x]" if $n1_x_ == 0;
    $n1x_ = myngram($x . '_');
    warn("kn0($s->[$i],$i,$n): n1x_ = $n1x_\n") if $debug;
    $n1_xy = myngram('_ ' . $x . $s->[$i]);
    warn("kn0($s->[$i],$i,$n): n1_xy = $n1_xy\n") if $debug;
    $D = 1/(1+exp(-$KN[$n]));
    warn("kn0($s->[$i],$i,$n): D = $D\n") if $debug;
    $n1_xy -= $D if $n1_xy > 0;
    $kn0 = $n1_xy / $n1_x_ + ($n1x_ * $D / $n1_x_) * kn0($s, $i, $n-1);

    # new formula: does not really work
#     my $n1x = myngram($x);
#     my $mc = $n1x - $n1_x_;
#     $kn0 = $n1_xy / $n1x + (($mc + $n1x_ * $D) / $n1x) * kn0($s, $i, $n-1);

    warn("kn0($s->[$i],$i,$n): kn0 = $kn0\n") if $debug;
    return $kn0;
}

sub myngram {
    my $str = shift;
    if ($patterns) {
	if (not defined $myngram{$str}) {
	    $myngram{$str} = 1;
	    print "$str\n";
	}
	return 1;
    } else {

	# BUG: temporary fix
	$str =~ s/^\s+//;
	$str =~ s/\s+$//;
	if ($str =~ /^<S>$/) {
	    return 95119665584;
	}
	# End of temporary fix

	return gngram($str);
    }
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

