#!/usr/bin/perl -w
warn q{$Id: model.pl,v 3.20 2008/01/29 22:09:24 dyuret Exp dyuret $ } ."\n";

use strict;
use Getopt::Long;
use Pod::Usage;
use PDL;
use PDL::Opt::Simplex;
use Data::Dumper;
require 'gtokenize.pl';
require 'fileio.pl';

my $oldfh = select(STDERR); $| = 1; select($oldfh); # fix pdl bug
my $log2 = log(2);
sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }
my $log10 = log(10);
#sub log10 { log($_[0])/$log10; } # defined in PDL
sub exp10 { exp($_[0]*$log10); }

# Discounting parameters:

my @C = (undef, undef, 0.12264181, 0.48531058, 0.73285371, 0.8485226); # baseline: 8.20830869973577
my @D = (undef, undef, 7.8440119, 6.6305746, 6.9277921, 7.3675802); # mc: 8.0477965326711
my @KN = (0.83456664, 0.80501932, 0.81691654, 0.90655321, 0.97261754, 0.96917268, 0.97422316); # kn: 8.23974660099736
my @E = (undef, undef, 0.58301752, 0.79111861, 0.92800359, 0.9840277); # cdiscount: 
my @KNMC = (0.83456664, 0.80501932, 0.81691654, 0.90655321, 0.97261754, 0.96917268, 0.97422316); # knmc: 
my @KNMOD = 			# knmod: 1k=7.80248865053247 full=7.85466149603376
    (
     1.4351794, 		# coef=(A+1)/(A+40) for kn2 (7.8536 if 0)
     1.3009039, 2.0350656, 2.4501021, # C2..C4 (x n1modx_) for kn2; D=1 for kn2 (8.0358 if 0)
     3.0263603, 2.9146123, 3.1023827, 3.5312749, # C2..C5 (x mc) for kn (8.8599 if 0)
     0.17869642, 3.3607215e-05, 0.23590772, 0.22686745 # D2..D5 (x40) for kn (7.8035 if 0)
     );

# Wildcard counts:
#                        n0              n1              n2
#         _              1024908267229   13588391        -
#         _ _            910884463583    314843401       12880878
#         _ _ _          739006848674    977069902       9329099
#         _ _ _ _        528435661704    1333820466      7934066
#         _ _ _ _ _      368402367169    1214460675      7015505

# Optimization functions:

my %init_fn = 
    ('mc' => \&init_mc, 
     'baseline' => \&init_baseline, 
     'kn' => \&init_kn, 
     'cdiscount' => \&init_cdiscount, 
     'knmc' => \&init_knmc,
     'knmod' => \&init_knmod,
     );

my %score_fn = 
    ('mc' => \&score_mc, 
     'baseline' => \&score_baseline, 
     'kn' => \&score_kn, 
     'cdiscount' => \&score_cdiscount,
     'knmc' => \&score_knmc,
     'knmod' => \&score_knmod,
     );

# Globals:

my $ZERO_WARNING;
my (%n0, %n1, %n2);
my $nline = 0;
my $nscore = 0;
my $infinity = 1E9;
my $epsilon = 1E-6;
my $MINNGRAM = 40;
my %myngram;
my @corpus;


# Configuration options:

my %config = 
    (
     counts => undef,
     debug => 0,
     dhc => 0,
     mincnt => 0,
     ngram => 5,
     optimize => 0,
     patterns => 0,
     random => 0,
     simplex => 0,
     smoothing => 'knmod', # { baseline, mc, mcmod, kn, cdiscount, knmc, wbdiscount, knmod }
     string => '',
     verbose => 0,
     verify => '',
     zeroes => 0,
     );

# Main:

main() if $0 =~ /\bmodel.pl$/;

sub main {

    GetOptions
	(\%config,
	 'counts=s',
	 'debug',
	 'dhc',
	 'mincnt=i',
	 'ngram=i',
	 'optimize',
	 'patterns',
	 'random',
	 'simplex',
	 'smoothing=s',
	 'string=s',
	 'verbose',
	 'verify=s',
	 'zeroes',
	 ) or pod2usage(2);

    if ($config{counts}) {
	warn sprintf("read_counts($config{counts})=>%d\n", read_counts($config{counts}));
    } elsif (not $config{patterns}) {
	warn "Count file not specified, turning on patterns\n";
	$config{patterns} = 1;
    }
    if ($config{string}) {
	my @s = split(' ', $config{string});
	my $i = $#s;
	print $s[$i];
	for (my $n = 1; $n <= $config{ngram}; $n++) {
	    printf "\t%.4f", bits(\@s, $i, $n);
	}
	printf "\n";
	exit;
    }

    if ($config{verify}) {
	my @vocab;
	my @prob;
	readfile('zcat vocab.gz|', sub { s/\t.*\n//; push @vocab, $_; });
	my @s = split(' ', $config{verify});
	my $i = $#s;
	for my $w (@vocab) {
	    print $w;
	    $s[$i] = $w;
	    for (my $n = 1; $n <= $config{ngram}; $n++) {
		my $p = exp2(-bits(\@s, $i, $n));
		$prob[$n-1] += $p;
		printf "\t%g", $p;
	    }
	    print "\n";
	}
	warn Dumper(\@prob);
	exit;
    }

    read_corpus();
    warn "corpus=" . scalar(@corpus) . " sentences\n";
    $config{optimize} ? optimize() : ngram();
}

# Subroutines:

sub config {
    my ($href) = @_;
    while (my ($opt, $val) = each %$href) {
	$config{$opt} = $val;
    }
}

sub optimize {
    die "Choose ngram >= 2 for optimization"
	unless $config{ngram} >= 2;
    my $init = &{$init_fn{$config{smoothing}}}();
#    warn "init=$init\n";
    warn "initial score:\n";
    my $initscore = &{$score_fn{$config{smoothing}}}($init);
    my $initsize = 1.0;
    my $minsize = 1E-4;
    my $maxiter = 1E6;
    if ($config{simplex}) {
	my ($optimum, $ssize) = simplex($init, $initsize, $minsize,
					$maxiter, \&score, \&display);
    } else {			# dhc is the default
	my ($optimum, $final) = dhc(\&score, $init);
	warn "$optimum <= $final\n";
    }
}

sub ngram {
    my ($nword, $nbits) = (0, 0);
    for my $s (@corpus) {
	my ($b, $w) = process_sentence($s);
	$nbits += $b; $nword += $w;
    }
    my $avgbits = $nbits / $nword;
    warn "$nword $avgbits\n"
	if not $config{optimize};
    return $avgbits;
}

sub init_mc {
    my $init;
    my $ndims = $config{ngram} - 1;
    $init = $config{zeroes} ? zeroes($ndims) 
	: $config{random} ? (10 * random($ndims))
	: pdl(@D[2 .. $config{ngram}]);
    return $init;
}

sub init_baseline {
    my $init;
    my $ndims = $config{ngram} - 1;
    $init = $config{zeroes} ? zeroes($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@C[2 .. $config{ngram}]);
    return $init;
}

sub init_cdiscount {
    my $init;
    my $ndims = $config{ngram} - 1;
    $init = $config{zeroes} ? zeroes($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@E[2 .. $config{ngram}]);
    return $init;
}

sub init_kn {
    my $init;
    my $ndims = 2 * $config{ngram} - 3;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KN[0 .. $ndims-1]);
    return $init;
}

sub init_knmc {
    my $init;
    my $ndims = 2 * $config{ngram} - 3;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KNMC[0 .. $ndims-1]);
    return $init;
}

sub init_knmod {
    my $init;
    my $ndims = 3 * $config{ngram} - 3;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KNMOD[0 .. $ndims-1]);
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
	    $score[$i] = &{$score_fn{$config{smoothing}}}($xi);
	}
	return pdl(@score);
    } else {
	return pdl(&{$score_fn{$config{smoothing}}}($x));
    }
}

sub score_mc {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $config{ngram}; $i++) {
	$D[$i] = nonnegative($x->at($i - 2));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_baseline {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $config{ngram}; $i++) {
	$C[$i] = zero_one($x->at($i - 2));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_cdiscount {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $config{ngram}; $i++) {
	$E[$i] = zero_one($x->at($i - 2));
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_kn {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	if ($x->at($i) < 0 or $x->at($i) > 1) {
	    return $infinity;
	}
	$KN[$i] = $x->at($i);
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_knmc {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	if ($x->at($i) < 0 or $x->at($i) > 1) {
	    return $infinity;
	}
	$KNMC[$i] = $x->at($i);
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_knmod {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	my $xi = $x->at($i);
	if ($xi < 0) {
	    return $infinity;
	} elsif ($i >= $x->nelem - 4 and $xi > 1) { 
# BUG: this wont work with order < 5
	    return $infinity;
	}
	$KNMOD[$i] = $xi;
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
    if (!$config{patterns}) {
	for (my $i = 1; $i < $#{$s}; $i++) {
	    if (n0($s->[$i]) == 0) { 
		warn "Warning[$nline]: [$s->[$i]] unknown, skipping sentence.\n"
		    if not $config{optimize};
		return (0, 0);
	    }
	}
    }
    for (my $i = 1; $i < $#{$s}; $i++) {
	$nword++;
	my $b = bits($s, $i, $config{ngram});
	$nbits += $b;
	if ($config{verbose}) {
	    print $s->[$i];
	    for (my $n = 1; $n < $config{ngram}; $n++) {
		printf "\t%.4f", bits($s, $i, $n);
	    }
	    printf "\t%.4f\n", $b;
	}
# The unknown word hack
  	if ($config{mincnt} and n0($s->[$i]) < $config{mincnt}) {
  	    $s->[$i] = '<UNK>';
  	}
    }
    print "\n" if $config{verbose};
    return ($nbits, $nword);
}

sub read_corpus {
    my ($file) = @_;
    readfile($file, sub {
	return unless /\S/;
	my @s = ('<S>', gtokenize($_), '</S>');
	push @corpus, \@s;
    });
}

sub bits {
    my ($s, $i, $n) = @_;
    $n = 5 if not defined $n;
    if ($config{smoothing} =~ /^kn/) {
	my $pb = kn($s, $i, $n);
	die "Kneser-Ney returned $pb" if $pb <= 0;
	return -log2($pb);
    }
    if ($n == 1) {
	my $g = n0($s->[$i]);

	# BUG: Any word below 200 count is excluded from the google
	# data.  We will give such words a count of 100.  Note that
	# this may happen because of tokenization errors, and it is
	# not probabilistically sound.

	# die "Unknown word [$s->[$i]]" if $g == 0;
	$g = 100 if $g == 0;

	my $n_ = n0('_');
	warn "p0($s->[$i],$i,$n): $g/$n_=".($g/$n_)."\n" if $config{debug};

	return log2($n_ / $g);
    } elsif ($n > $i + 1) {
	return bits($s, $i, $i+1);
    }

    my $a = join(' ', @{$s}[($i-$n+1) .. ($i-1)]); # a = n-1 word prefix
    my $ga = n0($a);	# ga = count of a
    my $x = bits($s, $i, $n-1);	# x = lower order model bits
    if ($ga == 0) { 
#	    warn "Warning: Zero a-count[$a]\n" 
#		if $n <= 3 and $a =~ /[^\w ]/;		# check what is going on with punctuation
	return $x;		# return lower order result
    }
    my $px = exp2(-$x);		# px = lower order model probability
    my $b = $a . " $s->[$i]";	# b = all n words
    my $gb = n0($b);	# gb = count of b
    my $c = $a . " _";	# c = all ngrams that start with a
    my $gc = n0($c);	# gc = count of c
    my $missing_count = $ga - $gc; # missing_count = occurances of (a) 
    #   that are missing from ngram data
    die "Negative missing count" if $missing_count < 0;
    my $pb;

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.  In fact the only three cases are [!], [?], [;].  I am
# going to add [.] to the list for the same reason even though there
# seem to be bigrams starting with [.].

    # BUG: Note that the following will fail if we try to score the
    # </S> tokens.

    my $bad_context = 0;
    for my $tok (@{$s}[($i-$n+1) .. ($i-1)]) {
	if ($tok =~ /^[;!?.]$/) {
	    $bad_context = 1; last;
	}
    }

    if ($bad_context) {
	$pb = $px;

    } elsif ($config{smoothing} =~ /^mc/) {
	warn "mc($s->[$i],$i,$n): ML=$gb/$gc=".($gc>0?$gb/$gc:0)."\n" if $config{debug};
	my $extra = $D[$n] * $missing_count;
	$extra = 1 if $extra == 0;
	warn "mc($s->[$i],$i,$n): mc=$missing_count x D=$D[$n] = $extra\n" if $config{debug};
	if ($config{smoothing} eq 'mcmod') {
	    $px = mc2($s, $i, $n-1);
	}
	$pb = ($gb + $px * $extra)
	    / ($gc + $extra);
	warn "mc($s->[$i],$i,$n): (nxy=$gb + px=$px * $extra) / (nx_=$gc + $extra) = $pb\n" if $config{debug};
    } elsif ($config{smoothing} eq 'baseline') {
	$pb = $C[$n] * $px;
	$pb += (1 - $C[$n]) * ($gb / $gc)
	    if $gc > 0;
    } elsif ($config{smoothing} eq 'cdiscount') {
	my $D = 40 * $E[$n];
 	if ($ga == 0) {
 	    $pb = $px;
 	} elsif ($gb == 0) {
 	    $pb = $px * ($D * n1($c) + $missing_count) / $ga;
 	} else {
 	    $gb -= $D;
 	    $pb = $gb / $ga;
 	    $pb += $px * ($D * n1($c) + $missing_count) / $ga;
 	}

# The following gives horrible results
# It seems important to take mc into account

# 	if ($gc == 0) {
# 	    $pb = $px;
# 	} elsif ($gb == 0) {
# 	    $pb = $px * ($D * n1($c)) / $gc;
# 	} else {
#  	    $gb -= $D;
#  	    $pb = $gb / $gc;
#  	    $pb += $px * ($D * n1($c)) / $gc;
# 	}

    } elsif ($config{smoothing} eq 'wbdiscount') {

# This one gave 8.5773 on the 1k data:
# 	if ($missing_count > 0) {
# 	    $pb = (1 - $missing_count/$ga) * ($gb/$ga) + ($missing_count/$ga) * $px;
# 	} else {
# 	    $pb = (1 - 1/$ga) * ($gb/$ga) + (1/$ga) * $px;
# 	}

# This one gives 11.9696
# 	my $n1c = n1($c);
# 	my $bow = $n1c/($n1c+$ga);
# 	$bow = 1 if $bow == 0;
# 	$pb = (1-$bow) * ($gb/$ga) + $bow * $px;

# This one gives 7.9996: it is equivalent to mc if we ignore n1($c)
 	if ($gc == 0) {
 	    $pb = $px;
 	} else {
 	    my $n1c = n1($c) + 7*($ga - $gc);
 	    my $bow = $n1c/($n1c+$gc);
 	    die "WB: bow=0" if $bow == 0;
 	    $pb = (1-$bow) * ($gb/$gc) + $bow * $px;
 	}

    } else {
	die "Unknown smoothing $config{smoothing}";
    }
    return ($pb > 0 ? -log2($pb) : $infinity);
}


sub kn {
    my ($s, $i, $n) = @_;
    my ($nxy, $x, $nx_, $n1x_, $kn, $D);
    my ($nx, $mc);
    $n = 5 if not defined $n;
    my $y = $s->[$i];

    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {

 	# BUG: no need for 0 order model, if the word is known, this
 	# is exactly equivalent to k/N.  If the word is not known it
 	# is cheating.
	
 	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $g = n0($y);

	# BUG: Any word below 200 count is excluded from the google
	# data.  We will give such words a count of 100.  Note that
	# this may happen because of tokenization errors, and it is
	# not probabilistically sound.

	# die "Unknown word [$y]" if $g == 0;
	$g = 100 if $g == 0;

	my $n_ = n0('_');
	warn "kn($y,$i,$n): $g/$n_=".($g/$n_)."\n" if $config{debug};

	return $g / $n_;
	
    } elsif ($n > 5) {
	return kn($s, $i, 5);
    } elsif ($n > $i + 1) {
	return kn($s, $i, $i+1);
    }    

    die if $n <= 1;

    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $nx = n0($x);
    $nx_ = n0($x . ' _');
    $mc = $nx - $nx_;
    if ($nx == 0) {		#DBG - this should be nx_
	warn("kn($y,$i,$n): nx=0 for [$x]\n") if $config{debug};
	return kn($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn($y,$i,$n): nx=$nx nx_=$nx_ bad [$x]\n") if $config{debug};
	return kn($s, $i, $n-1);
    }	

    my $kn0 = ($config{smoothing} eq 'kn') ? kn0($s, $i, $n-1) : 
	($config{smoothing} eq 'knmc') ? kn1($s, $i, $n-1) :
	($config{smoothing} eq 'knmod') ? kn2($s, $i, $n-1) :
	die "KN: Unknown smoothing $config{smoothing}";
    warn("kn($y,$i,$n): kn0 = $kn0\n") if $config{debug};

    $nxy = n0("$x $y");
    warn("kn($y,$i,$n): ML: nxy = $nxy / nx = $nx => ".($nxy/$nx)."\n") if $config{debug};
    warn("kn($y,$i,$n): nx=$nx nx_=$nx_ mc=$mc\n") if $config{debug};
    $D = $MINNGRAM * 
	(($config{smoothing} eq 'kn') ? $KN[2*$n-4] :
	 ($config{smoothing} eq 'knmc') ? $KNMC[2*$n-4] :
	 ($config{smoothing} eq 'knmod') ? $KNMOD[$n+6] :
	 die "Bad smoothing [$config{smoothing}]");
    warn("kn($y,$i,$n): D = $D\n") if $config{debug};
    $nxy -= $D if $nxy > 0;

    $n1x_ = n1("$x _");
    warn("kn($y,$i,$n): n1x_ = $n1x_\n") if $config{debug};
    warn("kn($y,$i,$n): x = [$x] nx = $nx\n") if $config{debug};
    warn("kn($y,$i,$n): mc = $mc\n") if $config{debug};

    # Original formulation:
    #$kn = $nxy / $nx + (($mc + $n1x_ * $D) / $nx) * $kn0;

    # Using nx_ for denominator instead of nx, makes it worse.
    #if ($nx_ == 0) { return kn($s, $i, $n-1); }
    #$kn = $nxy / $nx_ + $kn0 * $D * $n1x_ / $nx_;

    # mc formulation
    #my $extra = $KNMOD[$n+2] * $mc;  #7.8232
    #7.82316150800905 <= [0.34999882 0.058431153  0.3384098 0.71771877  3.1924666  3.0841902   3.838957  4.1892848]

    #my $extra = $KNMOD[$n+2] * $nx_; #7.90
    #7.90320385739 <= 0.262751247301 0.048686713198 0.27605715289 0.568193284213 0.088951456544 0.715625 1.73570024289 2.16562605758 ssize=-0.0123526471100327

    #my $extra = $KNMOD[$n+2] * 100 * $n1x_; #7.8341
    #7.8340994726598 <= [0.32673847 0.055480826 0.32197277 0.66516695  1.9684104  2.5096591  3.3671446  3.7403359]

    my $extra = $KNMOD[$n+2] * $mc;
    $extra = 1 if $extra == 0;

    #orig: $kn = $nxy / $nx + (($mc + $n1x_ * $D) / $nx) * $kn0;
    #modf: $kn = ($nxy  + $kn0 * $extra) / ($nx_ + $extra);
    #eqiv: $kn = ($nxy + $kn0 * ($extra + $n1x_ * $D)) / ($nx + $extra - $mc);

    $kn = ($nxy + $kn0 * ($extra + $n1x_ * $D)) / ($nx_ + $extra);

    die "nxy=$nxy kn0=$kn0 mc=$mc nx_=$nx_ nx=$nx kn=0" if $kn == 0;

    warn("kn($y,$i,$n): kn = $kn\n") if $config{debug};
    if ($config{patterns}) { 
	# if outputting patterns do the lower order models as well
	kn($s, $i, $n-1); 
    }
    return $kn;
}

sub kn0 {
    my ($s, $i, $n) = @_;
    my ($x, $n1_xy, $n1_x_, $D, $kn0);
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {
	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1('_ ' . $s->[$i]);
	if ($n1_y == 0) {
	    warn("kn0($s->[$i],$i,$n): n1_y($s->[$i]) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _');
	warn("kn0($s->[$i],$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }
    die if $n <= 1;
    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $n1_x_ = n1("_ $x _");
    if ($n1_x_ == 0) {
	warn("kn0($s->[$i],$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return kn0($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
    { # bad context
	warn("kn0($s->[$i],$i,$n): n1_x_ = $n1_x_ bad x=[$x]\n") if $config{debug};
	return kn0($s, $i, $n-1);
    }	
    die "KN0 zero denominator [$x]" if $n1_x_ == 0;
    my $px = kn0($s, $i, $n-1);

    $n1_xy = n1("_ $x $s->[$i]");
    warn("kn0($s->[$i],$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};
    $D = $KN[2*$n-3];
    warn("kn0($s->[$i],$i,$n): D = $D\n") if $config{debug};
    $n1_xy -= $D if $n1_xy > 0;
    my $n2_x_ = n2("_ $x _");
    warn("kn0($s->[$i],$i,$n): n2_x_ = $n2_x_\n") if $config{debug};
#    $kn0 = $n1_xy / $n1_x_ + ($n1x_ * $D / $n1_x_) * kn0($s, $i, $n-1);	#BUGGY
    $kn0 = $n1_xy / $n1_x_ + ($n2_x_ * $D / $n1_x_) * $px;

    # new formula: does not really work
#     my $n1x = myngram($x);
#     my $mc = $n1x - $n1_x_;
#     $kn0 = $n1_xy / $n1x + (($mc + $n1x_ * $D) / $n1x) * kn0($s, $i, $n-1);

    warn("kn0($s->[$i],$i,$n): kn0 = $kn0\n") if $config{debug};
    return $kn0;
}

sub kn1 {
    my ($s, $i, $n) = @_;
    warn("kn1($s->[$i],$i,$n): hello\n") if $config{debug};
    if ($n <= 0) {
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n0y = n0($s->[$i]);
	my $n0_y = n0("_ $s->[$i]");
	my $mc = $n0y - $n0_y;
	$mc = 1 if $mc == 0;	# BUG
	my $n0_ = n0('_');
	my $n0__ = n0('_ _');
	my $MC = $n0_ - $n0__;
	return $mc/$MC;
    }
    die if $n <= 1;
    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $n0x_ = n0("$x _");
    my $n0_x_ = n0("_ $x _");
    my $MC = $n0x_ - $n0_x_;
    
    if ($n0x_ == 0) {
	return kn1($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
    { # bad context
	return kn1($s, $i, $n-1);
    }	
    elsif ($MC == 0) {
	# BUG: don't know what to do here exactly
	# warn "KN1 zero denominator [$x]";
	return kn1($s, $i, $n-1);
    }
    
    my $n0xy = n0("$x $s->[$i]");
    my $n0_xy = n0("_ $x $s->[$i]");
    my $mc = $n0xy - $n0_xy;
    my $D = $KNMC[2*$n-3];
    $mc -= $D if $mc > 0;
    my $n1x_ = n1("$x _");
    
    my $kn1 = $mc / $MC + ($n1x_ * $D / $MC) * kn1($s, $i, $n-1);

    warn("kn1($s->[$i],$i,$n): kn1 = $kn1\n") if $config{debug};
    return $kn1;
}

sub kn2 {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = ($KNMOD[0] + 1)/($KNMOD[0] + 40);
    if ($n <= 0) {
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1mod("_ $y", $coef);
	if ($n1_y == 0) {
	    warn("kn2($y,$i,$n): n1_y($y) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1mod('_ _', $coef);

	# OK, let's not use the extension here, gives better results:
	$n1_y = (n1("_ $y") or 1);
	$n1__ = n1("_ _");

	warn("kn2($y,$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }

    my $p0 = kn2($s, $i, $n-1);

    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $n1_x_ = n1mod("_ $x _", $coef);
    if ($n1_x_ == 0) {
	warn("kn2($y,$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return $p0;
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn2($y,$i,$n): n1_x_ = $n1_x_ bad [$x]\n") if $config{debug};
	return $p0;
    }	

    my $n1_xy = n1mod("_ $x $y", $coef);
    die "n1_xy=$n1_xy" if ($n1_xy > 0 and $n1_xy < 1);
    warn("kn2($y,$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};
#    my $D = $KNMOD[$n+6];
#    $D = 1 if $D > 1;
    my $D = 1;			# optimization gives 1.
    warn("kn2($y,$i,$n): D = $D\n") if $config{debug};
    $n1_xy -= $D if $n1_xy > 0;

    my $n1x_ = n1("$x _");
    warn("kn2($y,$i,$n): n1x_ = $n1x_\n") if $config{debug};

    my $C = $KNMOD[$n-1];

    my $n1modx_ = n1mod("$x _");
    my $kn2 = ($n1_xy + $p0 * ($C * $n1modx_ + $D * $n1x_)) / ($C * $n1modx_ + $n1_x_);

#7.8232 my $kn2 = ($n1_xy + $p0 * ($C * $n1_x_ + $D * $n1x_)) / ($C * $n1_x_ + $n1_x_);
#7.82316150800905 <= [0.34999882 0.058431153  0.3384098 0.71771877  3.1924666  3.0841902   3.838957  4.1892848]

#7.8072 my $kn2 = ($n1_xy + $p0 * ($C * $n1x_ + $D * $n1x_)) / ($C * $n1x_ + $n1_x_); #*run*
#7.80719881238559 <= [ 1.1273339  2.5163229  4.8874927  6.8448645  3.1924666  2.9414249  3.5130246  3.9000424]

#    my $n1modx_ = n1mod("$x _");# *run4*
#7.8031 my $kn2 = ($n1_xy + $p0 * ($C * $n1modx_ + $D * $n1x_)) / ($C * $n1modx_ + $n1_x_);	# *run4*
#7.80310938385236 <= [ 1.4369137  1.3023956  2.0387025  2.4608132  3.1924666  2.9146141  3.4571246  3.9006402]

#    my $n1_x = n1("_ $x");	# *run2*
#7.8162 my $kn2 = ($n1_xy + $p0 * ($C * $n1_x + $D * $n1x_)) / ($C * $n1_x + $n1_x_); #*run2*
#7.81622506825492 <= [0.56228426   1.422944  3.5268414  5.2919063  3.1924666  3.0078189  3.6986144  4.2018432]

#    my $n1mod_x = n1mod("_ $x");	# *run3*
#7.8149 my $kn2 = ($n1_xy + $p0 * ($C * $n1mod_x + $D * $n1x_)) / ($C * $n1mod_x + $n1_x_); #*run3*
#7.81492091713032 <= [0.59317426 0.63644359  1.3401929  1.8563427  3.1924666  2.9976826  3.6738253  4.2409812]

    warn("kn2($y,$i,$n): kn2 = $kn2\n") if $config{debug};
    return $kn2;
}

sub n1mod {
    my ($pat, $coef) = @_;
    if (defined $coef) {
	die "Bad n1mod coefficient $coef" if $coef < 1/40;
    } else {
	$coef = 1/40;
    }
    my $n1 = n1($pat);
    my $n0 = n0($pat);
    $pat =~ s/\s*_\s*//;
    my $n0x = n0($pat);
    return $n1 + $coef * ($n0x - $n0);
}

sub mc2 {
    my ($s, $i, $n) = @_;
    if ($n <= 0) {
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1mod("_ $s->[$i]");
	if ($n1_y == 0) {
	    warn("mc2($s->[$i],$i,$n): n1_y($s->[$i]) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1mod('_ _');

	# BUG: OK, let's not use the extension here:
	$n1_y = (n1("_ $s->[$i]") or 1);
	$n1__ = n1("_ _");
	#$n1_y = n0($s->[$i]);
	#$n1__ = n0('_');

	warn("mc2($s->[$i],$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }

    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $px = mc2($s, $i, $n-1);

    my $n1_x_ = n1mod("_ $x _");
    if ($n1_x_ == 0) {
	warn("mc2($s->[$i],$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return $px;
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
    { # bad context
	warn("mc2($s->[$i],$i,$n): n1_x_ = $n1_x_ bad [$x]\n") if $config{debug};
	return $px;
    }	
    my $n1_xy = n1mod("_ $x $s->[$i]");
    die "n1_xy=$n1_xy" if ($n1_xy > 0 and $n1_xy < 1);
    warn("mc2($s->[$i],$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};
    #my $D = $KNMOD[2*$n-2];
    warn("mc2($s->[$i],$i,$n): D = 1\n") if $config{debug};
    $n1_xy -= 1 if $n1_xy > 0;

    my $n1x_ = n1("$x _");
    warn("mc2($s->[$i],$i,$n): n1x_ = $n1x_\n") if $config{debug};

    my $mc2 = $n1_xy / $n1_x_ + $px * $n1x_ / $n1_x_;

    warn("mc2($s->[$i],$i,$n): mc2 = $mc2\n") if $config{debug};
    return $mc2;
}

sub dhc {
    my ($f, $xpdl) = @_;
    my @x = list $xpdl;
    my (@u, @v, @xv);
    my ($fx, $fxv);
    my ($vi, $vvec);
    my $vr;
    my ($i, $iter, $maxiter);
    my $NDIM = scalar(@x);
    my $THRESHOLD = 1e-4;
    my $INIT_SIZE = 0.1;
    my $count = 0;

    for($i=0; $i<$NDIM; $i++) {
	$u[$i] = 0; 
	$v[$i] = 0;
    }
    $vi = -1; $vvec = 1;
    $vr = -$INIT_SIZE;
    $fx = &$f(pdl(@x));
    $fxv = 1E10;

    printf("%d. %.12g <= ", ++$count, $fx);
    for($i=0; $i<$NDIM; $i++) { printf("%.12g ", $x[$i]); }; 
    printf("ssize=$vr\n");

    while(abs($vr) >= $THRESHOLD) {
	$maxiter = ((abs($vr) < 2*$THRESHOLD) ? 2*$NDIM : 2);
	$iter = 0;
	while(($fxv >= $fx) && ($iter < $maxiter)) {
	    if($iter == 0) { for($i=0; $i<$NDIM; $i++) { $xv[$i] = $x[$i]; } }
	    else { $xv[$vi] -= $vr; }
	    if($vvec) { $vvec = 0; }
	    $vr = -$vr;
	    if($vr > 0) { $vi = (($vi+1) % $NDIM); }
	    $xv[$vi] += $vr;
	    $fxv = &$f(pdl(@xv));
	    $iter++;
	}
	if($fxv >= $fx) {
	    $fxv = 1E10;
	    $vr /= 2;
	} else {
	    $fx = $fxv; printf("%d. %.12g <= ", ++$count, $fx);
	    for($i=0; $i<$NDIM; $i++) { $x[$i] = $xv[$i]; printf("%.12g ", $x[$i]); }
	    printf("ssize=$vr\n");
	    if($iter == 0) {
		if($vvec) {
		    for($i=0; $i<$NDIM; $i++) {
			$u[$i] += $v[$i]; $v[$i] *= 2; $xv[$i] += $v[$i];
		    }
		    $vr *= 2;
		} else {
		    $u[$vi] += $vr; $vr *= 2; $xv[$vi] += $vr;
		}
		$fxv = &$f(pdl(@xv));
	    } else {
		for($i=0; $i<$NDIM; $i++) { $xv[$i] += $u[$i]; }
		$xv[$vi] += $vr;
		$fxv = &$f(pdl(@xv));
		if($fxv >= $fx) {
		    for($i=0; $i<$NDIM; $i++) { $u[$i] = 0; $xv[$i] = $x[$i]; }
		    $u[$vi] = $vr; $vr *= 2;
		    $xv[$vi] += $vr; $fxv = &$f(pdl(@xv));
		} else {
		    for($i=0; $i<$NDIM; $i++) { $x[$i] = $xv[$i]; } $fx = $fxv;
		    $u[$vi] += $vr;
		    for($i=0; $i<$NDIM; $i++) { $v[$i] = 2*$u[$i]; } $vvec = 1;
		    for($i=0; $i<$NDIM; $i++) { $xv[$i] += $v[$i]; } $fxv = &$f(pdl(@xv));
		    for($vr=0,$i=0; $i<$NDIM; $i++) { $vr += $v[$i]*$v[$i]; }
		    $vr = sqrt($vr);
		}
	    }
	}
    }
    return (&$f(pdl(@x)), pdl(@x));
}

sub read_counts {
    my $path = shift;
    readfile($path, sub {
	my ($pat, $n0, $n1, $n2) = @_;
	die if not defined $pat;
	die if not defined $n0;
	$n0{$pat} = $n0;
	$n1{$pat} = $n1 if defined $n1;
	$n2{$pat} = $n2 if defined $n2;
	print STDERR '.' unless $. % 100000;
    }, "\t");
    return $n0{_};
}

sub n0 {
    my $pat = shift;
    if ($config{patterns}) {
	print "$pat\n" if 0 == $myngram{$pat}++;
    } 
    if (defined $n0{$pat}) {
	return $n0{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return $config{patterns} ? 1 : 0;
    }
}

sub n1 {
    my $pat = shift;
    if ($config{patterns}) {
	print "$pat\n" if 0 == $myngram{$pat}++;
    } 
    if (defined $n1{$pat}) {
	return $n1{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return $config{patterns} ? 1 : 0;
    }
}

sub n2 {
    my $pat = shift;
    if ($config{patterns}) {
	print "$pat\n" if 0 == $myngram{$pat}++;
    } 
    if (defined $n2{$pat}) {
	return $n2{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return $config{patterns} ? 1 : 0;
    }
}

1;

__END__

=head1 NAME

model.pl - Optimize and test language models

=head1 SYNOPSIS

model.pl [options] [file ...]

 Options:
   -counts file		read counts from file
   -debug		detailed debug output
   -dhc			use DHC for optimization
   -mincnt n		use <UNK> for words in context with count < n
   -ngram n		use ngram order n, default=5
   -optimize		perform parameter optimization
   -patterns		output ngram patterns
   -random		start optimization at random point
   -simplex		use Simplex for optimization
   -smoothing str	use smoothing method "str", default=knmod
   -string str		model last word of "str"
   -verbose		output bits for each word and each ngram order
   -verify str		verify probability=1 for last position of "str"		
   -zeroes		start optimization at zero

B<This program> will read the given input file(s) and compute cross
entropy for a language model specified by -smoothing and the counts
given in the -counts file.  The -optimize and related options perform
parameter optimization for the given -smoothing option.  The -patterns
option outputs the patterns that the -counts file will need to contain
in order to compute cross entropy.

The valid smoothing options are: baseline, kn, knmc, knmod, mc, mcmod,
cdiscount, wbdiscount.  Default is knmod, which gives 7.85 bits per
word on the brown corpus.

=cut

=begin comment

=head1 OPTIONS

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> will read the given input file(s) and do someting
useful with the contents thereof.

=head1 NOTES

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

