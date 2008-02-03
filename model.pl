#!/usr/bin/perl -w
warn q{$Id: model.pl,v 3.24 2008/02/03 07:27:58 dyuret Exp dyuret $ } ."\n";

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

# baseline: 8.20830869973577 if we use C(abc)/C(ab*)
my @C = (undef, undef, 0.12264181, 0.48531058, 0.73285371, 0.8485226);

# baseline2: 8.36994337610673 if we use C(abc)/C(ab)
# my @C = (undef, undef, 0.11574071, 0.48960756, 0.77180766, 0.89962748);

# cdiscount: 8.74236257742229 using C(ab*) for denominator
#my @E = (undef, undef, 0.94914893, 0.99209517, 0.99951819, 0.99995625);
# 3-gram does best: 8.52998951451499 (1k=8.46934122015906)
my @E = (undef, undef, 0.94934424, 0.99209517);

# cdiscount2: 8.23573621678019 (1k) using C(ab) for denominator
# my @E = (undef, undef, 0.58301752, 0.79111861, 0.92800359, 0.9840277);

# kn: 8.46259231435588 (1k=8.40261215945863) Using unmodified n1 counts.
# my @KN = (0.98520067, 0.8046389, 0.99577209, 0.91786366, 0.99933982, 0.98042218, 0.99996885);
# ngram=4 gives better results 8.40070320128183 (1k=8.33130330762054)
my @KN = (0.98520067, 0.80424828, 0.99577209, 0.93954335, 0.99953513);

# kn1: 8.25961028092419 (1k=8.19656987210192) Using modified n1 counts.
# my @KN1 = (0.038962617, 0.9851726, 0.97914173, 0.99581479, 0.99992424, 0.99931418, 0.99985808, 0.99996568);
# ngram=4 gives better results 8.22729782692671 (1k=8.15830220985079)
my @KN1 = (0.095564851, 0.9851726, 0.99995823, 0.99581479, 0.99992424, 0.99949497);

# mc: 8.0477965326711 (1k=7.99946312536937)
my @D = (undef, undef, 7.8440119, 6.6305746, 6.9277921, 7.3675802); 

# kn3: 7.96264920962347 (1k=7.90877748764095) Use mc with kn backoff using modified n1 counts.
my @KN3 = (0.020662374, 3.1924203, 0.95571439, 3.2104425, 0.99992424, 4.8814022, 0.9999869, 6.8346343);

# kn4: 7.86311246273543 (1k=7.81034468881) Use dirichlet (A=K*mc) with no D discount with both kn and its backoff
# Reported in paper
my @KN4 = (0.119330390707, 3.192501, 0.156310135, 2.8096546125, 0.205351523914, 3.23767800788, 0.238121750327, 3.58515159429);

# knmod: 7.85546292193106 (1k=7.80310938385546)
# Best results
my @KNMOD = 			# Use mc with both kn and its backoff
    (
     # D = 1 for kn2, D=0 for kn

     # coef=(A+1)/(A+40) for kn2:
     #1.4368796,

     # coef for kn2:
     .05880943795777517957,

     # C2..C5 (x mc) for kn:	# C2..C4 (x n1modx_) for kn2:
     3.192501, 			1.3022398,
     2.9146123, 		2.0386943, 
     3.4570501, 		2.4608205,
     3.9005968,
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
     'kn1' => \&init_kn1,
     'kn3' => \&init_kn3,
     'kn4' => \&init_kn4,
     'knmod' => \&init_knmod,
     );

my %score_fn = 
    ('mc' => \&score_mc, 
     'baseline' => \&score_baseline, 
     'kn' => \&score_kn, 
     'cdiscount' => \&score_cdiscount,
     'kn1' => \&score_kn1,
     'kn3' => \&score_kn3,
     'kn4' => \&score_kn4,
     'knmod' => \&score_knmod,
     );

# Globals:

my $ZERO_WARNING;
my $LOWER_WARNING;
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
     smoothing => 'knmod', # { baseline, mc, kn, cdiscount, kn1, kn3, kn4, wbdiscount, knmod }
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

sub init_kn1 {
    my $init;
    my $ndims = 2 * $config{ngram} - 2;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KN1[0 .. $ndims-1]);
    return $init;
}

sub init_kn3 {
    my $init;
    my $ndims = 2 * $config{ngram} - 2;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KN3[0 .. $ndims-1]);
    return $init;
}

sub init_kn4 {
    my $init;
    my $ndims = 2 * $config{ngram} - 2;
    $init = $config{zeroes} ? ones($ndims) 
	: $config{random} ? random($ndims)
	: pdl(@KN4[0 .. $ndims-1]);
    return $init;
}

sub init_knmod {
    my $init;
    my $ndims = 2 * $config{ngram} - 2;
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
	my $xi = $x->at($i-2);
	if ($xi < 0 or $xi > 1) { return $infinity; }
	$E[$i] = $xi;
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

sub score_kn1 {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	if ($x->at($i) < 0 or ($i > 0 and $x->at($i) > 1)) {
	    return $infinity;
	}
	$KN1[$i] = $x->at($i);
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_kn3 {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	my $xi = $x->at($i);
	if ($xi < 0) {
	    return $infinity;
	} elsif ($i > 0 && $i%2 == 0 && $xi > 1) {
	    return $infinity;
	}
	$KN3[$i] = $xi;
    }
    my $bits = ngram();
    warn "score[$nscore]: $bits $x\n";
    return $bits;
}

sub score_kn4 {			
    my $x = shift;
    $nscore++;
    for (my $i = 0; $i < $x->nelem; $i++) {
	my $xi = $x->at($i);
	if ($xi < 0) {
	    return $infinity;
	} elsif ($i == 0 and ($xi < 1/40 or $xi > 1)) {
	    return $infinity;
	}
	$KN4[$i] = $xi;
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
	} elsif ($i == 0 and ($xi > 1 or $xi < 1/40)) {
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
	$pb = ($gb + $px * $extra)
	    / ($gc + $extra);
	warn "mc($s->[$i],$i,$n): (nxy=$gb + px=$px * $extra) / (nx_=$gc + $extra) = $pb\n" if $config{debug};
    } elsif ($config{smoothing} eq 'baseline') {
	$pb = $C[$n] * $px;
	$pb += (1 - $C[$n]) * ($gb / $gc)
	    if $gc > 0;
    } elsif ($config{smoothing} eq 'cdiscount') {
	if ($n > 3) {
	    warn "Warning: cdiscount does best with 3-grams, using lower order model\n"
		unless $LOWER_WARNING++;
	    $pb = $px;
	} else {
	    my $D = 40 * $E[$n];
#  	if ($ga == 0) {
#  	    $pb = $px;
# 	} else {
# 	    $pb = (($gb > 0) ? ($gb - $D) : 0) / $ga +
# 		$px * ($D * n1($c) + $missing_count) / $ga;
#  	}

# The following gives horrible results if n > 3.
# It seems important to take mc into account

	    if ($gc == 0) {
		$pb = $px;
	    } else {
		$pb = (($gb > 0) ? ($gb - $D) : 0) / $gc +
		    $px * ($D * n1($c)) / $gc;
	    }
	}
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
    if (($nx == 0) or
	(($nx_ == 0) and 
	 (($config{smoothing} eq 'kn') or
	  ($config{smoothing} eq 'kn1'))))
    {
	warn("kn($y,$i,$n): nx=0 for [$x]\n") if $config{debug};
	return kn($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn($y,$i,$n): nx=$nx nx_=$nx_ bad [$x]\n") if $config{debug};
	return kn($s, $i, $n-1);
    } elsif ($n == 5 and ($config{smoothing} eq 'kn' or $config{smoothing} eq 'kn1')) {
	warn("Warning: $config{smoothing} does best with 4-grams, using lower order model\n")
	    unless $LOWER_WARNING++;
	return kn($s, $i, $n-1);
    }

    my $kn0 = ($config{smoothing} eq 'kn') ? kn0($s, $i, $n-1) : 
	($config{smoothing} eq 'kn1') ? kn1($s, $i, $n-1) :
	($config{smoothing} eq 'knmod') ? kn2($s, $i, $n-1) :
	($config{smoothing} eq 'kn3') ? kn3($s, $i, $n-1) :
	($config{smoothing} eq 'kn4') ? kn4($s, $i, $n-1) :
	die "KN: Unknown smoothing $config{smoothing}";
    warn("kn($y,$i,$n): kn0 = $kn0\n") if $config{debug};

    $nxy = n0("$x $y");
    warn("kn($y,$i,$n): ML: nxy = $nxy / nx = $nx => ".($nxy/$nx)."\n") if $config{debug};
    warn("kn($y,$i,$n): nx=$nx nx_=$nx_ mc=$mc\n") if $config{debug};

    $n1x_ = n1("$x _");
    warn("kn($y,$i,$n): n1x_ = $n1x_\n") if $config{debug};
    warn("kn($y,$i,$n): x = [$x] nx = $nx\n") if $config{debug};
    warn("kn($y,$i,$n): mc = $mc\n") if $config{debug};

    if ($config{smoothing} eq 'kn') {
	$D = $MINNGRAM * $KN[2*$n-4];
	warn("kn($y,$i,$n): D = $D\n") if $config{debug};
	$kn = ($nxy == 0 ? 0 : $nxy - $D) / $nx_ +
	    $kn0 * $n1x_ * $D / $nx_;
	
    } elsif ($config{smoothing} eq 'kn1') {
	$D = $MINNGRAM * $KN1[2*$n-3];
	warn("kn($y,$i,$n): D = $D\n") if $config{debug};
	$kn = ($nxy == 0 ? 0 : $nxy - $D) / $nx_ +
	    $kn0 * $n1x_ * $D / $nx_;
	
    } elsif ($config{smoothing} eq 'kn3') {
	my $extra = $KN3[2*$n-3] * $mc;
	$extra = 1 if $extra == 0;
	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);

    } elsif ($config{smoothing} eq 'kn4') {
	my $extra = $KN4[2*$n-3] * $mc;
	$extra = 1 if $extra == 0;
	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);

    } elsif ($config{smoothing} eq 'knmod') {
	my $extra = $KNMOD[2*$n-3] * $mc;
	$extra = 1 if $extra == 0;
	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);

    } else {
	die "$config{smoothing} is not implemented.";
    }

    die "nxy=$nxy kn0=$kn0 mc=$mc nx_=$nx_ nx=$nx kn=0" if $kn == 0;
    warn("kn($y,$i,$n): kn = $kn\n") if $config{debug};
    if ($config{patterns}) { 
	# if outputting patterns do the lower order models as well
	kn($s, $i, $n-1); 
    }
    return $kn;
}

# kn experiments:
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

	#orig: $kn = $nxy / $nx + (($mc + $n1x_ * $D) / $nx) * $kn0;
	#modf: $kn = ($nxy  + $kn0 * $extra) / ($nx_ + $extra);
	#equiv: $kn = ($nxy + $kn0 * ($extra + $n1x_ * $D)) / ($nx + $extra - $mc);
	#withD: $kn = ($nxy + $kn0 * ($extra + $n1x_ * $D)) / ($nx_ + $extra);


# KN0: backoff function that uses unmodified n1 counts:
# (kn) 	$D = $MINNGRAM * $KN[2*$n-4];
# 	$kn = ($nxy == 0 ? 0 : $nxy - $D) / $nx_ +
# 	    $kn0 * $n1x_ * $D / $nx_;
# (kn0) $D = $KN[2*$n-3];
#       $kn0 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
# 	    $px * $n2_x_ * $D / $n1_x_;

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
    # $n1_xy -= $D if $n1_xy > 0;
    my $n2_x_ = n2("_ $x _");
    warn("kn0($s->[$i],$i,$n): n2_x_ = $n2_x_\n") if $config{debug};
#    $kn0 = $n1_xy / $n1_x_ + ($n1x_ * $D / $n1_x_) * kn0($s, $i, $n-1);	#BUGGY
    $kn0 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
	$px * $n2_x_ * $D / $n1_x_;

    # new formula: does not really work
#     my $n1x = myngram($x);
#     my $mc = $n1x - $n1_x_;
#     $kn0 = $n1_xy / $n1x + (($mc + $n1x_ * $D) / $n1x) * kn0($s, $i, $n-1);

    warn("kn0($s->[$i],$i,$n): kn0 = $kn0\n") if $config{debug};
    return $kn0;
}


# KN1: backoff function that uses modified n1 counts:
# (kn)	$D = $MINNGRAM * $KN1[2*$n-3];
# 	$kn = ($nxy == 0 ? 0 : $nxy - $D) / $nx_ +
# 	    $kn0 * $n1x_ * $D / $nx_;
# (kn1)	$D = $KN1[2*$n-2];
# 	$kn1 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
# 	    $p0 * $n1x_ * $D / $n1_x_;

sub kn1 {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = ($KN1[0] + 1)/($KN1[0] + 40);
    my ($x, $n1_xy, $n1_x_, $D, $kn1);
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {
	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1('_ ' . $y);
	if ($n1_y == 0) {
	    warn("kn1($y,$i,$n): n1_y($y) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _');
	warn("kn1($y,$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }
    die if $n <= 1;

    my $p0 = kn1($s, $i, $n-1);

    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $n1_x_ = n1mod("_ $x _", $coef);
    if ($n1_x_ == 0) {
	warn("kn1($y,$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return $p0;
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn1($y,$i,$n): n1_x_ = $n1_x_ bad x=[$x]\n") if $config{debug};
	return $p0;
    }	
    die "KN1 zero denominator [$x]" if $n1_x_ == 0;

    $n1_xy = n1mod("_ $x $y", $coef);
    warn("kn1($y,$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};
    $D = $KN1[2*$n-2];
    warn("kn1($y,$i,$n): D = $D\n") if $config{debug};
    # $n1_xy -= $D if $n1_xy > 0;
    my $n1x_ = n1("$x _");
    $kn1 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
	$p0 * $n1x_ * $D / $n1_x_;

    warn("kn1($y,$i,$n): kn1 = $kn1\n") if $config{debug};
    return $kn1;
}

# KN3: identical to kn1, but kn uses missing-count discounting
# (kn)	my $extra = $KN3[2*$n-3] * $mc;
# 	$extra = 1 if $extra == 0;
# 	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);
# (kn1)	$D = $KN3[2*$n-2];
# 	$kn3 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
# 	    $p0 * $n1x_ * $D / $n1_x_;

sub kn3 {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = ($KN3[0] + 1)/($KN3[0] + 40);
    my ($x, $n1_xy, $n1_x_, $D, $kn3);
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {
	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1('_ ' . $y);
	if ($n1_y == 0) {
	    warn("kn3($y,$i,$n): n1_y($y) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _');
	warn("kn3($y,$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }
    die if $n <= 1;

    my $p0 = kn3($s, $i, $n-1);

    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $n1_x_ = n1mod("_ $x _", $coef);
    if ($n1_x_ == 0) {
	warn("kn3($y,$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return $p0;
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn3($y,$i,$n): n1_x_ = $n1_x_ bad x=[$x]\n") if $config{debug};
	return $p0;
    }	
    die "KN3 zero denominator [$x]" if $n1_x_ == 0;

    $n1_xy = n1mod("_ $x $y", $coef);
    warn("kn3($y,$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};
    $D = $KN3[2*$n-2];
    warn("kn3($y,$i,$n): D = $D\n") if $config{debug};
    # $n1_xy -= $D if $n1_xy > 0;
    my $n1x_ = n1("$x _");
    $kn3 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / $n1_x_ + 
	$p0 * $n1x_ * $D / $n1_x_;

    warn("kn3($y,$i,$n): kn3 = $kn3\n") if $config{debug};
    return $kn3;
}

# KN2: both kn and kn2 uses missing-count discounting
# (kn)	my $extra = $KNMOD[2*$n-3] * $mc;
# 	$extra = 1 if $extra == 0;
# 	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);
# (kn2)	$D = 1;
#	$C = $KNMOD[2*$n-2];
#       $kn2 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / ($C * $n1modx_ + $n1_x_) + 
#	       $p0 * ($C * $n1modx_ + $n1x_ * $D) / ($C * $n1modx_ + $n1_x_);

sub kn2 {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = $KNMOD[0];
    if ($n <= 0) {
	die "Bad n [$n]";
    } elsif ($n == 1) {
	# we do not use the modified n1 count here.
	my $n1_y = n1("_ $y", $coef);
	if ($n1_y == 0) {
	    warn("kn2($y,$i,$n): n1_y($y) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _', $coef);

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

    my $n1x_ = n1("$x _");
    warn("kn2($y,$i,$n): n1x_ = $n1x_\n") if $config{debug};

    my $D = 1;			# optimization gives 1.
    my $C = $KNMOD[2*$n-2];
    my $n1modx_ = n1mod("$x _");
    warn("kn2($y,$i,$n): D = $D C = $C n1modx_ = $n1modx_\n") if $config{debug};

    my $kn2 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / ($C * $n1modx_ + $n1_x_) + 
	$p0 * ($C * $n1modx_ + $n1x_ * $D) / ($C * $n1modx_ + $n1_x_);
    warn("kn2($y,$i,$n): kn2 = $kn2\n") if $config{debug};
    return $kn2;
}

# kn2 experiments:
#
# Note that KNMOD[0] in some old experiments was interpreted not as
# the coef itself but as coef = (1+A)/(40+A)
#
#    my $n1modx_ = n1mod("$x _");# *run4*
#7.8031 my $kn2 = ($n1_xy + $p0 * ($C * $n1modx_ + $D * $n1x_)) / ($C * $n1modx_ + $n1_x_);	# *run4*
#(D=1): 7.80310938385236 <= [ 1.4369137  1.3023956  2.0387025  2.4608132  3.1924666  2.9146141  3.4571246  3.9006402]
#(D=0): 7.8093857188333 <= [ 3.5378598   3.192501  3.0752446  2.8224962  4.4838955  3.2652199   5.410972  3.6069656]

#7.8072 my $kn2 = ($n1_xy + $p0 * ($C * $n1x_ + $D * $n1x_)) / ($C * $n1x_ + $n1_x_); #*run*
#7.80719881238559 <= [ 1.1273339  2.5163229  4.8874927  6.8448645  3.1924666  2.9414249  3.5130246  3.9000424]

#    my $n1mod_x = n1mod("_ $x");	# *run3*
#7.8149 my $kn2 = ($n1_xy + $p0 * ($C * $n1mod_x + $D * $n1x_)) / ($C * $n1mod_x + $n1_x_); #*run3*
#7.81492091713032 <= [0.59317426 0.63644359  1.3401929  1.8563427  3.1924666  2.9976826  3.6738253  4.2409812]

#    my $n1_x = n1("_ $x");	# *run2*
#7.8162 my $kn2 = ($n1_xy + $p0 * ($C * $n1_x + $D * $n1x_)) / ($C * $n1_x + $n1_x_); #*run2*
#7.81622506825492 <= [0.56228426   1.422944  3.5268414  5.2919063  3.1924666  3.0078189  3.6986144  4.2018432]

#7.8232 my $kn2 = ($n1_xy + $p0 * ($C * $n1_x_ + $D * $n1x_)) / ($C * $n1_x_ + $n1_x_);
#7.82316150800905 <= [0.34999882 0.058431153  0.3384098 0.71771877  3.1924666  3.0841902   3.838957  4.1892848]


# KN4: both kn and kn4 uses missing-count discounting
# (kn)	my $extra = $KN4[2*$n-3] * $mc;
# 	$extra = 1 if $extra == 0;
# 	$kn = ($nxy + $kn0 * $extra) / ($nx_ + $extra);
# (kn4)	$extra = $KN4[2*$n-2] * $mc;
# (D=1) $kn4 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / ($extra + $n1_x_) + 
#	       $p0 * ($extra + $n1x_ * $D) / ($extra + $n1_x_);
# (D=0)	$kn4 = $n1_xy / ($extra + $n1_x_) + 
# 		$p0 * $extra / ($extra + $n1_x_);

#7.8038 my $kn2 = ($n1_xy + $p0 * ($extra + $D * $n1x_)) / ($extra + $n1_x_);
#(D=1): 7.80377791076084 <= [0.059170783  3.1925007 0.058437701  2.9088102 0.082038069  3.4459117 0.094556554  3.9088292]
#(D=0):	7.81034498858537 [0.11966381   3.192501 0.15670076  2.8094593 0.20582305   3.237459 0.23886949  3.5852088]

sub kn4 {
    my ($s, $i, $n) = @_;
    my $y = $s->[$i];
    my $coef = $KN4[0];
    if ($n <= 0) {
	die "Bad n [$n]";
    } elsif ($n == 1) {
	#my $n1_y = n1mod("_ $y", $coef);
	my $n1_y = n1("_ $y", $coef);
	if ($n1_y == 0) {
	    warn("kn4($y,$i,$n): n1_y($y) == 0\n") if $config{debug};
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	#my $n1__ = n1mod('_ _', $coef);
	my $n1__ = n1('_ _', $coef);

	warn("kn4($y,$i,$n): n1_y=$n1_y / n1__=$n1__ => ".($n1_y/$n1__)."\n") if $config{debug};
	return $n1_y / $n1__;
    }

    my $p0 = kn4($s, $i, $n-1);

    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $n1_x_ = n1mod("_ $x _", $coef);
    if ($n1_x_ == 0) {
	warn("kn4($y,$i,$n): n1_x_ = $n1_x_\n") if $config{debug};
	return $p0;
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $y eq '</S>'))
    { # bad context
	warn("kn4($y,$i,$n): n1_x_ = $n1_x_ bad [$x]\n") if $config{debug};
	return $p0;
    }	

    my $n1_xy = n1mod("_ $x $y", $coef);
    die "n1_xy=$n1_xy" if ($n1_xy > 0 and $n1_xy < 1);
    warn("kn4($y,$i,$n): ML: n1_xy = $n1_xy / n1_x_ = $n1_x_ => ".($n1_xy/$n1_x_)."\n") if $config{debug};

    my $n1x_ = n1("$x _");
    warn("kn4($y,$i,$n): n1x_ = $n1x_\n") if $config{debug};

    # my $D = 1;			# optimization gives 1.
    # warn("kn4($y,$i,$n): D = $D\n") if $config{debug};

    my $mc = n0($x) - n0("$x _");
    my $extra = $KN4[2*$n-2] * $mc;
    $extra = 1 if $extra == 0;

#     my $kn4 = ($n1_xy == 0 ? 0 : $n1_xy - $D) / ($extra + $n1_x_) + 
# 	$p0 * ($extra + $n1x_ * $D) / ($extra + $n1_x_);
    my $kn4 = $n1_xy / ($extra + $n1_x_) + 
 	$p0 * $extra / ($extra + $n1_x_);
    warn("kn4($y,$i,$n): kn4 = $kn4\n") if $config{debug};
    return $kn4;
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

Typical usage:

  model.pl -patterns < text > patterns
  glookup < patterns > counts
  model.pl -counts counts < text

B<This program> will read the given input file(s) and compute cross
entropy for a language model specified by -smoothing and the counts
given in the -counts file.  The -optimize and related options perform
parameter optimization for the given -smoothing option.  The -patterns
option outputs the patterns that the -counts file will need to contain
in order to compute cross entropy.

The valid smoothing options are: baseline, kn, knmod, mc, cdiscount,
wbdiscount.  Default is knmod, which gives 7.85 bits per word on the
brown corpus.

=head1 OPTIONS
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

