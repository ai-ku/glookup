#!/usr/bin/perl -w
warn q{$Id: model.pl,v 3.10 2008/01/28 10:56:44 dyuret Exp dyuret $ } ."\n";

use strict;
use Getopt::Long;
use PDL;
use PDL::Opt::Simplex;
my $oldfh = select(STDERR); $| = 1; select($oldfh);
require 'gtokenize.pl';
my $log2 = log(2);
sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }
my $log10 = log(10);
#sub log10 { log($_[0])/$log10; } # defined in PDL
sub exp10 { exp($_[0]*$log10); }

# Command line options:

my $verbose = 0;
my $debug = 0;
my $optimize = 0;
my $dhc = 0;
my $simplex = 0;
my $patterns;			# output patterns
my $cachefile;
my $ngram = 5;
my $random = 0;
my $zeroes = 0;
# smoothing = { mc, baseline, kn, cdiscount }
my $smoothing = 'mc';
my @C = (undef, undef, 0.12264181, 0.48531058, 0.73285371, 0.8485226); # baseline: 8.20830869973577
my @D = (undef, undef, 7.8440119, 6.6305746, 6.9277921, 7.3675802); # mc: 8.0477965326711
my @KN = (0.83456664, 0.80501932, 0.81691654, 0.90655321, 0.97261754, 0.96917268, 0.97422316); # kn: 8.23974660099736
my @E = (undef, undef, 0.58301752, 0.79111861, 0.92800359, 0.9840277); # cdiscount: 

GetOptions('cache=s' => \$cachefile,
           'verbose' => \$verbose,
           'debug' => \$debug,
           'optimize' => \$optimize,
	   'dhc' => \$dhc,
	   'simplex' => \$simplex,
	   'patterns' => \$patterns,
           'random' => \$random,
           'zeroes' => \$zeroes,
	   'ngram=i' => \$ngram,
           'smoothing=s' => \$smoothing,
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
 	   'e2=f' => \$E[2],
 	   'e3=f' => \$E[3],
 	   'e4=f' => \$E[4],
 	   'e5=f' => \$E[5],
);

my %init_fn = ('mc' => \&init_mc, 'baseline' => \&init_baseline, 'kn' => \&init_kn, 'cdiscount' => \&init_cdiscount);
my %score_fn = ('mc' => \&score_mc, 'baseline' => \&score_baseline, 'kn' => \&score_kn, 'cdiscount' => \&score_cdiscount);

# Main:

die "Please specify -cache cntfile or -patterns to output patterns\n"
    if not $cachefile and not $patterns;
my $ZERO_WARNING;
my (%n0, %n1, %n2);
my $GTotal = $patterns ? 1024908267229 : cnt_init($cachefile);
warn "\ngtotal=$GTotal\n";
my $corpus = read_corpus();
warn "corpus=" . scalar(@$corpus) . " sentences\n";
my $nline = 0;
my $nscore = 0;
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
    if ($simplex) {
	my ($optimum, $ssize) = simplex($init, $initsize, $minsize,
					$maxiter, \&score, \&display);
    } else {			# dhc is the default
	my ($optimum, $final) = dhc(\&score, $init);
	warn "$optimum <= $final\n";
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

sub init_mc {
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

sub init_cdiscount {
    my $init;
    my $ndims = $ngram - 1;
    $init = $zeroes ? zeroes($ndims) 
	: $random ? random($ndims)
	: pdl(@E[2 .. $ngram]);
    return $init;
}

sub init_kn {
    my $init;
    my $ndims = 2 * $ngram - 3;
    $init = $zeroes ? ones($ndims) 
	: $random ? random($ndims)
	: pdl(@KN[0 .. $ndims-1]);
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

sub score_mc {			
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

sub score_cdiscount {			
    my $x = shift;
    $nscore++;
    for (my $i = 2; $i <= $ngram; $i++) {
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
	    if (n0($s->[$i]) == 0) { 
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
    while(<>) {
	return unless /\S/;
	my @s = ('<S>', gtokenize($_), '</S>');
	push @corpus, \@s;
    }
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
	my $g = n0($s->[$i]);

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

    } elsif ($smoothing eq 'mc') {
	if ($missing_count > 0) { # apply our smoothing formula
	    my $extra = $D[$n] * $missing_count;
	    $pb = ($gb + $px * $extra)
		/ ($gc + $extra);
	} elsif ($missing_count == 0) {
	    #warn "Warning: Zero missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	    $pb = ($gb + $px) / ($gc + 1);
	    
	} elsif ($missing_count < 0) {
		
# This is a bug in gngram.  This means there are more instances of
# words following "A" than there are instances of "A" by itself.
		
	    die "Error: Negative missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	    #$pb = ($gb + $px) / ($gc + 1);
	}
    } elsif ($smoothing eq 'baseline') {
	$pb = $C[$n] * $px;
	$pb += (1 - $C[$n]) * ($gb / $gc)
	    if $gc > 0;
    } elsif ($smoothing eq 'cdiscount') {
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

    } else {
	die "Unknown smoothing $smoothing";
    }
    return ($pb > 0 ? -log2($pb) : $infinity);
}


sub kn {
    my ($s, $i, $n) = @_;
    my ($nxy, $x, $nx_, $n1x_, $kn, $D);
    my ($nx, $mc);
    $n = 5 if not defined $n;
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {

 	# BUG: no need for 0 order model, if the word is known, this
 	# is exactly equivalent to k/N.  If the word is not known it
 	# is cheating.
	
 	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $g = n0($s->[$i]);

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

    die if $n <= 1;
    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $nx = n0($x);
    $nx_ = n0($x . ' _');
    $mc = $nx - $nx_;
    warn("kn($s->[$i],$i,$n): nx_ = $nx_\n") if $debug;
    if ($nx == 0) {
	return kn($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
    { # bad context
	return kn($s, $i, $n-1);
    }	
    $n1x_ = n1("$x _");
    warn("kn($s->[$i],$i,$n): n1x_ = $n1x_\n") if $debug;
    $nxy = n0("$x $s->[$i]");
    warn("kn($s->[$i],$i,$n): nxy = $nxy\n") if $debug;
    $D = $MINNGRAM * $KN[2*$n-4];
    warn("kn($s->[$i],$i,$n): D = $D\n") if $debug;
    $nxy -= $D if $nxy > 0;
    #$kn = $nxy / $nx_ + ($n1x_ * $D / $nx_) * kn0($s, $i, $n-1);

    warn("kn($s->[$i],$i,$n): x = [$x] nx = $nx\n") if $debug;
    warn("kn($s->[$i],$i,$n): mc = $mc\n") if $debug;
    my $kn0 = kn0($s, $i, $n-1);
    warn("kn($s->[$i],$i,$n): kn0 = $kn0\n") if $debug;

    $kn = $nxy / $nx + (($mc + $n1x_ * $D) / $nx) * $kn0;

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
    warn("kn0($s->[$i],$i,$n): hello\n") if $debug;
    my ($x, $n1_xy, $n1_x_, $D, $kn0);
    if ($n < 0) {
	die "Bad n [$n]";
    } elsif ($n == 0) {
	# return 1/n1('_');
	die "Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1('_ ' . $s->[$i]);
	if ($n1_y == 0) {
	    warn("kn0($s->[$i],$i,$n): n1_y($s->[$i]) == 0\n") if $debug;
	    # BUG: verify that this is exactly what would happen if we smoothed with a 0-order 1/V model
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _');
	warn("kn0($s->[$i],$i,$n): n1_y=$n1_y n1__=$n1__\n") if $debug;
	return $n1_y / $n1__;
    }
    die if $n <= 1;
    $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    $n1_x_ = n1("_ $x _");
    warn("kn0($s->[$i],$i,$n): n1_x_ = $n1_x_\n") if $debug;
    if ($n1_x_ == 0) {
	return kn0($s, $i, $n-1);
    } elsif ((" $x " =~ / [;!?.] /)
	     and not ($x =~ /^[^;!?.]*[;!?.]$/ and $s->[$i] eq '</S>'))
    { # bad context
	return kn0($s, $i, $n-1);
    }	
    die "KN0 zero denominator [$x]" if $n1_x_ == 0;
    $n1_xy = n1("_ $x $s->[$i]");
    warn("kn0($s->[$i],$i,$n): n1_xy = $n1_xy\n") if $debug;
    $D = $KN[2*$n-3];
    warn("kn0($s->[$i],$i,$n): D = $D\n") if $debug;
    $n1_xy -= $D if $n1_xy > 0;
    my $n2_x_ = n2("_ $x _");
#    $kn0 = $n1_xy / $n1_x_ + ($n1x_ * $D / $n1_x_) * kn0($s, $i, $n-1);	#BUGGY
    $kn0 = $n1_xy / $n1_x_ + ($n2_x_ * $D / $n1_x_) * kn0($s, $i, $n-1);

    # new formula: does not really work
#     my $n1x = myngram($x);
#     my $mc = $n1x - $n1_x_;
#     $kn0 = $n1_xy / $n1x + (($mc + $n1x_ * $D) / $n1x) * kn0($s, $i, $n-1);

    warn("kn0($s->[$i],$i,$n): kn0 = $kn0\n") if $debug;
    return $kn0;
}

# sub myngram {
#     my $str = shift;
#     if ($patterns) {
# 	if (not defined $myngram{$str}) {
# 	    $myngram{$str} = 1;
# 	    print "$str\n";
# 	}
# 	return 1;
#     } else {

# 	# BUG: temporary fix
# 	$str =~ s/^\s+//;
# 	$str =~ s/\s+$//;
# 	if ($str =~ /^<S>$/) {
# 	    return 95119665584;
# 	}
# 	# End of temporary fix

# 	return gngram($str);
#     }
# }

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
    for($i=0; $i<$NDIM; $i++) { printf("%.12g ", $x[$i]); }; printf("\n");

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
	    printf("\n");
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

sub cnt_init {
    my $path = shift;
    return if not defined $path;
    open(FP, $path) or die $!;
    while(<FP>) {
	chomp;
	my ($pat, $n0, $n1, $n2) = split(/\t/);
	die if not defined $pat;
	die if not defined $n0;
	$n0{$pat} = $n0;
	$n1{$pat} = $n1 if defined $n1;
	$n2{$pat} = $n2 if defined $n2;
	print STDERR '.' unless $. % 100000;
    }
    close(FP);
    return $n0{_};
}

sub n0 {
    my $pat = shift;
    if ($patterns) {
	print "$pat\n" if 0 == $myngram{$pat}++;
	return 1;
    } elsif (defined $n0{$pat}) {
	return $n0{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return 0;
    }
}

sub n1 {
    my $pat = shift;
    if ($patterns) {
	print "$pat\n" if 0 == $myngram{$pat}++;
	return 1;
    } elsif (defined $n1{$pat}) {
	return $n1{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return 0;
    }
}

sub n2 {
    my $pat = shift;
    if ($patterns) {
	print "$pat\n" if 0 == $myngram{$pat}++;
	return 1;
    } elsif (defined $n2{$pat}) {
	return $n2{$pat};
    } else {
	warn "[$pat] some patterns not found, using zero\n" 
	    unless $ZERO_WARNING++;
	return 0;
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

