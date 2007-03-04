#!/usr/bin/perl -w
use strict;
warn '$Id: model.pl,v 1.10 2007/03/04 10:01:38 dyuret Exp dyuret $' ."\n";

my $log2 = log(2);
sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }
my $log10 = log(10);
sub log10 { log($_[0])/$log10; }
sub exp10 { exp($_[0]*$log10); }

use Getopt::Long;
my $verbose = 0;
my $cachefile;
my $ngram = 1;
my @A = (undef, undef, 4.5725, 5.2775, 4.9850, 3.9550);
my @B = (undef, undef, 2.5375, 2.1850, 2.1450, 2.2475);
GetOptions('cache=s' => \$cachefile,
	   'verbose' => \$verbose,
	   'ngram=i' => \$ngram,
	   'a2=f' => \$A[2],
	   'a3=f' => \$A[3],
	   'a4=f' => \$A[4],
	   'a5=f' => \$A[5],
	   'b2=f' => \$B[2],
	   'b3=f' => \$B[3],
	   'b4=f' => \$B[4],
	   'b5=f' => \$B[5]);

require 'gngram.pl';
my $gtotal = ginit($cachefile);
my $logtotal = log2($gtotal);

my $nword;
my $nbits;
my $nline;

while(<>) {
    $nline++;
    s/\bn\'t\b/not/g;
    s/(\w)([-\/])(\w)/$1 $2 $3/gi;
    s/(\w)([-;,\/\'\%\)]) /$1 $2 /gi;
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
	my $b = bits(\@s, $i, $ngram);
	$nbits += $b;
	if ($verbose) {
	    print $s[$i];
	    for (my $n = 1; $n < $ngram; $n++) {
		printf "\t%.4f", bits(\@s, $i, $n);
	    }
	    printf "\t%.4f\n", $b;
	}
    }
    print "\n" if $verbose;
}

print $nword . ' ' . ($nbits/$nword) . "\n";

# bits(s,i,n): returns the number of bits according to an n-gram model
# for the i'th word in sentence s.  n defaults to 1.  The assumption
# is that the first order model is complete, i.e. there are no unknown
# words.  The first token in the sentence array should be <S>.

sub bits {
    my ($s, $i, $n) = @_;
    $n = 1 if not defined $n;
    if ($n == 1) {
	my $g = gngram($s->[$i]);
	die "Unknown word [$s->[$i]]" if $g == 0;
	return $logtotal - log2($g);
    } elsif ($n > $i + 1) {
	return bits($s, $i, $i+1);
    } else {
	my $a = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);	# a = n-1 word prefix
	my $ga = gngram($a);				# ga = count of a
	my $x = bits($s, $i, $n-1);			# x = lower order model bits
	if ($ga == 0) { 
	    warn "Warning: Zero a-count[$a]\n" 
		if $n <= 3 and $a =~ /[^\w ]/;		# check what is going on with punctuation
	    return $x;					# return lower order result
	}
	my $px = exp2(-$x);				# px = lower order model probability
	my $b = $a . " $s->[$i]";			# b = all n words
	my $gb = gngram($b);				# gb = count of b
	my $c = $a . " :?:";				# c = all ngrams that start with a
	my $gc = gngram($c);				# gc = count of c
	my $missing_count = $ga - $gc;			# missing_count = occurances of (a) 
							#   that are missing from ngram data
	my $pb;
	if ($missing_count > 0) {			# apply our smoothing formula
	    my $extra = $A[$n] * $missing_count + exp10($B[$n]);
	    $pb = ($gb + $px * ($missing_count + $extra)) 
		/ ($ga + $extra);
	} elsif ($missing_count == 0) {

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.  In fact the only three cases are [!], [?], [;].
	    
	    if ($a =~ /[;!?]/) {
		$pb = $px;
	    } else {
		# In other cases we have a legitimate 100% ngram:
		warn "Warning: Zero missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
		$pb = ($gb + $px) / ($ga + 1);
	    }
	} elsif ($missing_count < 0) {

# This is a bug in gngram.  This means there are more instances of
# words following "A" than there are instances of "A" by itself.

	    warn "Warning: Negative missing_count [$b] ga=$ga gb=$gb gc=$gc\n"
		unless $a =~ /^A$/;
	    $pb = ($gb + $px) / ($gc + 1);
	}
	return -log2($pb);
    }
}

