#!/usr/bin/perl -w
use strict;
use Getopt::Std;
require 'gngram.pl';
my %opts = ('n' => 1, 'a' => 1, 'b' => 0);
getopt('nab', \%opts);

ginit();
my $nword;
my $nbits;
my $nline;
my $log2 = log(2);
my $logtotal = log2($::GTotal);

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
	$nbits += bits(\@s, $i, $opts{n});
    }
}

printf "%d %f\n", $nword, $nbits/$nword;

sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }

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
	    warn "Warning: Zero a-count[$a]\n";
	    return $x;				      # is this right?
	}
	my $px = exp2(-$x);				# px = lower order model probability
	my $b = $a . " $s->[$i]";			# b = all n words
	my $gb = gngram($b);				# gb = count of b
	my $c = $a . " :?:";				# c = all ngrams that start with a
	my $gc = gngram($c);				# gc = count of c
	my $missing_count = $ga - $gc;			# missing_count = occurances of (a) 
							#   that are missing from ngram data
	my $pb;
	if ($missing_count > 0) {
	    my $extra = $opts{a} * $missing_count + $opts{b};
	    $pb = ($gb + $px * ($missing_count + $extra)) # our smoothing formula
		/ ($ga + $extra);
	} else {

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.
	    
	    # $pb = ($gb + $ga * $opts{c} * $px) / ($ga * (1 + $opts{c})); # if no missing count
	    $pb = $px;
	    warn "Warning: Zero b-count[$b] and zero missing_count\n" if $gb == 0;
	}
	return -log2($pb);
    }
}

sub smooth {
    my ($a) = @_;
    my $ga = gngram($a);
    die if $ga == 0;
    my $gb = gngram($a . " :?:");
    my $h = $gb/$ga;
    die if $h == 1;
    return $h;
}
