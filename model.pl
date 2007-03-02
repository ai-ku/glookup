#!/usr/bin/perl -w
use strict;
my $log2 = log(2);

use Getopt::Long;
my $ngram = 1;
my @A = (1, 1,  1.2, 2.3, 1, 1);
my @B = (0, 0, 1000, 1, 1, 1);
GetOptions('ngram=i' => \$ngram,
	   'a2=f' => \$A[2],
	   'a3=f' => \$A[3],
	   'a4=f' => \$A[4],
	   'a5=f' => \$A[5],
	   'b2=f' => \$B[2],
	   'b3=f' => \$B[3],
	   'b4=f' => \$B[4],
	   'b5=f' => \$B[5]);

require 'gngram.pl';
my $gtotal = ginit();
my $logtotal = log2($gtotal);

my $nword;
my $nbits;
my $nline;

while(<>) {
    $nline++;
    s/\bn't\b/not/g;
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
	$nbits += bits(\@s, $i, $ngram);
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
	    my $extra = $A[$n] * $missing_count + $B[$n];
	    $pb = ($gb + $px * ($missing_count + $extra)) # our smoothing formula
		/ ($ga + $extra);
	} elsif ($missing_count == 0) {

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.  In fact the only three cases are [!], [?], [;].
	    
	    $pb = $px;
	    warn "Warning: New zero missing_count [$b] ga=$ga gb=$gb gc=$gc\n"
		unless $a =~ /^[;!?]$/;
	} else {

# This is a bug in gngram.  This means there are more instances of
# words following "A" than there are instances of "A" by itself.

	    warn "Warning: New negative missing_count [$b] ga=$ga gb=$gb gc=$gc\n"
		unless $a =~ /^A$/;
	    $pb = ($gb + $px) / ($gc + 1);
	}
	return -log2($pb);
    }
}

