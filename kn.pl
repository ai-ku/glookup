# kn(s,i,n): returns the number of bits according to an n-gram model
# for the i'th word in sentence s.  n defaults to 5.  For unknown
# words we assume a count of 100, which is half of the minimum count
# for google tokens.  The first token in the sentence array should be
# <S>.  The model needs to be initialized (kn_init) with a file of
# counts before use.  If this is not performed, then the program
# simply writes out the patterns that would be needed for its
# computation.

# Supporting stuff for kn
my @D40 = ( undef, undef, 28.11096928644257128068, 33.98733562666299928484, 39.16948427744310745367, 39.81172633310000762000 );
my @D1 = ( undef, undef, .85948550385118495566, .96034527604501211555, .99565930386648919162, undef );

my $log2 = log(2);
sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }

my $WILD = '_';
my %n0;
my %n1;
my %n2;
my %printed;

sub kn_init {
    my $path = shift;
    open(FP, $path) or die $!;
    while(<FP>) {
	@_ = split(/\t/);
	my $pat = shift;
	$n0{$pat} = shift;
	$n1{$pat} = shift if @_;
	$n2{$pat} = shift if @_;
    }
    close(FP);
}

sub kn {
    my ($s, $i, $n) = @_;

    # Handle $n = 1 and $n out of range:

    $n = 5 if not defined $n;
    if ($n <= 0) {
	die "kn: Bad n [$n]";
    } elsif ($n == 1) {
	my $ny = n0($s->[$i]);
	if ($ny == 0) {
	    # BUG: Any word below 200 count is excluded from the google
	    # data.  We will give such words a count of 100.  Note that
	    # this may happen because of tokenization errors.  It is
	    # not probabilistically sound.
	    $ny = 100;
	}
	my $n_ = n0($WILD);
	return -log2($ny / $n_);
    } elsif ($n > 5) {
	return kn($s, $i, 5);
    } elsif ($n > $i + 1) {
	return kn($s, $i, $i+1);
    }

    # See if the n-1 word context x is any good

    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $nx = n0($x);
    if ($nx == 0) {
	# BUG: instead of backing off here, we should try using the
	# <UNK> token for unknown words in context.
	return kn($s, $i, $n-1);
    } elsif (" $x " =~ / [;!?.] /) { 
        # bad context: google data does not have any regular tokens
	# following these four characters.
	return kn($s, $i, $n-1);
    }

    # Implement Kneser-Ney smoothing

    my $nxy = n0("$x $s->[$i]");
    my $D = $D40[$n];	    # parameter to subtract from the nxy count
    $nxy -= $D if $nxy > 0;	# we only subtract if it is positive
    my $p = $nxy / $nx;	# first approximation

    # Missing count: if we add up nxy for all y, the total will be
    # less than nx, and we will use the difference for smoothing.  The
    # first source of the missing count is the $D we subtract for each
    # y with nxy > 0.  The following gives us the missing count from
    # the D subtractions:

    my $n1x_ = n1("$x _");	# number of distinct words following x
    my $mc1 = $D * $n1x_;	# missing count from D subtractions

    # Second source of the missing count is that the google data has
    # filtered out all ngrams below the count of 40.  Thus there is a
    # missing count between all ngrams that start with x and the count
    # of the n-1 gram x itself:
    
    my $nx_ = n0("$x _");
    my $mc2 = $nx - $nx_;
    
    # Add the kn smoothing term:

    $p += (($mc1 + $mc2) / $nx) * kn2($s, $i, $n-1);

    return -log2($p);
}

sub kn2 {
    my ($s, $i, $n) = @_;

    # Handle out of range $n and the $n=1 special case:

    if (($n <= 0) or ($n >= 5)) {
	die "kn2: Bad n [$n]";
    } elsif ($n == 1) {
	my $n1_y = n1('_ ' . $s->[$i]);
	if ($n1_y == 0) {
	    # Making this 1 is exactly what would happen if we smoothed with a
	    # 0-order 1/V model.  It is not probabilistically sound.
	    $n1_y = 1;
	}
	my $n1__ = n1('_ _');
	return $n1_y / $n1__;
    }

    my $x = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);
    my $n1_x_ = n1("_ $x _");

    # See if the context has non-zero count:
    if ($n1_x_ == 0) {
	# BUG: instead of backing off here, we should try using the
	# <UNK> token for unknown words in context.
	return kn2($s, $i, $n-1);
    }    

    # Implement Kneser-Ney smoothing

    my $n1_xy = n1("_ $x $s->[$i]");
    my $D = $D1[$n];		# parameter to subtract
    $n1_xy -= $D if $n1_xy > 0;	# only subtract if positive count
    my $p = $n1_xy / $n1_x_;	# first approximation

    # Missing count: if we add up n1_xy for all y, the total will be
    # less than n1_x_, and we will use the difference for smoothing.
    # The source of the missing count is the $D we subtract for each y
    # with n1_xy > 0.  How many such y are there?  We can't use n1x_,
    # because of the filtering of ngrams with counts less than 40, nx_
    # is not equal to n_x_.  We can't use n1_x_ because that counts
    # each axb for distinct a and b.  We need a special counter for
    # the number of y with n1_xy > 0.  That is the n2_x_ count.

    my $n2_x_ = n2("_ $x _");
    my $mc = $D * $n2_x_;

    # Add the kn smoothing term:

    $p += ($mc / $n1_x_) * kn2($s, $i, $n-1);
    return $p;
}

sub n0 { cnt($_[0], \%n0) }
sub n1 { cnt($_[0], \%n1) }
sub n2 { cnt($_[0], \%n2) }
sub cnt {
    my ($pat, $cnt) = @_;
    if (not %n0) {
	# if kn not initialized, just print out the patterns needed
	if (not defined $printed{$pat}) {
	    $printed{$pat} = 1;
	    print "$pat\n";
	}
	return 1;
    } elsif (not defined $cnt->{$pat}) {
	die "Error: pattern not found [$pat]";
    } else {
	return $cnt->{$pat};
    }
}

1;
