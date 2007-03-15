warn '$Id$' . "\n";

# submatch($head, $word, [$pos]): finds the versions of the word that
# matches the head in terms of capitalization and morphology.
# Optional pos argument restricts which parts of speech to consider,
# should be one of 'n', 'v', or 'a'.

my %emw;
my %var = ('n' => 2, 'a' => 3, 'v' => 5);

# print join(' ', submatch(@ARGV)) . "\n";

sub submatch {
    my ($head, $word, $pos) = @_;
    my @wlst1 = mormatch($head, $word, $pos);
    my %wlst2 = ();
    for my $w (@wlst1) {
	my $wcap = capmatch($head, $w);
	$wlst2{$wcap} = 1;
    }
    # warn join('/', @_, '=>', keys %wlst2) . "\n";
    return keys %wlst2;
}

sub mormatch {
    my ($head, $word, $pos) = @_;
    if (not %emw) {
	open(FP, 'emw.dat') or return ($word);
	while(<FP>) {
	    my @a = split /[\t\n]/;
	    for my $w (@a) {
		push @{$emw{$w}}, \@a;
	    }
	}
	close(FP);
    }
    my $lchead = lc($head);
    my $lcword = lc($word);
    return ($word) unless (defined $emw{$lchead} and defined $emw{$lcword});
    my %answer;
    for my $a (@{$emw{$lchead}}) {
	my $n = scalar(@$a);
	next if defined $pos and $n != $var{$pos};
	for (my $i = 0; $i <= $#{$a}; $i++) {
	    next unless $a->[$i] eq $lchead;
	    for my $b (@{$emw{$lcword}}) {
		next unless $n == scalar(@$b); # wrong pos
		$answer{$b->[$i]} = 1;
	    }
	}
    }
    return keys %answer;
}

sub capmatch {
    my ($head, $word) = @_;
    if ($head =~ /^[a-z\W]+$/) {
	return lc($word);
    } elsif ($head =~ /^[A-Z][a-z\W]*$/) {
	return ucfirst($word);
    } elsif ($head =~ /^[A-Z\W]+$/) {
	return uc($word);
    } else {
	return $word;
    }
}

1;
