warn '$Id: submatch.pl,v 1.2 2007/03/27 10:25:26 dyuret Exp dyuret $' . "\n";

# submatch($head, $word, [$pos]): finds the versions of the word that
# matches the head in terms of capitalization and morphology.
# Optional pos argument restricts which parts of speech to consider,
# should be one of 'n', 'v', or 'a'.

my %emw;
my %var = ('n' => 2, 'v' => 5, 'a' => 3);

# print join(' ', submatch(@ARGV)) . "\n";

sub submatch {
    my ($head, $word, $pos) = @_;
    my @wlst1 = mormatch($head, $word, $pos);
#    warn join('/', 'mormatch', @_, '=>', @wlst1) . "\n";
    my %wlst2 = ();
    for my $w (@wlst1) {
	my $wcap = capmatch($head, $w);
	$wlst2{$wcap} = 1;
    }
#    warn join('/', 'submatch', @_, '=>', keys %wlst2) . "\n";
    return keys %wlst2;
}

sub mormatch {
    my ($head, $word, $pos) = @_;
    return ($word) unless (defined $pos and defined $var{$pos});
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
    if (%answer) {
	return keys %answer;
    } else {
	return ($word);
    }
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
