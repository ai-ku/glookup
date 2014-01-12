warn q{$Id: submatch.pl,v 1.9 2008/06/27 17:45:23 dyuret Exp $ }."\n";

# submatch($head, $word, [$pos]): finds the versions of the word that
# matches the head in terms of capitalization and morphology.
# Optional pos argument restricts which parts of speech to consider,
# should be one of 'n', 'v', 'a', 'r'.  If not specified, all parts of
# speech will be considered.  Head and word may be multiword phrases,
# in which case the words should be separated by spaces.

use strict;
require 'celex.pl';

my $debug = 0;
my %celexpos;			# celex pos numbers
my %wfids;			# wordform ids for a string
my @wfpos;			# wordform part of speech
my @wfmor;			# wordform inflection type
my @wflm;			# wordform lemma id
my @lmwfs;			# wordform ids for lemma
my @wfstrs;			# wordform strings

# print join(' ', submatch(@ARGV)) . "\n";

sub submatch {
    my ($h0, $w0, $pos) = @_;
#    warn join(',', 'submatch', @_)."\n";
    my $head = headword($h0, $pos);
    my $word = headword($w0, $pos);
    my @wlst1 = mormatch($head, $word, $pos);
#    warn join('/', 'mormatch', @_, '=>', @wlst1) . "\n";
    
    for (@wlst1) {
	my $x = $w0;
	$x =~ s/\Q$word\E/$_/;
	$_ = $x;
    }

    my %wlst2 = ();
    for my $w (@wlst1) {
	my $wcap = capmatch($h0, $w);
	$wlst2{$wcap} = 1;
    }
#    warn join('/', 'submatch', @_, '=>', keys %wlst2) . "\n";

    return keys %wlst2;
}

sub headword {
    my ($phrase, $pos) = @_;
    my @p = split(/\s+/, $phrase);
    if (@p == 1) {
	return $phrase;
    } elsif (defined $pos and $pos eq 'n' and $p[1] ne 'of') {
	return $p[$#p];
    } else {
	return $p[0];
    }
}

sub mormatch {
    my ($head, $word, $pos) = @_;

    # Find celex wordforms for head and word
    mormatch_init() if not %wfids;
    my $headids = $wfids{lc($head)};
    my $wordids = $wfids{lc($word)};
    if (not defined $headids) {
	warn "Warning: [$head] not found in celex\n" if $debug;
	return ($word);
    }
    if (not defined $wordids) {
	warn "Warning: [$word] not found in celex\n" if $debug;
	return ($word);
    }

    # Filter the wordforms for pos if given
    if (defined $pos) {
	die "Unknown pos character [$pos]" if not defined $celexpos{$pos};
	my $cpos = $celexpos{$pos};
	my @headids = grep { $wfpos[$_] eq $cpos } @$headids;
	if (not @headids) {
	    warn "Warning: [$head] not found with pos [$pos] in celex\n" if $debug;
	    return ($word);
	}
	$headids = \@headids;
	my @wordids = grep { $wfpos[$_] eq $cpos } @$wordids;
	if (not @wordids) {
	    warn "Warning: [$word] not found with pos [$pos] in celex\n" if $debug;
	    return ($word);
	}
	$wordids = \@wordids;
    }

    # Output matching forms for word
    my %answer;
    for my $hid (@$headids) {
	my $hpos = $wfpos[$hid];
	my $hmor = $wfmor[$hid];
	
	for my $wid (@$wordids) {
	    my $wpos = $wfpos[$wid];
	    next if $wpos ne $hpos;
	    my $lemma = $wflm[$wid];
	    
	    for my $xid (@{$lmwfs[$lemma]}) {
		my $xmor = $wfmor[$xid];
		if ($xmor eq $hmor) {
		    for my $str (@{$wfstrs[$xid]}) {
			$answer{$str}++;
		    }
		}
	    }
	}
    }

    if (%answer) {
	return keys %answer;
    } else {
	warn "Warning: no wordforms found for [$word] that match [$head] in celex\n" if $debug;
	return ($word);
    }
}

sub capmatch {
    my ($head, $word) = @_;
    if ($head eq lc($head)) {
	return lc($word);
    } elsif ($head eq uc($head) and $head =~ /[A-Z].*[A-Z]/) {
	return uc($word);
    } elsif ($head =~ /^[A-Z]/) {
	return join(' ', map { ucfirst } split(/\s+/, $word));
    } else {
	return $word;
    }
}

sub mormatch_init {
    warn "\nInitializing submatch.pl\n";
    %celexpos = ('n' => 1, 'a' => 2, 'v' => 4, 'r' => 7);
    my @lmpos; readCelexFiles('esl', sub { $lmpos[$_[0]] = $_[3]; });

    readCelexFiles('emw', sub { 
	$wflm[$_[0]] = $_[3]; 
	$wfmor[$_[0]] = $_[4]; 
	$wfpos[$_[0]] = $lmpos[$_[3]];
	push @{$lmwfs[$_[3]]}, $_[0];
    });
    readCelexFiles('eow', sub {
	my $id = $_[0];
	for (my $i = 8; $i <= $#_; $i += 5) {
	    my $str = lc(StripDia(StripSyl($_[$i])));
	    push @{$wfstrs[$id]}, $str;
	    push @{$wfids{$str}}, $id;
	}
    });
}

1;
