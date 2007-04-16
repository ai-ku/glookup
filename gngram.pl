#!/usr/bin/perl -w
warn q{$Id: gngram.pl,v 1.19 2007/04/14 11:10:25 dyuret Exp dyuret $ } . "\n";

use strict;
use IO::File;
use Search::Binary;
use Data::Dumper;
require 'fileio.pl';

my $GTotal;		       # total number of words from 1gms/total
my $GDataDir = '/mnt/sdc1/google-ngram'; # google data, subdirs 1gms .. 5gms
my $GCachePath;			# file to use as cache
my $GCacheHandle;		# file handle to append to the cache file
my @GIndex;			# $GIndex[n][k] = first entry in ngm file k
				
my %GCache =			# ngram => freq values from cache
    (				# initialize with some special values
     ':?:' => 1024908267229,	# number of tokens
     '_' => 13588391,		# number of unigrams
     '_ _' => 314843401,	# number of bigrams
     );
     

# gpair(a,b) gives the number of times a and b appear within a five
# word window.

sub gpair {
    my ($a, $b) = @_;
    my $gcount = gngram("$a $b") + gngram("$b $a");
    for (my $n = 3; $n <= 5; $n++) {
	$gcount += gngram("$a :$n: $b");
	$gcount += gngram("$b :$n: $a");
    }
    return $gcount;
}

# gngram gets the count of an ngram string consisting of n space
# separated words.
# if the query is of the form "$a :$n: $b", then it gets the total
# count of all ngrams that start with a and end with b.

sub gngram {
    my ($query) = @_;

    die 'No query' if not defined $query or $query eq '';

    # Initialize and search the cache:
    if (not $GTotal) { ginit(); }
    my $gcount = $GCache{$query};
    if (defined $gcount) { return $gcount; }

    # Could not find in cache, need to look for it in data file:
    # warn "query=$query\n";
    die "[$query] not found in cache" if not @GIndex;
    my $file = gfile($query);
    die "gfile error" if (not defined $file);
    # warn "file=$file\n";
    my $handle = new IO::File;
    $handle->open("< $file") or die "$file: (query=$query) $!";
    my @stat = stat $handle;
    my $size = $stat[7];
    # warn Dumper($handle);

    # "binary_search" implements a generic binary search algorithm returning
    # the position of the first record whose index value is greater than or
    # equal to $val.
    my $pos = binary_search(0, $size, $query, \&gread, $handle);

    # warn "pos=$pos\n";
    if (not defined $pos) {
	# warn "Warning: binary_search returned undef (0, $size, $query, $file)\n";
	return 0;
    }
    seek($handle, $pos, SEEK_SET)
	or die "Cannot seek to $pos";
    my $record = <$handle>;
    $handle->close;
    my ($ngram, $cnt) = split(/[\t\n]/, $record);
    $gcount = (($ngram eq $query) ? $cnt : 0);

    # Record the answer and return:
    $GCache{$query} = $gcount;
    print $GCacheHandle "$query\t$gcount\n" if $GCacheHandle;
    return $gcount;
}

# gread(): During the search the read function will be called with
# three arguments: the input parameters $handle and $val, and a
# position.  If the position is not "undef", the read function should
# read the first whole record starting at or after the position;
# otherwise, the read function should read the record immediately
# following the last record it read.  The search algorithm will
# guarantee that the first call to the read function will not be with
# a position of "undef".  The read function needs to return a two
# element array consisting of the result of comparing $val with the
# index value of the read record and the position of the read
# record. The comparison value must be positive if $val is strictly
# greater than the index value of the read record, 0 if equal, and
# negative if strictly less. Furthermore, the returned position value
# must be greater than or equal to the position the read function was
# called with.

sub gread {
    my ($handle, $val, $pos) = @_;
    # warn Dumper("gread", \@_);
    if (defined $pos) { 
	if ($pos > 0) {
	    seek($handle, $pos - 1, SEEK_SET)
		or die "Cannot seek to $pos";
	    my $ignore = <$handle>;
	} else {
	    seek($handle, $pos, SEEK_SET)
		or die "Cannot seek to $pos";
	}
    }
    my $readpos = tell($handle);
    die "Cannot tell position" if $readpos < 0;
    return if (not defined ($_ = <$handle>));
    my ($ngram, $cnt) = split(/\t/);
    my $readcmp = ($val cmp $ngram);
    # warn "pos=$readpos cmp=$readcmp $ngram\n";
    return ($readcmp, $readpos);
}

# ginit(): reads the cache file, the index (first ngram) files, and
# opens the cache file for appending.  Returns total count.

sub ginit {
    my ($cachefile) = @_;
    return $GTotal if $GTotal;
    if (-r "$GDataDir/1gms/total") {
	readfile("$GDataDir/1gms/total", sub {
	    $GTotal = 0 + $_;
	});
	for my $n (2, 3, 4, 5) {
	    readfile("$GDataDir/${n}gms/${n}gm.idx", sub {
		push @{$GIndex[$n]}, $_[1];
	    }, "\t");
	}
    } else {
	warn "Warning: Cannot find Google ngrams, relying on the cache\n";
	$GTotal = 1024908267229;
    }
    if (defined $cachefile) {
	$GCachePath = $cachefile;
	my $path = $GCachePath;
	if ($GCachePath =~ /\.gz$/) {
	    warn "Warning: Cache is read-only\n";
	    $path = "zcat $path |";
	} else {
	    $GCacheHandle = new IO::File ">>$GCachePath";
	    $GCacheHandle->autoflush(1);
	}
	readfile($path, sub {
	    unless (/^(.+?)\t(\d+)\n$/) {
		warn "Warning: Incomplete cache line [$_]\n";
		return;
	    }
	    $GCache{$1} = 0 + $2;
	});
    } else {
	warn "Warning: No cache file\n";
    }
    warn "ginit: gtotal = $GTotal\n";
    return $GTotal;
}

# gtokenize(): 's, 'd etc. split, n't not split. intra-word dash
# split. numbers with hyphens slashes etc. not split.

sub gtokenize {
    my $str = shift;
    
    # Here is a list of all ascii punctuation:
    # (33-47) !"#$%&'()*+,-./ (58-64) :;<=>?@ (91-96) [\]^_` (123-126) {|}~

    # The following punct characters are always split:
    # ([?] used in urls which we ignore)
    $str =~ s/(\S)([\!\"\$\%\(\)\<\>\?\@\[\\\]\^\_\{\|\}\~])/$1 $2/g;
    $str =~ s/([\!\"\$\%\(\)\<\>\?\@\[\\\]\^\_\{\|\}\~])(\S)/$1 $2/g;

    # Backtick split if not with another backtick:
    $str =~ s/([^\`\s])(\`)/$1 $2/g;
    $str =~ s/(\`)([^\`\s])/$1 $2/g;

    # [&#;] only allowed by themselves or as part of a &#120; escape:
    # (and url's and emoticons, but we'll ignore those for now)
    $str =~ s/\&\#(\S+);/_AB_$1_C_/g;
    $str =~ s/\&(\S+)\;/_A_$1_C_/g;
    $str =~ s/(\S)([\&\#\;])/$1 $2/g;
    $str =~ s/([\&\#\;])(\S)/$1 $2/g;
    $str =~ s/_A_(\S+)_C_/\&$1\;/g;
    $str =~ s/_AB_(\S+)_C_/\&\#$1\;/g;

    # [*=] only allowed by themselves or with other punctuation (except urls):
    $str =~ s/([\pL\pN])([\*\=])/$1 $2/g;
    $str =~ s/([\*\=])([\pL\pN])/$1 $2/g;

    # [+,:] only allowed by themselves or with numbers and other punctuation:
    # (except urls)
    $str =~ s/(\pL)([\+\,\:])/$1 $2/gi;
    $str =~ s/([\+\,\:])(\pL)/$1 $2/gi;

    # This leaves the following punctuation:
    # '-./

    # They seem to have converted contractions into not's:
    $str =~ s/\bn\'t\b/not/g;

    # Do not leave certain punctuation at the end of a word:
    # $str =~ s/(\w)([^\w\s]+)\s/$1 $2 /g;
    $str =~ s/(\w)([\'\-\/]) /$1 $2 /gi;

    # Split all contractions:
    $str =~ s/(\S)(\'(d|em|ll|m|re|s|ve))\b/$1 $2/gi;

    # Intra word dashes and slashes are split:
    # Reason for while: otherwise vis-a-vis does not get split
    while($str =~ s/(\w)([\-\/])(\w)/$1 $2 $3/gi) {};

    # Split periods at the end of a sentence:
    $str =~ s/(\w)(\.\W*)$/$1 $2/;

    return split(' ', $str);
}

# gfile(): Finds the file in which a given query may be found.

sub gfile {
    my ($query) = @_;
    if (not $GTotal) { ginit(); }

    my @query = split(' ', $query);
    my $nword = scalar(@query);
    my $type = '';

    # Handle special query types:
    if (scalar(@query) == 3 and $query[1] =~ /^:(\d):$/) {
	$type = 'pair';
	$nword = $1;
	die "nword=$nword out of range" unless $nword =~ /^[345]$/;
    } elsif ($query[$#query] eq ':?:') {
	$type = 'left';
    } elsif ($query[0] eq '_' and $query[$#query] eq '_') {
	$type = '_x_';
    } elsif ($query[0] eq '_') {
	$type = '_x';
    } elsif ($query[$#query] eq '_') {
	$type = 'x_';
    } elsif ($nword > 3) {
	$type = 'dvd7';
    }

    if ($nword == 1) {
	die unless $type eq '';
	return "$GDataDir/1gms/vocab";
    } elsif ($nword !~ /^[2345]$/) {
	die "Do not have $nword-grams";
    }
    my $search = $query;
    if ($type eq 'pair') { 
	$search = $query[0]; 
    } else {
	# ignore the wildcards for search comparison
	$search =~ s/\s*:\?:\s*//g;
	$search =~ s/\s*_\s*//g;
    }
    my $idx = $GIndex[$nword];
    my $file;
    # warn "query = [$query], nword = $nword, type = $type, search = [$search]\n";
    for (my $i = 0; $i <= $#{$idx}; $i++) {
	my $first = $idx->[$i];
	if ($search lt $first) {
	    # For file 0000 the search query can be less than the first entry:
	    $i-- if $i > 0;
	    $file = sprintf("%s/%dgms/%dgm-%04d", $GDataDir, $nword, $nword, $i);
	    last;
	}
    }
    # For the last file the search query will not be less than the first entry:
    $file = sprintf("%s/%dgms/%dgm-%04d", $GDataDir, $nword, $nword, $#{$idx})
	if not defined $file;
    $file .= ".$type";
    return $file;
}

# gbits(s,i,n): returns the number of bits according to an n-gram model
# for the i'th word in sentence s.  n defaults to 1.  The assumption
# is that the first order model is complete, i.e. there are no unknown
# words.  The first token in the sentence array should be <S>.

# Supporting stuff for gbits
my @A = (undef, undef, 6.7131, 5.9414, 6.5528, 5.7061);

my $log2 = log(2);
sub log2 { log($_[0])/$log2; }
sub exp2 { exp($_[0]*$log2); }
my $log10 = log(10);
sub log10 { log($_[0])/$log10; }
sub exp10 { exp($_[0]*$log10); }

sub gbits {
    my ($s, $i, $n) = @_;
    $n = 1 if not defined $n;
    if ($n == 1) {
	my $g = gngram($s->[$i]);

	# BUG: Any word below 200 count is excluded from the google
	# data.  We will give such words a count of 100.  Note that
	# this may happen because of tokenization errors, and it is
	# not probabilistically sound.

	# die "Unknown word [$s->[$i]]" if $g == 0;
	$g = 100 if $g == 0;

	return log2($GTotal / $g);
    } elsif ($n > $i + 1) {
	return gbits($s, $i, $i+1);
    } else {
	my $a = join(' ', @{$s}[($i-$n+1) .. ($i-1)]);	# a = n-1 word prefix
	my $ga = gngram($a);				# ga = count of a
	my $x = gbits($s, $i, $n-1);			# x = lower order model bits
	if ($ga == 0) { 
#	    warn "Warning: Zero a-count[$a]\n" 
#		if $n <= 3 and $a =~ /[^\w ]/;		# check what is going on with punctuation
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

# If there is no missing count, surprizingly it is better to ignore
# the context and just use the backed-off model.  This happens for
# punctuation marks which probably have buggy counts.  For example,
# semi-colon seems to end sentences, so there are no regular words
# following it.  In fact the only three cases are [!], [?], [;].  I am
# going to add [.] to the list for the same reason even though there
# seem to be bigrams starting with [.].
	    
	my $bad_context = 0;
	for my $tok (@{$s}[($i-$n+1) .. ($i-1)]) {
	    if ($tok =~ /^[;!?.]$/) {
		$bad_context = 1; last;
	    }
	}

	if ($bad_context) {
	    $pb = $px;

	} elsif ($missing_count > 0) {			# apply our smoothing formula
	    my $extra = $A[$n] * $missing_count;
	    $pb = ($gb + $px * ($missing_count + $extra)) 
		/ ($ga + $extra);

	} elsif ($missing_count == 0) {
#	    warn "Warning: Zero missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	    $pb = ($gb + $px) / ($ga + 1);

	} elsif ($missing_count < 0) {

# This is a bug in gngram.  This means there are more instances of
# words following "A" than there are instances of "A" by itself.

	    warn "Warning: Negative missing_count [$b] ga=$ga gb=$gb gc=$gc\n";
	    $pb = ($gb + $px) / ($gc + 1);
	}
	return -log2($pb);
    }
}

1;
