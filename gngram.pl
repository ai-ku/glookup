#!/usr/bin/perl -w
use strict;
use IO::File;
use Search::Binary;
use Data::Dumper;
require 'fileio.pl';

$main::GTotal = 0;	       # total number of words from 1gms/total
my $GDataDir = '/mnt/sdc1/google-ngram'; # google data, subdirs 1gms .. 5gms
my $GCachePath = 'gngram.cache'; # file to use as cache
my $GCacheHandle;		# file handle to append to the cache file
my %GCache;			# ngram => freq values from cache
my @GIndex;			# $GIndex[n][k] = first entry in ngm file k

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
    if (not %GCache) { 	ginit(); }
    my $gcount = $GCache{$query};
    if (defined $gcount) { return $gcount; }

    # Could not find in cache, need to look for it in data file:
    # warn "query=$query\n";
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
	warn "Warning: binary_search returned undef (0, $size, $query, $file)\n";
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
    print $GCacheHandle "$query\t$gcount\n";
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
	seek($handle, $pos, SEEK_SET)
	    or die "Cannot seek to $pos";
	if ($pos > 0) { 
	    # BUG?: must find the beginning of entry, but what if pos
	    # is already at the beginning of one.
	    my $ignore = <$handle>;
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
# opens the cache file for appending.

sub ginit {
    readfile("$GDataDir/1gms/total", sub {
	$main::GTotal = 0 + $_;
    });
    warn "ginit: gtotal = $main::GTotal\n";
    for my $n (2, 3, 4, 5) {
	readfile("$GDataDir/${n}gms/${n}gm.idx", sub {
	    push @{$GIndex[$n]}, $_[1];
	}, "\t");
    }
    $GCacheHandle = new IO::File ">>$GCachePath";
    readfile($GCachePath, sub {
	$GCache{$_[0]} = 0 + $_[1];
    }, "\t");
}

# gtokenize(): 's, 'd etc. split, n't not split. intra-word dash
# split. numbers with hyphens slashes etc. not split.

sub gtokenize {
    my $str = shift;
    $str =~ s/(\w)([^\w\s]+)\s/$1 $2 /g;
    $str =~ s/(\w)(\'(s|m|re|ve|ll|d))\b/$1 $2/g;
    $str =~ s/([a-z])-([a-z])/$1 - $2/gi;
    $str =~ s/(\w) n\'t\b/$1n\'t/gi;
    return split(' ', $str);
}

# gfile(): Finds the file in which a given query may be found.

sub gfile {
    my ($query) = @_;
    my @query = split(' ', $query);
    my ($nword, $pair);
    if (scalar(@query) == 3 and $query[1] =~ /^:(\d):$/) {
	$nword = $1;
	die "nword=$nword out of range" if ($nword < 3) or ($nword > 5);
	$pair = 1;
    } else {
	$nword = scalar(@query);
	$pair = 0;
    }
    if ($nword == 1) {
	return "$GDataDir/1gms/vocab";
    } elsif ($nword > 5) {
	die "Do not have $nword-grams";
    }
    my $file;
    my $idx = $GIndex[$nword];
    my $search = $pair ? $query[0] : $query;
    for (my $i = 0; $i <= $#{$idx}; $i++) {
	my $first = $idx->[$i];
	if (($first cmp $search) > 0) {
	    if ($i == 0) { return; }
	    $i--;
	    $file = sprintf("%s/%dgms/%dgm-%04d", $GDataDir, $nword, $nword, $i);
	    last;
	}
    }
    $file = sprintf("%s/%dgms/%dgm-%04d", $GDataDir, $nword, $nword, $#{$idx})
	if not defined $file;
    $file .= '.pair' if $pair;
    return $file;
}

1;
