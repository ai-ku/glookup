# gtokenize(): Tries to replicate the tokenization of Google ngram
# data.  's, 'd etc. split, n't not split. intra-word dash
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

1;
