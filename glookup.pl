#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Search::Binary;
use IO::File;

my $wfirst = "!";	# "!" is the first word in vocab
my $wlast = "\xc9\x81";	# "\xc981" is the last word in vocab
my $wild = "_";		# "_" does not exist in vocab, so we use it for wildcard

my @n0wild = (1024908267229, 910884463583, 739006848674, 528435661704, 368402367169);
my @n1wild = (     13588391,    314843401,    977069902,   1333820466,   1214460675);
my @n2wild = (     13588391,     12880878,      9329099,      7934066,      7015505);

my (@gindex1, @gindex2);

my $path = '/work/gngram';
my $ngram;
my $allout;
my $debug;
GetOptions('path=s' => \$path, 'ngram=s' => \$ngram, 'all' => \$allout, 'debug' => \$debug);

while(<DATA>) {
    my ($sortkey, $order, $fnum, $line) = split /\t/;
    if ($sortkey == 1) {
	$gindex1[$order][$fnum] = $line;
    } else {
	$gindex2[$order][$fnum] = $line;
    }
}

if (defined $ngram) {
    glookup($ngram);
} else {
    glookup($_) while(<>);
}

sub glookup {
    my ($ngram) = @_;
    my @ngram = split(' ', $ngram);
    $ngram = join(' ', @ngram);
    my $order = scalar(@ngram);
    if ($ngram =~ /^[$wild ]+$/) {
	glookup0($ngram, $order);
    } elsif ($ngram =~ /^[^$wild]/) {
	glookup1($ngram, $order);
    } else {
	glookup2($ngram, $order);
    }
}

sub glookup0 {
    my ($ngram, $n) = @_;
    printf("%s\t%d\t%d\t%d\n", $ngram,
	   $n0wild[$n-1], $n1wild[$n-1], $n2wild[$n-1]);
}

sub gfile1 {
    my ($query, $n) = @_;
    if ($n == 1) {
	return "$path/1gms/vocab";
    }
    my $pat = $query;
    $pat =~ s/$wild/$wfirst/g;
    my $idx = $gindex1[$n];
    my $fnum;
    for ($fnum = 0; $fnum <= $#{$idx}; $fnum++) {
	if ($pat lt $idx->[$fnum]) {
	    $fnum--; last;
	}
    }
    if ($fnum < 0) { $fnum = 0; }
    elsif ($fnum > $#{$idx}) { $fnum = $#{$idx}; }
    return sprintf("$path/${n}gms/${n}gm-%04d", $fnum);
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

sub gread1 {
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

sub glookup1 {
    my ($query, $n) = @_;
    warn "glookup1: query=[$query] n=$n\n" if $debug;
    my $pat = $query;
    $pat =~ s/$wild/$wfirst/g;
    my $file = gfile1($pat, $n);
    warn "file=$file\n" if $debug;
    my $handle = new IO::File;
    $handle->open("< $file") or die "$file: (query=$query) $!";
    my @stat = stat $handle;
    my $size = $stat[7];
    # "binary_search" implements a generic binary search algorithm returning
    # the position of the first record whose index value is greater than or
    # equal to $val.
    my $pos = binary_search(0, $size, $pat, \&gread1, $handle);
    if (not defined $pos) {
	warn "Warning: binary_search returned undef (0, $size, $query, $file)\n";
	return 0;
    }
    $handle->seek($pos, SEEK_SET)
	or die "Cannot seek to $pos";
    
    my $n0 = 0;
    my $n1 = 0;
    my %n2;
    while(<$handle>) {
	chop;
	my ($ngram, $cnt) = split(/\t/);
	my $cmp = ngram_cmp($ngram, $query);
	warn "ngram_cmp[$ngram][$query]=$cmp\n" if $debug;
	if ($cmp == 0) {
	    $n0 += $cnt;
	    $n1++;
	    my ($last) = ($ngram =~ /\S+$/);
	    $n2{$last}++;
	    print "$ngram\t$cnt\n" if $allout;
	} elsif ($cmp == 1) {
	    last;
	}
    }
    if ($handle->eof) {
	warn "Warning: EOF: don't know how to cross files yet.\n";
    }
    $handle->close;
    if ($query !~ /$wild/) {
	print "$query\t$n0\n";
    } else {
	my $nwild;
	$nwild++ while $query =~ /$wild/g;
	if ($nwild > 1 and $query =~ /$wild$/) {
	    my $n2 = scalar(keys(%n2));
	    print "$query\t$n0\t$n1\t$n2\n";
	} else {
	    print "$query\t$n0\t$n1\n";
	}
    }
}

# Binary "cmp" returns -1, 0, or 1 depending on whether the left argu-
# ment is stringwise less than, equal to, or greater than the right
# argument.  Should stop searching if cmp=1.

sub ngram_cmp {
    my ($ngram, $pat) = @_;
    my @ngram = split(' ', $ngram);
    my @pat = split(' ', $pat);
    return 1 unless @ngram == @pat;
    my $mpat = $pat;
    $mpat =~ s/$wild/$wlast/g;
    return 1 if $ngram gt $mpat;
    for (my $i = 0; $i <= $#ngram; $i++) {
	next if $pat[$i] eq $wild;
	return -1 if $pat[$i] ne $ngram[$i];
    }
    return 0;
}

sub glookup2 {
    die "Wildcard at the beginning not implemented yet\n";
}

# Use this to generate the table at the end of the file if data changes
sub init {
    my ($n) = @_;
    if (not defined $n) {
	init($_) for 2..5;
    } elsif ($n < 2 or $n > 5) {
	die "init: Bad ngram order [$n]\n";
    }
    for my $f (glob("$path/${n}gms/${n}gm-????")) {
	my ($fnum) = ($f =~ /gm-(\d+)/);
	open(FP, $f); $gindex1[$n][$fnum] = <FP>; close(FP);
	warn "1\t$n\t$fnum\t$gindex1[$n][$fnum]";
    }
    for my $f (glob("$path/${n}gms/${n}gm-????.sort2")) {
	my ($fnum) = ($f =~ /gm-(\d+).sort2/);
	open(FP, $f); $gindex2[$n][$fnum] = <FP>; close(FP);
	warn "2\t$n\t$fnum\t$gindex2[$n][$fnum]";
    }
}

__END__
1	2	0000	! "	27516372
1	2	0001	, djy	75
1	2	0002	04530 i	75
1	2	0003	17 GIBBS	56
1	2	0004	29425 Tel	124
1	2	0005	6th Path	54
1	2	0006	<S> Challamel	89
1	2	0007	ABCbodybuilding Co	88
1	2	0008	BRIDGEPORT /	736
1	2	0009	Captain Argos	67
1	2	0010	David Sisam	42
1	2	0011	FRAME DOLL	40
1	2	0012	HUSH IN	52
1	2	0013	Johann Sr.	196
1	2	0014	MISIÓN </S>	111
1	2	0015	Network Chevey	92
1	2	0016	Pixel Pro	664
1	2	0017	Russian HIll	66
1	2	0018	Star Sue	435
1	2	0019	U +5231	44
1	2	0020	Zebra Fish	2348
1	2	0021	at toh	41
1	2	0022	civIV forums	79
1	2	0023	doubt linger	188
1	2	0024	frist website	48
1	2	0025	influx takes	106
1	2	0026	measuring Time	85
1	2	0027	ordinate potranno	144
1	2	0028	readiness conference	504
1	2	0029	somehow illuminating	47
1	2	0030	to Hayashi	852
1	2	0031	woman crochet	66
2	2	0000	" !	6592917
2	2	0001	84.872 ,	91
2	2	0002	repines .	83
2	2	0003	125,000 11	79
2	2	0004	Westfield 21	107
2	2	0005	55260 55261	52
2	2	0006	3.0081 </S>	176
2	2	0007	v8.05.02.02 >	204
2	2	0008	TEAL BURGUNDY	109
2	2	0009	/ Cdkn1b	110
2	2	0010	Nathan Di	83
2	2	0011	flowers Florist	509
2	2	0012	Guarantee Hosted	87
2	2	0013	NEW LENSATIC	63
2	2	0014	With Mikhail	390
2	2	0015	Also PLOWS	44
2	2	0016	Ventilation Rats	55
2	2	0017	' ShowSource.pm	98
2	2	0018	<UNK> ToB	144
2	2	0019	value X2	101
2	2	0020	Paninis are	195
2	2	0021	<S> carnivals	3997
2	2	0022	shipment declines	190
2	2	0023	C3 extends	102
2	2	0024	Linear have	81
2	2	0025	awaiting keyboard	46
2	2	0026	by moviemogul27	280
2	2	0027	lash perming	272
2	2	0028	qualified repairs	516
2	2	0029	winter spearing	79
2	2	0030	web toshiba	393
2	2	0031	nicest wooden	44
1	3	0000	! " ''	508
1	3	0001	& UnifyPow "	122
1	3	0002	( Men who	1090
1	3	0003	+ * config.table	57
1	3	0004	, Mary Alvina	87
1	3	0005	, nude hentai	8009
1	3	0006	- Farr (	135
1	3	0007	- relational environment	84
1	3	0008	/ songs *	150
1	3	0009	10 January .	14826
1	3	0010	180 4128 |	91
1	3	0011	2006 Chocolate /	50
1	3	0012	397KB ] </S>	73
1	3	0013	7 û ;	145
1	3	0014	: DarkCounter </S>	43
1	3	0015	<S> 1006 KB	244
1	3	0016	<S> Bell often	41
1	3	0017	<S> Gloss medium	42
1	3	0018	<S> OUTLET ADS	69
1	3	0019	<S> Thinker >	40
1	3	0020	<S> nent model	346
1	3	0021	<UNK> Coastpath Nearby	114
1	3	0022	A MEMBER box	48
1	3	0023	Aitchison ( Sat	48
1	3	0024	Auckland Festival </S>	423
1	3	0025	Bingo - Here	49
1	3	0026	CEOS FORM COALITION	47
1	3	0027	Chicago Membership Posted	52
1	3	0028	Copán Ruinas (	124
1	3	0029	Deep & Skin	1090
1	3	0030	ESSENTIALS IN ONE	69
1	3	0031	FOAM - <UNK>	200
1	3	0032	From southeast star	80
1	3	0033	Gualeguaychú ( ar	122
1	3	0034	Hotel In Kusadasi	157
1	3	0035	Installing xforms .	74
1	3	0036	KIDS DESIGNER CLOTHES	83
1	3	0037	Lexapro 10mg (	132
1	3	0038	Making Noise (	464
1	3	0039	Moment With God	509
1	3	0040	Nice ( adamr.	40
1	3	0041	Our Honey &	41
1	3	0042	Phaser 350 Tektronix	84
1	3	0043	Pullets ; </S>	46
1	3	0044	Residency program at	672
1	3	0045	Sabatini , DA	449
1	3	0046	Simple Guestbook </S>	644
1	3	0047	Students realized that	200
1	3	0048	Tests module for	68
1	3	0049	Tool Dallas Garden	53
1	3	0050	VPD , and	447
1	3	0051	White Leather Gift	58
1	3	0052	] Logging remote	73
1	3	0053	academic and testing	232
1	3	0054	also need comprehensive	145
1	3	0055	and democratizing .	47
1	3	0056	approached Lenin after	78
1	3	0057	auctions are intended	71
1	3	0058	between , pure	81
1	3	0059	by : patriot	444
1	3	0060	ceramics chlorine free	207
1	3	0061	companion part ,	49
1	3	0062	cracking programs )	52
1	3	0063	deserved praise as	108
1	3	0064	driveway is plowed	49
1	3	0065	error executing wwmaster.initiate	70
1	3	0066	feet secretary stockings	48
1	3	0067	for nude lesbians	1955
1	3	0068	gas card 0	106
1	3	0069	hardcore XXX video	2185
1	3	0070	hotels , ideas	407
1	3	0071	in states do	60
1	3	0072	iorr.org , yesterday	93
1	3	0073	la presión a	86
1	3	0074	local buildings such	77
1	3	0075	means officers can	149
1	3	0076	moxie " and	76
1	3	0077	not have NC	115
1	3	0078	of emotionally loaded	185
1	3	0079	one cause at	254
1	3	0080	output is additionally	49
1	3	0081	photos acteress indian	63
1	3	0082	priced proposal and	70
1	3	0083	raised during 2001	44
1	3	0084	research materials *	102
1	3	0085	school right at	234
1	3	0086	shows old school	45
1	3	0087	spermine , and	737
1	3	0088	supplement and after	46
1	3	0089	that includes account	190
1	3	0090	the documents transferred	198
1	3	0091	then will refer	124
1	3	0092	to access Penn	149
1	3	0093	trials will make	224
1	3	0094	valid debts of	149
1	3	0095	well as brutal	720
1	3	0096	with two mechanized	77
1	3	0097	| Finite IT	794
2	3	0000	" ! "	344801
2	3	0001	9 % Tell	127
2	3	0002	daughter 's ultimate	61
2	3	0003	3 ( Trag.	62
2	3	0004	sun4u ) S	82
2	3	0005	<UNK> , 47.4	87
2	3	0006	Biga , Kay	73
2	3	0007	down , amidst	198
2	3	0008	March , life	55
2	3	0009	intrenchments , when	47
2	3	0010	Humanware - India	52
2	3	0011	www.newandusedgolfclubs.com - sponsored	47
2	3	0012	Plus / Ali	137
2	3	0013	NJ 08742 Category	68
2	3	0014	- 14 SEC	80
2	3	0015	NSSE 2001 </S>	515
2	3	0016	Presario 3000 Battery	138
2	3	0017	With 59 points	40
2	3	0018	hole : 0784A	580
2	3	0019	Hydroid : Unknown	110
2	3	0020	° <UNK> Main	44
2	3	0021	, A. Ganguli	40
2	3	0022	Apollo Alliance plan	51
2	3	0023	MULTIPLA BARCHETTA </S>	143
2	3	0024	Zippered Book /	67
2	3	0025	LETTER CUTTING MACHINES	44
2	3	0026	<S> Cleveland who	53
2	3	0027	Sell Cuckoo Clock	51
2	3	0028	and Dir QuickView	94
2	3	0029	and Emily never	53
2	3	0030	Timothy Files by	65
2	3	0031	| Gatorade Thirst	148
2	3	0032	3 Harp Boogie	158
2	3	0033	NURSERY INSPECTION </S>	97
2	3	0034	By Jan. ,	76
2	3	0035	& LOCATION »	231
2	3	0036	or MA student	123
2	3	0037	<S> Memorize my	1934
2	3	0038	<S> NPHS (	453
2	3	0039	WoodSongs Old -	624
2	3	0040	old Parliament Building	59
2	3	0041	Baker Printing ,	63
2	3	0042	Star Recycling </S>	235
2	3	0043	<S> SEWN -	136
2	3	0044	of Sewing Tips	69
2	3	0045	Jiki Statue PRE	232
2	3	0046	<S> Target A3	51
2	3	0047	's Training Technology	113
2	3	0048	Captures Video and	264
2	3	0049	Furnished Wooden Doll	70
2	3	0050	by a 365	320
2	3	0051	claimed actual damages	56
2	3	0052	<UNK> an indigenous	96
2	3	0053	Iranians and give	63
2	3	0054	basis any fees	71
2	3	0055	joint assembly ,	932
2	3	0056	grizzly bear are	260
2	3	0057	a brief annual	1453
2	3	0058	components can co	95
2	3	0059	patient clinics ,	2283
2	3	0060	The continuous design	141
2	3	0061	o de un	2507
2	3	0062	<UNK> div2011.appendChild (	375
2	3	0063	<S> employing strategies	152
2	3	0064	Sidran faced a	44
2	3	0065	prices for Brisas	115
2	3	0066	good friends stayed	49
2	3	0067	nofx goldfinger .	48
2	3	0068	with heart wrapping	470
2	3	0069	by illegal abortion	153
2	3	0070	flashlights in stock	223
2	3	0071	exercises involving several	43
2	3	0072	exercised jurisdiction in	853
2	3	0073	who likes sweet	138
2	3	0074	Ross martin and	49
2	3	0075	numbers more </S>	294
2	3	0076	support no <UNK>	66
2	3	0077	care of Tillie	99
2	3	0078	Cross okc roommates	123
2	3	0079	harmful or nuisance	139
2	3	0080	> past experience	177
2	3	0081	a poorly patched	85
2	3	0082	ZjStream protocol .	54
2	3	0083	10 regional representatives	54
2	3	0084	phone rumor mill	43
2	3	0085	in setup menu	501
2	3	0086	We sold tons	81
2	3	0087	first strike from	470
2	3	0088	lesbians teenage sister	66
2	3	0089	With the Raft	188
2	3	0090	<UNK> the the2	46
2	3	0091	puffy tit girls	45
2	3	0092	alternative to whining	58
2	3	0093	scope units .	354
2	3	0094	class wagon .	65
2	3	0095	singles who 've	60
2	3	0096	very workmanlike performance	44
2	3	0097	Resources | Switzerland	109
1	4	0000	! " '' </S>	500
1	4	0001	$ 2.19 Save :	1021
1	4	0002	's 213T 21 Inch	270
1	4	0003	( Cannon St )	48
1	4	0004	) 2006 Cost (	66
1	4	0005	+6 <UNK> 1.916667 +6	216
1	4	0006	, Erobin , Finucane	43
1	4	0007	, Thomas ( 1830	84
1	4	0008	, do not reconfigure	65
1	4	0009	, proactive problem solving	119
1	4	0010	- $ 172000 </S>	677
1	4	0011	- Jean Rhys :	91
1	4	0012	- high , 25	53
1	4	0013	... More theme of	48
1	4	0014	0.1250 0.0000 0.6250 2.6250	280
1	4	0015	10:09 pm jah shaka	40
1	4	0016	1997 ) : Effect	278
1	4	0017	23 - to 37	107
1	4	0018	40 - fm (	49
1	4	0019	8 PH08 - 19	475
1	4	0020	: Donetsk ; </S>	68
1	4	0021	<S> " That department	70
1	4	0022	<S> 8 pm Red	64
1	4	0023	<S> Banjo ( 26	57
1	4	0024	<S> Dhiraj Joshi ,	61
1	4	0025	<S> Heraldic queries may	65
1	4	0026	<S> Mansion del Valle	103
1	4	0027	<S> Poster : Charlie	435
1	4	0028	<S> Stuart Stephen </S>	58
1	4	0029	<S> VISHNU PURAN (	43
1	4	0030	<S> environmental procedures </S>	45
1	4	0031	<S> they rebelled against	85
1	4	0032	<UNK> <UNK> KIA4558P <UNK>	55
1	4	0033	= : ownerId "	60
1	4	0034	ATOM 386 CB SER	65
1	4	0035	And time doth weary	57
1	4	0036	BOUNCE IS BACK JH	131
1	4	0037	Bracelets - Welcome to	111
1	4	0038	Camp Birch Trails is	81
1	4	0039	Clocks Daylight Daylight Dont	92
1	4	0040	Crisp with Pine Nut	57
1	4	0041	Dialup Networking Installation Instructions	92
1	4	0042	Electronic devices are not	123
1	4	0043	Federal Government -- Judicial	46
1	4	0044	GIS data for USDA	73
1	4	0045	HIV infections among IDUs	217
1	4	0046	However , the DPA	117
1	4	0047	Inc. <UNK> 4601 East	356
1	4	0048	Jesus is relevant to	127
1	4	0049	Laser ML - 2250	89
1	4	0050	MISSED JUMPER by Theobald	97
1	4	0051	Michigan Indian tribal police	63
1	4	0052	Nathan , Lurie &	45
1	4	0053	Offers Broadband Regulatory Relief	55
1	4	0054	Pan Asian Growth Class	87
1	4	0055	Power | Phone and	48
1	4	0056	Random walks , </S>	52
1	4	0057	Rossmoor , California Education	43
1	4	0058	Security : Security Systems	62
1	4	0059	Soto - VICTORIA </S>	310
1	4	0060	Switzerland ( with German	44
1	4	0061	The OECD Convention on	507
1	4	0062	To An Order Already	65
1	4	0063	Upper / Middle <UNK>	102
1	4	0064	Web Ecommerce San Diego	43
1	4	0065	You seemed to know	517
1	4	0066	a big support group	166
1	4	0067	a spokesperson for San	45
1	4	0068	adults ) This is	110
1	4	0069	amid a deepening economic	50
1	4	0070	and abstract of up	45
1	4	0071	and orphan files .	49
1	4	0072	any four year business	46
1	4	0073	as a beautiful little	109
1	4	0074	available for quick review	65
1	4	0075	been , which have	52
1	4	0076	both google.com and google.co.uk	48
1	4	0077	by ignorant persons but	41
1	4	0078	cells in the diagonal	55
1	4	0079	colour and color )	53
1	4	0080	continues on through the	1305
1	4	0081	dancing lessons , a	1065
1	4	0082	dhcp server ) .	89
1	4	0083	driven Sites E -	180
1	4	0084	enjoy myself now ,	43
1	4	0085	extract text and metadata	52
1	4	0086	fits 20 " x	54
1	4	0087	for synchronous systems with	49
1	4	0088	from variations on the	159
1	4	0089	goofs , awards ,	54
1	4	0090	have been in mid	246
1	4	0091	his recommendation . </S>	6943
1	4	0092	implied that the US	713
1	4	0093	in the Seller </S>	48
1	4	0094	intended to provide care	210
1	4	0095	is the smallest 2	85
1	4	0096	kind of scenes .	140
1	4	0097	light by the action	40
1	4	0098	magnet for families </S>	67
1	4	0099	member account , otherwise	8583
1	4	0100	mostly defenseless country .	45
1	4	0101	new projects , such	1020
1	4	0102	nurse in hospital bisexual	616
1	4	0103	of both being the	194
1	4	0104	of the Queen )	867
1	4	0105	on foot and motorcycle	43
1	4	0106	or hereafter have for	49
1	4	0107	owned by an other	60
1	4	0108	pg - 432 Demanar	47
1	4	0109	powerful points about the	80
1	4	0110	protocols , QoS ,	65
1	4	0111	received by the Brigade	49
1	4	0112	response to SuggestedEnhancements .	48
1	4	0113	sbp2 hang when loaded	44
1	4	0114	sex naked pussy free	68
1	4	0115	smtpscan.nasl , v retrieving	50
1	4	0116	starts with 4 cities	79
1	4	0117	support needs are ,	119
1	4	0118	than in ( i	160
1	4	0119	the Enterprise team .	74
1	4	0120	the contact channels of	42
1	4	0121	the most valued area	99
1	4	0122	the tunnel they keep	53
1	4	0123	this article : comment	47
1	4	0124	to RABQSA website </S>	48
1	4	0125	to preempt any local	128
1	4	0126	training and seminars <UNK>	61
1	4	0127	update to 2.1.18 -	196
1	4	0128	vmware win95 win95 Windows	40
1	4	0129	well as a synthetic	401
1	4	0130	will in general perturb	42
1	4	0131	work have nothing to	138
1	4	0132	your business today ,	1547
1	4	0133	« 702 » by	81
2	4	0000	<S> ! " ''	48
2	4	0001	contractor ) with the	104
2	4	0002	stumps ) . </S>	533
2	4	0003	the African American elite	93
2	4	0004	short [ 1 ]	379
2	4	0005	the bazaar area .	150
2	4	0006	that discipline could be	74
2	4	0007	the handouts ) .	106
2	4	0008	simply linking to a	219
2	4	0009	set on input and	48
2	4	0010	the romance element ,	41
2	4	0011	than the usual day	133
2	4	0012	Seward 's Day </S>	120
2	4	0013	but I already posted	212
2	4	0014	command and run some	46
2	4	0015	capital crime in England	49
2	4	0016	calm headlands were And	108
2	4	0017	control mentioned above ,	66
2	4	0018	be orange , not	55
2	4	0019	care services , Mandatory	160
2	4	0020	are the key work	46
2	4	0021	annual vacation entitlement </S>	51
2	4	0022	photograph , comes in	250
2	4	0023	on Dec 5 to	52
2	4	0024	or advertiser asserts that	43
2	4	0025	proposal can be drafted	49
2	4	0026	registered for employment .	107
2	4	0027	pacha is arriving in	102
2	4	0028	one of openness to	161
2	4	0029	press release and distribute	224
2	4	0030	of the Schedule SE	50
2	4	0031	order to understand interactions	150
2	4	0032	State , ExpandedKey +	86
2	4	0033	TORONTO -- Creative Technology	40
2	4	0034	Run : [ BCMSMMSG	3284
2	4	0035	middot - 1150 join	249
2	4	0036	named Preston Michael Spears	80
2	4	0037	management and troubleshooting tasks	214
2	4	0038	in consultation with Tribes	90
2	4	0039	must give a ballot	97
2	4	0040	if legible , may	50
2	4	0041	most often left out	275
2	4	0042	more sense than these	159
2	4	0043	inside the sound of	150
2	4	0044	int width ) Creates	59
2	4	0045	Part C.1 . </S>	79
2	4	0046	Northrop Frye also bought	76
2	4	0047	Middletown Midland Midway Mifflinburg	147
2	4	0048	PTA SECTION FRIENDS '	138
2	4	0049	Safaris Wildwatch Vision Home	125
2	4	0050	Sun from $ 8,698.00	87
2	4	0051	Pap test results of	103
2	4	0052	2 ) Composition /	109
2	4	0053	2 - PC Riding	85
2	4	0054	1077968400 15 1077973200 15	71
2	4	0055	- 864 transexuales -	3199
2	4	0056	gangbangs , booobs penetration	89
2	4	0057	email Basi Subject :	78
2	4	0058	for an enriched quality	53
2	4	0059	distressed by the process	43
2	4	0060	except for a top	88
2	4	0061	easy installation , to	49
2	4	0062	enemy of Israel 's	78
2	4	0063	do proper research ,	147
2	4	0064	find that your enemies	60
2	4	0065	equal treat - ment	53
2	4	0066	Dept 's ltr dated	70
2	4	0067	Journal , and membership	53
2	4	0068	Euro 225,750 <UNK> £	46
2	4	0069	LIGHT BRIGHT AAA Battery	57
2	4	0070	FFF FFe Fe |	50
2	4	0071	LN London ; </S>	152
2	4	0072	Human Resources Job Hotline	69
2	4	0073	Inn Wentworth Plaza ,	88
2	4	0074	Free for Reservations </S>	43
2	4	0075	Linda real estate market	47
2	4	0076	Us ( * =	189
2	4	0077	to 2.2.3 , apply	53
2	4	0078	to a single provision	108
2	4	0079	would become clear :	54
2	4	0080	via email with or	117
2	4	0081	useful if ` xh	68
2	4	0082	was much rarer .	41
2	4	0083	to produce its accounts	203
2	4	0084	will take the sting	664
2	4	0085	time to distill the	110
2	4	0086	world without love is	55
2	4	0087	$ 1,000 minimum </S>	760
2	4	0088	, <UNK> , theonecalledh	189
2	4	0089	' Eddie is the	53
2	4	0090	' Mechanism , Fixed	253
2	4	0091	, TX ( 2nd	55
2	4	0092	, and the Tweeter	74
2	4	0093	, expertly arranged by	107
2	4	0094	's minimum , and	118
2	4	0095	( sic ) along	58
2	4	0096	, who ranks </S>	47
2	4	0097	[ 21 ] BH	68
2	4	0098	YOUR DETAILS / This	47
2	4	0099	Your Move - SMARTpages.com	85
2	4	0100	an Urn by Jan	60
2	4	0101	aid being offered .	58
2	4	0102	The entrance exam for	135
2	4	0103	and initial consultation </S>	40
2	4	0104	ambiguity on this question	50
2	4	0105	all taxes , tip	1282
2	4	0106	Index » While Salsa	52
2	4	0107	| The Presidency and	609
2	4	0108	ACM ** null **	175
2	4	0109	<UNK> - y #	151
2	4	0110	<UNK> <UNK> <UNK> 752400	59
2	4	0111	Benser Christopher Hall --	253
2	4	0112	Communications Inc. Marc Hyman	71
2	4	0113	Beth R. Hart <UNK>	57
2	4	0114	Boris Yeltsin in March	48
2	4	0115	Complainant in person </S>	73
2	4	0116	Corporate profile : TRICO	46
2	4	0117	<S> weighing , </S>	79
2	4	0118	- COURT ROLL -	41
2	4	0119	- Indigenous people [	99
2	4	0120	23 SCOTT Neil m	49
2	4	0121	- beginner - intermediate	115
2	4	0122	/ may not cause	77
2	4	0123	200,000 towards the project	47
2	4	0124	<S> 4050 . </S>	1227
2	4	0125	<S> At Weill Cornell	50
2	4	0126	<S> DRAPED TABLES -	46
2	4	0127	<S> HANES ® </S>	54
2	4	0128	<S> Listen to Grateful	59
2	4	0129	<S> Over 250 species	516
2	4	0130	<S> Scholarship includes traditional	47
2	4	0131	<S> The minimum IELTS	145
2	4	0132	<S> Yeeeeah . </S>	257
2	4	0133	: replicate - wild	83
1	5	0000	! " ) ) )	213
1	5	0001	$ 300 - $ 1,300	106
1	5	0002	's Tennis Headlines : </S>	104
1	5	0003	( I know it all	60
1	5	0004	) :383 - 9 </S>	42
1	5	0005	, " unto , "	76
1	5	0006	, Gloucestershire , Up to	99
1	5	0007	, Wigton , Cumbria .	147
1	5	0008	, cosmetic dentist , restorative	59
1	5	0009	, one person is "	41
1	5	0010	, truck , SUV bumper	554
1	5	0011	- AutoClass - 0.09.tar.gz 01	74
1	5	0012	- based Mental Health </S>	170
1	5	0013	--- [ 261 <UNK> Adder	109
1	5	0014	0 ENDBLK 5 6D 100	54
1	5	0015	1094676301 0.00 1094679902 0.00 1094683502	70
1	5	0016	2 5809 <UNK> <UNK> <UNK>	46
1	5	0017	29 7 00000 12 <UNK>	599
1	5	0018	5350 W Lovers Ln Ste	55
1	5	0019	: 1984 Sales Volume :	59
1	5	0020	: The power of silence	64
1	5	0021	<S> 11 Q. Tell us	93
1	5	0022	<S> All beads are handcrafted	147
1	5	0023	<S> Commissioned May 1943 .	41
1	5	0024	<S> Garmin GPSMAP 2010 GPS	816
1	5	0025	<S> Jon Alexander Furniture -	385
1	5	0026	<S> Once the fat and	48
1	5	0027	<S> Second time parents or	78
1	5	0028	<S> The practice contains a	105
1	5	0029	<S> With online classes ,	81
1	5	0030	<S> not an option ;	44
1	5	0031	<UNK> - Do not source	103
1	5	0032	<UNK> cos <UNK> 60 20	74
1	5	0033	A tree had fallen on	60
1	5	0034	American Academy of Dermatology 63rd	50
1	5	0035	BELIEFS OUR DONORS OUR FACILITIES	61
1	5	0036	Britons , Gaels , and	78
1	5	0037	Catalog of Research on Secondary	47
1	5	0038	ConcernArticle 4 Destruction of Organisms	637
1	5	0039	Day of the time appointed	56
1	5	0040	EX SHOES IN BOX <UNK>	294
1	5	0041	Fault - Tolerant Broadcast of	261
1	5	0042	GRADED MS 70 1981 -	52
1	5	0043	Has primary responsibility for the	172
1	5	0044	I dont know how fast	108
1	5	0045	Initial Available Appropriations for 2005	50
1	5	0046	Jungle - > ( 20	85
1	5	0047	Linux on Cygwin on Linux	44
1	5	0048	Matt Nielsen has signed a	41
1	5	0049	NC6000 - Centrino 1.4 GHz	126
1	5	0050	Oakdale , Louisiana Financial </S>	90
1	5	0051	Pantalla & Limpieza ( 23	54
1	5	0052	Price Match on Cardinals Tickets	56
1	5	0053	Red Ranger Gun Set </S>	86
1	5	0054	SUBJECT : RE : Oklahoma	83
1	5	0055	Similarly , if t </S>	90
1	5	0056	Subject : Business Company Name	57
1	5	0057	The Foregate ' . </S>	50
1	5	0058	Tires ( 3 ) Trade	50
1	5	0059	Utilities > Communication Tools </S>	68
1	5	0060	What new technology is available	48
1	5	0061	[ MXZ ] 174 C.	41
1	5	0062	a few other examples from	54
1	5	0063	a women 's ministry and	40
1	5	0064	alkaline plutonism along the Pacific	42
1	5	0065	and Counting - out Rhymes	84
1	5	0066	and insurance plans for employees	60
1	5	0067	and will not exceed 3	41
1	5	0068	arpa / CVS / Repository	104
1	5	0069	author a note if you	113
1	5	0070	becoming aware of : </S>	114
1	5	0071	broad - band filter .	44
1	5	0072	camping packages or coach tours	103
1	5	0073	claim and self - </S>	47
1	5	0074	contact the primary contact officer	109
1	5	0075	day to mess with me	51
1	5	0076	displayed on interactive seating maps	46
1	5	0077	elasticity & firms up the	84
1	5	0078	explosion or other release of	47
1	5	0079	follow the argument . </S>	708
1	5	0080	forgive us when we make	49
1	5	0081	general , abstract ideas that	42
1	5	0082	happiness , happy , infant	46
1	5	0083	here for the quest for	68
1	5	0084	if you really cared for	126
1	5	0085	in someone 's head ,	642
1	5	0086	insurance , customs duties and	46
1	5	0087	is revealing about this book	49
1	5	0088	kdm / index.docbook / usr	1298
1	5	0089	like to beg for your	69
1	5	0090	man who has learned not	84
1	5	0091	misleading campaign , the people	78
1	5	0092	neat stuff Like Dolls Toys	189
1	5	0093	note - taking , preparing	43
1	5	0094	of another organism of the	55
1	5	0095	of the book ' La	612
1	5	0096	on current development and humanitarian	83
1	5	0097	or loose abdominal skin that	717
1	5	0098	pardon was granted to the	73
1	5	0099	plug in to your TV	102
1	5	0100	project , is in a	47
1	5	0101	recessive dystrophic epidermolysis bullosa in	49
1	5	0102	rigidity , which translates into	45
1	5	0103	semi - detached house :	73
1	5	0104	sites run efficiently 24 hours	45
1	5	0105	stated as a sum certain	441
1	5	0106	tail end of a </S>	91
1	5	0107	that student should be allowed	61
1	5	0108	the Village of Jefferson to	53
1	5	0109	the framework of the 2001	166
1	5	0110	the project has played an	46
1	5	0111	them how to react when	60
1	5	0112	through the option byte .	48
1	5	0113	to dwell in a land	86
1	5	0114	to the education that I	79
1	5	0115	types of office equipment .	386
1	5	0116	very hard after he entered	210
1	5	0117	we think that the dresses	66
1	5	0118	why we are facing a	78
1	5	0119	women , ranging in ages	79
1	5	0120	you the ability to participate	136
1	5	0121	£ 9.95 £ 30.00 </S>	89
2	5	0000	" ! " ) )	108
2	5	0001	- " Baraka " August	67
2	5	0002	table - chases away chills	326
2	5	0003	to MySQL 4.1 from 4.0	82
2	5	0004	to any messages that have	63
2	5	0005	stud cheats games flush poker	90
2	5	0006	shall exercise the following rights	59
2	5	0007	to i - 1 .	190
2	5	0008	so many people denounced for	50
2	5	0009	the outboard effects ? </S>	40
2	5	0010	the significance of the built	68
2	5	0011	strangers to us , and	95
2	5	0012	] 28th May 2005 </S>	41
2	5	0013	and RTF formats are not	67
2	5	0014	area and temporarily remove two	62
2	5	0015	by country ... . </S>	214
2	5	0016	and garage ( supplement )	65
2	5	0017	and kill Kira , and	43
2	5	0018	arrival of the first members	64
2	5	0019	are reproductive and developmental toxicants	40
2	5	0020	and the abdominal skin )	46
2	5	0021	are ugly people everywhere .	77
2	5	0022	char * string , +	279
2	5	0023	cockpit -- the A2 Aviator	85
2	5	0024	give a heads up about	144
2	5	0025	d avions pas cher paris	41
2	5	0026	comprehensive educational , training and	64
2	5	0027	girls in g strings black	145
2	5	0028	deserves no stars . </S>	109
2	5	0029	carries out a multi -	46
2	5	0030	didn t have much time	79
2	5	0031	consider to be between bright	59
2	5	0032	everything would be hunky dory	194
2	5	0033	lha - 1.14i - 4.i386.rpm	80
2	5	0034	instructions Shipment will be made	507
2	5	0035	history and theology , and	256
2	5	0036	he could avoid most of	81
2	5	0037	man has a piece of	54
2	5	0038	is limited to 6 persons	104
2	5	0039	me on : i love	43
2	5	0040	is set when a signal	52
2	5	0041	lacks the polish and refinement	61
2	5	0042	he wants to go first	76
2	5	0043	scooters , gas skateboard ,	125
2	5	0044	of DDR RAM ( expandable	47
2	5	0045	once again , I just	411
2	5	0046	provided by well trained ,	41
2	5	0047	priorities for the scientific and	50
2	5	0048	principal issue at trial was	47
2	5	0049	number of the VISM -	118
2	5	0050	or revocation of licenses ,	314
2	5	0051	on the census of 24	130
2	5	0052	now to explain a little	45
2	5	0053	Games & Game Pieces Playing	52
2	5	0054	Input , Mono Output -	84
2	5	0055	L 01 - Aug -	216
2	5	0056	Lyrics > The Beatles >	82
2	5	0057	Luca Deidda DP26 Do We	43
2	5	0058	IN KEY STAGES 1 AND	3595
2	5	0059	Free Press , 2006 .	40
2	5	0060	FUNCTIONALLY UNCLASSIFIED <UNK> UNCLASSIFIED LOCALIZATION	40
2	5	0061	January deadline to phase out	44
2	5	0062	Large on / off switch	90
2	5	0063	Heritage | Hawaii Marketplace |	58
2	5	0064	Zhou , Chew Lim Tan	69
2	5	0065	Student / Family Reunification Team	42
2	5	0066	Review <UNK> Uptime : 100.0	126
2	5	0067	Toxic Crusaders Transformers Tranzor Z	42
2	5	0068	RU International Emergency Assistance </S>	64
2	5	0069	The Pope is the Catholic	107
2	5	0070	Summers Tits Amateur Showering Photos	7305
2	5	0071	( 1.1 million pound )	166
2	5	0072	( <UNK> ) with views	60
2	5	0073	, Ending Thursday , 16	57
2	5	0074	' Mo , Vita Lyrics	49
2	5	0075	" They gave the impression	129
2	5	0076	, and the newly incorporated	41
2	5	0077	( errors == '' )	92
2	5	0078	) located at 1883 Grant	64
2	5	0079	, sardonic killers who dispatch	42
2	5	0080	, version 1.93A or higher	55
2	5	0081	Posting date , Material ,	55
2	5	0082	School of Psychiatry at the	82
2	5	0083	The women fare wonderfully ,	88
2	5	0084	00 ) | 12/11/2005 (	74
2	5	0085	: - reading file /	53
2	5	0086	/ 2 , Kissimmee </S>	40
2	5	0087	7.16 <UNK> ... 1 /	65
2	5	0088	Crosley + Company | </S>	56
2	5	0089	CD - ROM Garmin MapSource	44
2	5	0090	Cast : Groucho Marx ,	182
2	5	0091	<UNK> BDA Shop Branches and	822
2	5	0092	CV FAQ 's and answers	68
2	5	0093	California Medi - Cal Long	48
2	5	0094	<UNK> Seventy ( 70 )	53
2	5	0095	<S> * All - Conference	205
2	5	0096	: Eros Therapy <UNK> </S>	675
2	5	0097	- NULL <UNK> < -	466
2	5	0098	- Then All Was Silent	1041
2	5	0099	45 days after the end	8034
2	5	0100	- patterning ) . </S>	97
2	5	0101	- z " $ prefix	81
2	5	0102	Check for Missing Five Element	234
2	5	0103	Animal sex with . </S>	301
2	5	0104	you - I used to	138
2	5	0105	| News | Feedback /	53
2	5	0106	would be apt to come	42
2	5	0107	very good and contributes significantly	100
2	5	0108	with non - <UNK> </S>	278
2	5	0109	to scale up production and	78
2	5	0110	to the subtle power of	49
2	5	0111	with worldwide sales . </S>	158
2	5	0112	<S> 4 ( Fall 1997	139
2	5	0113	<S> As with any caching	42
2	5	0114	<S> Day of month [	43
2	5	0115	<S> He will be starting	1078
2	5	0116	<S> Litigation Department of the	122
2	5	0117	<S> Peridot , Iolite &	108
2	5	0118	<S> Sometimes , the state	133
2	5	0119	<S> This forum has 2095	201
2	5	0120	<S> accredited online university degree	96
2	5	0121	<S> should be modified ,	212
