/* -*- mode: C; mode: Outline-minor; outline-regexp: "/[*][*]+"; -*- */
const char *rcsid = "$Id: glookup.c,v 1.5 2007/11/21 16:47:08 dyuret Exp dyuret $";
const char *help = "glookup [-p google-ngram-path] [-a] < patterns > counts";

/** glookup.c: This program will read ngram patterns from stdin,
 * lookup their counts in one pass from google ngram data, and print
 * the counts.  The -p option can be used to specify the location of
 * the google ngram data, otherwise the current directory will be
 * assumed.  The input should have a single pattern on each line
 * consisting of space separated tokens with '_' representing the
 * wildcard token that matches any word.  The output will have up to
 * three counts (tab separated) next to the pattern:
 *
 * n0: computes the total count of the ngrams matching a given pattern
 * n1: computes the number of distinct ngrams matching a given pattern
 * n2: computes the number of distinct words that appear as the last
 * word in a pattern that ends with a wildcard and has more than one
 * wildcard.  This is needed for kn smoothing.
 *
 * I added an -a option to glookup, which not only prints out counts
 * for the patterns, but also prints out the actual ngrams with their
 * counts that matched the patterns.  This is useful if you are trying
 * to figure out the most likely words that complete a given pattern.
 * The printing of actual ngrams is supressed for patterns that
 * consist of wildcards only like '_ _ _'.
 *
 * If needed, here are the counts for the all-wildcard patterns:
 * const guint64 n0wild[] = {0, 1024908267229, 910884463583, 739006848674, 528435661704, 368402367169 };
 * const guint32 n1wild[] = {0, 13588391, 314843401, 977069902, 1333820466, 1214460675};
 * const guint32 n2wild[] = {0, 13588391, 12880878, 9329099, 7934066, 7015505};
 *
 */

/** Kneser-Ney Discounting Primer
 *
 *  KN discounting generally gives superior results and glookup has
 *  been written to extract all the counts necessary.  Here is a brief
 *  overview of interpolated KN discounting as described in Chen and
 *  Goodman 1998.  (The original KN is backed-off not interpolated.
 *  Chen and Goodman suggest that interpolation gives better results.
 *  They also define a modified KN which uses different discount
 *  factors D for each order n.)
 *
 *  For the top level ngram:
 *
 *	      max(0, n0(1..n)-D) + D*n1(1.._)*Q(2..n)
 *  P(1..n) = ---------------------------------------
 *	                    n0(1.._)
 *
 *  (1..n) represents an n word ngram, and (1.._) represents the ngram
 *  with the last word replaced with a wildcard and n0, n1, n2 have
 *  the meanings defined above.  P is guaranteed to sum to 1 over all
 *  the nth words if Q adds up to 1.
 *
 *  For the lower level ngrams, we'd like to make their probability
 *  proportional to the number of unique words that precede them.
 *
 *	      max(0, n1(_2..n)-D) + D*n2(_2.._)*Q(3..n)
 *  Q(2..n) = -----------------------------------------
 *			     n1(_2.._)
 *
 *  And finally when we come to the unigram model:
 *
 *	   max(0, n1(_n)-D) + D*n2(__)*(1/V)
 *  Q(n) = ----------------------------------
 *		        n1(__)
 *
 *  That is if we wanted to interpolate it with the 0'th order model
 *  1/V.  Otherwise we can just take n1(_n)/n1(__) as our Q(n).
 *
 *  Of course if the context 1..(n-1) has zero count, we should start
 *  with P(2..n), not Q(2..n).
 *
 *  There are a couple of complications associated with the Google
 *  ngram counts.  The counts are incomplete, ngrams with counts less
 *  than 40 have been omitted.  Thus for example we cannot use n0(12)
 *  instead of n0(12_).  Secondly, the best D used for P is found to
 *  be close to 40 because the minimum n0 number it will be subtracted
 *  from is 40.  The best D for Q is found to be close to 1.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <glib.h>
#include "foreach.h"
#include "procinfo.h"
#include "ghashx.h"


/** Command line processing */

static char *path = ".";	/* google-ngram path */
static int all_ngrams = 0;	/* whether to print all ngrams matching a wildcard pattern */

static struct option long_options[] = {
  {"path", 1, 0, 'p'},
  {"all", 0, 0, 'a'},
  {0, 0, 0, 0}
};

void get_options(int argc, char *argv[]) {
  while(1) {
    int option_index = 0;
    int c = getopt_long_only(argc, argv, "p:a", long_options, &option_index);
    if (c == -1) break;
    switch(c) {
    case 'p': path = optarg; break;
    case 'a': all_ngrams = 1; break;
    default: g_message(help); exit(0);
    }
  }
}

/** Word tokens are represented as guint32 integers (GQuarks).  The
 *  data types like GQuark, GString etc. are defined in glib2.  The
 *  string "_" and its associated GQuark represent the wildcard
 *  token. */

typedef GQuark Token;
static Token ANY;		/* GQuark corresponding to wildcard */
#define WILDCARD_STR "_"	/* String corresponding to wildcard */
void wildcard_init() { ANY = g_quark_from_string(WILDCARD_STR); }


/** Ngrams and related structures are represented as guint32 arrays
 * where the first element of the array gives the length.  The other
 * elements are positive integers representing tokens, special token
 * ANY represents a wildcard.  The three types (Ngram, Pattern, Mask)
 * differ only by intent of usage: Ngram represents ngrams read from
 * data, Pattern represents patterns with wildcards that the user
 * supplies, and Mask represents different types of pattern determined
 * by the position of wildcards.
 */

#define NGRAM_ORDER 5		/* maximum ngram order */
typedef Token *Ngram, *Pattern, *Mask;
typedef Token EmptyPattern[NGRAM_ORDER + 1];

#define ngram_size(p) (*(p))
#define ngram_copy(p) g_memdup(p, (ngram_size(p)+1)*sizeof(Token))

#define foreach_tok(var, ngm)\
  for (register guint32 _i = 1, _l = ngram_size(ngm), var = 0;\
       (_i <= _l) && ((var = ngm[_i]) || 1); _i++)

void ngram_from_string(Ngram ngm, char *str) {
  int ntok = 0;
  foreach_token(word, str) {
    g_assert(ntok < NGRAM_ORDER);
    ngm[++ntok] = g_quark_from_string(word);
  }
  ngm[0] = ntok;
}

void ngram_to_gstring(GString *g, Ngram pat) {
  int n = 0;
  g_string_truncate(g, 0);
  foreach_tok(pi, pat) {
    if (n++) g_string_append_c(g, ' ');
    g_string_append(g, g_quark_to_string(pi));
  }
}

int ngram_count_wildcards(Ngram pat) {
  int nw = 0;
  foreach_tok(pi, pat) {
    if (pi == ANY) nw++;
  }
  return nw;
}

void ngram_free(gpointer p) { free(p); }


/** Counters: we keep track of up to three types of counts.
 * n0 is the total number of times matching ngrams occur.
 * n1 is the total number of distinct ngrams that match.
 * n2 is the hash of distinct last words that match.
 * (used in kn smoothing) */

typedef struct CounterS {
  guint64 n0;
  guint32 n1;
  Hash n2;
} *Counter;

Counter counter_new(Pattern pat) {
  Counter cnt = g_new0(struct CounterS, 1);
  int nwild = ngram_count_wildcards(pat);
  if ((nwild > 1) &&
      (pat[ngram_size(pat)] == ANY))
    cnt->n2 = g_hash_table_new(g_direct_hash, g_direct_equal);
  return cnt;
}

void counter_free(gpointer val) {
  Counter cnt = val;
  if (cnt->n2 != NULL) free_hash(cnt->n2);
  free(cnt);
}


/** Pattern hash is initialized with the patterns from the user and
    empty counters.  The counters are updated when going over the
    ngram files. */

static guint pattern_hash_rnd[NGRAM_ORDER + 1];

guint pattern_hash(Pattern pat) {
  guint hash = 0;
  for (int i = ngram_size(pat); i > 0; i--) {
    hash += (1+pat[i]) * pattern_hash_rnd[i];
  }
  return hash;
}

gboolean pattern_equal(Pattern pa, Pattern pb) {
  if (ngram_size(pa) != ngram_size(pb))
    return 0;
  for (int i = ngram_size(pa); i > 0; i--) {
    if (pa[i] != pb[i]) return 0;
  }
  return 1;
}

Hash pattern_hash_init() {
  for (int i = NGRAM_ORDER; i > 0; i--) {
    pattern_hash_rnd[i] = rand();
  }
  return g_hash_table_new_full((guint(*)(gconstpointer)) pattern_hash,
			       (gboolean(*)(gconstpointer,gconstpointer)) pattern_equal,
			       ngram_free, counter_free);
}

/* To print the pattern counts after going through google ngrams */

void pattern_hash_print_fn(Pattern pat, Counter cnt, GString *gstr) {
  ngram_to_gstring(gstr, pat);
  fputs(gstr->str, stdout);
  fputc('\t', stdout);
  printf("%" G_GUINT64_FORMAT, cnt->n0);
  int nwild = ngram_count_wildcards(pat);
  if (nwild > 0) printf("\t%u", cnt->n1);
  int n = ngram_size(pat);
  if (pat[n] == ANY) {
    printf("\t%u",
	   (cnt->n2 != NULL) ? hlen(cnt->n2)
	   : nwild == 1 ? cnt->n1
	   : cnt->n0 == 0 ? 0
	   : (g_message("Warning: No n2 for [%s]", gstr->str), 0));
  }
  fputc('\n', stdout);
}

void pattern_hash_print(Hash patterns) {
  GString *gstr = g_string_new(NULL);
  hmap(patterns, pattern_hash_print_fn, gstr);
  g_string_free(gstr, TRUE);
}


/** Masks represent types of patterns based on their wildcard
 * positions.  They are important for efficiency: while going over a
 * 5-gram data file one would need to do 32 pattern hash lookups for
 * each ngram to cover all wildcard combinations.  With masks we only
 * check the actual types that the user asked for. */

/* mask_match(msk, pat) returns true if the wildcard positions in pat
   matches ones in msk */

gboolean mask_match(Mask msk, Pattern pat) {
  if (ngram_size(msk) != ngram_size(pat)) return 0;
  for (int i = ngram_size(msk); i > 0; i--) {
    if ((pat[i] == ANY) != (msk[i] == ANY))
      return 0;
  }
  return 1;
}

/* mask_apply(msk, ngm, pat) copies ngm to pat, replacing certain
   positions indicated by msk with the wildcard */

void mask_apply(Mask msk, Ngram ngm, Pattern pat) {
  g_assert(ngram_size(msk) == ngram_size(ngm));
  ngram_size(pat) = ngram_size(ngm);
  for (int i = ngram_size(ngm); i > 0; i--) {
    pat[i] = ((msk[i] == ANY) ? msk[i] : ngm[i]);
  }
}

#define foreach_mask(m, a)\
  for (register guint32 *m = NULL, _i = 0; (m=a[_i])!=NULL; _i++)

/* mask_table_add(pat, cnt, table) adds a mask for the pattern pat to
   the table if it does not already exist.  The cnt variable is the
   value in the hash table iterated over and is ignored. */

void mask_table_add(Pattern pat, Counter cnt, Mask **table) {
  int n = ngram_size(pat);
  Mask *masks = table[n];
  int mask_count = 0;
  foreach_mask(m, masks) {
    if (mask_match(m, pat)) return;
    mask_count++;
  }
  g_assert(mask_count < (1<<n));
  masks[mask_count] = ngram_copy(pat);
}

/* The mask_table[n] is an array of masks for patterns of length n
   given by the user */

Mask **mask_table_init(Hash patterns) {
  Mask **mask_table = g_malloc0((1 + NGRAM_ORDER) * sizeof(Mask *));
  for (int n = NGRAM_ORDER; n > 0; n--) {
    mask_table[n] = g_malloc0((1+(1<<n)) * sizeof(Mask));
  }
  hmap(patterns, mask_table_add, mask_table);
  return mask_table;
}

void mask_table_free(Mask **mask_table) {
  for (int n = NGRAM_ORDER; n > 0; n--) {
    free(mask_table[n]);
  }
  free(mask_table);
}

void mask_table_dump(Mask **mask_table) {
  GString *g = g_string_new(NULL);
  for (int n = 1; n <= NGRAM_ORDER; n++) {
    Mask *masks = mask_table[n];
    foreach_mask(m, masks) {
      ngram_to_gstring(g, m);
      fputs(g->str, stderr);
      fputc('\n', stderr);
    }
  }
  g_string_free(g, TRUE);
}


/** read_patterns(): Reads ngram patterns from stdin and returns a
 * hash table where keys are the patterns and the values are empty
 * counters. */

Hash read_patterns() {
  EmptyPattern pat;
  wildcard_init();
  Hash patterns = pattern_hash_init();
  foreach_line(str, NULL) {
    ngram_from_string(pat, str);
    if (!hgot(patterns, pat)) {
      Counter cnt = counter_new(pat);
      hset(patterns, ngram_copy(pat), cnt);
    }
  }
  return patterns;
}

/** count_ngram(patterns, ngm, ngram_cnt, lastword, masks) increments
 *  the appropriate counts in the patterns hash looking at ngram ngm
 *  and its count ngram_cnt.  The parameter masks is an array of masks
 *  for the ngram size of ngm.
 */

void count_ngram(Hash patterns, Ngram ngm, guint64 ngram_cnt, Mask *masks) {
  EmptyPattern pat;
  int nowild = 0;	/* whether this ngram matches exactly */
  int nmatch = 0;	/* whether this ngram matches any patterns
			   (except all wildcards) */

  foreach_mask(m, masks) {
    mask_apply(m, ngm, pat);
    Counter cnt = hget(patterns, pat);
    if (cnt != NULL) {
      int nwild = ngram_count_wildcards(pat);
      if (nwild == 0) {
	g_assert(cnt->n0 == 0);
	cnt->n0 = ngram_cnt;
	nowild = 1;
      } else {
	int n = ngram_size(ngm);
	if (nwild < n) nmatch = 1;
	cnt->n0 += ngram_cnt;
	cnt->n1++;
	if (cnt->n2 != NULL) {
	  Token q = ngm[n];
	  hset(cnt->n2, GUINT_TO_POINTER(q), 1);
	}
      }
    }
  }

  /* if we have a match for this ngram (nmatch) and (1) the user wants
     us to count all matching ngrams as well (all_ngrams), (2) the
     ngram itself is not in the pattern list (!nowild), we add the
     ngram with its count to the pattern list.  The only exception is
     when the only pattern the ngram matched is an all wild pattern,
     in which case the user probably did not want us to dump the whole
     database.  The (nwild < n) condition above caters for that. */

  if (all_ngrams && nmatch && !nowild) {
    g_assert(hget(patterns, ngm) == NULL);
    Counter cnt = counter_new(ngm);
    cnt->n0 = ngram_cnt;
    hset(patterns, ngram_copy(ngm), cnt);
  }
}

/** count_ngrams(patterns, pattern_masks) iterates over each data file
 *  in the google ngram directory and calls count_ngram for each line.
 */

void count_ngrams(Hash patterns, Mask **pattern_masks) {
  GString *ngram_file = g_string_new(NULL);
  GString *ngram_files = g_string_new(NULL);
  EmptyPattern ngm;
  guint64 cnt;

  g_string_printf(ngram_files, "|find %s -follow -name \"?gm-????.gz\" -or -name vocab.gz", path);

  foreach_line(fname, ngram_files->str) {
    g_strchomp(fname);
    g_string_printf(ngram_file, "|zcat %s", fname);
    g_message("Reading ngrams from [%s]", ngram_file->str);
    Mask *masks = NULL;
    foreach_line(str, ngram_file->str) {
      char *tab = strchr(str, '\t');
      g_assert(tab != NULL);
      cnt = atoll(tab+1);
      *tab = 0;
      ngram_from_string(ngm, str);
      if (masks == NULL) {
	masks = pattern_masks[ngram_size(ngm)];
	/* If we do not have any masks for this size quit.
	   This assumes a single n in each ngram file. */
	if (masks[0] == NULL) break;
      }
      count_ngram(patterns, ngm, cnt, masks);
    }
  }
  g_string_free(ngram_file, TRUE);
  g_string_free(ngram_files, TRUE);
}

/** main */

int main(int argc, char *argv[]) {
  g_message_init();
  g_message(rcsid);
  get_options(argc, argv);
  g_message("Google ngram path = [%s], all_ngrams = %d", path, all_ngrams);

  g_message("Reading patterns from stdin...");
  Hash patterns = read_patterns();
  g_message("Read %d patterns.", hlen(patterns));
  g_hash_table_analyze(patterns);

  g_message("Creating masks...");
  Mask **masks = mask_table_init(patterns);
  mask_table_dump(masks);

  g_message("Counting ngrams...");
  count_ngrams(patterns, masks);

  g_message("Printing output...");
  pattern_hash_print(patterns);

  g_message("Cleaning up...");
  free_hash(patterns);
  mask_table_free(masks);

  return 0;
}
