/* -*- mode: C; mode: Outline-minor; outline-regexp: "/[*][*]+"; -*- */
const char *help = "glookup [-a20] [-p google-ngram-path] < patterns > counts";

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include "dlib.h"

/** Command line processing */

static char *path = ".";	/* google-ngram path */
static int all_ngrams = 0;	/* whether to print all ngrams matching a wildcard pattern */
static int output_n2 = 0; 	/* whether to compute n2 (by default we do not it takes memory) */
static int output_zero = 0;	/* whether to print out ngrams with zero count (default no) */

static struct option long_options[] = {
  {"path", 1, 0, 'p'},
  {"all", 0, 0, 'a'},
  {"output_n2", 0, 0, '2'},
  {"output_zero", 0, 0, '0'},
  {0, 0, 0, 0}
};

void get_options(int argc, char *argv[]) {
  while(1) {
    int option_index = 0;
    int c = getopt_long_only(argc, argv, "p:a20", long_options, &option_index);
    if (c == -1) break;
    switch(c) {
    case 'p': path = optarg; break;
    case 'a': all_ngrams = 1; break;
    case '2': output_n2 = 1; break;
    case '0': output_zero = 1; break;
    default: msg("%s", help); exit(0);
    }
  }
}

/** Word tokens are represented as u32 integers (sym_t).  The string
    "_" and its associated sym_t represent the wildcard token. */

typedef sym_t Token;
static Token ANY;		/* sym_t corresponding to wildcard */
#define WILDCARD_STR "_"	/* String corresponding to wildcard */
void wildcard_init() { ANY = str2sym(WILDCARD_STR, true); }

/** Ngrams and related structures are represented as u32 arrays
 * where the first element of the array gives the length.  The other
 * elements are positive integers representing tokens, special token
 * ANY represents a wildcard.  The three types (Ngram, Pattern, Mask)
 * differ only by intent of usage: Ngram represents ngrams read from
 * data, Pattern represents patterns with wildcards that the user
 * supplies, and Mask represents different types of pattern determined
 * by the position of wildcards.
 */

#define MAX_NGRAM_ORDER 5
typedef Token *Ngram, *Pattern, *Mask;
typedef Token EmptyPattern[MAX_NGRAM_ORDER + 1];

#define ngram_size(p) (*(p))

#define forngram(var, ngm)\
  for (register u32 _i = 1, _l = ngram_size(ngm), var = 0;\
       (_i <= _l) && ((var = ngm[_i]) || 1); _i++)

Ngram ngram_copy(Ngram p) {
  size_t len = sizeof(Token) * (1 + ngram_size(p));
  Ngram q = dalloc(len);
  memcpy(q, p, len);
  return q;
}

void ngram_from_string(Ngram ngm, char *str, bool create_new_symbols) {
  int ntok = 0;
  fortok (word, str) {
    assert(ntok < MAX_NGRAM_ORDER);
    ngm[++ntok] = str2sym(word, create_new_symbols);
  }
  ngm[0] = ntok;
}

str_t ngram_to_string(str_t *s, size_t *n, Ngram pat) {
  if (*s == NULL || *n == 0) {
    *n = 120;
    *s = _d_malloc(*n);
  }
  str_t str = *s;
  *str = '\0';
  size_t nstr = 0;
  for (u32 i = 1; i <= ngram_size(pat); i++) {
    str_t wi = sym2str(pat[i]);
    size_t ni = strlen(wi);
    while (nstr + ni + 2 > *n) {
      *n <<= 1;
      *s = _d_realloc(*s, *n);
      str = *s;
    }
    if (nstr > 0) str[nstr++] = ' ';
    strcpy(&str[nstr], wi);
    nstr += ni;
  }
  return str;
}

int ngram_count_wildcards(Ngram pat) {
  int nw = 0;
  forngram(pi, pat) {
    if (pi == ANY) nw++;
  }
  return nw;
}

void ngram_print_count(Ngram ngm, u64 n0) {
  static char *str = NULL;
  static size_t nstr = 0;
  printf("%s\t%lu\n", ngram_to_string(&str, &nstr, ngm), n0);
}


/** Counters: we keep track of up to three types of counts.
 * n0 is the total number of times matching ngrams occur.
 * n1 is the total number of distinct ngrams that match.
 * n2 is the hash of distinct last words that match.
 * (used in kn smoothing) */

/** Patterns with no wildcards have NULL counters, their counts do not
 *  need to be stored in memory and are directly printed out.  Only
 *  patterns that have two or more wildcards, one of which is the last
 *  token need n2 hash.
 */

D_HASH(thash_, Token, Token, d_eqmatch, d_ident, d_ident, d_ident, d_iszero, d_mkzero)

typedef struct CounterS {
  u64 n0;
  u32 n1;
  darr_t n2;
} *Counter;

Counter counter_new(Pattern pat) {
  Counter cnt = _d_calloc(1, sizeof(struct CounterS));
  if (output_n2) {
    int nwild = ngram_count_wildcards(pat);
    if ((nwild > 1) &&
	(pat[ngram_size(pat)] == ANY))
      cnt->n2 = thash_new(0);
  }
  return cnt;
}

void counter_free(Counter cnt) {
  if (cnt != NULL) {
    if (cnt->n2 != NULL) thash_free(cnt->n2);
    _d_free(cnt);
  }
}


/** Pattern hash is initialized with the patterns from the user and
    empty counters.  The counters are updated when going over the
    ngram files. */

u64 pattern_hash(const Pattern pat) {
  u64 hash = 14695981039346656037ULL;
  for (int i = ngram_size(pat); i >= 0; i--) {
    hash ^= pat[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

bool pattern_equal(Pattern pa, Pattern pb) {
  if (ngram_size(pa) != ngram_size(pb))
    return 0;
  for (int i = ngram_size(pa); i > 0; i--) {
    if (pa[i] != pb[i]) return 0;
  }
  return 1;
}

typedef struct PC_s { Pattern key; Counter val; } *PC_t;
#define _pc_init(p) ((struct PC_s) { ngram_copy(p), counter_new(p) })
D_HASH(phash_, struct PC_s, Pattern, pattern_equal, pattern_hash, d_keyof, _pc_init, d_keyisnull, d_keymknull)

/* To print the pattern counts after going through google ngrams */

void pattern_hash_print (darr_t patterns) {
  char *str = NULL;
  size_t nstr = 0;
  forhash (PC_t, pc, patterns, d_keyisnull) {
    Pattern pat = pc->key;
    Counter cnt = pc->val;
    assert (cnt != NULL);
    /* By default, do not print patterns with zero count */
    if (!output_zero) {
      if (cnt->n0 == 0) continue;
    }
    printf("%s\t%lu", ngram_to_string(&str, &nstr, pat), cnt->n0);
    int nwild = ngram_count_wildcards(pat);
    if (nwild > 0) {
      printf("\t%u", cnt->n1);
      if (cnt->n2 != NULL) 
	printf("\t%llu", len(cnt->n2));
    }
    fputc('\n', stdout);
  }
  free(str);
}


/** Masks represent types of patterns based on their wildcard
 * positions.  They are important for efficiency: while going over a
 * 5-gram data file one would need to do 32 pattern hash lookups for
 * each ngram to cover all wildcard combinations.  With masks we only
 * check the actual types that the user asked for. */

/* mask_match(msk, pat) returns true if the wildcard positions in pat
   matches ones in msk */

bool mask_match(Mask msk, Pattern pat) {
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
  assert(ngram_size(msk) == ngram_size(ngm));
  ngram_size(pat) = ngram_size(ngm);
  for (int i = ngram_size(ngm); i > 0; i--) {
    pat[i] = ((msk[i] == ANY) ? msk[i] : ngm[i]);
  }
}

#define foreach_mask(m, a)\
  for (register u32 *m = NULL, _i = 0; (m=a[_i])!=NULL; _i++)

/* mask_table_add(pat, cnt, table) adds a mask for the pattern pat to
   the table if it does not already exist.  The cnt variable is the
   value in the hash table iterated over and is ignored. */

void mask_table_add(Mask **table, darr_t patterns) {
  forhash (PC_t, pc, patterns, d_keyisnull) {
    Pattern pat = pc->key;
    int n = ngram_size(pat);
    Mask *masks = table[n];
    int mask_count = 0;
    foreach_mask(m, masks) {
      if (mask_match(m, pat)) {
	mask_count = -1; break;
      }
      mask_count++;
    }
    if (mask_count == -1) continue;
    assert(mask_count < (1<<n));
    masks[mask_count] = ngram_copy(pat);
  }
}

/* The mask_table[n] is an array of masks for patterns of length n
   given by the user */

Mask **mask_table_init(darr_t patterns) {
  Mask **mask_table = _d_calloc((1 + MAX_NGRAM_ORDER), sizeof(Mask *));
  for (int n = MAX_NGRAM_ORDER; n > 0; n--) {
    mask_table[n] = _d_calloc((1+(1<<n)), sizeof(Mask));
  }
  mask_table_add(mask_table, patterns);
  return mask_table;
}

void mask_table_free(Mask **mask_table) {
  for (int n = MAX_NGRAM_ORDER; n > 0; n--) {
    _d_free(mask_table[n]);
  }
  _d_free(mask_table);
}

void mask_table_dump(Mask **mask_table) {
  char *str = NULL;
  size_t nstr = 0;
  for (int n = 1; n <= MAX_NGRAM_ORDER; n++) {
    Mask *masks = mask_table[n];
    foreach_mask(m, masks) {
      fprintf(stderr, "%s\n", ngram_to_string(&str, &nstr, m));
    }
  }
  free(str);
}


/** read_patterns(): Reads ngram patterns from stdin and returns a
 * hash table where keys are the patterns and the values are empty
 * counters. 
 */

darr_t read_patterns() {
  EmptyPattern pat;
  wildcard_init();
  darr_t patterns = phash_new(0);
  forline (str, NULL) {
    ngram_from_string(pat, str, true);
    phash_get(patterns, pat, true);
  }
  return patterns;
}

/** count_ngram(patterns, ngm, ngram_cnt, lastword, masks) increments
 *  the appropriate counts in the patterns hash looking at ngram ngm
 *  and its count ngram_cnt.  The parameter masks is an array of masks
 *  for the ngram size of ngm.
 */

void count_ngram(darr_t patterns, Ngram ngm, u64 ngram_cnt, Mask *masks) {
  EmptyPattern pat;
  PC_t pc;
  Counter cnt;
  bool nowild = false;	/* whether this ngram matches exactly */
  bool nmatch = false;	/* whether this ngram matches any patterns
			   (except all wildcards) */

  foreach_mask(m, masks) {
    mask_apply(m, ngm, pat);
    if ((pc = phash_get(patterns, pat, false)) != NULL) {
      int nwild = ngram_count_wildcards(pat);
      if (nwild == 0) nowild = true;
      cnt = pc->val;
      int n = ngram_size(ngm);
      if (nwild < n) nmatch = true;
      cnt->n0 += ngram_cnt;
      cnt->n1++;
      if (cnt->n2 != NULL) {
	Token q = ngm[n];
	thash_get(cnt->n2, q, true);
      }
    }
  }

  /* if we have a match for this ngram (nmatch) and (1) the user wants
     us to print all matching ngrams as well (all_ngrams), (2) the
     ngram itself is not in the pattern list (!nowild), we print the
     ngram with its count right away.  The only exception is when the
     only pattern the ngram matched is an all wild pattern, in which
     case the user probably did not want us to dump the whole
     database.  The (nwild < n) condition above caters for that. */

  if (all_ngrams && nmatch && !nowild) {
    assert(phash_get(patterns, ngm, false) == NULL);
    ngram_print_count(ngm, ngram_cnt);
  }
}

/** count_ngrams(patterns, pattern_masks) iterates over each data file
 *  in the google ngram directory and calls count_ngram for each line.
 */

void count_ngrams(darr_t patterns, Mask **pattern_masks) {
  char *ngram_files = malloc(strlen(path) + 60);
  EmptyPattern ngm;
  u64 cnt;

  sprintf(ngram_files, "<find %s -follow -name \"?gm-????.gz\" -or -name vocab.gz", path);
  forline(fname, ngram_files) {
    fname[strlen(fname) - 1] = '\0';
    msg("Reading ngrams from [%s]", fname);
    Mask *masks = NULL;
    forline (str, fname) {
      char *tab = strchr(str, '\t');
      assert(tab != NULL);
      cnt = atoll(tab+1);
      *tab = 0;
      ngram_from_string(ngm, str, false);
      if (masks == NULL) {
	masks = pattern_masks[ngram_size(ngm)];
	/* If we do not have any masks for this size quit.
	   This assumes a single n in each ngram file. */
	if (masks[0] == NULL) break;
      }
      count_ngram(patterns, ngm, cnt, masks);
    }
  }
  free(ngram_files);
}

/** main */

int main(int argc, char *argv[]) {
  get_options(argc, argv);
  msg("Google ngram path = [%s], all_ngrams = %d", path, all_ngrams);

  msg("Reading patterns from stdin...");
  darr_t patterns = read_patterns();
  msg("Read %d patterns.", len(patterns));

  msg("Creating masks...");
  Mask **masks = mask_table_init(patterns);
  mask_table_dump(masks);

  msg("Counting ngrams...");
  count_ngrams(patterns, masks);

  msg("Printing output...");
  pattern_hash_print(patterns);

  msg("Cleaning up...");
  phash_free(patterns);
  mask_table_free(masks);

  msg("done");
  return 0;
}
