/* -*- mode: C; mode: Outline-minor; outline-regexp: "/[*][*]+"; -*- */
const char *rcsid = "$Id: glookup.c,v 1.1 2007/06/01 10:40:32 dyuret Exp dyuret $";
const char *help = "glookup [-p prefix] [-a] < patterns > counts";
 
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <libgen.h>
#include <glib.h>
#include "foreach.h"
#include "procinfo.h"
#include "ghashx.h"

 
/** Command line processing */

static char *prefix = ".";
static int all_ngrams = 0;
 
static struct option long_options[] = {
  {"prefix", 1, 0, 'p'},
  {"all", 0, 0, 'a'},
  {0, 0, 0, 0}
};

void get_options(int argc, char *argv[]) {
  while(1) {
    int option_index = 0;
    int c = getopt_long_only(argc, argv, "p:a", long_options, &option_index);
    if (c == -1) break;
    switch(c) {
    case 'p': prefix = optarg; break;
    case 'a': all_ngrams = 1; break;
    default: g_error(help); break;
    }
  }
}
 
/** Word tokens are represented as guint32 integers.  The integer 0
 * and the string "_" represent the wildcard token. */
 
typedef guint32 Token;
#define WILDCARD_TOK 0
#define WILDCARD_STR "_"
#define str2tok(s) (strcmp((s), WILDCARD_STR) ? g_quark_from_string(s) : WILDCARD_TOK)
#define tok2str(t) (((t) == WILDCARD_TOK) ? WILDCARD_STR : g_quark_to_string(t))


/** read_vocabulary() reads vocab.gz, initializes the tokens
 * (represented with GQuarks), and returns the counts.  Note that this
 * code sort of depends on WILDCARD being zero and GQuarks being
 * consecutive positive numbers. */

GArray *read_vocabulary(const char *path) {
  GArray *cnt = g_array_new(FALSE, FALSE, sizeof(guint64));
  guint64 total = 0;
  g_array_append_val(cnt, total); /* reserved for the total count */
  foreach_line(str, path) {
    char *tab = strchr(str, '\t');
    g_assert(tab != NULL);
    guint64 n = atoll(tab+1);
    g_assert(n > 0);
    total += n;
    g_array_append_val(cnt, n);
    *tab = 0;
    GQuark q = g_quark_from_string(str);
    g_assert(q == cnt->len - 1);
  }
  g_array_index(cnt, guint64, 0) = total;
  return cnt;
}
 
/** Ngrams and related structures are represented as guint32 arrays
 * where the first element of the array gives the length.  The other
 * elements are positive integers representing tokens or zero which
 * represents a wildcard.  The three types differ only by intent of
 * usage: Ngram represents ngrams read from data, Pattern represents
 * patterns with wildcards that the user supplies, and Mask represents
 * different types of pattern determined by the position of wildcards.
 */

#define NGRAM_ORDER 5
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
    ngm[++ntok] = str2tok(word);
  }
  ngm[0] = ntok;
}

void ngram_to_gstring(GString *g, Ngram pat) {
  int n = 0;
  g_string_truncate(g, 0);
  foreach_tok(pi, pat) {
    if (n++) g_string_append_c(g, ' ');
    g_string_append(g, tok2str(pi));
  }
}

int ngram_count_wildcards(Ngram pat) {
  int nw = 0;
  foreach_tok(pi, pat) {
    if (pi == WILDCARD_TOK) nw++;
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
  if ((ngram_count_wildcards(pat) >= 2) &&
      (pat[ngram_size(pat)] == WILDCARD_TOK))
    cnt->n2 = g_hash_table_new(g_direct_hash, g_direct_equal);
  return cnt;
}

void counter_free(gpointer val) {
  Counter cnt = val;
  if (cnt->n2 != NULL) free(cnt->n2);
  free(cnt);
}
 
 
/** Pattern hash is initialized with the patterns from the user and
    empty counters.  The counters are updated when going over the
    ngram files. */

static guint pattern_hash_rnd[NGRAM_ORDER + 1];

guint pattern_hash(gconstpointer key) {
  const guint32 *pat = key;
  guint hash = 0;
  for (int i = ngram_size(pat); i > 0; i--) {
    hash += (1+pat[i]) * pattern_hash_rnd[i];
  }
  return hash;
}

gboolean pattern_equal(gconstpointer a, gconstpointer b) {
  const guint32 *pa = a;
  const guint32 *pb = b;
  if (ngram_size(pa) != ngram_size(pb)) 
    return 0;
  for (int i = ngram_size(pa); i > 0; i--) {
    if (pa[i] != pb[i]) return 0;
  }
  return 1;
}

Hash pattern_hash_init() {
  for (int i = 1; i <= NGRAM_ORDER; i++) {
    pattern_hash_rnd[i] = rand();
  }
  return g_hash_table_new_full(pattern_hash, pattern_equal, 
			       ngram_free, counter_free);
}
 
void pattern_hash_print_fn(Pattern pat, Counter cnt, GString *gstr) {
  ngram_to_gstring(gstr, pat);
  fputs(gstr->str, stdout);
  fputc('\t', stdout);
  printf("%" G_GUINT64_FORMAT, cnt->n0);
  if (cnt->n1 > 0) printf("\t%u", cnt->n1);
  if (cnt->n2 != NULL) printf("\t%u", hlen(cnt->n2));
  fputc('\n', stdout);
}
 
void pattern_hash_print(Hash patterns) {
  GString *gstr = g_string_new(NULL);
  hmap(patterns, pattern_hash_print_fn, gstr);
  g_string_free(gstr, TRUE);
}

 
/** read_patterns(): Reads ngram patterns from stdin and returns a
 * hash table where keys are the patterns and the values are empty
 * counters. */

Hash read_patterns() {
  EmptyPattern pat;
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

 
/** Masks represent types of patterns based on their wildcard
 * positions.  They are important for efficiency: while going over a
 * 5-gram data file one would need to do 32 pattern hash lookups to
 * cover all wildcard combinations.  With masks we only check the
 * actual types that the user asked for. */

gboolean mask_match(Mask msk, Pattern pat) {
  if (ngram_size(msk) != ngram_size(pat)) return 0;
  for (int i = ngram_size(msk); i > 0; i--) {
    if ((pat[i] == WILDCARD_TOK) != (msk[i] == WILDCARD_TOK))
      return 0;
  }
  return 1;
}

void mask_apply(Mask msk, Ngram ngm, Pattern pat) {
  g_assert(ngram_size(msk) == ngram_size(ngm));
  ngram_size(pat) = ngram_size(ngm);
  for (int i = ngram_size(ngm); i > 0; i--) {
    pat[i] = ((msk[i] == WILDCARD_TOK) ? msk[i] : ngm[i]);
  }
}

void mask_table_add(Pattern pat, Counter cnt, Mask **table) {
  int n = ngram_size(pat);
  Mask *masks = table[n];
  int i;
  for (i = 0; masks[i] != NULL; i++) {
    if (mask_match(masks[i], pat)) return;
  }
  g_assert(i < (1<<n));
#if WILD_CACHE
  if (ngram_count_wildcards(pat) < ngram_size(pat))
#endif
  masks[i] = ngram_copy(pat);
}

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

void mask_table_dump(Mask **mask_table) {		/* for debugging */
  GString *g = g_string_new(NULL);
  for (int n = 1; n <= NGRAM_ORDER; n++) {
    Mask *m = mask_table[n];
    for (int i = 0; m[i] != NULL; i++) {
      ngram_to_gstring(g, m[i]);
      g_message(g->str);
    }
  }
  g_string_free(g, TRUE);
}
  
 
/** counting patterns */

#if WILD_CACHE

/* This is an optimization that prevents all-wildcard patterns that
   match every ngram from being calculated using the data */

const guint64 n0wild[] = {0, 1024908267229, 880784659306, 732213757736, 522447038907, 362967949887};
const guint32 n1wild[] = {0, 13588391, 314843401, 977069902, 1333820466, 1214460675};
const guint32 n2wild[] = {0, 13588391, 12880878, 9329099, 7934066, 7015505};

void count_all_wildcards(Hash patterns) {
  EmptyPattern pat;
  for (int i = NGRAM_ORDER; i > 0; i--) {
    pat[i] = WILDCARD_TOK;
  }
  for (int i = NGRAM_ORDER; i > 0; i--) {
    ngram_size(pat) = i;
    Counter cnt = hget(patterns, pat);
    if (cnt != NULL) {
      cnt->n0 = n0wild[i];
      cnt->n1 = n1wild[i];
    }
  }
}

#endif

void count_single_word(Ngram pat, Counter cnt, GArray *vocab_cnt) {
  if (ngram_size(pat) != 1) return;
  if (pat[1] < vocab_cnt->len) {
    cnt->n0 = g_array_index(vocab_cnt, guint64, pat[1]);
  } else { /* unknown word */
    cnt->n0 = 0;
  }
}

void count_ngram(Hash patterns, Ngram ngm, guint64 ngram_cnt, Mask *masks) {
  EmptyPattern pat;
  int n = ngram_size(ngm);
  int nmatch = 0;
  int nowild = 0;

  for (int i = 0; masks[i] != NULL; i++) {
    mask_apply(masks[i], ngm, pat);
    Counter cnt = hget(patterns, pat);
    if (cnt != NULL) {
      nmatch++;
      cnt->n0 += ngram_cnt;
      if (ngram_count_wildcards(pat)) { 
	cnt->n1++;
      } else {
	nowild = 1;
      }
      if (cnt->n2 != NULL) hset(cnt->n2, GUINT_TO_POINTER(ngm[n]), 1);
    }
  }

  if (all_ngrams && nmatch && !nowild) {
    g_assert(hget(patterns, ngm) == NULL);
    Counter cnt = counter_new(ngm);
    cnt->n0 = ngram_cnt;
    hset(patterns, ngram_copy(ngm), cnt);
  }
}

void count_ngrams(Hash patterns, Mask **pattern_masks) {
  GString *gstr = g_string_new(NULL);
  EmptyPattern ngm;
  guint64 cnt;

  char *ngram_files =
    g_strdup_printf("|find %s -follow -name \"?gm-????.gz\"", prefix);
  
  foreach_line(fname, ngram_files) {
    g_strchomp(fname);
    g_string_printf(gstr, "|zcat %s", fname);
    /* BUG: I don't like this hack: */
    guint32 n = atoi(basename(fname));
    g_assert(n >= 2 && n <= 5);
    if (pattern_masks[n][0] == NULL) continue;
    g_message("Reading %d-grams from [%s]", n, gstr->str);
    foreach_line(str, gstr->str) {
      char *tab = strchr(str, '\t');
      g_assert(tab != NULL);
      cnt = atoll(tab+1);
      *tab = 0;
      ngram_from_string(ngm, str);
      count_ngram(patterns, ngm, cnt, pattern_masks[n]);
    }
  }
  g_string_free(gstr, TRUE);
}

/** main */
int main(int argc, char *argv[]) {
  g_message_init();
  g_message(rcsid);
  get_options(argc, argv);
  g_message("Gindex path = [%s], all_ngrams = %d", prefix, all_ngrams);
  char *vocab_path = g_strdup_printf("|zcat %s/1gms/vocab.gz", prefix);
  g_message("Reading [%s]...", vocab_path);
  GArray *vocab_cnt = read_vocabulary(vocab_path);
  g_message("%d words found in vocabulary, total count = %llu...", 
	    vocab_cnt->len - 1, g_array_index(vocab_cnt, guint64, 0));
  g_message("Reading patterns from stdin...");
  Hash patterns = read_patterns();
  g_message("Read %d patterns.", hlen(patterns));
  g_hash_table_analyze(patterns);
#if WILD_CACHE
  g_assert(g_array_index(vocab_cnt, guint64, 0) == n0wild[1]);
  g_assert(vocab_cnt->len == n1wild[1] + 1);
  g_message("Counting all wildcards...");
  count_all_wildcards(patterns);
#endif
  g_message("Counting single words...");
  hmap(patterns, count_single_word, vocab_cnt);
  g_message("Creating masks...");
  Mask **masks = mask_table_init(patterns);
  mask_table_dump(masks);
  count_ngrams(patterns, masks);
  g_message("Printing output...");
  pattern_hash_print(patterns);
  g_message("Cleaning up...");
  free(vocab_path);
  g_array_free(vocab_cnt, TRUE);
  free_hash(patterns); 
  mask_table_free(masks);
  return 0;
}
