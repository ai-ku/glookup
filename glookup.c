/* -*- mode: C; mode: Outline-minor; outline-regexp: "/[*][*]+"; -*- */
/** glookup.c: count patterns in gindex */
const char *rcsid = "$Id$";
 
/** includes */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <getopt.h>
#include <glib.h>
#include "foreach.h"
#include "procinfo.h"
#include "ghashx.h"
 
/** defines */
#define NGRAM_ORDER 5
#define WILDCARD_INT 0
#define WILDCARD_STR "_"
#define ORIGINAL_DATA 1
 
/** type declarations and macros */

/* ngrams and related structures are represented as guint32 arrays
   where the first element of the array gives the length */

typedef guint32 *Pattern, *Ngram, *Mask;
#define plen(p) (*(p))
#define foreach_ptok(var, pat)\
  for (register guint32 _i = 1, _l = plen(pat), var = 0;\
       (_i <= _l) && ((var = pat[_i]) || 1); _i++)

/* We keep track of up to three types of counts.
   n0 is the total number of times matching ngrams occur.
   n1 is the total number of distinct ngrams that match.
   n2 is the hash of distinct last words that match.
   (used in kn smoothing) */

typedef struct CounterS {
  gint64 n0;
  guint32 n1;
  Hash n2;
} *Counter;
 
/** globals */
const gint64 n0vocab = 1024908267229;
const guint32 n1vocab = 13588391;
static gint64 *vocab_cnt;
static guint pattern_hash_rnd[NGRAM_ORDER + 1];
static Mask *pattern_masks[NGRAM_ORDER + 1];
static Hash bignum;
static char *prefix = ".";
static int all_ngrams = 0;
 
/** function declarations */
int main(int argc, char *argv[]);
void get_options(int argc, char *argv[]);
int read_vocabulary(const char *vocab);
Hash read_exception_files();
Hash read_patterns();
void count_all_wildcards(Hash patterns);
void count_single_words(Hash patterns);
void count_single_word(gpointer key, gpointer val, gpointer dat);
gint64 find_bignum(Pattern pat);
void count_ngrams(Hash patterns);
void count_ngram(Hash patterns, Ngram ngm, gint64 ngm_cnt);
void *mapf(const char *path, int *len);
void print_pattern_count(gpointer key, gpointer val, gpointer dat);
guint pattern_hash(gconstpointer key);
Hash pattern_hash_init();
void pattern_from_string(Pattern pat, char *str);
gboolean pattern_equal(gconstpointer a, gconstpointer b);
void free_pattern(gpointer pat);
void free_counter(gpointer cnt);
Pattern make_pattern(char *str);
Pattern copy_pattern(Pattern pat);
void pattern_string(GString *g, Pattern pat);
void pattern_mask_init();
void pattern_mask_add(Pattern pat);
void pattern_mask_apply(Pattern pat, Ngram ngm, Mask msk);
int pattern_mask_match(Pattern pat, Mask msk);
void pattern_mask_dump();
Counter make_counter(Pattern pat);
void print_pattern(Pattern pat);
gboolean n2_pattern(Pattern pat);
int wildcards(Pattern pat);
guint32 str2int(const char *str);
const char *int2str(guint32 tok);
void g_hash_table_analyze(Hash h);
 
/** main */
int main(int argc, char *argv[]) {
  g_message_init();
  g_message("$Id$");
  get_options(argc, argv);
  g_message("Gindex path = [%s], all_ngrams = %d", prefix, all_ngrams);
  g_message("Reading vocabulary...");
#if ORIGINAL_DATA
  char *vocab_file = g_strdup_printf("|zcat %s/1gms/vocab.gz", prefix);
#else  
  char *vocab_file = g_strdup_printf("%s/1gms/vocab.str", prefix);
#endif
  int nvocab = read_vocabulary(vocab_file);
  g_message("%d words found in vocabulary...", nvocab);
#if !ORIGINAL_DATA
  g_message("Reading bignum exception files...");
  bignum = read_exception_files();
  g_message("%d ngrams have bignum counts...", hlen(bignum));
#endif
  g_message("Reading patterns from stdin...");
  Hash patterns = read_patterns();
  g_hash_table_analyze(patterns);
  g_message("Read %d patterns, counting...", hlen(patterns));
  count_all_wildcards(patterns);
  count_single_words(patterns);
  count_ngrams(patterns);
  g_message("Printing output...");
  g_hash_table_foreach(patterns, print_pattern_count, NULL);
  g_message("Cleaning up...");
  free_hash(patterns);
  return 0;
}

/** command line processing */

static struct option long_options[] = {
  {"prefix", 1, 0, 'p'},
  {"all", 0, 0, 'a'},
  {0, 0, 0, 0}
};

const char *help = "glookup [-p prefix] [-a] < patterns > counts";

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
 
/** read_vocabulary */

#if ORIGINAL_DATA

int read_vocabulary(const char *vocab) {
  g_message(vocab);
  int nvocab = 0;
  vocab_cnt = g_malloc(n1vocab * sizeof(gint64));
  foreach_line(str, vocab) {
    char *tab = strchr(str, '\t');
    g_assert(tab != NULL);
    *tab = 0;
    g_assert(++nvocab == g_quark_from_string(str));
    vocab_cnt[nvocab] = atoll(tab+1);
  }
  return nvocab;
}

#else
int read_vocabulary(const char *vocab) {
  g_message(vocab);
  /* Read vocabulary and construct a string-int mapping */
  /* First word is 1, 0 is reserved for wildcard */
  int len = strlen(vocab);
  char *fname;
  if (!strcmp(vocab+len-3, ".gz") ||
      !strcmp(vocab+len-2, ".Z")) {
    /* If it is compressed, put zcat in the beginning */
    /* The '|' as a first character tells foreach.h to open a pipe */
    fname = g_strdup_printf("|zcat %s", vocab);
  } else {
    fname = g_strdup(vocab);
  }
  int nvocab = 0;
  foreach_line(str, fname) {
    int len = strcspn(str, "\t\r\n");
    str[len] = 0;
    g_assert(++nvocab == g_quark_from_string(str));
  }
  if (fname != vocab) free(fname);
  return nvocab;
}
#endif

/** read bignum exceptions */
Hash read_exception_files() {
  Hash hash = g_hash_table_new_full(pattern_hash, pattern_equal, 
				    free_pattern, free_counter);
  char *bignum_exception_files = 
    g_strdup_printf("|find %s -follow -name \"*.exc\"", prefix);
  foreach_line(fname, bignum_exception_files) {
    g_strchomp(fname);
    g_message(fname);
    foreach_line(str, fname) {
      char *tab = strchr(str, '\t');
      *tab = 0;
      Pattern pat = make_pattern(str);
      Counter cnt = make_counter(pat);
      cnt->n0 = atoll(tab+1);
      hset(hash, pat, cnt);
    }
  }
  return hash;
}
 
/** read_patterns */
Hash read_patterns() {
  Hash patterns = pattern_hash_init();
  pattern_mask_init();
  foreach_line(str, NULL) {
    Pattern pat = make_pattern(str);
    if (!hgot(patterns, pat)) {
      Counter cnt = make_counter(pat);
      hset(patterns, pat, cnt);
      if (wildcards(pat) < plen(pat))
	pattern_mask_add(pat);
    } else {
      free_pattern(pat);
    }
  }
  pattern_mask_dump();
  return patterns;
}
 
/** counting patterns */
void count_all_wildcards(Hash patterns) {
  g_message("Counting wildcards");
  const gint64  n0wild[] = {0, 1024908267229, 880784659306, 732213757736, 522447038907, 362967949887};
  const guint32 n1wild[] = {0, 13588391, 314843401, 977069902, 1333820466, 1214460675};
  // const guint32 n2wild[] = {0, 13588391, 12880878, 9329099, 7934066, 7015505};

  guint32 pat[NGRAM_ORDER+1];
  for (int i = 1; i <= NGRAM_ORDER; i++) {
    pat[i] = WILDCARD_INT;
  }
  for (int i = 1; i <= NGRAM_ORDER; i++) {
    plen(pat) = i;
    Counter cnt = hget(patterns, pat);
    if (cnt != NULL) {
      cnt->n0 = n0wild[i];
      cnt->n1 = n1wild[i];
    }
  }
}

#if ORIGINAL_DATA

void count_single_words(Hash patterns) {
  g_message("Counting single words");
  g_hash_table_foreach(patterns, count_single_word, vocab_cnt);
}

#else

void count_single_words(Hash patterns) {
  guint32 *cnt, len;
  char *vocab_cnt_file = 
    g_strdup_printf("%s/1gms/vocab.cnt", prefix);
  cnt = mapf(vocab_cnt_file, &len);
  g_hash_table_foreach(patterns, count_single_word, cnt);
  if (munmap(cnt, len)) g_error("Error in munmap");
}

#endif

void count_single_word(gpointer key, gpointer val, gpointer dat) {
  Pattern pat = key;
  if (plen(pat) != 1) return;
  Counter cnt = val;
  if (pat[1] == WILDCARD_INT) {
    /* these should be set by count_all_wildcards */
    g_assert(cnt->n0 == n0vocab);
    g_assert(cnt->n1 == n1vocab);
  } else if (pat[1] > n1vocab) {
    /* unknown word */
    cnt->n0 = 0;
  } else {
#if ORIGINAL_DATA    
    cnt->n0 = vocab_cnt[pat[1]];
#else
    guint32 *vocab_cnt = dat;
    cnt->n0 = vocab_cnt[pat[1]-1];
    if (cnt->n0 == 0) cnt->n0 = find_bignum(pat);
#endif
  }
}

gint64 find_bignum(Pattern pat) {
  Counter tmp = hget(bignum, pat);
  if (tmp == NULL) g_error("Bignum not found");
  return tmp->n0;
}

#if ORIGINAL_DATA

void count_ngrams(Hash patterns) {
  guint32 ngm[NGRAM_ORDER + 1];
  gint64 cnt;
  char *ngram_files =
    g_strdup_printf("|find %s -follow -name \"?gm-????.gz\"", prefix);
  foreach_line(fname, ngram_files) {
    g_strchomp(fname);
    g_message(fname);
    /* BUG: I don't like this hack: */
    guint32 n = atoi(fname + strlen(fname) - 11);
    g_assert(n >= 2 && n <= 5);
    if (pattern_masks[n][0] == NULL)
      continue;
    foreach_line(str, g_strdup_printf("|zcat %s", fname)) {
      char *tab = strchr(str, '\t');
      *tab = 0;
      pattern_from_string(ngm, str);
      cnt = atoll(tab+1);
      count_ngram(patterns, ngm, cnt);
    }
  }
}

#else

void count_ngrams(Hash patterns) {
  guint32 *map, len, ngm[NGRAM_ORDER + 1];
  gint64 cnt;
  char *binary_ngram_files =
    g_strdup_printf("|find %s -follow -name \"*.ngm-????\"", prefix);
  foreach_line(fname, binary_ngram_files) {
    g_strchomp(fname);
    /* BUG: I don't like this hack: */
    guint32 n = atoi(fname + strlen(fname) - 10);
    g_assert(n >= 2 && n <= 5);
    if (pattern_masks[n][0] == NULL)
      continue;
    map = mapf(fname, &len);
    guint32 imax = len >> 2;
    guint32 iinc = n + 1;
    ngm[0] = n;
    for (guint32 i = 0; i < imax; i += iinc) {
      for (guint32 j = 0; j < n; j++) {
	guint32 tok = map[i+j];
	g_assert(tok >= 1 && tok <= n1vocab);
	ngm[j+1] = tok;
      }
      cnt = map[i+n];
      if (cnt == 0) cnt = find_bignum(ngm);
      count_ngram(patterns, ngm, cnt);
    }
    if (munmap(map, len)) g_error("Error in munmap");
  }
}

#endif
 
void count_ngram(Hash patterns, Ngram ngm, gint64 ngm_cnt) {
  guint32 pat[NGRAM_ORDER + 1];
  guint32 n = plen(ngm);
  Mask *masks = pattern_masks[n];
  guint32 nmatch = 0;

  for (int i = 0; masks[i] != NULL; i++) {
    pattern_mask_apply(pat, ngm, masks[i]);
    Counter cnt = hget(patterns, pat);
    if (cnt != NULL) {
      nmatch++;
      cnt->n0 += ngm_cnt;
      if (wildcards(pat)) cnt->n1++;
      if (cnt->n2 != NULL) hset(cnt->n2, GUINT_TO_POINTER(ngm[n]), 1);
    }
  }

  if (nmatch && all_ngrams) {
    Counter cnt = hget(patterns, ngm);
    if (cnt == NULL) {
      guint32 *copy = copy_pattern(ngm);
      cnt = make_counter(copy);
      cnt->n0 = ngm_cnt;
      hset(patterns, copy, cnt);
    }
  }
}
 
/** printing pattern counts */
void print_pattern_count(gpointer key, gpointer val, gpointer dat) {
  Pattern pat = key;
  Counter cnt = val;
  print_pattern(pat);
  printf("\t%li", cnt->n0);
  if (cnt->n1 > 0) printf("\t%d", cnt->n1);
  if (cnt->n2 != NULL) printf("\t%d", hlen(cnt->n2));
  fputc('\n', stdout);
}
 
void print_pattern(Pattern pat) {
  int n = 0;
  foreach_ptok(pi, pat) {
    if (n++) fputc(' ', stdout);
    fputs(int2str(pi), stdout);
  }
}
 
void pattern_string(GString *g, Pattern pat) {
  int n = 0;
  g_string_truncate(g, 0);
  foreach_ptok(pi, pat) {
    if (n++) g_string_append_c(g, ' ');
    g_string_append(g, int2str(pi));
  }
}

/** Pattern functions */
Pattern make_pattern(char *str) {
  static guint32 pat[NGRAM_ORDER + 1];
  int ntok = 0;
  foreach_token(word, str) {
    pat[++ntok] = str2int(word);
  }
  pat[0] = ntok;
  return copy_pattern(pat);
}

Pattern copy_pattern(Pattern pat) {
  return g_memdup(pat, (plen(pat)+1)*sizeof(guint32));
}

void pattern_from_string(Pattern pat, char *str) {
  int ntok = 0;
  foreach_token(word, str) {
    pat[++ntok] = str2int(word);
  }
  pat[0] = ntok;
}

void free_pattern(gpointer pat) {
  free(pat);
}

guint pattern_hash(gconstpointer key) {
  const guint32 *pat = key;
  guint hash = 0;
  for (int i = plen(pat); i > 0; i--) {
    hash += (1+pat[i]) * pattern_hash_rnd[i];
  }
  return hash;
}

Hash pattern_hash_init() {
  *pattern_hash_rnd = NGRAM_ORDER;
  for (int i = 1; i <= NGRAM_ORDER; i++) {
    pattern_hash_rnd[i] = rand();
    /* g_message("pattern_hash_rnd[%d] = %d", i, pattern_hash_rnd[i]); */
  }
  return g_hash_table_new_full(pattern_hash, pattern_equal, 
			       free_pattern, free_counter);
}

gboolean pattern_equal(gconstpointer a, gconstpointer b) {
  const guint32 *pa = a;
  const guint32 *pb = b;
  if (plen(pa) != plen(pb)) return 0;
  for (int i = 1; i <= plen(pa); i++) {
    if (pa[i] != pb[i]) return 0;
  }
  return 1;
}

gboolean n2_pattern(Pattern pat) {
  if (wildcards(pat) < 2) return 0;
  if (pat[plen(pat)] != WILDCARD_INT) return 0;
  return 1;
}

int wildcards(Pattern pat) {
  int nw = 0;
  foreach_ptok(pi, pat) {
    if (pi == WILDCARD_INT) nw++;
  }
  return nw;
}
 
guint32 str2int(const char *str) {
  guint32 ans = g_quark_try_string(str);
  if ((ans == WILDCARD_INT) && strcmp(str, WILDCARD_STR))
    ans = g_quark_from_string(str);
  return ans;
}

const char *int2str(guint32 tok) {
  const char *str = (tok == WILDCARD_INT) ? WILDCARD_STR :
    g_quark_to_string(tok);
  if (str == NULL) {
    g_message("Integer [%d] not found in vocabulary", tok);
    return g_strdup_printf("[%u]", tok);
  }
  return str;
}
 
/** Mask functions */

void pattern_mask_init() {
  for (int i = 1; i <= NGRAM_ORDER; i++) {
    pattern_masks[i] = g_malloc0((1+(1<<i)) * sizeof(Mask));
  }
}

void pattern_mask_add(Pattern pat) {
  int n = plen(pat);
  Mask *masks = pattern_masks[n];
  int i;
  for (i = 0; masks[i] != NULL; i++) {
    if (pattern_mask_match(pat, masks[i])) return;
  }
  g_assert(i < (1<<n));
  masks[i] = copy_pattern(pat);
}

void pattern_mask_dump() {
  GString *g = g_string_new(NULL);
  g_message("Unique patterns:");
  for (int n = 1; n <= NGRAM_ORDER; n++) {
    Mask *m = pattern_masks[n];
    for (int i = 0; m[i] != NULL; i++) {
      pattern_string(g, m[i]);
      g_message(g->str);
    }
  }
}

int pattern_mask_match(Pattern pat, Mask msk) {
  if (msk[0] != pat[0]) return 0;
  for (int i = pat[0]; i > 0; i--) {
    if (((pat[i] == 0) && (msk[i] != 0)) ||
	((pat[i] != 0) && (msk[i] == 0)))
      return 0;
  }
  return 1;
}

void pattern_mask_apply(Pattern pat, Ngram ngm, Mask msk) {
  pat[0] = ngm[0];
  for (int i = ngm[0]; i > 0; i--) {
    pat[i] = msk[i] ? ngm[i] : 0;
  }
}

/** Counter functions */
Counter make_counter(Pattern pat) {
  Counter cnt = g_new0(struct CounterS, 1);
  if (n2_pattern(pat)) 
    cnt->n2 = g_hash_table_new(g_direct_hash, g_direct_equal);
  return cnt;
}

void free_counter(gpointer val) {
  Counter cnt = val;
  free(cnt->n2);
  free(cnt);
}
 
/** file i/o */
/* BUG: get this out of this file into its own header */
/* Also make procinfo.h self contained, without a need for procinfo.c */

void *mapf(const char *path, int *len) {
  g_message(path);
  int fd = open(path, O_RDONLY);
  if (fd == -1) g_error("Cannot open file %s", path);
  struct stat buf;
  int rval = fstat(fd, &buf);
  if (rval == -1) g_error("Cannot stat file %s", path);
  *len = buf.st_size;
  void *map = mmap(0, *len, PROT_READ, MAP_SHARED, fd, 0);
  if (map == MAP_FAILED) g_error("Cannot mmap file %s", path);
  return map;
}
 
/** debugging: move this to ghashx */

typedef struct _GHashNode      GHashNode;

struct _GHashNode
{
  gpointer   key;
  gpointer   value;
  GHashNode *next;
};

struct _GHashTable
{
  gint             size;
  gint             nnodes;
  GHashNode      **nodes;
  GHashFunc        hash_func;
  GEqualFunc       key_equal_func;
  GDestroyNotify   key_destroy_func;
  GDestroyNotify   value_destroy_func;
};

void g_hash_table_analyze(Hash h) {
  struct _GHashTable *ht = h;
  int sum=0, zero=0, max=0;
  for (int i = 0; i < ht->size; i++) {
    GHashNode **node = &ht->nodes[i];
    int cnt = 0;
    while (*node) {
      cnt++; node = &(*node)->next;
    }
    sum += cnt;
    if (cnt > max) max = cnt;
    if (cnt == 0) zero++;
  }
  g_message("size=%d nelem=%d full=%d avg=%g max=%d",
	    ht->size, sum, ht->size - zero, (float)sum/ht->size, max);
}
