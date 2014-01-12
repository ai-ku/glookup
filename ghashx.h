/* $Id: ghashx.h,v 1.1 2014/01/12 15:32:10 dyuret Exp $ */

#ifndef __GHASHX_H__
#define __GHASHX_H__

#include <glib.h>

/* shorthands */
typedef GHashTable *Hash;
#define hset(h,k,v) g_hash_table_insert((h),(k),(gpointer)(v))
#define hget(h,k) g_hash_table_lookup((h),(k))
#define hgot(h,k) g_hash_table_lookup_extended((h),(k),NULL,NULL)
#define hdel(h,k) g_hash_table_remove((h),(k))
#define hlen(h) g_hash_table_size(h)
#define hmap(h,f,d) g_hash_table_foreach((h),(GHFunc)(f),(d))
#define free_hash(h) g_hash_table_destroy(h)
#define hfree(h) g_hash_table_destroy(h)

/* string-to-int hash functions */
#define I2P(v) GINT_TO_POINTER(v)
#define s2i_new() g_hash_table_new_full(g_str_hash, g_str_equal, g_free, NULL)
#define s2i_get(h,k) GPOINTER_TO_INT(hget(h,k))
#define s2i_set(h,k,v) (hgot(h,k) ? hset(h,k,I2P(v)) : hset(h,g_strdup(k),I2P(v)))
#define s2i_add(h,k,v) (hgot(h,k) ? hset(h,k,I2P(s2i_get(h,k)+(v))) : hset(h,g_strdup(k),I2P(v)))
#define s2i_inc(h,k) s2i_add(h,k,1)
#define s2i_free(h) hfree(h)

#endif
