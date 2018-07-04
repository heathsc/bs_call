//
//  dbSNP_idx.c
//  dbSNP
//
//  Created by Simon Heath on 15/11/2017.
//  Copyright © 2017 Simon Heath. All rights reserved.
//

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>

static char *header;
static contig *contigs;
static prefix *prefixes;
static uint16_t n_prefix;

static void clear_bins(bin *b, int ct) {
  for(int i = 0; i < ct; i++) {
    b[i].entry_size = 0;
    b[i].entry = NULL;
    b[i].name_buf_size = b[i].name_buf_idx = 0;
    b[i].name_buf = NULL;
    b[i].mask = 0;
    b[i].n_entries = 0;
  }
}

static char dtab[256] =
  {
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0x0, 0x1, 0x2, 0x3, 0x4, 0x5, 0x6, 0x7, 0x8, 0x9, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf,
   0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf, 0xf
  };

static unsigned char dtab2[256] =
  {
   33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 0, 0, 0, 0, 0, 133,
   43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 0, 0, 0, 0, 0, 134,
   53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 0, 0, 0, 0, 0, 135,
   63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 0, 0, 0, 0, 0, 136,
   73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 0, 0, 0, 0, 0, 137,
   83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 0, 0, 0, 0, 0, 138,
   93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 0, 0, 0, 0, 0, 139,
   103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 0, 0, 0, 0, 0, 140,
   113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 0, 0, 0, 0, 0, 141,
   123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 0, 0, 0, 0, 0, 142,
  };

static void add_to_bin(bin *b, uint8_t off, char *name) {
  uint64_t msk = (uint64_t)1 << off;
  if(b->mask & msk) return;
  b->mask |= msk;
  if(!b->entry_size) {
    b->entry_size = INITIAL_ENTRY_SIZE;
    b->entry = malloc(sizeof(uint16_t) * INITIAL_ENTRY_SIZE * 2);
    b->name_buf_size = INITIAL_NAME_BUF_SIZE;
    b->name_buf = malloc((size_t)INITIAL_NAME_BUF_SIZE);
  }
  if(b->n_entries == b->entry_size) {
    b->entry_size *= 1.5;
    if(b->entry_size > 64) b->entry_size = 64;
    b->entry = realloc(b->entry, sizeof(uint16_t) * b->entry_size * 2);
  }
  size_t l = strlen(name);
  size_t k = l;
  for(; k > 0; k--) if(name[k - 1] < '0' || name[k-1] > '9') break;
  prefix *p = NULL;
   
  HASH_FIND(hh, prefixes, name, k, p);
  if(p != NULL) p->count++;
  else {
    if(n_prefix == 0xffff) {
      fprintf(stderr, "Too many SNP prefixes found (%d)\n", n_prefix);
      exit(-1);
    }
    p = malloc(sizeof(prefix));
    p->ix = n_prefix++;
    p->count = 1;
    p->pref = malloc((size_t)(k+1));
    memcpy(p->pref, name, k);
    p->pref[k] = 0;
    HASH_ADD_KEYPTR(hh, prefixes, p->pref, k, p);
    fprintf(stderr,"Added new prefix %d, '%s'\n", p->ix, p->pref);
  }
   
  //   fprintf(stderr,"name %s, prefix %.*s, suffix %s\n", name, (int)k, name, name + k);
  name += k;
  l -= k;
   
  if(l > 510) {
    fprintf(stderr,"%s:%d %s() Marker name exceeds 510 characters: %s\n", __FILE__, __LINE__, __func__, name);
    exit(-1);
  }
  uint16_t l1 = (uint16_t)(l + 1) >> 1;
  uint16_t sz = b->name_buf_idx;
  if(sz + l1 > b->name_buf_size) {
    b->name_buf_size = (b->name_buf_size + l1) * 1.25;
    b->name_buf = realloc(b->name_buf, (size_t)b->name_buf_size);
  }
  int x = (b->n_entries++) << 1;
  b->entry[x] = (l1 << 8) | off;
  b->entry[x + 1] = p->ix;
  for(size_t i = 0; i < l; i += 2) {
    b->name_buf[b->name_buf_idx++] = (dtab[(int)name[i]] << 4) | dtab[(int)name[i + 1]];
  }
}

static bool process_file(char *s) {
  bool error = false;
  FILE *fptr;
  bool filter = false;
  if(s == NULL) fptr = stdin;
  else fptr = open_readfile(s, &filter);
  if(fptr == NULL) return true;
  char *buf = NULL;
  size_t buf_size = 0;
  tokens *tok = NULL;
  char *curr_ctg = NULL;
  contig *ctg = NULL;
  while(!error) {
    ssize_t l = getline(&buf, &buf_size, fptr);
    if(l < 0) break;
    if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
    if(l > 0) {
      if(l > 6 && !strncmp(buf, "track ", 6)) {
	if(header == NULL) header = strdup(buf);
      } else {
	tok = tokenize(buf, '\t', tok);
	if(tok->n_tok >= 4) {
	  char *p;
	  uint32_t x = (uint32_t)strtoul(tok->toks[1], &p, 10);
	  uint32_t y = (uint32_t)strtoul(tok->toks[2], &p, 10);
	  if(y > x && y - x == 1) {
	    uint32_t binx = y >> 6;
	    if(curr_ctg == NULL || strcmp(curr_ctg, tok->toks[0])) {
	      if(curr_ctg != NULL) free(curr_ctg);
	      curr_ctg = strdup(tok->toks[0]);
	      HASH_FIND_STR(contigs, curr_ctg, ctg);
	      if(ctg == NULL) {
		ctg = malloc(sizeof(contig));
		ctg->min_bin = ctg->max_bin = binx;
		ctg->bins = malloc(sizeof(bin));
		ctg->name = strdup(curr_ctg);
		size_t sz = strlen(ctg->name);
		HASH_ADD_KEYPTR(hh, contigs, ctg->name, sz, ctg);
		clear_bins(ctg->bins, 1);
		//                        printf("New contig %s\n", curr_ctg);
	      }
	    }
	    if(binx > ctg->max_bin) {
	      ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(binx - ctg->min_bin + 1));
	      clear_bins(ctg->bins + (ctg->max_bin - ctg->min_bin + 1), binx - ctg->max_bin);
	      //                     printf("Expanding bin range for %s from %d - %d to %d - %d\n", curr_ctg, ctg->min_bin, ctg->max_bin, ctg->min_bin, binx);
	      ctg->max_bin = binx;
	    } else if(binx < ctg->min_bin) {
	      ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - binx + 1));
	      memmove(ctg->bins + (ctg->min_bin - binx), ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - ctg->min_bin + 1));
	      clear_bins(ctg->bins, ctg->min_bin - binx);
	      //                   printf("Expanding bin range for %s from %d - %d to %d - %d\n", curr_ctg, ctg->min_bin, ctg->max_bin, binx, ctg->max_bin);
	      ctg->min_bin = binx;
	    }
	    add_to_bin(ctg->bins + (binx - ctg->min_bin), y & 63, tok->toks[3]);
	  }
	}
      }
    }
  }
  if(curr_ctg != NULL) free(curr_ctg);
  if(buf != NULL) free(buf);
  if(tok != NULL) free_tokens(tok);
  if(s != NULL) {
    fclose(fptr);
    if(filter) {
      int i;
      while (waitpid(-1, &i, WNOHANG) > 0);
    }
  }
  return error;
}

static bin *cmp_bin_b;

static int cmp_bin_entries(const void *s1, const void *s2) {
  const int ix = *(const int *)s1;
  const int iy = *(const int *)s2;
  const bin *b = cmp_bin_b;
  return (b->entry[ix << 1] & 63) - (b->entry[iy << 1] & 63);
}

static int cmp_prefix(const prefix *p1, const prefix *p2) {
  if(p2->count < p1->count) return -1;
  if(p1->count < p2->count) return 1;
  return 0;
}

static void process_contigs(void) {
  int idx[64];
  uint16_t start_ix[64];
  if(header != NULL) printf("%s\n", header);
	else printf("track name = dbSNP_index description = \"dbSNP index produced by dbSNP_idx\"\n");
  prefix **pref_list = NULL;
  if(prefixes != NULL) {
    HASH_SORT(prefixes, cmp_prefix);
    pref_list = malloc(sizeof(void *) * n_prefix);
    int ix = 0;
    for(prefix *p = prefixes; p != NULL; p = p->hh.next) {
      pref_list[p->ix] = p;
      p->ix = ix++;
      printf("+%s\n", p->pref);
    }
  }
  for(contig *ctg = contigs; ctg != NULL; ctg = ctg->hh.next) {
    printf(">%s\t%u\t%u\n", ctg->name, ctg->min_bin, ctg->max_bin);
    bin *b = ctg->bins;
    uint32_t curr_bin = ctg->min_bin;
    for(uint32_t i = ctg->min_bin; i <= ctg->max_bin; i++, b++) {
      int ne = b->n_entries;
      if(!ne) continue;
      if(i > curr_bin) {
	if(i - curr_bin > 1) printf("+%d\n", i - curr_bin);
	else fputs("+\n", stdout);
      }
      curr_bin = i;
      // Sort entries in bin
      for(int j = 0; j < ne; j++) idx[j] = j;
			cmp_bin_b = b;
      qsort(idx, ne, sizeof(int), cmp_bin_entries);
      // Get starting positions for names
      uint16_t x = 0;
      for(int j = 0; j < ne; j++) {
	start_ix[j] = x;
	x += b->entry[j << 1] >> 8;
      }
      for(int j = 0; j < ne; j++) {
	int j1 = idx[j];
	uint16_t z = b->entry[j1 << 1];
	int l = z >> 8;
	int ix = pref_list[(int)b->entry[(j1 << 1) + 1]]->ix;
	if(ix < 3) z |= ((ix + 1) << 6);
	printf("%02x", z & 0xff);
	if(ix >= 3) printf("%04x", ix);
	unsigned char *np = b->name_buf + start_ix[j1];
	for(int k = 0; k < l; k++) fputc(dtab2[np[k]], stdout);
	fputc('\n',stdout);
      }
    }
  }
}

int main(int argc, char *argv[]) {
  bool st = false;
  if(argc > 1) {
    for(int i = 1; i < argc && !st; i++) st = process_file(argv[i]);
  } else {
    st = process_file(NULL);
  }
  process_contigs();
  return st;
}
