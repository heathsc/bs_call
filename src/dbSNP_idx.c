//
//  dbSNP_idx.c
//  dbSNP
//
//  Created by Simon Heath on 15/11/2017.
//  Copyright 2017 Simon Heath. All rights reserved.
//

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <zlib.h>

static char *header;
static contig *contigs;
static prefix *prefixes;
static uint16_t n_prefix;
static uint64_t n_snps;

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
	n_snps++;
}

static bool process_file(char *s) {
	fprintf(stderr,"Processing file %s\n", s ? s : "<stdin>");
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

#define wrt2bufn(a,sz) { \
	if(sz + ucomp_idx > ucomp_buf_size) { \
		ucomp_buf_size *= 1.5; \
		ucomp_buf = realloc(ucomp_buf, ucomp_buf_size); \
	} \
	memcpy(ucomp_buf + ucomp_idx, a, sz); \
	ucomp_idx += sz; \
}

static void *comp_buf;
static size_t comp_buf_size, max_buf_size;

static uLongf compress_buf(const int * const buf, const size_t size) {
	uLong req_size = compressBound((uLong) size);
	if(req_size > comp_buf_size) {
		comp_buf_size = req_size * 1.2;
		comp_buf = realloc(comp_buf, comp_buf_size);
	}
	uLongf compress_size = comp_buf_size;
	if(size > max_buf_size) max_buf_size = size;
//	memcpy(comp_buf, buf, size);
//	compress_size = size;
	int ret = compress((Bytef *)comp_buf, &compress_size, (Bytef *)buf, (uLong)size);
	if(ret != 0) {
		fprintf(stderr, "Failed to compress data block\n");
		exit(-3);
	}
	return compress_size;
}

#define write_buf(f, sz) fwrite(comp_buf, 1, sz, f)
#define wrt2buf(a) wrt2bufn(&a, sizeof(a))

static void process_contigs(FILE * const fp) {
	size_t ucomp_buf_size = 4096, ucomp_idx = 0;
	void *ucomp_buf = malloc(ucomp_buf_size);
	int idx[64];
	uint16_t start_ix[64];

	const uint8_t vers = 1;
	const uint8_t zerob = 0;
	const uint8_t oneb = 1;
	wrt2buf(vers);
	wrt2buf(zerob);
	wrt2buf(n_prefix);
	uint32_t n_ctgs = HASH_COUNT(contigs);
	wrt2buf(n_ctgs);
	uint64_t * const ctg_off = malloc(sizeof(uint64_t) * n_ctgs);
	if(header != NULL) {
		wrt2bufn(header, 1 + strlen(header));
	} else {
		const char *hd = "track name = dbSNP_index description = \"dbSNP index produced by dbSNP_idx\"";
		wrt2bufn(hd, 1 + strlen(hd));
	}
	prefix **pref_list = NULL;
	if(prefixes != NULL) {
		HASH_SORT(prefixes, cmp_prefix);
		pref_list = malloc(sizeof(void *) * n_prefix);
		int ix = 0;
		for(prefix *p = prefixes; p != NULL; p = p->hh.next) {
			pref_list[p->ix] = p;
			p->ix = ix++;
			wrt2bufn(p->pref, 1 + strlen(p->pref));
		}
	}
	for(contig *ctg = contigs; ctg != NULL; ctg = ctg->hh.next) {
		wrt2buf(ctg->min_bin);
		wrt2buf(ctg->max_bin);
		wrt2bufn(ctg->name, 1 + strlen(ctg->name));
	}
	uint64_t header_size = compress_buf(ucomp_buf, ucomp_idx);
	fseek(fp, 24, SEEK_SET);
	fwrite(&header_size, 1, sizeof(header_size), fp);
	write_buf(fp, header_size);
	ucomp_idx = 0;
	int ctg_idx = 0;
	for(contig *ctg = contigs; ctg != NULL; ctg = ctg->hh.next) {
		int n_items = 0;
		ctg_off[ctg_idx++] = ftell(fp);
		fprintf(stderr, "Writing data for contig %s (%d of %d)\n", ctg->name, ctg_idx, n_ctgs);
		bin *b = ctg->bins;
		uint32_t curr_bin = ctg->min_bin;
		for(uint32_t i = ctg->min_bin; i <= ctg->max_bin; i++, b++) {
			int ne = b->n_entries;
			if(!ne) continue;
			uint32_t k = i - curr_bin;
			uint8_t x = 0;
			if(k < 64) {
				x = k << 2;
				wrt2buf(x);
			} else if(k < 256) {
				x = 1;
				wrt2buf(x);
				x = k;
				wrt2buf(x);
			} else if(k < 65536) {
				x = 2;
				uint16_t k1 = k;
				wrt2buf(x);
				wrt2buf(k1);
			} else {
				x = 3;
				wrt2buf(x);
				wrt2buf(k);
			}
			curr_bin = i;
			// Sort entries in bin
			for(int j = 0; j < ne; j++) idx[j] = j;
			cmp_bin_b = b;
			qsort(idx, ne, sizeof(int), cmp_bin_entries);
			// Get starting positions for names
			uint16_t x2 = 0;
			for(int j = 0; j < ne; j++) {
				start_ix[j] = x2;
				x2 += b->entry[j << 1] >> 8;
			}
			for(int j = 0; j < ne; j++) {
				if(j) wrt2buf(zerob);
				int j1 = idx[j];
				uint16_t z = b->entry[j1 << 1];
				int l = z >> 8;
				uint16_t ix = pref_list[(int)b->entry[(j1 << 1) + 1]]->ix;
				if(ix < 3) z |= ((ix + 1) << 6);
				uint8_t xb = z & 0xff;
				wrt2buf(xb);
				if(ix >= 3) wrt2buf(ix);
				unsigned char *np = b->name_buf + start_ix[j1];
				for(int k = 0; k < l; k++) wrt2bufn(&dtab2[np[k]], 1);
			}
			wrt2buf(oneb);
			if(++n_items == ITEMS_PER_BLOCK) {
				uint64_t sz = compress_buf(ucomp_buf, ucomp_idx);
				fwrite(&sz, 1, sizeof(sz), fp);
				write_buf(fp, sz);
				ucomp_idx = 0;
				n_items = 0;
			}
		}
		if(n_items > 0) {
			uint64_t sz = compress_buf(ucomp_buf, ucomp_idx);
			fwrite(&sz, 1, sizeof(sz), fp);
			write_buf(fp, sz);
			ucomp_idx = 0;
			n_items = 0;
		}
		// End of contig marker
		uint64_t zero = 0;
		fwrite(&zero, 1, sizeof(zero), fp);
	}
	uint64_t off = ftell(fp);
	fwrite(ctg_off, n_ctgs, sizeof(uint64_t), fp);
	fseek(fp, 0, SEEK_SET);
	const uint32_t magic = IDX_MAGIC;
	const uint32_t reserve = 0;
	const uint64_t max = max_buf_size;
	fwrite(&magic, 1, sizeof(uint32_t), fp);
	fwrite(&reserve, 1, sizeof(uint32_t), fp);
	fwrite(&off, 1, sizeof(uint64_t), fp);
	fwrite(&max, 1, sizeof(uint64_t), fp);
	fclose(fp);
}

int main(int argc, char *argv[]) {
	if(argc < 2) {
		fprintf(stderr, "usage: dbSNP_idx outfile [infile] [infile] ...\n");
		return -1;
	}
	char *const outfile = argv[1];
	FILE *fp = fopen(outfile, "wb");
	if(!fp) {
		fprintf(stderr, "Couldn't open file %s for output: %s\n", outfile, strerror(errno));
		return -2;
	}
	bool st = false;
	if(argc > 2) {
		for(int i = 2; i < argc && !st; i++) st = process_file(argv[i]);
	} else {
		st = process_file(NULL);
	}
	if(!st) process_contigs(fp);
	fprintf(stderr, "Index file %s created: %lu snps processed\n", outfile, n_snps);
	return st;
}
