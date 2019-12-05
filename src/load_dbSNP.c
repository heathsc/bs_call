/*
 * load_dbSNP.c
 *
 *  Created on: Nov 24, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

static void store_dbsnp_entries(dbsnp_bin *bin, int n_entries, int name_buf_sz, uint16_t *entries, uint8_t *name_buf) {
	bin->entries = malloc(sizeof(uint16_t) * n_entries);
	bin->name_buf = malloc((size_t)name_buf_sz);
	bin->n_entries = n_entries;
	uint64_t msk = (uint64_t)0;
	for(int i = 0; i < n_entries; i++) {
		bin->entries[i] = entries[i];
		msk |= ((uint64_t)1 << (entries[i] & 63));
	}
	bin->mask = msk;
	memcpy(bin->name_buf, name_buf, name_buf_sz);
}

void *read_dbSNP_file(void *arg) {
	bool ok = true;
	sr_param *par = arg;

	static uint8_t db_tab[] = {
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x10, 0x11, 0x12, 0x13, 0x14,
			0x15, 0x16, 0x17, 0x18, 0x19, 0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, 0x30,
			0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x40, 0x41, 0x42, 0x43, 0x44, 0x45, 0x46,
			0x47, 0x48, 0x49, 0x50, 0x51, 0x52, 0x53, 0x54, 0x55, 0x56, 0x57, 0x58, 0x59, 0x60, 0x61, 0x62,
			0x63, 0x64, 0x65, 0x66, 0x67, 0x68, 0x69, 0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77, 0x78,
			0x79, 0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87, 0x88, 0x89, 0x90, 0x91, 0x92, 0x93, 0x94,
			0x95, 0x96, 0x97, 0x98, 0x99, 0x0f, 0x1f, 0x2f, 0x3f, 0x4f, 0x5f, 0x6f, 0x7f, 0x8f, 0x9f, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff,
			0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff
	};

	fprintf(stderr,"Loading dbSNP from %s\n", par->dbSNP_name);
	gt_input_file *file = gt_input_file_generic_open(par->dbSNP_name, false);
	gt_string *sbuf = gt_string_new(128);
	while(gt_expect_true(!file->eof && GT_INPUT_FILE_CURRENT_CHAR(file) != EOL && GT_INPUT_FILE_CURRENT_CHAR(file) != DOS_EOL)) {
		gt_string_append_char(sbuf,GT_INPUT_FILE_CURRENT_CHAR(file));
		GT_INPUT_FILE_NEXT_CHAR(file);
	}
	GT_INPUT_FILE_SKIP_EOL(file);
	gt_string *ctg_name = gt_string_new(128);
	dbsnp_ctg *ctg = NULL;
	uint16_t *entries = malloc(sizeof(uint16_t) * 64);
	uint8_t *name_buf = malloc(sizeof(uint8_t) * 256 * 64);
	int n_snps = 0, n_bins = 0, n_ctgs = 0;
	work_t * const work = &par->work;
	if(sbuf->length > 6 && !strncmp(sbuf->buffer, "track ", 6)) {
		work->dbSNP_header = malloc(sbuf->length - 5);
		memcpy(work->dbSNP_header, sbuf->buffer + 6, sbuf->length - 6);
		work->dbSNP_header[sbuf->length - 6] = 0;
		int n_p_store = 8;
		char **p_store = malloc(sizeof(void *) * n_p_store);
		while(gt_expect_true(!file->eof)) {
			if(GT_INPUT_FILE_CURRENT_CHAR(file) == '+') {
				gt_string_clear(sbuf);
				GT_INPUT_FILE_NEXT_CHAR(file);
				while(gt_expect_true(!file->eof && GT_INPUT_FILE_CURRENT_CHAR(file) != EOL && GT_INPUT_FILE_CURRENT_CHAR(file) != DOS_EOL)) {
					gt_string_append_char(sbuf,GT_INPUT_FILE_CURRENT_CHAR(file));
					GT_INPUT_FILE_NEXT_CHAR(file);
				}
				GT_INPUT_FILE_SKIP_EOL(file);
				if(sbuf->length > 0) {
					if(work->n_dbSNP_prefixes == 0xffff) {
						fprintf(stderr, "Error in dbSNP file: too many prefixes\n");
						exit(-1);
					}
					if(work->n_dbSNP_prefixes == n_p_store) {
						n_p_store *= 1.5;
						p_store = realloc(p_store, sizeof(void *) * n_p_store);
					}
					char *tp = p_store[work->n_dbSNP_prefixes++] = malloc(sbuf->length + 1);
					memcpy(tp, sbuf->buffer, sbuf->length);
					tp[sbuf->length] = 0;
				}
			} else break;
		}
		if(!work->n_dbSNP_prefixes) {
			fprintf(stderr, "Error in dbSNP file: no prefix information\n");
			ok = false;
		} else {
			work->dbSNP_prefix = malloc(sizeof(void *) * work->n_dbSNP_prefixes);
			memcpy(work->dbSNP_prefix, p_store, sizeof(void *) * work->n_dbSNP_prefixes);
			free(p_store);
		}
		int cbin = 0, max_bin = 0, n_entries = 0, name_buf_ptr = 0, ix = -1;
		dbsnp_bin *bins = NULL;
		while(gt_expect_true(!file->eof && ok)) {
			if(GT_INPUT_FILE_CURRENT_CHAR(file) == '>') {
				if(n_entries) {
					store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
					n_bins++;
					n_snps += n_entries;
				}
				GT_INPUT_FILE_NEXT_CHAR(file);
				gt_string_clear(sbuf);
				while(gt_expect_true(!file->eof && GT_INPUT_FILE_CURRENT_CHAR(file) != EOL && GT_INPUT_FILE_CURRENT_CHAR(file) != DOS_EOL && GT_INPUT_FILE_CURRENT_CHAR(file) != '\t')) {
					gt_string_append_char(sbuf,GT_INPUT_FILE_CURRENT_CHAR(file));
					GT_INPUT_FILE_NEXT_CHAR(file);
				}
				if(file->eof || GT_INPUT_FILE_CURRENT_CHAR(file) != '\t') {
					fprintf(stderr,"Error in dbSNP file - bad chromosome header\n");
					ok = false;
					break;
				}
				GT_INPUT_FILE_NEXT_CHAR(file);
				if(cbin != max_bin) {
					fprintf(stderr, "Error in dbSNP file - wrong number of bins (expected %d, saw %d\n", max_bin, cbin);
					ok = false;
					break;
				}
				cbin = 0, max_bin = 0;
				while(gt_expect_true(!file->eof && gt_is_number(GT_INPUT_FILE_CURRENT_CHAR(file)))) {
					cbin = cbin * 10 + gt_get_cipher(GT_INPUT_FILE_CURRENT_CHAR(file));
					GT_INPUT_FILE_NEXT_CHAR(file);
				}
				if(file->eof || GT_INPUT_FILE_CURRENT_CHAR(file) != '\t') {
					fprintf(stderr,"Error in dbSNP file - bad chromosome header\n");
					ok = false;
					break;
				}
				GT_INPUT_FILE_NEXT_CHAR(file);
				while(gt_expect_true(!file->eof && gt_is_number(GT_INPUT_FILE_CURRENT_CHAR(file)))) {
					max_bin = max_bin * 10 + gt_get_cipher(GT_INPUT_FILE_CURRENT_CHAR(file));
					GT_INPUT_FILE_NEXT_CHAR(file);
				}
				if(max_bin < cbin) {
					fprintf(stderr,"Error in dbSNP file - bad chromosome header\n");
					ok = false;
					break;
				}
				gt_string_copy(ctg_name, sbuf);
				HASH_FIND(hh, work->dbSNP, ctg_name->buffer, sbuf->length, ctg);
				if(gt_expect_false(ctg != NULL)) {
					fprintf(stderr,"Error in dbSNP file - duplicate contigs (%*.s)\n", PRIgts_content(ctg_name));
					ok = false;
					break;
				}
				ctg = malloc(sizeof(dbsnp_ctg));
				ctg->name = malloc(sbuf->length + 1);
				ctg->min_bin = cbin;
				ctg->max_bin = max_bin;
				ctg->bins = malloc(sizeof(dbsnp_bin) * (max_bin - cbin + 1));
				memcpy(ctg->name, sbuf->buffer, sbuf->length);
				ctg->name[sbuf->length] = 0;
				HASH_ADD_KEYPTR(hh, work->dbSNP, ctg->name, sbuf->length, ctg);
				n_entries = name_buf_ptr = 0;
				ix = -1;
				n_ctgs++;
				bins = ctg->bins;
			} else {
				if(ctg == NULL) {
					fprintf(stderr,"Error in dbSNP file - missing contig header\n");
					ok = false;
					break;
				} else if(GT_INPUT_FILE_CURRENT_CHAR(file) == '+') {
					if(n_entries) {
						store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
						n_bins++;
						n_snps += n_entries;
					}
					int d = 0;
					GT_INPUT_FILE_NEXT_CHAR(file);
					while(gt_expect_true(!file->eof && gt_is_number(GT_INPUT_FILE_CURRENT_CHAR(file)))) {
						d = d * 10 + gt_get_cipher(GT_INPUT_FILE_CURRENT_CHAR(file));
						GT_INPUT_FILE_NEXT_CHAR(file);
					}
					if(!d) d = 1;
					cbin += d;
					bins += d;
					if(cbin > max_bin) {
						fprintf(stderr,"Error in dbSNP file - too many bins for chromosome\n");
						ok = false;
						break;
					}
					n_entries = name_buf_ptr = 0;
					ix = -1;
				} else {
					if(n_entries == 64) {
						fprintf(stderr,"Error in dbSNP file - too many entries for bin (max 64)\n");
						ok = false;
						break;
					}
					gt_string_clear(sbuf);
					while(gt_expect_true(!file->eof && GT_INPUT_FILE_CURRENT_CHAR(file) != EOL && GT_INPUT_FILE_CURRENT_CHAR(file) != DOS_EOL)) {
						gt_string_append_char(sbuf,GT_INPUT_FILE_CURRENT_CHAR(file));
						GT_INPUT_FILE_NEXT_CHAR(file);
					}
					if(sbuf->length < 3 || sbuf->length > 256) {
						fprintf(stderr,"Error in dbSNP file - bad line length: %" PRIu64 "\n", sbuf->length);
						ok = false;
						break;
					}
					int ix1 = (gt_get_hex_cipher(sbuf->buffer[0]) << 4) | gt_get_hex_cipher(sbuf->buffer[1]);
					int prefix_ix = ix1 >> 6;
					if(prefix_ix > work->n_dbSNP_prefixes) {
						fprintf(stderr,"Error in dbSNP file - invalid prefix\n");
						ok = false;
						break;
					}
					ix1 &= 63;
					if(ix1 <= ix) {
						fprintf(stderr,"Error in dbSNP file - entries out of order or invalid\n");
						ok = false;
						break;
					}
					ix = ix1;
					int kx = 2;
					if(prefix_ix == 0) {
						if(sbuf->length < 7) {
							ok = false;
							fprintf(stderr,"Error in dbSNP file - bad line length: %" PRIu64 "\n", sbuf->length);
							break;
						}
						uint8_t ix_high = (gt_get_hex_cipher(sbuf->buffer[2]) << 4) | gt_get_hex_cipher(sbuf->buffer[3]);
						uint8_t ix_low = (gt_get_hex_cipher(sbuf->buffer[4]) << 4) | gt_get_hex_cipher(sbuf->buffer[5]);
						name_buf[name_buf_ptr++] = ix_high;
						name_buf[name_buf_ptr++] = ix_low;
						kx = 6;
					}
					int k = sbuf->length - kx;
					entries[n_entries++] = (k << 8) | (prefix_ix << 6) | ix;
					uint8_t *tp = (uint8_t *)(sbuf->buffer + kx);
					for(int j = 0; j < k; j++) name_buf[name_buf_ptr++] =  db_tab[(int)tp[j]];
				}
			}
			GT_INPUT_FILE_SKIP_EOL(file);
		}
		if(n_entries) {
			store_dbsnp_entries(bins, n_entries, name_buf_ptr, entries, name_buf);
			n_bins++;
			n_snps += n_entries;
		}
	} else {
		fprintf(stderr,"Error reading track line\n");
		ok = false;
	}

	gt_string_delete(sbuf);
	gt_string_delete(ctg_name);
	gt_input_file_close(file);
	free(entries);
	free(name_buf);
	if(ok) fprintf(stderr,"Completed loading dbSNP (no. contigs %d, no. bins %d, no. SNPs %d\n", n_ctgs, n_bins, n_snps);
	else fprintf(stderr,"Error loading dbSNP\n");
	return NULL;
}

