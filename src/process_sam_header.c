/*
 * process_sam_header.c
 *
 *  Created on: Dec 2, 2019
 *      Author: heath
 */



#include <stdio.h>
#include <sys/wait.h>
#include <pthread.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "gem_tools.h"
#include "bs_call.h"

typedef struct {
	uint32_t x, y;
} pair_t;

typedef struct {
	ctg_t *ctg;
	char *name;
	gt_vector *coords;
	UT_hash_handle hh;
} ctg_hash;

static ctg_hash *new_ctg_hash(ctg_t * const ctg) {
	ctg_hash *tp = gt_alloc(ctg_hash);
	tp->coords = NULL;
	tp->ctg = ctg;
	tp->name = tp->ctg->name;
	return tp;
}

static ctg_t *new_contig(char * const name) {
	ctg_t * const ctg = gt_alloc(ctg_t);
	ctg->bam_tid = ctg->fai_id = ctg->vcf_rid = -1;
	ctg->seq = NULL;
	ctg->ctg_stats = NULL;
	ctg->seq_len = 0;
	ctg->start_pos = ctg->end_pos = 0;
	ctg->name = name;
	ctg->curr_reg = NULL;
	return ctg;
}

// Make a common list of contigs that are in the supplied bed file (if present), and in the
// sam/bam header and in the reference sequence
bool process_sam_header(sr_param * const par, bam_hdr_t * hdr) {
	bool err = false;

	ctg_hash *contigs = NULL;
	// First we store the contigs in the supplied bed file
	if(par->contig_bed != NULL) {
		FILE *cfp = fopen(par->contig_bed, "r");
		if(!cfp) fprintf(stderr,"Could not open contig bed file '%s' for input\n", par->contig_bed);
		else {
			char *buf = NULL;
			size_t buf_size = 0;
			while(true) {
				ssize_t l = getline(&buf, &buf_size, cfp);
				if(l < 0) break;
				char *p = strchr(buf, '\t');
				if(p) *p = 0;
				else continue;
				char *p1;
				uint32_t x = strtol(p + 1, &p1, 10);
				if(p == p1 || !isspace((int)(*p1))) continue;
				p = strchr(p1, '\t');
				if(!p) continue;
				uint32_t y = strtol(p + 1, &p1, 10);
				if(p == p1 || !(isspace((int)(*p1)) || *p1 == 0)) continue;
				if(x >= y) {
					fprintf(stderr, "Invalid region in contig bed file - x >= y\n");
					continue;
				}
				size_t sz = strlen(buf);
				if(sz > 0) {
					ctg_hash *tp = NULL;
					HASH_FIND_STR(contigs, buf, tp);
					if (!tp) {
						char *cname = gt_malloc(sz + 1);
						memcpy(cname, buf, sz + 1);
						tp = new_ctg_hash(new_contig(cname));
						tp->coords = gt_vector_new(1, sizeof(pair_t));
						pair_t *pair = gt_vector_get_mem(tp->coords, pair_t);
						pair->x = x + 1;
						pair->y = y;
						gt_vector_set_used(tp->coords, 1);
						HASH_ADD_KEYPTR(hh, contigs, cname, sz, tp);
					} else {
						int n = gt_vector_get_used(tp->coords);
						pair_t *pair = gt_vector_get_mem(tp->coords, pair_t);
						int k;
						for(k = 0; k < n; k++) if(y >= pair[k].x && x<= pair[k].y) break;
						if(k < n) {
							fprintf(stderr,"Warning - Region %s:%u-%u overlaps with previous region and will be ignored\n", buf, x + 1, y);
							continue;
						}
						gt_vector_reserve_additional(tp->coords, 1);
						pair = gt_vector_get_free_elm(tp->coords, pair_t);
						pair->x = x + 1;
						pair->y = y;
						gt_vector_inc_used(tp->coords);
					}
				}
			}
			fclose(cfp);
			if(buf != NULL) free(buf);
		}
	}
	// Now add contigs from reference sequence
	faidx_t *idx = par->work.seq_idx;
	int nsq = faidx_nseq(idx);
	for(int i = 0; i < nsq; i++) {
		const char * const name = faidx_iseq(idx, i);
		ctg_hash *tp = NULL;
		HASH_FIND_STR(contigs, name, tp);
		if(tp == NULL) {
			if(par->contig_bed == NULL) {
				size_t sz = strlen(name);
				char *cname = gt_malloc(sz + 1);
				memcpy(cname, name, sz + 1);
				tp = new_ctg_hash(new_contig(cname));
				HASH_ADD_KEYPTR(hh, contigs, cname, sz, tp);
			}
		}
		if(tp != NULL) {
			tp->ctg->fai_id = i;
			tp->ctg->seq_len = faidx_seq_len(idx, name);
			if(tp->coords) {
				int n = gt_vector_get_used(tp->coords);
				pair_t *pair = gt_vector_get_mem(tp->coords, pair_t);
				for(int k = 0; k < n; k++, pair++) {
					if(pair->y > tp->ctg->seq_len) {
						fprintf(stderr, "Warning - requested region (%u - %d) for %s outside of reference sequence (%u - %d) - trimming region\n", pair->x, pair->y, name, 1, tp->ctg->seq_len);
					}
				}
			}
		}
	}
	// And contigs from the SAM/BAM header
	if(hdr->n_targets) {
		par->work.tid2id = gt_malloc(hdr->n_targets * sizeof(int));
		for(int k = 0; k < hdr->n_targets; k++) par->work.tid2id[k] = -1;
		for(int i = 0; i < hdr->n_targets; i++) {
			const char * const name = hdr->target_name[i];
			ctg_hash *tp = NULL;
			HASH_FIND_STR(contigs, name, tp);
			if(tp == NULL) {
				if(par->contig_bed == NULL) {
					size_t sz = strlen(name);
					char *cname = gt_malloc(sz + 1);
					memcpy(cname, name, sz + 1);
					tp = new_ctg_hash(new_contig(cname));
					HASH_ADD_KEYPTR(hh, contigs, cname, sz, tp);
				}
			}
			if(tp != NULL) {
				ctg_t *ctg = tp->ctg;
				ctg->bam_tid = i;
				if(ctg->fai_id >= 0) {
					if(ctg->seq_len != hdr->target_len[i]) {
						fprintf(stderr, "Warning: mismatch in sequence length for contig %s between reference sequence and SAM/BAM header\n", name);
					}
				} else ctg->seq_len = hdr->target_len[i];
			}
		}
	}
	// Count retained contigs
	int n_ctg = 0;
	for(ctg_hash *tp = contigs; tp != NULL; tp = tp->hh.next) {
		if(par->contig_bed == NULL) {
			if(tp->ctg->fai_id >= 0 && tp->ctg->bam_tid >= 0) n_ctg++;
		} else if(tp->coords != NULL) {
			if(tp->ctg->fai_id < 0) fprintf(stderr, "Requested contig %s not in reference file - region omitted\n", tp->name);
			else if(tp->ctg->bam_tid < 0) fprintf(stderr, "Requested contig %s not in SAM/BAM header - region omitted\n", tp->name);
			else n_ctg++;
		}
	}
	if(n_ctg) {
		par->work.n_contigs = n_ctg;
		par->work.contigs = gt_malloc(n_ctg * sizeof(ctg_t *));
		int n_reg = 0;
		int k = 0;
		for(ctg_hash *tp = contigs; tp != NULL; tp = tp->hh.next) {
			if(tp->ctg->fai_id >= 0 && tp->ctg->bam_tid >= 0) {
				assert(k < n_ctg);
				par->work.tid2id[tp->ctg->bam_tid] = k;
				if(par->report_file != NULL) tp->ctg->ctg_stats = gt_calloc((size_t)1, gt_ctg_stats, true);
				par->work.contigs[k++] = tp->ctg;
				if(tp->coords != NULL) n_reg += gt_vector_get_used(tp->coords);
			}
		}
		assert(k == n_ctg);
		par->work.n_regions = n_reg;
		region_t *reg = NULL;
		if(n_reg > 0) reg = par->work.regions = gt_malloc(n_reg * sizeof(region_t));
		k = 0;
		for(ctg_hash *tp = contigs; tp != NULL; tp = tp->hh.next) {
			if(tp->ctg->fai_id >= 0 && tp->ctg->bam_tid >= 0) {
				if(tp->coords != NULL) {
					int n = gt_vector_get_used(tp->coords);
					pair_t *pair = gt_vector_get_mem(tp->coords, pair_t);
					for(int j = 0; j < n; j++, pair++) {
						assert(k < n_reg);
						reg[k].ctg = tp->ctg;
						reg[k].start = pair->x;
						reg[k++].stop = pair->y;
					}
				}
				tp->ctg = NULL;
			}
		}
		assert(k == n_reg);
	}
	ctg_hash *tp = contigs;
	while(tp) {
		ctg_hash *tp1 = tp->hh.next;
		if(tp->ctg != NULL) {
			free(tp->ctg->name);
			free(tp->ctg);
		}
		HASH_DEL(contigs, tp);
		tp = tp1;
	}
	return err;
}
