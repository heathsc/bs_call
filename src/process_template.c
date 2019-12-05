/*
 * process_template.c
 *
 *  Created on: Nov 27, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <sys/wait.h>
#include <ctype.h>
#include <pthread.h>

#include "gem_tools.h"

#include "bs_call.h"

gt_status process_template_vector(gt_vector *align_list, ctg_t * const ctg, uint32_t y, sr_param *param) {
	const work_t * work = &param->work;
	bs_stats * const stats = work->stats;
	int ix = gt_vector_get_used(align_list);
	assert(ix);
	align_details **al_p = gt_vector_get_mem(align_list, align_details *);
	uint32_t x = (*al_p)->forward_position;
	if(x == 0) x = (*al_p)->reverse_position;
	assert(x > 0 && x <= y);
	x = x > 2 ? x - 2 : 1;
	uint32_t sz = y - x + 1;
	gt_string_resize(param->work.ref1, sz + 3);
	bool ret = get_sequence_string(ctg, x, sz + 2, work->vcf_ctg, work->ref1, param);
	if (ret) {
		char *ctgname = ctg->name;
		fprintf(stderr,"Problem loading reference sequence for contig '%s' %" PRIu32 " %" PRIu32 "\n", ctgname, x, sz);
		return GT_STATUS_FAIL;
	}
	for (int i = 0; i < ix; i++, al_p++) {
		if (param->left_trim || param->right_trim) {
			trim_read((*al_p)->read[0], param->left_trim, param->right_trim);
			trim_read((*al_p)->read[1], param->left_trim, param->right_trim);
		}
		uint32_t trim_left[2] = {0, 0};
		uint32_t trim_right[2] = {0, 0};
		trim_soft_clips(*al_p, stats, trim_left, trim_right);
		handle_overlap(*al_p, stats, trim_left, trim_right);
		uint32_t rdl[2]={0,0};
		for(int k = 0; k < 2; k++) {
			if((*al_p)->read[k] == NULL) continue;
			rdl[k]=gt_vector_get_used((*al_p)->read[k]);
			if(stats != NULL) {
				uint8_t *sp = gt_vector_get_mem((*al_p)->read[k], uint8_t);
				for(int k1 = 0; k1 < rdl[k]; k1++) {
					uint8_t c = *sp++;
					uint8_t q = GET_QUAL(c);
					if(q == FLT_QUAL) stats->base_filter[base_trim]++;
					else if(q < param->min_qual) stats->base_filter[base_lowqual]++;
					else stats->base_filter[base_none]++;
				}
				stats->filter_cts[gt_flt_none] ++;
				stats->filter_bases[gt_flt_none] += rdl[k];
			}
			// Here we simply normalize the reads w.r.t. indels by
			// adding N's in the place of deletions and removing insertions
			uint32_t num_misms = gt_vector_get_used((*al_p)->mismatches[k]);
			gt_vector_clear(work->orig_pos);
			uint32_t del_size = 0;
			for (int z = 0; z < num_misms; z++) {
				gt_misms *misms = gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
				if (misms->misms_type == INS) del_size += misms->size;
			}
			if (del_size > 0) gt_vector_reserve_additional((*al_p)->read[k], del_size);
			uint8_t * const sp = gt_vector_get_mem((*al_p)->read[k], uint8_t);
			gt_vector_reserve(work->orig_pos, rdl[k] + del_size, true);
			// *Orig tracks position in original read (taking into account the original strand)
			int * const orig = gt_vector_get_mem(work->orig_pos, int);
			int posx;
			if(k) {
				posx = rdl[k] + trim_right[k] - 1;
				for(int k1 = 0; k1 < rdl[k]; k1++) orig[k1] = posx - k1;
			} else {
				posx = trim_left[k];
				for(int k1 = 0; k1 < rdl[k]; k1++) orig[k1] = posx + k1;
			}
			uint32_t adj=0, ix1;
			for (int z = 0; z < num_misms; z++) {
				gt_misms *misms =  gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
				ix1 = misms->position + adj;
				if (misms->misms_type == INS) {
					memmove(sp + ix1 + misms->size, sp + ix1, rdl[k] + adj - ix1);
					memmove(&orig[ix1 + misms->size], &orig[ix1], sizeof(int) * (rdl[k] + adj - ix1));
					for (int k1 = 0; k1 < misms->size; k1++) {
						sp[ix1 + k1] = 0;
						orig[ix1 + k1] = -1;
					}
					adj += misms->size;
				} else if (misms->misms_type == DEL) {
					memmove(sp + ix1, sp + ix1 + misms->size, rdl[k] + adj - ix1 - misms->size);
					memmove(&orig[ix1], &orig[ix1 + misms->size], sizeof(int) * (rdl[k] + adj - ix1 - misms->size));
					adj -= misms->size;
				}
			}
			ix1 = rdl[k] + adj;
			gt_vector_set_used((*al_p)->read[k], ix1);
			gt_vector_set_used(work->orig_pos, ix1);
			meth_profile(*al_p, k, x, work->orig_pos, param);
		}

	}
	if (ix) call_genotypes_ML(ctg, align_list, x, y, param);
	return GT_STATUS_OK;
}

