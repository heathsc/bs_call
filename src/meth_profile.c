/*
 * meth_profile.c
 *
 *  Created on: Nov 27, 2019
 *      Author: heath
 */

#include <stdio.h>

#include "gem_tools.h"
#include "bs_call.h"

typedef struct {
	sr_param *par;
	uint32_t x;
	gt_vector *align_list;
} mprofile_par;

//static const char base_tab_st[3][256] = {
//		{['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4},
//		{['A'] = 1, ['C'] = 6, ['G'] = 3, ['T'] = 8},
//		{['A'] = 5, ['C'] = 2, ['G'] = 7, ['T'] = 4}};

static const char base_tab_st[3][4] = {
		{1, 2, 3, 4}, {1, 6, 3, 8}, {5, 2, 7, 4}
};

static int mprof_tab[4][4] = {
		{0, 0, 0, 0}, // A in reference
		{0, 1, 0, 2}, // C in reference, informative C is non-converted, informative T is converted
		{2, 0, 1, 0}, // G in reference, informative G is non-converted, informative A is converted
		{0, 0, 0, 0} // T in reference
};

void meth_profile(const align_details * const al, const int k, const uint32_t x, gt_vector * const orig_pos, sr_param * const par) {
	if(!al->read[k]) return;
	uint32_t rl = gt_vector_get_used(al->read[k]);
	if(rl == 0) return;
	gt_vector *mprof = par->work.stats->meth_profile;
	int min_qual = par->min_qual;
	const char *ref_st = gt_string_get_string(par->work.ref1);
	const char * const base_tab = par->defs.base_tab;
	int mx = gt_vector_get_used(mprof);
	int st = al->bs_strand;
	uint8_t *sp = gt_vector_get_mem(al->read[k], uint8_t);
	uint32_t pos = k ? al->reverse_position : al->forward_position;
	const char *rf = ref_st + pos - x;
	uint32_t j;
	int *posx = gt_vector_get_mem(orig_pos, int);
	char prev_r = pos > x ? base_tab[(int)ref_st[pos - x - 1]] : 0;
	for(j = 0; j < rl; j++) {
		uint8_t q = GET_QUAL(sp[j]);
		if (q >= min_qual && q != FLT_QUAL) {
			int c = base_tab_st[st][GET_BASE(sp[j])];
			int r = base_tab[(int)rf[j]];
			if(r) {
				int xx = (c > 4) ? mprof_tab[r - 1][c - 5] : mprof_tab[r - 1][c - 1] + 2;
				if(xx) {
					// Exclude CpG positions
					if((r == 3 && prev_r != 2) || (r == 2 && rf[j + 1] != 'G')) {
						int p = posx[j];
						assert(p >= 0);
						if(p >= mx) {
							gt_vector_reserve(mprof, p + 1, true);
							mx = p + 1;
						}
						meth_cts *mc = gt_vector_get_elm(mprof, p, meth_cts);
						mc->conv_cts[xx - 1]++;
					}
				}
			}
			prev_r = r;
		}
	}
	if(mx > gt_vector_get_used(mprof)) gt_vector_set_used(mprof, mx);
}
