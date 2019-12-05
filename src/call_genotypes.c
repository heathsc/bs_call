/*
 * call_genotypes.c
 *
 *  Created on: Nov 27, 2019
 *      Author: heath
 */

#include <stdio.h>
#include <sys/wait.h>
#include <ctype.h>
#include <pthread.h>
#include <gsl/gsl_multimin.h>
#include "gem_tools.h"

#include "bs_call.h"

static const char base_tab_st[3][4] = {
		{1, 2, 3, 4}, {1, 6, 3, 8}, {5, 2, 7, 4}
};

typedef struct {
	pthread_t thr;
	sr_param *param;
	gt_meth *gtm;
	pileup *cts;
	uint32_t x;
	uint32_t y;
	int step;
} cthread_par;

void *call_thread(void *arg) {
	cthread_par *cpar = arg;
	uint32_t x = cpar->x, y = cpar->y;
	int step = cpar->step;
	gt_meth *tg = cpar->gtm;
	pileup *tp = cpar->cts;
	work_t * const work = &cpar->param->work;
	gt_vcf *v = work->vcf + x - work->vcf_x;
	const sr_param * const par = cpar->param;
	const char *ref_st = gt_string_get_string(par->work.ref);
	for (uint32_t i = x; i <= y; i += step, tp += step, tg += step, v += step) {
		if (tp->n) {
			float tot_qual = 0.0;
			for(int j = 0; j < 8; j++) {
				float nn = (float)(tp->counts[0][j] + tp->counts[1][j]);
				if(nn > 0) {
					tot_qual += tp->quality[j];
					tg->qual[j] = (int)floorf(0.5 + tp->quality[j] / nn);
				} else tg->qual[j] = 0;
			}
			tg->aq = (int)floorf(0.5 + tot_qual / (float)tp->n);
			tg->mq = (int)(0.5 + sqrt(tp->mapq2 / (float)tp->n));
			for (int j = 0; j < 8; j++) {
				if (tp->counts[0][j] + tp->counts[1][j]) {
					tg->counts[j] = tp->counts[0][j] + tp->counts[1][j];
				}
			}
			calc_gt_prob(tg, par, ref_st[i - work->vcf_x]);
			double fs = 0.0;
			if(par->defs.gt_het[tg->max_gt]) {
				int ftab[4];
				switch(tg->max_gt) {
				case 1: // AC
					ftab[0] = tp->counts[0][0] + tp->counts[0][4];
					ftab[1] = tp->counts[0][1] + tp->counts[0][5] + tp->counts[0][7];
					ftab[2] = tp->counts[1][0] + tp->counts[1][4];
					ftab[1] = tp->counts[1][1] + tp->counts[1][5] + tp->counts[1][7];
					break;
				case 2: // AG
					ftab[0] = tp->counts[0][0];
					ftab[1] = tp->counts[0][2] + tp->counts[0][6];
					ftab[2] = tp->counts[1][0];
					ftab[1] = tp->counts[1][2] + tp->counts[1][6];
					break;
				case 3: // AT
					ftab[0] = tp->counts[0][0] + tp->counts[0][4];
					ftab[1] = tp->counts[0][3] + tp->counts[0][7];
					ftab[2] = tp->counts[1][0] + tp->counts[1][4];
					ftab[1] = tp->counts[1][3] + tp->counts[1][7];
					break;
				case 5: // CG
					ftab[0] = tp->counts[0][1] + tp->counts[0][5] + tp->counts[0][7];
					ftab[1] = tp->counts[0][2] + tp->counts[0][4] + tp->counts[0][6];
					ftab[2] = tp->counts[1][1] + tp->counts[1][5] + tp->counts[1][7];
					ftab[3] = tp->counts[1][2] + tp->counts[1][4] + tp->counts[1][6];
					break;
				case 6: // CT
					ftab[0] = tp->counts[0][1] + tp->counts[0][5];
					ftab[1] = tp->counts[0][3];
					ftab[2] = tp->counts[1][1] + tp->counts[1][5];
					ftab[3] = tp->counts[1][3];
					break;
				case 8: // GT
					ftab[0] = tp->counts[0][2] + tp->counts[0][4] + tp->counts[0][6];
					ftab[1] = tp->counts[0][3] + tp->counts[0][7];
					ftab[2] = tp->counts[1][2] + tp->counts[1][4] + tp->counts[0][6];
					ftab[3] = tp->counts[1][3] + tp->counts[1][7];
					break;
				default: // Should not happen
					fprintf(stderr, "Internal error: illegal option in call_thread()\n");
					break;
				}
				double z = fisher(ftab, par->defs.lfact_store);
				if(z < 1.0e-20) z = 1.0e-20;
				fs = log(z) / LOG10;
			}
			tg->fisher_strand = fs;
			memcpy(&v->gtm, tg, sizeof(gt_meth));
			v->skip = false;
		} else v->skip = true;
		v->ready = true;
		pthread_cond_signal(&work->vcf_cond);
	}
	return NULL;
}

static gt_vector *pileupv, *gt_resv;

void call_genotypes_ML(ctg_t * const ctg, gt_vector * const align_list, const uint32_t x, const uint32_t y, sr_param * const param) {
	uint32_t nr = gt_vector_get_used(align_list);
	assert(y >= x);
	uint32_t sz = y - x + 1;
	if(pileupv == NULL) pileupv = gt_vector_new(sz, sizeof(pileup));
	else gt_vector_reserve(pileupv, sz, false);
	if(gt_resv == NULL) gt_resv = gt_vector_new(sz, sizeof(gt_meth));
	else gt_vector_reserve(gt_resv, sz, false);
	pileup *counts = gt_vector_get_mem(pileupv, pileup);
	gt_meth *gt_res = gt_vector_get_mem(gt_resv, gt_meth);
	memset(counts, 0, sizeof(pileup) * sz);
	memset(gt_res, 0, sizeof(gt_meth) * sz);
	align_details **al_p = gt_vector_get_mem(align_list, align_details *);
	for (uint32_t ix = 0; ix < nr; ix++, al_p++) {
		align_details *al = *al_p;
		uint32_t x1 = al->forward_position;
		if(x1 == 0) x1 = al->reverse_position;
		else if (al->reverse_position > 0 && al->reverse_position < x1) x1 = al->reverse_position;
		assert(x1 >= x);
		int ori = al->orientation;
		assert(ori < 2);
		int st = al->bs_strand;
		for (int k = 0; k < 2; k++) {
			if(al->read[k] == NULL) continue;
			uint32_t rl = gt_vector_get_used(al->read[k]);
			if(rl == 0) continue;
			float mapq2 = al->mapq[k] * al->mapq[k];
			uint8_t *sp = gt_vector_get_mem(al->read[k], uint8_t);
			uint32_t pos = k ? al->reverse_position : al->forward_position;
			uint32_t j;
			for(j = 0; j < rl; j++) {
				uint8_t q = GET_QUAL(sp[j]);
				if(q > 0 && q != FLT_QUAL) break;
			}
			uint32_t read_start;
			if(j < rl) read_start = j;
			else continue;
			for(j = rl; j > 0; j--) {
				uint8_t q = GET_QUAL(sp[j - 1]);
				if(q > 0 && q != FLT_QUAL) break;
			}
			uint32_t read_end;
			if(j > 0) read_end = j - 1;
			else continue;
			pos += read_start;
			pileup *curr_loc = counts + pos - x;
			for (j = read_start; j <= read_end && pos <= y; j++, pos++, curr_loc++) {
				int c = base_tab_st[st][(int)GET_BASE(sp[j])] - 1;
				uint8_t q = GET_QUAL(sp[j]);
				if (q >= param->min_qual && q != FLT_QUAL) {
					curr_loc->n++;
					curr_loc->quality[c] += (float)q;
					curr_loc->mapq2 += mapq2;
					curr_loc->counts[ori][c]++;
				}
			}
			ori ^= 1;
		}
	}
	work_t * const work = &param->work;
	// Prepare printing
//	fprintf(stderr, "call_genotypes() preparing for calc\n");
	pthread_mutex_lock(&work->print_mutex);
	while(param->work.vcf_n) {
//		fprintf(stderr, "call_genotypes() waiting to dump new batch\n");
		pthread_cond_wait(&work->print_cond, &work->print_mutex);
	}
	pthread_mutex_unlock(&work->print_mutex);

	if(sz > work->vcf_size) {
		work->vcf = realloc(work->vcf, sizeof(gt_vcf) * sz);
		work->vcf_size = sz;
	}
	for(int i = 0; i < sz; i++) work->vcf[i].ready = false;
	work->vcf_x = x;
	work->vcf_ctg = ctg;
	gt_string *tp = work->ref;
	work->ref = work->ref1;
	work->ref1 = tp;
	work->vcf_n = sz;
	pthread_mutex_lock(&work->print_mutex);
//	fprintf(stderr,"call_genotypes(), deblocking print_cond\n");
	pthread_cond_signal(&param->work.print_cond);
	pthread_mutex_unlock(&work->print_mutex);
	// Set up calculation threads
	int nthr = 1 + param->num_threads[CALC_THREADS];
	cthread_par *cpar = malloc(sizeof(cthread_par) * nthr);
	for(int i = 0; i < nthr; i++) {
		cpar[i].param = param;
		cpar[i].gtm = gt_res + i;
		cpar[i].cts = counts + i;
		cpar[i].x = x + i;
		cpar[i].y = y;
		cpar[i].step = nthr;
		pthread_create(&cpar[i].thr, NULL, call_thread, cpar + i);
	}
	for(int i = 0; i < nthr; i++) pthread_join(cpar[i].thr, 0);
	pthread_mutex_lock(&work->vcf_mutex);
	pthread_cond_signal(&work->vcf_cond);
	pthread_mutex_unlock(&work->vcf_mutex);
	free(cpar);
}
