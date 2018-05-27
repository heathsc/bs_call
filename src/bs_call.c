/*
 * bs_call.c
 *
 *  Created on: 30 Sep 2014
 *      Author: heath
 */

#define BS_CALL "bs_call"
#define BS_CALL_VERSION "2.01"

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <pthread.h>
#include <gsl/gsl_multimin.h>
#include "gem_tools.h"
#include "gt_pipe_io.h"

#include "bs_call.h"
#include "bs_call_options.h"

#define LFACT_STORE_SIZE 256

static double lfact_store[LFACT_STORE_SIZE];

static double logp[100]; // logp[i] = log((i+1)/100)

sr_param param = {
  .input_file = NULL,
  .name_reference_file = NULL,
  //    .name_gem_index_file = NULL,
  .output_prefix = NULL,
  .sample_name = NULL,
  .species_filter = NULL,
  .dbSNP_name = NULL,
  .report_file = NULL,
  .species_hash = NULL,
  .mmap_input = false,
  .compress = NONE,
  .verbose = false,
  .haploid = false,
  .blank_trim = false,
  //	  .pileup = false,
  .caller = maximum_likelihood,
  .left_trim = 0,
  .right_trim = 0,
  //	  .pileup_file_name = NULL,
  .output_file = NULL,
  .json_file = NULL,
  //	  .pileup_file = NULL,
  .sam_headers = NULL,
  .mapq_thresh = DEFAULT_MAPQ_THRESH,
  .min_qual = MIN_QUAL,
  .max_template_len = DEFAULT_MAX_TEMPLATE_LEN,
  .realign_tol = DEFAULT_REALIGN_TOL,
  .under_conv = DEFAULT_UNDER_CONVERSION,
  .over_conv = DEFAULT_OVER_CONVERSION,
  .ref_bias = DEFAULT_REF_BIAS,
  .no_split = false,
  .extra_stats = false,
  .keep_duplicates = false,

  .all_positions = false,
  .num_threads = 1,
  .sequence_archive = NULL,
  .dbSNP = NULL,
  .dbSNP_header = NULL,
  .dbSNP_prefix = NULL,
  .n_dbSNP_prefixes = 0,
  .stats = NULL,
  .ctg_stats = NULL
};

void lfact_store_init(void) {
  lfact_store[0] = lfact_store[1] = 0.0;
  double l = 0.0;
  for (int i = 2; i < LFACT_STORE_SIZE; i++) {
    l += log((double)i);
    lfact_store[i] = l;
  }
}

static inline double lfact(int x) {
  return (x < LFACT_STORE_SIZE ? lfact_store[x] : lgamma((double)(x + 1)));
}

double fisher(int *c) {
  int row[2], col[2];
  row[0] = c[0] + c[1];
  row[1] = c[2] + c[3];
  col[0] = c[0] + c[2];
  col[1] = c[1] + c[3];
  int n = row[0] + row[1];
  if(n == 0) return 1.0;
  double delta = (double)c[0] - (double)(row[0] * col[0]) / (double)n;
  double knst =
    lfact(col[0]) + lfact(col[1]) + lfact(row[0]) + lfact(row[1]) - lfact(n);
  double l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
  double p = l;
  if (delta > 0.0) {
    // Decrease counter diagonal elements until zero (this will increase delta)
    int mn = c[1] < c[2] ? c[1] : c[2];
    for (int i = 0; i < mn; i++) {
      l *= (double)((c[1] - i) * (c[2] - i)) /
	(double)((c[0] + i + 1) * (c[3] + i + 1));
      p += l;
    }
    mn = c[0] < c[3] ? c[0] : c[3];
    // Calculate amount required to increase delta by decreasing leading
    // diagonal elements
    int k = ceil(2.0 * delta);
    if (k <= mn) {
      c[0] -= k;
      c[3] -= k;
      c[1] += k;
      c[2] += k;
      l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
      p += l;
      for (int i = 0; i < mn - k; i++) {
        l *= (double)((c[0] - i) * (c[3] - i)) /
	  (double)((c[1] + i + 1) * (c[2] + i + 1));
        p += l;
      }
    }
  } else {
    // Decrease leading diagonal elements until zero (this will increase delta)
    int mn = c[0] < c[3] ? c[0] : c[3];
    for (int i = 0; i < mn; i++) {
      l *= (double)((c[0] - i) * (c[3] - i)) /
	(double)((c[1] + i + 1) * (c[2] + i + 1));
      p += l;
    }
    mn = c[1] < c[2] ? c[1] : c[2];
    // Calculate amount required to increase delta by decreasing counter
    // diagonal elements
    int k = ceil(-2.0 * delta);
    if (!k)
      k = 1;
    if (k <= mn) {
      c[0] += k;
      c[3] += k;
      c[1] -= k;
      c[2] -= k;
      l = exp(knst - lfact(c[0]) - lfact(c[1]) - lfact(c[2]) - lfact(c[3]));
      p += l;
      for (int i = 0; i < mn - k; i++) {
        l *= (double)((c[1] - i) * (c[2] - i)) /
	  (double)((c[0] + i + 1) * (c[3] + i + 1));
        p += l;
      }
    }
  }
  return p;
}

void usage(const gt_option *const options, char *groups[],
           const bool print_inactive) {
  fprintf(stderr, "USE: ./bs_call [ARGS]...\n");
  gt_options_fprint_menu(stderr, options, groups, true, print_inactive);
}

static int cmp_al(const void *s1, const void *s2) {
  const align_details *al1 = *(align_details **)s1,
    *al2 = *(align_details **)s2;
  int a = 0;
  int x1, x2, y1, y2;
  if(al1->forward_position > 0) {
    x1 = al1->forward_position;
    y1 = al1->reverse_position > 0 ? al1->reverse_position : x1;
  } else x1 = y1 = al1->reverse_position;
  if(al2->forward_position > 0) {
    x2 = al2->forward_position;
    y2 = al2->reverse_position > 0 ? al2->reverse_position : x1;
  } else x2 = y2 = al2->reverse_position;
	
  if (x1 < x2) {
    a = -1;
  } else if (x1 > x2) {
    a = 1;
  } else if (y1 < y2) {
    a = -1;
  } else if (y1 > y2) {
    a = 1;
  } else if (al1->bs_strand < al2->bs_strand) {
    a = -1;
  } else if (al1->bs_strand > al2->bs_strand) {
    a = 1;
  }
  return a;
}

static uint64_t get_al_qual(align_details *al) {
  uint64_t qual = 0;
  for (int k = 0; k < 2; k++) {
    if(al->qualities[k]) {
      uint64_t rl = gt_string_get_length(al->qualities[k]);
      char *sq = gt_string_get_string(al->qualities[k]);
      for (int j = 0; j < rl; j++)
	qual += sq[j];
    }
  }
  return qual;
}

// Every entry will be set to zero be default
static const char base_tab[256] = {['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4};

static const char base_tab_st[3][256] = {
  {['A'] = 1, ['C'] = 2, ['G'] = 3, ['T'] = 4},
  {['A'] = 1, ['C'] = 6, ['G'] = 3, ['T'] = 8},
  {['A'] = 5, ['C'] = 2, ['G'] = 7, ['T'] = 4}};

static void add_flt_counts(gt_vector *v, int ct, bool var) {
  if(ct >= v->elements_allocated) gt_vector_reserve(v, ct + 1, true);
  fstats_cts *c = gt_vector_get_elm(v, ct, fstats_cts);
  c->cts[var ? 1 : 0]++;
  if(ct >= v->used) v->used = ct + 1;
}

static dbsnp_ctg *dbSNP_ctg;
static gt_ctg_stats *ctg_stats;

static bool gt_het[10] = { false, true, true, true, false, true, true, false, true, false };
static char *flt_name[] = { "q20", "qd2", "fs60", "mq40", "gof20" };

void _print_vcf_entry(FILE *fp, char *ctg, gt_meth *gtm, const char *rf_ctxt,
                      const uint64_t x, char *gt_store) {
  static const char *ref_alt[10][5] = {
    {"A", ".", "A", "A", "A"},       // AA
    {"A,C", "C", "A", "A,C", "A,C"}, // AC
    {"A,G", "G", "A,G", "A", "A,G"}, // AG
    {"A,T", "T", "A,T", "A,T", "A"}, // AT
    {"C", "C", ".", "C", "C"},       // CC
    {"C,G", "C,G", "G", "C", "C,G"}, // CG
    {"C,T", "C,T", "T", "C,T", "C"}, // CT
    {"G", "G", "G", ".", "G"},       // GG
    {"G,T", "G,T", "G,T", "T", "G"}, // GT
    {"T", "T", "T", "T", "."}        // TT
  };
  static const stats_mut mut_type[10][5] = {
    {mut_no, mut_no, mut_CA, mut_GA, mut_TA},   // AA
    {mut_no, mut_AC, mut_CA, mut_no, mut_no},       // AC
    {mut_no, mut_AG, mut_no, mut_GA, mut_no},       // AG
    {mut_no, mut_AT, mut_no, mut_no, mut_TA},       // AT
    {mut_no, mut_AC, mut_no, mut_GC, mut_TC},   // CC
    {mut_no, mut_no, mut_CG, mut_GC, mut_no},       // CG
    {mut_no, mut_no, mut_CT, mut_no, mut_TC},       // CT
    {mut_no, mut_AG, mut_CG, mut_no, mut_TG},   // GG
    {mut_no, mut_no, mut_no, mut_GT, mut_TG},       // GT
    {mut_no, mut_AT, mut_CT, mut_GT, mut_no},   // TT
  };
  static char *cs_str[10] = {"NA", "+", "-", "NA", "+",
                             "+-", "+", "-", "-",  "NA"};
  static const int all_idx[10][5][2] = {
    {{1, 0}, {0, 0}, {1, 0}, {1, 0}, {1, 0}}, // AA
    {{1, 2}, {2, 0}, {1, 0}, {1, 2}, {1, 2}}, // AC
    {{1, 3}, {3, 0}, {1, 3}, {1, 0}, {1, 3}}, // AG
    {{1, 4}, {4, 0}, {1, 4}, {1, 4}, {1, 0}}, // AT
    {{2, 0}, {2, 0}, {0, 0}, {2, 0}, {2, 0}}, // CC
    {{2, 3}, {2, 3}, {3, 0}, {2, 0}, {2, 3}}, // CG
    {{2, 4}, {2, 4}, {4, 0}, {2, 4}, {2, 0}}, // CT
    {{3, 0}, {3, 0}, {3, 0}, {0, 0}, {3, 0}}, // GG
    {{3, 4}, {3, 4}, {3, 4}, {4, 0}, {3, 0}}, // GT
    {{4, 0}, {4, 0}, {4, 0}, {4, 0}, {0, 0}}  // TT
  };
  static const char *gt_str[10][5] = {
    {"1/1", "0/0", "1/1", "1/1", "1/1"}, // AA
    {"1/2", "0/1", "0/1", "1/2", "1/2"}, // AC
    {"1/2", "0/1", "1/2", "0/1", "1/2"}, // AG
    {"1/2", "0/1", "1/2", "1/2", "0/1"}, // AT
    {"1/1", "1/1", "0/0", "1/1", "1/1"}, // CC
    {"1/2", "1/2", "0/1", "0/1", "1/2"}, // CG
    {"1/2", "1/2", "0/1", "1/2", "0/1"}, // CT
    {"1/1", "1/1", "1/1", "0/0", "1/1"}, // GG
    {"1/2", "1/2", "1/2", "0/1", "0/1"}, // GT
    {"1/1", "1/1", "1/1", "1/1", "0/0"}  // TT
  };
  static const char gt_flag[10][5] = {
    {0, 1, 0, 0, 0}, // AA
    {0, 0, 0, 0, 0}, // AC
    {0, 0, 0, 0, 0}, // AG
    {0, 0, 0, 0, 0}, // AT
    {0, 0, 0, 0, 0}, // CC
    {0, 0, 0, 0, 0}, // CG
    {0, 0, 0, 0, 0}, // CT
    {0, 0, 0, 0, 0}, // GG
    {0, 0, 0, 0, 0}, // GT
    {0, 0, 0, 0, 1}, // TT
  };

  static char *iupac = "NAMRWCSYGKT";
  static int cflag[] = {0, 1, 0, 0, 1, 1, 1, 0, 0, 0};
  static int gflag[] = {0, 0, 1, 0, 0, 1, 0, 1, 1, 0};
  static char dtab[16] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 0, 0, 0, 0, 0, 0 };
  static char *old_ctg;
  static uint64_t old_x;
  char rs[256];
  char *db_prefix;
  static uint64_t prev_cpg_x;
  static bool prev_cpg_flt;
	
  if(x == 0) return;
  if(old_ctg != NULL && !strcmp(old_ctg, ctg)) {
    if(x <= old_x) return;
    if(old_ctg != NULL) free(old_ctg);
    old_ctg = strdup(ctg);
  }
  old_x = x;
  uint64_t *counts = gtm->counts;
  uint64_t dp = 0, d_inf = 0, dp1 = 0;
  for (int i = 0; i < 4; i++) dp1 += counts[i];
  for (int i = 4; i < 8; i++) d_inf += counts[i];
  dp = dp1 + d_inf;
  if (!dp) return;
  rs[0] = 0;
  db_prefix = NULL;
  if(dbSNP_ctg != NULL) {
    int bn = x >> 6;
    if(bn >= dbSNP_ctg->min_bin && bn <= dbSNP_ctg->max_bin) {
      dbsnp_bin *b = dbSNP_ctg->bins + bn - dbSNP_ctg->min_bin;
      int ix = x & 63;
      uint64_t mk = (uint64_t)1 << ix;
      if(b->mask & mk) {
	uint64_t mk1 = b->mask & (mk - (uint64_t)1);
	int i = 0, j = 0;
	while(mk1) {
	  if(mk1 & (uint64_t)1) {
	    uint16_t en = b->entries[i++];
	    j += en >> 8;
	    if(!((en >> 6) & 3)) j += 2;
	  }
	  mk1 >>= 1;
	}
	char *tp = rs;
	int prefix_id = (b->entries[i] >> 6) & 3;
	unsigned char *tp1 = b->name_buf + j;
	if((prefix_id--) == 0) {
	  prefix_id = (tp1[0] << 8) | tp1[1];
	  tp1+=2;
	}
        db_prefix = param.dbSNP_prefix[prefix_id];
	j = b->entries[i] >> 8;
	for(int k = 0; k < j; k++) {
	  unsigned char z = *tp1++;
	  *tp++ = dtab[z >> 4];
	  *tp++ = dtab[z & 15];
	}
	*tp = 0;
      }
    }
  }
  int rfc = toupper((int)rf_ctxt[2]);
  int rfix = base_tab[rfc];
  if(rfix < 1 || rfix > 4) {
    rfc = 'N';
    rfix = 0;
  }
  int gt = gt_store[2] - 1;
  // Skip homozygous reference if AA or TT
  bool skip = (!param.all_positions && db_prefix == NULL && gt_flag[gt][rfix]);
  double z = gtm->gt_prob[gt];
  int phred;
  double z1 = exp(z * LOG10);
  if (z1 >= 1.0)
    phred = 255;
  else {
    phred = (int)(-10.0 * log(1.0 - z1) / LOG10);
    if (phred > 255)
      phred = 255;
  }
  const char *alt = ref_alt[gt][rfix];
  const stats_mut mut = mut_type[gt][rfix];
  const int fs = (int)(-gtm->fisher_strand * 10.0 + 0.5);
  const int gof = (int)(gtm->gt_gof * 10.0 + 0.5);
  const uint64_t qd = dp1 > 0 ? phred / dp1 : phred;
  uint32_t flt = 0;
  if(!skip) {
    fprintf(fp, "%s\t%" PRIu64 "\t%s%s\t%c\t%s\t%d", ctg, x, db_prefix == NULL ? "." : db_prefix, rs, rfc, alt, phred);
    // FILTER fields
    if (phred < 20) flt |= 1;
    if (qd < 2) flt |= 2;
    if (fs > 60) flt |= 4;
    if (gtm->mq < 40) flt |= 8;
    if (gof > 20) flt |= 16;
    if(!flt) {
      bool mac1 = false;
      switch(gt) {
      case 1: // AC
	mac1 = (counts[1] + counts[5] + counts[7] <= 1 || counts[0] + counts[4] <= 1);
	break;
      case 2: // AG
	mac1 = (counts[2] + counts[6] <= 1 || counts[0] <= 1);
	break;
      case 3: // AT
	mac1 = (counts[3] + counts[7] <= 1 || counts[0] + counts[4] <= 1);
	break;
      case 5: // CG
	mac1 = (counts[2] + counts[6] + counts[4] <= 1 || counts[1] + counts[5] + counts[7] <= 1);
	break;
      case 6: // CT
	mac1 = (counts[3] <= 1 || counts[1] + counts[5] <= 1);
	break;
      case 8: // GT
	mac1 = (counts[3] + counts[7] <= 1 || counts[2] + counts[6] + counts[4] <= 1);
	break;
      }
      if(mac1) {
	flt |= 128;
	fputs("\tmac1", fp);
      } else fputs("\tPASS", fp);
    } else fputs("\tfail", fp);
    // INFO field
    fprintf(fp, "\tCX=%.5s", rf_ctxt);
    // FORMAT field
    fputs("\tGT:FT:DP:MQ:GQ:QD:GL:MC8:AMQ:CS:CG:CX:GOF", fp);
    if(gt_het[gt]) fputs(":FS", fp);
  }
  // Genotype fields
  char ctxt[5];
  for (int i = 0; i < 5; i++)
    ctxt[i] = iupac[(int)gt_store[i]];
  char *cpg = ".";
  if ((gt_store[2] == 5 && gt_store[3] == 8) ||
      (gt_store[2] == 8 && gt_store[1] == 5))
    cpg = "CG";
  else if (gt_store[2] == 5) {
    if (gt_store[3]) {
      if (gflag[(int)gt_store[3] - 1])
	cpg = "H";
      else
	cpg = "N";
    } else
      cpg = "?";
  } else if (gt_store[2] == 8) {
    if (gt_store[1]) {
      if (cflag[(int)gt_store[1] - 1])
	cpg = "H";
      else
        cpg = "N";
    } else
      cpg = "?";
  } else if (cflag[(int)gt_store[2] - 1]) {
    if (gt_store[3]) {
      if (gflag[(int)gt_store[3] - 1])
        cpg = "H";
      else
        cpg = "N";
    } else
      cpg = "?";
  } else if (gflag[(int)gt_store[2] - 1]) {
    if (gt_store[1]) {
      if (cflag[(int)gt_store[1] - 1])
        cpg = "H";
      else
        cpg = "N";
    } else
      cpg = ".";
  }
  if(!skip) {
    fprintf(fp,"\t%s", gt_str[gt][rfix]);
    if(flt & 31) {
      char tmp = ':';
      int f_ix = 0;
      uint32_t flt1 = flt & 31;
      while(flt1) {
	if(flt1 & 1) {
	  fprintf(fp, "%c%s", tmp, flt_name[f_ix]);
	  tmp = ';';
	}
	flt1 >>= 1;
	f_ix++;
      }
    } else fputs(":PASS", fp);
    fprintf(fp, ":%" PRIu64 ":%d:%d:%" PRIu64, dp1, gtm->mq, phred, qd);
    const int *aix = all_idx[gt][rfix];
    if (rfix) {
      int j = rfix * (9 - rfix) / 2 + rfix - 5;
      z = gtm->gt_prob[j];
    } else
      z = -99.999;
    fprintf(fp, ":%.3f", z < -99.999 ? -99.999 : z);
    for (int i = 0; i < 2 && aix[i] > 0; i++) {
      int j;
      if (rfix) {
	if (rfix < aix[i])
	  j = rfix * (9 - rfix) / 2 + aix[i] - 5;
	else
	  j = aix[i] * (9 - aix[i]) / 2 + rfix - 5;
	z = gtm->gt_prob[j];
      }
      fprintf(fp, ",%.3f", z < -99.999 ? -99.999 : z);
      for (int k = 0; k < i; k++) {
	if (aix[k] < aix[i])
	  j = aix[k] * (9 - aix[k]) / 2 + aix[i] - 5;
	else
	  j = aix[i] * (9 - aix[i]) / 2 + aix[k] - 5;
	z = gtm->gt_prob[j];
	fprintf(fp, ",%.3f", z < -99.999 ? -99.999 : z);
      }
      j = aix[i] * (9 - aix[i]) / 2 + aix[i] - 5;
      z = gtm->gt_prob[j];
      fprintf(fp, ",%.3f", z < -99.999 ? -99.999 : z);
    }
    fprintf(fp, ":%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64 ",%" PRIu64,
	    counts[0], counts[1], counts[2], counts[3], counts[4], counts[5], counts[6], counts[7]);
    
    int k = 0;
    for(int j = 0; j < 8; j++) {
      if(counts[j] > 0) fprintf(fp, "%c%d", k++ ? ',' : ':', gtm->qual[j]);
    }
    fprintf(fp, ":%s:%s:%.5s:%d", cs_str[gt], cpg, ctxt, gof); 
    if(gt_het[gt]) fprintf(fp, ":%d", fs);
    fputc('\n', fp);
  }
  bs_stats *stats = param.stats;
  if(stats != NULL) {
    bool snp = false;
    bool multi = false;
    gt_cov_stats *gcov;
    HASH_FIND(hh, stats->cov_stats, &dp, sizeof(uint64_t), gcov);
    if(gcov == NULL) {
      gcov = calloc((size_t)1, sizeof(gt_cov_stats));
      gcov->coverage = dp;
      HASH_ADD(hh, stats->cov_stats, coverage, sizeof(uint64_t), gcov);
    }
    gcov->all++;
    int bn = x / 100;
    if(bn < ctg_stats->nbins) {
      int gc = ctg_stats->gc[bn];
      if(gc <= 100) gcov->gc_pcent[gc]++;
    }
    if(!skip) {
      if(alt[0] != '.') { //Variant site
	if(alt[1] == ',') multi = true;
	else snp = true;
	if(snp) {
	  stats->snps[stats_all]++;
	  ctg_stats->snps[stats_all]++;
	  if(!flt) {
	    stats->snps[stats_passed]++;
	    ctg_stats->snps[stats_passed]++;
	  }
	} else {
	  stats->multi[stats_all]++;
	  ctg_stats->multi[stats_all]++;
	  if(!flt) {
	    stats->multi[stats_passed]++;
	    ctg_stats->multi[stats_passed]++;
	  }
	}
	stats->qual[variant_sites][phred]++;
	gcov->var++;
      }
      add_flt_counts(stats->qd_stats, qd, gt_het[gt]);
      add_flt_counts(stats->fs_stats, fs, gt_het[gt]);
      add_flt_counts(stats->mq_stats, gtm->mq, gt_het[gt]);
      add_flt_counts(stats->gof_stats, gof, gt_het[gt]);
      stats->filter_counts[gt_het[gt] ? 1 : 0][flt & 31]++;
      stats->qual[all_sites][phred]++;			
      if(db_prefix != NULL) {
	stats->dbSNP_sites[stats_all]++;
	ctg_stats->dbSNP_sites[stats_all]++;
	if(snp || multi) {
	  stats->dbSNP_var[stats_all]++;
	  ctg_stats->dbSNP_var[stats_all]++;
	}
	if(!flt) {
	  stats->dbSNP_sites[stats_passed]++;
	  ctg_stats->dbSNP_sites[stats_passed]++;
	  if(snp || multi) {
	    stats->dbSNP_var[stats_passed]++;
	    ctg_stats->dbSNP_var[stats_passed]++;
	  }
	}
      }
      if(!strcmp(cpg, "CG")) {
	bool ref_cpg = false;
	bool cpg_ok = false;
	uint64_t a,b;
	if(!strcmp(cs_str[gt], "+")) {
	  prev_cpg_x = x;
	  prev_cpg_flt = (flt != 0);
	  if(!strncmp(rf_ctxt+2, "CG", 2)) ref_cpg = true;
	  a = counts[5];
	  b = counts[7];
	  cpg_ok = true;
	} else if(!strcmp(cs_str[gt], "-")) {
	  if(!strncmp(rf_ctxt+1, "CG", 2)) ref_cpg = true;
	  if(x - prev_cpg_x == 1) {
	    if(ref_cpg) {
	      stats->CpG_ref[stats_all]++;
	      ctg_stats->CpG_ref[stats_all]++;
	      if(!(prev_cpg_flt || flt)) {
		stats->CpG_ref[stats_passed]++;
		ctg_stats->CpG_ref[stats_passed]++;
	      }
	    } else {
	      ctg_stats->CpG_nonref[stats_all]++;
	      stats->CpG_nonref[stats_all]++;
	      if(!(prev_cpg_flt || flt)) {
		stats->CpG_nonref[stats_passed]++;
		ctg_stats->CpG_nonref[stats_passed]++;
	      }
	    }
	  }
	  a = counts[6];
	  b = counts[4];
	  cpg_ok = true;
	}
	if(cpg_ok) {
	  if(ref_cpg) {
	    stats->qual[CpG_ref_sites][phred]++;
	  } else {
	    stats->qual[CpG_nonref_sites][phred]++;
	  }
	  gcov->CpG[ref_cpg ? 0 : 1]++;
	  gt_cov_stats *gcov1;
	  HASH_FIND(hh, stats->cov_stats, &d_inf, sizeof(uint64_t), gcov1);
	  if(gcov1 == NULL) {
	    gcov1 = calloc((size_t)1, sizeof(gt_cov_stats));
	    gcov1->coverage = d_inf;
	    HASH_ADD(hh, stats->cov_stats, coverage, sizeof(uint64_t), gcov1);
	  }
	  gcov1->CpG_inf[ref_cpg ? 0 : 1]++;
	  if(a + b) {
	    double meth[101];
	    double konst = lfact(a + b + 1) - lfact(a) - lfact(b);
	    double sum = 0.0;
	    if(a) meth[0] = 0.0;
	    else sum = meth[0] = exp(konst);
	    if(b) meth[100] = 0.0;
	    else sum = (meth[100] = exp(konst));
	    double da = (double)a;
	    double db = (double)b;
	    for(int i = 1; i < 100; i++) {
	      sum += (meth[i] = exp(konst + logp[i - 1] * da + logp[99 - i] * db));
	    }
	    for(int i = 0; i < 101; i++) {
	      double z = meth[i] / sum;
	      if(ref_cpg) {
		stats->CpG_ref_meth[stats_all][i] += z;
		if(!flt) stats->CpG_ref_meth[stats_passed][i] += z;
	      } else {
		stats->CpG_nonref_meth[stats_all][i] += z;
		if(!flt) stats->CpG_nonref_meth[stats_passed][i] += z;
	      }
	    }
	  }
	}
      }
      if(mut != mut_no) {
	stats->mut_counts[mut][stats_all]++;
	if(!flt) stats->mut_counts[mut][stats_passed]++;
	if(db_prefix != NULL) { // dbSNP
	  stats->dbSNP_mut_counts[mut][stats_all]++;
	  if(!flt) stats->dbSNP_mut_counts[mut][stats_passed]++;
	}
      }
    }
  }
}

static gt_vcf *vcf;
static uint64_t vcf_x;
static char *vcf_ctg;
static gt_string *ref;
static int vcf_size, vcf_n;
static bool process_end;
static bool print_end;

static char gt_store[5];
uint64_t store_x = 0;
static gt_meth gtm_store[5];
static char *curr_ctg;
static char rf_ctxt[7];

// Print the last 2 entries in gt_store
void flush_vcf_entries(FILE *fp) {
  if (curr_ctg && store_x) {
    for (int i = 0; i < 2; i++) {
      memmove(gt_store, gt_store + 1, 4);
      memmove(gtm_store, gtm_store + 1, 4 * sizeof(gt_meth));
      memmove(rf_ctxt, rf_ctxt + 1, 6);
      if (gt_store[2]) _print_vcf_entry(fp, curr_ctg, gtm_store + 2, rf_ctxt, store_x - 1 + i, gt_store);
    }
    store_x = 0;
  }
}

void print_vcf_entry(FILE *fp, char *ctg, gt_meth *gtm, const char *rf,
                     const uint64_t x, const uint64_t xstart, bool skip) {
  if (curr_ctg != NULL) {
    if (strcmp(curr_ctg, ctg)) {
      //      flush_vcf_entries(fp);
      free(curr_ctg);
      curr_ctg = NULL;
    }
  }
  if (curr_ctg == NULL) {
    uint64_t l = strlen(ctg);
    curr_ctg = gt_malloc(l + 1);
    memcpy(curr_ctg, ctg, l);
    curr_ctg[l] = 0;
    if(param.dbSNP != NULL) {
      HASH_FIND(hh, param.dbSNP, curr_ctg, l, dbSNP_ctg);
    }
    if(param.ctg_stats != NULL) {
      HASH_FIND(hh, param.ctg_stats, curr_ctg, l, ctg_stats);
    }
  }
  uint64_t l = x - store_x;
  if (l < 5) {
    memmove(gt_store, gt_store + l, 5 - l);
    memmove(gtm_store, gtm_store + l, (5 - l) * sizeof(gt_meth));
    for (uint64_t i = 4; i >= 5 - l; i--) {
      gt_store[i] = 0;
    }
  } else memset(gt_store, 0, 5);
  assert(x > store_x);
  store_x = x;
  memcpy(gtm_store + 4, gtm, sizeof(gt_meth));
  if (x - xstart >= 4)
    strncpy(rf_ctxt, rf + x - xstart - 4, 7);
  else {
    uint64_t l = x - xstart;
    for (uint64_t i = 0; i < 4 - l; i++)
      rf_ctxt[i] = 'N';
    strncpy(rf_ctxt + 4 - l, rf, 3 + l);
  }
  if(skip) gt_store[4] = 0;
  else {
    double z = gtm->gt_prob[0];
    int gt = 0;
    for (int i = 1; i < 10; i++)
      if (gtm->gt_prob[i] > z) {
	z = gtm->gt_prob[i];
	gt = i;
      }
    gt_store[4] = gt + 1;
  }
  if (gt_store[2]) _print_vcf_entry(fp, curr_ctg, gtm_store + 2, rf_ctxt, x - 2, gt_store);
}

void *print_thread(void *arg) {
  sr_param *par = arg;
  FILE *fp = par->output_file;
  //	fprintf(stderr,"print_thread()\n");
  while(1) {
    while(!vcf_n && !print_end) usleep(25);
    if(vcf_n) {
      const char *ref_st = gt_string_get_string(ref);
      //			fprintf(stderr,"print_thread() vcf_n = %d, x = %lu, y = %lu\n", vcf_n, vcf_x, vcf_x + vcf_n - 1);
      for(int i = 0; i < vcf_n; i++) {
	while(!vcf[i].ready) usleep(25);
	//				fprintf(stderr,"print_thread() %d / %d, %lu, %c, %c\n", i, vcf_n, vcf_x + i, ref_st[i], vcf[i].skip ? '0' : '1');
	print_vcf_entry(fp, vcf_ctg, &vcf[i].gtm, ref_st, i + vcf_x, vcf_x, vcf[i].skip);
      }
      flush_vcf_entries(fp);
      vcf_n = 0;
      //			fprintf(stderr,"print_thread() finished block\n");
    } else break;
  }
  return NULL;
}

// static char stop_base[256] = { ['A'] = 'a', ['C'] = 'c', ['G'] = 'g', ['T'] = 't', ['.'] = ',' };

typedef struct {
  sr_param *par;
  uint64_t x;
  gt_vector *align_list;
} mprofile_par;

static int mprof_tab[4][4] = {
  {0, 0, 0, 0}, // A in reference
  {0, 1, 0, 2}, // C in reference, informative C is non-converted, informative T is converted
  {2, 0, 1, 0}, // G in reference, informative G is non-converted, informative A is converted
  {0, 0, 0, 0} // T in reference
};

static void *meth_profile_thread(void *arg) {
  mprofile_par  *mpar = arg;
  gt_vector *mprof = mpar->par->stats->meth_profile;
  int min_qual = mpar->par->min_qual;
  uint64_t nr = gt_vector_get_used(mpar->align_list);
  uint64_t x = mpar->x;
  const char *ref_st = gt_string_get_string(ref);
  align_details **al_p = gt_vector_get_mem(mpar->align_list, align_details *);
  int mx = gt_vector_get_used(mprof);
  for (uint64_t ix = 0; ix < nr; ix++, al_p++) {
    align_details *al = *al_p;
    int st = al->bs_strand;
    for (int k = 0; k < 2; k++) {
      if(!al->read[k]) continue;
      uint64_t rl = gt_string_get_length(al->read[k]);
      if(rl == 0) continue;
      // fprintf(stderr,"%d\t%lu\t%lu\n",k,k?al->reverse_position:al->forward_position,rl);
      char *sp = gt_string_get_string(al->read[k]);
      char *sq = gt_string_get_string(al->qualities[k]);
      uint64_t pos = k ? al->reverse_position : al->forward_position;
      const char *rf = ref_st + pos - x;
      uint64_t j;
      int *posx = gt_vector_get_mem(al->orig_pos[k], int);
      char prev_r = pos > x ? base_tab[(int)ref_st[pos - x - 1]] : 0;
      for(j = 0; j < rl; j++) {
        int c = base_tab_st[st][(int)sp[j]];
        int q = sq[j] - QUAL_CONV;
        if (c > 0 && q >= min_qual) {
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
    }
  }
  if(mx > gt_vector_get_used(mprof)) gt_vector_set_used(mprof, mx);
  return NULL;
}

typedef struct {
  pthread_t thr;
  sr_param *param;
  gt_meth *gtm;
  pileup *cts;
  uint64_t x;
  uint64_t y;
  int step;
} cthread_par; 

void *call_thread(void *arg) {
  cthread_par *cpar = arg;
  uint64_t x = cpar->x, y = cpar->y;
  int step = cpar->step;
  gt_meth *tg = cpar->gtm;
  pileup *tp = cpar->cts;
  gt_vcf *v = vcf + x - vcf_x;
  sr_param *par = cpar->param;
  const char *ref_st = gt_string_get_string(ref);
  for (uint64_t i = x; i <= y; i += step, tp += step, tg += step, v += step) {
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
      calc_gt_prob(tg, par, ref_st[i - vcf_x]);
      double fs = 0.0;
      if(gt_het[tg->max_gt]) {
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
	default:
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
	}
	double z = fisher(ftab);
	if(z < 1.0e-20) z = 1.0e-20;
	fs = log(z) / LOG10;
      }
      tg->fisher_strand = fs;			
      memcpy(&v->gtm, tg, sizeof(gt_meth));
      v->skip = false;
    } else v->skip = true;
    v->ready = true;
  }
  return NULL;
}

static gt_vector *pileupv, *gt_resv;

void call_genotypes_ML(char *ctg, gt_vector *align_list, uint64_t x, uint64_t y, sr_param *param) {
  //  char *gt_string[] = {"AA", "AC", "AG", "AT", "CC",
  //                       "CG", "CT", "GG", "GT", "TT"};

  //	FILE *fp_pile = param->pileup ? param->pileup_file : NULL;
  uint64_t nr = gt_vector_get_used(align_list);
  assert(y >= x);
  uint64_t sz = y - x + 1;
  if(pileupv == NULL) pileupv = gt_vector_new(sz, sizeof(pileup));
  else gt_vector_reserve(pileupv, sz, false);
  if(gt_resv == NULL) gt_resv = gt_vector_new(sz, sizeof(gt_meth));
  else gt_vector_reserve(gt_resv, sz, false);
  pileup *counts = gt_vector_get_mem(pileupv, pileup);
  gt_meth *gt_res = gt_vector_get_mem(gt_resv, gt_meth);
  memset(counts, 0, sizeof(pileup) * sz);
  memset(gt_res, 0, sizeof(gt_meth) * sz);
  align_details **al_p = gt_vector_get_mem(align_list, align_details *);
  for (uint64_t ix = 0; ix < nr; ix++, al_p++) {
    align_details *al = *al_p;
    uint64_t x1 = al->forward_position;
    if(x1 == 0) x1 = al->reverse_position;
    else if (al->reverse_position > 0 && al->reverse_position < x1) x1 = al->reverse_position;
    assert(x1 >= x);
    int ori = al->orientation;
    assert(ori < 2);
    int st = al->bs_strand;
    for (int k = 0; k < 2; k++) {
      if(!al->read[k]) continue;
      uint64_t rl = gt_string_get_length(al->read[k]);
      if(rl == 0) continue;
      // fprintf(stderr,"%d\t%lu\t%lu\n",k,k?al->reverse_position:al->forward_position,rl);
      float mapq2 = al->mapq[k] * al->mapq[k];
      char *sp = gt_string_get_string(al->read[k]);
      char *sq = gt_string_get_string(al->qualities[k]);
      uint64_t pos = k ? al->reverse_position : al->forward_position;
      uint64_t j;
      for(j = 0; j < rl; j++) {
        int c = base_tab_st[st][(int)sp[j]];
	if(c) break;
      }
      uint64_t read_start;
      if(j < rl) read_start = j;
      else continue;
      for(j = rl; j > 0; j--) {
        int c = base_tab_st[st][(int)sp[j - 1]];
	if(c) break;
      }
      uint64_t read_end;
      if(j > 0) read_end = j - 1;
      else continue;
      pos += read_start;
      pileup *curr_loc = counts + pos - x;
      for (j = read_start; j <= read_end && pos <= y; j++, pos++, curr_loc++) {
        int c = base_tab_st[st][(int)sp[j]];
        int q = sq[j] - QUAL_CONV;
        if (c-- && q >= param->min_qual) {
          if (q > MAX_QUAL) q = MAX_QUAL;
          curr_loc->n++;
          curr_loc->quality[c] += (float)q;
          curr_loc->mapq2 += mapq2;
          curr_loc->counts[ori][c]++;
        }
      }
      ori ^= 1;
    }
  }
  // Prepare printing
  while(vcf_n) usleep(25);
  if(sz > vcf_size) {
    vcf = realloc(vcf, sizeof(gt_vcf) * sz);
    vcf_size = sz;
  }
  for(int i = 0; i < sz; i++) vcf[i].ready = false;
  vcf_x = x;
  if(vcf_ctg == NULL || strcmp(vcf_ctg, ctg)) {
    if(vcf_ctg != NULL) free(vcf_ctg);
    vcf_ctg = strdup(ctg);
  }
  gt_string_resize(ref, sz + 2);
  //  gt_string_delete(ref);
  //	ref = gt_string_new(sz);
  gt_status gts = gt_sequence_archive_get_sequence_string(param->sequence_archive, ctg, FORWARD, x - 1, sz + 2, ref);
  if (gts) {
    if(gts == GT_SEQUENCE_NOT_FOUND) fprintf(stderr,"Could not find reference sequence for contig '%s'\n", ctg);
    else fprintf(stderr,"No good  - aborting for contig '%s' %" PRIu64 " %" PRIu64 "\n", ctg, x, sz);
    return;
  }
  vcf_n = sz;
  pthread_t meth_profile_thr;
  if(param->stats) {
    mprofile_par mpar = { param, x, align_list };
    pthread_create(&meth_profile_thr, NULL, meth_profile_thread, &mpar);
  }
  // Set up calculation threads
  int nthr = param->num_threads;
  //	int nthr = 1;
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
  if(param->stats) pthread_join(meth_profile_thr, 0);
  free(cpar);
  //	flush_vcf_entries(fp);
  //	if(buf != NULL) free(buf);
}

#define call_genotypes(ctg, align_list, x, y, par) call_genotypes_ML(ctg, align_list, x, y, par);

gt_status process_template_vector(gt_vector *align_list, char *ctg, uint64_t y,
                                  sr_param *param) {
  align_details **al_p = gt_vector_get_mem(align_list, align_details *);
  bs_stats *stats = param->stats;
  int ix = gt_vector_get_used(align_list);
  int orig_ix = ix;
  uint64_t x = 0;
  uint8_t thresh = param->mapq_thresh;
  assert(ix);
  // We sort by both forward and reverse reads as well as by bs_strand to
  // detect duplicates
  qsort(al_p, ix, sizeof(void *), cmp_al);
  // Now we can collapse duplicates, picking the 'best' read pair for analysis
  // in each duplicates set
  // Best being the pair with the highest total quality score for all
  // non-trimmed bases
  x = (*al_p)->forward_position;
  if(x == 0) x = (*al_p)->reverse_position;
  int j = 0, k = 0;
  int curr_ct = 0;
  uint64_t x1 = 0;
  uint64_t x2 = 0;
  uint64_t mx_qual = 0;
  uint8_t mx_mapq = 0;
  gt_bs_strand strand = NON_CONVERTED;
  uint64_t blk_length = 0;
  bool paired = false;
  //  fprintf(stderr, "Processing %" PRIu64 " - %" PRIu64 " (%d pairs)\n", x, y, ix);
  //  fprintf(stderr, "Processing:\t%s\t%lu\t%lu\t%d\n", ctg, x, y, ix);
  for (int i = 0; i < ix; i++) {
    align_details *al =
      *(align_details **)gt_vector_get_elm(align_list, i, align_details *);
    uint64_t len1 = al->read[0] ? al->read[0]->length : 0;
    uint64_t len2 = al->read[1] ? al->read[1]->length : 0;
    if(len1 + len2 == 0) continue;
    bool frag_paired = len1 && len2;		
    if (al->reverse_position > 0 && al->reverse_position < x)
      x = al->reverse_position;
    if (!param->keep_duplicates && al->forward_position == x1 && al->reverse_position == x2 && al->bs_strand == strand) {
      //			fprintf(stderr,"AA: %" PRIu64 ", %" PRIu64 " i = %d\n", al->forward_position, al->reverse_position, i);
      curr_ct++;
      blk_length += len1 + len2;
      uint8_t mapq=0;
      int kn = 0;
      for(int ix1 = 0; ix1 < 2; ix1++) {
	if(al->read[ix1] && al->read[ix1]->length) {
	  mapq += al->mapq[ix1];
	  kn++;
	}
      }
      if(kn) mapq /= kn;
      uint64_t qual = get_al_qual(al);
      if (mapq > mx_mapq) {
        mx_mapq = mapq;
        mx_qual = qual;
        k = i;
      } else if (mapq == mx_mapq && qual > mx_qual) {
        mx_qual = qual;
        k = i;
      }
    } else {
      if (curr_ct > 0) {
	if(mx_mapq >= thresh) {
	  align_details *al1 = *(align_details **)gt_vector_get_elm(align_list, j, align_details *);
	  if (j < k) {
	    align_details *al2 = *(align_details **)gt_vector_get_elm(align_list, k, align_details *);
	    gt_vector_set_elm(align_list, j, align_details *, al2);
	    gt_vector_set_elm(align_list, k, align_details *, al1);
	  }
	  j++;
	  if(stats != NULL) {
	    uint64_t l = (al1->read[0] ? al1->read[0]->length : 0) + (al1->read[1] ? al1->read[1]->length : 0);
	    if(curr_ct > 1) {
	      stats->filter_cts[gt_flt_duplicate] += (curr_ct - 1) * (paired ? 2 : 1);
	      assert(l <= blk_length);
	      stats->filter_bases[gt_flt_duplicate] += blk_length - l;
	    }
	    stats->filter_bases[gt_flt_none] += l;
	    stats->filter_cts[gt_flt_none] += paired ? 2 : 1;
	  }
	} else if(stats != NULL) {
	  stats->filter_cts[gt_flt_mapq] += curr_ct * (paired ? 2 : 1);
	  stats->filter_bases[gt_flt_mapq] += blk_length;
	}
      }
      x1 = al->forward_position;
      x2 = al->reverse_position;
      strand = al->bs_strand;
      curr_ct = 1;
      blk_length = len1 + len2;
      paired = frag_paired;
      mx_qual = get_al_qual(al);
      mx_mapq=0;
      int kn = 0;
      for(int ix1 = 0; ix1 < 2; ix1++) {
	if(al->read[ix1] && al->read[ix1]->length) {
	  mx_mapq += al->mapq[ix1];
	  kn++;
	}
      }
      mx_mapq/= kn;
      k = i;
    }
  }
  if (curr_ct) {
    if(mx_mapq >= thresh) {
      align_details *al1 = *(align_details **)gt_vector_get_elm(align_list, j, align_details *);
      if (j < ix) {
	align_details *al2 =*(align_details **)gt_vector_get_elm(align_list, k, align_details *);
	gt_vector_set_elm(align_list, j++, align_details *, al2);
	gt_vector_set_elm(align_list, k, align_details *, al1);
      }
      if(stats != NULL) {
	uint64_t l = (al1->read[0] ? al1->read[0]->length : 0) + (al1->read[1] ? al1->read[1]->length : 0);
	if(curr_ct > 1) {
	  stats->filter_cts[gt_flt_duplicate] += (curr_ct - 1) * (paired ? 2 : 1);
	  assert(l <= blk_length);
	  stats->filter_bases[gt_flt_duplicate] += blk_length - l;
	}
	stats->filter_bases[gt_flt_none] += l;
	stats->filter_cts[gt_flt_none] += paired ? 2 : 1;
      }
    } else if(stats != NULL) {
      stats->filter_cts[gt_flt_mapq] += curr_ct * (paired ? 2 : 1);
      stats->filter_bases[gt_flt_mapq] += blk_length;
    }
  }
  gt_vector_set_used(align_list, j);
  ix = j;
  // We do trimming by setting to lower case the trimmed bases
  if (param->left_trim || param->right_trim) {
    for (int i = 0; i < ix; i++) {
      align_details *al = *(align_details **)gt_vector_get_elm(align_list, i, align_details *);
      if(al->read[0] && al->read[0]->length) {
	char *sp = gt_string_get_string(al->read[0]);
	uint64_t rl = gt_string_get_length(al->read[0]);
	for (int k1 = 0; k1 < param->left_trim && k1 < rl; k1++) sp[k1] = tolower((int)sp[k1]);
	for (int k1 = 0; k1 < param->right_trim && k1 < rl; k1++)	sp[rl - k1 - 1] = tolower((int)sp[rl - k1 - 1]);
	if(param->right_trim) {
	  int kk = param->right_trim;
	  if(al->read[1] && al->read[1]->length) {
	    kk = 0;
	    if(al->forward_position + al->reference_span[0] >= al->reverse_position + al->reference_span[1]) kk = param->right_trim;
	    else if(al->reverse_position + al->reference_span[1] - al->forward_position - al->reference_span[0] < param->right_trim)
	      kk = param->right_trim - (al->reverse_position + al->reference_span[1] - al->forward_position - al->reference_span[0]);
	  }
	  for (int k1 = 0; k1 < kk && k1 < rl; k1++) sp[rl - k1 - 1] = tolower((int)sp[rl - k1 - 1]);
	}
      }
      if(al->read[1] && al->read[1]->length) {
	char *sp = gt_string_get_string(al->read[1]);
	int rl = gt_string_get_length(al->read[1]);
	for (int k1 = 0; k1 < param->left_trim && k1 < rl; k1++) sp[rl - k1 - 1] = tolower((int)sp[rl - k1 - 1]);
	for (int k1 = 0; k1 < param->right_trim && k1 < rl; k1++) sp[k1] = tolower((int)sp[k1]);
	if(param->left_trim) {
	  int kk = param->left_trim;
	  if(al->read[0] && al->read[0]->length) {
	    kk = 0;
	    if(al->reverse_position <= al->forward_position) kk = param->left_trim;
	    else if(al->reverse_position - al->forward_position < param->left_trim)
	      kk = param->left_trim - (al->reverse_position - al->forward_position);
	  }
	  for (int k1 = 0; k1 < kk && k1 < rl; k1++) sp[k1] = tolower((int)sp[k1]);
	}
      }
    }
  }

  assert(x > 0 && x <= y);
  x1 = x;
  uint64_t y1 = y;
  if(x1 > 2) x1 -= 2;
  else x1 = 1;
  y1 += 2;
  x=x1;
  y=y1;
  al_p = gt_vector_get_mem(align_list, align_details *);
  for (int i = 0; i < ix; i++, al_p++) {
    for (int k = 0; k < 2; k++) {
      if((*al_p)->read[k]==NULL) continue;
      uint64_t rl = gt_string_get_length((*al_p)->read[k]);
      if(rl == 0) continue;
      uint64_t num_misms = gt_vector_get_used((*al_p)->mismatches[k]);
      // Trim Soft clips if present
      int nclip = 0;
      uint64_t adj = 0;
      for (int z = 0; z < num_misms; z++) {
        gt_misms *misms =
	  gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
        if (misms->misms_type == SOFT) {
          if (z && z != num_misms - 1)
            gt_fatal_error_msg(
			       "Error in CIGAR - Soft clip not at extremity of read\n");
          nclip++;
          if (!misms->position) {
            if (misms->size >= rl)
              gt_fatal_error_msg(
				 "Error in CIGAR - Illegal soft clip (%d %" PRIu64 " %" PRIu64 " %" PRIu64 ")\n", z,
				 misms->position, misms->size, rl);
            adj = misms->size;
	    if(stats) stats->base_filter[base_clip] += adj;
            gt_string_trim_left((*al_p)->read[k], adj);
            gt_string_trim_left((*al_p)->qualities[k], adj);
	    (*al_p)->trim_left[k] = adj;
          } else {
            if (misms->position + misms->size != rl)
              gt_fatal_error_msg(
				 "Error in CIGAR - Illegal soft clip (%d %" PRIu64 " %" PRIu64 " %" PRIu64 ")\n", z,
				 misms->position, misms->size, rl);
            gt_string_trim_right((*al_p)->read[k], misms->size);
            gt_string_trim_right((*al_p)->qualities[k], misms->size);
	    (*al_p)->trim_right[k] = misms->size;
	    if(stats) stats->base_filter[base_clip] += misms->size;
          }
        } else if (nclip) {
          misms->position -= adj;
          gt_vector_set_elm((*al_p)->mismatches[k], z - nclip, gt_misms, *misms);
        }
      }
      if (nclip) {
        num_misms -= nclip;
        gt_vector_set_used((*al_p)->mismatches[k], num_misms);
        rl = gt_string_get_length((*al_p)->read[k]);
      }
    }
    uint64_t rdl[2]={0,0};
    if((*al_p)->read[0]) rdl[0]=gt_string_get_length((*al_p)->read[0]);
    if((*al_p)->read[1]) rdl[1]=gt_string_get_length((*al_p)->read[1]);
    bool rev;
    int64_t overlap;
    if(rdl[0] > 0 && rdl[1] > 0) {
      if ((*al_p)->forward_position <= (*al_p)->reverse_position) {   
	overlap = (*al_p)->reference_span[0] - (*al_p)->reverse_position + (*al_p)->forward_position;
	rev = false;
      } else {
	overlap = (*al_p)->reference_span[1] + (*al_p)->reverse_position - (*al_p)->forward_position;
	rev = true;
      }
      // Look for overlapping reads - keep only best part
      if ((*al_p)->forward_position + (*al_p)->reference_span[0] >=	(*al_p)->reverse_position) {
	// Find the overlapping parts
	uint64_t *rspan;
	// Find best quality read (normally read 1)
	rspan = (*al_p)->reference_span;
	int trim_read;
	if (rspan[0] > rspan[1])
	  trim_read = 1;
	else if (rspan[0] < rspan[1])
	  trim_read = 0;
	else {
	  uint64_t tot[2] = {0, 0};
	  for (int k = 0; k < 2; k++) {
	    char *sq = gt_string_get_string((*al_p)->qualities[k]);
	    for (int i = 0; i < rdl[k]; i++)
	      tot[k] += sq[i];
	  }
	  trim_read = tot[0] < tot[1] ? 0 : 1;
	}
	// Adjust starting position if a left trim is used
	if (!(rev && trim_read) && (rev || trim_read)) {
	  if (trim_read)
	    (*al_p)->reverse_position += overlap;
	  else
	    (*al_p)->forward_position += overlap;
	}
	// Find overlap point for read
	uint64_t num_misms = gt_vector_get_used((*al_p)->mismatches[trim_read]);
	if (!num_misms) {
	  if ((rev && trim_read) || !(rev || trim_read)) {
	    gt_string_trim_right((*al_p)->read[trim_read], overlap);
	    gt_string_trim_right((*al_p)->qualities[trim_read], overlap);
	  } else {
	    gt_string_trim_left((*al_p)->read[trim_read], overlap);
	    gt_string_trim_left((*al_p)->qualities[trim_read], overlap);	
	  }
	} else {
	  bool trimmed = false;
	  if ((rev && trim_read) || !(rev || trim_read)) {
	    uint64_t xx = (*al_p)->reference_span[trim_read] - overlap;
	    int64_t adj = 0;
	    int z;
	    for (z = 0; z < num_misms; z++) {
	      gt_misms *misms =
                gt_vector_get_elm((*al_p)->mismatches[trim_read], z, gt_misms);
	      if (misms->position + adj >= xx) {
		int64_t trim = rdl[trim_read] - xx + adj;
		gt_string_trim_right((*al_p)->read[trim_read], trim);
		gt_string_trim_right((*al_p)->qualities[trim_read], trim);	
		num_misms = z;
		gt_vector_set_used((*al_p)->mismatches[trim_read], num_misms);
		trimmed = true;
		break;
	      }
	      if (misms->misms_type == INS) {
		if (misms->position + adj + misms->size >= xx) {
		  int64_t trim = rdl[trim_read] - misms->position;
		  misms->size = xx - (misms->position + adj);
		  gt_string_trim_right((*al_p)->read[trim_read], trim);
		  gt_string_trim_right((*al_p)->qualities[trim_read], trim);
		  num_misms = z + 1;
		  gt_vector_set_used((*al_p)->mismatches[trim_read], num_misms);
		  trimmed = true;
		}
		adj += misms->size;
	      } else if (misms->misms_type == DEL) {
		adj -= misms->size;
	      }
	    }
	    if (trimmed == false) {
	      gt_string_trim_right((*al_p)->read[trim_read], overlap);
	      gt_string_trim_right((*al_p)->qualities[trim_read], overlap);
	    }
	  } else {
	    uint64_t xx = overlap;
	    int z;
	    int64_t adj = 0;
	    for (z = 0; z < num_misms; z++) {
	      gt_misms *misms =
                gt_vector_get_elm((*al_p)->mismatches[trim_read], z, gt_misms);
	      if (misms->position + adj >= xx) {
		uint64_t trim = overlap - adj;
		gt_string_trim_left((*al_p)->read[trim_read], trim);
		gt_string_trim_left((*al_p)->qualities[trim_read], trim);	
		trimmed = true;
		if (z) {
		  for (int z1 = z; z1 < num_misms; z1++) {
		    gt_misms *misms = gt_vector_get_elm((*al_p)->mismatches[trim_read], z1, gt_misms);
		    misms->position -= trim;
		    gt_misms *misms1 = gt_vector_get_elm((*al_p)->mismatches[trim_read], z1 - z, gt_misms);
		    gt_vector_set_elm((*al_p)->mismatches[trim_read], z1 - z,gt_misms, *misms);
		    gt_vector_set_elm((*al_p)->mismatches[trim_read], z1,gt_misms, *misms1);
		  }
		  num_misms -= z;
		  gt_vector_set_used((*al_p)->mismatches[trim_read], num_misms);
		  break;
		} else {
		  for (int z1 = 0; z1 < num_misms; z1++) {
		    gt_misms *misms = gt_vector_get_elm((*al_p)->mismatches[trim_read], z1, gt_misms);
		    misms->position -= trim;
		  }
		}
		break;
	      }
	      if (misms->misms_type == INS) {
		if (misms->position + adj + misms->size >= xx) {
		  misms->size = misms->position + misms->size + adj - xx;
		  uint64_t trim = misms->position;
		  gt_string_trim_left((*al_p)->read[trim_read], trim);
		  gt_string_trim_left((*al_p)->qualities[trim_read], trim);	
		  trimmed = true;
		  int z2 = misms->size ? z : z + 1;
		  for (int z1 = z2; z1 < num_misms; z1++) {
		    gt_misms *misms = gt_vector_get_elm((*al_p)->mismatches[trim_read], z1, gt_misms);
		    misms->position -= trim;
		    if (z2) {
		      gt_misms *misms1 = gt_vector_get_elm((*al_p)->mismatches[trim_read], z1 - z2, gt_misms);
		      gt_vector_set_elm((*al_p)->mismatches[trim_read], z1 - z2,gt_misms, *misms);
		      gt_vector_set_elm((*al_p)->mismatches[trim_read], z1, gt_misms, *misms1);
		    }
		  }
		  num_misms -= z2;
		  gt_vector_set_used((*al_p)->mismatches[trim_read], num_misms);
		  break;
		}
		adj += misms->size;
	      } else if (misms->misms_type == DEL) {
		adj -= misms->size;
	      }
	    }
	    if (trimmed == false) {
	      gt_string_trim_left((*al_p)->read[trim_read], overlap - adj);
	      gt_string_trim_left((*al_p)->qualities[trim_read], overlap - adj);
	      gt_vector_set_used((*al_p)->mismatches[trim_read], 0);
	    }
	  }
	}
	uint64_t rdl1[2];
	rdl1[0] = (*al_p)->read[0] ? gt_string_get_length((*al_p)->read[0]) : 0;
	rdl1[1] = (*al_p)->read[1] ? gt_string_get_length((*al_p)->read[1]) : 0;
	if(stats) stats->base_filter[base_overlap] += (rdl[0] - rdl1[0] + rdl[1] - rdl1[1]);
	if ((rev && trim_read) || !(rev || trim_read)) { // Right trim
	  (*al_p)->trim_right[trim_read] += (rdl[trim_read] - rdl1[trim_read]);
	} else { // Left trim
	  (*al_p)->trim_left[trim_read] += (rdl[trim_read] - rdl1[trim_read]);
	}
	rdl[0] = rdl1[0];
	rdl[1] = rdl1[1];
      }
    }
    if(stats) {
      for(int k = 0; k < 2; k++) {
	if((*al_p)->read[k]) {
	  char *sp = gt_string_get_string((*al_p)->read[k]);
	  char *sq = gt_string_get_string((*al_p)->qualities[k]);
	  for(int k1 = 0; k1 < rdl[k]; k1++) {
	    char c = *sp++;
	    char q = (*sq++) - QUAL_CONV;
	    if(islower((int)c)) stats->base_filter[base_trim]++;
	    else if(c == 'N' || q < param->min_qual) stats->base_filter[base_lowqual]++;
	    else stats->base_filter[base_none]++;
	  }
	}
      }
    }
	
    // Quick and dirty - to be replaced in next version!
    // Here we simply normalize the reads w.r.t. indels by
    // adding N's in the place of deletions and removing insertions
    for (int k = 0; k < 2; k++) {
      if(rdl[k] == 0) continue;
      uint64_t num_misms = gt_vector_get_used((*al_p)->mismatches[k]);
      if (num_misms) {
        uint64_t del_size = 0;
        for (int z = 0; z < num_misms; z++) {
          gt_misms *misms =
	    gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
          if (misms->misms_type == INS) {
            del_size += misms->size;
          }
        }
        char *tp, *tq, *sp, *sq;
        if (del_size) {
          gt_string_resize((*al_p)->read[k], rdl[k] + del_size);
          gt_string_resize((*al_p)->qualities[k], rdl[k] + del_size);
          sp = gt_string_get_string((*al_p)->read[k]);
          sq = gt_string_get_string((*al_p)->qualities[k]);
          tp = gt_malloc(rdl[k] * 2);
          tq = tp + rdl[k];
          memcpy(tp, sp, rdl[k]);
          memcpy(tq, sq, rdl[k]);
        } else {
          tp = sp = gt_string_get_string((*al_p)->read[k]);
          tq = sq = gt_string_get_string((*al_p)->qualities[k]);
        }
	if((*al_p)->orig_pos[k]) gt_vector_reserve((*al_p)->orig_pos[k], rdl[k] + del_size, false);
	else (*al_p)->orig_pos[k] = gt_vector_new(rdl[k] + del_size, sizeof(int));
        uint64_t ix = 0, ix1 = 0;
	int posx = k ? rdl[k] + (*al_p)->trim_right[k] - 1 : (*al_p)->trim_left[k];
	int *orig = gt_vector_get_mem((*al_p)->orig_pos[k], int);
        for (int z = 0; z < num_misms; z++) {
          gt_misms *misms =  gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
          uint64_t s = misms->position - ix;
          memcpy(sp + ix1, tp + ix, s);
          memcpy(sq + ix1, tq + ix, s);
	  if(k) {
	    for(int k1 = 0; k1 < s; k1++) orig[ix1 + k1] = posx - ix - k1;
	  } else {
	    for(int k1 = 0; k1 < s; k1++) orig[ix1 + k1] = posx + ix + k1;
	  }
          if (misms->misms_type == INS) {
            ix = misms->position;
            ix1 += s;
            for (int k1 = 0; k1 < misms->size; k1++, ix1++) {
              sp[ix1] = 'N';
              sq[ix1] = 33;
	      orig[ix1] = -1;
            }
          } else if (misms->misms_type == DEL) {
            ix = misms->position + misms->size;
            ix1 += s;
          }
        }
        uint64_t s = rdl[k] - ix;
        memcpy(sp + ix1, tp + ix, s);
        memcpy(sq + ix1, tq + ix, s);
	if(k) {
	  for(int k1 = 0; k1 < s; k1++) orig[ix1 + k1] = posx - ix - k1;
	} else {
	  for(int k1 = 0; k1 < s; k1++) orig[ix1 + k1] = posx + ix + k1;
	}
        ix1 += s;
        gt_string_set_length((*al_p)->read[k], ix1);
        gt_string_set_length((*al_p)->qualities[k], ix1);
	gt_vector_set_used((*al_p)->orig_pos[k], ix1);
        if (del_size) gt_free(tp);
      } else {
	if((*al_p)->orig_pos[k]) gt_vector_reserve((*al_p)->orig_pos[k], rdl[k], false);
	else (*al_p)->orig_pos[k] = gt_vector_new(rdl[k], sizeof(int));
	int posx = k ? rdl[k] + (*al_p)->trim_right[k] - 1 : (*al_p)->trim_left[k];
	int *orig = gt_vector_get_mem((*al_p)->orig_pos[k], int);
	if(k) {
	  for(int k1 = 0; k1 < rdl[k]; k1++) orig[k1] = posx - k1;
	} else {
	  for(int k1 = 0; k1 < rdl[k]; k1++) orig[k1] = posx + k1;
	}				
      }
    }
  }
  if (ix) {
    call_genotypes(ctg, align_list, x, y, param);
  }
  //  gt_string_delete(ref);
  gt_vector_set_used(align_list, orig_ix);
  return GT_STATUS_OK;
}


static gt_vector *align_list_waiting;
static gt_vector *free_list_waiting;
static char *ctg_waiting;
static uint64_t y_waiting;

void *process_thread(void *arg) {
  sr_param *par = arg;
  gt_vector *prev_align = NULL;
  while(true) {
    while(align_list_waiting == NULL && !process_end) usleep(25);
    gt_vector *alist = align_list_waiting;
    if(alist != NULL) {
      char *ctg = ctg_waiting;
      uint64_t pos = y_waiting;
      //      			fprintf(stderr, "Process_thread(): ctg_waiting = %s (%p), alist = %p (%lu, %lu), prev_align = %p (%lu, %lu), pos = %lu\n", ctg_waiting, ctg_waiting, alist, alist ? alist->elements_allocated: 0,
      // alist ? alist->used : 0, prev_align, prev_align ? prev_align->elements_allocated : 0, prev_align ? prev_align->used : 0, pos);
      ctg_waiting = NULL;
      free_list_waiting = prev_align;
      align_list_waiting = NULL;
      process_template_vector(alist, ctg, pos, par);
      free(ctg);
      prev_align = alist;
    } else break;
  }
  //	flush_vcf_entries(par->output_file);
  return NULL;
}
	
// Expects a SAM file sorted on coordinates.  Selects mapped pairs with TQ and
// TQ >= thresh, and make up vector of overlapping templates
gt_status input_sam_parser_get_template_vector(
					       gt_buffered_input_file *const buffered_sam_input, gt_vector *align_list,
					       gt_sam_parser_attributes *const sam_parser_attr, sr_param *param) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_VECTOR_CHECK(align_list);
  gt_vector_clear(align_list);
  gt_status error_code;
  char *curr_ctg = NULL, *old_ctg = NULL;
  bool chr_skip = false;
  bool new_contig = false;
  uint64_t max_pos = 0; // Position of righmost end of current pileup
  uint64_t start_pos = 0; // Position of leftmost end of current pileup
  uint64_t read_idx = 0;
  //	uint64_t n_align = 0;
  gt_string *tag = gt_string_new(128);
  align_details *al = 0, *align_hash = 0;
  gt_vector *free_list = gt_vector_new(256, sizeof(align_details *));
  // Cycle through input lines until next read does not overlap with current
  // pileup
  do {
    // Check the end_of_block. Reload buffer if needed
    if (gt_buffered_input_file_eob(buffered_sam_input)) {
      if ((error_code = gt_input_sam_parser_reload_buffer(buffered_sam_input)) != GT_ISP_OK) {
	// Insert unmatched reads into vector
	align_details *thash = align_hash;
	while (thash) {
	  gt_vector_reserve(align_list, thash->idx + 1, false);
	  if (gt_vector_get_used(align_list) <= thash->idx)
	    gt_vector_set_used(align_list, thash->idx + 1);
	  gt_vector_set_elm(align_list, thash->idx, align_details *, thash);
	  align_details *thash1 = thash->hh.next;
	  HASH_DEL(align_hash, thash);
	  thash = thash1;
	}
        int ix = gt_vector_get_used(align_list);
        if (ix) {
	  while(align_list_waiting != NULL) usleep(25);
	  y_waiting = max_pos;
	  ctg_waiting = strdup(curr_ctg);
	  align_list_waiting = align_list;
	  process_end = true;
	  //          process_template_vector(align_list, curr_ctg, max_pos, param);
        }
        return error_code;
      }
    }
    // Check file format
    gt_input_file *input_file = buffered_sam_input->input_file;
    if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
      gt_error(PARSE_SAM_BAD_FILE_FORMAT, input_file->file_name,
               buffered_sam_input->current_line_num, (uint64_t)0);
      return GT_ISP_FAIL;
    }
    char *const line_start = buffered_sam_input->cursor;
    const uint64_t line_num = buffered_sam_input->current_line_num;
    // Parse alignment
    const char **const text_line =
      (const char **const) & (buffered_sam_input->cursor);
    // Read initial TAG (QNAME := Query template)
    if ((error_code = gt_isp_read_tag(text_line, text_line, tag)))
      return error_code;
    // Parse SAM Alignment
    bool clear_al = true;
    if (!al) {
      if (gt_vector_get_used(free_list)) {
        al = *(gt_vector_get_last_elm(free_list, align_details *));
        gt_vector_dec_used(free_list);
      } else {
	//				n_align++;
        al = gt_malloc(sizeof(align_details));
        al->seq_name = al->read[0] = al->read[1] = al->qualities[0] = al->qualities[1] = NULL;
	al->orig_pos[0] = al->orig_pos[1] = NULL;
        al->mismatches[0] = gt_vector_new(8, sizeof(gt_misms));
        al->mismatches[1] = gt_vector_new(8, sizeof(gt_misms));
	al->tag = gt_string_new(32);
	clear_al = false;
      }
    }
    if(clear_al) {
      gt_string_clear(al->tag);
      gt_string_clear(al->seq_name);
      if (al->read[0]) {
	gt_string_clear(al->read[0]);
	gt_string_clear(al->qualities[0]);
	gt_vector_clear(al->mismatches[0]);
      }
      if (al->read[1]) {
	gt_string_clear(al->read[1]);
	gt_string_clear(al->qualities[1]);
	gt_vector_clear(al->mismatches[1]);
      }
      if(al->orig_pos[0]) gt_vector_clear(al->orig_pos[0]);
      if(al->orig_pos[1]) gt_vector_clear(al->orig_pos[1]);
      gt_vector_clear(al->mismatches[0]);
      gt_vector_clear(al->mismatches[1]);
    }
    al->mapq[0] = al->mapq[1] = 0;
    al->trim_left[0] = al->trim_left[1] = al->trim_right[0] = al->trim_right[1] = 0;
    al->forward_position = al->reverse_position = 0;
    bool reverse;
    error_code = gt_isp_quick_parse_bs_sam_alignment(
						     text_line, al, param->mapq_thresh, param->max_template_len, param->keep_unmatched, &reverse);
    gt_input_sam_parser_next_record(buffered_sam_input);
    if (error_code) {
      if (error_code == GT_ISP_SAM_FILTERED) {
	if(param->stats != NULL) {
	  param->stats->filter_cts[al->filtered]++;
	  int ix = reverse ? 1 : 0;
	  uint64_t l = al->read[ix] ? al->read[ix]->length : 0;
	  param->stats->filter_bases[al->filtered] += l;
	  //				fprintf(stderr,"A: FILTERED!\n");]
	}
        continue;
      }
      gt_input_sam_parser_prompt_error(buffered_sam_input, line_num,
                                       buffered_sam_input->cursor - line_start,
                                       error_code);
      return GT_ISP_FAIL;
    }
    gt_input_parse_tag_chomp_pairend_info(tag);
    gt_string_copy(al->tag, tag);
    bool new_block = false;
    new_contig = false;
    if (!curr_ctg || strncmp(curr_ctg, gt_string_get_string(al->seq_name),
                             gt_string_get_length(al->seq_name))) {
      chr_skip = false;
      old_ctg = curr_ctg;
      new_contig = true;
      uint64_t l = gt_string_get_length(al->seq_name);
      curr_ctg = gt_malloc(l + 1);
      memcpy(curr_ctg, gt_string_get_string(al->seq_name), l);
      curr_ctg[l] = 0;
      if (param->species_hash) {
        ctg_hash *tp;
        HASH_FIND_STR(param->species_hash, curr_ctg, tp);
        if (!tp)
          chr_skip = true;
      }
      fprintf(stderr, "Processing chromosome %s (%s)\n", curr_ctg, chr_skip ? "SKIP" : "OK");
      new_block = true;
    }
    if (chr_skip) {
      if(old_ctg) {
	gt_free(old_ctg);
	old_ctg = NULL;
      }
      continue;
    }
    bool insert;
    if((al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) && al->forward_position > 0 && al->reverse_position > 0) {
      if (reverse)
	insert = al->forward_position > al->reverse_position;
      else
	insert = al->forward_position < al->reverse_position;
      if (al->forward_position == al->reverse_position) {
	insert = false;
	if (!new_block) {
	  align_details *thash;
	  HASH_FIND(hh, align_hash, gt_string_get_string(al->tag), gt_string_get_length(al->tag), thash);
	  if (!thash)
	    insert = true;
	} else
	  insert = true;
      }
    } else {
      insert = true;
    }
    if (new_block == false && insert == true) {
      if(start_pos > 0) {
	if(al->forward_position > 0) {
	  if(al->forward_position > max_pos && (al->reverse_position > max_pos || al->reverse_position == 0)) {
	    if(al->forward_position - max_pos > 8) new_block = true;
	  }
	} else if(al->reverse_position > max_pos && al->reverse_position - max_pos > 8) new_block = true;
      }
    }
    if (new_block == true) {
      //			fprintf(stderr,"New Block: max_pos = %lu, fp = %lu, rp = %lu\n", max_pos, al->forward_position, al->reverse_position);
      // Insert unmatched reads into vector
      align_details *thash = align_hash;
      while (thash) {
        gt_vector_reserve(align_list, thash->idx + 1, false);
        if (gt_vector_get_used(align_list) <= thash->idx)
          gt_vector_set_used(align_list, thash->idx + 1);
        gt_vector_set_elm(align_list, thash->idx, align_details *, thash);
        align_details *thash1 = thash->hh.next;
        HASH_DEL(align_hash, thash);
        thash = thash1;
      }
      insert = true;
      read_idx = 0;
      int ix = gt_vector_get_used(align_list);
      if (ix) {
	align_details **al_p = gt_vector_get_mem(align_list, align_details *);
	uint64_t xx = (*al_p)->forward_position;
	if(xx == 0) xx = (*al_p)->reverse_position;
	while(align_list_waiting != NULL) usleep(25);
	ctg_waiting = strdup(new_contig ? old_ctg : curr_ctg);
	y_waiting = max_pos;
	gt_vector *new_free_list = free_list_waiting;
	//				uint64_t t_align = ix + gt_vector_get_used(free_list) + (new_free_list != NULL ? gt_vector_get_used(new_free_list) : 0);
	//				fprintf(stderr, "Parser(): n_align = %lu, t_align = %lu, ix = %d, ctg_waiting = %s (%p), align_list = %p, free_list_waiting = %p, max_pos = %lu\n",
	//								n_align, t_align, ix, ctg_waiting, ctg_waiting, align_list, free_list_waiting, max_pos);
	free_list_waiting = NULL;
	align_list_waiting = align_list;
	if(new_free_list != NULL) {
	  //        process_template_vector(align_list, new_contig ? old_ctg : curr_ctg, max_pos, param);
	  uint64_t used = gt_vector_get_used(free_list);
	  ix = gt_vector_get_used(new_free_list);
	  gt_vector_reserve(free_list, ix + used, false);
	  memcpy(gt_vector_get_mem(free_list, align_details *)+used,
		 gt_vector_get_mem(new_free_list, align_details *),
		 sizeof(align_details *) * ix);
	  gt_vector_set_used(free_list, ix + used);
	  align_list = new_free_list;
	  gt_vector_clear(align_list);
	} else align_list =  gt_vector_new(32, sizeof(align_details *));
      }
      if(new_contig && old_ctg) {
	gt_free(old_ctg);
	old_ctg = NULL;
	new_contig = false;
      }
      //     fprintf(stderr, "Reading new block\n");
      max_pos = start_pos = 0;
    }
    uint64_t x = al->forward_position;
    uint64_t x1 = al->reverse_position;
    if(x >= x1) {
      assert(x > 0);
      if(x + al->align_length > max_pos) max_pos = x + al->align_length;
      if(x1 > 0) {
	if(start_pos == 0 || x1 < start_pos) start_pos = x1;
      } else if(start_pos == 0 || x < start_pos) start_pos = x;
    } else {
      if(x1 + al->align_length > max_pos) max_pos = x1 + al->align_length;
      if(x > 0) {
	if(start_pos == 0 || x < start_pos) start_pos = x;
      } else if(start_pos == 0 || x1 < start_pos) start_pos = x1;
    }
    if(al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) {
      //			fprintf(stderr,"OOOK: %lu %lu %lu %lu\n", al->forward_position, al->reverse_position, start_pos, max_pos);
      if (insert == false) {
	align_details *thash;
	HASH_FIND(hh, align_hash, gt_string_get_string(al->tag), gt_string_get_length(al->tag), thash);
	if (thash) {
	  gt_vector_reserve(align_list, thash->idx + 1, false);
	  if (gt_vector_get_used(align_list) <= thash->idx)
	    gt_vector_set_used(align_list, thash->idx + 1);
	  gt_vector_set_elm(align_list, thash->idx, align_details *, thash);
	  HASH_DEL(align_hash, thash);
	  int ix = reverse ? 1 : 0;
	  gt_string *ts = thash->read[ix];
	  thash->read[ix] = al->read[ix];
	  al->read[ix] = ts;
	  ts = thash->qualities[ix];
	  thash->qualities[ix] = al->qualities[ix];
	  al->qualities[ix] = ts;
	  thash->mapq[ix] = al->mapq[ix];
	  thash->reference_span[ix] = al->reference_span[ix];
	  gt_vector *tv = thash->mismatches[ix];
	  thash->mismatches[ix] = al->mismatches[ix];
	  al->mismatches[ix] = tv;
	  assert(al->forward_position = thash->forward_position && al->reverse_position = thash->reverse_position);
	} else {
	  if(param->stats) param->stats->filter_cts[14]++;
	  if(param->keep_unmatched) {
	    if(al->forward_position > 0) x = al->forward_position + al->align_length;
	    else x = al->reverse_position + al->align_length;
	    if(x > max_pos) max_pos = x;
	    al->idx = read_idx++; // Preserve read order after matching
	    gt_vector_reserve(align_list, al->idx + 1, false);
	    if (gt_vector_get_used(align_list) <= al->idx)
	      gt_vector_set_used(align_list, al->idx + 1);
	    gt_vector_set_elm(align_list, al->idx, align_details *, al);
	    al = 0;
	  } else {
		  fprintf(stdout, "Warning not found: " PRIgts " %" PRIu64 " %" PRIu64 " %c\n",
					 PRIgts_content(tag), al->forward_position, al->reverse_position,
					 al->orientation == FORWARD ? '+' : '-');
	  }
	}
      } else {
	// Here we have a forward facing pair, so we need to store end to be
	// matched up later
	al->idx = read_idx++; // Preserve read order after matching
	align_details *thash;
	HASH_FIND(hh, align_hash, gt_string_get_string(al->tag), gt_string_get_length(al->tag), thash);
	gt_cond_fatal_error(thash != NULL, PARSE_SAM_DUPLICATE_SEQUENCE_TAG,
			    PRIgts_content(tag));
	HASH_ADD_KEYPTR(hh, align_hash, gt_string_get_string(al->tag), gt_string_get_length(al->tag), al);
	al = 0;
      }
    } else { // Single (non-paired) reads
      al->idx = read_idx++; // Preserve read order after matching
      gt_vector_reserve(align_list, al->idx + 1, false);
      if (gt_vector_get_used(align_list) <= al->idx)
	gt_vector_set_used(align_list, al->idx + 1);
      gt_vector_set_elm(align_list, al->idx, align_details *, al);
      al = 0;
    }
  } while (1);
  return GT_ISP_OK;
}

#define GT_INPUT_MULTIFASTA_RETURN_ERROR(error_code)	\
  gt_string_delete(buffer);				\
  gt_segmented_sequence_delete(seg_seq);		\
  return GT_IFP_PE_TAG_BAD_BEGINNING


GT_INLINE gt_status gt_input_multifasta(gt_input_file* const input_multifasta_file,gt_sequence_archive* const sequence_archive) {
  static int tbin[256] = { ['A'] = 1, ['C'] = 2, ['G'] = 2, ['T'] = 1, ['a'] = 1, ['c'] = 2, ['g'] = 2, ['t'] = 1};			
  GT_INPUT_FILE_CHECK(input_multifasta_file);
  GT_SEQUENCE_ARCHIVE_CHECK(sequence_archive);
  if (input_multifasta_file->eof) return GT_IFP_OK;
  GT_INPUT_FILE_CHECK(input_multifasta_file);
  // Read all sequences
  gt_string* const buffer = gt_string_new(200); // TODO: Should be done in terms of dna_string
  bool calc_gc = param.report_file == NULL ? false : true;
  gt_vector *bins = calc_gc ? gt_vector_new(16374, sizeof(uint8_t)) : NULL;
  gt_ctg_stats *gs = NULL;
  int ix = 0, ct[2] = {0, 0};
  while (!input_multifasta_file->eof) {
    gt_segmented_sequence* seg_seq = gt_segmented_sequence_new();
    // Parse TAG
    if (GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)!='>') {
      GT_INPUT_MULTIFASTA_RETURN_ERROR(GT_IFP_PE_TAG_BAD_BEGINNING);
    }
    GT_INPUT_FILE_NEXT_CHAR(input_multifasta_file);
    bool store_str = true;
    while (gt_expect_true(!input_multifasta_file->eof &&
			  GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)!='*' &&
			  GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)!='=' &&
			  GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)!=EOL &&
			  GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)!=DOS_EOL)) {
      if(store_str) {
	if(isspace(GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file))) store_str = false;
	else gt_string_append_char(seg_seq->seq_name,GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file));
      }
      GT_INPUT_FILE_NEXT_CHAR(input_multifasta_file);
    }
    GT_INPUT_FILE_SKIP_EOL(input_multifasta_file);
    gt_string_append_eos(seg_seq->seq_name);
    if(calc_gc) {
      HASH_FIND(hh, param.ctg_stats, seg_seq->seq_name->buffer, seg_seq->seq_name->length, gs);
      if(gs != NULL) {
	fprintf(stderr,"Duplicate contig found in reference\n");
	return GT_STATUS_FAIL;
      }
      gs = calloc((size_t)1, sizeof(gt_ctg_stats));
      gs->contig = malloc(seg_seq->seq_name->length + 1);
      gs->nbins = 0;
      memcpy(gs->contig, seg_seq->seq_name->buffer, seg_seq->seq_name->length);
      gs->contig[seg_seq->seq_name->length] = 0;
      HASH_ADD_KEYPTR(hh, param.ctg_stats, gs->contig, seg_seq->seq_name->length, gs);
      gt_vector_clear(bins);
      ix = ct[0] = ct[1] = 0;
    }
    // Parse all the sequence lines
    while (gt_expect_true(!input_multifasta_file->eof)) {
      if (gt_expect_true(gt_is_iupac_code(GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)))) {
	// TODO: gt_get_dna_normalized should be transparent
	char b = gt_get_dna_normalized(GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file));
	gt_string_append_char(buffer, b);
	if(calc_gc) {
	  b = tbin[(int)b];
	  if(b) ct[b - 1]++;
	  ix++;
	  if(ix == 100) {
	    *(gt_vector_get_free_elm(bins, uint8_t)) = (ct[0] + ct[1] == 100 ? ct[1] : 255);
	    gt_vector_reserve_additional(bins, 1);
	    gt_vector_inc_used(bins);
	    ix = ct[0] = ct[1] = 0;
	  }
	}
	GT_INPUT_FILE_NEXT_CHAR(input_multifasta_file);
      } else if (GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)==EOL ||
		 GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)==DOS_EOL) {
	// Append string and clear buffer (next line!)
	gt_segmented_sequence_append_string(seg_seq, gt_string_get_string(buffer),gt_string_get_length(buffer));
	gt_string_clear(buffer);
	GT_INPUT_FILE_SKIP_EOL(input_multifasta_file);
      } else if (GT_INPUT_FILE_CURRENT_CHAR(input_multifasta_file)=='>') {
	break; // Next sequence !
      }
    } 
    // Store the parsed sequence
    gt_sequence_archive_add_segmented_sequence(sequence_archive,seg_seq);
    if(calc_gc) {
      gs->nbins = gt_vector_get_used(bins);
      if(gs->nbins) {
	gs->gc = malloc((size_t)gs->nbins);
	memcpy(gs->gc, bins->memory, (size_t)gs->nbins);
      }
      //			fprintf(stderr,"Adding ctg '%s' (nbins = %d)\n", gs->contig, gs->nbins);			
    }
  }
  gt_string_delete(buffer);
  if(calc_gc) gt_vector_delete(bins);
  return GT_IFP_OK;
}

gt_sequence_archive *gt_open_sequence_archive(void) {
  gt_sequence_archive *sequence_archive = NULL;
  gt_input_file *const reference_file = 	gt_input_file_open(param.name_reference_file, false);
  sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
  if (gt_input_multifasta(reference_file, sequence_archive) != GT_IFP_OK) {
    gt_fatal_error_msg("Error parsing reference file '%s'\n", param.name_reference_file);
  }
  gt_input_file_close(reference_file);
  return sequence_archive;
}

void print_vcf_header(sr_param *param) {
  FILE *fp = param->output_file;
  fputs("##fileformat=VCFv4.2x\n", fp);
  time_t cl = time(0);
  struct tm *tt = localtime(&cl);
  fprintf(fp, "##fileDate(dd/mm/yyyy)=%02d/%02d/%04d\n", tt->tm_mday,
          tt->tm_mon + 1, tt->tm_year + 1900);
  fprintf(fp, "##source=bs_call_v%s,under_conversion=%g,over_conversion=%g,mapq_thresh=%d,bq_thresh=%d\n", BS_CALL_VERSION, param->under_conv, param->over_conv, param->mapq_thresh, param->min_qual);
  if(param->dbSNP_header != NULL) {
    fprintf(fp, "##dbsnp=<%s>\n", param->dbSNP_header);
  }
  GT_VECTOR_ITERATE(param->sam_headers->sequence_dictionary, header_record_p,
                    line_num, gt_sam_header_record *) {
    gt_sam_header_record *header_record = *header_record_p;
    gt_string *ctg = gt_shash_get(header_record, "SN", gt_string);
    gt_segmented_sequence* seg_seq = gt_sequence_archive_get_segmented_sequence(param->sequence_archive, gt_string_get_string(ctg));
    if (seg_seq == NULL) continue;
    gt_string *len = gt_shash_get(header_record, "LN", gt_string);
    if (ctg && len) {
      gt_string *sp = gt_shash_get(header_record, "SP", gt_string);
      if (param->species_filter) {
        if (!sp) continue;
        uint64_t l = gt_string_get_length(sp);
        if (strncasecmp(gt_string_get_string(sp), param->species_filter, l))
          continue;
        l = gt_string_get_length(ctg);
        char *ctg1 = gt_malloc(l + 1);
        memcpy(ctg1, gt_string_get_string(ctg), l);
        ctg1[l + 1] = 0;
        ctg_hash *tp = NULL;
        HASH_FIND_STR(param->species_hash, ctg1, tp);
        if (!tp) {
          tp = gt_alloc(ctg_hash);
          tp->flag = true;
          tp->ctg = ctg1;
          HASH_ADD_KEYPTR(hh, param->species_hash, ctg1, l, tp);
        }
      }
      fprintf(fp, "##contig=<ID=" PRIgts ",length=" PRIgts, PRIgts_content(ctg),
              PRIgts_content(len));
      gt_string *str = gt_shash_get(header_record, "AS", gt_string);
      if (str)
        fprintf(fp, ",assembly=" PRIgts, PRIgts_content(str));
      str = gt_shash_get(header_record, "M5", gt_string);
      if (str)
        fprintf(fp, ",md5=" PRIgts, PRIgts_content(str));
      if (sp)
        fprintf(fp, ",sp=\"" PRIgts "\"", PRIgts_content(sp));
      fputs(">\n", fp);
    }
  }
  fputs("##INFO=<ID=CX,Number=1,Type=String,Description=\"5 base sequence context (from position -2 to +2 on the positive strand) determined from the reference\">\n", fp);
  fputs("##FILTER=<ID=fail,Description=\"No sample passed filters\">\n", fp);
  fputs("##FILTER=<ID=q20,Description=\"Genotype Quality below 20\">\n", fp);
  fputs("##FILTER=<ID=qd2,Description=\"Quality By Depth below 2\">\n", fp);
  fputs("##FILTER=<ID=fs60,Description=\"Fisher Strand above 60\">\n", fp);
  fputs("##FILTER=<ID=mq40,Description=\"RMS Mapping Quality below 40\">\n", fp);
  fputs("##FILTER=<ID=gof20,Description=\"Goodness of fit above 20\">\n", fp);
  fputs("##FILTER=<ID=mac1,Description=\"Minor allele count <= 1\">\n", fp);
  fputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", fp);
  fputs("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample Genotype Filter\">\n", fp);	
  fputs("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype Likelihood\">\n", fp);
  fputs("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Phred scaled conditioanl genotype quality\">\n", fp);
  fputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth (non converted reads only)\">\n", fp);
  fputs("##FORMAT=<ID=MQ,Number=1,Type=Integer,Description=\"RMS Mapping Quality\">\n", fp);
  //  fputs("##FORMAT=<ID=AQ,Number=1,Type=Integer,Description=\"Average base quality\">\n", fp);
  fputs("##FORMAT=<ID=QD,Number=1,Type=Integer,Description=\"Quality By Depth (Variant quality / read depth (non-converted reads only))\">\n", fp);
  fputs("##FORMAT=<ID=MC8,Number=8,Type=Integer,Description=\"Base counts: non-informative for methylation (ACGT) followed by informative for methylation (ACGT)\">\n", fp);
  fputs("##FORMAT=<ID=AMQ,Number=.,Type=Integer,Description=\"Average base quailty for where MC8 base count non-zero\">\n", fp);
  fputs("##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Strand of Cytosine relative to reference sequence (+/-/+-/NA)\">\n", fp);
  fputs("##FORMAT=<ID=CG,Number=1,Type=String,Description=\"CpG Status (from genotype calls: Y/N/H/?)\">\n", fp);
  fputs("##FORMAT=<ID=CX,Number=1,Type=String,Description=\"5 base sequence context (from position -2 to +2 on the positive strand) determined from genotype call\">\n", fp);
  fputs("##FORMAT=<ID=FS,Number=1,Type=Integer,Description=\"Phred scaled log p-value from Fishers exact test of strand bias\"\n", fp);
  fputs("##FORMAT=<ID=GOF,Number=1,Type=Integer,Description=\"Phred scaled goodness of fit LR test for best genotype call against best call with free allele frequency parameter\n", fp);
  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", fp);
  if (param->sample_name)
    fprintf(fp, "\t%s\n", param->sample_name);
  else
    fputs("\tSAMPLE\n", fp);
}

gt_status bs_call_process(sr_param *param) {
  gt_status err = GT_STATUS_OK;
  param->output_file = stdout;
  /*
    if(param->pileup) {
    bool add_suffix = true;
    char *p = strrchr(param->pileup_file_name, '.');
    if(p && !strcmp(p, ".gz")) add_suffix = false;
    char *pcom;
    if(add_suffix) asprintf(&pcom, "bgzip > %s.gz",param->pileup_file_name);
    else asprintf(&pcom, "bgzip > %s",param->pileup_file_name);
    param->pileup_file = popen(pcom, "w");
    free(pcom);
    }
  */
  gt_input_file *input_file =
    param->input_file
    ? gt_input_file_sam_open(param->input_file, param->mmap_input)
    : gt_input_stream_sam_open(stdin);

  gt_sam_headers *sam_headers = gt_sam_header_new();
  uint64_t characters_read = 0, lines_read = 0;
  gt_status error_code = gt_input_file_sam_read_headers(
							(char *)input_file->file_buffer, input_file->buffer_size, sam_headers,
							&characters_read, &lines_read);
  if (error_code)
    gt_error(PARSE_SAM_HEADER_NOT_SAM, input_file->file_name);
  param->sam_headers = sam_headers;
  print_vcf_header(param);
  ref = gt_string_new(16384);

  gt_buffered_input_file *buffered_input =
    gt_buffered_input_file_new(input_file);
  gt_sam_parser_attributes *input_sam_attributes =
    gt_input_sam_parser_attributes_new();

  fill_base_prob_table();
  process_end = print_end = false;
  pthread_t process_thr, print_thr;
  pthread_create(&print_thr, NULL, print_thread, param);
  pthread_create(&process_thr, NULL, process_thread, param);
	
  gt_vector *templates = gt_vector_new(32, sizeof(gt_template *));
  error_code = input_sam_parser_get_template_vector(buffered_input, templates, input_sam_attributes, param);

  process_end = true;
  pthread_join(process_thr, 0);
  print_end = true;
  pthread_join(print_thr, 0);
  gt_input_sam_parser_attributes_delete(input_sam_attributes);
  gt_buffered_input_file_close(buffered_input);
  gt_input_file_close(input_file);
  //	if(param->pileup) fclose(param->pileup_file);
  return err;
}

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
  
static void *read_dbSNP_file(void *arg) {
  bool ok = true;
  sr_param *par = arg;
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
  if(sbuf->length > 6 && !strncmp(sbuf->buffer, "track ", 6)) {
    par->dbSNP_header = malloc(sbuf->length - 5);
    memcpy(par->dbSNP_header, sbuf->buffer + 6, sbuf->length - 6);
    par->dbSNP_header[sbuf->length - 6] = 0;
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
	  if(par->n_dbSNP_prefixes == 0xffff) {
	    fprintf(stderr, "Error in dbSNP file: too many prefixes\n");
	    exit(-1);
	  }
	  if(par->n_dbSNP_prefixes == n_p_store) {
	    n_p_store *= 1.5;
	    p_store = realloc(p_store, sizeof(void *) * n_p_store);
	  }
	  char *tp = p_store[par->n_dbSNP_prefixes++] = malloc(sbuf->length + 1);
	  memcpy(tp, sbuf->buffer, sbuf->length);
	  tp[sbuf->length] = 0;
	}
      } else break;
    }
    if(!par->n_dbSNP_prefixes) {
      fprintf(stderr, "Error in dbSNP file: no prefix information\n");
      ok = false;
    } else {
      par->dbSNP_prefix = malloc(sizeof(void *) * par->n_dbSNP_prefixes);
      memcpy(par->dbSNP_prefix, p_store, sizeof(void *) * par->n_dbSNP_prefixes);
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
	HASH_FIND(hh, par->dbSNP, ctg_name->buffer, sbuf->length, ctg);
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
	HASH_ADD_KEYPTR(hh, par->dbSNP, ctg->name, sbuf->length, ctg);
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
	  if(prefix_ix > par->n_dbSNP_prefixes) {
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
  } else ok = false;
	
  gt_string_delete(sbuf);
  gt_string_delete(ctg_name);
  gt_input_file_close(file);
  free(entries);
  free(name_buf);
  if(ok) fprintf(stderr,"Completed loading dbSNP (no. contigs %d, no. bins %d, no. SNPs %d\n", n_ctgs, n_bins, n_snps);
  else fprintf(stderr,"Error loading dbSNP\n");
  return NULL;
}

static int sort_cov_stats(gt_cov_stats *a, gt_cov_stats *b) {
  if(a->coverage < b->coverage) return -1;
  else if(a->coverage > b->coverage) return 1;
  return 0;
}

static void output_stats(sr_param *par) {
  static char *mut_type[] = {"A>C", "A>G", "A>T", "C>A", "C>G", "C>T", "G>A", "G>C", "G>T", "T>A", "T>C", "T>G"};		
  static char *filter_names[] = {
    "Passed", "Unmapped", "QC_Flags", "SecondaryAlignment", "MateUnmapped", "Duplicate", "NoPosition", "NoMatePosition", "MismatchContig", "BadOrientation", "LargeInsertSize", "NoSequence", "LowMAPQ", "NotCorrectlyAligned", "PairNotFound"
  };
  static char *base_filters[] = {
    "Passed", "Trimmed", "Clipped", "Overlapping", "LowQuality"
  };
  FILE *fp = par->json_file;
  bs_stats *stats = par->stats;
  fprintf(fp, "{\n\t\"source\": \"bs_call_v2.1, under_conversion=%g, over_conversion=%g, mapq_thresh=%d, bq_thread=%d\",\n",
          par->under_conv, par->over_conv, par->mapq_thresh, par->min_qual);
  time_t cl = time(0);
  struct tm *tt = localtime(&cl);
  fprintf(fp, "\t\"date\": \"%02d/%02d/%04d\",\n", tt->tm_mday, tt->tm_mon + 1, tt->tm_year + 1900);
  fputs("\t\"filterStats\": {\n\t\t\"ReadLevel\": {\n", fp);
  fprintf(fp, "\t\t\t\"%s\": {\n\t\t\t\t\"Reads\": %" PRIu64 ",\n\t\t\t\t\"Bases\": %" PRIu64 "\n\t\t\t}", filter_names[0], stats->filter_cts[0], stats->filter_bases[0]);
  for(int i = 1; i < 15; i++) {
    if(stats->filter_cts[i] > 0) {
      fprintf(fp, ",\n\t\t\t\"%s\": {\n\t\t\t\t\"Reads\": %" PRIu64 ",\n\t\t\t\t\"Bases\": %" PRIu64 "\n\t\t\t}", filter_names[i], stats->filter_cts[i], stats->filter_bases[i]);
    }
  }
  fputs("\n\t\t},\n\t\t\"BaseLevel\": {\n", fp);
  fprintf(fp, "\t\t\t\"%s\": %" PRIu64, base_filters[0], stats->base_filter[0]);
  for(int i = 1; i < 5; i++) {
    if(stats->base_filter[i] > 0) {
      fprintf(fp, ",\n\t\t\t\"%s\": %" PRIu64, base_filters[i], stats->base_filter[i]);
    }
  }
  fputs("\n\t\t}\n\t},\n\t\"totalStats\": {\n", fp);
  fprintf(fp, "\t\t\"SNPS\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->snps[stats_all], stats->snps[stats_passed]);
  fprintf(fp, "\t\t\"Indels\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->indels[stats_all], stats->indels[stats_passed]);
  fprintf(fp, "\t\t\"Multiallelic\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->multi[stats_all], stats->multi[stats_passed]);
  if(par->dbSNP != NULL) {
    fprintf(fp, "\t\t\"dbSNPSites\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->dbSNP_sites[stats_all], stats->dbSNP_sites[stats_passed]);
    fprintf(fp, "\t\t\"dbSNPVariantSites\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->dbSNP_var[stats_all], stats->dbSNP_var[stats_passed]);
  }
  fprintf(fp, "\t\t\"RefCpG\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->CpG_ref[stats_all], stats->CpG_ref[stats_passed]);
  fprintf(fp, "\t\t\"NonRefCpG\": {\n\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\"Passed\": %" PRIu64 "\n\t\t},\n", stats->CpG_nonref[stats_all], stats->CpG_nonref[stats_passed]);
  fputs("\t\t\"QCDistributions\": {\n", fp);
  fputs("\t\t\t\"FisherStrand\": ", fp);
  char term = '{';
  for(int i = 0; i < stats->fs_stats->used; i++) {
    fstats_cts *c = gt_vector_get_elm(stats->fs_stats, i, fstats_cts);
    if(c->cts[1] > 0) {
      fprintf(fp,"%c\n\t\t\t\t\"%d\": %" PRIu64, term, i, c->cts[1]);
      term = ',';
    }
  }
  if(term == '{') fputc(term, fp);
  fputs("\n\t\t\t},\n", fp);
  fputs("\t\t\t\"QualityByDepth\": ", fp);
  term = '{';
  for(int i = 0; i < stats->qd_stats->used; i++) {
    fstats_cts *c = gt_vector_get_elm(stats->qd_stats, i, fstats_cts);
    if(c->cts[0] + c->cts[1] > 0) {
      fprintf(fp,"%c\n\t\t\t\t\"%d\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", term, i, c->cts[0], c->cts[1]);
      term = ',';
    }
  }
  if(term == '{') fputc(term, fp);
  fputs("\n\t\t\t},\n", fp);
  fputs("\t\t\t\"RMSMappingQuality\": ", fp);
  term = '{';
  for(int i = 0; i < stats->mq_stats->used; i++) {
    fstats_cts *c = gt_vector_get_elm(stats->mq_stats, i, fstats_cts);
    if(c->cts[0] + c->cts[1] > 0) {
      fprintf(fp,"%c\n\t\t\t\t\"%d\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", term, i, c->cts[0], c->cts[1]);
      term = ',';
    }
  }
  if(term == '{') fputc(term, fp);
  fputs("\n\t\t\t},\n", fp);
  fputs("\t\t\t\"GoodnessOfFit\": ", fp);
  term = '{';
  for(int i = 0; i < stats->gof_stats->used; i++) {
    fstats_cts *c = gt_vector_get_elm(stats->gof_stats, i, fstats_cts);
    if(c->cts[0] + c->cts[1] > 0) {
      fprintf(fp,"%c\n\t\t\t\t\"%d\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", term, i, c->cts[0], c->cts[1]);
      term = ',';
    }
  }
  if(term == '{') fputc(term, fp);
  fputs("\n\t\t\t}\n\t\t},\t\t\"VCFFilterStats\": {\n", fp);
  fprintf(fp,"\t\t\t\"PASS\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", stats->filter_counts[0][0], stats->filter_counts[1][0]);
  for(int i = 1; i < 32; i++) {
    fputs(",\n\t\t\t", fp);
    int k = i;
    int f_ix = 0;
    char tmp = '\"';
    while(k) {
      if(k & 1) {
	fprintf(fp, "%c%s", tmp, flt_name[f_ix]);
	tmp = ',';
      }
      k >>= 1;
      f_ix++;
    }
    fprintf(fp,"\": {\"NonVariant\": %" PRIu64 ", \"Variant\": %" PRIu64 "}", stats->filter_counts[0][i], stats->filter_counts[1][i]);
  }
  fputs("\n\t\t},\n", fp);
  HASH_SORT(stats->cov_stats, sort_cov_stats);
  fputs("\t\t\"coverage\": {\n",fp);
  fputs("\t\t\t\"All\": ", fp);
  int ix = 0;
  term = '{';
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->all != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->all);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"Variant\": ", fp);
  term = '{';
  ix = 0;
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->var != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->var);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"RefCpG\": ", fp);
  term = '{';
  ix = 0;
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->CpG[0] != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG[0]);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"RefCpGInf\": ", fp);
  term = '{';
  ix = 0;
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->CpG_inf[0] != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG_inf[0]);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"NonRefCpG\": ", fp);
  term = '{';
  ix = 0;
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->CpG[1] != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG[1]);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"NonRefCpGInf\": ", fp);
  term = '{';
  ix = 0;
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(gcov->CpG_inf[1] != 0) {
      if(!(ix++)) {
	fprintf(fp, "%c\n\t\t\t\t", term);
	term = ',';
      } else fputs(", ", fp);
      fprintf(fp,"\"%" PRIu64 "\": %" PRIu64, gcov->coverage, gcov->CpG_inf[1]);
      ix %= 12;
    }
  }
  fputs("\n\t\t\t},\n\t\t\t\"GC\": ", fp);
  term = '{';
  for(gt_cov_stats *gcov = stats->cov_stats; gcov != NULL; gcov = gcov->hh.next) {
    if(!gcov->all) continue;
    fprintf(fp, "%c\n\t\t\t\t\"%" PRIu64 "\": [\n\t\t\t\t\t", term, gcov->coverage);
    term = ',';
    for(int i = 0; i < 100; i++) {
      fprintf(fp, "%" PRIu64 ",", gcov->gc_pcent[i]);
      if((i & 15) == 15) fputs("\n\t\t\t\t\t", fp);
      else fputc(' ', fp);
    }
    fprintf(fp, "%" PRIu64 "\n\t\t\t\t]", gcov->gc_pcent[100]);
  }
  fputs("\n\t\t\t}\n\t\t},\n\t\t\"quality\": {\n", fp);
  fputs("\t\t\t\"All\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 255; i++) {
    fprintf(fp, "%" PRIu64 ", ", stats->qual[all_sites][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
  }
  fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[all_sites][255]);
  fputs("\t\t\t\"Variant\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 255; i++) {
    fprintf(fp, "%" PRIu64 ",", stats->qual[variant_sites][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
    else fputc(' ', fp);
  }
  fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[variant_sites][255]);
  fputs("\t\t\t\"RefCpG\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 255; i++) {
    fprintf(fp, "%" PRIu64 ",", stats->qual[CpG_ref_sites][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
    else fputc(' ', fp);
  }
  fprintf(fp, "%" PRIu64 "\n\t\t\t],\n", stats->qual[CpG_ref_sites][255]);
  fputs("\t\t\t\"NonRefCpG\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 255; i++) {
    fprintf(fp, "%" PRIu64 ",", stats->qual[CpG_nonref_sites][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
    else fputc(' ', fp);
  }
  fprintf(fp, "%" PRIu64 "\n\t\t\t]\n", stats->qual[CpG_nonref_sites][255]);
  fputs("\t\t},\n\t\t\"mutations\": {\n", fp);
  for(int mut = 0; mut < 11; mut++) {
    fprintf(fp, "\t\t\t\"%s\": { \"All\": %" PRIu64 ", \"Passed\": %" PRIu64 ", \"dbSNPAll\": %" PRIu64 ", \"dbSNPPassed\": %" PRIu64 " },\n",
	    mut_type[mut], stats->mut_counts[mut][stats_all], stats->mut_counts[mut][stats_passed], 
	    stats->dbSNP_mut_counts[mut][stats_all], stats->dbSNP_mut_counts[mut][stats_passed]);
  }
  fprintf(fp, "\t\t\t\"%s\": { \"All\": %" PRIu64 ", \"Passed\": %" PRIu64 ", \"dbSNPAll\": %" PRIu64 ", \"dbSNPPassed\": %" PRIu64 " }\n",
	  mut_type[11], stats->mut_counts[11][stats_all], stats->mut_counts[11][stats_passed], 
	  stats->dbSNP_mut_counts[11][stats_all], stats->dbSNP_mut_counts[11][stats_passed]);
  fputs("\t\t},\n\t\t\"methylation\": {\n", fp);
  fputs("\t\t\t\"AllRefCpg\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 100; i++) {
    fprintf(fp, "%.8g, ", stats->CpG_ref_meth[stats_all][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
  }
  fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_ref_meth[stats_all][100]);
  fputs("\t\t\t\"PassedRefCpg\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 100; i++) {
    fprintf(fp, "%.8g, ", stats->CpG_ref_meth[stats_passed][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
  }
  fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_ref_meth[stats_passed][100]);
  fputs("\t\t\t\"AllNonRefCpg\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 100; i++) {
    fprintf(fp, "%.8g, ", stats->CpG_nonref_meth[stats_all][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
  }
  fprintf(fp, "%.8g\n\t\t\t],\n", stats->CpG_nonref_meth[stats_all][100]);
  fputs("\t\t\t\"PassedNonRefCpg\": [\n\t\t\t\t", fp);
  for(int i = 0; i < 100; i++) {
    fprintf(fp, "%.8g, ", stats->CpG_nonref_meth[stats_passed][i]);
    if((i & 15) == 15) fputs("\n\t\t\t\t", fp);
  }
  fprintf(fp, "%.8g\n\t\t\t]", stats->CpG_nonref_meth[stats_passed][100]);
  int nr = gt_vector_get_used(stats->meth_profile);
  if(nr) {
    fputs(",\n\t\t\t\"NonCpGreadProfile\": ", fp);
    term ='[';
    for(int i = 0; i < nr; i++) {
      meth_cts *mc = gt_vector_get_elm(stats->meth_profile, i, meth_cts);
      fprintf(fp, "%c\n\t\t\t\t[ %" PRIu64 ", %" PRIu64 ", %" PRIu64 ", %" PRIu64 " ]", term, mc->conv_cts[0], mc->conv_cts[1], mc->conv_cts[2], mc->conv_cts[3]);
      term = ',';
    }
    fputs("\n\t\t\t]", fp);
  }
  fputs("\n\t\t}\n\t},\n\t\"contigStats\": ", fp);
  term = '{';
  for(gt_ctg_stats *gs = par->ctg_stats; gs != NULL; gs = gs->hh.next) {
    if(gs->snps[stats_all] == 0) continue;
    fprintf(fp, "%c\n\t\t\"%s\": {\n", term, gs->contig);
    term = ',';
    fprintf(fp, "\t\t\t\"SNPS\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->snps[stats_all], gs->snps[stats_passed]);
    fprintf(fp, "\t\t\t\"Indels\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->indels[stats_all], gs->indels[stats_passed]);
    fprintf(fp, "\t\t\t\"Multiallelic\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->multi[stats_all], gs->multi[stats_passed]);
    if(par->dbSNP != NULL) {
      fprintf(fp, "\t\t\t\"dbSNPSites\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->dbSNP_sites[stats_all], gs->dbSNP_sites[stats_passed]);
      fprintf(fp, "\t\t\t\"dbSNPVariantSites\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->dbSNP_var[stats_all], gs->dbSNP_var[stats_passed]);
    }
    fprintf(fp, "\t\t\t\"RefCpG\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t},\n", gs->CpG_ref[stats_all], gs->CpG_ref[stats_passed]);
    fprintf(fp, "\t\t\t\"NonRefCpG\": {\n\t\t\t\t\"All\": %" PRIu64 ",\n\t\t\t\t\"Passed\": %" PRIu64 "\n\t\t\t}\n\t\t}", gs->CpG_nonref[stats_all], gs->CpG_nonref[stats_passed]);
  }
  fputs("\n\t}\n}\n", fp);
}

static void init_stats(sr_param *par) {
  par->json_file = fopen(par->report_file, "w");
  if(par->json_file == NULL) fprintf(stderr,"Could not open report file '%s' for output: %s\n", par->report_file, strerror(errno));
  else par->stats = calloc((size_t)1, sizeof(bs_stats));
  for(int i = 0; i < 100; i++) logp[i] = log(0.01 * (double)(i + 1));
  par->stats->meth_profile = gt_vector_new(256, sizeof(meth_cts));
  par->stats->qd_stats = gt_vector_new(256, sizeof(fstats_cts));
  par->stats->fs_stats = gt_vector_new(256, sizeof(fstats_cts));
  par->stats->mq_stats = gt_vector_new(256, sizeof(fstats_cts));
  par->stats->gof_stats = gt_vector_new(256, sizeof(fstats_cts));
  memset(par->stats->meth_profile->memory, 0, sizeof(meth_cts) * par->stats->meth_profile->elements_allocated);
  memset(par->stats->qd_stats->memory, 0, sizeof(fstats_cts) * par->stats->qd_stats->elements_allocated);
  memset(par->stats->fs_stats->memory, 0, sizeof(fstats_cts) * par->stats->fs_stats->elements_allocated);
  memset(par->stats->mq_stats->memory, 0, sizeof(fstats_cts) * par->stats->mq_stats->elements_allocated);
  memset(par->stats->gof_stats->memory, 0, sizeof(fstats_cts) * par->stats->gof_stats->elements_allocated);
}

gt_status parse_arguments(int argc, char **argv) {
  gt_status err = GT_STATUS_OK;
  struct option *bs_call_getopt = gt_options_adaptor_getopt(bs_call_options);
  gt_string *const bs_call_short_getopt =  gt_options_adaptor_getopt_short(bs_call_options);
  int option, option_index;
  bool load_seq = false;
  param.num_threads = sysconf(_SC_NPROCESSORS_ONLN);
  while (true) {
    // Get option &  Select case
    if ((option =
	 getopt_long(argc, argv, gt_string_get_string(bs_call_short_getopt),
		     bs_call_getopt, &option_index)) == -1)
      break;
    switch (option) {
      /* Operations */
    case 'N':
      param.no_split = true;
      break;
    case '1':
      param.haploid = true;
      break;
    case 'd':
      param.keep_duplicates = true;
      break;
    case 'k':
      param.keep_unmatched = true;
      break;
    case 's':
      param.extra_stats = true;
      break;
      //		 case 'P':
      //			if(param.pileup_file_name != NULL) free(param.pileup_file_name);
      //			param.pileup_file_name = strdup(optarg);
      //	 	  param.pileup = true;
      //			break;
    case 'R':
      param.right_trim = atol(optarg);
      break;
    case 'B':
      param.blank_trim = true;
      break;
    case 'L':
      param.left_trim = atol(optarg);
      break;
    case 'q':
      param.mapq_thresh = (uint8_t)atoi(optarg);
      break;
    case 'Q':
      param.min_qual = (uint8_t)atoi(optarg);
      break;
    case 'l':
      param.max_template_len = atol(optarg);
    case 'T':
      param.realign_tol = atol(optarg);
      break;
      /* IO */
    case 'o':
      param.output_prefix = optarg;
      break;
    case 'n':
      param.sample_name = optarg;
      break;
    case 'x':
      param.species_filter = optarg;
      break;
    case 'r':
      param.name_reference_file = optarg;
      load_seq = 1;
      break;
      //    case 'I':
      //      param.name_gem_index_file = optarg;
      //      load_seq = 1;
      //      break;
    case 'D':
      param.dbSNP_name = optarg;
      break;
    case 'A':
      param.all_positions = true;
      break;
    case 'z':
#ifdef HAVE_ZLIB
      param.compress = GZIP;
#endif
      break;
    case 'j':
#ifdef HAVE_BZLIB
      param.compress = BZIP2;
#endif
      break;
    case 'Z':
      param.compress = NONE;
      break;
    case 201:
      param.mmap_input = true;
      break;
    case 202:
      param.report_file = optarg;
      break;
      /* Model */
    case 'c':
      if (sscanf(optarg, "%lf,%lf", &param.under_conv, &param.over_conv) != 2)
        gt_fatal_error_msg(
			   "c (conversion) option requires two comma separated arguments");
      break;
      //    case 301:
      //      param.caller = graphical_model;
      //      break;
      //    case 302:
      //      param.caller = maximum_likelihood;
      //      break;
    case 303:
      param.ref_bias = atof(optarg);
      break;
      /* Misc */
    case 'v':
      param.verbose = true;
      break;
    case 't':
      param.num_threads = atol(optarg);
      break;
    case 'h':
    case '?':
      usage(bs_call_options, bs_call_groups, false);
      exit(1);
      break;
    case 'H':
      usage(bs_call_options, bs_call_groups, true);
      exit(1);
      break;
    default:
      usage(bs_call_options, bs_call_groups, false);
      gt_fatal_error_msg("Option '%c' %d not recognized", option, option);
      break;
    }
  }
  if (optind < argc) param.input_file = argv[optind];
  pthread_t dbsnp_thread;
  if(param.dbSNP_name != NULL) {
    pthread_create(&dbsnp_thread, NULL, read_dbSNP_file, &param);
  }
  // Sanity
  if (param.min_qual < 1)
    param.min_qual = 1;
  else if (param.min_qual > MAX_QUAL)
    param.min_qual = MAX_QUAL;
  if (param.num_threads < 1)
    param.num_threads = 1;
  if (load_seq) {
    fprintf(stderr, "Loading reference sequences\n");
    param.sequence_archive = gt_open_sequence_archive();
    fprintf(stderr, "Completed loading reference sequences\n");
  } else {
    fputs("Error in bs_calls(): a sequence archive is mandatory\n",
          stderr);
    err = GT_STATUS_FAIL;
  }
  // Free
  gt_string_delete(bs_call_short_getopt);
  gt_free(bs_call_getopt);
  // JSON Stats
  if(param.report_file != NULL) init_stats(&param);
  if(param.dbSNP_name != NULL && err == GT_STATUS_OK) pthread_join(dbsnp_thread, 0);
  return err;
}

int main(int argc, char *argv[]) {
  gt_status err = 0;
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  err = parse_arguments(argc, argv);
  if (err == GT_STATUS_OK) {
    lfact_store_init();
    err = bs_call_process(&param);
  }
  if (param.sequence_archive) gt_sequence_archive_delete(param.sequence_archive);
  if(param.json_file != NULL) {
    output_stats(&param);
    fclose(param.json_file);
  }
  return err == GT_STATUS_OK ? 0 : -1;
}
