/*
* bs_call.c
*
*  Created on: 30 Sep 2014
*      Author: heath
*/

#define BS_CALL "bs_call"
#define BS_CALL_VERSION "2.0"

#include <stdio.h>
#include <getopt.h>
#include <ctype.h>
#include <pthread.h>
#include <gsl/gsl_integration.h>
#include "gem_tools.h"
#include "gt_pipe_io.h"

#include "bs_call.h"

#define LOG10 (2.30258509299404568402)
#define LOG2 (0.69314718055994530942)

void fexact_(int *, int *, double *, int *, double *, double *, double *,
             double *, double *);

#define LFACT_STORE_SIZE 256

static double lfact_store[LFACT_STORE_SIZE];

typedef struct {
  char *ctg;
  bool flag;
  UT_hash_handle hh;
} ctg_hash;

typedef struct {
  char *input_file;
  char *name_reference_file;
  char *name_gem_index_file;
  char *output_prefix;
  char *sample_name;
  char *species_filter;
  ctg_hash *species_hash;
  bool mmap_input;
  /* Control flags */
  bool is_paired;
  bool no_split;
  bool extra_stats;
  bool keep_duplicates;
  bool haploid;
  bool verbose;
  bool blank_trim;
	bool pileup;
  BS_Caller caller;
  int left_trim;
  int right_trim;
	char *pileup_file_name;
  FILE *output_file;
  FILE *pileup_file;
  gt_sam_headers *sam_headers;
  uint8_t mapq_thresh;
  uint8_t min_qual;
  uint64_t max_template_len;
  uint64_t realign_tol;
  double under_conv, over_conv;
  gt_output_file_compression compress;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output;
  int num_threads;
  gt_sequence_archive *sequence_archive;
} sr_param;

sr_param param = {
    .input_file = NULL,
    .name_reference_file = NULL,
    .name_gem_index_file = NULL,
    .output_prefix = NULL,
    .sample_name = NULL,
    .species_filter = NULL,
    .species_hash = NULL,
    .mmap_input = false,
    .compress = NONE,
    .verbose = false,
    .haploid = false,
    .blank_trim = false,
	  .pileup = false,
    .caller = maximum_likelihood,
    .left_trim = 0,
    .right_trim = 0,
	  .pileup_file_name = NULL,
    .output_file = NULL,
	  .pileup_file = NULL,
    .sam_headers = NULL,
    .mapq_thresh = DEFAULT_MAPQ_THRESH,
    .min_qual = MIN_QUAL,
    .max_template_len = DEFAULT_MAX_TEMPLATE_LEN,
    .realign_tol = DEFAULT_REALIGN_TOL,
    .under_conv = DEFAULT_UNDER_CONVERSION,
    .over_conv = DEFAULT_OVER_CONVERSION,
    .is_paired = false,
    .no_split = false,
    .extra_stats = false,
    .keep_duplicates = false,
    .num_threads = 1,
    .sequence_archive = NULL,
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

double fisher(double *tab) {
  int c[4], row[2], col[2];
  for (int i = 0; i < 4; i++)
    c[i] = (int)tab[i];
  row[0] = c[0] + c[1];
  row[1] = c[2] + c[3];
  col[0] = c[0] + c[2];
  col[1] = c[1] + c[3];
  int n = row[0] + row[1];
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

typedef struct {
  double match, match1, no_match;
  double a1, a2, e, b, z;
} qual_prob;

typedef struct {
  uint64_t counts[8];
  double gt_prob[10];
  double gt_gof; /* Goodness of fit LR */
  uint8_t max_gt;
  uint8_t base_qual;
  uint8_t mapq_qual;
} gt_meth;

typedef struct _base_counts {
  struct _base_counts *next;
  int8_t idx[MAX_QUAL + 1];
  uint32_t counts[MAX_QUAL + 1];
} base_counts;

typedef struct {
  base_counts *base[8];
  uint32_t mask;
} locus_counts;

typedef struct {
  uint32_t counts[2][8];
  uint32_t n;
  float quality;
  float mapq;
	uint8_t *seq;
	size_t seq_size;
	size_t seq_idx;
} pileup;

static gsl_integration_glfixed_table *gauss_tab[MAX_GAUSS_N];
static double gauss_y_tab[8][MAX_GAUSS_N];
double *gauss_x[MAX_GAUSS_N];

// Optimized 8 level Quality binning scheme from Illumina
static char qual_code[256] = {
	 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2,
	 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5,
	 5, 5, 5, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
	 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};
static int qual_bin[8] = {0, 6, 15, 22, 27, 33, 37, 40};

static double *gauss_y[2][10] = {
    {0, gauss_y_tab[0], 0, 0, gauss_y_tab[1], gauss_y_tab[2], gauss_y_tab[3], 0,
     0, 0},
    {0, 0, gauss_y_tab[4], 0, 0, gauss_y_tab[5], 0, gauss_y_tab[6],
     gauss_y_tab[7], 0},
};
static qual_prob q_prob[MAX_QUAL + 1];

static base_counts *base_count_free_list;

static base_counts *new_base_counts(void) {
  base_counts *bc = base_count_free_list;
  if (!bc) {
    bc = gt_alloc(base_counts);
  } else {
    base_count_free_list = bc->next;
  }
  bc->next = NULL;
  memset(bc->idx + 1, 0, MAX_QUAL + 1);
  bc->idx[0] = -1;
  return bc;
}

static void fill_base_prob_table(double under, double over) {
  const double ln10 = LOG10;
  double lambda = 1.0 - under;
  double tau = over;

  for (int q = 0; q <= MAX_QUAL; q++) {
    double e = exp(-.1 * (double)q * ln10);
    double no_match = log(e / 3.0);
    double match1 = log(.5 - e / 3.0);
    double match = log(1.0 - e);
    double z = 1.0 - 4.0 * e / 3.0;
    double a1 = 1.0 - e - lambda * z;
    double a2 = e / 3.0 + lambda * z;
    double b = (lambda - tau) * z;
    q_prob[q].match = match;
    q_prob[q].match1 = match1;
    q_prob[q].no_match = no_match;
    q_prob[q].a1 = a1;
    q_prob[q].a2 = a2;
    q_prob[q].e = e;
    q_prob[q].b = b;
    q_prob[q].z = z;
  }
}

static void fill_x(gsl_integration_glfixed_table *tab, double *x) {
  int n = tab->n;
  int m = (n + 1) >> 1;
  double *tx = tab->x;
  int i = 0;
  int j = 0;
  if (n & 1) {
    x[j++] = .5;
    i++;
  }
  for (; i < m; i++) {
    x[j++] = (1.0 + tx[i]) * .5;
    x[j++] = (1.0 - tx[i]) * .5;
  }
}

static gsl_integration_glfixed_table *get_gl_tab(int n, double *xx) {
  if (n) {
    int nt = (n + 4) >> 1;
    if (nt > MAX_GAUSS_N)
      nt = MAX_GAUSS_N;
    if (!gauss_tab[nt - 1]) {
      gauss_tab[nt - 1] = gsl_integration_glfixed_table_alloc(nt);
      gauss_x[nt - 1] = gt_malloc(sizeof(double) * nt);
      fill_x(gauss_tab[nt - 1], gauss_x[nt - 1]);
    }
    memcpy(xx, gauss_x[nt - 1], sizeof(double) * nt);
    return gauss_tab[nt - 1];
  } else
    return 0;
}

void do_calling(FILE *fp, locus_counts *curr_loc, gt_meth *gtmeth,
                sr_param *param) {

  static int vtab[] = {0, 1, 1, 1, 0, 1, 1, 0, 1, 0}; // Heterozygote genotypes

  double gt_prob[10] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  uint64_t ncounts[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  uint64_t nn = 0;
  for (int i = 0; i < 8; i++) {
    base_counts *bc = curr_loc->base[i];
    if (bc) {
      for (int j = bc->idx[0]; j > 0; j = bc->idx[j]) {
        uint32_t n = bc->counts[j];
        ncounts[i] += n;
        nn += n;
      }
    }
  }
  if (!nn) {
    for (int i = 0; i < 8; i++)
      gtmeth->counts[i] = 0;
    for (int i = 0; i < 10; i++)
      gtmeth->gt_prob[i] = 0.0;
    return;
  }
  uint64_t n_c2t = ncounts[5] + ncounts[7]; // C and T on C2T strand
  uint64_t n_g2a = ncounts[4] + ncounts[6]; // A and G on G2A strand
  gsl_integration_glfixed_table *tab1, *tab2;
  double x1[MAX_GAUSS_N], x2[MAX_GAUSS_N];
  if (n_c2t) {
    tab1 = get_gl_tab(n_c2t, x1);
    for (int j = 0; j < 4; j++)
      for (size_t i = 0; i < tab1->n; i++)
        gauss_y_tab[j][i] = 0.0;
  } else
    tab1 = 0;
  if (n_g2a) {
    tab2 = get_gl_tab(n_g2a, x2);
    for (int j = 4; j < 8; j++)
      for (size_t i = 0; i < tab2->n; i++)
        gauss_y_tab[j][i] = 0.0;
  } else
    tab2 = 0;
  // A - non-informative
  base_counts *bc = curr_loc->base[0];
  if (bc) {
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[0] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[1] += zz;
      gt_prob[2] += zz;
      gt_prob[3] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[4] += zz;
      gt_prob[5] += zz;
      gt_prob[6] += zz;
      gt_prob[7] += zz;
      gt_prob[8] += zz;
      gt_prob[9] += zz;
    }
  }
  // C - non-informative
  bc = curr_loc->base[1];
  if (bc) {
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[4] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[1] += zz;
      gt_prob[5] += zz;
      gt_prob[6] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[2] += zz;
      gt_prob[3] += zz;
      gt_prob[7] += zz;
      gt_prob[8] += zz;
      gt_prob[9] += zz;
    }
  }
  // G - non-informative
  bc = curr_loc->base[2];
  if (bc) {
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[7] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[2] += zz;
      gt_prob[5] += zz;
      gt_prob[8] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[1] += zz;
      gt_prob[3] += zz;
      gt_prob[4] += zz;
      gt_prob[6] += zz;
      gt_prob[9] += zz;
    }
  }
  // T - non-informative
  bc = curr_loc->base[3];
  if (bc) {
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[9] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[3] += zz;
      gt_prob[6] += zz;
      gt_prob[8] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[1] += zz;
      gt_prob[2] += zz;
      gt_prob[4] += zz;
      gt_prob[5] += zz;
      gt_prob[7] += zz;
    }
  }
  // A - informative
  bc = curr_loc->base[4];
  if (bc) {
    int n = tab2->n;
    double *restrict y2 = gauss_y[1][2];
    double *restrict y5 = gauss_y[1][5];
    double *restrict y7 = gauss_y[1][7];
    double *restrict y8 = gauss_y[1][8];
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[0] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[1] += zz;
      gt_prob[3] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[4] += zz;
      gt_prob[6] += zz;
      gt_prob[9] += zz;
      gt_prob[2] -= LOG2;
      gt_prob[5] -= LOG2;
      gt_prob[8] -= LOG2;
      double a2 = q_prob[q].a2;
      double b = q_prob[q].b;
      double e = q_prob[q].e;
      double ax = a2 + 1.0 - e;
      double ay = a2 + e / 3.0;
      for (int j = 0; j < n; j++) {
        zz = b * x2[j];
        y2[j] += zn * log(ax - zz);
        y7[j] += zn * log(a2 - zz);
        zz = zn * log(ay - zz);
        y5[j] += zz;
        y8[j] += zz;
      }
    }
  }
  // C - informative
  bc = curr_loc->base[5];
  if (bc) {
    int n = tab1->n;
    double *restrict y1 = gauss_y[0][1];
    double *restrict y4 = gauss_y[0][4];
    double *restrict y5 = gauss_y[0][5];
    double *restrict y6 = gauss_y[0][6];
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      double zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[2] += zz;
      gt_prob[3] += zz;
      gt_prob[7] += zz;
      gt_prob[8] += zz;
      gt_prob[9] += zz;
      gt_prob[1] -= LOG2;
      gt_prob[5] -= LOG2;
      gt_prob[6] -= LOG2;
      double a1 = q_prob[q].a1;
      double b = q_prob[q].b;
      double e = q_prob[q].e;
      double ax = a1 + e / 3.0;
      for (int j = 0; j < n; j++) {
        zz = b * x1[j];
        y4[j] += zn * log(a1 + zz);
        zz = zn * log(ax + zz);
        y1[j] += zz;
        y5[j] += zz;
        y6[j] += zz;
      }
    }
  }
  // G - informative
  bc = curr_loc->base[6];
  if (bc) {
    int n = tab2->n;
    double *restrict y2 = gauss_y[1][2];
    double *restrict y5 = gauss_y[1][5];
    double *restrict y7 = gauss_y[1][7];
    double *restrict y8 = gauss_y[1][8];
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      double zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[1] += zz;
      gt_prob[3] += zz;
      gt_prob[4] += zz;
      gt_prob[6] += zz;
      gt_prob[9] += zz;
      gt_prob[2] -= LOG2;
      gt_prob[5] -= LOG2;
      gt_prob[8] -= LOG2;
      double a1 = q_prob[q].a1;
      double b = q_prob[q].b;
      double e = q_prob[q].e;
      double ax = a1 + e / 3.0;
      for (int j = 0; j < n; j++) {
        zz = b * x2[j];
        y7[j] += zn * log(a1 + zz);
        zz = zn * log(ax + zz);
        y2[j] += zz;
        y5[j] += zz;
        y8[j] += zz;
      }
    }
  }
  // T - informative
  bc = curr_loc->base[7];
  if (bc) {
    int n = tab1->n;
    double *restrict y1 = gauss_y[0][1];
    double *restrict y4 = gauss_y[0][4];
    double *restrict y5 = gauss_y[0][5];
    double *restrict y6 = gauss_y[0][6];
    for (int q = bc->idx[0]; q > 0; q = bc->idx[q]) {
      double zn = (double)bc->counts[q];
      gt_prob[9] += zn * q_prob[q].match;
      double zz = zn * q_prob[q].match1;
      gt_prob[3] += zz;
      gt_prob[8] += zz;
      zz = zn * q_prob[q].no_match;
      gt_prob[0] += zz;
      gt_prob[2] += zz;
      gt_prob[7] += zz;
      gt_prob[1] -= LOG2;
      gt_prob[5] -= LOG2;
      gt_prob[6] -= LOG2;
      double a2 = q_prob[q].a2;
      double b = q_prob[q].b;
      double e = q_prob[q].e;
      double ax = a2 + 1.0 - e;
      double ay = a2 + e / 3.0;
      for (int j = 0; j < n; j++) {
        zz = b * x1[j];
        y6[j] += zn * log(ax - zz);
        y4[j] += zn * log(a2 - zz);
        zz = zn * log(ay - zz);
        y1[j] += zz;
        y5[j] += zz;
      }
    }
  }
  for (int i = 0; i < 10; i++) {
    if (param->haploid && vtab[i])
      continue;
    double *y = gauss_y[0][i];
    if (y && tab1) {
      int n = tab1->n;
      int m = (n + 1) >> 1;
      double ymax = y[0];
      for (int j = 1; j < n; j++)
        if (y[j] > ymax)
          ymax = y[j];
      for (int j = 0; j < n; j++)
        y[j] = exp(y[j] - ymax);
      int k, k1;
      double s;
      double *w = tab1->w;
      if (n & 1) {
        s = w[0] * y[0];
        k = k1 = 1;
      } else {
        s = 0.0;
        k = k1 = 0;
      }
      for (; k < m; k++, k1 += 2)
        s += w[k] * (y[k1] + y[k1 + 1]);
      gt_prob[i] += ymax + log(s * .5);
    }
    y = gauss_y[1][i];
    if (y && tab2) {
      int n = tab2->n;
      int m = (n + 1) >> 1;
      double ymax = y[0];
      for (int j = 1; j < n; j++)
        if (y[j] > ymax)
          ymax = y[j];
      for (int j = 0; j < n; j++)
        y[j] = exp(y[j] - ymax);
      int k, k1;
      double s;
      double *w = tab2->w;
      if (n & 1) {
        s = w[0] * y[0];
        k = k1 = 1;
      } else {
        s = 0.0;
        k = k1 = 0;
      }
      for (; k < m; k++, k1 += 2)
        s += w[k] * (y[k1] + y[k1 + 1]);
      gt_prob[i] += ymax + log(s * .5);
    }
  }
  double zm = gt_prob[0];
  double sum = 0.0;
  int best_gt = 0;
  if (param->haploid) {
    for (int i = 1; i < 10; i++) {
      if (vtab[i])
        gt_prob[i] = -DBL_MAX;
      else if (gt_prob[i] > zm) {
        zm = gt_prob[i];
        best_gt = i;
      }
    }
    for (int i = 0; i < 10; i++) {
      if (!vtab[i]) {
        gt_prob[i] = exp(gt_prob[i] - zm);
        sum += gt_prob[i];
      }
      double *restrict gp = gtmeth->gt_prob;
      for (int i = 0; i < 10; i++)
        if (!vtab[i])
          gp[i] = log(gt_prob[i] / sum) / LOG10;
    }
  } else {
    for (int i = 1; i < 10; i++)
      if (gt_prob[i] > zm) {
        zm = gt_prob[i];
        best_gt = i;
      }
    for (int i = 0; i < 10; i++) {
      gt_prob[i] = exp(gt_prob[i] - zm);
      sum += gt_prob[i];
    }
    double *restrict gp = gtmeth->gt_prob;
    for (int i = 0; i < 10; i++)
      gp[i] = log(gt_prob[i] / sum) / LOG10;
  }
  memcpy(gtmeth->counts, ncounts, sizeof(uint64_t) * 8);
}

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
  uint64_t *counts = gtm->counts;
  uint64_t dp = 0;
  for (int i = 0; i < 8; i++)
    dp += counts[i];
  if (!dp)
    return;
  const int rfc = toupper((int)rf_ctxt[2]);
  int rfix = base_tab[rfc];
  int gt = gt_store[2] - 1;
  // Skip homozygous reference if AA or TT
  if (gt_flag[gt][rfix])
    return;
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
  fprintf(fp, "%s\t%" PRIu64 "\t.\t%c\t%s\t%d", ctg, x, rfc, ref_alt[gt][rfix], phred);
  // FILTER field
  if (phred < 20)
    fputs("\tq20", fp);
  else
    fputs("\tPASS", fp);
  // INFO field
  fprintf(fp, "\tCX=%.5s", rf_ctxt);
  // FORMAT field
  fputs("\tGT:DP:GL:MC8:CS:CG:CX", fp);
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
  fprintf(fp, "\t%s:%" PRIu64, gt_str[gt][rfix], dp);
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
  fprintf(fp, ":%s:%s:%.5s\n", cs_str[gt], cpg, ctxt);
}

static char gt_store[5];
uint64_t store_x = 0;
static gt_meth gtm_store[5];
static char *curr_ctg;
static char rf_ctxt[7];

// Print the last 2 entries in gt_store
void flush_vcf_entries(FILE *fp) {
  if (curr_ctg && store_x) {
    for (int i = 0; i < 2; i++) {
      memcpy(gt_store, gt_store + 1, 4);
      memcpy(gtm_store, gtm_store + 1, 4 * sizeof(gt_meth));
      memcpy(rf_ctxt, rf_ctxt + 1, 6);
      if (gt_store[2])
        _print_vcf_entry(fp, curr_ctg, gtm_store + 2, rf_ctxt, store_x - 1 + i,
                         gt_store);
    }
    gt_free(curr_ctg);
    curr_ctg = 0;
    store_x = 0;
  }
}

void print_vcf_entry(FILE *fp, char *ctg, gt_meth *gtm, const char *rf,
                     const uint64_t x, const uint64_t xstart) {
  if (curr_ctg) {
    if (strcmp(curr_ctg, ctg)) {
      flush_vcf_entries(fp);
    }
  }

  if (!curr_ctg) {
    uint64_t l = strlen(ctg);
    curr_ctg = gt_malloc(l + 1);
    memcpy(curr_ctg, ctg, l);
    curr_ctg[l] = 0;
  }
  uint64_t l = x - store_x;
  if (l < 5) {
    memcpy(gt_store, gt_store + l, 5 - l);
    memcpy(gtm_store, gtm_store + l, (5 - l) * sizeof(gt_meth));
    for (uint64_t i = 4; i >= 5 - l; i--) {
      gt_store[i] = 0;
    }
  } else
    memset(gt_store, 0, 5);
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
  double z = gtm->gt_prob[0];
  int gt = 0;
  for (int i = 1; i < 10; i++)
    if (gtm->gt_prob[i] > z) {
      z = gtm->gt_prob[i];
      gt = i;
    }
  gt_store[4] = gt + 1;
  if (gt_store[2])
    _print_vcf_entry(fp, curr_ctg, gtm_store + 2, rf_ctxt, x - 2, gt_store);
}

inline static void calc_meth(double n1, double n2, double n3, double n4,
                             double e, double l, double t, double w, double *p,
                             double *m, double *pz) {
  double q = 3.0 - 4.0 * e;
  if (n3 + n4 == 0) { // No converted reads
    double u = n1 / (n1 + n2);
    *p = (u * (w * q + 2.0 * e) - e) / q;
    //		printf("u=%g, p=%g\n",u,*p);
    if (*p < 0.0)
      *p = 0.0;
    else if (*p > w)
      *p = w;
  } else {
    double r = e / q;
    if (n1 + n2 == 0) { // No non-converted reads
      double Q = (n3 + n4) / (n3 * (w * q + e) - n4 * e);
      *pz = 1.0 / (q * Q);
      if (*pz > w)
        *pz = w;
      else if (*pz < 0) {
        *pz = 0;
      }
    } else {
      // Both converted and non-converted reads
      double s = w + e / q;
      // First check whether the gradient at c=0 is 0 or negative.  If so, set c
      // to the boundary
      if ((n1 + n3) / r - (n2 + n4) / s <= 0.0) {
        *p = *pz = 0.0;
      } else {
        double Q = (n3 + n4) / (n3 * (w * q + e) - n4 * e);
        double K = (n3 > 0.0 ? n3 / (1.0 + e * Q) : 0.0) +
                   (n4 > 0.0 ? n4 / (1.0 + (q * w + 3.0 * e - 3.0) * Q) : 0.0);
        double A = n1 + n2 + K;
        double B = (n2 + K) * r - (n1 + K) * s;
        //        printf("Q=%g, K=%.10g, A=%.10g, B=%.10g\n", Q, K, A, B);

        *p = (-B + sqrt(B * B + 4.0 * A * K * r * s)) / (2.0 * A);

        if (*p > w)
          *p = w;
        double z = 1.0 / (Q * q * (*p));
        *m = (z - 1.0 + l) / (l - t);
        //        printf("z=%g, p=%g, m=%g\n", z, *p, *m);
        if (*m > 1.0 || *m < 0.0) {
          if (*m > 1.0) {
            *m = 1.0;
            z = 1.0 - t;
            *p = w * (n1 + n3) / (n1 + n2 + n3 + n4);
          } else {
            *m = 0.0;
            z = 1.0 - l;
            *p = w * n1 / (n1 + n2);
          }
          double K1, K2, K3, K4;
          // Evaluate gradient at w to see if still positive
          double gd =
              n1 / (w + r) + n2 / (w - s) + n3 / (w + r / z) + n4 / (w - s / z);
          //          printf("gd = %.10g\n", gd);
          if (gd >= 0.0) {
            *p = w;
          } else {
            for (int it = 0; it < 100; it++) {
              K1 = (*p + r);
              K2 = (*p - s);
              K3 = (*p + r / z);
              K4 = (*p - s / z);
              gd = n1 / K1 + n2 / K2 + n3 / K3 + n4 / K4;
              double h = -(n1 / (K1 * K1) + n2 / (K2 * K2) + n3 / (K3 * K3) +
                           n4 / (K4 * K4));
              //              printf("%g %g %g %g\n", *p, gd, h, gd / h);

              *p -= gd / h;
              if (fabs(gd) < 1.0e-6)
                break;
            }
          }
        }
        *pz = *p * z;
      }
    }
  }
}

#define TEST_M(z, p, l, t)                                                     \
  {                                                                            \
    double m = ((z) / (p)-1.0 + (l)) / ((l) - (t));                            \
    if ((m) < 0.0)                                                             \
      (z) = (p) * (1.0 - (l));                                                 \
    else if ((m) > 1.0)                                                        \
      (z) = (p) * (1.0 - (t));                                                 \
  }

static inline void get_zp(double x1, double x2, double e, double l, double t,
                          double *z) {
  double q = 3.0 - 4.0 * e;
  // w = 1/2, p = 1/2
  double Q = (x1 + x2) / (x1 * (0.5 * q + e) - x2 * e);
  z[0] = 1.0 / (Q * q);
  TEST_M(z[0], 0.5, l, t);
  // w = 1, p = 1/2
  Q = (x1 + x2) / (x1 * (q + e) - x2 * e);
  z[1] = z[2] = 1.0 / (Q * q);
  TEST_M(z[1], 0.5, l, t);
  // w = 1, p = 1
  TEST_M(z[2], 1.0, l, t);
}

void calc_gt_prob(gt_meth *gt, pileup *tp, sr_param *param, char rf) {
  double e = exp((float)tp->quality * -0.1 * LOG10);
  double l = 1.0 - param->under_conv;
  double t = param->over_conv;
  double a = -1.0, c = -1.0, g = -1.0, cz = -1.0, gz = -1.0;
  double m1 = -1.0, m2 = -1.0;
  double n[8];
  for (int i = 0; i < 8; i++)
    n[i] = (double)gt->counts[i];
  // Estimate frequency of A+G
  double x1 = n[1] + n[3] + n[5] + n[7];
  double x2 = n[0] + n[2] + n[4] + n[6];
  double u = x2 / (x1 + x2);
  double w = (3.0 * u - 2.0 * e) / (3.0 - 4.0 * e);
  if (w > 1.0) {
    w = 1.0;
    c = 0.0;
  } else {
    if (w < 0.0)
      w = 0.0;
    // Find c & m1
    calc_meth(n[1], n[3], n[5], n[7], e, l, t, 1.0 - w, &c, &m1, &cz);
  }
  if (w > 0.0) { // Find g & m2
    calc_meth(n[2], n[0], n[6], n[4], e, l, t, w, &g, &m2, &gz);
  } else {
    g = a = 0.0;
  }
  if (g >= 0.0)
    a = w - g;
  double q = (1.0 - 4.0 * e / 3.0);
  double ll[11];
  // Add in prior from reference
  // We just slightly bias the reference homozygote as twice as likely as the
  // others
  {
    int i;
    switch (rf) {
    case 'A':
      ll[0] = ll[1] = 0.0;
      for (i = 2; i < 11; i++)
        ll[i] = -LOG2;
      break;
    case 'C':
      for (i = 1; i < 11; i++)
        ll[i] = -LOG2;
      ll[0] = ll[5] = 0.0;
      break;
    case 'G':
      for (i = 1; i < 11; i++)
        ll[i] = -LOG2;
      ll[0] = ll[8] = 0.0;
      break;
    case 'T':
      for (i = 1; i < 11; i++)
        ll[i] = -LOG2;
      ll[0] = ll[10] = 0.0;
      break;
    default:
      for (i = 0; i < 11; i++)
        ll[i] = 0.0;
      break;
    }
  }
  double k1 = log(e / 3.0);
  double k2 = log(1.0 - e);
  double k3 = log(0.5 - e / 3.0);
  if (n[0]) {
    ll[0] += log(a * q + e / 3.0) * n[0]; // ML
    ll[1] += k2 * n[0];                   // AA
    double tz = k3 * n[0];
    ll[2] += tz; // AC
    ll[3] += tz; // AG
    ll[4] += tz; // AT
    tz = k1 * n[0];
    ll[5] += tz;  // CC
    ll[6] += tz;  // CG
    ll[7] += tz;  // CT
    ll[8] += tz;  // GG
    ll[9] += tz;  // GT
    ll[10] += tz; // TT
  }
  if (n[1]) {
    ll[0] += log(c * q + e / 3.0) * n[1];
    ll[5] += k2 * n[1]; // CC
    double tz = k3 * n[1];
    ll[2] += tz; // AC
    ll[6] += tz; // CG
    ll[7] += tz; // CT
    tz = k1 * n[1];
    ll[1] += tz;  // AA
    ll[3] += tz;  // AG
    ll[4] += tz;  // AT
    ll[8] += tz;  // GG
    ll[9] += tz;  // GT
    ll[10] += tz; // TT
  }
  if (n[2]) {
    ll[0] += log(g * q + e / 3.0) * n[2];
    ll[8] += k2 * n[2]; // GG
    double tz = k3 * n[2];
    ll[3] += tz; // AG
    ll[6] += tz; // CG
    ll[9] += tz; // TG
    tz = k1 * n[2];
    ll[1] += tz;  // AA
    ll[2] += tz;  // AC
    ll[4] += tz;  // AT
    ll[5] += tz;  // CC
    ll[7] += tz;  // CT
    ll[10] += tz; // TT
  }
  if (n[3]) {
    ll[0] += log((1.0 - w - c) * q + e / 3.0) * n[3];
    ll[10] += k2 * n[3]; // TT
    double tz = k3 * n[3];
    ll[4] += tz; // AT
    ll[7] += tz; // CT
    ll[9] += tz; // GT
    tz = k1 * n[3];
    ll[1] += tz; // AA
    ll[2] += tz; // AC
    ll[3] += tz; // AG
    ll[5] += tz; // CC
    ll[6] += tz; // CG
    ll[8] += tz; // GG
  }
  double zp[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  if (n[5] + n[7] > 0.0)
    get_zp(n[5], n[7], e, l, t, zp);
  if (n[4] + n[6] > 0.0)
    get_zp(n[6], n[4], e, l, t, zp + 3);
  for (int i = 0; i < 6; i++) {
    double m = -1.0;
    if (zp[i] >= 0.0) {
      if (zp[i] >= 0.0) {
        m = (zp[i] / (i % 3 == 2 ? 1.0 : 0.5) - 1.0 + l) / (l - t);
      }
    }
  }
  if (n[4]) {
    ll[0] += log(1.0 - q * (1.0 - w + gz) - e) * n[4];
    ll[1] += k2 * n[4];                       // AA
    ll[3] += log(1.0 - q * zp[4] - e) * n[4]; // AG
    ll[8] += log(1.0 - q * zp[5] - e) * n[4]; // GG
    double tz = log(1.0 - q * (0.5 + zp[3]) - e) * n[4];
    ll[6] += tz; // CG
    ll[9] += tz; // GT
    tz = k3 * n[4];
    ll[2] += tz; // AC
    ll[4] += tz; // AT
    tz = k1 * n[4];
    ll[5] += tz;  // CC
    ll[7] += tz;  // CT
    ll[10] += tz; // TT
                  //		printf("AA: %c %g %g\n",rf,ll[0],ll[8]);
  }
  if (n[5]) {
    ll[0] += log(cz * q + e / 3.0) * n[5];
    ll[5] += log(zp[2] * q + e / 3.0) * n[5]; // CC
    double tz = log(zp[0] * q + e / 3.0) * n[5];
    ll[2] += tz;                              // AC
    ll[6] += tz;                              // CG
    ll[7] += log(zp[1] * q + e / 3.0) * n[5]; // CT
    tz = k1 * n[5];
    ll[1] += tz;  // AA
    ll[3] += tz;  // AG
    ll[4] += tz;  // AT
    ll[8] += tz;  // GG
    ll[9] += tz;  // GT
    ll[10] += tz; // TT
  }
  if (n[6]) {
    ll[0] += log(gz * q + e / 3.0) * n[6];
    ll[8] += log(zp[5] * q + e / 3.0) * n[6]; // GG
    double tz = log(zp[3] * q + e / 3.0) * n[6];
    ll[6] += tz; // CG
    ll[9] += tz; // TG
    ll[3] += log(zp[4] * q + e / 3.0) * n[6];
    tz = k1 * n[6];
    ll[1] += tz;  // AA
    ll[2] += tz;  // AC
    ll[4] += tz;  // AT
    ll[5] += tz;  // CC
    ll[7] += tz;  // CT
    ll[10] += tz; // TT
  }
  if (n[7]) {
    ll[0] += log(1.0 - q * (w + cz) - e) * n[7];
    ll[10] += k2 * n[7];                      // TT
    ll[5] += log(1.0 - q * zp[2] - e) * n[7]; // CC
    ll[7] += log(1.0 - q * zp[1] - e) * n[7]; // CT
    double tz = log(1.0 - q * (0.5 + zp[0]) - e) * n[7];
    ll[2] += tz; // AC
    ll[6] += tz; // CG
    tz = k3 * n[7];
    ll[4] += tz; // AT
    ll[9] += tz; // GT
    tz = k1 * n[7];
    ll[1] += tz; // AA
    ll[3] += tz; // AG
    ll[8] += tz; // GG
  }
  double max = ll[1];
  int mx = 1;
  for (int i = 2; i < 11; i++) {
    if (ll[i] > max) {
      max = ll[i];
      mx = i;
    }
  }
  gt->max_gt = mx - 1;
  double sum = 0.0;
  for (int i = 0; i < 10; i++) {
    sum += exp(ll[i + 1] - max);
  }
  sum = log(sum);
  for (int i = 0; i < 10; i++) {
    gt->gt_prob[i] = (ll[i + 1] - max - sum) / LOG10;
  }
  gt->gt_gof = ll[0] - max;
}

static char stop_base[256] = { ['A'] = 'a', ['C'] = 'c', ['G'] = 'g', ['T'] = 't', ['.'] = ',' };

void call_genotypes_ML(char *ctg, gt_vector *align_list, uint64_t x, uint64_t y,
                       gt_string *ref, sr_param *param) {
  char *gt_string[] = {"AA", "AC", "AG", "AT", "CC",
                       "CG", "CT", "GG", "GT", "TT"};

  FILE *fp = param->output_file;
	FILE *fp_pile = param->pileup ? param->pileup_file : NULL;
  uint64_t nr = gt_vector_get_used(align_list);
//  fprintf(stderr,"Getting pileup %" PRIu64 " - %" PRIu64 " no. read pairs: %" PRIu64 "\n", x, y, nr);
  assert(y >= x);
  uint64_t sz = y - x + 1;
  pileup *counts = gt_calloc(sz, pileup, true);
  gt_meth *gt_res = gt_calloc(sz, gt_meth, true);
  const char *ref_st = gt_string_get_string(ref);
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
		uint64_t old_pos = 0;
    for (int k = 0; k < 2; k++) {
			if(!al->read[k]) continue;
      uint64_t rl = gt_string_get_length(al->read[k]);
			if(rl == 0) continue;
			// fprintf(stderr,"%d\t%lu\t%lu\n",k,k?al->reverse_position:al->forward_position,rl);
      float mapq = al->mapq[k];
      if (!rl)
        continue;
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
			for(j = rl - 1; j >= 0; j--) {
        int c = base_tab_st[st][(int)sp[j]];
				if(c) break;
			}
			uint64_t read_end;
			if(j >= 0) read_end = j;
			else continue;
			pos += read_start;
			if(param->pileup) {
				// Fill in gap between read pair (if any) with '_'
				if(k == 1 && old_pos > 0 && pos > old_pos) {
//					fprintf(stderr,"%lu\t%lu\n",old_pos,pos);
					pileup *cl = counts + old_pos - x;
					for(;old_pos < pos; old_pos++, cl++) {
//						fprintf(stderr," : %lu\t%lu\n",old_pos,pos);
						if(cl->seq_size == 0) {
							cl->seq_size = 32;
							cl->seq = malloc(cl->seq_size);
						}
						if(cl->seq_idx >= cl->seq_size) {
							cl->seq_size <<= 1;
							cl->seq = realloc(cl->seq, cl->seq_size);
						}
						cl->seq[cl->seq_idx++] = 3 | (st << 2);
					}
				}
			}
      pileup *curr_loc = counts + pos - x;
      for (j = read_start; j <= read_end && pos <= y; j++, pos++, curr_loc++) {
        int c = base_tab_st[st][(int)sp[j]];
        int q = sq[j] - QUAL_CONV;
				if(param->pileup) {
					bool first = (j == read_start && k == 0) ? true : false;
					bool last = ((j == read_end && k == 1) || pos == y) ? true : false;
					 if(curr_loc->seq_size == 0) {
							curr_loc->seq_size = 32;
							curr_loc->seq = malloc(curr_loc->seq_size);
					 }
					size_t req = (first || last) ? 2 : 1;
					if(curr_loc->seq_idx + req > curr_loc->seq_size) {
						curr_loc->seq_size <<= 1;
						curr_loc->seq = realloc(curr_loc->seq, curr_loc->seq_size);
					}
					if(first) { // Start of read
						curr_loc->seq[curr_loc->seq_idx++] = 12 | (st << 4);
					}
					int qbin = qual_code[q];
					int c1 = base_tab[(int)sp[j]];
					uint8_t pbyte = (c1 > 0 && qbin > 0) ? (c1 - 1) | (st << 2) | (qbin << 4) : st << 2;
					curr_loc->seq[curr_loc->seq_idx++] = pbyte;
					if(last) { // End of read
						curr_loc->seq[curr_loc->seq_idx++] = 13 | (st << 4);
					}
				}
        if (c-- && q >= param->min_qual) {
          if (q > MAX_QUAL)
            q = MAX_QUAL;
          curr_loc->n++;
          curr_loc->quality += (float)q;
          curr_loc->mapq += (float)mapq;
          curr_loc->counts[ori][c]++;
        }
      }
			old_pos = pos;
      ori ^= 1;
    }
  }
  pileup *tp = counts;
  gt_meth *tg = gt_res;
//  printf("Calling\n");
  	char *base_str = "NACGT";
	size_t buf_size = 256;
	char *buf = NULL;
	if(param->pileup) {
		size_t buf_size = 256;
		buf = malloc(buf_size);
	}
  for (uint64_t i = x; i <= y; i++, tp++, tg++) {
    if (tp->n) {
      tp->quality = floor(0.5 + tp->quality / (float)tp->n);
      tp->mapq = floor(0.5 + tp->mapq / (float)tp->n);
      double table[16];
      double *tab = table;
      int nr = 2, nc = 0;
      double fisher_p = -1.0;
      for (int j = 0; j < 8; j++) {
        if (tp->counts[0][j] + tp->counts[1][j]) {
          *tab++ = tp->counts[0][j];
          *tab++ = tp->counts[1][j];
          nc++;
          tg->counts[j] = tp->counts[0][j] + tp->counts[1][j];
        }
      }
      calc_gt_prob(tg, tp, param, ref_st[i - x]);
      print_vcf_entry(fp, ctg, tg, ref_st, i, x);
			if(param->pileup) {
				int n[3] = {0, 0, 0};
				for(size_t j = 0; j < tp->seq_idx; j++) {
					int pbyte = tp->seq[j];
					int st = (pbyte >> 2) & 3;
					if(st < 3) n[st]++;
				}
				size_t sz1 = 2 * (n[0] + n[1] + n[2]) + 6;
				if(sz1 > buf_size) {
					buf_size = sz1;
					buf = realloc(buf, buf_size);
				}
				int gt = tg->max_gt;
				double z = exp(tg->gt_prob[gt] * LOG10);
				int phred;
				if (z >= 1.0)
					phred = 255;
				else {
					phred = (int)(-10.0 * log(1.0 - z) / LOG10);
					if (phred > 255)
						phred = 255;
				}
				int ix[6];
				ix[0] = 0;
				ix[1] = n[0] + 1;
				ix[2] = ix[1] + n[1] + 1;
				ix[3] = ix[2] + n[2] + 1;
				ix[4] = ix[3] + n[0] + 1;
				ix[5] = ix[4] + n[1] + 1;
				char *sbp[3] = {buf, buf + ix[1], buf + ix[2]};
				char *qbp[3] = {buf + ix[3], buf + ix[4], buf + ix[5]};
				char ref = ref_st[i - x];
				for(size_t j = 0; j < tp->seq_idx; j++) {
					int pbyte = tp->seq[j];
					int st = (pbyte >> 2) & 3;
					if(st == 3) {
						st = (pbyte >> 4) & 3;
						assert(st < 3);
						if(pbyte & 1) {
							int ixp = ix[st] - 1;
							assert(ixp >= 0);
							buf[ixp] = stop_base[buf[ixp]];
						}
					} else {
						int qual = qual_bin[(pbyte >> 4) & 7];
						char base;
						if(qual == 0) {
							base = ((pbyte & 3) == 0) ? 'N' : '_';
						} else {
							base = base_str[1 + (pbyte & 3)];
						}
						if(base == ref) base = '.';
						buf[ix[st]++] = base;
						buf[ix[st + 3]++] = qual + QUAL_CONV;
					}
				}
				for(int k = 0; k < 6; k++) buf[ix[k]] = 0;
				fprintf(fp_pile, "%s\t%" PRIu64 "\t%c\t%s\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", ctg, i, ref, gt_string[gt], phred, sbp[1], sbp[2], sbp[0], qbp[1], qbp[2], qbp[0]);
			}
#if 0						
      if (nc == 2) {
        fisher_p = fisher(table);
      } else if (nc > 2) {
        double tt = 0.0, zp;
        //        fexact_(&nr, &nc, table, &nr, &tt, &tt, &tt, &zp, &fisher_p);
      }
      int mx = tg->max_gt;
      int phred = (int)(0.5 - 10.0 * log(1.0 - tg->gt_prob[mx]) / LOG10);
      int gof = (int)(0.5 + 10.0 * tg->gt_gof / LOG10);
      printf("%lu\t%c\t%g\t%g\t%d\t\t%d\t%g\t%s\t", i, ref_st[i - x],
             tp->quality, tp->mapq, phred, gof, tg->gt_gof, gt_string[mx]);
      for (int j = 0; j < 2; j++) {
        printf("%u", tp->counts[j][0]);
        for (int k = 1; k < 8; k++)
          printf(",%u", tp->counts[j][k]);
        if (!j)
          fputc('\t', stdout);
      }
      if (fisher_p > 0.0) {
        int fish = (int)(0.5 - 10.0 * log(fisher_p) / LOG10);
        printf("\t%d", fish);
      } else if (fisher_p < 0.0)
        fputs("\t-", stdout);
      else
        fputs("\t-99999", stdout);
      fputc('\n', stdout);
#endif
    }
  }
  gt_free(gt_res);
	if(param->pileup) {
		for(size_t i = 0; i < sz; i++) {
			if(counts[i].seq_size) free(counts[i].seq);
		}
	}
  gt_free(counts);
	if(buf != NULL) free(buf);
}

void call_genotypes_GM(char *ctg, gt_vector *align_list, uint64_t x, uint64_t y,
                       gt_string *ref, sr_param *param) {

  FILE *fp = param->output_file;
  uint64_t nr = gt_vector_get_used(align_list);
  assert(y >= x);
  uint64_t sz = y - x + 1;
  locus_counts *counts = gt_calloc(sz, locus_counts, true);
  gt_meth gtm;
  const char *ref_st = gt_string_get_string(ref);
  align_details **al_p = gt_vector_get_mem(align_list, align_details *);
  uint64_t curr_x = x;
  for (uint64_t ix = 0; ix < nr; ix++, al_p++) {
    align_details *al = *al_p;
    uint64_t x1 = al->forward_position;
    if (al->reverse_position < x1)
      x1 = al->reverse_position;
    assert(x1 >= curr_x);
    for (; curr_x < x1; curr_x++) {
      locus_counts *curr_loc = counts + curr_x - x;
      if (curr_loc->mask) {
        do_calling(fp, curr_loc, &gtm, param);
        print_vcf_entry(fp, ctg, &gtm, ref_st, curr_x, x);
        for (int i = 0; i < 8; i++) {
          if (curr_loc->base[i]) {
            curr_loc->base[i]->next = base_count_free_list;
            base_count_free_list = curr_loc->base[i];
            curr_loc->base[i] = NULL;
          }
        }
        curr_loc->mask = 0;
      }
    }
    int st = al->bs_strand;
    for (int k = 0; k < 2; k++) {
			if(!al->read[k]) continue;
      uint64_t rl = gt_string_get_length(al->read[k]);
      if (rl == 0) continue;
      char *sp = gt_string_get_string(al->read[k]);
      char *sq = gt_string_get_string(al->qualities[k]);
      uint64_t pos = k ? al->reverse_position : al->forward_position;
      locus_counts *curr_loc = counts + pos - x;
      for (uint64_t j = 0; j < rl && pos <= y; j++, pos++, curr_loc++) {
        int c = base_tab_st[st][(int)sp[j]];
        int q = sq[j] - QUAL_CONV;
        if (c-- && q >= param->min_qual) {
          if (q > MAX_QUAL)
            q = MAX_QUAL;
          if (!(curr_loc->mask & (1 << c))) {
            curr_loc->base[c] = new_base_counts();
            curr_loc->mask |= (1 << c);
          }
          base_counts *bp = curr_loc->base[c];
          if (!bp->idx[q]) {
            bp->idx[q] = bp->idx[0];
            bp->idx[0] = q;
            bp->counts[q] = 1;
          } else
            bp->counts[q]++;
        }
      }
    }
  }
  for (; curr_x <= y; curr_x++) {
    locus_counts *curr_loc = counts + curr_x - x;
    if (curr_loc->mask) {
      do_calling(fp, curr_loc, &gtm, param);
      print_vcf_entry(fp, ctg, &gtm, ref_st, curr_x, x);
      for (int i = 0; i < 8; i++) {
        if (curr_loc->base[i]) {
          curr_loc->base[i]->next = base_count_free_list;
          base_count_free_list = curr_loc->base[i];
          curr_loc->base[i] = NULL;
        }
      }
      curr_loc->mask = 0;
    }
  }
  gt_free(counts);
}

void call_genotypes(char *ctg, gt_vector *align_list, uint64_t x, uint64_t y,
                    gt_string *ref, sr_param *param) {

  switch (param->caller) {
  case graphical_model:
    call_genotypes_GM(ctg, align_list, x, y, ref, param);
    break;
  case maximum_likelihood:
    call_genotypes_ML(ctg, align_list, x, y, ref, param);
    break;
  default:
    fprintf(stderr, "Unknown caller\n");
    exit(-1);
  }
}

typedef struct {
  gt_vector *align_list;
  char *ctg;
  uint64_t y;
} align_block;

gt_status process_template_vector(gt_vector *align_list, char *ctg, uint64_t y,
                                  sr_param *param) {
  align_details **al_p = gt_vector_get_mem(align_list, align_details *);
  int ix = gt_vector_get_used(align_list);
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
//  fprintf(stderr, "Processing %" PRIu64 " - %" PRIu64 " (%d pairs)\n", x, y, ix);
  for (int i = 0; i < ix; i++) {
    align_details *al =
        *(align_details **)gt_vector_get_elm(align_list, i, align_details *);
    if (!(al->read[0] || al->read[1]))
      continue;
    if (al->reverse_position > 0 && al->reverse_position < x)
      x = al->reverse_position;
    if (al->forward_position == x1 && al->reverse_position == x2 &&
        al->bs_strand == strand) {
//			fprintf(stderr,"AA: %" PRIu64 ", %" PRIu64 " i = %d\n", al->forward_position, al->reverse_position, i);
      curr_ct++;
			uint8_t mapq=0;
			int kn = 0;
			for(int ix1 = 0; ix1 < 2; ix1++) {
				if(al->read[ix1]) {
					mapq += al->mapq[ix1];
					kn++;
				}
			}
			mapq/=kn;
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
      if (curr_ct > 0 && mx_mapq >= thresh) {
        if (j < k) {
          align_details *al2 = *(align_details **)gt_vector_get_elm(
              align_list, k, align_details *);
          align_details *al1 = *(align_details **)gt_vector_get_elm(
              align_list, j, align_details *);
          gt_vector_set_elm(align_list, j, align_details *, al2);
          gt_vector_set_elm(align_list, k, align_details *, al1);
        }
        j++;
      }
      x1 = al->forward_position;
      x2 = al->reverse_position;
      strand = al->bs_strand;
      curr_ct = 1;
      mx_qual = get_al_qual(al);
			mx_mapq=0;
			int kn = 0;
			for(int ix1 = 0; ix1 < 2; ix1++) {
				if(al->read[ix1]) {
					mx_mapq += al->mapq[ix1];
					kn++;
				}
			}
      mx_mapq/= kn;
			k = i;
    }
  }
  if (curr_ct && mx_mapq >= thresh) {
    if (j < ix) {
      align_details *al2 =
          *(align_details **)gt_vector_get_elm(align_list, k, align_details *);
      align_details *al1 =
          *(align_details **)gt_vector_get_elm(align_list, j, align_details *);
      gt_vector_set_elm(align_list, j++, align_details *, al2);
      gt_vector_set_elm(align_list, k, align_details *, al1);
    }
  }
  gt_vector_set_used(align_list, j);
  ix = j;
  // We do trimming by setting to lower case the trimmed bases
  if (param->left_trim || param->right_trim) {
    for (int i = 0; i < ix; i++) {
      align_details *al =
          *(align_details **)gt_vector_get_elm(align_list, i, align_details *);
			if(al->read[0]) {
				char *sp = gt_string_get_string(al->read[0]);
				uint64_t rl = gt_string_get_length(al->read[0]);
				for (int k1 = 0; k1 < param->left_trim && k1 < rl; k1++)
					sp[k1] = tolower((int)sp[k1]);
				for (int k1 = 0; k1 < param->right_trim && k1 < rl; k1++)
					sp[rl - k1 - 1] = tolower((int)sp[rl - k1 - 1]);
			}
			if(al->read[1]) {
				char *sp = gt_string_get_string(al->read[1]);
				int rl = gt_string_get_length(al->read[1]);
				for (int k1 = 0; k1 < param->left_trim && k1 < rl; k1++)
					sp[rl - k1 - 1] = tolower((int)sp[rl - k1 - 1]);
				for (int k1 = 0; k1 < param->right_trim && k1 < rl; k1++)
					sp[k1] = tolower((int)sp[k1]);
			}
    }
  }
  assert(x > 0 && x <= y);
//  fprintf(stderr,"%lu -> %lu\n",x,y);
  x1 = (x > 2) ? x - 2 : 1;
  uint64_t y1 = y + 4;
  uint64_t length = y1 - x1 + 1;
  gt_string *ref = gt_string_new(length);
  gt_status gts = gt_sequence_archive_get_sequence_string(
      param->sequence_archive, ctg, FORWARD, x1 - 1, length, ref);
  if (gts) {
		if(gts == GT_SEQUENCE_NOT_FOUND) {
			fprintf(stderr,"Could not find reference sequence for contig '%s'\n", ctg);
			return GT_STATUS_OK;
		}
		fprintf(stderr,"Ref problem - retrying\n");
    x1 = x;
    y1 = y;
    length = y1 - x1 + 1;
    gt_status gts = gt_sequence_archive_get_sequence_string(
        param->sequence_archive, ctg, FORWARD, x1 - 1, length, ref);
    if (gts) {
			fprintf(stderr,"No good  - aborting\n");
      return GT_STATUS_OK;
		}
  }
  x = x1;
  y = y1;
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
            gt_string_trim_left((*al_p)->read[k], adj);
            gt_string_trim_left((*al_p)->qualities[k], adj);
          } else {
            if (misms->position + misms->size != rl)
              gt_fatal_error_msg(
                  "Error in CIGAR - Illegal soft clip (%d %" PRIu64 " %" PRIu64 " %" PRIu64 ")\n", z,
                  misms->position, misms->size, rl);
            gt_string_trim_right((*al_p)->read[k], misms->size);
            gt_string_trim_right((*al_p)->qualities[k], misms->size);
          }
        } else if (nclip) {
          misms->position -= adj;
          gt_vector_set_elm((*al_p)->mismatches[k], z - nclip, gt_misms,
                            *misms);
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
				overlap = (*al_p)->reference_span[0] - (*al_p)->reverse_position +
					(*al_p)->forward_position;
				rev = false;
			} else {
				overlap = (*al_p)->reference_span[1] + (*al_p)->reverse_position -
					(*al_p)->forward_position;
				rev = true;
			}
			// Look for overlapping reads - keep only best part
			if ((*al_p)->forward_position + (*al_p)->reference_span[0] >=
					(*al_p)->reverse_position) {
				// Find the overlapping parts
				uint64_t *rspan;
				// Find best quality read (normally read 1)
				int trim_read;
				rspan = (*al_p)->reference_span;
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
										gt_misms *misms = gt_vector_get_elm(
																												(*al_p)->mismatches[trim_read], z1, gt_misms);
										misms->position -= trim;
										gt_misms *misms1 = gt_vector_get_elm(
																												 (*al_p)->mismatches[trim_read], z1 - z, gt_misms);
										gt_vector_set_elm((*al_p)->mismatches[trim_read], z1 - z,
																			gt_misms, *misms);
										gt_vector_set_elm((*al_p)->mismatches[trim_read], z1,
																			gt_misms, *misms1);
									}
									num_misms -= z;
									gt_vector_set_used((*al_p)->mismatches[trim_read], num_misms);
									break;
								} else {
									for (int z1 = 0; z1 < num_misms; z1++) {
										gt_misms *misms = gt_vector_get_elm(
																												(*al_p)->mismatches[trim_read], z1, gt_misms);
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
										gt_misms *misms = gt_vector_get_elm(
																												(*al_p)->mismatches[trim_read], z1, gt_misms);
										misms->position -= trim;
										if (z2) {
											gt_misms *misms1 = gt_vector_get_elm(
																													 (*al_p)->mismatches[trim_read], z1 - z2, gt_misms);
											gt_vector_set_elm((*al_p)->mismatches[trim_read], z1 - z2,
																				gt_misms, *misms);
											gt_vector_set_elm((*al_p)->mismatches[trim_read], z1,
																				gt_misms, *misms1);
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
			}
    }
    rdl[0] = (*al_p)->read[0] ? gt_string_get_length((*al_p)->read[0]) : 0;
    rdl[1] = (*al_p)->read[1] ? gt_string_get_length((*al_p)->read[1]) : 0;
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
        uint64_t ix = 0, ix1 = 0;
        for (int z = 0; z < num_misms; z++) {
          gt_misms *misms =
              gt_vector_get_elm((*al_p)->mismatches[k], z, gt_misms);
          uint64_t s = misms->position - ix;
          memcpy(sp + ix1, tp + ix, s);
          memcpy(sq + ix1, tq + ix, s);

          if (misms->misms_type == INS) {
            ix = misms->position;
            ix1 += s;
            for (int k1 = 0; k1 < misms->size; k1++, ix1++) {
              sp[ix1] = 'N';
              sq[ix1] = 33;
            }
          } else if (misms->misms_type == DEL) {
            ix = misms->position + misms->size;
            ix1 += s;
          }
        }
        uint64_t s = rdl[k] - ix;
        memcpy(sp + ix1, tp + ix, s);
        memcpy(sq + ix1, tq + ix, s);
        ix1 += s;
        gt_string_set_length((*al_p)->read[k], ix1);
        gt_string_set_length((*al_p)->qualities[k], ix1);
        if (del_size)
          gt_free(tp);
      }
    }
  }
  if (ix) {
    call_genotypes(ctg, align_list, x, y, ref, param);
  }
  gt_string_delete(ref);
  return GT_STATUS_OK;
}

#if 0
// Expects a SAM file sorted on coordinates.  Selects mapped pairs with TQ and
// TQ >= thresh, and make up vector of overlapping templates
gt_status input_sam_parser_get_template_vector(
    gt_buffered_input_file *const buffered_sam_input, gt_vector *align_list,
    gt_sam_parser_attributes *const sam_parser_attr, sr_param *param) {

  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  gt_status error_code;
  char *curr_ctg = NULL;
  bool chr_skip = false;
  uint64_t max_pos = 0; // Position of righmost end of current pileup
  uint64_t read_idx = 0;
  gt_string *tag = gt_string_new(128);
  align_details *al = 0, *align_hash = 0;
  gt_vector *free_list = gt_vector_new(256, sizeof(align_details *));
  // Cycle through input lines until next read does not overlap with current
  // pileup
  do {
    // Check the end_of_block. Reload buffer if needed
    if (gt_buffered_input_file_eob(buffered_sam_input)) {
      if ((error_code = gt_input_sam_parser_reload_buffer(
               buffered_sam_input)) != GT_ISP_OK) {
        int ix = gt_vector_get_used(align_list);
        if (ix) {
          // Clean up
        }
        return error_code;
      }
    }
		fprintf(stderr,"And here\n");
    // Check file format
    if (gt_input_sam_parser_check_sam_file_format(buffered_sam_input)) {
      gt_error(PARSE_SAM_BAD_FILE_FORMAT,
               buffered_sam_input->input_file->file_name,
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
    if (!al) {
      if (gt_vector_get_used(free_list)) {
        al = *(gt_vector_get_last_elm(free_list, align_details *));
        gt_vector_dec_used(free_list);
        if (al->tag)
          free(al->tag);
        al->tag = 0;
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
      } else {
        al = gt_malloc(sizeof(align_details));
        al->seq_name = al->read[0] = al->read[1] = al->qualities[0] =
            al->qualities[1] = 0;
        al->tag = 0;
        al->mismatches[0] = gt_vector_new(8, sizeof(gt_misms));
        al->mismatches[1] = gt_vector_new(8, sizeof(gt_misms));
      }
    }
    bool reverse;
    error_code = gt_isp_quick_parse_bs_sam_alignment(
        text_line, al, param->mapq_thresh, param->max_template_len, &reverse);
    gt_input_sam_parser_next_record(buffered_sam_input);
    if (error_code) {
      if (error_code == GT_ISP_SAM_FILTERED) {
				fprintf(stderr,"FILTERED!\n");
        continue;
      }
      gt_input_sam_parser_prompt_error(buffered_sam_input, line_num,
                                       buffered_sam_input->cursor - line_start,
                                       error_code);
      return GT_ISP_FAIL;
    }
		fprintf(stderr,"Read OK!\n");
    gt_input_parse_tag_chomp_pairend_info(tag);
    if (al->tag)
      free(al->tag);
    al->tag = strdup(gt_string_get_string(tag));
    bool new_block = false;
    if (!curr_ctg || strncmp(curr_ctg, gt_string_get_string(al->seq_name),
                             gt_string_get_length(al->seq_name))) {
      chr_skip = false;
      if (curr_ctg)
        gt_free(curr_ctg);
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
      fprintf(stderr, "Processing chromosome %s (%s)\n", curr_ctg,
              chr_skip ? "SKIP" : "OK");
      new_block = true;
    }
    if (chr_skip)
      continue;
    bool insert;
		if(al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) {
			if (reverse)
				insert = al->forward_position > al->reverse_position ? true : false;
			else
				insert = al->forward_position < al->reverse_position ? true : false;
			if (al->forward_position == al->reverse_position) {
				insert = false;
				if (!new_block) {
					align_details *thash;
					HASH_FIND(hh, align_hash, al->tag, gt_string_get_length(tag), thash);
					if (!thash)
						insert = true;
				} else
					insert = true;
			}
		}
    if (new_block == false && insert == true) {
      if (al->forward_position > max_pos && al->reverse_position > max_pos) {
        new_block = true;
      }
    }
  } while(1);
  return GT_ISP_OK;
}

#else
// Expects a SAM file sorted on coordinates.  Selects mapped pairs with TQ and
// TQ >= thresh, and make up vector of overlapping templates
gt_status input_sam_parser_get_template_vector(
    gt_buffered_input_file *const buffered_sam_input, gt_vector *align_list,
    gt_sam_parser_attributes *const sam_parser_attr, sr_param *param) {
  GT_BUFFERED_INPUT_FILE_CHECK(buffered_sam_input);
  GT_VECTOR_CHECK(align_list);
  gt_vector_clear(align_list);
  gt_status error_code;
  char *curr_ctg = NULL;
  bool chr_skip = false;
  uint64_t max_pos = 0; // Position of righmost end of current pileup
  uint64_t read_idx = 0;
  gt_string *tag = gt_string_new(128);
  align_details *al = 0, *align_hash = 0;
  gt_vector *free_list = gt_vector_new(256, sizeof(align_details *));

  // Cycle through input lines until next read does not overlap with current
  // pileup
  do {
    // Check the end_of_block. Reload buffer if needed
    if (gt_buffered_input_file_eob(buffered_sam_input)) {
      if ((error_code = gt_input_sam_parser_reload_buffer(
               buffered_sam_input)) != GT_ISP_OK) {
        int ix = gt_vector_get_used(align_list);
        if (ix) {
          process_template_vector(align_list, curr_ctg, max_pos, param);
          flush_vcf_entries(param->output_file);
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
    if (!al) {
      if (gt_vector_get_used(free_list)) {
        al = *(gt_vector_get_last_elm(free_list, align_details *));
        gt_vector_dec_used(free_list);
        if (al->tag)
          free(al->tag);
        al->tag = 0;
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
      } else {
        al = gt_malloc(sizeof(align_details));
        al->seq_name = al->read[0] = al->read[1] = al->qualities[0] =
            al->qualities[1] = 0;
        al->tag = 0;
        al->mismatches[0] = gt_vector_new(8, sizeof(gt_misms));
        al->mismatches[1] = gt_vector_new(8, sizeof(gt_misms));
      }
    }
    bool reverse;
    error_code = gt_isp_quick_parse_bs_sam_alignment(
        text_line, al, param->mapq_thresh, param->max_template_len, &reverse);
    gt_input_sam_parser_next_record(buffered_sam_input);
    if (error_code) {
      if (error_code == GT_ISP_SAM_FILTERED) {
//				fprintf(stderr,"A: FILTERED!\n");
        continue;
      }
      gt_input_sam_parser_prompt_error(buffered_sam_input, line_num,
                                       buffered_sam_input->cursor - line_start,
                                       error_code);
      return GT_ISP_FAIL;
    }
    gt_input_parse_tag_chomp_pairend_info(tag);
    if (al->tag)
      free(al->tag);
    al->tag = strdup(gt_string_get_string(tag));
    bool new_block = false;
    if (!curr_ctg || strncmp(curr_ctg, gt_string_get_string(al->seq_name),
                             gt_string_get_length(al->seq_name))) {
      chr_skip = false;
      if (curr_ctg)
        gt_free(curr_ctg);
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
      fprintf(stderr, "Processing chromosome %s (%s)\n", curr_ctg,
              chr_skip ? "SKIP" : "OK");
      new_block = true;
    }
    if (chr_skip)
      continue;
    bool insert;
		if(al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) {
			if (reverse)
				insert = al->forward_position > al->reverse_position ? true : false;
			else
				insert = al->forward_position < al->reverse_position ? true : false;
			if (al->forward_position == al->reverse_position) {
				insert = false;
				if (!new_block) {
					align_details *thash;
					HASH_FIND(hh, align_hash, al->tag, gt_string_get_length(tag), thash);
					if (!thash)
						insert = true;
				} else
					insert = true;
			}
		} else {
			insert = true;
		}
		if (new_block == false && insert == true) {
			if(al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) {
				if (al->forward_position > max_pos && al->reverse_position > max_pos) {
					new_block = true;
				}
			} else {
				if (al->forward_position > max_pos || al->reverse_position > max_pos) {
					new_block = true;
				}
			}
		}
    if (new_block == true) {
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
        process_template_vector(align_list, curr_ctg, max_pos, param);
        uint64_t used = gt_vector_get_used(free_list);
        gt_vector_reserve(free_list, ix + used, false);
        memcpy(gt_vector_get_mem(free_list, align_details *)+used,
               gt_vector_get_mem(align_list, align_details *),
               sizeof(align_details *) * ix);
        gt_vector_set_used(free_list, ix + used);
        gt_vector_clear(align_list);
      }
//      printf("Reading new block\n");
      max_pos = 0;
    }
    uint64_t x;
		if(al->alignment_flag & GT_SAM_FLAG_MULTIPLE_SEGMENTS) {
			x = al->forward_position + llabs(al->template_len) - 1;
			if (x > max_pos) max_pos = x;
			if (insert == false) {
				align_details *thash;
				HASH_FIND(hh, align_hash, al->tag, gt_string_get_length(tag), thash);
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
				} else {
					fprintf(stdout, "Warning not found: " PRIgts " %" PRIu64 " %" PRIu64 " %c\n",
									PRIgts_content(tag), al->forward_position, al->reverse_position,
									al->orientation == FORWARD ? '+' : '-');
				}
			} else {
				// Here we have a forward facing pair, so we need to store end to be
				// matched up later
				al->idx = read_idx++; // Preserve read order after matching
				align_details *thash;
				HASH_FIND(hh, align_hash, al->tag, gt_string_get_length(tag), thash);
				gt_cond_fatal_error(thash != NULL, PARSE_SAM_DUPLICATE_SEQUENCE_TAG,
														PRIgts_content(tag));
				HASH_ADD_KEYPTR(hh, align_hash, al->tag, gt_string_get_length(tag), al);
				al = 0;
			}
		} else { // Single (non-paired) reads
			if(al->forward_position > 0) x = al->forward_position + al->align_length;
			else x = al->reverse_position + al->align_length;
//			fprintf(stderr,"%lu %lu %lu %lu %lu\n", al->forward_position, al->reverse_position, al->align_length, x, max_pos);
			if(x > max_pos) max_pos = x;
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
#endif
gt_sequence_archive *gt_open_sequence_archive(void) {
  gt_sequence_archive *sequence_archive = NULL;
  if (param.name_gem_index_file != NULL) { // Load GEM-IDX
    sequence_archive = gt_sequence_archive_new(GT_BED_ARCHIVE);
    gt_gemIdx_load_archive(param.name_gem_index_file, sequence_archive, true);
  } else {
    gt_input_file *const reference_file =
        gt_input_file_open(param.name_reference_file, false);
    sequence_archive = gt_sequence_archive_new(GT_CDNA_ARCHIVE);
    if (gt_input_multifasta_parser_get_archive(reference_file,
                                               sequence_archive) != GT_IFP_OK) {
      gt_fatal_error_msg("Error parsing reference file '%s'\n",
                         param.name_reference_file);
    }
    gt_input_file_close(reference_file);
  }
  return sequence_archive;
}

void print_vcf_header(sr_param *param) {
  FILE *fp = param->output_file;
  fputs("##fileformat=VCFv4.1x\n", fp);
  time_t cl = time(0);
  struct tm *tt = localtime(&cl);
  fprintf(fp, "##fileDate(dd/mm/yyyy)=%02d/%02d/%04d\n", tt->tm_mday,
          tt->tm_mon + 1, tt->tm_year + 1900);
  fprintf(fp, "##source=bs_call_v2.0,under_conversion=%g,over_conversion=%g,"
              "mapq_thresh=%d,bq_thresh=%d\n",
          param->under_conv, param->over_conv, param->mapq_thresh,
          param->min_qual);
  GT_VECTOR_ITERATE(param->sam_headers->sequence_dictionary, header_record_p,
                    line_num, gt_sam_header_record *) {
    gt_sam_header_record *header_record = *header_record_p;
    gt_string *ctg = gt_shash_get(header_record, "SN", gt_string);
    gt_string *len = gt_shash_get(header_record, "LN", gt_string);
    if (ctg && len) {
      gt_string *sp = gt_shash_get(header_record, "SP", gt_string);
      if (param->species_filter) {
        if (!sp)
          continue;
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
  fputs("##INFO=<ID=CX,Number=1,Type=String,Description=\"5 base sequence "
        "context (from position -2 to +2 on the positive strand) determined "
        "from the reference\">\n",
        fp);
  fputs("##FILTER=<ID=q20,Description=\"Quality below 20\">\n", fp);
  fputs("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n", fp);
  fputs("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype "
        "Likelihood\">\n",
        fp);
  fputs("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">\n",
        fp);
  fputs("##FORMAT=<ID=MC8,Number=8,Type=Integer,Description=\"Base counts "
        "non-informative for methylation (ACGT) followed by informative for "
        "methylation (ACGT)\">\n",
        fp);
  fputs("##FORMAT=<ID=CS,Number=1,Type=String,Description=\"Strand of Cytosine "
        "relative to reference sequence (+/-/+-/NA)\">\n",
        fp);
  fputs("##FORMAT=<ID=CG,Number=1,Type=String,Description=\"CpG Status (from "
        "genotype calls: Y/N/H/?)\">\n",
        fp);
  fputs("##FORMAT=<ID=CX,Number=1,Type=String,Description=\"5 base sequence "
        "context (from position -2 to +2 on the positive strand) determined "
        "from genotype call\">\n",
        fp);
  fputs("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT", fp);
  if (param->sample_name)
    fprintf(fp, "\t%s\n", param->sample_name);
  else
    fputs("\tSAMPLE\n", fp);
}

gt_status bs_call_process(sr_param *param) {
  gt_status err = GT_STATUS_OK;
  param->output_file = stdout;
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

  gt_buffered_input_file *buffered_input =
      gt_buffered_input_file_new(input_file);
  gt_sam_parser_attributes *input_sam_attributes =
      gt_input_sam_parser_attributes_new(param->is_paired);

  fill_base_prob_table(param->under_conv, param->over_conv);

    gt_vector *templates = gt_vector_new(32, sizeof(gt_template *));
    error_code = input_sam_parser_get_template_vector(
        buffered_input, templates, input_sam_attributes, param);

#if 0
  if (param->is_paired) {
    gt_vector *templates = gt_vector_new(32, sizeof(gt_template *));
    error_code = input_sam_parser_get_template_vector(
        buffered_input, templates, input_sam_attributes, param);
//    printf("error_code=%d\n", error_code);
  } else {
    gt_alignment *alignment = gt_alignment_new();
    while ((error_code = gt_input_sam_parser_get_single_alignment(
                buffered_input, alignment, input_sam_attributes))) {
      if (error_code != GT_IMP_OK)
        break;
    }
  }
#endif
  gt_input_sam_parser_attributes_delete(input_sam_attributes);
  gt_buffered_input_file_close(buffered_input);
  gt_input_file_close(input_file);
	if(param->pileup) fclose(param->pileup_file);
  return err;
}

gt_status parse_arguments(int argc, char **argv) {
  gt_status err = GT_STATUS_OK;
  struct option *bs_call_getopt = gt_options_adaptor_getopt(bs_call_options);
  gt_string *const bs_call_short_getopt =
      gt_options_adaptor_getopt_short(bs_call_options);
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
    case 's':
      param.extra_stats = true;
      break;
		 case 'P':
			if(param.pileup_file_name != NULL) free(param.pileup_file_name);
			param.pileup_file_name = strdup(optarg);
	 	  param.pileup = true;
			break;
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
    case 'p':
      param.is_paired = true;
      break;
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
    case 'I':
      param.name_gem_index_file = optarg;
      load_seq = 1;
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
    /* Model */
    case 'c':
      if (sscanf(optarg, "%lf,%lf", &param.under_conv, &param.over_conv) != 2)
        gt_fatal_error_msg(
            "c (conversion) option requires two comma separated arguments");
      break;
    case 301:
      param.caller = graphical_model;
      break;
    case 302:
      param.caller = maximum_likelihood;
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
  if (optind < argc)
    param.input_file = argv[optind];
  if (load_seq) {
    if (param.verbose)
      printf("Loading reference sequences:");
    fflush(stdout);
    param.sequence_archive = gt_open_sequence_archive();
    if (param.verbose)
      printf(" Done\n");
  } else {
    fputs("Error in bs_calls(): a sequence archive or GEM index is mandatory\n",
          stderr);
    err = GT_STATUS_FAIL;
  }
  // Free
  gt_string_delete(bs_call_short_getopt);
  // Sanity
  if (param.min_qual < 1)
    param.min_qual = 1;
  else if (param.min_qual > MAX_QUAL)
    param.min_qual = MAX_QUAL;
  if (param.num_threads < 1)
    param.num_threads = 1;
  return err;
}

int main(int argc, char *argv[]) {
  gt_status err = 0;
  // GT error handler
  gt_handle_error_signals();

  // Parsing command-line options
  err = parse_arguments(argc, argv);
  if (err == GT_STATUS_OK)
    lfact_store_init();
  err = bs_call_process(&param);
  if (param.sequence_archive)
    gt_sequence_archive_delete(param.sequence_archive);
  return err == GT_STATUS_OK ? 0 : -1;
}
