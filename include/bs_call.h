#ifndef BS_CALL_H

#define STRING_EXP(tok) #tok
#define STRING(tok) STRING_EXP(tok)

#define DEFAULT_MAPQ_THRESH 20
#define DEFAULT_MAX_TEMPLATE_LEN 1000
#define DEFAULT_REALIGN_TOL 8 // Allow reads to align within this many bp of reported position
#define DEFAULT_UNDER_CONVERSION 0.01
#define DEFAULT_OVER_CONVERSION 0.05
#define DEFAULT_REF_BIAS 2

#define MAX_GAUSS_N 256
#define MAX_QUAL 43
#define MIN_QUAL 20
#define QUAL_CONV 33
#define MAX_ITER 15
#define ITER_FIN (1.0e-8)

#define LOG10 (2.30258509299404568402)
#define LOG2 (0.69314718055994530942)

typedef struct {
  uint64_t mask;
  int n_entries;
  uint16_t *entries;
  uint8_t *name_buf;
} dbsnp_bin;

typedef struct {
  char *name;
  int min_bin;
  int max_bin;
  dbsnp_bin *bins;
  UT_hash_handle hh;
} dbsnp_ctg;

typedef enum {stats_all = 0, stats_passed} stats_cat;
typedef enum {all_sites = 0, variant_sites, CpG_ref_sites, CpG_nonref_sites} qual_cat;
typedef enum {mut_AC = 0, mut_AG, mut_AT, mut_CA, mut_CG, mut_CT, mut_GA, mut_GC, mut_GT, mut_TA, mut_TC, mut_TG, mut_no} stats_mut;
typedef enum {graphical_model, maximum_likelihood} BS_Caller;
typedef enum {base_none = 0, base_trim, base_clip, base_overlap, base_lowqual} base_filter_types;

typedef struct {
  char *ctg;
  bool flag;
  UT_hash_handle hh;
} ctg_hash;

typedef struct {
  char *contig;
  uint64_t snps[2];
  uint64_t indels[2];
  uint64_t multi[2];
  uint64_t dbSNP_sites[2]; // How many dbSNP sites are covered
  uint64_t dbSNP_var[2]; // How many dbSNP sites are variant
  uint64_t CpG_ref[2];
  uint64_t CpG_nonref[2];
  int nbins;
  uint8_t *gc;
  UT_hash_handle hh;
} gt_ctg_stats;

typedef struct {
  uint64_t coverage;
  uint64_t var;
  uint64_t CpG[2];
  uint64_t CpG_inf[2];
  uint64_t all;
  uint64_t gc_pcent[101];
  UT_hash_handle hh;
} gt_cov_stats;

typedef struct {
  uint64_t conv_cts[4];
} meth_cts;

typedef struct {
  uint64_t cts[2];
} fstats_cts;

typedef struct {
  uint64_t snps[2];
  uint64_t indels[2];
  uint64_t multi[2];
  uint64_t dbSNP_sites[2]; // How many dbSNP sites are covered
  uint64_t dbSNP_var[2]; // How many dbSNP sites are variant
  uint64_t CpG_ref[2];
  uint64_t CpG_nonref[2];
  uint64_t mut_counts[12][2];
  uint64_t dbSNP_mut_counts[12][2];
  uint64_t qual[4][256];
  uint64_t filter_cts[15];
  uint64_t filter_bases[14];
  uint64_t base_filter[5];
  gt_vector *meth_profile;
  gt_vector *fs_stats;
  gt_vector *qd_stats;
  gt_vector *mq_stats;
  gt_vector *gof_stats;
  uint64_t filter_counts[2][32];
  double CpG_ref_meth[2][101];
  double CpG_nonref_meth[2][101];
  gt_cov_stats *cov_stats;
} bs_stats;

typedef struct {
  char *input_file;
  char *name_reference_file;
  //  char *name_gem_index_file;
  char *output_prefix;
  char *sample_name;
  char *dbSNP_name;
  char *report_file;
  char *contig_bed;
  ctg_hash *contig_hash;
  bool mmap_input;
  /* Control flags */
  bool no_split;
  bool extra_stats;
  bool keep_duplicates;
  bool ignore_duplicates;
  bool keep_unmatched;
  bool haploid;
  bool verbose;
  bool blank_trim;
  //	bool pileup;
  bool all_positions;
  BS_Caller caller;
  int left_trim;
  int right_trim;
  //	char *pileup_file_name;
  FILE *output_file;
  FILE *json_file;
  //  FILE *pileup_file;
  gt_sam_headers *sam_headers;
  uint8_t mapq_thresh;
  uint8_t min_qual;
  uint64_t max_template_len;
  uint64_t realign_tol;
  double under_conv, over_conv;
  double ref_bias;
  gt_output_file_compression compress;
  gt_generic_printer_attributes *printer_attr;
  gt_buffered_output_file *buf_output;
  int num_threads;
  gt_sequence_archive *sequence_archive;
  dbsnp_ctg *dbSNP;
  uint16_t n_dbSNP_prefixes;
  char **dbSNP_prefix;
  char *dbSNP_header;
  bs_stats *stats;
  gt_ctg_stats *ctg_stats;
} sr_param;

typedef struct {
  double e, k, ln_k, ln_k_half, ln_k_one;
} qual_prob;

typedef struct {
  uint64_t counts[8];
  int qual[8]; // Average quality per base type
  double gt_prob[10]; // Genotype log probabilities (Log10)
  double gt_gof; // Goodness of fit LR (Log10) 
  double fisher_strand; // Allele strand bias LR (Log10) 
  int mq; // RMS Mapping quality
  int aq; // Average base quality
  uint8_t max_gt;
} gt_meth;

typedef struct {
  gt_meth gtm;
  bool ready;
  bool skip;
} gt_vcf;

typedef struct _base_counts {
  struct _base_counts *next;
  int8_t idx[MAX_QUAL + 1];
  uint32_t counts[MAX_QUAL + 1];
} base_counts;

typedef struct {
  uint32_t counts[2][8];
  uint32_t n;
  float quality[8];
  float mapq2;
  uint8_t *seq;
  size_t seq_size;
  size_t seq_idx;
} pileup;

void fill_base_prob_table(void);
void calc_gt_prob(gt_meth *gt, sr_param *param, char rf);

#define BS_CALL_H 1
#endif
