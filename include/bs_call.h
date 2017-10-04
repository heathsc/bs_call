#ifndef BS_CALL_H

#define STRING_EXP(tok) #tok
#define STRING(tok) STRING_EXP(tok)

#define DEFAULT_MAPQ_THRESH 20
#define DEFAULT_MAX_TEMPLATE_LEN 1000
#define DEFAULT_REALIGN_TOL 8 // Allow reads to align within this many bp of reported position
#define DEFAULT_UNDER_CONVERSION 0.01
#define DEFAULT_OVER_CONVERSION 0.05

#define MAX_GAUSS_N 256
#define MAX_QUAL 43
#define MIN_QUAL 20
#define QUAL_CONV 33
#define MAX_ITER 15
#define ITER_FIN (1.0e-8)

typedef enum {graphical_model, maximum_likelihood} BS_Caller;

/*
* bs_call options
*/
gt_option bs_call_options[] = {
  /* Operations */
  { 'N', "no-split", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Do not split output on contig"},
  { '1', "haploid", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Assume genome is haploid"},
  { 'd', "keep-duplicates", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Don't merge duplicate reads"},
  { 's', "extra-stats", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Generate extra stats files"},
  { 'R', "right-trim", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "", "Bases to trim from right of read pair"},
  { 'L', "left-trim", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "", "Bases to trim from left of read pair"},
  { 'B', "blank-trim", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 1, true, "", "Don't use trimmed bases for genotype estimation"},
  { 'q', "mapq-threshold", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set MAPQ threshold for selecting reads (default "STRING(DEFAULT_MAPQ_THRESH)")"},
  { 'Q', "bq-threshold", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set base quality threshold for calling (default "STRING(MIN_QUAL)")"},
	{ 'l', "max-template-length", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Set maximum template length for a pair (default "STRING(DEFAULT_MAX_TEMPLATE_LEN)")"},
  { 'T', "realign-tolerance", GT_OPT_REQUIRED, GT_OPT_INT, 1, true, "<int> ","Tolerance for realignment positions (default "STRING(DEFAULT_REALIGN_TOL)")"},
  /* I/O */
  { 'p', "paired-end", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
  { 'z', "gzip", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
  { 'j', "bzip2", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
  { 'Z', "no-compress", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, true, "" , "" },
  { 201, "mmap-input", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 2, false, "" , "" },
  { 'o', "output", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<output prefix>" , "" },
  { 'P', "pileup", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<pileup file name>" , "" },
  { 'n', "sample", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<sample name>" , "SAMPLE" },
  { 'x', "species", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<species filter" , "" },	
  { 'r', "reference", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file> (MultiFASTA/FASTA)" , "" },
  { 'I', "gem-index", GT_OPT_REQUIRED, GT_OPT_STRING, 2, true, "<file> (GEM2-Index)" , "" },
	/* Model */
  { 'c', "conversion", GT_OPT_REQUIRED, GT_OPT_INT, 3, true, "<float>,<float>","Set under and over conversion rates (default "STRING(DEFAULT_UNDER_CONVERSION)","STRING(DEFAULT_OVER_CONVERSION)")"},	
  { 301, "graphical_model", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "" },
  { 302, "maximum_likelihood", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 3, true, "", "" },
	/* Misc */
  { 'v', "verbose", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  { 't', "threads", GT_OPT_REQUIRED, GT_OPT_INT, 4, false, "", ""},
  { 'h', "help", GT_OPT_NO_ARGUMENT, GT_OPT_NONE, 4, true, "", ""},
  {  0, 0, 0, 0, 0, false, "", ""}
};

char* bs_call_options_short = "N1dsBR:L:pzjZo:r:I:vt:h";

char* bs_call_groups[] = {
  /*  0 */ "Null",
  /*  1 */ "Operations",
  /*  2 */ "I/O",
  /*  3 */ "Model",
  /*  4 */ "Misc",
};

// Weights for SW alignment (same as BWA)

#define BS_CALL_MATCH 1
#define BS_CALL_MISM -4
#define BS_CALL_GAP_OPEN 6
#define BS_CALL_GAP_EXTEND 1
#define BS_CALL_QUAL_CUTOFF 26 // For realignment

int align_ss2(int16_t **tref,int rl,gt_string *query,int mm,int v,int u,int mat,int x,int tol);
int band_align_ss2(const char *ref,int rl,gt_string *query,int mm,int v,int u,int mat,int x);

#define BS_CALL_H 1
#endif
