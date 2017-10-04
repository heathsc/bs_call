/*
 * PROJECT: GEM-Tools library
 * FILE: gt.coverage.c
 * DATE: 10/03/2014
 * AUTHOR(S): Simon Heath <simon.heath@gmail.com>
 * DESCRIPTION: Calculate coverage from unsorted BAMs of GEM map files
 */

#include <getopt.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "gem_tools.h"

typedef enum {GT_UNIQUE_FIRST, GT_UNIQUE_GEM, GT_UNIQUE_XT, GT_UNIQUE_MQ} gt_unique_t;
typedef enum {GT_DUPLICATE_IGNORE, GT_DUPLICATE_REMOVE} gt_duplicate_t;

typedef struct {
  char *ranges_file;
  char *target_file;
  char *output_file;
  char *prefix;
  bool detailed_output;
  bool combine;
  bool paired;
  bool verbose;
  gt_unique_t unique;
  gt_duplicate_t duplicate;
  uint64_t min_quality;
  uint64_t MQ_limit;
  uint64_t extend;
  uint64_t block_size;
  uint64_t num_threads;
} gt_coverage_args;

gt_coverage_args parameters = {
  .ranges_file=NULL,
  .target_file=NULL,
  .output_file=NULL,
  .prefix=NULL,
  .detailed_output=false,
  .combine=false,
  .paired=false,
  .unique=GT_UNIQUE_GEM,
  .duplicate=GT_DUPLICATE_REMOVE,
  .min_quality=0,
  .MQ_limit=0,
  .verbose=false,
  .extend=0,
  .block_size=1,
  .num_threads=1
};

void usage(const gt_option* const options,char* groups[],const bool print_inactive) {
  fprintf(stderr, "USE: ./gt_coverage [ARGS] [FILES]...\n");
  gt_options_fprint_menu(stderr,options,groups,false,print_inactive);
}

gt_status parse_arguments(int argc,char** argv) {
  gt_status err=GT_STATUS_OK;
  struct option* gt_coverage_getopt = gt_options_adaptor_getopt(gt_coverage_options);
  gt_string* const gt_coverage_short_getopt = gt_options_adaptor_getopt_short(gt_coverage_options);
  int option, option_index;

  while (true) {
    if ((option=getopt_long(argc,argv,
                            gt_string_get_string(gt_coverage_short_getopt),gt_coverage_getopt,&option_index))==-1) break;
    switch (option) {
    case 'c':
      parameters.combine=true;
      break;
    case 'd':
      parameters.detailed_output=true;
      break;
    case 'b':
      parameters.block_size = atol(optarg);
      break;
    case 'o':
      parameters.output_file = optarg;
      break;
    case 'p':
      parameters.prefix = optarg;
      break;
    case 'r':
      parameters.ranges_file = optarg;
      break;
    case 't':
      parameters.target_file = optarg;
      break;
    case 'x':
      parameters.extend = atol(optarg);
      break;
    case 'P':
      parameters.paired = true;
      break;
    case 'U':
      if(!strcasecmp("first",optarg)) parameters.unique=GT_UNIQUE_FIRST;
      else if(!strcasecmp("gem",optarg)) parameters.unique=GT_UNIQUE_GEM;
      else if(!strcasecmp("xt",optarg)) parameters.unique=GT_UNIQUE_XT;
      else if(!strncasecmp("mq:",optarg,3)) {
	parameters.unique=GT_UNIQUE_MQ;
	parameters.MQ_limit=atol(optarg+3);
      } else {
	gt_fatal_error_msg("Unique mode not recognized: '%s'",optarg);
      }
      break;
    case 'D':
      if(!strcasecmp("ignore",optarg)) parameters.duplicate=GT_DUPLICATE_IGNORE;
      else if(!strcasecmp("remove",optarg)) parameters.duplicate=GT_DUPLICATE_REMOVE;
      else {
	gt_fatal_error_msg("Duplicate mode not recognized: '%s'",optarg);
      }
      break;
    case 'q':
      parameters.min_quality = atol(optarg);
      break;
    case 'v':
      parameters.verbose=true;
      break;
    case 'T':
#ifdef HAVE_OPENMP
      parameters.num_threads = atol(optarg);
#endif
      break;
    case 'h':
      usage(gt_coverage_options,gt_coverage_groups,false);
      exit(1);
      break;
    case 'H':
      usage(gt_coverage_options,gt_coverage_groups,true);
      exit(1);
    case 'J':
      gt_options_fprint_json_menu(stderr,gt_coverage_options,gt_coverage_groups,true,false);
      exit(1);
      break;
    case '?':
    default:
      usage(gt_coverage_options,gt_coverage_groups,false);
      gt_fatal_error_msg("Option '%c' %d not recognized",option,option);
      break;
    }
  }
  if(!parameters.ranges_file) {
    fputs("No ranges file specified.  Use -r or --ranges_file option\n",stderr);
  }
  return err;
}

int main(int argc,char** argv) {
  parse_arguments(argc,argv);
  return 0;
}


