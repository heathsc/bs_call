/*
 * PROJECT: GEM-Tools library
 * FILE: gt_input_sam_parser.h
 * DATE: 17/07/2012
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 * DESCRIPTION: // TODO
 */


#ifndef GT_INPUT_SAM_PARSER_H_
#define GT_INPUT_SAM_PARSER_H_

#include "gt_commons.h"
#include "gt_dna_string.h"
#include "gt_alignment_utils.h"
#include "gt_template_utils.h"

#include "gt_sequence_archive.h"

#include "gt_input_file.h"
#include "gt_buffered_input_file.h"
#include "gt_input_parser.h"
#include "gt_input_fasta_parser.h"

#include "gt_sam_attributes.h"

// Codes gt_status
#define GT_ISP_OK   GT_STATUS_OK
#define GT_ISP_FAIL GT_STATUS_FAIL
#define GT_ISP_EOF  0

/*
 * Parsing error/state codes
 */
#define GT_ISP_PE_WRONG_FILE_FORMAT 10
#define GT_ISP_PE_PREMATURE_EOL 11
#define GT_ISP_PE_EXPECTED_NUMBER 12
#define GT_ISP_PE_BAD_CHARACTER 14
#define GT_ISP_PE_WRONG_READ_CONTENT 15
/* CIGAR */
#define GT_ISP_PE_CIGAR_PREMATURE_END 20
#define GT_ISP_PE_SAM_UNMAPPED_XA 21
/* PairedEnd Parsing */
#define GT_ISP_PE_WRONG_NUM_XA 30
#define GT_ISP_PE_UNSOLVED_PENDING_MAPS 32

#define GT_ISP_SAM_FILTERED 40

/*
 * SAM file format constants
 */
#define GT_SAM_HEADER_BEGIN '@'

typedef struct {
  bool parse_optional_fields; // TODO
  bool sam_soap_style;
} gt_sam_parser_attributes;
#define GT_SAM_PARSER_ATTR_DEFAULT { .sam_soap_style=false }

GT_INLINE gt_sam_parser_attributes* gt_input_sam_parser_attributes_new();
GT_INLINE void gt_input_sam_parser_attributes_delete(gt_sam_parser_attributes* const attributes);
GT_INLINE void gt_input_sam_parser_attributes_reset_defaults(gt_sam_parser_attributes* const attributes);
GT_INLINE void gt_input_sam_parser_attributes_set_soap_compilant(gt_sam_parser_attributes* const attributes);

/*
 * SAM File basics
 */
GT_INLINE bool gt_input_file_test_sam(
    gt_input_file* const input_file,gt_sam_headers* const sam_headers,const bool show_errors);
GT_INLINE void gt_input_sam_parser_prompt_error(
    gt_buffered_input_file* const buffered_map_input,
    uint64_t line_num,uint64_t column_pos,const gt_status error_code);
GT_INLINE void gt_input_sam_parser_next_record(gt_buffered_input_file* const buffered_map_input);
GT_INLINE gt_status gt_input_file_sam_read_headers(
    char* const buffer,const uint64_t buffer_size,gt_sam_headers* const sam_headers,
    uint64_t* const characters_read,uint64_t* const lines_read);
GT_INLINE gt_status gt_input_sam_parser_reload_buffer(gt_buffered_input_file* const buffered_sam_input);
GT_INLINE gt_status gt_input_sam_parser_check_sam_file_format(gt_buffered_input_file* const buffered_sam_input);
GT_INLINE uint64_t gt_isp_read_tag(const char** const init_text_line,const char** const end_text_line,gt_string* const tag);


/*
 * High Level Parsers
 */
GT_INLINE gt_status gt_input_sam_parser_get_template(
    gt_buffered_input_file* const buffered_map_input,gt_template* const template,gt_sam_parser_attributes* const attributes);
GT_INLINE gt_status gt_input_sam_parser_get_alignment(
						      
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_sam_parser_attributes* const attributes);
GT_INLINE gt_status gt_input_sam_parser_get_single_alignment(
    gt_buffered_input_file* const buffered_map_input,gt_alignment* const alignment,gt_sam_parser_attributes* const attributes);

/*
 * These are for the Bisulphite specific routines
 */

typedef enum { NON_CONVERTED, STRAND_C2T, STRAND_G2A } gt_bs_strand;

typedef enum { gt_flt_none = 0, gt_flt_unmapped, gt_flt_qc, gt_flt_secondary, gt_flt_mate_unmapped, gt_flt_duplicate, gt_flt_nopos, gt_flt_nomatepos, gt_flt_mismatch_chr, gt_flt_orientation, gt_flt_insert_size, gt_flt_noseq, gt_flt_mapq,gt_flt_not_correctly_aligned } gt_filter_reason;

typedef struct {
  uint64_t alignment_flag;
  uint64_t forward_position;
  uint64_t reverse_position;
	uint64_t reference_span[2];
  uint64_t idx;
	uint64_t align_length;
  int64_t template_len;
	uint64_t trim_left[2];    // Note by how much the reads have been trimmed from it's original length
	uint64_t trim_right[2];
  gt_string* tag;
  gt_string* seq_name;
  gt_string* read[2];
  gt_string* qualities[2];
  gt_vector* orig_pos[2];
  gt_vector* mismatches[2];
  uint8_t mapq[2];
  uint8_t sam_tp;
  uint8_t sam_tq;
	gt_filter_reason filtered;
  gt_strand orientation;
  gt_bs_strand bs_strand;
  UT_hash_handle hh;  
} align_details;

GT_INLINE gt_status gt_isp_quick_parse_bs_sam_alignment(const char** const text_line, align_details *al, const uint64_t thresh, const uint64_t max_template_len, bool keep_unmatched, bool* reverse);

#endif /* GT_INPUT_SAM_PARSER_H_ */
