//
//  dbSNP_idx.h
//  dbSNP
//
//  Created by Simon Heath on 15/11/2017.
//  Copyright © 2017 Simon Heath. All rights reserved.
//

#ifndef dbSNP_idx_h
#define dbSNP_idx_h

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>
#include <pthread.h>
#include "uthash.h"

#define INITIAL_ENTRY_SIZE 4
#define INITIAL_NAME_BUF_SIZE 32

typedef struct {
   uint64_t mask;
   uint16_t name_buf_size;
   uint16_t name_buf_idx;
   uint8_t n_entries;
   uint8_t entry_size;
   uint16_t *entry;
   unsigned char *name_buf;
} bin;

typedef struct {
   char *name;
   bin *bins;
   uint32_t min_bin;
   uint32_t max_bin;
   UT_hash_handle hh;
} contig;

typedef struct {
  uint16_t ix;
  uint64_t count;
  char *pref;
  UT_hash_handle hh;
} prefix;
  
#define INITIAL
#endif /* dbSNP_idx_h */
