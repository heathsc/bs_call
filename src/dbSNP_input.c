/*
 * dbSNP_input.c
 *
 *  Created on: Feb 2, 2020
 *      Author: heath
 */

#include "dbSNP_idx.h"
#include "dbSNP_utils.h"
#include <sys/types.h>
#include <sys/wait.h>

static file_t *next_file(dbsnp_param_t * const par) {
	if(par->abort_flag) return NULL;
	pthread_mutex_lock(&par->param_mut);
	file_t *file = par->unsorted;
	if(file) HASH_DEL(par->unsorted, file);
	else {
		while(par->unsorted_flag && par->sorted) pthread_cond_wait(&par->param_cond, &par->param_mut);
		file = par->sorted;
		if(file) HASH_DEL(par->sorted, file);
	}
	if(file) {
		if(!file->sorted) par->n_unsorted_being_processed++;
		file->hh.next = par->used_files;
		par->used_files = file;
		pthread_mutex_unlock(&par->param_mut);
		if(strcmp(file->name, "-")) file->fp = open_readfile(file->name, &file->filter);
		else {
			file->fp = stdin;
			file->filter = false;
		}
	} else 	pthread_mutex_unlock(&par->param_mut);
	return file;
}

void free_used_files(dbsnp_param_t * const par) {
	file_t *f = par->used_files, *tmp;
	while(f) {
		tmp = f->hh.next;
		free(f);
		f = tmp;
	}
}

static void set_header(dbsnp_param_t * const par, const char * const hd) {
	pthread_mutex_lock(&par->param_mut);
	if(!par->header) par->header = strdup(hd);
	pthread_mutex_unlock(&par->param_mut);
}

static contig *new_contig(const char * const name, file_t * const file, const uint32_t binx) {
	contig * const ctg = malloc(sizeof(contig));
	ctg->min_bin = ctg->max_bin = binx;
	ctg->bins = malloc(sizeof(bin));
	ctg->name = strdup(name);
	ctg->first_file = file;
	ctg->in_queue = false;
	pthread_mutex_init(&ctg->mut, NULL);
	clear_bins(ctg->bins, 1);
	return ctg;
}

static contig * find_or_create_contig(dbsnp_param_t * const par, const char * const name, file_t * const file, const uint32_t binx) {
	contig *ctg;
	pthread_mutex_lock(&par->param_mut);
	HASH_FIND_STR(par->contigs, name, ctg);
	if(!ctg) {
		ctg = new_contig(name, file, binx);
		HASH_ADD_KEYPTR(hh, par->contigs, ctg->name, strlen(ctg->name), ctg);
	} else if(!par->unsorted_flag) {
		if(file == ctg->first_file) fprintf(stderr, "Contig %s found in multiple places in input '%s' - verify usage of --sorted option\n", ctg->name, file->name);
		else {
			file_t *file1= ctg->first_file;
			if(file1->sorted || !file1->read) {
				fprintf(stderr, "Contig %s found in multiple files ('%s' and '%s') - verify usage of --sorted option\n", ctg->name, file->name, ctg->first_file->name);
				ctg = NULL;
			} else ctg->first_file = file;
		}
	}
	pthread_mutex_unlock(&par->param_mut);
	return ctg;
}
static void add_contig_to_queue(dbsnp_param_t * const par, contig * const ctg) {
	if(!par->unsorted_flag) {
		pthread_mutex_lock(&par->contig_queue_mut);
		ctg->next = par->contig_queue;
		par->contig_queue = ctg;
		ctg->in_queue = true;
		pthread_cond_signal(&par->contig_queue_cond);
		pthread_mutex_unlock(&par->contig_queue_mut);
	}
}

void add_remaining_contigs_to_queue(dbsnp_param_t * const par) {
	for(contig *ctg = par->contigs; ctg; ctg = ctg->hh.next) {
		if(!ctg->in_queue) add_contig_to_queue(par, ctg);
	}
}

void *input_thread(void *pt) {
	dbsnp_param_t * const par = pt;
	bool *st = malloc(sizeof(bool));
	*st = false;
	char *buf = NULL;
	size_t buf_size = 0;
	tokens *tok = NULL;
	contig *ctg = NULL;
	prefix *pref = NULL;
	uint64_t n_snps = 0;
	while(!*st) {
		file_t * const file = next_file(par);
		if(!file) break;
		FILE *fin = file->fp;
		fprintf(stderr, "Reading from %s\n", fin != stdin ? file->name : "<stdin>");
		while(!*st) {
			if(par->abort_flag) {
				*st = true;
				break;
			}
			ssize_t l = getline(&buf, &buf_size, fin);
			if(l < 0) break;
			if(l > 0 && buf[l - 1] == '\n') buf[--l] = 0;
			if(l > 0) {
				if(l > 6 && !strncmp(buf, "track ", 6)) {
					if(par->header == NULL) set_header(par, buf);
				} else {
					tok = tokenize(buf, '\t', tok);
					if(tok->n_tok >= 4) {
						char *p;
						uint32_t x = (uint32_t)strtoul(tok->toks[1], &p, 10);
						uint32_t y = (uint32_t)strtoul(tok->toks[2], &p, 10);
						if(y > x && y - x == 1) {
							uint32_t binx = y >> 6;
							if(ctg == NULL || strcmp(ctg->name, tok->toks[0])) {
								if(ctg) add_contig_to_queue(par, ctg);
								ctg = find_or_create_contig(par, tok->toks[0], file, binx);
								if(!ctg) {
									*st = true;
									par->abort_flag = true;
									break;
								}
							}
							pthread_mutex_lock(&ctg->mut);
							if(binx > ctg->max_bin) {
								ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(binx - ctg->min_bin + 1));
								clear_bins(ctg->bins + (ctg->max_bin - ctg->min_bin + 1), binx - ctg->max_bin);
								ctg->max_bin = binx;
							} else if(binx < ctg->min_bin) {
								ctg->bins = realloc(ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - binx + 1));
								memmove(ctg->bins + (ctg->min_bin - binx), ctg->bins, sizeof(bin) * (size_t)(ctg->max_bin - ctg->min_bin + 1));
								clear_bins(ctg->bins, ctg->min_bin - binx);
								ctg->min_bin = binx;
							}
							n_snps += add_to_bin(ctg->bins + (binx - ctg->min_bin), y & 63, tok->toks[3], &pref, par);
							pthread_mutex_unlock(&ctg->mut);
						}
					}
				}
			}
		}
		if(fin != stdin) {
			fclose(fin);
			if(file->filter) {
				int i;
				while (waitpid(-1, &i, WNOHANG) > 0);
			}
		}
		file->read = true;
		if(!file->sorted) {
			pthread_mutex_lock(&par->param_mut);
			assert(par->unsorted_flag &&& par->n_unsorted_being_processed > 0);
			if(!(--par->n_unsorted_being_processed) && !par->unsorted) {
				par->unsorted_flag = false;
				fprintf(stderr, "Switching to sorted mode (ctg = %s)\n", ctg->name);
				pthread_cond_broadcast(&par->param_cond);
			}
			pthread_mutex_unlock(&par->param_mut);
			ctg = NULL;
		}
	}
	if(ctg) add_contig_to_queue(par, ctg);
	if(buf != NULL) free(buf);
	if(tok != NULL) free_tokens(tok);
	pthread_mutex_lock(&par->param_mut);
	par->n_snps += n_snps;
	pthread_mutex_unlock(&par->param_mut);
	return st;
}
