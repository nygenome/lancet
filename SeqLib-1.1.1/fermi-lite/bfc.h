#ifndef AC_BFC_H__
#define AC_BFC_H__

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "htab.h"
#include "kmer.h"
#include "internal.h"
#include "fml.h"
#include "khash.h"

#define _cnt_eq(a, b) ((a)>>14 == (b)>>14)
#define _cnt_hash(a) ((a)>>14)
KHASH_INIT(cnt, uint64_t, char, 0, _cnt_hash, _cnt_eq)
typedef khash_t(cnt) cnthash_t;

/********************
 * Correct one read *
 ********************/

#include "ksort.h"

#define ECCODE_MISC      1
#define ECCODE_MANY_N    2
#define ECCODE_NO_SOLID  3
#define ECCODE_UNCORR_N  4
#define ECCODE_MANY_FAIL 5

typedef struct {
	uint32_t ec_code:3, brute:1, n_ec:14, n_ec_high:14;
	uint32_t n_absent:24, max_heap:8;
} ecstat_t;

typedef struct {
	uint8_t ec:1, ec_high:1, absent:1, absent_high:1, b:4;
} bfc_penalty_t;



struct bfc_ch_s {
  int k;
  cnthash_t **h;
  // private
  int l_pre;
};

typedef struct {
	int n_threads, q, k, l_pre;
	int min_cov; // a k-mer is considered solid if the count is no less than this

	int max_end_ext;
	int win_multi_ec;
	float min_trim_frac;

	// these ec options cannot be changed on the command line
	int w_ec, w_ec_high, w_absent, w_absent_high;
	int max_path_diff, max_heap;
} bfc_opt_t;

/**********************
 *** K-mer counting ***
 **********************/

/***************
 *** Correct ***
 ***************/

#define BFC_MAX_KMER     63
#define BFC_MAX_BF_SHIFT 37

#define BFC_MAX_PATHS 4
#define BFC_EC_HIST 5
#define BFC_EC_HIST_HIGH 2

#define BFC_EC_MIN_COV_COEF .1

/**************************
 * Sequence struct for ec *
 **************************/

#include "kvec.h"

typedef struct { // NOTE: unaligned memory
	uint8_t b:3, q:1, ob:3, oq:1;
	uint8_t dummy;
	uint16_t lcov:6, hcov:6, solid_end:1, high_end:1, ec:1, absent:1;
	int i;
} ecbase_t;

typedef kvec_t(ecbase_t) ecseq_t;

#define CNT_BUF_SIZE 256

typedef struct { // cache to reduce locking
	uint64_t y[2];
	int is_high;
} insbuf_t;

typedef struct {
	int k, q;
	int n_seqs;
	const fseq1_t *seqs;
	bfc_ch_t *ch;
	int *n_buf;
	insbuf_t **buf;
} cnt_step_t;

float fml_correct_core(const fml_opt_t *opt, int flt_uniq, int n, fseq1_t *seq);



typedef struct {
	int tot_pen;
	int i; // base position
	int k; // position in the stack
	int32_t ecpos_high[BFC_EC_HIST_HIGH];
	int32_t ecpos[BFC_EC_HIST];
	bfc_kmer_t x;
} echeap1_t;

typedef struct {
	int parent, i, tot_pen;
	uint8_t b;
	bfc_penalty_t pen;
	uint16_t cnt;
} ecstack1_t;

typedef struct {
	const bfc_opt_t *opt;
	const bfc_ch_t *ch;
	kvec_t(echeap1_t) heap;
	kvec_t(ecstack1_t) stack;
	ecseq_t seq, tmp, ec[2];
	int mode;
	ecstat_t ori_st;
} bfc_ec1buf_t;


/********************
 * Error correction *
 ********************/

typedef struct {
	const bfc_opt_t *opt;
	const bfc_ch_t *ch;
	bfc_ec1buf_t **e;
	int64_t n_processed;
	int n_seqs, flt_uniq;
	fseq1_t *seqs;
} ec_step_t;

void kmer_correct(ec_step_t * es, int mode, bfc_ch_t * ch);
void bfc_opt_init(bfc_opt_t *opt);

#endif
