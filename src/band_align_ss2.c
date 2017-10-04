#include "gem_tools.h"
#include <x86intrin.h>

#define UPDATE_R(r, p, q, r1, bts)                                             \
  {                                                                            \
    r1 = _mm_max_epi8(r, p);                                                   \
    r1 = _mm_max_epi8(r1, q);                                                  \
    __m128i _tt = _mm_cmpeq_epi8(p, r1);                                       \
    bts = _mm_and_si128(_tt, msk_a);                                           \
    _tt = _mm_cmpeq_epi8(q, r1);                                               \
    _tt = _mm_and_si128(_tt, msk_b);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
    _tt = _mm_cmpeq_epi8(r, r1);                                               \
    _tt = _mm_and_si128(_tt, msk_c);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
  }

#define _UPDATE_PQ_a(r, p, q, bts)                                             \
  {                                                                            \
    _tR = _mm_subs_epi8(r, gap_open);                                          \
    __m128i _mx = _mm_max_epi8(p, _tR);                                        \
    __m128i _tt = _mm_cmpeq_epi8(_mx, p);                                      \
    _tt = _mm_and_si128(_tt, msk_d);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
    _tt = _mm_cmpeq_epi8(_mx, _tR);                                            \
    _tt = _mm_and_si128(_tt, msk_e);                                           \
    p = _mm_subs_epi8(_mx, gap_ext);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
  }

#define _UPDATE_PQ_b(r, p, q, bts, ge)                                         \
  {                                                                            \
    __m128i _mx = _mm_max_epi8(q, _tR);                                        \
    __m128i _tt = _mm_cmpeq_epi8(_mx, q);                                      \
    _tt = _mm_and_si128(_tt, msk_f);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
    _tt = _mm_cmpeq_epi8(_mx, _tR);                                            \
    _tt = _mm_and_si128(_tt, msk_g);                                           \
    bts = _mm_or_si128(bts, _tt);                                              \
    q = _mm_subs_epi8(_mx, ge);                                                \
  }

#define UPDATE_PQ(r, p, q, bts)                                                \
  {                                                                            \
    __m128i _tR;                                                               \
    _UPDATE_PQ_a(r, p, q, bts);                                                \
    _UPDATE_PQ_b(r, p, q, bts, gap_ext);                                       \
  }

#define UPDATE_PQ1(r, p, q, bts, go1, ge1)                                     \
  {                                                                            \
    __m128i _tR;                                                               \
    _UPDATE_PQ_a(r, p, q, bts);                                                \
    _tR = _mm_subs_epi8(r, go1);                                               \
    _UPDATE_PQ_b(r, p, q, bts, ge1);                                           \
  }

#define _UPDATE_EDGES(B, B1, B1a, B2)                                          \
  {                                                                            \
    __m128i t1 = _mm_cmpeq_epi8(B, msk_0x7f);                                  \
    B = _mm_andnot_si128(t1, B);                                               \
    t1 = _mm_and_si128(B1a, msk_a);                                            \
    __m128i t2 = _mm_and_si128(B, msk_e);                                      \
    t1 = _mm_or_si128(t1, t2);                                                 \
    t1 = _mm_cmpeq_epi8(t1, msk_0x11);                                         \
    t2 = _mm_and_si128(B1, msk_b);                                             \
    __m128i t3 = _mm_and_si128(B, msk_g);                                      \
    t2 = _mm_or_si128(t2, t3);                                                 \
    t2 = _mm_cmpeq_epi8(t2, msk_0x42);                                         \
    t1 = _mm_or_si128(t1, t2);                                                 \
    t2 = _mm_and_si128(B2, msk_c);                                             \
    t2 = _mm_cmpeq_epi8(t2, msk_c);                                            \
    t1 = _mm_or_si128(t1, t2);                                                 \
    t1 = _mm_or_si128(t1, msk_0x78);                                           \
    B = _mm_and_si128(B, t1);                                                  \
    t1 = _mm_and_si128(B1a, msk_a);                                            \
    t2 = _mm_and_si128(B, msk_d);                                              \
    t1 = _mm_or_si128(t1, t2);                                                 \
    t1 = _mm_cmpeq_epi8(t1, msk_0x09);                                         \
    t2 = _mm_and_si128(B, msk_e);                                              \
    t2 = _mm_srli_epi32(t2, 1);                                                \
    t2 = _mm_xor_si128(t2, msk_d);                                             \
    new_d = _mm_and_si128(t2, t1);                                             \
    t2 = _mm_and_si128(B, msk_a);                                              \
    t2 = _mm_slli_epi32(t2, 4);                                                \
    t2 = _mm_xor_si128(t2, msk_0x11);                                          \
    __m128i new_ae = _mm_and_si128(t2, t1);                                    \
    t1 = _mm_and_si128(B1, msk_b);                                             \
    t2 = _mm_and_si128(B, msk_f);                                              \
    t1 = _mm_or_si128(t1, t2);                                                 \
    t1 = _mm_cmpeq_epi8(t1, msk_0x22);                                         \
    t2 = _mm_and_si128(B, msk_g);                                              \
    t2 = _mm_srli_epi32(t2, 1);                                                \
    t2 = _mm_xor_si128(t2, msk_f);                                             \
    new_f = _mm_and_si128(t2, t1);                                             \
    t2 = _mm_and_si128(B, msk_b);                                              \
    t2 = _mm_slli_epi32(t2, 5);                                                \
    t2 = _mm_xor_si128(t2, msk_0x42);                                          \
    t2 = _mm_and_si128(t2, t1);                                                \
    t2 = _mm_or_si128(t2, new_ae);                                             \
    B = _mm_and_si128(B, msk_0xaf);                                            \
    B = _mm_or_si128(B, t2);                                                   \
  }

#define UPDATE_EDGES(B, B1, B1a, B2)                                           \
  {                                                                            \
    __m128i new_d, new_f;                                                      \
    _UPDATE_EDGES(B, B1, B1a, B2);                                             \
    B1 = _mm_and_si128(B1, msk_0xd7);                                          \
    B1 = _mm_or_si128(B1, new_f);                                              \
    new_d = _mm_slli_si128(new_d, 1);                                          \
    B1 = _mm_or_si128(B1, new_d);                                              \
  }

#define UPDATE_EDGES1(B, B1, B1a, B2)                                          \
  {                                                                            \
    __m128i new_d, new_f;                                                      \
    _UPDATE_EDGES(B, B1, B1a, B2);                                             \
    B1 = _mm_and_si128(B1a, msk_0xd7);                                         \
    B1 = _mm_or_si128(B1, new_d);                                              \
    new_f = _mm_srli_si128(new_f, 1);                                          \
    B1 = _mm_or_si128(B1, new_f);                                              \
  }

#define COMPUTE_SCORE(S, Sx, Sy)                                               \
  {                                                                            \
    __m128i tt = _mm_and_si128(Sx, Sy);                                        \
    tt = _mm_cmpeq_epi8(tt, zero);                                             \
    __m128i tx = _mm_and_si128(tt, mismatch);                                  \
    S = _mm_andnot_si128(tt, match);                                           \
    S = _mm_or_si128(S, tx);                                                   \
    tt = _mm_and_si128(Sx, msk_g);                                             \
    tt = _mm_cmpeq_epi8(tt, zero);                                             \
    S = _mm_and_si128(S, tt);                                                  \
  }

#define CHECK_STATE(S, ret)                                                    \
  {                                                                            \
    __m128i tt = _mm_and_si128(S, msk_0x07);                                   \
    int count = popcnt128(tt);                                                 \
    if (count == 1) {                                                          \
      if (!(_mm_testz_si128(S, msk_b)))                                        \
        ret = 1;                                                               \
      else if (!(_mm_testz_si128(S, msk_c)))                                   \
        ret = 2;                                                               \
      else                                                                     \
        ret = 3;                                                               \
    } else if (count > 1)                                                      \
      ret = 4;                                                                 \
    else                                                                       \
      ret = 0;                                                                 \
  }

static inline int popcnt128(__m128i n) {
  const __m128i n_hi = _mm_unpackhi_epi64(n, n);
  return __builtin_popcountll(_mm_cvtsi128_si64(n)) +
         __builtin_popcountll(_mm_cvtsi128_si64(n_hi));
}

static inline int update_state(__m128i *S, int st, int rt, int cx, int cy,
                               int *cx1, int *cy1, const char *X,
                               const char *Y) {
  static int sx, sy;

  int8_t sts[16];
  _mm_store_si128((__m128i *)sts, *S);
  int k = -1;
  if (rt) {
    for (k = 0; k < 16 && !sts[k]; k++)
      ;
  }
  if (!st) {
    if (rt == 2) {
      printf("US: Starting match: cx=%d, cy=%d, k=%d\n", cx + k, cy - k, k);
      st = 1;
      *cx1 = cx + k;
      *cy1 = cy - k;
    } else if (rt > 2) {
      printf("US: Starting indel: cx=%d, cy=%d\n", *cx1, *cy1);
      sx = *cx1;
      sy = *cy1;
      st = 5;
    } else if (rt == 1) {
      *cx1 = cx + k;
      *cy1 = cy - k;
    }
  } else if (st == 1 || st == 3) {
    if (!rt)
      st++;
    else {
      printf("US: Starting indel: cx=%d, cy=%d\n", *cx1, *cy1);
      sx = *cx1;
      sy = *cy1;
      st = 5;
    }
  } else if (st == 2 || st == 4) {
    if (rt == 2) {
      *cx1 = cx + k;
      *cy1 = cy - k;
      if (!(X[*cx1] & Y[*cy1]))
        printf("US: Mismatch: %d/%d cx=%d, cy=%d\n", Y[*cy1], X[*cx1], *cx1,
               *cy1);
      st--;
    } else {
      if (*cx1 || rt > 1) {
        printf("US: Starting indel: cx=%d, cy=%d\n", *cx1, *cy1);
        sx = *cx1;
        sy = *cy1;
      }
      st = 5;
    }
  } else if (st == 5) {
    if (rt == 2) {
      *cx1 = cx + k;
      *cy1 = cy - k;
      st = 6;
    }
  } else if (st == 6) {
    if (!rt)
      st = 7;
    else
      st = 5;
  } else if (st == 7) {
    if (rt == 2) {
      printf("US: Finishing indel: [%d,%d] cx=%d, cy=%d\n", (sy - *cy1),
             (sx - *cx1), *cx1, *cy1);
      *cx1 = cx + k;
      *cy1 = cy - k;
      if (!(X[*cx1] & Y[*cy1]))
        printf("US: Mismatch: %d/%d cx=%d, cy=%d\n", Y[*cy1], X[*cx1], *cx1,
               *cy1);
      st = 3;
    } else
      st = 5;
  }
  return st;
}

static char *ss2_store;
static __m128i *bits;
static int ss2_store_len, bits_len;

int band_align_ss2(const char *tref, int rl, gt_string *query, int mm, int v,
                   int u, int mat, int x) {
  const char *X = gt_string_get_string(query);
  int m = gt_string_get_length(query);

  // The reference string is 30 bases larger than the query, 15 bases in each
  // direction
  // If the supplied reference is too short it is padded
  // The query string is padded to the same size

  int n = m + 30;

  // Set up padded read (X1) and reference (Y1) strings
  if (ss2_store_len < 2 * n) {
    ss2_store_len = 2 * n;
    if (ss2_store)
      realloc(ss2_store, (size_t)ss2_store_len);
    else
      ss2_store = malloc((size_t)ss2_store_len);
  }
  char *X1 = ss2_store;
  char *Y1 = X1 + n;
  int x1 = x - 15;
  if (x1 < 0)
    x1 = 0;
  int x2 = x + m + 15 - 1;
  if (x2 >= rl)
    x2 = rl - 1;

  const char *Y = tref + x1;
  int i = 0;
  for (i = 0; i < 15 - x; i++)
    Y1[i] = 32;
  memcpy(Y1 + i, Y, x2 - x1 + 1);
  for (i = x2 - x1 + 1; i < n; i++)
    Y1[i] = 32;
  for (i = 0; i < 15; i++)
    X1[i] = X1[m + 15 + i] = 16;
  memcpy(X1 + 15, X, m);

  char *bktrans = "-AC-Ga--T-t----N?---------------!";
  fputs("REF: ", stdout);
  for (int i = 0; i < n; i++)
    fputc(bktrans[Y1[i] & 63], stdout);
  fputc('\n', stdout);
  fputs("RED: ", stdout);
  for (int i = 0; i < n; i++)
    fputc(bktrans[X1[i] & 63], stdout);
  fputc('\n', stdout);

  // Set up bit array
  if ((m + 15) * 2 + 1 > bits_len) {
    if (bits)
      free(bits);
    bits_len = (m * 15) * 2 + 1;
    bits = calloc((size_t)bits_len, sizeof(__m128i));
  }
  memset(bits, 0, ((m + 15) * 2 + 1) * sizeof(__m128i));

  // Required constants
  __m128i zero = _mm_set1_epi32(0);
  __m128i gap_open = _mm_set1_epi8(v);
  __m128i gap_ext = _mm_set1_epi8(u);
  __m128i msk_a = _mm_set1_epi8(1);
  __m128i msk_b = _mm_set1_epi8(2);
  __m128i msk_c = _mm_set1_epi8(4);
  __m128i msk_d = _mm_set1_epi8(8);
  __m128i msk_e = _mm_set1_epi8(16);
  __m128i msk_f = _mm_set1_epi8(32);
  __m128i msk_g = _mm_set1_epi8(64);
  __m128i match = _mm_set1_epi8(mat);
  __m128i mismatch = _mm_set1_epi8(mm);
  __m128i msk_0x07 = _mm_set1_epi8(0x07);
  __m128i msk_0x09 = _mm_set1_epi8(0x09);
  __m128i msk_0x11 = _mm_set1_epi8(0x11);
  __m128i msk_0x22 = _mm_set1_epi8(0x22);
  __m128i msk_0x42 = _mm_set1_epi8(0x42);
  __m128i msk_0x78 = _mm_set1_epi8(0x78);
  __m128i msk_0x7f = _mm_set1_epi8(0x7f);
  __m128i msk_0xd7 = _mm_set1_epi8(0xd7);
  __m128i msk_0xaf = _mm_set1_epi8(0xaf);
  __m128i last_inf =
      _mm_set_epi8(-128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  __m128i first_inf =
      _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -128);
  __m128i last_ones =
      _mm_set_epi8(0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

  // Initialize arrays

  __m128i R, R1, R2, Q, P, S, Sx, Sy;

#if 0
  R = _mm_set_epi8(-128, -128, -(u + u + v), X[0] & Y[0] ? mat : mm, 0, -128,
                   -128, -128, -128, -128, -128, -128, -128, -128, -128, -128);
	R2 = _mm_set_epi8(-128, -128, -128, -(u + v), 0, -128, -128, -128, -128, -128,
                    -128, -128, -128, -128, -128, -128);
  P = _mm_set_epi8(-128, -128, -(u + u + v), -(u + v), -128, -128, -128, -128,
                   -128, -128, -128, -128, -128, -128, -128, -128);
  Q = _mm_set_epi8(-128, -128, -128, -2 * (u + v), -(u + v), -128, -128, -128,
                   -128, -128, -128, -128, -128, -128, -128, -128);
#endif
  R = _mm_set_epi8(X[0] & Y[0] ? mat : mm, 0, -128, -128, -128, -128, -128,
                   -128, -128, -128, -128, -128, -128, -128, -128, -128);
  P = _mm_set_epi8(-(u + v), -128, -128, -128, -128, -128, -128, -128, -128,
                   -128, -128, -128, -128, -128, -128, -128);
  Q = _mm_set_epi8(-2 * (u + v), -(u + v), -128, -128, -128, -128, -128, -128,
                   -128, -128, -128, -128, -128, -128, -128, -128);
  R2 = _mm_set_epi8(-(u + v), 0, -128, -128, -128, -128, -128, -128, -128, -128,
                    -128, -128, -128, -128, -128, -128);

  __m128i edge_mask = _mm_set_epi8(-1, -1, -1, -1, -1, 0, -1, -1, -1, -1, -1,
                                   -1, -1, -1, -1, -1);
  __m128i edge_mask1 = first_inf;
  __m128i edge_mask2 = last_inf;
  __m128i gap_open1 = gap_open;
  __m128i gap_ext1 = gap_ext;

  // First diagonal element
  __m128i bts;
  UPDATE_R(R, P, Q, R1, bts);
  UPDATE_PQ(R1, P, Q, bts);
  _mm_store_si128(bits, bts);

  // Set up initial sequence vectors
  int8_t Sxs[16], Sys[16];
  for (int k = 0; k < 16; k++) {
    Sys[k] = Y1[15 - k];
    Sxs[k] = X1[k + 1];
  }
  Sx = _mm_load_si128((__m128i *)Sxs);
  Sy = _mm_load_si128((__m128i *)Sys);
  __m128i *bitp = bits + 1;

  int n2 = m+15;
  int nx = 0;
  int ny = 0;
  int flag = 0;
  int ix;
  // Remaining diagonal elements (two updates per element)
  for (ix = 1; ix <= n2; ix++) {
    // Shift Q vector right 1 byte
    Q = _mm_srli_si128(Q, 1);
    Q = _mm_or_si128(Q, last_inf);
    COMPUTE_SCORE(S, Sx, Sy);
    Sy = _mm_slli_si128(Sy, 1);
    __m128i msk =
        _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, Y1[ix + 15]);
    Sy = _mm_or_si128(Sy, msk);
    R = _mm_adds_epi8(R2, S);
    if (ix + 3 == m) {
      flag |= 2;
      gap_open1 = _mm_srli_si128(gap_open1, 1);
      gap_ext1 = _mm_srli_si128(gap_ext1, 1);
    }
    if (flag & 2) {
      nx++;
      gap_open1 = _mm_srli_si128(gap_open1, 1);
      gap_ext1 = _mm_srli_si128(gap_ext1, 1);
      __m128i tt = _mm_cmpeq_epi8(edge_mask2, zero);
      R = _mm_and_si128(R, tt);
      R = _mm_or_si128(R, edge_mask2);
      P = _mm_and_si128(P, tt);
      P = _mm_or_si128(P, edge_mask2);
      Q = _mm_and_si128(Q, tt);
      Q = _mm_or_si128(Q, edge_mask2);
    }
    if (flag & 1) {
      __m128i tt = _mm_cmpeq_epi8(edge_mask1, zero);
      R = _mm_and_si128(R, tt);
      R = _mm_or_si128(R, edge_mask1);
      P = _mm_and_si128(P, tt);
      P = _mm_or_si128(P, edge_mask1);
      Q = _mm_and_si128(Q, tt);
      Q = _mm_or_si128(Q, edge_mask1);
      edge_mask1 = _mm_slli_si128(edge_mask1, 1);
      edge_mask1 = _mm_or_si128(edge_mask1, first_inf);
    }

    R2 = R1;
    if (ix < 12) {
      R = _mm_and_si128(R, edge_mask);
      UPDATE_R(R, P, Q, R1, bts);
      bts = _mm_and_si128(bts, edge_mask);
    } else
      UPDATE_R(R, P, Q, R1, bts);

    UPDATE_PQ1(R1, P, Q, bts, gap_open1, gap_ext1);
    _mm_store_si128(bitp++, bts);

    if (nx + ny == 15) {
      flag |= 4;
      break;
    }
    flag |= (ix + 12 == n);
    // Shift P vector left 1 byte
    P = _mm_slli_si128(P, 1);
    P = _mm_or_si128(P, first_inf);
    // Calculate score vector
    COMPUTE_SCORE(S, Sx, Sy);
    Sx = _mm_srli_si128(Sx, 1);
    msk =
        _mm_set_epi8(X1[ix + 16], 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    Sx = _mm_or_si128(Sx, msk);
    R = _mm_adds_epi8(R2, S);
    if (flag & 1) {
      ny++;
      __m128i tt = _mm_cmpeq_epi8(edge_mask1, zero);
      R = _mm_and_si128(R, tt);
      R = _mm_or_si128(R, edge_mask1);
      P = _mm_and_si128(P, tt);
      P = _mm_or_si128(P, edge_mask1);
      Q = _mm_and_si128(Q, tt);
      Q = _mm_or_si128(Q, edge_mask1);
    }
    if (flag & 2) {
      __m128i tt = _mm_cmpeq_epi8(edge_mask2, zero);
      R = _mm_and_si128(R, tt);
      R = _mm_or_si128(R, edge_mask2);
      P = _mm_and_si128(P, tt);
      P = _mm_or_si128(P, edge_mask2);
      Q = _mm_and_si128(Q, tt);
      Q = _mm_or_si128(Q, edge_mask2);
      edge_mask2 = _mm_srli_si128(edge_mask2, 1);
      edge_mask2 = _mm_or_si128(edge_mask2, last_inf);
    }

    R2 = R1;
    if (ix < 12) {
      R = _mm_and_si128(R, edge_mask);
      UPDATE_R(R, P, Q, R1, bts);
      bts = _mm_and_si128(bts, edge_mask);
      edge_mask = _mm_srli_si128(edge_mask, 1);
      edge_mask = _mm_or_si128(edge_mask, last_ones);
    } else
      UPDATE_R(R, P, Q, R1, bts);

    UPDATE_PQ1(R1, P, Q, bts, gap_open1, gap_ext1);
    _mm_store_si128(bitp++, bts);
    if (nx + ny == 15) {
      flag |= 8;
      break;
    }
  }

  // Reconstruct edges
  __m128i B, B1, B2, B1a;

  // First element
  B2 = zero;
  B1 = _mm_load_si128(--bitp);
  __m128i tt = _mm_cmpeq_epi8(B1, msk_0x7f);
  tt = _mm_andnot_si128(tt, msk_0x07);
  int cx, cy, cx1, cy1, cy2;
  int state = 0;
  cx = cx1 = m-1;
  cy2 = m + 15;

  if (flag & 4) {
    B1a = _mm_and_si128(B1, tt);
    _mm_store_si128(bits + 1, B1a);
    B1 = _mm_slli_si128(B1a, 1);
    cy = cy1 = m+13+x-x1;
  } else {
    cy = cy1 = m+14+x-x1;
    B1 = _mm_and_si128(B1, tt);
    _mm_store_si128(bits + 1, B1);
    B1a = _mm_srli_si128(B1, 1);
    B = _mm_load_si128(--bitp);
    UPDATE_EDGES(B, B1, B1a, B2);
    int k = _mm_testz_si128(B1, msk_b);
    if (k) {
      printf("Starting at %d,%d\n", m + 11, cy2);
    }
    _mm_store_si128(bits + 1, B1);
    int rt;
    CHECK_STATE(B1, rt);
    state = update_state(&B1, state, rt, cx, cy, &cx1, &cy1, X, Y);
    B2 = B1;
    B1a = B;
    B1 = _mm_slli_si128(B1a, 1);
    cy--;
    cy2--;
  }
  ix--;
  for (; ix >= 0; ix--) {
    B = _mm_load_si128(--bitp);
    UPDATE_EDGES1(B, B1, B1a, B2);
    _mm_store_si128(bits + 1, B1);
    if (!state) {
      int k = _mm_testz_si128(B1, msk_b);
      if (k) {
        printf("Starting at %d,%d\n", m + 11, cy2);
      }
    }
    int rt;
    CHECK_STATE(B1, rt);
    state = update_state(&B1, state, rt, cx, cy, &cx1, &cy1, X, Y);
    cx--;
    B2 = B1;
    B1 = B;
    B1a = _mm_srli_si128(B1, 1);
    B = _mm_load_si128(--bitp);
    UPDATE_EDGES(B, B1, B1a, B2);
    _mm_store_si128(bits + 1, B1);
    if (!state) {
      int k = _mm_testz_si128(B1, msk_b);
      if (k) {
        printf("Starting at %d,%d\n", m + 11, cy2);
      }
    }
    CHECK_STATE(B1, rt);
    state = update_state(&B1, state, rt, cx, cy, &cx1, &cy1, X, Y);
    B2 = B1;
    B1a = B;
    B1 = _mm_slli_si128(B1a, 1);
    cy--;
    cy2--;
  }
  return 0;
}
