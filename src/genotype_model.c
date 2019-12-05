#include <stdio.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_multimin.h>

#include "gem_tools.h"
#include "bs_call.h"

static qual_prob q_prob[MAX_QUAL + 1];

void fill_base_prob_table(void) {
  for (int q = 0; q <= MAX_QUAL; q++) {
    double e = exp(-.1 * (double)q * LOG10);
    if(e > .5) e = .5;
    double k = e / (3.0 - 4.0 * e);
    q_prob[q].e = e;
    q_prob[q].k = k;
    q_prob[q].ln_k = log(k);
    q_prob[q].ln_k_half = log(0.5 + k);
    q_prob[q].ln_k_one = log(1.0 + k);
  }
}

static inline void get_Z(double x1, double x2, double k1, double k2, double l, double t, double *Z) {
  double lpt = l + t;
  double lmt = l - t;
  double d = (x1 + x2) * lmt;
  // w = 1, p = 1
  double sinm = (x1 * (lpt + 2.0 * k2) - x2 * (2.0 - lpt + 2.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[0] = 0.5 * (lmt * sinm + 2.0 - lpt);
  // w = 1, p = 1/2
  sinm = (x1 * (2.0 + lpt + 4.0 * k2) - x2 * (2.0 - lpt + 4.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[1] = 0.5 * (lmt * sinm + 2.0 - lpt);
  // w = 1/2, p = 1
  sinm = (x1 * (lpt + 4.0 * k2) - x2 * (2.0 - lpt + 4.0 * k1)) / d;
  if(sinm < -1.0) sinm = -1.0;
  else if(sinm > 1.0) sinm = 1.0;
  Z[2] = 0.5 * (lmt * sinm + 2.0 - lpt);
}

typedef struct {
  double l,t;

  double bias;  
  char ref;
  double k[8];
  double n[8];
} like_par;

double calc_full_model_like(const gsl_vector *v, void *par_p) {
  like_par *par = par_p;			
  double l = par->l;				
  double t = par->t;				
  double *n = par->n;				
  double *k = par->k;
  double lk = 0.0;
  double w1 = gsl_vector_get(v, 0);
  double w = 0.5 * (1.0 + sin(w1));
  double p1 = gsl_vector_get(v, 1);
  double p = 0.5 * (1.0 + sin(p1));
  double q1 = gsl_vector_get(v, 2);
  double q = 0.5 * (1.0 + sin(q1));
  double bias = par->bias;
  switch(par->ref) {
  case 'A' :
    lk += log(1.0 + (1.0 - w) * (1.0 - q) * (bias - 1));
    break;
  case 'C':
    lk += log(1.0 + w * p * (bias - 1));
    break;
  case 'G':
    lk += log(1.0 + (1.0 - w) * q * (bias - 1));
    break;
  case 'T':
    lk += log(1.0 + w * (1.0 - p) * (bias - 1));
    break;
  }
  if(n[0] + n[2] + n[4] + n[6] > 0.0) {
    lk += n[0] * log((1.0 - w) * (1.0 - q) + k[0]) + n[2] * log((1.0 - w) * q + k[2]);
    if(n[4] + n[6] > 0) {
      double mg1 = gsl_vector_get(v, 4);
      double zg = 0.5 * (1.0 + sin(mg1)) * (l - t) + 1.0 - l;
      lk += n[4] * log((1.0 - w) * (1.0 - q * zg) + k[4]) + n[6] * log((1.0 - w) * q * zg + k[6]);
    }
  }
  if(n[1] + n[3] + n[5] + n[7] > 0.0) {
    lk += n[1] * log(w * p + k[1]) + n[3] * log(w * (1.0 - p) + k[3]);
    if(n[5] + n[7] > 0) {
      double mc1 = gsl_vector_get(v, 3);
      double zc = 0.5 * (1.0 + sin(mc1)) * (l - t) + 1.0 - l;
      lk += n[5] * log(w * p * zc + k[5]) + n[7] * log(w * (1.0 - p * zc) + k[7]);
    }
  }
  /*  fprintf(stderr, "Like():");
  for(int j = 0; j < 5; j++) {
    double z = gsl_vector_get(v, j);
    double z1 = (1.0 + sin(z)) / 2.0;
    fprintf(stderr, "\t(%g, %g)", z, z1);
  }
  fprintf(stderr, " -> %.15g\n",lk); */
  return -lk;
}

void calc_full_model_like_gradient(const gsl_vector *v, void *par_p, double *like, gsl_vector *df) {
  like_par *par = par_p;			
  double l = par->l;				
  double t = par->t;				
  double *n = par->n;				
  double *k = par->k;
  double lk = 0.0;
  double sin_v[5], cos_v[5];
  double w1 = gsl_vector_get(v, 0);
  sin_v[0] = sin(w1);
  cos_v[0] = cos(w1);
  double w = 0.5 * (1.0 + sin_v[0]);
  double g[] = {0, 0, 0, 0, 0};
  double p1 = gsl_vector_get(v, 1);
  sin_v[1] = sin(p1);
  cos_v[1] = cos(p1);
  double p = 0.5 * (1.0 + sin_v[1]);
  double q1 = gsl_vector_get(v, 2);
  sin_v[2] = sin(q1);
  cos_v[2] = cos(q1);
  double q = 0.5 * (1.0 + sin_v[2]);
  double bias = par->bias;
  double z;
  switch(par->ref) {
  case 'A' :
    z = 1.0 + (1.0 - w) * (1.0 - q) * (bias - 1);
    g[0] += -(cos_v[0] * (1.0 - q) * (bias - 1.0)) / (2.0 * z);
    g[2] += -(cos_v[2] * (1.0 - w) * (bias - 1.0)) / (2.0 * z);
    break;
  case 'C':
    z = 1.0 + w * p * (bias - 1);
    g[0] += (cos_v[0] * p * (bias - 1.0)) / (2.0 * z);
    g[1] += (cos_v[1] * w * (bias - 1.0)) / (2.0 * z);
    break;
  case 'G':
    z = 1.0 + (1.0 - w) * q * (bias - 1);
    g[0] += -(cos_v[0] * q * (bias - 1.0)) / (2.0 * z);
    g[2] += (cos_v[2] * (1.0 - w) * (bias - 1.0)) / (2.0 * z);
    break;
  case 'T':
    z = 1.0 + w * (1.0 - p) * (bias - 1);
    g[0] += (cos_v[0] * (1.0 - p) * (bias - 1.0)) / (2.0 * z);
    g[1] += -(cos_v[1] * w * (bias - 1.0)) / (2.0 * z);
    break;
  default:
    z = 1.0;
    break;
  }  
  lk += log(z);
  if(n[0] + n[2] + n[4] + n[6] > 0.0) {
    if(n[0] > 0.0) {
      lk += n[0] * log((1.0 - w) * (1.0 - q) + k[0]);
      g[0] += -n[0] * (1.0 - q) * cos_v[0] / ((1.0 - sin_v[0]) * (1.0 - q) + 2.0 * k[0]);
      g[2] += -n[0] * (1.0 - w) * cos_v[2] / ((1.0 - sin_v[2]) * (1.0 - w) + 2.0 * k[0]);
    }
    if(n[2] > 0.0) {
      lk += n[2] * log((1.0 - w) * q + k[2]);
      g[0] += -n[2] * q * cos_v[0] / ((1.0 - sin_v[0]) * q + 2.0 * k[2]);
      g[2] += n[2] * (1.0 - w) * cos_v[2] / ((1.0 + sin_v[2]) * (1.0 - w) + 2.0 * k[2]);
    }
    if(n[4] + n[6] > 0.0) {
      double mg1 = gsl_vector_get(v, 4);
      sin_v[4] = sin(mg1);
      cos_v[4] = cos(mg1);
      double zg = 0.5 * (1.0 + sin_v[4]) * (l - t) + 1.0 - l;
      if(n[4] > 0.0) {
	lk += n[4] * log((1.0 - w) * (1.0 - q * zg) + k[4]);
	g[0] += -n[4] * (1.0 - q * zg) * cos_v[0] / ((1.0 - sin_v[0]) * (1.0 - q * zg) + 2.0 * k[4]);
	g[2] += -n[4] * (1.0 - w) * zg * cos_v[2]  / ((1.0 - w) * (2.0 - zg * (sin_v[2] + 1.0)) + 2 * k[4]);
	g[4] += -n[4] * (1.0 - w) * q * (l - t) * cos_v[4] / (2 * (k[4] + (1.0 - w) * (1.0 - q * zg)));
      }
      if(n[6] > 0.0) {
	lk += n[6] * log((1.0 - w) * q * zg + k[6]);
	g[0] += -n[6] * q * zg * cos_v[0] / ((1.0 - sin_v[0]) * q * zg + 2.0 * k[6]);
	g[2] += n[6] * (1.0 - w) * zg * cos_v[2] / ((1.0 + sin_v[2]) * (1.0 - w) * zg + 2.0 * k[6]);
	g[4] += n[6] * (1.0 - w) * q * (l - t) * cos_v[4] / (2 * (k[6] + (1.0 - w) * q * zg));
      }			
    }
  }
  if(n[1] + n[3] + n[5] + n[7] > 0.0) {
    if(n[1] > 0.0) {
      lk += n[1] * log(w * p + k[1]);
      g[0] += n[1] * p * cos_v[0] / ((1.0 + sin_v[0]) * p + 2.0 * k[1]);
      g[1] += n[1] * w * cos_v[1] / ((1.0 + sin_v[1]) * w + 2.0 * k[1]);
    }
    if(n[3] > 0.0) {
      lk += n[3] * log(w * (1.0 - p) + k[3]);
      g[0] += n[3] * (1.0 - p) * cos_v[0] / ((1.0 + sin_v[0]) * (1.0 - p) + 2.0 * k[3]);
      g[1] += -n[3] * w * cos_v[1] / ((1.0 - sin_v[1]) * w + 2.0 * k[3]);
    }
    if(n[5] + n[7] > 0.0) {
      double mc1 = gsl_vector_get(v, 3);
      sin_v[3] = sin(mc1);
      cos_v[3] = cos(mc1);
      double zc = 0.5 * (1.0 + sin_v[3]) * (l - t) + 1.0 - l;
      if(n[5] > 0.0) {
	lk += n[5] * log(w * p * zc + k[5]);
	g[0] += n[5] * p * zc * cos_v[0] / ((1.0 + sin_v[0]) * p * zc + 2.0 * k[5]);
	g[1] += n[5] * w * zc * cos_v[1] / ((1.0 + sin_v[1]) * w * zc + 2.0 * k[5]);
	g[3] += n[5] * w * p * (l - t) * cos_v[3] / (2 * (k[5] + w * p * zc));
      }
      if(n[7] > 0.0) {
	lk += n[7] * log(w * (1.0 - p * zc) + k[7]);
	g[0] += n[7] * (1.0 - p * zc) * cos_v[0] / ((1.0 + sin_v[0]) * (1.0 - p * zc) + 2.0 * k[7]);
	g[1] += -n[7] * w * zc * cos_v[1] / (w * (2.0 - zc * (sin_v[1] + 1.0)) + 2 * k[7]);
	g[3] += -n[7] * w * p * (l - t) * cos_v[3] / (2 * (k[7] + w * (1.0 - p * zc)));
      }
    }
  }
  for(int i = 0; i < 5; i++) gsl_vector_set(df, i, -g[i]);
  /*  fprintf(stderr, "Like_Gradient():");
  for(int j = 0; j < 5; j++) {
    double z = gsl_vector_get(v, j);
    double z1 = (1.0 + sin(z)) / 2.0;
    fprintf(stderr, "\t(%g, %g)", z, z1);
  }
  fputc('\n',stderr);
  fprintf(stderr, "   -> %.15g",lk);
  for(int j = 0; j < 5; j++) fprintf(stderr,"\t%g", -g[j]);
  fputc('\n',stderr); */
  *like = -lk;
}

void calc_full_model_gradient(const gsl_vector *v, void *par_p, gsl_vector *df) {
  like_par *par = par_p;			
  double l = par->l;				
  double t = par->t;				
  double *n = par->n;				
  double *k = par->k;
  double sin_v[5], cos_v[5];
  double w1 = gsl_vector_get(v, 0);
  sin_v[0] = sin(w1);
  cos_v[0] = cos(w1);
  double w = 0.5 * (1.0 + sin_v[0]);
  double g[] = {0, 0, 0, 0, 0};
  double p1 = gsl_vector_get(v, 1);
  sin_v[1] = sin(p1);
  cos_v[1] = cos(p1);
  double p = 0.5 * (1.0 + sin_v[1]);
  double q1 = gsl_vector_get(v, 2);
  sin_v[2] = sin(q1);
  cos_v[2] = cos(q1);
  double q = 0.5 * (1.0 + sin_v[2]);
  double bias = par->bias;
  double z;
  switch(par->ref) {
  case 'A' :
    z = 1.0 + (1.0 - w) * (1.0 - q) * (bias - 1);
    g[0] += -(cos_v[0] * (1.0 - q) * (bias - 1.0)) / (2.0 * z);
    g[2] += -(cos_v[2] * (1.0 - w) * (bias - 1.0)) / (2.0 * z);
    break;
  case 'C':
    z = 1.0 + w * p * (bias - 1);
    g[0] += (cos_v[0] * p * (bias - 1.0)) / (2.0 * z);
    g[1] += (cos_v[1] * w * (bias - 1.0)) / (2.0 * z);
    break;
  case 'G':
    z = 1.0 + (1.0 - w) * q * (bias - 1);
    g[0] += -(cos_v[0] * q * (bias - 1.0)) / (2.0 * z);
    g[2] += (cos_v[2] * (1.0 - w) * (bias - 1.0)) / (2.0 * z);
    break;
  case 'T':
    z = 1.0 + w * (1.0 - p) * (bias - 1);
    g[0] += (cos_v[0] * (1.0 - p) * (bias - 1.0)) / (2.0 * z);
    g[1] += -(cos_v[1] * w * (bias - 1.0)) / (2.0 * z);
    break;
  default:
    z = 1.0;
    break;
  }  
  if(n[0] + n[2] + n[4] + n[6] > 0.0) {
    if(n[0] > 0.0) {
      g[0] += -n[0] * (1.0 - q) * cos_v[0] / ((1.0 - sin_v[0]) * (1.0 - q) + 2.0 * k[0]);
      g[2] += -n[0] * (1.0 - w) * cos_v[2] / ((1.0 - sin_v[2]) * (1.0 - w) + 2.0 * k[0]);
    }
    if(n[2] > 0.0) {
      g[0] += -n[2] * q * cos_v[0] / ((1.0 - sin_v[0]) * q + 2.0 * k[2]);
      g[2] += n[2] * (1.0 - w) * cos_v[2] / ((1.0 + sin_v[2]) * (1.0 - w) + 2.0 * k[2]);
    }
    if(n[4] + n[6] > 0.0) {
      double mg1 = gsl_vector_get(v, 4);
      sin_v[4] = sin(mg1);
      cos_v[4] = cos(mg1);
      double zg = 0.5 * (1.0 + sin_v[4]) * (l - t) + 1.0 - l;
      if(n[4] > 0.0) {
	g[0] += -n[4] * (1.0 - q * zg) * cos_v[0] / ((1.0 - sin_v[0]) * (1.0 - q * zg) + 2.0 * k[4]);
	g[2] += -n[4] * (1.0 - w) * zg * cos_v[2]  / ((1.0 - w) * (2.0 - zg * (sin_v[2] + 1.0)) + 2 * k[4]);
	g[4] += -n[4] * (1.0 - w) * q * (l - t) * cos_v[4] / (2 * (k[4] + (1.0 - w) * (1.0 - q * zg)));
      }
      if(n[6] > 0.0) {
	g[0] += -n[6] * q * zg * cos_v[0] / ((1.0 - sin_v[0]) * q * zg + 2.0 * k[6]);
	g[2] += n[6] * (1.0 - w) * zg * cos_v[2] / ((1.0 + sin_v[2]) * (1.0 - w) * zg + 2.0 * k[6]);
	g[4] += n[6] * (1.0 - w) * q * (l - t) * cos_v[4] / (2 * (k[6] + (1.0 - w) * q * zg));
      }			
    }
  }
  if(n[1] + n[3] + n[5] + n[7] > 0.0) {
    if(n[1] > 0.0) {
      g[0] += n[1] * p * cos_v[0] / ((1.0 + sin_v[0]) * p + 2.0 * k[1]);
      g[1] += n[1] * w * cos_v[1] / ((1.0 + sin_v[1]) * w + 2.0 * k[1]);
    }
    if(n[3] > 0.0) {
      g[0] += n[3] * (1.0 - p) * cos_v[0] / ((1.0 + sin_v[0]) * (1.0 - p) + 2.0 * k[3]);
      g[1] += -n[3] * w * cos_v[1] / ((1.0 - sin_v[1]) * w + 2.0 * k[3]);
    }
    if(n[5] + n[7] > 0.0) {
      double mc1 = gsl_vector_get(v, 3);
      sin_v[3] = sin(mc1);
      cos_v[3] = cos(mc1);
      double zc = 0.5 * (1.0 + sin_v[3]) * (l - t) + 1.0 - l;
      if(n[5] > 0.0) {
	g[0] += n[5] * p * zc * cos_v[0] / ((1.0 + sin_v[0]) * p * zc + 2.0 * k[5]);
	g[1] += n[5] * w * zc * cos_v[1] / ((1.0 + sin_v[1]) * w * zc + 2.0 * k[5]);
	g[3] += n[5] * w * p * (l - t) * cos_v[3] / (2 * (k[5] + w * p * zc));
      }
      if(n[7] > 0.0) {
	g[0] += n[7] * (1.0 - p * zc) * cos_v[0] / ((1.0 + sin_v[0]) * (1.0 - p * zc) + 2.0 * k[7]);
	g[1] += -n[7] * w * zc * cos_v[1] / (w * (2.0 - zc * (sin_v[1] + 1.0)) + 2 * k[7]);
	g[3] += -n[7] * w * p * (l - t) * cos_v[3] / (2 * (k[7] + w * (1.0 - p * zc)));
      }
    }
  }
  for(int i = 0; i < 5; i++) gsl_vector_set(df, i, -g[i]);
  /*  fprintf(stderr, "Gradient():");
  for(int j = 0; j < 5; j++) {
    double z = gsl_vector_get(v, j);
    double z1 = (1.0 + sin(z)) / 2.0;
    fprintf(stderr, "\t(%g, %g)", z, z1);
  }
  fputc('\n',stderr);
  for(int j = 0; j < 5; j++) fprintf(stderr, "\t%g", -g[j]);
  fputc('\n',stderr); */
}

double maximize_full_model(gt_meth *gt, qual_prob *qp, double l, double t, double bias, char rf, double *beta) {
  like_par par = {
		  .l = l,
		  .t = t,
		  .bias = bias,
		  .ref = rf
  };
  gsl_multimin_function_fdf func = {
				    .n = 5,
				    .f = calc_full_model_like,
				    .df = calc_full_model_gradient,
				    .fdf = calc_full_model_like_gradient,
				    .params = &par
  };
  for(int i = 0; i < 8; i++) {
    par.n[i] = (double)gt->counts[i];
    par.k[i] = qp[i].k;
  }
  gsl_vector *x = gsl_vector_alloc(5);
  const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_conjugate_fr;
  // const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
  gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc (T, 5);
  // Initially set all estimates to 0.0 (corresponds to a ratio of 0.5)
  for(int i = 0; i < 5; i++) gsl_vector_set(x, i, 0.0);
  // Where we have information we can make some guesses to improve the starting position
  double nn = 0.0;
  for(int i = 0; i < 8; i++) nn += par.n[i];
  double z = (par.n[1] + par.n[3] + par.n[5] + par.n[7]) / nn;
  gsl_vector_set(x, 0, asin(2.0 * z - 1.0));
  bool est_p = false, est_q = false;
  double p = 0.0, q = 0.0;
  if(par.n[1] + par.n[3] > 0.0) {
    p = par.n[1] / (par.n[1] + par.n[3]);
    gsl_vector_set(x, 1, asin(2.0 * p - 1.0));
    est_p = true;
  }
  if(par.n[0] + par.n[2] > 0.0) {
    q = par.n[2] / (par.n[0] + par.n[2]);
    gsl_vector_set(x, 2, asin(2.0 * q - 1.0));		
    est_q = true;
  }
  if(par.n[5] + par.n[7] > 0.0) {
    z = par.n[5] / (par.n[5] + par.n[7]);
    if(est_p && p > 0.0) {
      z /= p;
      if(z > 1.0) z = 1.0;
    }
    gsl_vector_set(x, 3, asin(2.0 * z - 1.0));
  }
  if(par.n[4] + par.n[6] > 0.0) {
    z = par.n[6] / (par.n[4] + par.n[6]);		
    if(est_q && q > 0.0) {
      z /= q;
      if(z > 1.0) z = 1.0;
    }
    gsl_vector_set(x, 4, asin(2.0 * z - 1.0));		
  }
  gsl_multimin_fdfminimizer_set (s, &func, x, 0.1, 0.1);
  int status = 0;
  for(int it = 0; it < 100; it++) {
    status = gsl_multimin_fdfminimizer_iterate(s);
    status = gsl_multimin_test_gradient (s->gradient, 1e-4);
    if (status == GSL_SUCCESS) break;
  }
  for(int i = 0; i < 5; i++) beta[i] = 0.5 * (1.0 + sin(gsl_vector_get(s->x, i)));
  double lk = -s->f;
  /*  fprintf(stderr,"FM:(");
  for(int i = 0; i < 8; i++) fprintf(stderr," %g", par.n[i]);
  fputc(')', stderr);
  for(int i = 0; i < 5; i++) fprintf(stderr," %g", beta[i]);
  fprintf(stderr,"\t(%c, %g)\t%g\n", rf, bias, lk); */
  gsl_multimin_fdfminimizer_free(s);
  gsl_vector_free(x);
  return lk;
}

void calc_gt_prob(gt_meth *gt, const sr_param * const param, char rf) {
  qual_prob qp[8];
  for(int i = 0; i < 8; i++) qp[i] = q_prob[gt->qual[i]];
  double l = 1.0 - param->under_conv;
  double t = param->over_conv;

  /***********************************************************************************
   * Base and methylation frequencies are described by 5 parameters: w, p, q, mc, mg
   * 
   * Let n(X) be the count for base X, and N the total number of bases seen
   * w = (n(C) + n(T)) / N
   * p = n(C) / (n(C) + n(T))
   * q = n(G) / (n(A) + n(G))
   * mc is the proportion of methylated Cs on the top strand
   * mg is the proportion of methylated Cs on the bottom strand
   *
   * Base frequencies are therefore:
   *  f(A) = (1 - w) * (1 - q)
   *  f(C) = w * p
   *  f(G) = (1 - w) * q
   *  f(T) = w * (1 - p)
   *
   * All 5 parameters are ratios are are therefore independently constrained 
   * to be between 0 and 1.
   * 
   * We first maximize the full model, allowing w, p, q, mc and mg to take 
   * any legal value.  The log likelihood of this model is l_full.
   *
   * We then calculate the marginal likelihood for the 10 possible genetic models compatible
   * with a diploid state (thereby fixing w, p, q) and maximizing the likelihood over (mc, mg).
   *
   * The called genotype is that with the highest likelihood 
   * The phred score is calculated as the phred scaled posterior genotype probability (considering
   * only the 10 possible diploid genotypes)
   * The goodness of fit score is the phred scaled likelihood ratio between l_full and the likelihood
   * of the called genotype.
   *
   ***********************************************************************************/

  double beta[5];
  double ref_bias = param->ref_bias;
  double full_lk = maximize_full_model(gt, qp, l, t, ref_bias, rf, beta);
  double n[8];
  for (int i = 0; i < 8; i++) n[i] = (double)gt->counts[i];
  double ll[10];
  memset(ll, 0, sizeof(double) * 10);
  // Add in prior from reference
  {
    double lrb = log(ref_bias);
    double lrb1 = log(0.5 * (1.0 + ref_bias));
    switch (rf) {
    case 'A':
      ll[0] = lrb;
      ll[1] = ll[2] = ll[3] = lrb1;
      break;
    case 'C':
      ll[4] = lrb;
      ll[1] = ll[5] = ll[6] = lrb1;
      break;
    case 'G':
      ll[7] = lrb;
      ll[2] = ll[5] = ll[8] = lrb1;
      break;
    case 'T':
      ll[9] = lrb;
      ll[3] = ll[6] = ll[8] = lrb1;
      break;
    }
  }
  if (n[0]) {
    ll[0] += n[0] * qp[0].ln_k_one; // AA
    double tz = n[0] * qp[0].ln_k_half; 
    ll[1] += tz; // AC
    ll[2] += tz; // AG
    ll[3] += tz; // AT
    tz = n[0] * qp[0].ln_k;
    ll[4] += tz;  // CC
    ll[5] += tz;  // CG
    ll[6] += tz;  // CT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[1]) {
    ll[4] += n[1] * qp[1].ln_k_one; // CC
    double tz = n[1] * qp[1].ln_k_half;
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    ll[6] += tz; // CT
    tz = n[1] * qp[1].ln_k;
    ll[0] += tz;  // AA
    ll[2] += tz;  // AG
    ll[3] += tz;  // AT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[2]) {
    ll[7] += n[2] * qp[2].ln_k_one; // GG
    double tz = n[2] * qp[2].ln_k_half;
    ll[2] += tz; // AG
    ll[5] += tz; // CG
    ll[8] += tz; // TG
    tz = n[2] * qp[2].ln_k;
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[3]) {
    ll[9] += n[3] * qp[3].ln_k_one; // TT
    double tz = n[3] * qp[3].ln_k_half;
    ll[3] += tz; // AT
    ll[6] += tz; // CT
    ll[8] += tz; // GT
    tz = n[3] * qp[3].ln_k;
    ll[0] += tz; // AA
    ll[1] += tz; // AC
    ll[2] += tz; // AG
    ll[4] += tz; // CC
    ll[5] += tz; // CG
    ll[7] += tz; // GG
  }
  double Z[6] = {-1.0, -1.0, -1.0, -1.0, -1.0, -1.0};
  if (n[5] + n[7] > 0.0) {
    get_Z(n[5], n[7], qp[5].k, qp[7].k, l, t, Z);
  }
  if (n[4] + n[6] > 0.0) {
    get_Z(n[6], n[4], qp[6].k, qp[4].k, l, t, Z+3);
  }

  if (n[4]) {
    ll[0] += n[4] * qp[4].ln_k_one; // AA
    ll[2] += log(1.0 - 0.5 * Z[4] + qp[4].k) * n[4]; // AG
    ll[7] += log(1.0 - Z[3] + qp[4].k) * n[4]; // GG
    double tz = log(0.5 * (1.0 - Z[5]) + qp[4].k) * n[4];
    ll[5] += tz; // CG
    ll[8] += tz; // GT
    tz = n[4] * qp[4].ln_k_half;
    ll[1] += tz; // AC
    ll[3] += tz; // AT
    tz = n[4] * qp[4].ln_k;
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[5]) {
    ll[4] += log(Z[0] + qp[5].k) * n[5]; // CC
    double tz = log(0.5 * Z[2] + qp[5].k) * n[5];
    ll[1] += tz;                              // AC
    ll[5] += tz;                              // CG
    ll[6] += log(0.5 * Z[1] + qp[5].k) * n[5]; // CT
    tz = n[5] * qp[5].ln_k;
    ll[0] += tz;  // AA
    ll[2] += tz;  // AG
    ll[3] += tz;  // AT
    ll[7] += tz;  // GG
    ll[8] += tz;  // GT
    ll[9] += tz; // TT
  }
  if (n[6]) {
    ll[7] += log(Z[3] + qp[6].k) * n[6]; // GG
    double tz = log(0.5 * Z[5] + qp[6].k) * n[6];
    ll[5] += tz; // CG
    ll[8] += tz; // TG
    ll[2] += log(0.5 * Z[4] + qp[6].k) * n[6]; // AG
    tz = n[6] * qp[6].ln_k;
    ll[0] += tz;  // AA
    ll[1] += tz;  // AC
    ll[3] += tz;  // AT
    ll[4] += tz;  // CC
    ll[6] += tz;  // CT
    ll[9] += tz; // TT
  }
  if (n[7]) {
    ll[9] += n[7] * qp[7].ln_k_one;  // TT
    ll[4] += log(1.0 - Z[0] + qp[7].k) * n[7]; // CC
    ll[6] += log(1.0 - 0.5 * Z[1] + qp[7].k) * n[7]; // CT
    double tz = log(0.5 * (1.0 - Z[2]) + qp[7].k) * n[7];
    ll[1] += tz; // AC
    ll[5] += tz; // CG
    tz = n[7] * qp[7].ln_k_half;
    ll[3] += tz; // AT
    ll[8] += tz; // GT
    tz = n[7] * qp[7].ln_k;
    ll[0] += tz; // AA
    ll[2] += tz; // AG
    ll[7] += tz; // GG
  }
  double max = ll[0];
  int mx = 0;
  for (int i = 1; i < 10; i++) {
    if (ll[i] > max) {
      max = ll[i];
      mx = i;
    }
  }
  gt->max_gt = mx;
  double sum = 0.0;
  for (int i = 0; i < 10; i++) sum += exp(ll[i] - max);
  sum = log(sum);
  for (int i = 0; i < 10; i++) {
    gt->gt_prob[i] = (ll[i] - max - sum) / LOG10;
  }
  gt->gt_gof = full_lk > max ? (full_lk - max) / LOG10 : 0.0;
}


