/* ----------------------------------------------- */
/*        The MGC Vocoder (part of MGC_engine)     */
/*        developed by HTS Working Group           */
/*        maintained by Xingyu Na                  */
/* ----------------------------------------------- */
/* Copyright (c)                                   */
/*   2011-2012  Beijing Institute of Technology    */
/*              RCDCT, @RICT                       */
/*   2001-2011  Nagoya Institute of Technology     */
/*              Department of Computer Science     */
/* ----------------------------------------------- */

#include <stdlib.h>
#include <math.h>               /* for sqrt(),log(),exp(),pow(),cos() */

#include "Synthesizer.h"

static const double MGC_pade[21] = {
   1.00000000000,
   1.00000000000,
   0.00000000000,
   1.00000000000,
   0.00000000000,
   0.00000000000,
   1.00000000000,
   0.00000000000,
   0.00000000000,
   0.00000000000,
   1.00000000000,
   0.49992730000,
   0.10670050000,
   0.01170221000,
   0.00056562790,
   1.00000000000,
   0.49993910000,
   0.11070980000,
   0.01369984000,
   0.00095648530,
   0.00003041721
};

/* MGC_movem: move memory */
static void MGC_movem(double *a, double *b, const int nitem)
{
   long i = (long) nitem;

   if (a > b)
      while (i--)
         *b++ = *a++;
   else {
      a += i;
      b += i;
      while (i--)
         *--b = *--a;
   }
}

/* MGC_mlsafir: sub functions for MLSA filter */
static double MGC_mlsafir(const double x, const double *b, const int m, const double a, const double aa, double *d)
{
   double y = 0.0;
   int i;

   d[0] = x;
   d[1] = aa * d[0] + a * d[1];

   for (i = 2; i <= m; i++)
      d[i] += a * (d[i + 1] - d[i - 1]);

   for (i = 2; i <= m; i++)
      y += d[i] * b[i];

   for (i = m + 1; i > 1; i--)
      d[i] = d[i - 1];

   return (y);
}

/* MGC_mlsadf1: sub functions for MLSA filter */
static double MGC_mlsadf1(double x, const double *b, const int m, const double a, const double aa, const int pd, double *d, const double *ppade)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd + 1];

   for (i = pd; i >= 1; i--) {
      d[i] = aa * pt[i - 1] + a * d[i];
      pt[i] = d[i] * b[1];
      v = pt[i] * ppade[i];
      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

/* MGC_mlsadf2: sub functions for MLSA filter */
static double MGC_mlsadf2(double x, const double *b, const int m, const double a, const double aa, const int pd, double *d, const double *ppade)
{
   double v, out = 0.0, *pt;
   int i;

   pt = &d[pd * (m + 2)];

   for (i = pd; i >= 1; i--) {
      pt[i] = MGC_mlsafir(pt[i - 1], b, m, a, aa, &d[(i - 1) * (m + 2)]);
      v = pt[i] * ppade[i];

      x += (1 & i) ? v : -v;
      out += v;
   }

   pt[0] = x;
   out += x;

   return (out);
}

/* MGC_mlsadf: functions for MLSA filter */
static double MGC_mlsadf(double x, const double *b, const int m, const double a, const int pd, double *d)
{
   const double aa = 1 - a * a;
   const double *ppade = &(MGC_pade[pd * (pd + 1) / 2]);

   x = MGC_mlsadf1(x, b, m, a, aa, pd, d, ppade);
   x = MGC_mlsadf2(x, b, m, a, aa, pd, &d[2 * (pd + 1)], ppade);

   return (x);
}

/* MGC_rnd: functions for random noise generation */
static double MGC_rnd(unsigned long *next)
{
   double r;

   *next = *next * 1103515245L + 12345;
   r = (*next / 65536L) % 32768L;

   return (r / RANDMAX);
}

/* MGC_nrandom: functions for gaussian random noise generation */
static double MGC_nrandom(MGC_Vocoder * v)
{
   if (v->sw == 0) {
      v->sw = 1;
      do {
         v->r1 = 2 * MGC_rnd(&v->next) - 1;
         v->r2 = 2 * MGC_rnd(&v->next) - 1;
         v->s = v->r1 * v->r1 + v->r2 * v->r2;
      } while (v->s > 1 || v->s == 0);
      v->s = sqrt(-2 * log(v->s) / v->s);
      return (v->r1 * v->s);
   } else {
      v->sw = 0;
      return (v->r2 * v->s);
   }
}

/* MGC_srnd: functions for gaussian random noise generation */
static unsigned long MGC_srnd(unsigned long seed)
{
   return (seed);
}

/* MGC_mceq: function for M-sequence random noise generation */
static int MGC_mseq(MGC_Vocoder * v)
{
   int x0, x28;

   v->x >>= 1;
   if (v->x & B0)
      x0 = 1;
   else
      x0 = -1;
   if (v->x & B28)
      x28 = 1;
   else
      x28 = -1;
   if (x0 + x28)
      v->x &= B31_;
   else
      v->x |= B31;

   return (x0);
}

/* MGC_mc2b: transform mel-cepstrum to MLSA digital fillter coefficients */
static void MGC_mc2b(double *mc, double *b, int m, const double a)
{
   if (mc != b) {
      if (a != 0.0) {
         b[m] = mc[m];
         for (m--; m >= 0; m--)
            b[m] = mc[m] - a * b[m + 1];
      } else
         MGC_movem(mc, b, m + 1);
   } else if (a != 0.0)
      for (m--; m >= 0; m--)
         b[m] -= a * b[m + 1];
}

/* MGC_b2bc: transform MLSA digital filter coefficients to mel-cepstrum */
static void MGC_b2mc(const double *b, double *mc, int m, const double a)
{
   double d, o;

   d = mc[m] = b[m];
   for (m--; m >= 0; m--) {
      o = b[m] + a * d;
      d = b[m];
      mc[m] = o;
   }
}

/* MGC_freqt: frequency transformation */
static void MGC_freqt(MGC_Vocoder * v, const double *c1, const int m1, double *c2, const int m2, const double a)
{
   int i, j;
   const double b = 1 - a * a;
   double *g;

   if (m2 > v->freqt_size) {
      if (v->freqt_buff != NULL)
         free(v->freqt_buff);
      v->freqt_buff = (double *) calloc(m2 + m2 + 2, sizeof(double));
      v->freqt_size = m2;
   }
   g = v->freqt_buff + v->freqt_size + 1;

   for (i = 0; i < m2 + 1; i++)
      g[i] = 0.0;

   for (i = -m1; i <= 0; i++) {
      if (0 <= m2)
         g[0] = c1[-i] + a * (v->freqt_buff[0] = g[0]);
      if (1 <= m2)
         g[1] = b * v->freqt_buff[0] + a * (v->freqt_buff[1] = g[1]);
      for (j = 2; j <= m2; j++)
         g[j] = v->freqt_buff[j - 1] + a * ((v->freqt_buff[j] = g[j]) - g[j - 1]);
   }

   MGC_movem(g, c2, m2 + 1);
}

/* MGC_c2ir: The minimum phase impulse response is evaluated from the minimum phase cepstrum */
static void MGC_c2ir(const double *c, const int nc, double *h, const int leng)
{
   int n, k, upl;
   double d;

   h[0] = exp(c[0]);
   for (n = 1; n < leng; n++) {
      d = 0;
      upl = (n >= nc) ? nc - 1 : n;
      for (k = 1; k <= upl; k++)
         d += k * c[k] * h[n - k];
      h[n] = d / n;
   }
}

/* MGC_b2en: calculate frame energy */
static double MGC_b2en(MGC_Vocoder * v, const double *b, const int m, const double a)
{
   int i;
   double en = 0.0;
   double *cep;
   double *ir;

   if (v->spectrum2en_size < m) {
      if (v->spectrum2en_buff != NULL)
         free(v->spectrum2en_buff);
      v->spectrum2en_buff = (double *) calloc((m + 1) + 2 * IRLENG, sizeof(double));
      v->spectrum2en_size = m;
   }
   cep = v->spectrum2en_buff + m + 1;
   ir = cep + IRLENG;

   MGC_b2mc(b, v->spectrum2en_buff, m, a);
   MGC_freqt(v, v->spectrum2en_buff, m, cep, IRLENG - 1, -a);
   MGC_c2ir(cep, IRLENG, ir, IRLENG);

   for (i = 0; i < IRLENG; i++)
      en += ir[i] * ir[i];

   return (en);
}

/* MGC_ignorm: inverse gain normalization */
static void MGC_ignorm(double *c1, double *c2, int m, const double g)
{
   double k;
   if (g != 0.0) {
      k = pow(c1[0], g);
      for (; m >= 1; m--)
         c2[m] = k * c1[m];
      c2[0] = (k - 1.0) / g;
   } else {
      MGC_movem(&c1[1], &c2[1], m);
      c2[0] = log(c1[0]);
   }
}

/* MGC_gnorm: gain normalization */
static void MGC_gnorm(double *c1, double *c2, int m, const double g)
{
   double k;
   if (g != 0.0) {
      k = 1.0 + g * c1[0];
      for (; m >= 1; m--)
         c2[m] = c1[m] / k;
      c2[0] = pow(k, 1.0 / g);
   } else {
      MGC_movem(&c1[1], &c2[1], m);
      c2[0] = exp(c1[0]);
   }
}

/* MGC_lsp2lpc: transform LSP to LPC */
static void MGC_lsp2lpc(MGC_Vocoder * v, double *lsp, double *a, const int m)
{
   int i, k, mh1, mh2, flag_odd;
   double xx, xf, xff;
   double *p, *q;
   double *a0, *a1, *a2, *b0, *b1, *b2;

   flag_odd = 0;
   if (m % 2 == 0)
      mh1 = mh2 = m / 2;
   else {
      mh1 = (m + 1) / 2;
      mh2 = (m - 1) / 2;
      flag_odd = 1;
   }

   if (m > v->lsp2lpc_size) {
      if (v->lsp2lpc_buff != NULL)
         free(v->lsp2lpc_buff);
      v->lsp2lpc_buff = (double *) calloc(5 * m + 6, sizeof(double));
      v->lsp2lpc_size = m;
   }
   p = v->lsp2lpc_buff + m;
   q = p + mh1;
   a0 = q + mh2;
   a1 = a0 + (mh1 + 1);
   a2 = a1 + (mh1 + 1);
   b0 = a2 + (mh1 + 1);
   b1 = b0 + (mh2 + 1);
   b2 = b1 + (mh2 + 1);

   MGC_movem(lsp, v->lsp2lpc_buff, m);

   for (i = 0; i < mh1 + 1; i++)
      a0[i] = 0.0;
   for (i = 0; i < mh1 + 1; i++)
      a1[i] = 0.0;
   for (i = 0; i < mh1 + 1; i++)
      a2[i] = 0.0;
   for (i = 0; i < mh2 + 1; i++)
      b0[i] = 0.0;
   for (i = 0; i < mh2 + 1; i++)
      b1[i] = 0.0;
   for (i = 0; i < mh2 + 1; i++)
      b2[i] = 0.0;

   /* lsp filter parameters */
   for (i = k = 0; i < mh1; i++, k += 2)
      p[i] = -2.0 * cos(v->lsp2lpc_buff[k]);
   for (i = k = 0; i < mh2; i++, k += 2)
      q[i] = -2.0 * cos(v->lsp2lpc_buff[k + 1]);

   /* impulse response of analysis filter */
   xx = 1.0;
   xf = xff = 0.0;

   for (k = 0; k <= m; k++) {
      if (flag_odd) {
         a0[0] = xx;
         b0[0] = xx - xff;
         xff = xf;
         xf = xx;
      } else {
         a0[0] = xx + xf;
         b0[0] = xx - xf;
         xf = xx;
      }

      for (i = 0; i < mh1; i++) {
         a0[i + 1] = a0[i] + p[i] * a1[i] + a2[i];
         a2[i] = a1[i];
         a1[i] = a0[i];
      }

      for (i = 0; i < mh2; i++) {
         b0[i + 1] = b0[i] + q[i] * b1[i] + b2[i];
         b2[i] = b1[i];
         b1[i] = b0[i];
      }

      if (k != 0)
         a[k - 1] = -0.5 * (a0[mh1] + b0[mh2]);
      xx = 0.0;
   }

   for (i = m - 1; i >= 0; i--)
      a[i + 1] = -a[i];
   a[0] = 1.0;
}

/* MGC_gc2gc: generalized cepstral transformation */
static void MGC_gc2gc(MGC_Vocoder * v, double *c1, const int m1, const double g1, double *c2, const int m2, const double g2)
{
   int i, min, k, mk;
   double ss1, ss2, cc;

   if (m1 > v->gc2gc_size) {
      if (v->gc2gc_buff != NULL)
         free(v->gc2gc_buff);
      v->gc2gc_buff = (double *) calloc(m1 + 1, sizeof(double));
      v->gc2gc_size = m1;
   }

   MGC_movem(c1, v->gc2gc_buff, m1 + 1);

   c2[0] = v->gc2gc_buff[0];
   for (i = 1; i <= m2; i++) {
      ss1 = ss2 = 0.0;
      min = m1 < i ? m1 : i - 1;
      for (k = 1; k <= min; k++) {
         mk = i - k;
         cc = v->gc2gc_buff[k] * c2[mk];
         ss2 += k * cc;
         ss1 += mk * cc;
      }

      if (i <= m1)
         c2[i] = v->gc2gc_buff[i] + (g2 * ss2 - g1 * ss1) / i;
      else
         c2[i] = (g2 * ss2 - g1 * ss1) / i;
   }
}

/* MGC_mgc2mgc: frequency and generalized cepstral transformation */
static void MGC_mgc2mgc(MGC_Vocoder * v, double *c1, const int m1, const double a1, const double g1, double *c2, const int m2, const double a2, const double g2)
{
   double a;

   if (a1 == a2) {
      MGC_gnorm(c1, c1, m1, g1);
      MGC_gc2gc(v, c1, m1, g1, c2, m2, g2);
      MGC_ignorm(c2, c2, m2, g2);
   } else {
      a = (a2 - a1) / (1 - a1 * a2);
      MGC_freqt(v, c1, m1, c2, m2, a);
      MGC_gnorm(c2, c2, m2, g1);
      MGC_gc2gc(v, c2, m2, g1, c2, m2, g2);
      MGC_ignorm(c2, c2, m2, g2);
   }
}

/* MGC_lsp2mgc: transform LSP to MGC */
static void MGC_lsp2mgc(MGC_Vocoder * v, double *lsp, double *mgc, const int m, const double alpha)
{
   int i;
   /* lsp2lpc */
   MGC_lsp2lpc(v, lsp + 1, mgc, m);
   if (v->use_log_gain)
      mgc[0] = exp(lsp[0]);
   else
      mgc[0] = lsp[0];

   /* mgc2mgc */
   if (NORMFLG1)
      MGC_ignorm(mgc, mgc, m, v->gamma);
   else if (MULGFLG1)
      mgc[0] = (1.0 - mgc[0]) * v->stage;
   if (MULGFLG1)
      for (i = m; i >= 1; i--)
         mgc[i] *= -v->stage;
   MGC_mgc2mgc(v, mgc, m, alpha, v->gamma, mgc, m, alpha, v->gamma);
   if (NORMFLG2)
      MGC_gnorm(mgc, mgc, m, v->gamma);
   else if (MULGFLG2)
      mgc[0] = mgc[0] * v->gamma + 1.0;
   if (MULGFLG2)
      for (i = m; i >= 1; i--)
         mgc[i] *= v->gamma;
}

/* MGC_mglsadff: sub functions for MGLSA filter */
static double MGC_mglsadff(double x, const double *b, const int m, const double a, double *d)
{
   int i;

   double y;
   y = d[0] * b[1];
   for (i = 1; i < m; i++) {
      d[i] += a * (d[i + 1] - d[i - 1]);
      y += d[i] * b[i + 1];
   }
   x -= y;

   for (i = m; i > 0; i--)
      d[i] = d[i - 1];
   d[0] = a * d[0] + (1 - a * a) * x;
   return x;
}

/* MGC_mglsadf: sub functions for MGLSA filter */
static double MGC_mglsadf(double x, const double *b, const int m, const double a, const int n, double *d)
{
   int i;

   for (i = 0; i < n; i++)
      x = MGC_mglsadff(x, b, m, a, &d[i * (m + 1)]);

   return x;
}

/* MGC_white_noise: return white noise */
static double MGC_white_noise(MGC_Vocoder * v)
{
   if (v->gauss)
      return (double) MGC_nrandom(v);
   else
      return MGC_mseq(v);
}

/* MGC_Vocoder_initialize_excitation: initialize excitation */
static void MGC_Vocoder_initialize_excitation(MGC_Vocoder * v)
{
   v->p1 = 0.0;
   v->pc = 0.0;
   v->p = 0.0;
   v->inc = 0.0;
}

/* MGC_Vocoder_start_excitation: start excitation of each frame */
static void MGC_Vocoder_start_excitation(MGC_Vocoder * v, const double pitch)
{
   if (v->p1 != 0.0 && pitch != 0.0)
      v->inc = (pitch - v->p1) * v->iprd / v->fprd;
   else {
      v->inc = 0.0;
      v->p1 = 0.0;
      v->pc = pitch;
   }
   v->p = pitch;
}

/* MGC_Vocoder_get_excitation: get excitation of each sample */
static double MGC_Vocoder_get_excitation(MGC_Vocoder * v, const int fprd_index, const int iprd_index)
{
   double x;

    if (v->p1 == 0.0)
        x = MGC_white_noise(v);
    else {
        if ((v->pc += 1.0) >= v->p1) {
        x = sqrt(v->p1);
        v->pc -= v->p1;
        } else
        x = 0.0;
    }
    if (iprd_index <= 1)
        v->p1 += v->inc;

   return x;
}

/* MGC_Vocoder_end_excitation: end excitation of each frame */
static void MGC_Vocoder_end_excitation(MGC_Vocoder * v)
{
   v->p1 = v->p;
}

/* MGC_Vocoder_initialize: initialize vocoder */
void MGC_Vocoder_initialize(MGC_Vocoder * v, const int m, const int stage, Boolean use_log_gain, const int rate, const int fperiod)
{
   /* set parameter */
   v->stage = stage;
   if (stage != 0)
      v->gamma = -1.0 / v->stage;
   else
      v->gamma = 0.0;
   v->use_log_gain = use_log_gain;
   v->fprd = fperiod;
   v->iprd = IPERIOD;
   v->seed = SEED;
   v->next = SEED;
   v->gauss = GAUSS;
   v->rate = rate;
   v->p1 = -1.0;
   v->sw = 0;
   v->x = 0x55555555;
   /* init buffer */
   v->freqt_buff = NULL;
   v->freqt_size = 0;
   v->gc2gc_buff = NULL;
   v->gc2gc_size = 0;
   v->lsp2lpc_buff = NULL;
   v->lsp2lpc_size = 0;
   v->postfilter_buff = NULL;
   v->postfilter_size = 0;
   v->spectrum2en_buff = NULL;
   v->spectrum2en_size = 0;
   if (v->stage == 0) {         /* for MCP */
      v->c = (double *) calloc(m * (3 + PADEORDER) + 5 * PADEORDER + 6, sizeof(double));
      v->cc = v->c + m + 1;
      v->cinc = v->cc + m + 1;
      v->d1 = v->cinc + m + 1;
   } else {                     /* for LSP */
      v->c = (double *) calloc((m + 1) * (v->stage + 3), sizeof(double));
      v->cc = v->c + m + 1;
      v->cinc = v->cc + m + 1;
      v->d1 = v->cinc + m + 1;
   }
}

/* MGC_Vocoder_synthesize: pulse/noise excitation and MLSA/MGLSA filster based waveform synthesis */
void MGC_Vocoder_synthesize(MGC_Vocoder * v, const int m, double lf0, double *spectrum, double alpha, double beta, short *rawdata)
{
   double x;
   int i, j;
   short xs;
   int rawidx = 0;
   double p;

   /* lf0 -> pitch */
   if (lf0 == LZERO)
      p = 0.0;
   else
      p = v->rate / exp(lf0);

   /* first time */
   if (v->p1 < 0.0) {
      if (v->gauss & (v->seed != 1))
         v->next = MGC_srnd((unsigned) v->seed);
      MGC_Vocoder_initialize_excitation(v);
      if (v->stage == 0) {      /* for MCP */
         MGC_mc2b(spectrum, v->c, m, alpha);
      } else {                  /* for LSP */
         if (v->use_log_gain)
            v->c[0] = LZERO;
         else
            v->c[0] = ZERO;
         for (i = 1; i <= m; i++)
            v->c[i] = i * PI / (m + 1);
         MGC_lsp2mgc(v, v->c, v->c, m, alpha);
         MGC_mc2b(v->c, v->c, m, alpha);
         MGC_gnorm(v->c, v->c, m, v->gamma);
         for (i = 1; i <= m; i++)
            v->c[i] *= v->gamma;
      }
   }

   MGC_Vocoder_start_excitation(v, p);
   if (v->stage == 0) {         /* for MCP */
      MGC_Vocoder_postfilter_mcp(v, spectrum, m, alpha, beta);
      MGC_mc2b(spectrum, v->cc, m, alpha);
      for (i = 0; i <= m; i++)
         v->cinc[i] = (v->cc[i] - v->c[i]) * v->iprd / v->fprd;
   } else {                     /* for LSP */
      MGC_lsp2mgc(v, spectrum, v->cc, m, alpha);
      MGC_mc2b(v->cc, v->cc, m, alpha);
      MGC_gnorm(v->cc, v->cc, m, v->gamma);
      for (i = 1; i <= m; i++)
         v->cc[i] *= v->gamma;
      for (i = 0; i <= m; i++)
         v->cinc[i] = (v->cc[i] - v->c[i]) * v->iprd / v->fprd;
   }

   for (j = 0, i = (v->iprd + 1) / 2; j < v->fprd; j++) {
      x = MGC_Vocoder_get_excitation(v, j, i);
      if (v->stage == 0) {      /* for MCP */
         if (x != 0.0)
            x *= exp(v->c[0]);
         x = MGC_mlsadf(x, v->c, m, alpha, PADEORDER, v->d1);
      } else {                  /* for LSP */
         if (!NGAIN)
            x *= v->c[0];
         x = MGC_mglsadf(x, v->c, m, alpha, v->stage, v->d1);
      }

      /* output */
      if (x > 32767.0)
         xs = 32767;
      else if (x < -32768.0)
         xs = -32768;
      else
         xs = (short) x;
      if (rawdata)
         rawdata[rawidx++] = xs;

      if (!--i) {
         for (i = 0; i <= m; i++)
            v->c[i] += v->cinc[i];
         i = v->iprd;
      }
   }

   MGC_Vocoder_end_excitation(v);
   MGC_movem(v->cc, v->c, m + 1);
}

/* MGC_Vocoder_postfilter_mcp: postfilter for MCP */
void MGC_Vocoder_postfilter_mcp(MGC_Vocoder * v, double *mcp, const int m, double alpha, double beta)
{
   double e1, e2;
   int k;

   if (beta > 0.0 && m > 1) {
      if (v->postfilter_size < m) {
         if (v->postfilter_buff != NULL)
            free(v->postfilter_buff);
         v->postfilter_buff = (double *) calloc(m + 1, sizeof(double));
         v->postfilter_size = m;
      }
      MGC_mc2b(mcp, v->postfilter_buff, m, alpha);
      e1 = MGC_b2en(v, v->postfilter_buff, m, alpha);

      v->postfilter_buff[1] -= beta * alpha * mcp[2];
      for (k = 2; k <= m; k++)
         v->postfilter_buff[k] *= (1.0 + beta);

      e2 = MGC_b2en(v, v->postfilter_buff, m, alpha);
      v->postfilter_buff[0] += log(e1 / e2) / 2;
      MGC_b2mc(v->postfilter_buff, mcp, m, alpha);
   }
}

/* MGC_Vocoder_clear: clear vocoder */
void MGC_Vocoder_clear(MGC_Vocoder * v)
{
   if (v != NULL) {
      /* free buffer */
      if (v->freqt_buff != NULL) {
         free(v->freqt_buff);
         v->freqt_buff = NULL;
      }
      v->freqt_size = 0;
      if (v->gc2gc_buff != NULL) {
         free(v->gc2gc_buff);
         v->gc2gc_buff = NULL;
      }
      v->gc2gc_size = 0;
      if (v->lsp2lpc_buff != NULL) {
         free(v->lsp2lpc_buff);
         v->lsp2lpc_buff = NULL;
      }
      v->lsp2lpc_size = 0;
      if (v->postfilter_buff != NULL) {
         free(v->postfilter_buff);
         v->postfilter_buff = NULL;
      }
      v->postfilter_size = 0;
      if (v->spectrum2en_buff != NULL) {
         free(v->spectrum2en_buff);
         v->spectrum2en_buff = NULL;
      }
      v->spectrum2en_size = 0;
      if (v->c != NULL) {
         free(v->c);
         v->c = NULL;
      }
   }
}