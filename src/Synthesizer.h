/* ----------------------------------------------- */
/*        The MGC Vocoder (part of hts_engine)     */
/*        developed by HTS Working Group           */
/*        maintained by Xingyu Na                  */
/* ----------------------------------------------- */
/* Copyright (c)                                   */
/*   2011-2012  Beijing Institute of Technology    */
/*              RCDCT, @RICT                       */
/*   2001-2011  Nagoya Institute of Technology     */
/*              Department of Computer Science     */
/* ----------------------------------------------- */

#include "MGCVocoder.h"

/*  -------------------------- vocoder ----------------------------  */

#ifndef PI
#define PI  3.14159265358979323846
#endif                          /* !PI */
#ifndef PI2
#define PI2 6.28318530717958647692
#endif                          /* !PI2 */

#define RANDMAX 32767

#define IPERIOD 1
#define SEED    1
#define B0      0x00000001
#define B28     0x10000000
#define B31     0x80000000
#define B31_    0x7fffffff
#define Z       0x00000000
#define GAUSS     TRUE
#define PADEORDER 5
#define IRLENG    96

/* for MGLSA filter */
#define NORMFLG1 TRUE
#define NORMFLG2 FALSE
#define MULGFLG1 TRUE
#define MULGFLG2 FALSE
#define NGAIN    FALSE
