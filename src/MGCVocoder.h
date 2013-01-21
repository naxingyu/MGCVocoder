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

#include <stdio.h>

/*  -------------------------- common -----------------------------  */

typedef int Boolean;
#ifndef TRUE
#define TRUE  1
#endif				/* !TRUE */
#ifndef FALSE
#define FALSE 0
#endif				/* !FALSE */

#define ZERO  1.0e-10		/* ~(0) */
#define LZERO (-1.0e+10)	/* ~log(0) */
#define LTPI  1.83787706640935	/* log(2*PI) */

/*  -------------------------- vocoder ----------------------------  */

/* MGC_Vocoder: structure for setting of vocoder */
typedef struct _MGC_Vocoder {
   int stage;			/* Gamma=-1/stage: if stage=0 then Gamma=0 */
   double gamma;		/* Gamma */
   Boolean use_log_gain;	/* log gain flag (for LSP) */
   int fprd;			/* frame shift */
   int iprd;			/* interpolation period */
   int seed;			/* seed of random generator */
   unsigned long next;		/* temporary variable for random generator */
   Boolean gauss;		/* flag to use Gaussian noise */
   double rate;			/* sampling rate */
   double p1;			/* used in excitation generation */
   double pc;			/* used in excitation generation */
   double p;			/* used in excitation generation */
   double inc;			/* used in excitation generation */
   int sw;			/* switch used in random generator */
   int x;			/* excitation signal */
   double *freqt_buff;		/* used in freqt */
   int freqt_size;		/* buffer size for freqt */
   double *spectrum2en_buff;	/* used in spectrum2en */
   int spectrum2en_size;	/* buffer size for spectrum2en */
   double r1, r2, s;		/* used in random generator */
   double *postfilter_buff;	/* used in postfiltering */
   int postfilter_size;		/* buffer size for postfiltering */
   double *c, *cc, *cinc, *d1;	/* used in the MLSA/MGLSA filter */
   double *lsp2lpc_buff;	/* used in lsp2lpc */
   int lsp2lpc_size;		/* buffer size of lsp2lpc */
   double *gc2gc_buff;		/* used in gc2gc */
   int gc2gc_size;		/* buffer size for gc2gc */
} MGC_Vocoder;

/*  ----------------------- vocoder method ------------------------  */

/* MGC_Vocoder_initialize: initialize vocoder */
void MGC_Vocoder_initialize(MGC_Vocoder * v, const int m, const int stage,
			    Boolean use_log_gain, const int rate,
			    const int fperiod);

/* MGC_Vocoder_synthesize: pulse/noise excitation and MLSA/MGLSA filster based waveform synthesis */
void MGC_Vocoder_synthesize(MGC_Vocoder * v, const int m, double lf0,
			    double *spectrum, double alpha, double beta,
			    short *rawdata);

/* MGC_Vocoder_postfilter_mcp: postfilter for MCP */
void MGC_Vocoder_postfilter_mcp(MGC_Vocoder * v, double *mcp, const int m,
				double alpha, double beta);

/* MGC_Vocoder_clear: clear vocoder */
void MGC_Vocoder_clear(MGC_Vocoder * v);
