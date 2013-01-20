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

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "MGCVocoder.h"

/* Usage: output usage */
void Usage(void)
{
   fprintf(stderr, "\n");
   fprintf(stderr, "MGCVocoder - The MGC Vocoder (part of hts_engine)\n");
   fprintf(stderr, "\n");
   fprintf(stderr, "  usage:\n");
   fprintf(stderr, "       MGCVocoder [ options ] \n");
   fprintf(stderr, "  options:                                                                   [  def][ min--max]\n");
   fprintf(stderr, "    -im mgc        : input file for spectrum                                 [  N/A]\n");
   fprintf(stderr, "    -if lf0        : input file for Log F0                                   [  N/A]\n");
   fprintf(stderr, "    -or raw        : filename of output raw audio (generated speech)         [  N/A]\n");
   fprintf(stderr, "    -ow wav        : filename of output wav audio (generated speech)         [  N/A]\n");
   fprintf(stderr, "    -m  i          : order of MGC/LSP                                        [   24][   1--]\n");
   fprintf(stderr, "    -s  i          : sampling frequency                                      [16000][   1--48000]\n");
   fprintf(stderr, "    -p  i          : frame period (point)                                    [   80][   1--]\n");
   fprintf(stderr, "    -a  f          : all-pass constant                                       [ 0.42][ 0.0--1.0]\n");
   fprintf(stderr, "    -g  i          : gamma = -1 / i (if i=0 then gamma=0)                    [    0][   0-- ]\n");
   fprintf(stderr, "    -b  f          : postfiltering coefficient                               [  0.0][-0.8--0.8]\n");
   fprintf(stderr, "    -l             : regard input as log gain (LSP)                          [  N/A]\n");
   fprintf(stderr, "    -t             : regard input as STRAIGHT spectrum                       [  N/A]\n");
   fprintf(stderr, "    -h             : print this help                                         [  N/A]\n");
   fprintf(stderr, "  note:\n");
   fprintf(stderr, "    Both MGC and LSP spectrum are supported.\n");
   fprintf(stderr, "    Input spectrum and log F0 sequences are\n");
   fprintf(stderr, "    saved in natural endian, binary (float) format.\n");
   fprintf(stderr, "\n");

   exit(0);
}

/* Error: output error message */
void Error(const int error, char *message, ...)
{
   va_list arg;

   fflush(stdout);
   fflush(stderr);

   if (error > 0)
      fprintf(stderr, "\nError: ");
   else
      fprintf(stderr, "\nWarning: ");

   va_start(arg, message);
   vfprintf(stderr, message, arg);
   va_end(arg);

   fflush(stderr);

   if (error > 0)
      exit(error);
}

/* MGC_fopen: wrapper for fopen */
FILE *MGC_fopen(const char *name, const char *opt)
{
   FILE *fp = fopen(name, opt);

   if (fp == NULL) {
      Error(1, "MGC_fopen: Cannot open %s.\n", name);
      return NULL;
   }

   return fp;
}

void WAVE_save_header(int total_nsample, int sampling_rate, FILE *wavfp)
{
	int i;
	short temp;

	char data_01_04[] = { 'R', 'I', 'F', 'F' };
	int data_05_08 = total_nsample * sizeof(short) + 36;
	char data_09_12[] = { 'W', 'A', 'V', 'E' };
	char data_13_16[] = { 'f', 'm', 't', ' ' };
	int data_17_20 = 16;
	short data_21_22 = 1;        /* PCM */
	short data_23_24 = 1;        /* monoral */
	int data_25_28 = sampling_rate;
	int data_29_32 = sampling_rate * sizeof(short);
	short data_33_34 = sizeof(short);
	short data_35_36 = (short) (sizeof(short) * 8);
	char data_37_40[] = { 'd', 'a', 't', 'a' };
	int data_41_44 = total_nsample * sizeof(short);

	/* write header */
	fwrite(data_01_04, sizeof(char), 4, wavfp);
	fwrite(&data_05_08, sizeof(int), 1, wavfp);
	fwrite(data_09_12, sizeof(char), 4, wavfp);
	fwrite(data_13_16, sizeof(char), 4, wavfp);
	fwrite(&data_17_20, sizeof(int), 1, wavfp);
	fwrite(&data_21_22, sizeof(short), 1, wavfp);
	fwrite(&data_23_24, sizeof(short), 1, wavfp);
	fwrite(&data_25_28, sizeof(int), 1, wavfp);
	fwrite(&data_29_32, sizeof(int), 1, wavfp);
	fwrite(&data_33_34, sizeof(short), 1, wavfp);
	fwrite(&data_35_36, sizeof(short), 1, wavfp);
	fwrite(data_37_40, sizeof(char), 4, wavfp);
	fwrite(&data_41_44, sizeof(int), 1, wavfp);
}

int main(int argc, char **argv)
{
	FILE *mgcfp = NULL, *lf0fp = NULL, *wavfp = NULL, *rawfp = NULL;
	int total_frame;
	double lf0, *spectrum;
	MGC_Vocoder v;
	short *rawdata;              /* generated speech */

	/* global parameter */
	int order = 24;
	int sampling_rate = 16000;
	int fperiod = 80;
	double alpha = 0.42;
	int stage = 0;               /* Gamma=-1/stage: if stage=0 then Gamma=0 */
	double beta = 0.0;
	Boolean use_log_gain = FALSE;
	Boolean use_straight = FALSE;

	/* temp variables */
	int i, j, itemp;
	float ftemp;

	/* parse command line */
   if (argc == 1)
      Usage();
	/* read command */
	while (--argc) {
		if (**++argv == '-') {
			switch (*(*argv + 1)) {
			case 'i':
				switch (*(*argv + 2)) {
				case 'm':
					mgcfp = MGC_fopen(*++argv, "rb");
					break;
				case 'f':
					lf0fp = MGC_fopen(*++argv, "rb");
					break;
				default:
					Error(1, "MGCVocoder: Invalid option '-i%c'.\n", *(*argv + 2));
				}
				--argc;
				break;
			case 'o':
				switch (*(*argv + 2)) {
				case 'r':
					rawfp = MGC_fopen(*++argv, "wb");
					break;
				case 'w':
					wavfp = MGC_fopen(*++argv, "wb");
					break;
				default:
					Error(1, "MGCVocoder: Invalid option '-o%c'.\n", *(*argv + 2));
				}
				--argc;
				break;
			case 'm':
				order = atoi(*++argv);
				--argc;
				break;
			case 's':
				sampling_rate = atoi(*++argv);
				--argc;
				break;
			case 'p':
				fperiod = atoi(*++argv);
				--argc;
				break;
			case 'a':
				alpha = atof(*++argv);
				--argc;
				break;
			case 'g':
				stage = atoi(*++argv);
				--argc;
				break;
			case 'b':
				beta = atof(*++argv);
				--argc;
				break;
			case 'l':
				use_log_gain = TRUE;
				break;
			case 't':
				use_straight = TRUE;
				break;
			case 'h':
				Usage();
				break;
			default:
				Error(1, "MGCVocoder: Invalid option '-%c'.\n", *(*argv + 1));
			}
		}
	}

	/* initialize */
	fseek(mgcfp, 0, SEEK_END);
	total_frame = ftell(mgcfp)/(order*sizeof(float));
	fseek(mgcfp, 0, SEEK_SET);

	fseek(lf0fp, 0, SEEK_END);
	itemp = ftell(lf0fp)/sizeof(float);
	fseek(lf0fp, 0, SEEK_SET);

	if(itemp < total_frame)
		total_frame = itemp;

	fprintf(stderr, "\nTotal length: %dframes, %d ms.\n", total_frame, total_frame*5);
	rawdata = (short *)calloc(total_frame * fperiod, sizeof(short));
	spectrum = (double *)calloc(order + 1, sizeof(double));

	/* synthesize speech waveform */
	MGC_Vocoder_initialize(&v, order, stage, use_log_gain, sampling_rate, fperiod);

	for (i = 0; i < total_frame; i++) {
		/* copy parameter */
		fread(&ftemp, sizeof(float), 1, lf0fp);
		lf0 = (double)ftemp;
		for (j = 0; j < order + 1; j++) {
			fread(&ftemp, sizeof(float), 1, mgcfp);
			spectrum[j] = (double)ftemp;
		}
		if (use_straight)
			spectrum[0] += 10.0;
		MGC_Vocoder_synthesize(&v, order, lf0, spectrum, alpha, beta, &rawdata[i * fperiod]);
	}
	MGC_Vocoder_clear(&v);
	if (rawfp)
		for (i = 0; i < total_frame * fperiod; i++)
			fwrite(&rawdata[i], sizeof(short), 1, rawfp);
	if (wavfp) {
		WAVE_save_header(total_frame * fperiod, sampling_rate, wavfp);
		for (i = 0; i < total_frame * fperiod; i++)
			fwrite(&rawdata[i], sizeof(short), 1, wavfp);
	}

	/* free memory */
	free(spectrum);
	free(rawdata);

	/* close files */
	if (mgcfp != NULL)
		fclose(mgcfp);
	if (lf0fp != NULL)
		fclose(lf0fp);
	if (wavfp != NULL)
		fclose(wavfp);
	if (rawfp != NULL)
		fclose(rawfp);
}