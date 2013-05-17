/*
 * This material contains unpublished, proprietary software of 
 * Entropic Research Laboratory, Inc. Any reproduction, distribution, 
 * or publication of this work must be authorized in writing by Entropic 
 * Research Laboratory, Inc., and must bear the notice: 
 *
 *    "Copyright (c) 1990-1996 Entropic Research Laboratory, Inc. 
 *                   All rights reserved"
 *
 * The copyright notice above does not evidence any actual or intended 
 * publication of this source code.     
 *
 * Written by:  Derek Lin
 * Checked by:
 * Revised by:  David Talkin
 *
 * Brief description:  Estimates F0 using normalized cross correlation and
 *   dynamic programming.
 *
 */

static char *sccs_id = "@(#)get_f0.c	1.14	10/21/96	ERL";

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
/* #include <malloc.h> */

#include "dpwe_add.h"

#include <snd.h>

#include "f0.h"

#define USAGE(s)  fprintf(stderr, "%s\n", (s)); exit(1)

#define SYNTAX USAGE("get_f0 [-P param_file][-{pr} range][-s range][-S frame_step_samples]\n     [-i frame_step_seconds][-x debug_level] infile outfile")

char	    *ProgName = "get_f0";
static char *Version = "1.14";
static char *Date = "10/21/96";

int	    debug_level = 0;

static int check_f0_params();
int init_dp_f0();
int dp_f0();
static void get_range();

main(ac, av)
    int     ac;
    char    **av;
{
  extern char *optarg;
  extern int optind, getopt();
  char *get_cmd_line();
  void get_range();
  float *fdata;
  char c, *ifname, *ofname, *range = NULL;
  FILE *ofile;
  short n_cands;
  int done;
  long buff_size, actsize, s_rec, e_rec;
  double sf, start_time, output_starts, frame_rate;
  F0_params *par, *read_f0_params();
  char *param_file = NULL;
  float *f0p, *vuvp, *rms_speech, *acpkp;
  float tt = 0;
#ifdef ESPS_INPUT
  struct header *ihd;
  FILE *ifile;
#else
  SOUND *ifile;
#endif
#ifdef ESPS_OUTPUT
  struct header *ohd;
  struct feasd *sd_rec;
  struct fea_data *fea_rec;
  double *rec_F0, *rec_pv, *rec_rms, *rec_acp;
#endif
  int i, ii, vecsize;
  int framestep = -1, rflag = 0,
      sflag = 0, iflag = 0;
  long sdstep = 0, total_samps;

  par = (F0_params *) malloc(sizeof(F0_params));
  par->cand_thresh = 0.3;
  par->lag_weight = 0.3;
  par->freq_weight = 0.02;
  par->trans_cost = 0.005;
  par->trans_amp = 0.5;
  par->trans_spec = 0.5;
  par->voice_bias = 0.0;
  par->double_cost = 0.35;
  par->min_f0 = 50;
  par->max_f0 = 550;
  par->frame_step = 0.01;
  par->wind_dur = 0.0075;
  par->n_cands = 20;
  par->mean_f0 = 200;     /* unused */
  par->mean_f0_weight = 0.0;  /* unused */
  par->conditioning = 0;    /*unused */

  while((c = getopt(ac,av,"x:P:p:r:s:S:i:")) != EOF){
    switch(c){
    case 'P':
      param_file = optarg;
      break;
    case 'p':
    case 'r':
      if( range ){
	Fprintf(stderr, "%s: error: -r should not be used with -s.\n",
		ProgName);
	exit(1);
      }
      range = optarg;
      rflag++;
      break;
    case 's':
      if( range ){
	Fprintf(stderr, "%s: error: -s should not be used with -r.\n",
		ProgName);
	exit(1);
      }
      range = optarg;
      sflag++;
      break;
    case 'S':
      if( iflag ){
	Fprintf(stderr, "%s: error: -S should not be used with -i.\n",
		ProgName);
	exit(1);
      }
      framestep = atoi(optarg);
      break;
    case 'i':
      if(framestep > 0){
	Fprintf(stderr, "%s: error: -i should not be used with -S.\n",
		ProgName);
	exit(1);
      }
      par->frame_step = atof(optarg);
      iflag++;
      break;
    case 'x':
      debug_level = atoi(optarg);
      break;
    default:
      SYNTAX;
      exit(1);
    }
  }
  
  if((ac - optind) != 2){
    SYNTAX;
    exit(1);
  }
  
  (void) read_params(param_file, SC_NOCOMMON, (char *)NULL);
  
  if( framestep < 0 && !iflag && symtype("frame_step") != ST_UNDEF)
    par->frame_step = getsym_d("frame_step");

  if( symtype("cand_thresh") != ST_UNDEF)
    par->cand_thresh = getsym_d("cand_thresh");
  if( symtype("lag_weight") != ST_UNDEF)
    par->lag_weight = getsym_d("lag_weight");
  if( symtype("freq_weight") != ST_UNDEF)
    par->freq_weight = getsym_d("freq_weight");
  if( symtype("trans_cost") != ST_UNDEF)
    par->trans_cost = getsym_d("trans_cost");
  if( symtype("trans_amp") != ST_UNDEF)
    par->trans_amp = getsym_d("trans_amp");
  if( symtype("trans_spec") != ST_UNDEF)
    par->trans_spec = getsym_d("trans_spec");
  if( symtype("voice_bias") != ST_UNDEF)
    par->voice_bias = getsym_d("voice_bias");
  if( symtype("double_cost") != ST_UNDEF)
    par->double_cost = getsym_d("double_cost");
  if( symtype("min_f0") != ST_UNDEF)
    par->min_f0 = getsym_d("min_f0");
  if( symtype("max_f0") != ST_UNDEF)
    par->max_f0 = getsym_d("max_f0");
  if( symtype("wind_dur") != ST_UNDEF)
    par->wind_dur = getsym_d("wind_dur");
  if( symtype("n_cands") != ST_UNDEF)
    par->n_cands = getsym_i("n_cands");

#ifdef ESPS_INPUT
  ifname = eopen(ProgName, av[optind], "r", FT_FEA, FEA_SD, &ihd, &ifile);
#else
  ifile = sndOpen(NULL, av[optind], "r");
#endif

#ifdef ESPS_OUTPUT
  ofname = eopen(ProgName, av[optind+1], "w", NONE, NONE, &ohd, &ofile);
#else
  ofile = fopen(av[optind+1], "w");
  assert(ofile != NULL);
#endif

#ifdef ESPS_INPUT
  sf = get_genhd_val("record_freq", ihd, 0.0);
#else
  sf = sndGetSrate(ifile);
#endif
  if (sf == 0.0) {
    Fprintf(stderr, "%s: no sampling frequency---exiting.\n", ProgName);
    exit(1);
  }
  if (framestep > 0)  /* If a value was specified with -S, use it. */
    par->frame_step = framestep / sf;
#ifdef ESPS_INPUT
  start_time = get_genhd_val("start_time", ihd, 0.0);
#else
  start_time = 0;
#endif
  if(check_f0_params(par, sf)){
    Fprintf(stderr, "%s: invalid/inconsistent parameters -- exiting.\n",
	    ProgName);
    exit(1);
  }

#ifdef ESPS_INPUT
  get_range( &s_rec, &e_rec, range, rflag, sflag, ihd );
#else
  s_rec = 0;
  e_rec = sndGetFrames(ifile);
#endif
  if(debug_level) 
    Fprintf(stderr, "%s: input sample range: %ld - %ld.\n",
	    ProgName, s_rec, e_rec);
  total_samps = e_rec - s_rec + 1;
  if(total_samps < ((par->frame_step * 2.0) + par->wind_dur) * sf) {
    Fprintf(stderr, "%s: input range too small for analysis by get_f0.\n",
	    ProgName);
    SYNTAX;
    exit(1);
  }
#ifdef ESPS_OUTPUT
  ohd = new_header(FT_FEA);
  if (ohd == NULL) {
    Fprintf(stderr, "%s: failed to create output header---exiting.\n",
	    ProgName);
    exit(1);
  }
  (void) strcpy (ohd->common.prog, ProgName);
  (void) strcpy (ohd->common.vers, Version);
  (void) strcpy (ohd->common.progdate, Date);
  ohd->common.tag = NO;
  add_source_file(ohd,ifname,ihd);
  add_comment(ohd,get_cmd_line(ac,av));

  add_fea_fld("F0", 1L, 0, (long *) NULL, DOUBLE, (char **) NULL, ohd);
  add_fea_fld("prob_voice", 1L, 0, (long *) NULL, DOUBLE, (char **) NULL, ohd);
  add_fea_fld("rms", 1L, 0, (long *) NULL, DOUBLE, (char **) NULL, ohd);
  add_fea_fld("ac_peak", 1L, 0, (long *) NULL, DOUBLE, (char **) NULL, ohd);
  fea_rec = allo_fea_rec(ohd);
  rec_F0 = (double *) get_fea_ptr(fea_rec,"F0", ohd);
  rec_pv = (double *) get_fea_ptr(fea_rec,"prob_voice", ohd);
  rec_rms = (double *) get_fea_ptr(fea_rec,"rms", ohd);
  rec_acp = (double *) get_fea_ptr(fea_rec,"ac_peak", ohd);
#endif

  output_starts = ((s_rec-1) / sf) + start_time + par->wind_dur/2.0; 
  /* Average delay due to loc. of ref. window center. */
  frame_rate = 1.0 / par->frame_step;

#ifdef ESPS_OUTPUT
  (void)add_genhd_d("record_freq", &frame_rate, 1, ohd);
  (void)add_genhd_d("start_time", &output_starts, 1, ohd);
  (void)add_genhd_d("src_sf", &sf, 1, ohd);
  (void)add_genhd_f("cand_thresh", &par->cand_thresh, 1, ohd);
  (void)add_genhd_f("lag_weight", &par->lag_weight, 1, ohd);
  (void)add_genhd_f("freq_weight", &par->freq_weight, 1, ohd);
  (void)add_genhd_f("trans_cost", &par->trans_cost, 1, ohd);
  (void)add_genhd_f("trans_amp", &par->trans_amp, 1, ohd);
  (void)add_genhd_f("trans_spec", &par->trans_spec, 1, ohd);
  (void)add_genhd_f("voice_bias", &par->voice_bias, 1, ohd);
  (void)add_genhd_f("double_cost", &par->double_cost, 1, ohd);
  (void)add_genhd_f("min_f0", &par->min_f0, 1, ohd);
  (void)add_genhd_f("max_f0", &par->max_f0, 1, ohd);
  (void)add_genhd_f("frame_step", &par->frame_step, 1, ohd);
  (void)add_genhd_f("wind_dur", &par->wind_dur, 1, ohd);
  n_cands = (short) par->n_cands;
  (void)add_genhd_s("n_cands", &n_cands, 1, ohd);

  write_header(ohd, ofile);
#endif

  /* Initialize variables in get_f0.c; allocate data structures;
   * determine length and overlap of input frames to read.
   */
  if (init_dp_f0(sf, par, &buff_size, &sdstep)
      || buff_size > INT_MAX || sdstep > INT_MAX)
  {
    Fprintf(stderr, "%s: problem in init_dp_f0().\n", ProgName);
    exit(1);
  }

  if (debug_level)
    Fprintf(stderr, "%s: init_dp_f0 returned buff_size %ld, sdstep %ld.\n",
	    ProgName, buff_size, sdstep);

#ifdef ESPS_INPUT
  sd_rec = allo_feasd_recs(ihd, FLOAT, buff_size, NULL, NO);
  fdata = (float *) sd_rec->data;
  
  fea_skiprec(ifile, s_rec - 1, ihd);
#else
  fdata = (float *)malloc(buff_size*sizeof(float));
#endif

  if (buff_size > total_samps)
    buff_size = total_samps;

#if ESPS_INPUT
  actsize = get_feasd_recs(sd_rec, 0L, buff_size, ihd, ifile);
#else
  actsize = sndReadFrames(ifile, fdata, buff_size);
  for (ii = 0; ii < actsize; ++ii) fdata[ii] *= 32768.0;
#endif

  while (TRUE) {

    done = (actsize < buff_size) || (total_samps == buff_size);

    if (dp_f0(fdata, (int) actsize, (int) sdstep, sf, par,
	      &f0p, &vuvp, &rms_speech, &acpkp, &vecsize, done)) {
      Fprintf(stderr, "%s: problem in dp_f0().\n", ProgName);
      exit(1);
    }

    for (i = vecsize - 1; i >= 0; i--) {
#ifdef ESPS_OUTPUT
      *rec_F0 = f0p[i];
      *rec_pv = vuvp[i];
      *rec_rms = rms_speech[i];
      *rec_acp = acpkp[i];
      put_fea_rec(fea_rec, ohd, ofile);
#else /* ascii output */
      fprintf(ofile, "%.3f\t%.1f\t%.3f\t%.1f\t%.3f\n", tt, f0p[i], vuvp[i], rms_speech[i], acpkp[i]);
      tt += par->frame_step;
#endif
    }

    if (done)
      break;
    
#if ESPS_INPUT
    actsize = get_feasd_orecs( sd_rec, 0L, buff_size, sdstep, ihd, ifile);
#else
    /* fill buff_size frames, reading in sdstep new ones */
    int keeps = buff_size - sdstep;
    for (ii = 0; ii < keeps; ++ii)
	fdata[ii] = fdata[sdstep + ii];
    int newread = sndReadFrames(ifile, fdata + keeps, sdstep);
    for (ii = 0; ii < newread; ++ii) fdata[keeps+ii] *= 32768.0;
    actsize = keeps + newread;
#endif
    total_samps -= sdstep;

    if (actsize > total_samps)
      actsize = total_samps;
  }

#ifndef ESPS_OUTPUT
  fclose(ofile);
#endif 

  exit(0);
}


/*
 * Some consistency checks on parameter values.
 * Return a positive integer if any errors detected, 0 if none.
 */

static int
check_f0_params(par, sample_freq)
    F0_params   *par;
    double      sample_freq;
{
  int	  error = 0;
  double  dstep;

  if((par->cand_thresh < 0.01) || (par->cand_thresh > 0.99)) {
    Fprintf(stderr,
	    "%s: ERROR: cand_thresh parameter must be between [0.01, 0.99].\n",
	    ProgName);
    error++;
  }
  if((par->wind_dur > .1) || (par->wind_dur < .0001)) {
    Fprintf(stderr,
	    "ERROR: wind_dur parameter must be between [0.0001, 0.1].\n",
	    ProgName);
    error++;
  }
  if((par->n_cands > 100) || (par->n_cands < 3)){
    Fprintf(stderr,
	    "%s: ERROR: n_cands parameter must be between [3,100].\n",
	    ProgName); 
    error++;
  }
  if((par->max_f0 <= par->min_f0) || (par->max_f0 >= (sample_freq/2.0)) ||
     (par->min_f0 < (sample_freq/10000.0))){
    Fprintf(stderr,
	    "%s: ERROR: min(max)_f0 parameter inconsistent with sampling frequency.\n",
	    ProgName); 
    error++;
  }
  dstep = ((double)((int)(0.5 + (sample_freq * par->frame_step))))/sample_freq;
  if(dstep != par->frame_step) {
    if(debug_level)
      Fprintf(stderr,
	      "%s: Frame step set to %f to exactly match signal sample rate.\n",
	      ProgName, dstep);
    par->frame_step = dstep;
  }
  if((par->frame_step > 0.1) || (par->frame_step < (1.0/sample_freq))){
    Fprintf(stderr,
	    "%s: ERROR: frame_step parameter must be between [1/sampling rate, 0.1].\n",
	    ProgName); 
    error++;
  }

  return(error);
}
  

/*
 * This function facilitates ESPS range processing.  It sets srec and erec
 * to their parameter/common values unless a range option has been used, in
 * which case it uses the range specification to set srec and erec.  If
 * there is no range option and if start and nan do not appear in the
 * parameter/common file, then srec and erec are set to 1 and LONG_MAX.
 * Get_range assumes that read_params has been called; If a command-line
 * range option (e.g., -r range) has been used, get_range should be called
 * with positive pflag and with rng equal to the range specification.
 */

#ifdef ESPS_INPUT
static void
get_range(srec, erec, rng, pflag, Sflag, hd)
    long	    *srec;	/* starting record */
    long	    *erec;	/* end record */
    char	    *rng;	/* range string from range option */
    int		    pflag;	/* flag for whether -r or -p used */
    int		    Sflag;	/* flag for whether -S used */
    struct header   *hd;	/* input file header */
{
    long common_nan;

    if (!srec || !erec)
    {
	Fprintf(stderr, "get_range: NULL argument.\n");
	return;
    }

    *srec = 1;
    *erec = LONG_MAX;
    if (pflag)
        lrange_switch (rng, srec, erec, 1);
    else if (Sflag)
        trange_switch (rng, hd, srec, erec);
    else {
        if(symtype("start") == ST_INT) {
            *srec = getsym_i("start");
        }
        if(symtype("nan") == ST_INT) {
            common_nan = getsym_i("nan");
            if (common_nan != 0)
                *erec = *srec + common_nan - 1;
        }
    }
}

#endif
