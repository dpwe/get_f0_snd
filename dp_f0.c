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
 * Written by:  David Talkin
 * Checked by:
 * Revised by:  Derek Lin, David Talkin
 *
 * Brief description:  Estimate speech fundamental frequency.
 *
 */

static char *sccs_id = "@(#)dp_f0.c	1.14	10/21/96	ERL";

/* A fundamental frequency estimation algorithm using the normalized
   cross correlation function and dynamic programming.  The algorithm
   implemented here is similar to that presented by B. Secrest and
   G. Doddington, "An integrated pitch tracking algorithm for speech
   systems", Proc. ICASSP-83, pp.1352-1355.  It is fully described
   by D. Talkin, "A robust algorithm for ptich tracking (RAPT)", in
   W. B. Kleijn & K. K. Paliwal (eds.) Speech Coding and Synthesis,
   (New York: Elsevier, 1995). */

/* For each frame, up to par->n_cands cross correlation peaks are
   considered as F0 intervals.  Each is scored according to its within-
   frame properties (relative amplitude, relative location), and
   according to its connectivity with each of the candidates in the
   previous frame.  An unvoiced hypothesis is also generated at each
   frame and is considered in the light of voicing state change cost,
   the quality of the cross correlation peak, and frequency continuity. */

/* At each frame, each candidate has associated with it the following
   items:
	its peak value
	its peak value modified by its within-frame properties
	its location
	the candidate # in the previous frame yielding the min. err.
		(this is the optimum path pointer!)
	its cumulative cost: (local cost + connectivity cost +
		cumulative cost of its best-previous-frame-match). */

/* Dynamic programming is then used to pick the best F0 trajectory and voicing
   state given the local and transition costs for the entire utterance. */

/* To avoid the necessity of computing the full crosscorrelation at
   the input sample rate, the signal is downsampled; a full ccf is
   computed at the lower frequency; interpolation is used to estimate the
   location of the peaks at the higher sample rate; and the fine-grained
   ccf is computed only in the vicinity of these estimated peak
   locations. */

#include <stdio.h>
#include <math.h>
/* #include <malloc.h> */
#include <stdlib.h>
/* #include <esps/esps.h> */
#include <assert.h>
#include "f0.h"
#include "f0_structs.h"

#ifndef FLT_MAX
#define FLT_MAX         ((float) 3.40282347E+38)
#endif

#ifndef spsassert
#define spsassert(ex,ms) {if (!(ex)){fprintf(stderr,"SPS assertion failed: %s\n",(ms));assert(0);exit(1);}}
#endif


extern int  debug_level;
extern char *ProgName;
  
/*
 * READ_SIZE: length of input data frame in sec to read
 * DP_CIRCULAR: determines the initial size of DP circular buffer in sec
 * DP_HIST: stored frame history in second before checking for common path 
 *      DP_CIRCULAR > READ_SIZE, DP_CIRCULAR at least 2 times of DP_HIST 
 * DP_LIMIT: in case no convergence is found, DP frames of DP_LIMIT secs
 *      are kept before output is forced by simply picking the lowest cost
 *      path
 */

#define READ_SIZE 0.2
#define DP_CIRCULAR 1.5
#define DP_HIST 0.5
#define DP_LIMIT 1.0

/* 
 * stationarity parameters -
 * STAT_WSIZE: window size in sec used in measuring frame energy/stationarity
 * STAT_AINT: analysis interval in sec in measuring frame energy/stationarity
 */
#define STAT_WSIZE 0.030
#define STAT_AINT 0.020

/*
 * headF points to current frame in the circular buffer, 
 * tailF points to the frame where tracks start
 * cmpthF points to starting frame of converged path to backtrack
 */

static Frame *headF = NULL, *tailF = NULL, *cmpthF = NULL;

static  int *pcands = NULL;	/* array for backtracking in convergence check */
static int cir_buff_growth_count = 0;

static int size_cir_buffer,	/* # of frames in circular DP buffer */
           size_frame_hist,	/* # of frames required before convergence test */
           size_frame_out,	/* # of frames before forcing output */
           num_active_frames,	/* # of frames from tailF to headF */
           output_buf_size;	/* # of frames allocated to output buffers */

/* 
 * DP parameters
 */
static float tcost, tfact_a, tfact_s, frame_int, vbias, fdouble, wdur, ln2,
             freqwt, lagwt;
static int step, size, nlags, start, stop, ncomp, *locs = NULL;
static short maxpeaks;

static int wReuse = 0;  /* number of windows seen before resued */
static Windstat *windstat;

static float *f0p = NULL, *vuvp = NULL, *rms_speech = NULL, 
             *acpkp = NULL, *peaks = NULL;
static int first_time = 1, pad;


static Stat* get_stationarity();

/*--------------------------------------------------------------------*/
int
get_Nframes(buffsize, pad, step)
    long    buffsize;
    int     pad, step;
{
  if (buffsize < pad)
    return (0);
  else
    return ((buffsize - pad)/step);
}


/*--------------------------------------------------------------------*/
int
init_dp_f0(freq, par, buffsize, sdstep)
    double	freq;
    F0_params	*par;
    long	*buffsize, *sdstep;
{
  int nframes;
  int i;
  int stat_wsize, agap, ind, downpatch;

/*
 * reassigning some constants 
 */

  tcost = par->trans_cost;
  tfact_a = par->trans_amp;
  tfact_s = par->trans_spec;
  vbias = par->voice_bias;
  fdouble = par->double_cost;
  frame_int = par->frame_step;
  
  step = round(frame_int * freq);
  size = round(par->wind_dur * freq);
  frame_int = ((float)step)/freq;
  wdur = ((float)size)/freq;
  start = round(freq / par->max_f0);
  stop = round(freq / par->min_f0);
  nlags = stop - start + 1;
  ncomp = size + stop + 1; /* # of samples required by xcorr
			      comp. per fr. */
  maxpeaks = 2 + (nlags/2);	/* maximum number of "peaks" findable in ccf */
  ln2 = log(2.0);
  size_frame_hist = (int) (DP_HIST / frame_int);
  size_frame_out = (int) (DP_LIMIT / frame_int);

/*
 * SET UP THE D.P. WEIGHTING FACTORS:
 *      The intent is to make the effectiveness of the various fudge factors
 *      independent of frame rate or sampling frequency.                
 */
  
  /* Lag-dependent weighting factor to emphasize early peaks (higher freqs)*/
  lagwt = par->lag_weight/stop;
  
  /* Penalty for a frequency skip in F0 per frame */
  freqwt = par->freq_weight/frame_int;
  
  i = (int) (READ_SIZE *freq);
  if(ncomp >= step) nframes = ((i-ncomp)/step ) + 1;
  else nframes = i / step;

  /* *buffsize is the number of samples needed to make F0 computation
     of nframes DP frames possible.  The last DP frame is patched with
     enough points so that F0 computation on it can be carried.  F0
     computaion on each frame needs enough points to do

     1) xcross or cross correlation measure:
           enough points to do xcross - ncomp

     2) stationarity measure:
           enough to make 30 msec windowing possible - ind

     3) downsampling:
           enough to make filtering possible -- downpatch
 
     So there are nframes whole DP frames, padded with pad points
     to make the last frame F0 computation ok.

  */

  /* last point in data frame needs points of 1/2 downsampler filter length 
     long, 0.005 is the filter length used in downsampler */
  downpatch = (((int) (freq * 0.005))+1) / 2;

  stat_wsize = (int) (STAT_WSIZE * freq);
  agap = (int) (STAT_AINT * freq);
  ind = ( agap - stat_wsize ) / 2;
  i = stat_wsize + ind;
  pad = downpatch + ((i>ncomp) ? i:ncomp);
  *buffsize = nframes * step + pad;
  *sdstep = nframes * step;
  
  /* Allocate space for the DP storage circularly linked data structure */

  size_cir_buffer = (int) (DP_CIRCULAR / frame_int);

  /* creating circularly linked data structures */
  tailF = alloc_frame(nlags, par->n_cands);
  headF = tailF;

  /* link them up */
  for(i=1; i<size_cir_buffer; i++){
    headF->next = alloc_frame(nlags, par->n_cands);
    headF->next->prev = headF;
    headF = headF->next;
  }
  headF->next = tailF;
  tailF->prev = headF;

  headF = tailF;

  /* Allocate sscratch array to use during backtrack convergence test. */
  if( ! pcands ) {
    pcands = (int *) malloc( par->n_cands * sizeof(int));
    spsassert(pcands,"can't allocate pathcands");
  }

  /* Allocate arrays to return F0 and related signals. */

  /* Note: remember to compare *vecsize with size_frame_out, because
     size_cir_buffer is not constant */
  output_buf_size = size_cir_buffer;
  rms_speech = (float*)malloc(sizeof(float) * output_buf_size);
  spsassert(rms_speech,"rms_speech malloc failed");
  f0p = (float*)malloc(sizeof(float) * output_buf_size);
  spsassert(f0p,"f0p malloc failed");
  vuvp = (float*)malloc(sizeof(float)* output_buf_size);
  spsassert(vuvp,"vuvp malloc failed");
  acpkp = (float*)malloc(sizeof(float) * output_buf_size);
  spsassert(acpkp,"acpkp malloc failed");

  /* Allocate space for peak location and amplitude scratch arrays. */
  peaks = (float*)malloc(sizeof(float) * maxpeaks);
  spsassert(peaks,"peaks malloc failed");
  locs = (int*)malloc(sizeof(int) * maxpeaks);
  spsassert(locs, "locs malloc failed");
  
  /* Initialise the retrieval/saving scheme of window statistic measures */
  wReuse = agap / step;
  if (wReuse){
      windstat = (Windstat *) malloc( wReuse * sizeof(Windstat));
      spsassert(windstat, "windstat malloc failed");
      for(i=0; i<wReuse; i++){
	  windstat[i].err = 0;
	  windstat[i].rms = 0;
      }
  }

  if(debug_level){
    Fprintf(stderr, "%s: done with initialization:\n", ProgName);
    Fprintf(stderr,
	    " size_cir_buffer:%d  xcorr frame size:%d start lag:%d nlags:%d\n",
	    size_cir_buffer, size, start, nlags);
  }

  num_active_frames = 0;
  first_time = 1;

  return(0);
}
  

/*--------------------------------------------------------------------*/
int
dp_f0(fdata, buff_size, sdstep, freq,
      par, f0p_pt, vuvp_pt, rms_speech_pt, acpkp_pt, vecsize, last_time)
    float	*fdata;
    int		buff_size, sdstep;
    double	freq;
    F0_params	*par;		/* analysis control parameters */
    float	**f0p_pt, **vuvp_pt, **rms_speech_pt, **acpkp_pt;
    int		*vecsize, last_time;
{
  float  maxval, engref, *sta, *rms_ratio, *dsdata, *downsample();
  register float ttemp, ftemp, ft1, ferr, err, errmin;
  register int  i, j, k, loc1, loc2;
  int   nframes, maxloc, ncand, ncandp, minloc,
        decimate, samsds;

  Stat *stat = NULL;

  nframes = get_Nframes((long) buff_size, pad, step); /* # of whole frames */

  if(debug_level)
    Fprintf(stderr,
	    "%s: ******* Computing %d dp frames ******** from %d points\n",
	    ProgName, nframes, buff_size);

  /* Now downsample the signal for coarse peak estimates. */

  decimate = freq/2000.0;	/* downsample to about 2kHz */
  if (decimate <= 1)
    dsdata = fdata;
  else {
    samsds = ((nframes-1) * step + ncomp) / decimate;
    dsdata = downsample(fdata, buff_size, sdstep, freq, &samsds, decimate, 
			first_time, last_time);
    if (!dsdata) {
      Fprintf(stderr, "%s: can't get downsampled data.\n", ProgName);
      return 1;
    }
  }

  /* Get a function of the "stationarity" of the speech signal. */

  stat = get_stationarity(fdata, freq, buff_size, nframes, step, first_time);
  if (!stat) { 
    Fprintf(stderr, "%s: can't get stationarity\n", ProgName);
    return(1);
  }
  sta = stat->stat;
  rms_ratio = stat->rms_ratio;

  /***********************************************************************/
  /* MAIN FUNDAMENTAL FREQUENCY ESTIMATION LOOP */
  /***********************************************************************/
  if(!first_time && nframes > 0) headF = headF->next;

  for(i = 0; i < nframes; i++) {
 
    /* NOTE: This buffer growth provision is probably not necessary.
       It was put in (with errors) by Derek Lin and apparently never
       tested.  My tests and analysis suggest it is completely
       superfluous. DT 9/5/96 */
    /* Dynamically allocating more space for the circular buffer */
    if(headF == tailF->prev){
      Frame *frm;

      if(cir_buff_growth_count > 5){
	Fprintf(stderr,
		"%s: too many requests (%d) for dynamically allocating space.\n   There may be a problem in finding converged path.\n",
		ProgName, cir_buff_growth_count);
	return(1);
      }
      if(debug_level) 
	Fprintf(stderr, "%s: allocating %d more frames for DP circ. buffer.\n",
		ProgName, size_cir_buffer);
      frm = alloc_frame(nlags, par->n_cands);
      headF->next = frm;
      frm->prev = headF;
      for(k=1; k<size_cir_buffer; k++){
	frm->next = alloc_frame(nlags, par->n_cands);
	frm->next->prev = frm;
	frm = frm->next;
      }
      frm->next = tailF;
      tailF->prev = frm;
      cir_buff_growth_count++;
    }

    headF->rms = stat->rms[i];
    get_fast_cands(fdata, dsdata, i, step, size, decimate, start,
		   nlags, &engref, &maxloc,
		   &maxval, headF->cp, peaks, locs, &ncand, par);
    

    /*    Move the peak value and location arrays into the dp structure */
    {
      register float *ftp1, *ftp2;
      register short *sp1;
      register int *sp2;
      
      for(ftp1 = headF->dp->pvals, ftp2 = peaks,
	  sp1 = headF->dp->locs, sp2 = locs, j=ncand; j--; ) {
	*ftp1++ = *ftp2++;
	*sp1++ = *sp2++;
      }
      *sp1 = -1;		/* distinguish the UNVOICED candidate */
      *ftp1 = maxval;
      headF->dp->mpvals[ncand] = vbias+maxval; /* (high cost if cor. is high)*/
    }

    /* Apply a lag-dependent weight to the peaks to encourage the selection
       of the first major peak.  Translate the modified peak values into
       costs (high peak ==> low cost). */
    for(j=0; j < ncand; j++){
      ftemp = 1.0 - ((float)locs[j] * lagwt);
      headF->dp->mpvals[j] = 1.0 - (peaks[j] * ftemp);
    }
    ncand++;			/* include the unvoiced candidate */
    headF->dp->ncands = ncand;

    /*********************************************************************/
    /*    COMPUTE THE DISTANCE MEASURES AND ACCUMULATE THE COSTS.       */
    /*********************************************************************/

    ncandp = headF->prev->dp->ncands;
    for(k=0; k<ncand; k++){	/* for each of the current candidates... */
      minloc = 0;
      errmin = FLT_MAX;
      if((loc2 = headF->dp->locs[k]) > 0) { /* current cand. is voiced */
	for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */
	  /*    Get cost due to inter-frame period change. */
	  loc1 = headF->prev->dp->locs[j];
	  if (loc1 > 0) { /* prev. was voiced */
	    ftemp = log(((double) loc2) / loc1);
	    ttemp = fabs(ftemp);
	    ft1 = fdouble + fabs(ftemp + ln2);
	    if (ttemp > ft1)
	      ttemp = ft1;
	    ft1 = fdouble + fabs(ftemp - ln2);
	    if (ttemp > ft1)
	      ttemp = ft1;
	    ferr = ttemp * freqwt;
	  } else {		/* prev. was unvoiced */
	    ferr = tcost + (tfact_s * sta[i]) + (tfact_a / rms_ratio[i]);
	  }
	  /*    Add in cumulative cost associated with previous peak. */
	  err = ferr + headF->prev->dp->dpvals[j];
	  if(err < errmin){	/* find min. cost */
	    errmin = err;
	    minloc = j;
	  }
	}
      } else {			/* this is the unvoiced candidate */
	for(j=0; j<ncandp; j++){ /* for each PREVIOUS candidate... */
	  
	  /*    Get voicing transition cost. */
	  if (headF->prev->dp->locs[j] > 0) { /* previous was voiced */
	    ferr = tcost + (tfact_s * sta[i]) + (tfact_a * rms_ratio[i]);
	  }
	  else
	    ferr = 0.0;
	  /*    Add in cumulative cost associated with previous peak. */
	  err = ferr + headF->prev->dp->dpvals[j];
	  if(err < errmin){	/* find min. cost */
	    errmin = err;
	    minloc = j;
	  }
	}
      }
      /* Now have found the best path from this cand. to prev. frame */
      if (first_time && i==0) {		/* this is the first frame */
	headF->dp->dpvals[k] = headF->dp->mpvals[k];
	headF->dp->prept[k] = 0;
      } else {
	headF->dp->dpvals[k] = errmin + headF->dp->mpvals[k];
	headF->dp->prept[k] = minloc;
      }
    } /*    END OF THIS DP FRAME */

    if (i < nframes - 1)
      headF = headF->next;
    
    if (debug_level >= 2) {
      Fprintf(stderr,"%d engref:%10.0f max:%7.5f loc:%4d\n",
	      i,engref,maxval,maxloc);
    }
    
  } /* end for (i ...) */

  /***************************************************************/
  /* DONE WITH FILLING DP STRUCTURES FOR THE SET OF SAMPLED DATA */
  /*    NOW FIND A CONVERGED DP PATH                             */
  /***************************************************************/

  *vecsize = 0;			/* # of output frames returned */

  num_active_frames += nframes;

  if( num_active_frames >= size_frame_hist  || last_time ){
    Frame *frm;
    int  num_paths, best_cand, frmcnt, checkpath_done = 1;
    float patherrmin;
      
    if(debug_level)
      Fprintf(stderr, "%s: available frames for backtracking: %d\n", 
	      ProgName, num_active_frames);
      
    patherrmin = FLT_MAX;
    best_cand = 0;
    num_paths = headF->dp->ncands;

    /* Get the best candidate for the final frame and initialize the
       paths' backpointers. */
    frm = headF;
    for(k=0; k < num_paths; k++) {
      if (patherrmin > headF->dp->dpvals[k]){
	patherrmin = headF->dp->dpvals[k];
	best_cand = k;	/* index indicating the best candidate at a path */
      }
      pcands[k] = frm->dp->prept[k];
    }

    if(last_time){     /* Input data was exhausted. force final outputs. */
      cmpthF = headF;		/* Use the current frame as starting point. */
    } else {
      /* Starting from the most recent frame, trace back each candidate's
	 best path until reaching a common candidate at some past frame. */
      frmcnt = 0;
      while (1) {
	frm = frm->prev;
	frmcnt++;
	checkpath_done = 1;
	for(k=1; k < num_paths; k++){ /* Check for convergence. */
	  if(pcands[0] != pcands[k])
	    checkpath_done = 0;
	}
	if( ! checkpath_done) { /* Prepare for checking at prev. frame. */
	  for(k=0; k < num_paths; k++){
	    pcands[k] = frm->dp->prept[pcands[k]];
	  }
	} else {	/* All paths have converged. */
	  cmpthF = frm;
	  best_cand = pcands[0];
	  if(debug_level)
	    Fprintf(stderr,
		    "%s: paths went back %d frames before converging\n",
		    ProgName, frmcnt);
	  break;
	}
	if(frm == tailF){	/* Used all available data? */
	  if( num_active_frames < size_frame_out) { /* Delay some more? */
	    checkpath_done = 0; /* Yes, don't backtrack at this time. */
	    cmpthF = NULL;
	  } else {		/* No more delay! Force best-guess output. */
	    checkpath_done = 1;
	    cmpthF = headF;
	    Fprintf(stderr,
		    "%s: WARNING: no converging path found after going back %d frames, will use the lowest cost path\n",
		    ProgName, num_active_frames);
	  }
	  break;
	} /* end if (frm ...) */
      }	/* end while (1) */
    } /* end if (last_time) ... else */

    /*************************************************************/
    /* BACKTRACKING FROM cmpthF (best_cand) ALL THE WAY TO tailF    */
    /*************************************************************/
    i = 0;
    frm = cmpthF;	/* Start where convergence was found (or faked). */
    while( frm != tailF->prev && checkpath_done){
      if( i == output_buf_size ){ /* Need more room for outputs? */
	output_buf_size *= 2;
	if(debug_level)
	  Fprintf(stderr,
		  "%s: reallocating space for output frames: %d\n",
		  ProgName, output_buf_size);
	rms_speech = (float *) realloc((char *) rms_speech,
				       sizeof(float) * output_buf_size);
	spsassert(rms_speech, "rms_speech realloc failed in dp_f0()");
	f0p = (float *) realloc((char *) f0p,
				sizeof(float) * output_buf_size);
	spsassert(f0p, "f0p realloc failed in dp_f0()");
	vuvp = (float *) realloc(vuvp, sizeof(float) * output_buf_size);
	spsassert(vuvp, "vuvp realloc failed in dp_f0()");
	acpkp = (float *) realloc(acpkp, sizeof(float) * output_buf_size);
	spsassert(acpkp, "acpkp realloc failed in dp_f0()");
      }
      rms_speech[i] = frm->rms;
      acpkp[i] =  frm->dp->pvals[best_cand];
      loc1 = frm->dp->locs[best_cand];
      vuvp[i] = 1.0;
      best_cand = frm->dp->prept[best_cand];
      ftemp = loc1;
      if(loc1 > 0) {		/* Was f0 actually estimated for this frame? */
	if (loc1 > start && loc1 < stop) { /* loc1 must be a local maximum. */
	  float cormax, cprev, cnext, den;
		  
	  j = loc1 - start;
	  cormax = frm->cp->correl[j];
	  cprev = frm->cp->correl[j+1];
	  cnext = frm->cp->correl[j-1];
	  den = (2.0 * ( cprev + cnext - (2.0 * cormax) ));
	  /*
	   * Only parabolic interpolate if cormax is indeed a local 
	   * turning point. Find peak of curve that goes though the 3 points
	   */
		  
	  if (fabs(den) > 0.000001)
	    ftemp += 2.0 - ((((5.0*cprev)+(3.0*cnext)-(8.0*cormax))/den));
	}
	f0p[i] = freq/ftemp;
      } else {		/* No valid estimate; just fake some arbitrary F0. */
	f0p[i] = 0;
	vuvp[i] = 0.0;
      }
      frm = frm->prev;
	  
      if (debug_level >= 2)
	Fprintf(stderr," i:%4d%8.1f%8.1f\n",i,f0p[i],vuvp[i]);
      /* f0p[i] starts from the most recent one */ 
      /* Need to reverse the order in the calling function */
      i++;
    } /* end while() */
    if (checkpath_done){
      *vecsize = i;
      tailF = cmpthF->next;
      num_active_frames -= *vecsize;
    }
  } /* end if() */

  if (debug_level)
    Fprintf(stderr, "%s: writing out %d frames.\n", ProgName, *vecsize);
  
  *f0p_pt = f0p;
  *vuvp_pt = vuvp;
  *acpkp_pt = acpkp;
  *rms_speech_pt = rms_speech;
  *acpkp_pt = acpkp;
  
  if(first_time) first_time = 0;
  return(0);
}


/*--------------------------------------------------------------------*/
Frame *
alloc_frame(nlags, ncands)
    int     nlags, ncands;
{
  Frame *frm;
  int j;

  frm = (Frame*)malloc(sizeof(Frame));
  frm->dp = (Dprec *) malloc(sizeof(Dprec));
  spsassert(frm->dp,"frm->dp malloc failed in alloc_frame");
  frm->dp->ncands = 0;
  frm->cp = (Cross *) malloc(sizeof(Cross));
  spsassert(frm->cp,"frm->cp malloc failed in alloc_frame");
  frm->cp->correl = (float *) malloc(sizeof(float) * nlags);
  spsassert(frm->cp->correl, "frm->cp->correl malloc failed");
  /* Allocate space for candidates and working arrays. */
  frm->dp->locs = (short*)malloc(sizeof(short) * ncands);
  spsassert(frm->dp->locs,"frm->dp->locs malloc failed in alloc_frame()");
  frm->dp->pvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->pvals,"frm->dp->pvals malloc failed in alloc_frame()");
  frm->dp->mpvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->mpvals,"frm->dp->mpvals malloc failed in alloc_frame()");
  frm->dp->prept = (short*)malloc(sizeof(short) * ncands);
  spsassert(frm->dp->prept,"frm->dp->prept malloc failed in alloc_frame()");
  frm->dp->dpvals = (float*)malloc(sizeof(float) * ncands);
  spsassert(frm->dp->dpvals,"frm->dp->dpvals malloc failed in alloc_frame()");
    
  /*  Initialize the cumulative DP costs to zero */
  for(j = ncands-1; j >= 0; j--)
    frm->dp->dpvals[j] = 0.0;

  return(frm);
}


/*--------------------------------------------------------------------*/
/* push window stat to stack, and pop the oldest one */

static int
save_windstat(rho, order, err, rms)
    float   *rho;
    int     order;
    float   err;
    float   rms;
{
    int i,j;

    if(wReuse > 1){               /* push down the stack */
	for(j=1; j<wReuse; j++){
	    for(i=0;i<=order; i++) windstat[j-1].rho[i] = windstat[j].rho[i];
	    windstat[j-1].err = windstat[j].err;
	    windstat[j-1].rms = windstat[j].rms;
	}
	for(i=0;i<=order; i++) windstat[wReuse-1].rho[i] = rho[i]; /*save*/
	windstat[wReuse-1].err = err;
	windstat[wReuse-1].rms = rms;
	return 1;
    } else if (wReuse == 1) {
	for(i=0;i<=order; i++) windstat[0].rho[i] = rho[i];  /* save */
	windstat[0].err = err;
	windstat[0].rms = rms;
	return 1;
    } else 
	return 0;
}


/*--------------------------------------------------------------------*/
static int
retrieve_windstat(rho, order, err, rms)
    float   *rho;
    int     order;
    float   *err;
    float   *rms;
{
    Windstat wstat;
    int i;
	
    if(wReuse){
	wstat = windstat[0];
	for(i=0; i<=order; i++) rho[i] = wstat.rho[i];
	*err = wstat.err;
	*rms = wstat.rms;
	return 1;
    }
    else return 0;
}


/*--------------------------------------------------------------------*/
static float
get_similarity(order, size, pdata, cdata,
	       rmsa, rms_ratio, pre, stab, w_type, init)
    int     order, size;
    float   *pdata, *cdata;
    float   *rmsa, *rms_ratio, pre, stab;
    int     w_type, init;
{
  float rho3[BIGSORD+1], err3, rms3, rmsd3, b0, t, a2[BIGSORD+1], 
      rho1[BIGSORD+1], a1[BIGSORD+1], b[BIGSORD+1], err1, rms1, rmsd1;
  float itakura(), wind_energy();

/* (In the lpc() calls below, size-1 is used, since the windowing and
   preemphasis function assumes an extra point is available in the
   input data array.  This condition is apparently no longer met after
   Derek's modifications.) */

  /* get current window stat */
  lpc(order, stab, size-1, cdata,
      a2, rho3, (float *) NULL, &err3, &rmsd3, pre, w_type);
  rms3 = wind_energy(cdata, size, w_type);
  
  if(!init) {
      /* get previous window stat */
      if( !retrieve_windstat(rho1, order, &err1, &rms1)){
	  lpc(order, stab, size-1, pdata,
	      a1, rho1, (float *) NULL, &err1, &rmsd1, pre, w_type);
	  rms1 = wind_energy(pdata, size, w_type);
      }
      a_to_aca(a2+1,b,&b0,order);
      t = itakura(order,b,&b0,rho1+1,&err1) - .8;
      if(rms1 > 0.0)
	  *rms_ratio = (0.001 + rms3)/rms1;
      else
	  if(rms3 > 0.0)
	      *rms_ratio = 2.0;	/* indicate some energy increase */
	  else
	      *rms_ratio = 1.0;	/* no change */
  } else {
      *rms_ratio = 1.0;
      t = 10.0;
  }
  *rmsa = rms3;
  save_windstat( rho3, order, err3, rms3);
  return((float)(0.2/t));
}


/* -------------------------------------------------------------------- */
/* This is an ad hoc signal stationarity function based on Itakura
 * distance and relative amplitudes.
 */
/* 
  This illustrates the window locations when the very first frame is read.
  It shows an example where each frame step |  .  | is 10 msec.  The
  frame step size is variable.  The window size is always 30 msec.
  The window centers '*' is always 20 msec apart.
  The windows cross each other right at the center of the DP frame, or
  where the '.' is.

                          ---------*---------   current window

              ---------*---------  previous window

  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |  .  |
              ^           ^  ^
              ^           ^  ^
              ^           ^  fdata
              ^           ^
              ^           q
	      p

                          ---
                          ind

  fdata, q, p, ind, are variables used below.
   
*/

static Stat*
get_stationarity(fdata, freq, buff_size, nframes, frame_step, first_time)
    float   *fdata;
    double  freq;
    int     buff_size, nframes, frame_step, first_time;
{
  static Stat *stat;
  static int nframes_old = 0, memsize;
  static float *mem;
  float preemp = 0.4, stab = 30.0;
  float *p, *q, *r, *datend;
  int ind, i, j, m, size, order, agap, w_type = 3;

  agap = (int) (STAT_AINT *freq);
  size = (int) (STAT_WSIZE * freq);
  ind = (agap - size) / 2;

  if( nframes_old < nframes || !stat || first_time){
    /* move this to init_dp_f0() later */
    nframes_old = nframes;
    if(stat){ 
      free((char *) stat->stat);
      free((char *) stat->rms);
      free((char *) stat->rms_ratio);
      free((char *) stat);
    }
    stat = (Stat *) malloc(nframes *sizeof(Stat));
    spsassert(stat,"stat malloc failed in get_stationarity");
    stat->stat = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->stat,"stat->stat malloc failed in get_stationarity");
    stat->rms = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->rms,"stat->rms malloc failed in get_stationarity");
    stat->rms_ratio = (float*)malloc(sizeof(float)*nframes);
    spsassert(stat->rms_ratio,"stat->ratio malloc failed in get_stationarity");
    memsize = (int) (STAT_WSIZE * freq) + (int) (STAT_AINT * freq);
    mem = (float *) malloc( sizeof(float) * memsize);
    spsassert(mem, "mem malloc failed in get_stationarity()");
    for(j=0; j<memsize; j++) mem[j] = 0;
  }
  
  if(nframes == 0) return(stat);

  q = fdata + ind;
  datend = fdata + buff_size;

  if((order = 2.0 + (freq/1000.0)) > BIGSORD) {
    Fprintf(stderr,
	    "%s: Optimim order (%d) exceeds that allowable (%d); reduce Fs\n",
	    ProgName, order, BIGSORD);
    order = BIGSORD;
  }

  /* prepare for the first frame */
  for(j=memsize/2, i=0; j<memsize; j++, i++) mem[j] = fdata[i];

  /* never run over end of frame, should already taken care of when read */

  for(j=0, p = q - agap; j < nframes; j++, p += frame_step, q += frame_step){
      if( (p >= fdata) && (q >= fdata) && ( q + size <= datend) )
	  stat->stat[j] = get_similarity(order,size, p, q, 
					     &(stat->rms[j]),
					     &(stat->rms_ratio[j]),preemp,
					     stab,w_type, 0);
      else {
	  if(first_time) {
	      if( (p < fdata) && (q >= fdata) && (q+size <=datend) )
		  stat->stat[j] = get_similarity(order,size, NULL, q,
						     &(stat->rms[j]),
						     &(stat->rms_ratio[j]),
						     preemp,stab,w_type, 1);
	      else{
		  stat->rms[j] = 0.0;
		  stat->stat[j] = 0.01 * 0.2;   /* a big transition */
		  stat->rms_ratio[j] = 1.0;   /* no amplitude change */
	      }
	  } else {
	      if( (p<fdata) && (q+size <=datend) ){
		  stat->stat[j] = get_similarity(order,size, mem, 
						     mem + (memsize/2) + ind,
						     &(stat->rms[j]),
						     &(stat->rms_ratio[j]),
						     preemp, stab,w_type, 0);
		  /* prepare for the next frame_step if needed */
		  if(p + frame_step < fdata ){
		      for( m=0; m<(memsize-frame_step); m++) 
			  mem[m] = mem[m+frame_step];
		      r = q + size;
		      for( m=0; m<frame_step; m++) 
			  mem[memsize-frame_step+m] = *r++;
		  }
	      }
	  }
      }
  }

  /* last frame, prepare for next call */
  for(j=(memsize/2)-1, p=fdata + (nframes * frame_step)-1; j>=0; j-- )
    mem[j] = *p--;
  return(stat);
}


/* -------------------------------------------------------------------- */
/*	Round the argument to the nearest integer.			*/

#ifdef NOTDEF  /* use built-in round on modern systems */
int
round(flnum)
    double  flnum;
{
  return((flnum >= 0.0) ? (int)(flnum + 0.5) : (int)(flnum - 0.5));
}
#endif
