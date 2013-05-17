/* dpwe_add.c 
 *
 * Stuff to add to get_f0 to make it compile stand-alone.
 *
 * Dan Ellis dpwe@ee.columbia.edu 2012-09-19
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dpwe_add.h"

typedef struct keyval {
    char *key;
    char *value;
} keyval_t;

#define NUMPARAMS 1024

static keyval_t params[NUMPARAMS];
static int nparams = 0;

#define BUFLEN 1024

int find_ws(char *s)
{
    int i = 0;
    /* if no ws find, return pointer to EOS */
    while(s[i] != ' ' && s[i] != '\t' && s[i] != '\n' && s[i] != '\0') {
	/* if (s[i] == 0) { return -1; } */ /* .. or return special? */
	++i;
    }
    return i;
}

int find_nonws(char *s)
{
    int i = 0;
    /* if EOS hit before finding nonws, return pointer to EOS */
    while( (s[i] == ' ' || s[i] == '\t' || s[i] == '\n') && s[i] != '\0') {
	++i;
    }
    return i;
}

/* read a param file, return number of parameters read */
int read_params(char *filename, int flag, void *ptr)
{
    /* read in a name val param file and store in private structure */
    FILE *fp;
    char *s;
    char buf[BUFLEN];
    int nparams;

    /* null or empty filename means do nothing */

    if (filename && filename[0] != '\0') {

	/* flush the buffer (should free strings?) */
	nparams = 0;

	fp = fopen(filename, "r");
	if (fp == NULL) {
	    fprintf(stderr, "Could not read param file %s\n", filename);
	    exit(1);
	}
	while (!feof(fp)) {
	    /* get a line */
	    s = fgets(buf, BUFLEN, fp);
	    if (s != NULL) {
		/* skip any leading blank */
		s += find_nonws(s);
		size_t firstws = find_ws(s);
		if (firstws > 0) {
		    params[nparams].key = strndup(s,firstws);
		    ++nparams;
		    assert(nparams < NUMPARAMS);
		    s += firstws;
		    size_t firstnonws = find_nonws(s);
		    params[nparams].value = strdup(s+firstnonws);
		}
	    }
	}
    }
    return nparams;
}

double getsym_d(char *name)
{
    int i;
    double rslt = 0;

    for (i=0; i<nparams; ++i) {
	if (strcmp(name, params[nparams].key) == 0) {
	    rslt = atof(params[nparams].value);
	    break;
	}
    }
    return rslt;
}

int symtype(char *name)
{
    int i;
    int rslt = ST_UNDEF;

    for (i=0; i<nparams; ++i) {
	if (strcmp(name, params[nparams].key) == 0) {
	    rslt = ST_STRING;
	    break;
	}
    }
    return rslt;
}

int getsym_i(char *name)
{
    int i;
    int rslt = 0;

    for (i=0; i<nparams; ++i) {
	if (strcmp(name, params[nparams].key) == 0) {
	    rslt = atoi(params[nparams].value);
	    break;
	}
    }
    return rslt;
}

char *getsym_s(char *name)
{
    int i;
    char *rslt = 0;

    for (i=0; i<nparams; ++i) {
	if (strcmp(name, params[nparams].key) == 0) {
	    rslt = strdup(params[nparams].value);
	    break;
	}
    }
    return rslt;
}

char *eopen(char *ProgName, char *namein, char *mode, int flag1, int flag2, struct header **hdr, FILE **ofile) {
    return namein;
}

double get_genhd_val(char *name, struct header *ihd, double dflt) {
    return dflt;
}

struct header *new_header(int type) {
    return NULL;
}

void add_fea_fld(char *name, long len, int val, long *dflt, int type, char **charpp, struct header *hdr) {
}

struct fea_data *allo_fea_rec(struct header *hdr) {
    return NULL;
}

void *get_fea_ptr(struct fea_data *fea_rec, char *name, struct header *hdr) {
    return NULL;
}

struct feasd *allo_feasd_recs(struct header *hdr, int TYPE, size_t buff_size, void *ptr, int flag) {
    return NULL;
}

// Add a free-text comment to the header structure
void add_comment(struct header *ohdr, const char *cmt) {
}

//USAGE
