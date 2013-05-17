/* dpwe_add.h */

#include <snd.h>

#define ST_UNDEF 0
#define ST_STRING 1
#define ST_INT 2
#define ST_FLOAT 3

#define SC_NOCOMMON 1

#define FT_FEA 2
#define FEA_SD 3

#define NONE 0

#define YES (1)
#define NO (0)
//#define TRUE (1)
//#define FALSE (0)

#define DOUBLE (8)
#define FLOAT (4)

#define INT_MAX 0x7FFFFFFF
#define LONG_MAX 0x7FFFFFFFFFFFFFFF

typedef struct commonheader {
    char *prog;
    char *vers;
    char *progdate;
    int tag;
} COMMONHEADER;

typedef struct header {
    struct commonheader common;
    SOUND *snd;
} HEADER;

typedef struct fea_data {
} FEA_DATA;

typedef struct feasd {
    void *data;
} FEASD;

struct feasd *sd_rec;
struct fea_data *fea_rec;

int read_params(char *filename, int flag, void *ptr);
double getsym_d(char *name);
int getsym_i(char *name);
char *getsym_s(char *name);
int symtype(char *name);
// Open a data file of specified type and format
char *eopen(char *ProgName, char *namein, char *mode, int content_type, int format_type, struct header **hdr, FILE **file);
// Read a named header field from a header record
// "record_freq", "start_time", ...
double get_genhd_val(char *name, struct header *ihd, double dflt);


// Create a new data file header
struct header *new_header(int type);
// Add a named feature record to the data file header
void add_fea_fld(char *name, long chans, int size, long *val, int type, char **buf, struct header *hdr);
// Record the name of the file that this record is based on
void add_source_file(struct header *ohdr, const char *ifname, struct header *ihdr);
// Add a free-text comment to the header structure
void add_comment(struct header *ohdr, const char *cmt);
// Add double to header struct
void add_genhd_d(const char* name, double *pval, int len, struct header *ohdr);
// Add float to header struct
void add_genhd_f(const char* name, float *pval, int len, struct header *ohdr);
// Add short to header struct
void add_genhd_s(const char* name, short *pval, int len, struct header *ohdr);




struct fea_data *allo_fea_rec(struct header *hdr);
void *get_fea_ptr(struct fea_data *fea_rec, char *name, struct header *hdr);
struct feasd *allo_feasd_recs(struct header *hdr, int TYPE, size_t buff_size, void *ptr, int flag);

// from link errors
/* 
USAGE
*/
