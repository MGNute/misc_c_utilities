//
// Created by miken on 7/20/2018.
//

#ifndef MISC_C_UTILITIES_ALIGNMENT_H
#define MISC_C_UTILITIES_ALIGNMENT_H

//
//  Struct Definitions and Prototypes
//

typedef struct {
    int pos;
    char *name;
} listentry_t;

typedef struct {
    int nPos;
    int nSeq;
    char **names;
    char **seqs;
    int nSaved; /* actual allocated size of names and seqs */
    listentry_t **larr;
} alignment_t;

// These three functions were copied and pasted from the source code for FastTree
//   by Morgan Price
alignment_t *ReadAlignment(/*IN*/FILE *fp);
void *myfree(void *p, size_t sz);
void *myrealloc(void *data, size_t szOld, size_t szNew, bool bCopy);
void *mymemdup(void *data, size_t sz);
void *mymalloc(size_t sz);
int alignment_get_seq_position_by_name(alignment_t *aln, char *seqname);
int l_entry_compare_qsort(const void *a, const void *b);
int l_entry_compare_bsearch(const void *a, const void *b);

#endif //MISC_C_UTILITIES_ALIGNMENT_H
