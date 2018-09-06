//
// Created by miken on 7/20/2018.
//

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>
#include "alignment.h"



int main(int argc, char *argv[])
{
    char *aln_fi;
    int i;
    aln_fi = "pdist_test/subset_aln_0_n50.fasta";
    char *init_qseq1;
    char *init_qseq2;
    char **qseq = (char**)malloc(10* sizeof(char*));
//    qseq = "s_14645";
    qseq[0]="s_14645";
    qseq[1]="s_9624";
    qseq[2]="s_16368";
    qseq[3]="s_13922";
    qseq[4]="s_4681";
    qseq[5]="s_11371";
    qseq[6]="s_8151";
    qseq[7]="s_11128";
    qseq[8]="s_11405";
    qseq[9]="s_16904";

    init_qseq1 = "s_14645";
    init_qseq2 = qseq[2];

    alignment_t *test_aln;
    FILE *test_aln_f = fopen(aln_fi,"r");
    test_aln = ReadAlignment(test_aln_f);
    int test_pos;

    test_pos = alignment_get_seq_position_by_name(test_aln,init_qseq1);
    printf("name: %s\t\tposition: %d\n",init_qseq1,test_pos);
    test_pos = alignment_get_seq_position_by_name(test_aln,init_qseq2);
    printf("name: %s\t\tposition: %d\n",init_qseq2,test_pos);

    printf("init_qseq1: %X\tinit_qseq2: %X\t\n",init_qseq1,init_qseq2);
    printf("qseq[0]: %X\tqseq[1]: %X\tqseq[2]: %X\n",qseq[0],qseq[1],qseq[2]);
    printf("qseq %X\n",0,qseq);

    test_pos=0;
    for (i=0; i<10; i++) {
//        printf("qseq[%d] %s\t%X\n",i,qseq[i])  ;//, qseq[i]);
        test_pos = alignment_get_seq_position_by_name(test_aln,qseq[i]);
        printf("name: %s\t\tposition: %d\n",qseq[i],test_pos);
    }
    qseq[i];
//    test_pos = alignment_get_seq_position_by_name(test_aln,qseq);
//    printf("sequence name: %s\t\tposition: %d\n",qseq,test_pos);
    return 0;
}