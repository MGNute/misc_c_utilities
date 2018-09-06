//
// Created by miken on 7/18/2018.
//

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <ctype.h>


typedef struct {
    int pos;
    char *name;
} listentry_t;


char *usage ="";
int l_entry_compare_search(const void *a, const void *b) {

    char *a_le;

    listentry_t *b_le;
    a_le = (char*) (a);
//    printf("a = %X   ",a);
//    printf("\ta_le = %s",a_le);

    b_le = *(listentry_t **) (b);
//    printf("\t\tb_le->name: %s\t\t%X\t%X\n",b_le->name,b_le,b_le->name);
//    printf("\t\tcomparison: %d\n",strcmp(a_le, b_le->name));
    return strcmp(a_le, b_le->name);
}


int l_entry_compare(const void *a, const void *b)
{
//    int j;
//    printf("comparing: %X - %X\n",a,b);
    listentry_t* a_le;
    listentry_t* b_le;
    a_le = *(listentry_t**)(a);
    b_le = *(listentry_t**)(b);
//    printf("            %X - %X\n",*a_le, *b_le);

    // debugging
//    for (j=0; j<strlen(a_le->name);j++) {
//        if (a_le->name[j]==10) {
//            break;
//        } else {
//            printf("%c",a_le->name[j]);
//        }
//    }
//    printf("  --  ");
//    for (j=0; j<strlen(b_le->name);j++) {
//        if (b_le->name[j]==10) {
//            break;
//        } else {
//            printf("%c",b_le->name[j]);
//        }
//    }
//    printf("\t%s\t%s\t*** %d ***\n",a_le->name,b_le->name,strcmp(a_le->name, b_le->name));


    return strcmp(a_le->name, b_le->name);
}

int main2(int argc, char *argv[]) {
    char str1[15];
    char str2[15];

    int i;
//    char** carr = (char**)malloc(5 * sizeof(char*));
    listentry_t** larr = (listentry_t**)malloc(30 * sizeof(listentry_t*));

    for (i=0; i<30; i++) {
//        carr[i] = (char*)malloc(200);
        larr[i] = malloc(sizeof(listentry_t));
        larr[i]->name = (char*)malloc(200);
//        fgets(carr[i],199,stdin);
        fgets(larr[i]->name,199,stdin);
        larr[i]->pos = i;
    }
    int pos;
    char *myc;
    myc = strchr(larr[0]->name,10);
    printf("myc = %d",*myc);
    printf("myc = %d",myc-larr[0]->name);
    pos = myc-larr[0]->name;

    printf("\n");
    printf("%s",larr[0]->name);
    larr[0]->name[pos] = 0;
    printf("%s",larr[0]->name);

//    printf("\n\n");
//    for (i=0; i<30; i++) {
//        printf("%s",larr[i]->name);
//    }
//    printf("\npointers: larr[i]\n");
//    for (i=0; i<6; i++) {
//        printf("%d: %X\t%d: %X\t%d: %X\t%d: %X\t%d: %X\n",0+i,larr[0+i],6+i,larr[6+i],12+i,larr[12+i],18+i,larr[18+i],24+i,larr[24+i]);
//    }
//    printf("\n\naddresses: &larr[i]\n");
//    for (i=0; i<6; i++) {
//        printf("%d: %X\t%d: %X\t%d: %X\t%d: %X\t%d: %X\n",0+i,&larr[0+i],6+i,&larr[6+i],12+i,&larr[12+i],18+i,&larr[18+i],24+i,&larr[24+i]);
//    }
//    printf("\n\n");


//    printf("%d\n %s %s",l_entry_compare((void*)&larr[4],(void*)&larr[5]),larr[4]->name, larr[5]->name);
//    printf("sizeof(listentry_t): %d\n", sizeof(listentry_t));
//    printf("sizeof(listentry_t*): %d\n", sizeof(listentry_t*));
//    printf("sizeof(larr[1]): %d\n", sizeof(larr[1]));
//    printf("sizeof((void*)larr[1]): %d\n\n", sizeof((void*)larr[1]));

//    printf("larr:\t%X\n",larr);
//    printf("*larr:\t%X\n",*larr);
//    printf("larr[0]:\t%X\n",larr[0]);
//    printf("&larr[0]:\t%X\n",&larr[0]);
//    printf("larr[1]:\t%X\n",larr[1]);
//    printf("&larr[1]:\t%X\n",&larr[1]);



    qsort(&larr[0], 30, sizeof(listentry_t*),l_entry_compare);
    printf("done sorting...\n\n");
//    qsort(&larr[0], 30, sizeof(listentry_t*),l_entry_compare);


    char *qseq = (char*)malloc(20);
    qseq = "dnstrk\n";
    printf("qseq = %X   %s\n",qseq,qseq);
    void *bsres;
    bsres = bsearch(qseq,&larr[0],30, sizeof(listentry_t*),l_entry_compare_search);
    printf("bsres = %X\n",bsres);
    listentry_t *newle;
    newle = *(listentry_t**)bsres;
    printf("name: %s pos: %d\n",newle->name,newle->pos);


//    printf("\n\n");
//    for (i=0; i<30; i++) {
//        printf("%s",larr[i]->name);
//    }



//    for (i=0; i<5; i++) {
//        printf("%s\n",carr[i]);
//        printf("name: %s\t\ti: %d\n",larr[i]->name, larr[i]->pos);
//    }


//    strcpy(str1, "b");
//    strcpy(str2, "c");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "a");
//    strcpy(str2, "A");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "bc");
//    strcpy(str2, "ac");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "bb");
//    strcpy(str2, "ac");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "ab");
//    strcpy(str2, "ac");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "abc");
//    strcpy(str2, "acb");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);
//
//    strcpy(str1, "a");
//    strcpy(str2, "Z");
//    printf("%d\t%s vs %s\n",strcmp(str1,str2), str1, str2);

    return(0);
}


//    printf("runtot: %f\t\trunct: %d\n",runtot, runct);

//   i=22;
//   j=7;
//   pd = (double)i / (double)j;
//   printf("%f\n",pd);

//    char *a = (char*)malloc(10);
//    a = "-";
//    if (a[0]==45) printf("a is 45\n");
//    printf("a = %d\n",a[0]);

//    printf("args:\n");
//    for (i=0; i < argc; i++) {
//        printf("arg %d:\t%s\n",i,argv[i]);
//    }
//    if (strcmp(argv[1],"max")==0) {
//        printf("max\n");
//    } else {
//        printf("not max\n");
//    }
//    for (i = 0; i < aln->nSeq; i++)
//    {
//        printf("\t%s\n",aln->names[i]);
//    }

int main(int argc, char *argv[]) {

}