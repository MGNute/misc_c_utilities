#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/* 
    This program is a quick and dirty implementation of the needleman-wunsch algorithm. Currently it only
    returns the best score rather than the implied alignment, although that can be added by retracing
    the steps back through the matrix F. 
 */

 void inplace_reverse(char * str)
{
  if (str)
  {
    char * end = str + strlen(str) - 1;

    // swap the values in the two given variables
    // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a,b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
        } while (0)

    // walk inwards from both ends of the string, 
    // swapping until we get to the middle
    while (str < end)
    {
      XOR_SWAP(*str, *end);
      str++;
      end--;
    }
#   undef XOR_SWAP
  }
}
void zeroLongLongArray(long long *arr, long dim)
{
    int i;
    for (i=0; i<dim; i++) arr[i]=0;
}

/*
 * Function:    needleScore
 * Description: Computes the score of the best Needleman-Wunsch alignment and populates a DP matrix
 *              along with a matrix of directions of travel through the DP matrix.
 * */
long long needleScore(int len1, char* seq1, int len2, char* seq2, long long* F, long long d, long long m, long long g, 
						long long a, int* dirmat)
{
    // assumes that Fmat is big enough to handle the scores
    // d = difference score, m = match score, g = gap
    int i, j, p, r, c;
    long long mat, del, ins, tmp;
	
	int dir, dirtmp;
	// char *ast = "*";
    
    // fill the edges of F:
    for (i=0; i<len1+1;i++)
    {
        F[i*(len2+1)]=i*g; //should this be i*(len2+1)+1?
    }
    for (j=0; j<len2+1;j++)
    {
        F[j]=j*g;
    }
	// fill the middle of F:
    for (i=1; i<len1+1; i++)
    {
        for (j=1; j<len2+1; j++)
        {
            mat = F[(i-1)*(len2+1)+(j-1)] +
                  (seq1[i-1]==42 || seq2[j-1]==42 ? a : (seq1[i-1]==seq2[j-1] ? m : d)); // (dir = 1)
            del = F[(i-1)*(len2+1)+(j)] + g; // from above (dir = 2)
            ins = F[(i)*(len2+1)+(j-1)] + g; // from the left (dir = 0)
            tmp = (mat >= del ? mat : del);
			dirtmp = (mat >= del ? 1 : 2);
            F[i*(len2+1)+(j)] = (tmp > ins ? tmp : ins);
			dir = (tmp > ins ? dirtmp : 0);
			dirmat[i*(len2+1)+(j)] = dir;
        }
    }

    return F[(len2+1)*(len1+1)-1];
}

/*
 * Function:    getNeedleAlignment
 * Description: Computes the acutal Needleman-Wunsch alignment for two sequences, given
 *              penalties d/g/m/a. Aligned strings are placed in pre-allocated pointers.
 *
 * */
long long getNeedleAlignment(int len1, char* seq1, int len2, char* seq2, long long* F, long long d, long long m, long long g, 
						long long a, char* aln1, char* aln2, int* dirmat, int *aln_len, int *gapct)
{
	int i, j, p, r, c, n_gaps;
    long long mat, del, ins, tmp, best_score;
	int dir, dirtmp;
	
	best_score = needleScore(len1, seq1, len2, seq2, F, d, m, g, a, dirmat);
	
	p = 0;   n_gaps = 0;
	r = len1; c = len2;
	// get the sequences
	while (r > 0 || c > 0)
	{
		dir = dirmat[r*(len2+1)+c];
		if (r==0) dir = 0;
		if (c==0) dir = 2;
		if (dir == 0)
		{
			aln1[p] = 45;
			aln2[p] = seq2[c-1];
			c--;
			n_gaps++;
		} else if (dir == 2) {
			aln1[p] = seq1[r-1];
			aln2[p] = 45;
			r--;
			n_gaps++;
		} else {
			aln1[p] = seq1[r-1];
			aln2[p] = seq2[c-1];
			r--;
			c--;
		}
		p++;
	}

	inplace_reverse(aln1);
	inplace_reverse(aln2);
//    printf("aln1 len: %d\n", strlen(aln1));
//    printf("aln2 len: %d\n", strlen(aln2));
	aln1[p] = 0;
	aln2[p] = 0;
    aln_len[0] = p;
    gapct[0]=n_gaps;
	return best_score;
}

/*
 * Function:    intArrayMax
 * Desc:        Computes the max value of an integer array. Simple loop implementation.
 *
 * Inputs:      intarr: Array of integers to calc the max of.
 *              len:    length of the array.
 *
 * Returns:     maximum value of the array intarr.
 * */
int intArrayMax(int* intarr, int len)
{
    int i, tmx;
    tmx = -1;
    for (i=0; i < len; i++)
    {
        tmx = (tmx > intarr[i] ? tmx : intarr[i]);
    }
    return tmx;
}

void print_f_matrix(long long *F, int l1, int l2)
{
	printf("\n\n");
	int r, c;
	for (r=0; r< l1+1; r++) 
	{
		for (c=0; c< l2+1; c++) 
		{
			printf("%lld, ", F[r*(l2+1)+c]);
		}
		printf("\n");
	}
	printf("\n");
}

void print_help()
{
    printf("help menu: (not implemented yet)\n");
}

int main(int argc, char *argv[])
{
    /*
     * There are two usage options, one that does a big batch of these with files in a particular structure
     * (to avoid counting sequence lengths and the like).
     *
     * To use this option, the first argument must be "-f"
     *
     *  usage: -f "<num_seqs_file1> <seq_len_file1> <sequences_file1> ...(same for 2)..."
        where:
            <num_seqs_file1> is the number of sequences in the first set of files.
            <seq_len_file1> is a binary file with consecutive 32-bit integers representing the
                length of each sequence.
            <sequences_file1> is the path to a file with the sequences packed together via
                concatenation and no other characters anywhere in it.
     
     *
     * The other just takes two sequences from STDIN (each newline terminated) and prints the score back to 
     * STDOUT. This one requires no command line arguments so as long as the user provides more than two
     * arguments we assume the first one is what they want.
     *
     *  other options:
     *      -p  Prints output in a more readable form. If not set, output is tab-delimited (for easy
     *          machine reading).
 *          -3seq    Takes four sequences from STDIN instead of just two. Compares each of the first
     *          three against the fourth. Additionally returns which of thre three had the best score,
     *          along with the relevant stats.
     *      Cost Score Options:
     *      -d  Difference  (default: -1)
     *      -m  Match       (default: 1)
     *      -g  Gap         (default: -2)
     *      -a  Asterix (unknown character) (default: 0)
     *
     *
     */
    int use_case, pretty_output, i, j;
    use_case = 0;   // 0 --> 2 sequences, newline separated, from stdin
                    // 1 --> "-f" option. Several sequences packed together in a file
    pretty_output = 0; // set to 1 if the "-p" flag is used. Makes output readable.

    // COST SCORES
    int d, m, g, a;
    long long res;
    // d = mismatch penalty (*D*ifferent), m = Match score, g = Gap penalty.
    d = -1; m = 1; g = -2; a = 0;


    for (i=0; i<argc; i++) {
        if (strcmp(argv[i],"-h")==0) {
            print_help();
        } else if (strcmp(argv[i],"-p")==0) {
            pretty_output = 1;
        } else if (strcmp(argv[i],"-f")==0) {
            use_case = 1;
        } else if (strcmp(argv[i],"-3seq")==0) {
            use_case = 2;
        } else if (strcmp(argv[i],"-d")==0) {
            if (argc>i+1) {
                d = atoi(argv[i+1]);
            }
        } else if (strcmp(argv[i],"-m")==0) {
            if (argc>i+1) {
                m = atoi(argv[i+1]);
            }
        } else if (strcmp(argv[i],"-g")==0) {
            if (argc>i+1) {
                g = atoi(argv[i+1]);
            }
        } else if (strcmp(argv[i],"-a")==0) {
            if (argc>i+1) {
                a = atoi(argv[i+1]);
            }
        }

    }




        
    if (use_case==1) {
        char* ns1_str = argv[1];
        char* sl1_fn = argv[2];
        char* seq1_fn = argv[3];
        char* ns2_str = argv[4];
        char* sl2_fn = argv[5];
        char* seq2_fn = argv[6];
        char* outfile_path = argv[7];

        int ns1, ns2, sl1, sl2, max_sl1, max_sl2;
        ns1 = atoi(ns1_str);
        ns2 = atoi(ns2_str);
        int* seq_lens1 = (int *) malloc((ns1+2) * sizeof(int));
        int* seq_lens2 = (int *) malloc((ns2+2) * sizeof(int));
        long long* seq_pos2 = (long long *) malloc((ns2) * sizeof(long long));
        
        // this will be a matrix of the results that we would write to a binary file
        // ...but we're just going to do tab-delimited instead...
        //
        // long long* allResults = (long long *) malloc(ns2 * ns1 * sizeof(long long));

        
        FILE *ns1_f = fopen(sl1_fn,"r");
        fread((void*)(seq_lens1), 4, ns1, ns1_f);
        fclose(ns1_f);
        
        FILE *ns2_f = fopen(sl2_fn,"r");
        fread((void*)(seq_lens2), 4, ns2, ns2_f);
        fclose(ns2_f);
        
        // printf("ns1: %d\tns2: %d\n", ns1, ns2);
        // printf("seq_lens1: %d\tseq_lens2: %d\n", seq_lens1[0], seq_lens2[0]);
        
        // get the length of the longest sequence in each group, to malloc only once.
        max_sl1 = intArrayMax(seq_lens1,ns1);
        max_sl2 = intArrayMax(seq_lens2,ns2);
        long long seq2_totalsize=0;
        for (i=0;i<ns2;i++){
            seq_pos2[i] =seq2_totalsize;
            seq2_totalsize+=seq_lens2[i];
        }
        
        
        char* c1 = (char *)malloc(max_sl1+20);
        char** c2 = (char **)malloc(ns2 * sizeof(char *));
        
        FILE *seq1_f = fopen(seq1_fn,"r");
        FILE *seq2_f = fopen(seq2_fn,"r");
        for (i=0; i<ns2; i++) {
            c2[i] = (char *)malloc((seq_lens2[i]+1));
            fgets(c2[i], seq_lens2[i]+1, seq2_f);
        }
        fclose(seq2_f);
        
            
        long long* myFmat = (long long *)malloc((max_sl1+10) * (max_sl2+10) * sizeof(long long));
		int* dirmat = (int*)malloc((max_sl1+10) * (max_sl2+10) * sizeof(int));
        
        FILE *results_f = fopen(outfile_path,"w");
        
        
        for (i=0; i<ns1; i++) {
            fgets(c1, seq_lens1[i]+1, seq1_f);
            for (j=0; j<ns2; j++) {
                res = needleScore(seq_lens1[i], c1, seq_lens2[j], c2[j], myFmat, d, m , g, a, dirmat);
                fprintf(results_f,"%lld\t",res);
            }
            fprintf(results_f,"\n");
        }
        
        // clean up after ourselves.
        fclose(seq1_f);
        fclose(results_f);
        
        for (i=0; i<ns2; i++) {
            free(c2[i]);
        }
        
        free(myFmat);
		free(dirmat);
        free(c1);
        free(c2);
        free(seq_lens1);
        free(seq_lens2);
        free(seq_pos2);
    } else if (use_case==2){  //should probably try to make this avoid duplicate ccode, but not for today.
        char* s1 = (char *)malloc(5000);
        char* s2 = (char *)malloc(5000);
        char* s3 = (char *)malloc(5000);
        char* s4 = (char *)malloc(5000);
        size_t s1len, s2len, s3len, s4len, totlen, maxlen;
        long s1size, s2size, s3size, s4size;
        s1size = 5000; s2size = 5000; s3size = 5000; s4size = 5000;
        s1len = 0; s2len = 0; s3len = 0; s4len = 0;

        while (1)
        {
            fgets(&s1[s1len],5000,stdin);
            s1len = strlen(s1);
            if (s1len+1 < s1size) {
                break;
            } else {
                s1 = (char *)realloc(s1, s1size + 5000 - 1);
            }
            s1size += 4999;
        }
        while (1)
        {
            fgets(&s2[s2len],5000,stdin);
            s2len = strlen(s2);
            if (s2len+1 < s2size) {
                break;
            } else {
                s2 = (char *)realloc(s2, s2size + 5000 - 1);
            }
            s2size += 4999;
        }
        while (1)
        {
            fgets(&s3[s3len],5000,stdin);
            s3len = strlen(s3);
            if (s3len+1 < s3size) {
                break;
            } else {
                s3 = (char *)realloc(s3, s3size + 5000 - 1);
            }
            s2size += 4999;
        }
        while (1)
        {
            fgets(&s4[s4len],5000,stdin);
            s4len = strlen(s4);
            if (s4len+1 < s4size) {
                break;
            } else {
                s4 = (char *)realloc(s4, s4size + 5000 - 1);
            }
            s4size += 4999;
        }
        if (s1[s1len-1]=='\n')
        {
            s1[s1len-1]='\0';
            s1len = strlen(s1);
        }
        if (s2[s2len-1]=='\n')
        {
            s2[s2len-1]='\0';
            s2len = strlen(s2);
        }
        if (s3[s3len-1]=='\n')
        {
            s3[s3len-1]='\0';
            s3len = strlen(s3);
        }
        if (s4[s4len-1]=='\n')
        {
            s4[s4len-1]='\0';
            s4len = strlen(s4);
        }

        maxlen=(s2len > s3len ? (s1len > s2len ? s1len : s2len) : (s1len > s3len ? s1len : s3len));
        totlen = s4len + maxlen;
        char* aln1 = (char *)malloc(totlen + 1);
        char* aln2 = (char *)malloc(totlen + 1);
        char* aln3 = (char *)malloc(totlen + 1);
        char* aln4 = (char *)malloc(totlen + 1);
        char* aln5 = (char *)malloc(totlen + 1);
        char* aln6 = (char *)malloc(totlen + 1);

        // For Debugging:
        // printf("s1len: %d\ts2len: %d\n",s1len,s2len);
        // printf("%s\n",s1);
        // printf("%s\n",s2);

        long long* myFmat = (long long *)malloc((maxlen+2) * (s4len+2) * sizeof(long long));
        int* dirmat = (int *)malloc((maxlen+2) * (s4len+2) * sizeof(int));
        int n_gaps, n_gaps2, n_gaps3, aln_length, aln_length2, aln_length3, match_ct, diff_ct, gap_ct, maxres;
        match_ct = 0;  diff_ct = 0; gap_ct = 0;
        // long long d, m, g, a, res;
//        d = -1; m = 2; g = -4; a = 1;
//        d = 0; m = 1; g = 0; a = 0;
        long long res1, res2, res3;
        res1 = getNeedleAlignment(s1len, s1, s4len, s4, myFmat, d, m, g, a, aln1, aln4, dirmat, &aln_length, &n_gaps);
        zeroLongLongArray(myFmat, (maxlen+2) * (s4len+2));
        res2 = getNeedleAlignment(s2len, s2, s4len, s4, myFmat, d, m, g, a, aln2, aln5, dirmat, &aln_length2, &n_gaps2);
        zeroLongLongArray(myFmat, (maxlen+2) * (s4len+2));
        res3 = getNeedleAlignment(s3len, s3, s4len, s4, myFmat, d, m, g, a, aln3, aln6, dirmat, &aln_length3, &n_gaps3);

        maxres = (res1 > res2 ?  (res1 > res3 ? 0 : 2) : (res2 > res3 ? 1 : 2));

        if (res2 > res1) {
            if (res3 > res2) {
                strcpy(aln1, aln3);
                strcpy(aln4, aln6);
                aln_length = aln_length3;
                n_gaps = n_gaps3;
                res1 = res3;
            } else {
                strcpy(aln1,aln2);
                strcpy(aln4,aln5);
                aln_length = aln_length2;
                n_gaps = n_gaps2;
                res1 = res2;
            }
        }

        /*if (res1 > res2) {
            if (res1 > res3) {
                for (i=0; i<aln_length; i++) {
                    if (aln1[i]==45 || aln4[i]==45) {
                        gap_ct++;
                    } else if (aln1[i]==aln4[i]) {
                        match_ct++;
                    } else {
                        diff_ct++;
                    }
                }
            } else {
                for (i=0; i<aln_length; i++) {
                    if (aln3[i]==45 || aln6[i]==45) {
                        gap_ct++;
                    } else if (aln3[i]==aln6[i]) {
                        match_ct++;
                    } else {
                        diff_ct++;
                    }
                }
            }
        } else {
            if (res2 > res3) {
                for (i=0; i<aln_length; i++) {
                    if (aln2[i]==45 || aln5[i]==45) {
                        gap_ct++;
                    } else if (aln2[i]==aln5[i]) {
                        match_ct++;
                    } else {
                        diff_ct++;
                    }
                }
            } else {
                for (i=0; i<aln_length; i++) {
                    if (aln3[i]==45 || aln6[i]==45) {
                        gap_ct++;
                    } else if (aln3[i]==aln6[i]) {
                        match_ct++;
                    } else {
                        diff_ct++;
                    }
                }
            }
        }*/
        for (i=0; i<aln_length; i++) {
            if (aln1[i]==45 || aln4[i]==45) {
                gap_ct++;
            } else if (aln1[i]==aln4[i]) {
                match_ct++;
            } else {
                diff_ct++;
            }
        }
        // printf("%lld\n",res);
        // getNeedleAlignment(s1len, s1, s2len, s2, myFmat, d, m, g, a, aln1, aln2);
        // output (seq_1 len, seq_2 len, match ct, diff ct, alignment length, gap ct, score
        if (pretty_output==1) {
            printf("Seq 1:           %s\n", s1);
            printf("Seq 2:           %s\n", s2);
            printf("Seq 3:           %s\n", s3);
            printf("Seq 4:           %s\n", s4);
            printf("Seq 1 Length:    %d\n", s1len);
            printf("Seq 2 Length:    %d\n", s2len);
            printf("Seq 3 Length:    %d\n", s3len);
            printf("Seq 4 Length:    %d\n", s4len);
            printf("Best (0,1,2):    %d\n", maxres);
            printf("Aligned Length:  %d\n", aln_length);
            printf("# Matches:       %d\n", match_ct);
            printf("# Differences:   %d\n", diff_ct);
            printf("# Gaps:          %d\n", n_gaps);
            printf("Score:           %d\n", res1);
            printf("Aligned Seq 1:   %s\n", aln1);
            printf("Aligned Seq 2:   %s\n", aln4);
        } else {
            printf("%d\t%d\t%d\t%d\t", s1len, s4len, match_ct, diff_ct);
            printf("%d\t%d\t%lld\t", aln_length, n_gaps, res1);  // add: match ct, diff ct, seq1_len, seq2_len,
            printf(aln1);
            printf("\t");
            printf(aln4);
            printf("\t%d\n",maxres);
        }
        // print_f_matrix(myFmat,s1len, s2len);

        free(s1);
        free(s2);
        free(s3);
        free(s4);
        free(aln1);
        free(aln2);
        free(aln3);
        free(aln4);
        free(aln5);
        free(aln6);
        free(myFmat);
        free(dirmat);
    } else {
        // take two newline-terminated sequences from standard input, 2000 characters at a time
        // Right now there is a bug in this when the input is longer than 2k sequences_file1
        char* s1 = (char *)malloc(5000);
        char* s2 = (char *)malloc(5000);
        size_t s1len, s2len, totlen;
        long s1size, s2size;
        s1size = 5000; s2size = 5000;
        s1len = 0; s2len = 0;
        
        while (1)
        {
            fgets(&s1[s1len],5000,stdin);
            s1len = strlen(s1);
            if (s1len+1 < s1size) {
                break;
            } else {
                s1 = (char *)realloc(s1, s1size + 5000 - 1);
            }
            s1size += 4999;
        }
        while (1)
        {
            fgets(&s2[s2len],5000,stdin);
            s2len = strlen(s2);
            if (s2len+1 < s2size) {
                break;
            } else {
                s2 = (char *)realloc(s2, s2size + 5000 - 1);
            }
            s2size += 4999;
        }
        if (s1[s1len-1]=='\n')
        {
            s1[s1len-1]='\0';
            s1len = strlen(s1);
        }
        if (s2[s2len-1]=='\n')
        {
            s2[s2len-1]='\0';
            s2len = strlen(s2);
        }
		
		totlen = s1len + s2len;
		char* aln1 = (char *)malloc(totlen + 1);
		char* aln2 = (char *)malloc(totlen + 1);
        
        // For Debugging:
        // printf("s1len: %d\ts2len: %d\n",s1len,s2len);
        // printf("%s\n",s1);
        // printf("%s\n",s2);
        
        long long* myFmat = (long long *)malloc((s1len+2) * (s2len+2) * sizeof(long long));
		int* dirmat = (int *)malloc((s1len+2) * (s2len+2) * sizeof(int));
		int n_gaps, aln_length, match_ct, diff_ct, gap_ct;
		match_ct = 0;  diff_ct = 0; gap_ct = 0;
        // long long d, m, g, a, res;
//        d = -1; m = 2; g = -4; a = 1;
//        d = 0; m = 1; g = 0; a = 0;
        res = getNeedleAlignment(s1len, s1, s2len, s2, myFmat, d, m, g, a, aln1, aln2, dirmat, &aln_length, &n_gaps);

        for (i=0; i<aln_length; i++) {
            if (aln1[i]==45 || aln2[i]==45) {
                gap_ct++;
            } else if (aln1[i]==aln2[i]) {
                match_ct++;
            } else {
                diff_ct++;
            }
        }
        // printf("%lld\n",res);
		// getNeedleAlignment(s1len, s1, s2len, s2, myFmat, d, m, g, a, aln1, aln2);
        // output (seq_1 len, seq_2 len, match ct, diff ct, alignment length, gap ct, score
        if (pretty_output==1) {
            printf("Seq 1 Length:    %d\n", s1len);
            printf("Seq 2 Length:    %d\n", s2len);
            printf("Aligned Length:  %d\n", aln_length);
            printf("# Matches:       %d\n", match_ct);
            printf("# Differences:   %d\n", diff_ct);
            printf("# Gaps:          %d\n", n_gaps);
            printf("Score:           %d\n", res);
            printf("Aligned Seq 1:   %s\n", aln1);
            printf("Aligned Seq 2:   %s\n", aln2);
        } else {
            printf("%d\t%d\t%d\t%d\t", s1len, s2len, match_ct, diff_ct);
            printf("%d\t%d\t%lld\t", aln_length, n_gaps, res);  // add: match ct, diff ct, seq1_len, seq2_len,
            printf(aln1);
            printf("\t");
            printf(aln2);
            printf("\n");
        }
		// print_f_matrix(myFmat,s1len, s2len);
		
		free(myFmat);
		free(dirmat);
    }
	exit(0);
}
