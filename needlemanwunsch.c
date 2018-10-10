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
            mat = F[(i-1)*(len2+1)+(j-1)] + (seq1[i-1]==42 || seq2[j-1]==42 ? a : (seq1[i-1]==seq2[j-1] ? m : d)); // (dir = 1)
            del = F[(i-1)*(len2+1)+(j)] + g; // from above (dir = 2)
            ins = F[(i)*(len2+1)+(j-1)] + g; // from the left (dir = 0)
            tmp = (mat > del ? mat : del);
			dirtmp = (mat > del ? 1 : 2);
            F[i*(len2+1)+(j)] = (tmp > ins ? tmp : ins);
			dir = (tmp > ins ? dirtmp : 0);
			dirmat[i*(len2+1)+(j)] = dir;
        }
    }

    return F[(len2+1)*(len1+1)-1];
}

long long getNeedleAlignment(int len1, char* seq1, int len2, char* seq2, long long* F, long long d, long long m, long long g, 
						long long a, char* aln1, char* aln2, int* dirmat)
{
	int i, j, p, r, c;
    long long mat, del, ins, tmp, best_score;
	int dir, dirtmp;
	
	best_score = needleScore(len1, seq1, len2, seq2, F, d, m, g, a, dirmat);
	
	p = 0;
	r = len1; c = len2;
	// get the sequences
	while (r > 0 && c > 0)
	{
		dir = dirmat[r*(len2+1)+c];
		if (dir == 0)
		{
			aln1[p] = 45;
			aln2[p] = seq2[c-1];
			c--;
		} else if (dir == 2) {
			aln1[p] = seq1[r-1];
			aln2[p] = 45;
			r--;
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
	aln1[p] = 0;
	aln2[p] = 0;
	return best_score;
}

/* void getNeedleAlignment(int len1, char* seq1, int len2, char* seq2, long long* F, long long d, long long m, long long g, char* aln1, char* aln2, int* dirmat)
{
    // printf("The function to collect and return the alignment has not been implemented yet.\n");
	int r, c;
	long long d1, ma, d2;
	int dir; // just 0, 1 or 2. 0 = go left, 1 = go up-left, 2 = go up.
	int p;
	r=len1; c=len2; p=0;
	// printf("len1=%d\tlen2=%d\n",len1,len2);
	
	while(r >0 || c >0)
	{
		if (r>0) { // left value
			d1 = F[(r)*(len2+1)+c-1];
		} else {
			d1 = -99999999;
		}
		if (r>0 && c > 0) { // match value (top left)
			ma = F[(r-1)*(len2+1)+(c-1)];
		} else {
			ma = -99999999;
		}
		if (c > 0) { // top value
			d2 = F[(r-1)*(len2+1)+(c)];
		} else {
			d2 = -99999999;
		}
		
		// determin which direction to go
		if (d1 > ma) {
			if (d1 > d2) {
				dir = 0;
			} else {
				dir = 2;
			}
		} else {
			if (ma > d2) {
				dir = 1;
			} else {
				dir = 2;
			}
		}
		
		// advance and add to aligned sequences
		
		if (dir == 0) { // left
			// printf("p=%d\ta1= %c\ta2= %c, (r,c)=(%d,%d)\n",p,45,seq2[len2-c],r,c);
			aln1[p] = 45; //hyphen
			aln2[p] = seq2[c-1]; 
			c--;
		} else if (dir == 1) { //diagonal
			// printf("p=%d\ta1= %c\ta2= %c, (r,c)=(%d,%d)\n",p,seq1[len1-r],seq2[len2-c],r,c);
			aln1[p] = seq1[r-1];
			aln2[p] = seq2[c-1];
			r--;
			c--;
		} else { // up
			// printf("p=%d\ta1= %c\ta2= %c, (r,c)=(%d,%d)\n",p,seq1[len1-r],45,r,c);
			aln1[p] = seq1[r-1];
			aln2[p] = 45;
			r--;
		}
		p++;
	}
	inplace_reverse(aln1);
	inplace_reverse(aln2);
	aln1[p]=0;
	aln2[p]=0;
} */

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
     *
     *
     */

        
    // COST SCORES (change before compiling if desired): 
    //      (SHOULD PROBABLY MAKE THESE A COMMAND LINE ARGUMENT OR FROM A SETTINGS FILE OR SOMETHING)
    long long d, m, g, a, res;
    // d = mismatch penalty (*D*ifferent), m = Match score, g = Gap penalty.
    d = -1; m = 1; g = -1; a = 1;
        
    if (argc>2) {
        char* ns1_str = argv[1];
        char* sl1_fn = argv[2];
        char* seq1_fn = argv[3];
        char* ns2_str = argv[4];
        char* sl2_fn = argv[5];
        char* seq2_fn = argv[6];
        char* outfile_path = argv[7];
        int i, j;

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
    } else {
        // take two newline-terminated sequences from standard input, 2000 characters at a time
        // Right now there is a bug in this when the input is longer than 2k sequences_file1
        char* s1 = (char *)malloc(2000);
        char* s2 = (char *)malloc(2000);
        size_t s1len, s2len, totlen;
        long s1size, s2size;
        s1size = 2000; s2size = 2000;
        s1len = 0; s2len = 0;
        
        while (1)
        {
            fgets(&s1[s1len],2000,stdin);
            s1len = strlen(s1);
            if (s1len+1 < s1size) {
                break;
            } else {
                s1 = (char *)realloc(s1, s1size + 2000 - 1);
            }
            s1size += 1999;
        }
        while (1)
        {
            fgets(&s2[s2len],2000,stdin);
            s2len = strlen(s2);
            if (s2len+1 < s2size) {
                break;
            } else {
                s2 = (char *)realloc(s2, s2size + 2000 - 1);
            }
            s2size += 1999;
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
        // long long d, m, g, a, res;
        d = -1; m = 2; g = -4; a = 1;
        res = getNeedleAlignment(s1len, s1, s2len, s2, myFmat, d, m, g, a, aln1, aln2, dirmat);
		
        // printf("%lld\n",res);
		// getNeedleAlignment(s1len, s1, s2len, s2, myFmat, d, m, g, a, aln1, aln2);
		printf(aln1);
		printf("\t");
		printf(aln2);
		printf("\n");
		// print_f_matrix(myFmat,s1len, s2len);
		
		free(myFmat);
		free(dirmat);
    }
	exit(0);
}
