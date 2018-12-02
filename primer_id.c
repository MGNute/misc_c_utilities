#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MP 50
#define maxlines 1000000
#define nuc_cutoff 0.05
#define takecomp 0

char linebuf[10000];
int nlines;
int slen, slentemp;

long maxprimer[5*MP];
float maxprimer_f[5*MP];
// A=0, C=1, G=2, T=3, other=4
// ACTGactg = 65, 67, 84, 71, 97, 99, 116, 103


int main(int argc, char *argv[])
{
   /*  // while(!e
    fgets(linebuf, 200, stdin);
    slen = strlen(linebuf);
    printf("\tString Length: %d\n",slen);
    // printf("\tthe values of the last 8 characters are: %d, %d, %d, %d, %d, %d, %d,    %d\n",linebuf[0],linebuf[1],linebuf[2],linebuf[3],linebuf[4],linebuf[5],linebuf[6],linebuf[7]);
    // printf("\tlast three: %d, %d, %d, %d\n",linebuf[slen-2],linebuf[slen-1],linebuf[slen], linebuf[slen+1]);
    
    fgets(&linebuf[slen], 200, stdin);
    slen += strlen(linebuf);
    slentemp = strlen(linebuf);
    printf("\tString Length: %d\nslen: %d\n",slentemp,slen); */
    
    int i, j, r, tot, a;
    char c, nuc;
    float q;
    for (i=0; i<5*MP; i++)
    {
        maxprimer[i]=0;
    }
    
    
    while(feof(stdin)==0 && nlines < maxlines) {
        fgets(linebuf, 10000, stdin);
        slen = strlen(linebuf);
        if (slen==0) break;
        for (i=0; i<MP; i++) {
            if (i >= slen) break;
            
            c = linebuf[i];
            switch(c) {
                case 'A':
                    r = 0; break;
                case 'a':
                    r = 0; break;
                case 'C':
                    r = 1; break;
                case 'c':
                    r = 1; break;
                case 'G':
                    r = 2; break;
                case 'g':
                    r = 2; break;
                case 'T':
                    r = 3; break;
                case 't':
                    r = 3; break;
                default:
                    r=4;
            }
            
            maxprimer[i*5+r]++;
        }
        nlines++;
    }
    printf("# lines: %d\n",nlines);
    
    for (i=0; i<MP; i++)
    {
        tot = maxprimer[i*5] + maxprimer[i*5+1] + maxprimer[i*5+2]+maxprimer[i*5+3]+maxprimer[i*5+4];
        for (j=0; j<5; j++)
        {
            a=maxprimer[i*5+j]; 
            q = (float)a / (float)tot;
            maxprimer_f[i*5+j] = q;
        }
    }
    
    printf("first 20 cols:\n");
    printf("A: ");
    for (i=0;i<20;i++) printf("%.1f ",maxprimer_f[i*5]);
    printf("\nC: ");
    for (i=0;i<20;i++) printf("%.1f ",maxprimer_f[i*5+1]);
    printf("\nT: ");
    for (i=0;i<20;i++) printf("%.1f ",maxprimer_f[i*5+2]);
    printf("\nG: ");
    for (i=0;i<20;i++) printf("%.1f ",maxprimer_f[i*5+3]);
    printf("\n?: ");
    for (i=0;i<20;i++) printf("%.1f ",maxprimer_f[i*5+4]);
    printf("\n\nConsensus:\n");
    
    for (i=0; i<20; i++) {
        c=0;
        if (maxprimer_f[i*5] > nuc_cutoff) c = c | 1;
        if (maxprimer_f[i*5+1] > nuc_cutoff) c = c | 2;
        if (maxprimer_f[i*5+2] > nuc_cutoff) c = c | 4;
        if (maxprimer_f[i*5+3] > nuc_cutoff) c = c | 8;
        switch((int)c) {
            case 1:
                nuc = (takecomp ? 'A' : 'T' ); break;
            case 2:
                nuc = (takecomp ? 'C' : 'G' ); break;
            case 3:
                nuc = (takecomp ? 'M' : 'K' ); break;
            case 4:
                nuc = (takecomp ? 'T' : 'A' ); break;
            case 5:
                nuc = (takecomp ? 'W' : 'W' ); break;
            case 6:
                nuc = (takecomp ? 'Y' : 'R' ); break;
            case 7:
                nuc = (takecomp ? 'H' : 'D' ); break;
            case 8:
                nuc = (takecomp ? 'G' : 'C' ); break;
            case 9:
                nuc = (takecomp ? 'R' : 'Y' ); break;
            case 10:
                nuc = (takecomp ? 'S' : 'S' ); break;
            case 11:
                nuc = (takecomp ? 'V' : 'B' ); break;
            case 12:
                nuc = (takecomp ? 'K' : 'M' ); break;
            case 13:
                nuc = (takecomp ? 'D' : 'H' ); break;
            case 14:
                nuc = (takecomp ? 'B' : 'V' ); break;
            case 15:
                nuc = (takecomp ? 'N' : 'N' ); break;
        }
        printf("%c",nuc);
        // printf("%d-%c-",c,nuc);
        // if (c==2) printf("T");
        // printf("  ");
    }
    
    
    
    
    // fgets(linebuf, 200, stdin);
    // slen = strlen(linebuf);
    // printf("\tString Length: %d\n",slen);
    
    // if (eof(stdin)) {
        // printf("eof is true\n");
    // }
}