#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"
#include "parson.h"
#include "parson.c"
#include <stdio.h>
#include <time.h>


extern LIFneuron *cells;
extern ExpSyn *syns;
extern Pinput *stims;
extern double t;
extern int NumofCells, NumofSyns, NumofPinputs;


void print_v(FILE *fid);


int main(){

    double tmax=10000;
    FILE *fid = fopen("test_out.csv", "w");
    clock_t tic, toc;

    init_genrand64(1000);

    char fname[] = "test.json";
    int n=0;
    LIFneuron *cells;
    ExpSyn *syns;
    Pinput *stims;
    int NumofCells, NumofSyns, NumofPinputs;

    readNetworkInfo(fname);
    printf("running ... ");
    tic = clock();

    while (t <= tmax){
        print_v(fid);
        updateAll();
    }
    print_v(fid);

    fclose(fid);
    freeObjs();

    toc = clock();
    printf("Done, execution time = %.5fs\n", (double) ((toc-tic) / CLOCKS_PER_SEC));

    return 0;

}

void print_v(FILE *fid){
    int n;

    fprintf(fid, "%f,", t);
    for (n=0; n<NumofCells; n++){
        fprintf(fid, "%f,", cells[n].v);
    }
    fprintf(fid, "\n");
}