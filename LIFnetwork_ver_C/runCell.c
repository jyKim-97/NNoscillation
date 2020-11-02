#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"


int main(){
    // testing
    LIFneuron cell;
    double *vs;
    FILE *fid = fopen("./v_test.csv", "w");
    int nitr;

    _tmax = 500;
    _dt = 0.01;
    nitr = (int) (_tmax/_dt) + 1;
    vs = (double*) malloc(sizeof(double) * nitr);

    // init vcells
    cell.v = -65;
    cell.i = 0;
    cell.t0 = -100;
    cell.is_refrac = 0;
    cell.flag = 0;
    cell.tau = 20;
    cell.r = 100;
    cell.e = -65;
    cell.vth = -40;
    cell.vahp = -80;
    cell.vmax = 30;
    cell.t_refrac = 5;

    runOnecell_w_current(&cell, 100, 400, 1, &vs);

    // write
    for (int i=0; i<nitr; i++) {
        fprintf(fid, "%lf,", vs[i]);
    }

    free(vs);
    fclose(fid);
    printf("Done\n");
}


















