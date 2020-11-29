#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"
#include "parson.h"
#include "parson.c"
#include <stdio.h>

void print_cell_info(LIFneuron *cells, int n);
void print_syn_info(ExpSyn *syns, int n);
void print_stim_info(Pinput *stims, int n);

extern LIFneuron *cells;
extern ExpSyn *syns;
extern Pinput *stims;

extern int NumofCells, NumofSyns, NumofPinputs;

int main(){

    char fname[] = "test.json";
    int n=0;

    readNetworkInfo(fname);

    printf("%d, %d, %d\n", NumofCells, NumofSyns, NumofPinputs);

    printf("Check Cell info\n");
    print_cell_info(cells, 2);

    printf("Check Synapse info\n");
    print_syn_info(syns, 2);
    print_syn_info(syns, 30);

    printf("Check Stim info\n");
    print_stim_info(stims, 2);


    printf("Done\n");


}

void print_cell_info(LIFneuron *cells, int n){
    printf("tau=%f, r=%f, e=%f, v=%f, vth=%f, vahp=%f, vmax=%f, t_refrac=%f\n",
           cells[n].tau, cells[n].r, cells[n].e, cells[n].v, cells[n].vth, cells[n].vahp, cells[n].vmax, cells[n].t_refrac);
}

void print_syn_info(ExpSyn *syns, int n){
    printf("tau1=%f, tau2=%f, A=%f, gmax=%f, e=%f, d=%f, id_pre=%d, id_post=%d\n",
           syns[n].tau1, syns[n].tau2, syns[n].A, syns[n].gmax, syns[n].e, syns[n].d, syns[n].id_pre, syns[n].id_post);
}

void print_stim_info(Pinput *stims, int n){
    printf("tstart=%f, tend=%f, p=%f\n", stims[n].tstart, stims[n].tend, stims[n].p);
}
