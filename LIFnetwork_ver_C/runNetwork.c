#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"

void ProgressBar(char *label, int step, int maxstep, float etime);

extern double _tmax, _dt, _s;
double _t=0;

int main(int argc, char **argv){
    char *fdir, *prefix, *fsave, fname[1024];
    FILE *fidv, *fidi;
    int n_cells, n_syns, n_stims, nitr, is_current=0;
    LIFneuron *cells;
    ExpSyn *syns;
    Stim *stims;
    IClamp *ics;
    double **spike_times, *vcells;
    char *label="Prgress simulation...";
    unsigned long long seed = 0;
    clock_t tic, toc;


    for (int i=0; i<argc; i++){
        if (strstr(argv[i], "--")){
            if (!strcmp("--fdir", argv[i])){
                fdir = argv[i+1];
            } else if (!strcmp("--prefix", argv[i])){
                prefix = argv[i+1];
            } else if (!strcmp("--fsave", argv[i])){
                // sprintf(fsave, "%s", argv[i+1]);
                fsave = argv[i+1];
            } else if (!strcmp("--tmax", argv[i])){
                _tmax = atof(argv[i+1]);
            } else if (!strcmp("--dt", argv[i])){
                _dt = atof(argv[i+1]);
            } else if (!strcmp("--s", argv[i])){
                _s = atof(argv[i+1]);
            } else if (!strcmp("-c", argv[i])){
                is_current = 1;
            } else if (!strcmp("--seed", argv[i])){
                seed = atoi(argv[i]);
            } else {
                printf("Usage: %s [--fdir, --prefix, --fsave] filename", argv[0]);
            }
        }
    }

    if (seed == 0){
        seed = time(NULL);
    }

    // printf("Target file:\n");
    // printf("%s%s_info.csv\n", fdir, prefix);
    // printf("%s%s_cell.csv\n", fdir, prefix);
    // printf("%s%s_syn.csv\n", fdir, prefix);
    // printf("%s%s_t_spike.csv\n", fdir, prefix);
    // if (is_current == 1){
    //     printf("%s%s_IClamp.csv\n", fdir, prefix);
    // }
    // printf("--> result will save to %s\n", fsave);

    init_genrand64(seed);

    // read params
    readParams(fdir, prefix, &cells, &syns, &stims, &spike_times, &n_cells, &n_syns, &n_stims, 1);

    // run Network
    sprintf(fname, "%s_v.csv", fsave);
    fidv = fopen(fname, "w");
    // printf("%s\n", fname);

    sprintf(fname, "%s_i.csv", fsave);
    fidi = fopen(fname, "w");
    
    // save init v, i
    fprintf(fidv, "%f", _t);
    fprintf(fidi, "%f", _t);
    for (int i=0; i<n_cells; i++){
        fprintf(fidv, ",%f", cells[i].v);
        fprintf(fidi, ",%f", cells[i].i);
    }
    fprintf(fidv, "\n");
    fprintf(fidi, "\n");

    nitr = (int) (_tmax / _dt);
    tic = clock();
    for (int n=0; n<nitr; n++){

        _t += _dt;
        // update stim counter
        for (int i=0; i<n_stims; i++){
            updateStim(stims+i, spike_times[i], _t);
        }
        // update synapse
        for (int i=0; i<n_syns+n_stims; i++) {
            updateSyn(syns+i, _t);

        }
        // update cell & save data
        fprintf(fidv, "%f", _t);
        fprintf(fidi, "%f", _t);
        for (int i=0; i<n_cells; i++) {
            fprintf(fidi, ",%f", cells[i].i);            
            updateLIF(cells+i, _t);
            // save data
            fprintf(fidv, ",%f", cells[i].v);
        }
        fprintf(fidv, "\n");
        fprintf(fidi, "\n");

        // if (n - (n/100)*100 == 0){
        //     ProgressBar(label, n, nitr, (float) ((clock() - tic)/CLOCKS_PER_SEC));
        // }

    }

    // printf("\nExecution Done, %.3fs\n", (float) ((clock() - tic)/CLOCKS_PER_SEC));
    

    free(syns);
    for (int i=0; i<n_stims; i++){
        if (stims[i].len != 0){
            free(spike_times[i]);
        }
    }
    free(spike_times);
    free(stims);
    free(cells);
    fclose(fidv);
    fclose(fidi);

    return 0;

}

void ProgressBar(char *label, int step, int maxstep, float etime){

    const int maxw=100, len=strlen(label);
    int w=maxw-len, n=(step*w)/maxstep;

    // printf("\r%*c", 100, '*');
    printf("\r%s [", label);
    for (int i = 0;i < n; i++)  printf("%c", '=');
    printf("%*c %.2fs left", w-n, ']', etime/step * (maxstep-step));
}