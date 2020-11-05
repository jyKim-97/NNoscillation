#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"

void ProgressBar(char *label, int step, int maxstep, float etime);

int main(int argc, char **argv){
    char *fdir, *prefix, *fsave, fname[1024];
    FILE *fidv, *fidi;
    int n_cells, n_syns, n_stims, nitr, is_current=0;
    LIFneuron *cells;
    ExpSyn *syns;
    Stim *stims;
    IClamp *ics;
    double **spike_times, *vcells, t=0;
    double t0, t1, amp;
    char *label="Progress simulation...";
    clock_t tic, toc;


    for (int i=0; i<argc; i++){
        if (strstr(argv[i], "--")){
            if (!strcmp("--fdir", argv[i])){
                fdir = argv[i+1];
            } else if (!strcmp("--prefix", argv[i])){
                prefix = argv[i+1];
            } else if (!strcmp("--fsave", argv[i])){
                fsave = argv[i+1];
            } else if (!strcmp("--tmax", argv[i])){
                _tmax = atof(argv[i+1]);
            } else if (!strcmp("--dt", argv[i])){
                _dt = atof(argv[i+1]);
            } else if (!strcmp("--s", argv[i])){
                _s = atof(argv[i+1]);
            } else if (!strcmp("-c", argv[i])){
                is_current = 1;
            } else if (!strcmp("--t0", argv[i])){
                t0 = atof(argv[i+1]);
            } else if (!strcmp("--t1", argv[i])){
                t1 = atof(argv[i+1]);
            } else if (!strcmp("--amp", argv[i])){
                amp = atof(argv[i+1]);
            } else {
                printf("Usage: %s [--fdir, --prefix, --savename] filename", argv[0]);
            }
        }
    }

    printf("Target file:\n");
    printf("%s%s_info.csv\n", fdir, prefix);
    printf("%s%s_cell.csv\n", fdir, prefix);
    printf("%s%s_syn.csv\n", fdir, prefix);
    if (is_current == 1){
        printf("%s%s_IClamp.csv\n", fdir, prefix);
    }
    printf("--> result will save to %s\n", fsave);


    // read params
    readParams(fdir, prefix, &cells, &syns, &stims, &spike_times, &n_cells, &n_syns, &n_stims, 0);

    // printf("t0=%lf, t1=%lf, amp=%lf\n", t0, t1, amp);
    // run network

    sprintf(fname, "%s_v.csv", fsave);
    fidv = fopen(fname, "w");
    // printf("%s\n", fname);

    sprintf(fname, "%s_i.csv", fsave);
    fidi = fopen(fname, "w");
    // printf("%s\n", fname);
    
    // fidv = fopen("./out_v.csv", "w");
    // fidi = fopen("./out_i.csv", "w");

    // save init v, i
    fprintf(fidv, "%f", t);
    fprintf(fidi, "%f", t);
    for (int i=0; i<n_cells; i++){
        fprintf(fidv, ",%f", cells[i].v);
        fprintf(fidi, ",%f", cells[i].i);
    }
    fprintf(fidv, "\n");
    fprintf(fidi, "\n");


    nitr = (int) (_tmax / _dt);
    tic = clock();
    for (int n=0; n<nitr; n++){


        t += _dt;
        // update synapse
        for (int i=0; i<n_syns; i++) {
            updateSyn(syns+i, t);
        }
        // current
        if ((t >= t0) && (t <= t1)){
            for (int i=0; i<225; i++){ //// temporally
                if (cells[i].is_refrac == 0){
                    cells[i].i += amp;
                }
            }
        }
        // update cell & save data
        fprintf(fidv, "%f", t);
        fprintf(fidi, "%f", t);
        for (int i=0; i<n_cells; i++) {
            fprintf(fidi, ",%f", cells[i].i);
            updateLIF(cells+i, t);
            // save data
            fprintf(fidv, ",%f", cells[i].v);
            
        }
        fprintf(fidv, "\n");
        fprintf(fidi, "\n");

        if (n - (n/100)*100 == 0){
            ProgressBar(label, n, nitr, (float) ((clock() - tic)/CLOCKS_PER_SEC));
        }
    }
    
    printf("\nExecution Done, %.3fs\n", ((float) (clock() - tic)/(float) CLOCKS_PER_SEC));
    
    free(syns);
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