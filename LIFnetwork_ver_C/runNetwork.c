#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ArtNetSimulTools.h"
#include "ArtNetSimulTools.c"
// #include <unistd.h>

void ProgressBar(char *label, int step, int maxstep, float etime);

int main(int argc, char **argv){
    char *fdir, *prefix, *savename, fname[128], buffer[128];//="./";
    // char ;//="test";
    // char *savename="./test_v.csv";
    // char fname[128], buffer[128];
    FILE *fid_obj, *fid_info;
    int n_cells, n_syns, n_stims, tmp_int, id_pre, id_post, nitr;
    LIFneuron *cells;
    ExpSyn *syns;
    Stim *stims;
    double **spike_times, *vcells, t=0;
    char *label="Prgress simulation...";
    clock_t tic, toc;


    for (int i=0; i<argc; i++){
        if (strstr(argv[i], "--")){
            if (!strcmp("--fdir", argv[i])){
                fdir = argv[i+1];
            } else if (!strcmp("--prefix", argv[i])){
                prefix = argv[i+1];

            } else if (!strcmp("--savename", argv[i])){
                // sprintf(savename, "%s", argv[i+1]);
                savename = argv[i+1];

            } else if (!strcmp("--tmax", argv[i])){
                _tmax = atof(argv[i+1]);

            } else if (!strcmp("--dt", argv[i])){
                _dt = atof(argv[i+1]);
            } else if (!strcmp("--s", argv[i])){
                _s = atof(argv[i+1]);
            } else {
                printf("Usage: %s [--fdir, --prefix, --savename] filename", argv[0]);
            }
        }
    }
    printf("Target file:\n");
    printf("%s%s_info.csv\n", fdir, prefix);
    printf("%s%s_cell.csv\n", fdir, prefix);
    printf("%s%s_syn.csv\n", fdir, prefix);
    printf("%s%s_t_spike.csv\n", fdir, prefix);
    printf("--> result will save to %s\n", savename);

    // open file
    sprintf(fname, "%s%s_info.csv", fdir, prefix);
    fid_info = fopen(fname, "r");

    fscanf(fid_info, "%s\n", buffer); // time
    fscanf(fid_info, "%s %d\n", buffer, &n_cells); // n_cells
    fscanf(fid_info, "%s %d", buffer, &n_syns); // n_syns
    fscanf(fid_info, "%s %d", buffer, &n_stims); // n_stims    
    
    fgets(buffer, 128, fid_info);
    fgets(buffer, 128, fid_info);

    cells = (LIFneuron*) malloc(sizeof(LIFneuron) * n_cells);
    stims = (Stim*) malloc(sizeof(Stim) * n_stims);
    syns = (ExpSyn*) malloc(sizeof(ExpSyn) * (n_syns+n_stims));
    spike_times = (double**) malloc(sizeof(double*) * n_stims);

    // cell obj
    sprintf(fname, "%s%s_cell.csv", fdir, prefix);

    fid_obj = fopen(fname, "r");
    
    fscanf(fid_obj, "%s\n", buffer);
    for (int i=0; i<n_cells; i++){
        readCellParams(fid_obj, cells+i);
    }
    fclose(fid_obj);

    // read stim times
    sprintf(fname, "%s%s_t_spike.csv", fdir, prefix);
    fid_obj = fopen(fname, "r");
    fscanf(fid_obj, "%s\n", buffer);

    for (int i=0; i<n_stims; i++){
        // readStimParams(fid_obj, stims+i);
        fscanf(fid_obj, "%d,%d,", &tmp_int, &(stims[i].len));

        spike_times[i] = (double*) malloc(sizeof(double) * stims[i].len);

        for (int j=0; j<stims[i].len-1; j++){
            fscanf(fid_obj, "%lf,", &(spike_times[i][j]));
        }
        fscanf(fid_obj, "%lf\n", &(spike_times[i][stims[i].len-1]));

        stims[i].ref_spk = spike_times[i];
        stims[i].n = 0;

    }

    // read syns
    sprintf(fname, "%s%s_syn.csv", fdir, prefix);
    fid_obj = fopen(fname, "r");
    fscanf(fid_obj, "%s\n", buffer);
    for (int i=0; i<n_syns+n_stims; i++){

        readSynParams(fid_obj, syns+i);
        if (i < n_syns) {
            fscanf(fid_info, "%d,%d,%d\n", &tmp_int, &id_pre, &id_post);
            syns[i].ref_t0 = &(cells[id_pre].t0);
        } else {
            fscanf(fid_info, "%d,%d\n", &tmp_int, &id_post);
            syns[i].ref_t0 = &(stims[i-n_syns].t0);
        }

        syns[i].ref_v = &(cells[id_post].v);
        syns[i].ref_i = &(cells[id_post].i);
        syns[i].ref_is_refrac = &(cells[id_post].is_refrac);

    }
    fclose(fid_obj);
    fclose(fid_info);

    // run network
    fid_obj = fopen(savename, "w");
    
    // save init v
    fprintf(fid_obj, "%f", t);
    for (int i=0; i<n_cells; i++){
        fprintf(fid_obj, ",%f", cells[i].v);
    }; fprintf(fid_obj, "\n");

    // for (int i=0; i<n_stims+n_syns; i++){
    //     printf("%lf, %lf, %lf\n", syns[i].gmax, syns[i].e, syns[i].d);
    // }
    // for (int i=0; i<n_cells; i++){
    //     printf("%lf, %lf\n", cells[i].tau, cells[i].r);
    // }

    nitr = (int) (_tmax / _dt);
    tic = clock();
    for (int i=0; i<nitr; i++){
        // printf("\naddress = %p - %p\n", syns[2].ref_t0, stims[2-n_syns].ref_t0);

        t += _dt;
        // update stim counter
        for (int i=0; i<n_stims; i++){
            updateStim(stims+i, spike_times[i], t);
        }
        // update synapse
        for (int i=0; i<n_syns+n_stims; i++) {
            updateSyn(syns+i, t);

        }
        // update cell & save data
        fprintf(fid_obj, "%f", t);
        for (int i=0; i<n_cells; i++) {
            updateLIF(cells+i, t);
            // save data
            fprintf(fid_obj, ",%f", cells[i].v);
        }
        fprintf(fid_obj, "\n");

        // if (i - (i/100)*100 == 0){
            // usleep(1);
            ProgressBar(label, i, nitr, (float) ((clock() - tic)/CLOCKS_PER_SEC));
        // }

    }
    
    printf("\nExecution Done, %.3fs\n", (float) ((clock() - tic)/CLOCKS_PER_SEC));
    
    free(syns);
    for (int i=0; i<n_stims; i++){
        // stims[i].ref_t0 = NULL;

        free(spike_times[i]);
    }
    free(spike_times);
    free(stims);
    free(cells);
    fclose(fid_obj);

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