#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define _GNU_SOURCE

#ifndef _ARTNET_H_
#define _ARTNET_H_

// global variable
extern double _dt;
extern double _tmax;
extern double _s;

// structs;
typedef struct {
    // dynamic variables
    double v, i;
    // detect firing
    double t0;
    // flag 1: go to vmax, flag 2: return to vahp
    int is_refrac, flag;
    // parameters for each neuron
    double tau, r, e, vth, vahp, vmax, t_refrac;
} LIFneuron;


typedef struct {
    // presynaptic parameter 
    double *ref_t0;
    // postsynaptic parameter
    double *ref_v, *ref_i;
    int id_pre, id_post;
    int *ref_is_refrac;
    // dynamic variable
    double g;
    // parameters for each synapse
    double tau1, tau2, A, gmax, e, d;
} ExpSyn;


typedef struct {
    /*
    spike times: dynamic array, store spike timees
    ref_t0: spiking time indicator (event occuring time)
    is_refrac: always 0
    len
    n: counter
    id_target: target cell number
    */
    double t0, *ref_spk; // allocate dynamic array
    int len, n, id_target;
} Stim;


// functions
// void readNetwork(char *fdir, char *prefix, LIFneuron **cells, ExpSyn **syns, Stim **stims, int *n_cells, int *n_syns, int *n_stims);
// void runAll_w_save(LIFneuron **cells, ExpSyn **syns, Stim **stims, int n_cells, int n_syns, int n_stims, FILE *fid);

// update each parts
void updateLIF(LIFneuron *cell, double t);
// void updateSyn(ExpSyn *syn, double t);
void updateSyn(ExpSyn *syn, double t);
void updateStim(Stim *stim, double *spk_times, double t);

// solve 
double solveLIF(LIFneuron cell);
double fLIF(LIFneuron cell, double v);

// read input params
void readCellParams(FILE *fid, LIFneuron *cell);
void readSynParams(FILE *fid, ExpSyn *syn);
void readStimParams(FILE *fid, Stim *ext);

// test function
void runOnecell_w_current(LIFneuron *cell, double t0, double t1, double amp, double **vs);

// subfunctions
void saveData(double *vs, int n, FILE *fid);

#endif

/* Parameter information
Neuron: LIF neuron

Synspae: ExpSyn

Stim: 


*/