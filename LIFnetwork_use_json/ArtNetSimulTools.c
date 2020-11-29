
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parson.h"
#include "mt64.h"
#include "mt64.c"
#include "ArtNetSimulTools.h"

double _tmax = 1000, _dt = 0.01, _s = 0.1, t=0;

LIFneuron *cells;
ExpSyn *syns;
Pinput *stims;
int NumofCells, NumofSyns, NumofPinputs;


/* set neuron, synapse, stim object */
struct lif_neuron {
    // dynamic variables
    double v, i;
    // detect firing
    double t0;
    // flag 1: go to vmax, flag 2: return to vahp
    int is_refrac, flag;
    // parameters for each neuron
    double tau, r, e, vth, vahp, vmax, t_refrac;
};

struct exp_syn {
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
};

struct poisson_input {

    double tstart, tend, p, t;

};

struct i_clamp {
    double t_start, t_end, amp, *ref_i;
};


/* Parse Parameters */
// TODO : lifneuron type만 받을 수 있게 하는 것을 어떻게 처리할 것인가?, synapse / cell type 변경 방법 생각 필요
// void readNetworkInfo(char JSON_fname[], LIFneuron **cells, ExpSyn **syns, Pinput **stims, int *NumofCells, int *NumofSyns, int *NumofPinputs){
void readNetworkInfo(char JSON_fname[]){

    int n, seed, id_pre, id_post, NumOfPinputSynapses;
    JSON_Value *rootValue;
    JSON_Object *rootObject, *childObject;
    JSON_Array *arr;

    /* init JSON file */
    rootValue = json_parse_file(JSON_fname);
    rootObject = json_value_get_object(rootValue);

    /* read network information */
    NumofCells    = (int) json_object_get_number(rootObject, "NumofCells");
    NumofSyns     = (int) json_object_get_number(rootObject, "NumofSyns");
    NumofPinputs  = (int) json_object_get_number(rootObject, "NumofPinputs");
    NumOfPinputSynapses = (int) json_object_get_number(rootObject, "NumOfPinputSynapses");

    /* read seed number & initializing */
    // seed = json_object_get_number(rootObject, "seed");
    // init_genrand64(seed);

    /* read Cell information */
    arr = json_object_get_array(rootObject, "Cells");
    cells = (LIFneuron*) malloc(sizeof(LIFneuron) * NumofCells);
    for (n=0; n<NumofCells; n++){
        childObject = json_array_get_object(arr, n);
        cells[n] = readCellParams(childObject);
    }

    /* read Synapse information */
    arr = json_object_get_array(rootObject, "Synapses");
    syns = (ExpSyn*) malloc(sizeof(ExpSyn) * (NumofSyns + NumOfPinputSynapses));
    for (n=0; n<NumofSyns; n++){
        childObject = json_array_get_object(arr, n);
        syns[n] = readSynParams(childObject);

        // connect
        id_pre = syns[n].id_pre;
        syns[n].ref_t0 = &(cells[id_pre].t0);

    }

    /* read Stim information */
    arr = json_object_get_array(rootObject, "Pinputs");
    stims = (Pinput*) malloc(sizeof(Pinput) * NumofPinputs);
    for (n=0; n<NumofPinputs; n++){
        childObject = json_array_get_object(arr, n);
        stims[n] = readPinputParams(childObject);

    }

    /* read Poisson input synapse information */
    arr = json_object_get_array(rootObject, "PinputSynapses");
    for (n=0; n<NumOfPinputSynapses; n++){
        childObject = json_array_get_object(arr, n);
        syns[n+NumofSyns] = readSynParams(childObject);

        // connect
        id_pre = syns[n].id_pre;
        syns[n+NumofSyns].ref_t0 = &(stims[id_pre].t);

    }

    NumofSyns += NumOfPinputSynapses;
    /* connect synapse */ 
    for (n=0; n<NumofSyns; n++){

        id_post = syns[n].id_post;

        syns[n].ref_v = &(cells[id_post].v);
        syns[n].ref_i = &(cells[id_post].i);
        syns[n].ref_is_refrac = &(cells[id_post].is_refrac);
        
    }
    
}

LIFneuron readCellParams(JSON_Object *childObject){
    
    LIFneuron cell;

    cell.tau      = json_object_get_number(childObject, "tau");
    cell.r        = json_object_get_number(childObject, "r");
    cell.e        = json_object_get_number(childObject, "e");
    cell.v        = json_object_get_number(childObject, "v");
    cell.vth      = json_object_get_number(childObject, "vth");
    cell.vahp     = json_object_get_number(childObject, "vahp");
    cell.vmax     = json_object_get_number(childObject, "vmax");
    cell.t_refrac = json_object_get_number(childObject, "t_refrac");
    
    // initial parameters
    cell.t0 = -100;
    cell.is_refrac = 0;
    cell.flag = 0;
    cell.i = 0;

    return cell;

}

ExpSyn readSynParams(JSON_Object *childObject){

    ExpSyn syn;

    syn.tau1    = json_object_get_number(childObject, "tau1");
    syn.tau2    = json_object_get_number(childObject, "tau2");
    syn.A       = json_object_get_number(childObject, "A");
    syn.gmax    = json_object_get_number(childObject, "gmax");
    syn.e       = json_object_get_number(childObject, "e");
    syn.d       = json_object_get_number(childObject, "d");
    syn.id_pre  = (int) json_object_get_number(childObject, "id_pre");
    syn.id_post = (int) json_object_get_number(childObject, "id_post");

    // initial parameters
    syn.g = 0;

    return syn;

}

// read Poisson input
Pinput readPinputParams(JSON_Object *childObject){

    Pinput stim;
    double f;

    stim.tstart = json_object_get_number(childObject, "tstart");
    stim.tend   = json_object_get_number(childObject, "tend");

    // calc probability
    f = json_object_get_number(childObject, "f");
    stim.p = f*_dt / 1000;

    // init
    stim.t = -100;

    return stim;

}

/* Run Network */
void updateLIF(LIFneuron *cell){
    // update
    if (cell->is_refrac == 1) {
        // if flag is 1 or 2 -> end update
        if (cell->flag == 1) {
            cell->v = cell->vmax;
            cell->flag = 2;
            return;
        } else if (cell->flag == 2) {
            cell->v = cell->vahp;
            cell->flag = 0;
            return;
        }
        if ((t - cell->t0) > cell->t_refrac) {
            cell->is_refrac = 0;
        }
    }
    cell->v += solveLIF(*cell) + genrand64_normal(0, _s);
    // detect isfire
    if (cell->v > cell->vth) {
        cell->t0 = t;
        cell->is_refrac = 1;
        cell->flag = 1;
    }
    // empty current
    cell->i = 0;
}


void updateSyn(ExpSyn *syn){
    /*
    equation: Isyn = gmax * P * (v - e)
    P = A * (exp(-((t-t0-d) / tau1)) - exp(-((t-t0-d) / tau2)))
    */
    double P, del_t;

    if (*(syn->ref_is_refrac) == 0) {
        // printf("\n%lf, %lf\n", t, *(syn->ref_t0));
        del_t = t - *(syn->ref_t0) - syn->d;
        if ((del_t > 0) & (del_t < 50)) {
            // calculate P
            P = (syn->A) * (exp(-del_t/syn->tau1) - exp(-del_t/syn->tau2));
            *(syn->ref_i) -= (syn->gmax) * P * (*(syn->ref_v) - syn->e);
        }
    }
}

void updateStim(Pinput *stim){
    double p;

    if ((t>=stim->tstart) & (t<=stim->tend)){
        p = genrand64_real2();
        if (p < stim->p){
            stim->t = t;
        }
    }
}

void updateAll(){

    int n;

    t += _dt;

    for (n=0; n<NumofPinputs; n++){
        updateStim(stims+n);
    }

    for (n=0; n<NumofSyns; n++){
        updateSyn(syns+n);
    }

    for (n=0; n<NumofCells; n++){
        updateLIF(cells+n);
    }

}


double solveLIF(LIFneuron cell){
    double dv;
    // use RK4 method
    // double dv1, dv2, dv3, dv4;

    // dv1 = fLIF(cell, cell.v)*_dt;
    // dv2 = fLIF(cell, cell.v+dv1/2)*_dt;
    // dv3 = fLIF(cell, cell.v+dv2/2)*_dt;
    // dv4 = fLIF(cell, cell.v+dv3)*_dt;

    // return (dv1+2*dv2+2*dv3+dv4)/6;

    // user Euler method
    dv = fLIF(cell, cell.v) * _dt;
    return dv;
}


double fLIF(LIFneuron cell, double v){ // update LIF neuron
    // equation: V' = [e - v + r * i] / tau
    return (cell.e - v + cell.r*cell.i) / cell.tau;
}

void freeObjs(){
    free(cells);
    free(syns);
    free(stims);
}
