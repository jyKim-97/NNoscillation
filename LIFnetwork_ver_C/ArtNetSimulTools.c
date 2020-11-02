#include "ArtNetSimulTools.h"
#include "./mt64.h"
#include "./mt64.c"

double _tmax = 300;
double _dt = 0.01;
double _s = 0.1;


void readCellParams(FILE *fid, LIFneuron *cell){
    int tmp_int;

    fscanf(fid, "%d,", &tmp_int);
    // tau, r, e, vth, vahp, v0, vmax, t_refrac
    fscanf(fid, "%lf,", &(cell->tau));
    fscanf(fid, "%lf,", &(cell->r));
    fscanf(fid, "%lf,", &(cell->e));
    fscanf(fid, "%lf,", &(cell->vth));
    fscanf(fid, "%lf,", &(cell->vahp));
    fscanf(fid, "%lf,", &(cell->v));
    fscanf(fid, "%lf,", &(cell->vmax));
    fscanf(fid, "%lf", &(cell->t_refrac));

    // init params
    cell->t0 = -100;
    cell->is_refrac = 0;
    cell->flag = 0;
    cell->i = 0;

}


void readSynParams(FILE *fid, ExpSyn *syn){
    int tmp_int;

    fscanf(fid, "%d,", &tmp_int);
    // tau1, tau2, A, gmax, e, d, id_pre, id_post
    fscanf(fid, "%lf,", &(syn->tau1));
    fscanf(fid, "%lf,", &(syn->tau2));
    fscanf(fid, "%lf,", &(syn->A));
    fscanf(fid, "%lf,", &(syn->gmax));
    fscanf(fid, "%lf,", &(syn->e));
    fscanf(fid, "%lf", &(syn->d));

    // init params
    syn->g = 0;

}


void updateLIF(LIFneuron *cell, double t){
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


void updateSyn(ExpSyn *syn, double t){
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


void updateStim(Stim *stim, double *spk_times, double t){
    // printf("stim t =%lf\n", *(stim->ref_t0));
    if (stim->n < stim->len && spk_times[stim->n+1] <= t){
        stim->n++;
        stim->t0 = spk_times[stim->n];
    }
}

double solveLIF(LIFneuron cell){
    // use RK4 method
    double dv1, dv2, dv3, dv4;

    dv1 = fLIF(cell, cell.v)*_dt;
    dv2 = fLIF(cell, cell.v+dv1/2)*_dt;
    dv3 = fLIF(cell, cell.v+dv2/2)*_dt;
    dv4 = fLIF(cell, cell.v+dv3)*_dt;

    return (dv1+2*dv2+2*dv3+dv4)/6;
}


double fLIF(LIFneuron cell, double v){ // update LIF neuron
    // equation: V' = [e - v + r * i] / tau
    return (cell.e - v + cell.r*cell.i) / cell.tau;
}


// test function
void runOnecell_w_current(LIFneuron *cell, double t0, double t1, double amp, double **vs){
    // for testing
    double t = 0;
    int n = 0;

    *(vs[0]) = cell->v;
    while (t < _tmax) {
        if ((t>t0) & (t<t1)){
            if (cell->is_refrac == 0) {
                cell->i = amp;
            }
        }
        updateLIF(cell, t);
        (*vs)[n+1] = cell->v;
        t += _dt;
        n++;
    }
}

// subfunctions
void saveData(double *vs, int n, FILE *fid){
    for (int i=0; i<n; i++){
        fprintf(fid, "%f,", vs[i]);
    }
}
