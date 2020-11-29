#ifndef _ARTNET
#define _ARTNET

#include "parson.h"

typedef struct lif_neuron LIFneuron;
typedef struct exp_syn ExpSyn;
typedef struct poisson_input Pinput;
typedef struct i_clamp IClamp;


/* read parameters */
void readNetworkInfo(char JSON_fname[]);
LIFneuron readCellParams(JSON_Object *childObject);
ExpSyn readSynParams(JSON_Object *childObject);
Pinput readPinputParams(JSON_Object *childObject);

/* running Network */
void updateLIF(LIFneuron *cell);
void updateSyn(ExpSyn *syn);
void updateStim(Pinput *stim);
void updateAll();
double solveLIF(LIFneuron cell);
double fLIF(LIFneuron cell, double v);
void freeObjs();

#endif
