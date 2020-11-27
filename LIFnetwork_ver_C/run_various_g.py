import numpy as np
import PrintNNinfo as out
import os
import argparse
import pickle
import subprocess


def CreateERnetwork(cell_types, ps):
    n = len(ps)
    syn_types = np.arange(n*n).reshape(n, n)
    cnt_map = np.ones([len(cell_types), len(cell_types)], dtype=int) * (-1)
    
    for i, cpre in enumerate(cell_types):
        for j, cpost in enumerate(cell_types):
            if i != j:
                if np.random.uniform() < ps[cpre][cpost]:
                    cnt_map[i][j] = syn_types[cpre][cpost]
    return cnt_map

def get_Poisson_t(f, N, tmax, dt, t0=0):
    prob = f*dt
    times = np.arange(t0, tmax, dt)
    ts = []
    for i in range(N):
        rand = np.random.uniform(low=0, high=1e3, size=len(times))
        ids = rand < prob
        ts.append(times[ids])
    return ts

def select_targets(n_stims, n_cells, n_overlap):
    targets = []
    for i in range(n_overlap):
        targets = np.concatenate((targets, [i for i in range(n_cells)]))
    targets = list(targets.astype(np.int))
        
    target_ids = []
    for i in range(n_stims):
        tid = [np.random.choice(targets)]
        targets.remove(tid)
        for j in range(n_overlap-1):
            n = np.random.choice(targets)
            while n in tid:
                n = np.random.choice(targets)
            tid.append(n)
            targets.remove(n)
        target_ids.append(tid)
    return target_ids

def set_ext_types(target_ids, n_ext_exc, cell_types, ext_syn_types):
    n_ext_inh = len(target_ids) - n_ext_exc
#     ext_syn_types = [[0, 1], [2, 3]] # ext_PN -> ? / ext_Inh -> ?
    ext_types = []
    for i, tid in enumerate(target_ids):
        ext_types.append([])
        if i < n_ext_exc:
            pre = 0
        else:
            pre = 1
        for n in tid:
            ext_types[-1].append(ext_syn_types[pre][cell_types[n]])
    return ext_types

# params
params_pn = {'tau':20, 'r':500, 'e':-70, 'vahp':-80, 'vth':-40, 't_refrac':5, 'v0':-65, 'vmax':30}
params_pv = {'tau':2, 'r':500, 'e':-70, 'vahp':-80, 'vth':-40, 't_refrac':0.2, 'v0':-65, 'vmax':30}

# set synapse type & parameter
syn_pn2pn = {'gmax':5e-3, 'tau_r':0.3, 'tau_d':6.9, 'e':0, 'd':1}
syn_pn2pv = {'gmax':7e-3, 'tau_r':0.1, 'tau_d':2.4, 'e':0, 'd':1}
syn_pv2pn = {'gmax':7.2e-3, 'tau_r':0.5, 'tau_d':6.8, 'e':-80, 'd':1}
syn_pv2pv = {'gmax':4e-3, 'tau_r':0.5, 'tau_d':6.8, 'e':-80, 'd':1}
syn_type_params = [syn_pn2pn, syn_pn2pv, syn_pv2pn, syn_pv2pv]

# stimulation
ext_pn2pn = {'gmax':0.05, 'tau_r':0.3, 'tau_d':6, 'e':0, 'd':0}
ext_pn2pv = {'gmax':0.05, 'tau_r':0.3, 'tau_d':6, 'e':0, 'd':0}
ext_pv2pn = {'gmax':0.05, 'tau_r':0.5, 'tau_d':6, 'e':-80, 'd':0}
ext_pv2pv = {'gmax':0.05, 'tau_r':0.5, 'tau_d':6, 'e':-80, 'd':0}
ext_syn_type_params = [ext_pn2pn, ext_pn2pv, ext_pv2pn, ext_pv2pv]


parser = argparse.ArgumentParser(description='get seed')
parser.add_argument('--seed', required=True, type=int)
args = parser.parse_args()


if __name__ == '__main__':

    # set time parameter
    dt = 0.01
    tmax = 5000

    n_exc = 1000
    n_inh = 250
    n_stims = 1250
    n_ext_exc = 1000
    n_cells = n_exc + n_inh

    ps = [[0.2, 0.2], [0.3, 0.3]] # PN->PN/Inh, Inh->PN/Inh
    f_pos = 1

    fdir = './params_g/'
    prefix = 'ntk_%d'%(args.seed)
    fsave ='./result_g/out_%d'%(args.seed)

    print('')

    # set network parameters
    np.random.seed(args.seed)
    
    cell_types = []
    for i in range(n_cells):
        if i < n_exc:
            cell_types.append(0)
        else:
            cell_types.append(1)

    cnt_map = CreateERnetwork(cell_types, ps)

    t_stims = get_Poisson_t(f_pos, n_stims, tmax, dt, t0=0)
    target_ids = select_targets(n_stims, n_cells, 3)
    ext_types = set_ext_types(target_ids, n_ext_exc, cell_types, [[0, 1], [2, 3]])

    # run network
    subprocess.call(['./runNetwork_no_pbar',
                    '--fdir', './parameter/',
                    '--prefix', 'ntk_test',
                    '--fsave', './out',
                    '--tmax', '500',
                    '--dt', '0.01',
                    '--s', '0.15',
                    '--seed', '1000'])

    # load ntk
    times, vcells = out.readOut("./vcells_all_v.csv")
    _, icells = out.readOut("./vcells_all_i.csv")
    
    with open(fsave+'.pickle', 'wb') as fid:
        pickle.dump(times, fid)
        pickle.dump(vcells, fid)
        pickle.dump(icells, fid)
