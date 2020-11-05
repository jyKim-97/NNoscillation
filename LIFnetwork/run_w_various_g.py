import numpy as np
import matplotlib.pyplot as plt
import SimpleModels as nrn
import argparse
import pickle


def get_gext(f, p_pos, w_rand=100, t0=None, t1=None, tau_r=0.1, tau_d=5, d=0):
    if t0 is None:
        t0 = 0
    if t1 is None:
        t1 = nrn._tmax
    # f Hz
    i0 = int(t0 / nrn._dt)
    i1 = int(t1 / nrn._dt)
    id_onsets = []
    while i0 < i1:
        id_onsets.append(i0 + np.random.randint(-w_rand//2, w_rand//2))
        i0 += 1e3/f/nrn._dt
    # random values
    p_bd = p_pos * nrn._dt
    rands = np.random.uniform(low=0, high=1, size=nrn._nitr)
    for i in np.where(rands < p_bd)[0]:
        if (nrn._times[i]<t0) | (nrn._times[i]>t1):
            id_onsets.append(i)
    id_onsets = np.array(id_onsets, dtype=np.int)
    id_onsets = np.sort(id_onsets)
    # get ghat
    n = 0
    onset = 100
    tau1, tau2 = nrn.convert_tau_rdto12(tau_r, tau_d)
    B = nrn.get_syn_norm(tau_r, tau_d)
    g = np.zeros(nrn._nitr)
    for i, t in enumerate(nrn._times[1:]):
        if n < len(id_onsets):
            if i == id_onsets[n]:
                onset = nrn._times[i]
                n += 1
        if ((t - onset) >= d) and ((t - onset) < tau_d*5):
            g[i] = B * (np.exp(-(t-onset-d)/tau1) - np.exp(-(t-onset-d)/tau2))
    return g


def set_ext_params(n_excs, n_pfcs_e, n_pfcs_i, n_pick, n_overlap, cell_types, gbar_m, gbar_s):
    
    ######## gbar_exc # gbar_pfc_exc # gbar_pfc_inh ######
    # to PN           #              #              ######
    ######################################################
    # to PV           #              #              ######
    ######################################################
    
    target_ids_all = select_targets(n_excs+n_pfcs_e+n_pfcs_i, len(cell_types), n_pick, n_overlap)
    gbar_all = []
    es_all = []
    
    # region X
    for i in range(n_excs):
        gbar_all.append([])
        for j in target_ids_all[i]:
            if cell_types[j] == 0:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[0][0], scale=gbar_s[0][0])))
            else:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[1][0], scale=gbar_s[1][0])))
        es_all.append(0)
    
    # PFC signal to exc
    for i in range(n_pfcs_e):
        gbar_all.append([])
        for j in target_ids_all[i]:
            if cell_types[j] == 0:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[0][1], scale=gbar_s[0][1])))
            else:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[1][1], scale=gbar_s[1][1])))
        es_all.append(0)
        
    # PFC signal to inh
    for i in range(n_pfcs_i):
        gbar_all.append([])
        for j in target_ids_all[i]:
            if cell_types[j] == 0:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[0][2], scale=gbar_s[0][2])))
            else:
                gbar_all[-1].append(abs(np.random.normal(loc=gbar_m[1][2], scale=gbar_s[1][2])))
        es_all.append(-80)
    
    return target_ids_all, gbar_all, es_all


def create_random_ntk(n_exc, n_inh, probs):
    n_cells = n_exc + n_inh
    cell_types = []
    cnt_map = np.ones([n_cells, n_cells], dtype=np.int) * (-1)
    gtype = [[0, 1], [2, 3]]
    for i in range(n_cells):
        if i < n_exc:
            cell_types.append(0)
        else:
            cell_types.append(1)
            
    for i in range(n_cells):
        cpre = cell_types[i]
        for j in range(n_cells):
            if i != j:
                cpost = cell_types[j]
                p = np.random.uniform(0, 1)
                if p <= probs[cpre][cpost]:
                    cnt_map[i][j] = gtype[cpre][cpost]
    return cell_types, cnt_map


def select_targets(n_ext, n_target, n_pick, n_overlap):
    
    targets = []
    for i in range(n_target):
        for j in range(n_overlap):
            targets.append(i)
            targets.append(i)
        
    target_ids = []
    for i in range(n_ext):
        # pick n_pick cells
        target_ids.append([np.random.choice(targets)])
        targets.remove(target_ids[-1][0])
        for j in range(1, n_pick):
            while True:
                n = np.random.choice(targets)
                if n != target_ids[-1][j-1]:
                    target_ids[-1].append(n)
                    targets.remove(n)
                    break                    
    return target_ids


params_pn = {'tau':20, 'r':100, 'vth':-40, 'v0':-65, 'vahp':-80, 'vmax':30, 'ev':-65, 'tahp':10}
params_pv = {'tau':5, 'r':100, 'vth':-50, 'v0':-65, 'vahp':-80, 'vmax':30, 'ev':-65, 'tahp':0}


syn_pn2pn = {'gbar_syn':0.004, 'tau_r':0.1, 'tau_d':8, 'es':0}
syn_pn2pv = {'gbar_syn':0.002, 'tau_r':0.1, 'tau_d':2, 'es':0}
syn_pv2pn = {'gbar_syn':0.004, 'tau_r':0.1, 'tau_d':4, 'es':-80}
syn_pv2pv = {'gbar_syn':0.001, 'tau_r':0.1, 'tau_d':4, 'es':-80}


parser = argparse.ArgumentParser(description='get seed')
parser.add_argument('--seed', required=True, type=int)

args = parser.parse_args()

if __name__ == '__main__':
    nrn.set_times(tmax=500, dt=0.01)
    
    nrn.set_seed(args.seed)
    fdir = './data/w_pfc/'

    n_exc = 80
    n_inh = 20
    n_cells = n_exc + n_inh
    probs = [[0.2, 0.2], [0.8, 0.2]]
    # create cell types
    cell_types, cnt_map = create_random_ntk(n_exc, n_inh, probs)
    setting = nrn.get_params([params_pn, params_pv],
                            [syn_pn2pn, syn_pn2pv, syn_pv2pn, syn_pv2pv],
                            cell_types, cnt_map, delay_m=1, delay_s=0.1)

    n_alls = 100
    n_pfcs_exc = 50
    n_pfcs_inh = 50
    # create background excitatory input; g_excs, es
    g_excs_all = np.zeros([nrn._nitr, n_alls+n_pfcs_exc+n_pfcs_inh])
    for i in range(n_alls):
        g_excs_all[:, i] = nrn.gPoisson(0.1, 0.1, 5, delay=0, t0=0, t1=500)
    for i in range(n_pfcs_exc + n_pfcs_inh):
        g_excs_all[:, n_alls+i] = get_gext(30, 0.01, 1000, t0=50, t1=200, tau_r=0.1, tau_d=5)

    gbar_pfc_exc = np.linspace(0, 0.001, 20)
    gbar_pfc_inh = np.linspace(0, 0.001, 20)
    gbar_x = 0.004

    # save info file
    with open(fdir+f'{args.seed}_ntk_info.pickle', 'wb') as fid:
        pickle.dump(cnt_map, fid)
        pickle.dump(g_excs_all, fid)
    
    # ntk_w_pfcs
    for i, g_exc in enumerate(gbar_pfc_exc): # row
        for j, g_inh in enumerate(gbar_pfc_inh): # columnm
            print(args.seed, i, j)

            gbar_m = [[gbar_x, g_exc, g_inh], [gbar_x/10, g_exc/10, g_inh/10]]
            gbar_s = [[gbar_x/10, g_exc/10, g_inh/10], [gbar_x/100, g_exc/100, g_inh/100]]

            target_ids_all, gbar_all, es_all = set_ext_params(n_alls, n_pfcs_exc, n_pfcs_inh, 2, 4, cell_types, gbar_m, gbar_s)

            ntk_w_pfcs = nrn.CellNetwork(**setting, g_ext=g_excs_all, e_ext=es_all, gbar_ext=gbar_all, target_id=target_ids_all, std=0.1)
            ntk_w_pfcs.run()

            with open(fdir+f'{args.seed}_ntk_calc_vals_%d_%d.pickle'%(i, j), 'wb') as fid:
                pickle.dump(ntk_w_pfcs.vcells, fid)
                pickle.dump(ntk_w_pfcs.spks, fid)

