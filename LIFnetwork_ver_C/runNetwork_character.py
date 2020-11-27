import numpy as np
import matplotlib.pyplot as plt
import os, subprocess
import PrintNNinfo as out
from time import time
import pickle as pkl


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


def read_all(prefix):
    t, v = out.readOut(prefix+'_v.csv')
    _, ii = out.readOut(prefix+'_i.csv')
    tspks = []
    for i in range(v.shape[1]):
        tspks.append(t[v[:, i] == 30])
    return t, v, ii, tspks


def merge_lists(a1, a2):
    arrs = []
    for a in a1:
        arrs.append(a)
    for a in a2:
        arrs.append(a)
    return arrs


def auto_corr(x, tmax, srate, ind=None):
    if ind is not None:
        x = x[ind]
    kmax = int(tmax * srate)
    corr = []
    ks = np.arange(0, kmax+0.1, 1) / srate
    xm = np.average(x)
    for k in range(kmax+1):
        if k == 0:
            x1 = x
            x2 = x
        else: # k > 0
            x1 = x[:-k]
            x2 = x[k:]
        s1 = np.std(x1)
        s2 = np.std(x2)
        if ((s1==0) or (s2==0)):
            corr.append(0)
        else:
            corr.append(np.sum((x1-xm)*(x2-xm))/len(x)/(np.std(x1)*np.std(x2)))
    return ks, corr

def get_sync_index(vall, ind=None):
    # vall (:, n_cells)
    if ind is None:
        ind = np.ones(vall.shape[0], dtype=bool)
    
    v_avg = np.average(vall, axis=1)
    s_all = np.std(v_avg[ind])**2
    s_cell = np.std(vall[ind, :], axis=0)
    
    numer = np.average(s_all)
    denom = np.average(s_cell**2)
    
    return numer/denom

def runSimulation(info):

    tic = time()
    seed = info['seed']

    np.random.seed(seed)

    # set parameters
    tmax = 200
    dt = 0.01

    n_exc = 100
    n_inh = 25

    n_stims_pos = 125
    p_stims_pos = [0.2, 0.2]

    # Poisson input
    f_pos = info['f']

    out_prefix = 'out%d'%(info['id'])
    fname = 'id_%04d_f_%03d_ext_%2d-5_seed_%d'%(info['id'], info['f'], info['gmax'], seed)

    # set PN cell type & cells parameter
    params_pn = {'tau':20, 'r':500, 'e':-70, 'vahp':-80, 'vth':-40, 't_refrac':5, 'v0':-65, 'vmax':30}
    params_pv = {'tau':2, 'r':500, 'e':-70, 'vahp':-80, 'vth':-40, 't_refrac':0.2, 'v0':-65, 'vmax':30}

    # set synapse type & parameter
    syn_pn2pn = {'gmax':1e-3, 'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':1}
    syn_pn2pv = {'gmax':1e-4, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':1}
    syn_pv2pn = {'gmax':1e-2, 'tau_r':0.5, 'tau_d':7.6, 'e':-80, 'd':1}
    syn_pv2pv = {'gmax':1e-3, 'tau_r':0.5, 'tau_d':7.6, 'e':-80, 'd':1}
    syn_type_params = [syn_pn2pn, syn_pn2pv, syn_pv2pn, syn_pv2pv]

    # set network parameters
    ps = [[0.2, 0.2], [0.3, 0.3]] # PN->PN/Inh, Inh->PN/Inh

    # background input
    ext_pn2pn = {'gmax':info['gmax']*1e-5, 'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':0}
    ext_pn2pv = {'gmax':info['gmax']*1e-5, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':0}

    ext_syn_type_params = [ext_pn2pn, ext_pn2pv]

    # create network
    n_cells = n_exc + n_inh
    cell_types = []
    for i in range(n_cells):
        if i < n_exc:
            cell_types.append(0)
        else:
            cell_types.append(1)
    cnt_map = CreateERnetwork(cell_types, ps)

    # get ids
    ids_targets = []
    for i in range(n_stims_pos):
        ids_targets.append([])
        for j in range(n_cells):
            c = cell_types[j]
            p = np.random.uniform(low=0, high=1)
            if p < p_stims_pos[c]:
                ids_targets[i].append(j)

    # external types
    ext_types = []
    for i in range(n_stims_pos):
        ext_types.append([])
        for ind in ids_targets[i]:
            c = cell_types[ind]
            ext_types[i].append(c)

    t_stims = get_Poisson_t(f_pos, n_cells, tmax, dt, t0=0)


    out.print_nn_params('./parameter/', out_prefix, cell_types, cnt_map, [params_pn, params_pv], syn_type_params,
                        ids_targets, ext_types, t_stims, ext_syn_type_params, overwrite=True, ischeck=False,
                        d_noise=0.3, g_noise_syn=1e-5, g_noise_ext=0)
    
    subprocess.call(['./runNetwork_w_no_print',
                    '--fdir', './parameter/', '--prefix', out_prefix, '--fsave', f'./vcells_{out_prefix}', '--tmax', '%d'%(tmax),
                     '--dt', '0.01', '--s', '0.2', '--seed', '%d'%(seed)])

    times, vall, iall, tspks = read_all(f'./vcells_{out_prefix}')

    si = get_sync_index(vall)

    with open('./result_ch/data/'+fname+'.pkl', 'wb') as f:
            pkl.dump(vall, f)
            pkl.dump(times, f)

            pkl.dump(tspks, f)
            pkl.dump(si, f)
    
    os.remove(f'./parameter/{out_prefix}_info.csv')
    os.remove(f'./parameter/{out_prefix}_cell.csv')
    os.remove(f'./parameter/{out_prefix}_syn.csv')
    os.remove(f'./parameter/{out_prefix}_t_spike.csv')

    os.remove(f'./vcells_{out_prefix}_v.csv')
    os.remove(f'./vcells_{out_prefix}_i.csv')

    print('%d Done, t=%.3fs'%(info['id'], time()-tic))

    return


if __name__ == '__main__':

    import multiprocessing

    plt.ioff()

    f_range = np.linspace(1, 201, 20)
    gmax_range = np.linspace(5, 50, 20) #x1e-5
    n_mc = 5
    seeds = np.random.randint(10, 100000, n_mc)
    np.save('./result_ch/seed.npy', seeds)

    infos = []
    n = 0
    for seed in seeds:
        for f in f_range:
            for g in gmax_range:
                    infos.append({'id':n, 'f':f, 'gmax':g, 'seed':seed})
                    n += 1

    pool = multiprocessing.Pool(processes=18)

    tic = time()
    pool.map(runSimulation, infos)
    pool.close()
    pool.join()

    print('Execution Done, %.5fs'%(time()-tic))


        
        
        

    
        
        

        