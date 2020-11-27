import numpy as np
import matplotlib.pyplot as plt
import os, subprocess
import PrintNNinfo as out
from tqdm import tqdm
from time import time
import pickle as pkl
from scipy.ndimage import gaussian_filter1d
from scipy.signal import resample
from scipy.interpolate import interp1d


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


def getFFT(x, dt, fmax=None, idt=None):
    if idt is not None:
        x = x[idt]
    n = len(x)
    fx = np.fft.fft(x, n) / n
    f = np.fft.fftfreq(n, dt)
    fx = abs(fx[:int(n//2)])
    f = f[:int(n//2)]
    if fmax is not None:
        idf = f < fmax
        f = f[idf]
        fx = fx[idf]
    return f, fx


def runSimulation(info):

    tic = time()
    # info, 'f', 'id', 'gmax_to_pn', 'gmax_to_pv'
    # print('%d Start'%(info['id']))
    seed = info['seed']

    np.random.seed(seed)

    out_prefix = 'out%d'%(info['id'])
    fname = 'id_%03d_pfc_%02d_%02d_-5_ratio_%d'%(info['id'], info['gmax_pfc_pn'], info['gmax_pfc_pv'], info['ratio'])
    fnames = os.listdir('./results/data/')

    if fname+'.pkl' in fnames:    
        return
    

    # set parameters
    tmax = 500
    dt = 0.01

    n_exc = 100
    n_inh = 25

    n_stims_pfc_pn = 100
    n_stims_pfc_pv = 25
    p_stims_pfc = [0.1, 0.1, 0.1, 0.1]

    n_stims_pos = 125
    p_stims_pos = [0.2, 0.2]

    # Poisson input
    f_pos = 80

    # oscillatory input
    t0 = 200
    t1 = 500

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

    # PFC input 3, 1
    pfc_pn2pn = {'gmax':info['gmax_pfc_pn']*1e-5,               'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':0}
    pfc_pn2pv = {'gmax':info['gmax_pfc_pn']/info['ratio']*1e-5, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':0}
    pfc_pv2pn = {'gmax':info['gmax_pfc_pv']*1e-5,               'tau_r':0.5, 'tau_d':7.6, 'e':-80, 'd':0}
    pfc_pv2pv = {'gmax':info['gmax_pfc_pv']/info['ratio']*1e-5, 'tau_r':0.5, 'tau_d':7.6, 'e':-80, 'd':0}

    # background input
    ext_pn2pn = {'gmax':3e-4, 'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':0}
    ext_pn2pv = {'gmax':3e-4, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':0}

    ext_syn_type_params = [pfc_pn2pn, pfc_pn2pv, pfc_pv2pn, pfc_pv2pv, ext_pn2pn, ext_pn2pv]

    

    # ---------------------------------------- #
    # -------- Calculate network info -------- #
    # ---------------------------------------- #

    # create network
    n_cells = n_exc + n_inh
    cell_types = []
    for i in range(n_cells):
        if i < n_exc:
            cell_types.append(0)
        else:
            cell_types.append(1)
    cnt_map = CreateERnetwork(cell_types, ps)

    n_stims_pfc = n_stims_pfc_pn+n_stims_pfc_pv
    # create t_stims & id at 28 Hz
    t_stims_pfc = get_Poisson_t(1, n_stims_pfc, tmax, dt, t0=0)
    for i in range(n_stims_pfc):
        t = t0
        tmp = []
        while t < t1:
            val = t + np.random.normal(loc=0, scale=1)
            if val < t1:
                tmp.append(val)
            t += 1e3/28
        t_stims_pfc[i] = np.concatenate((t_stims_pfc[i], tmp))
        t_stims_pfc[i] = np.sort(t_stims_pfc[i])
    t_stims_pos = get_Poisson_t(f_pos, n_stims_pos, tmax, dt, t0=0)

    t_stims = merge_lists(t_stims_pfc, t_stims_pos)

    # get ids
    ids_targets = []
    for i in range(n_stims_pfc+n_stims_pos):
        if i < n_stims_pfc_pn:
            prob = p_stims_pfc[:2]
        elif i < n_stims_pfc:
            prob = p_stims_pfc[2:]
        else:
            prob = p_stims_pos
        ids_targets.append([])
        for j in range(n_cells):
            c = cell_types[j]
            p = np.random.uniform(low=0, high=1)
            if p < prob[c]:
                ids_targets[i].append(j)

    # external types
    ext_types = []
    for i in range(n_stims_pfc + n_stims_pos):
        if i < n_stims_pfc_pn:
            n0 = 0
        elif i < n_stims_pfc:
            n0 = 1
        else:
            n0 = 2
        ext_types.append([])
        for ind in ids_targets[i]:
            c = cell_types[ind]
            ext_types[i].append(n0*2+c)
    
    out.print_nn_params('./parameter/', out_prefix, cell_types, cnt_map, [params_pn, params_pv], syn_type_params,
                        ids_targets, ext_types, t_stims, ext_syn_type_params, overwrite=True, ischeck=False,
                        d_noise=0.3, g_noise_syn=1e-5, g_noise_ext=0)
        
    
    subprocess.call(['./runNetwork_w_no_print',
                    '--fdir', './parameter/', '--prefix', out_prefix, '--fsave', f'./vcells_{out_prefix}', '--tmax', '%d'%(tmax),
                     '--dt', '0.01', '--s', '0.2', '--seed', '%d'%(seed)])

    times, vall, iall, tspks = read_all(f'./vcells_{out_prefix}')

    v_avg = np.average(vall, axis=1)
        
    freq_bef, fx_bef = getFFT(v_avg, dt/1e3, fmax=100, idt=(times<200))
    freq_aft, fx_aft = getFFT(v_avg, dt/1e3, fmax=100, idt=(times>200)&(times<500))

    f0 = 5
    f_new = np.linspace(10, 80, 1000)

    interp_bef = interp1d(freq_bef[freq_bef>f0], fx_bef[freq_bef>f0], kind='linear')
    interp_aft = interp1d(freq_aft[freq_aft>f0], fx_aft[freq_aft>f0], kind='linear')
    
    with open('./results/data/'+fname+'.pkl', 'wb') as f:
        pkl.dump(vall, f)
        pkl.dump(times, f)
        
        pkl.dump(interp_bef(f_new), f)
        pkl.dump(interp_aft(f_new), f)

        pkl.dump(f_new, f)

    os.remove(f'./parameter/{out_prefix}_info.csv')
    os.remove(f'./parameter/{out_prefix}_cell.csv')
    os.remove(f'./parameter/{out_prefix}_syn.csv')
    os.remove(f'./parameter/{out_prefix}_t_spike.csv')

    os.remove(f'./vcells_{out_prefix}_v.csv')
    os.remove(f'./vcells_{out_prefix}_i.csv')

    # plot dat
    # fig = plt.figure(dpi=200, figsize=(6, 6))

    # plt.subplot(411)

    # for i in range(n_cells):
    #     if i < n_exc:
    #         c = 'r'
    #     else:
    #         c = 'b'
    #     plt.vlines(tspks[i], i-0.5, i+0.5, color=c, lw=1)

    # plt.vlines(np.arange(0, 600, 50), 0, n_cells, color='k', lw=0.5)

    # plt.xticks([])
    # plt.yticks([])
    # plt.xlim([0, tmax])
    # plt.ylim([0, n_cells])

    # plt.subplot(412)

    # v_avg = np.average(vall, axis=1)
    # lfp = gaussian_filter1d(v_avg, 2)
    # plt.plot(times, lfp, lw=0.5, c='k')
    # plt.yticks([])
    # # plt.xticks(np.arange(0, 2100, 250))
    # plt.xticks(np.arange(0, tmax+50, 100))
    # plt.xlim([0, tmax])

    # plt.subplot(413)

    # srate = 2000
    # lfp_new = resample(lfp, srate//2+1)
    # t_new = np.arange(0, tmax+1e3/2/srate, 1e3/srate)
    # ks_aft, corr_aft = auto_corr(lfp_new, tmax=0.1, srate=srate, ind=(t_new>=t0)&(t_new<t1))
    # ks_bef, corr_bef = auto_corr(lfp_new, tmax=0.1, srate=srate, ind=(t_new<t0))
    # plt.plot(ks_bef*1e3, corr_bef, 'k', lw=1)
    # plt.plot(ks_aft*1e3, corr_aft, 'r', lw=1)
    # plt.vlines([1e3/28, 1e3/50, 1e3/60], -1, 1, 'g', linestyle='--', lw=1)

    # plt.xlim([0, 100])
    # plt.ylim([-1, 1])

    # plt.subplot(414)

    # plt.plot(f_new, interp_bef(f_new), 'k', lw=1)
    # plt.plot(f_new, interp_aft(f_new), 'r', lw=1)

    # plt.xlim([10, 80])
    # plt.ylim([0, 2])

    # plt.tight_layout()
    # plt.savefig(f'./results/figs/{fname}.png')
    
    # plt.close(fig)

    print('%d Done, t=%.3fs'%(info['id'], time()-tic))

    return

if __name__ == '__main__':

    import multiprocessing

    plt.ioff()

    # gmax_to_pn_range = np.linspace(5e-5, 5e-4, 10)
    # gmax_to_pv_range = np.linspace(5e-5, 5e-4, 10)
    gmax_to_pn_range = np.linspace(5, 50, 20) # x 1e-5
    gmax_to_pv_range = np.linspace(5, 50, 20) # x 1e-5
    # ratio = np.linspace(1, 10, 10) # PN/PV
    n_mc = 5
    seeds = np.random.randint(10, 100000, n_mc)
    np.save('./results/seed.npy', seeds)
    # seeds = np.load('./results/seed.npy')

    infos = []
    n = 0
    for i in range(n_mc):
        # seed = /
        for gmax_to_pn in gmax_to_pn_range:
            for gmax_to_pv in gmax_to_pv_range:
                infos.append({'id':n, 'gmax_pfc_pn':gmax_to_pn, 'gmax_pfc_pv':gmax_to_pv, 'ratio': 1, 'seed': seeds[i]})
                n += 1

    pool = multiprocessing.Pool(processes=18)

    tic = time()
    pool.map(runSimulation, infos)
    pool.close()
    pool.join()

    print('Execution Done, %.5fm'%((time()-tic)/60))
        