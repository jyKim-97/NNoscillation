import numpy as np
import pandas as pd
import os
import time


def get_syn_type_params(syn_type):
    tau1, tau2 = convert_tau_rdto12(syn_type['tau_r'], syn_type['tau_d'])
    A = get_syn_norm(syn_type['tau_r'], syn_type['tau_d'])
    syn_type['tau1'] = tau1
    syn_type['tau2'] = tau2
    syn_type['A'] = A


def convert_tau_rdto12(tau_r, tau_d):
    tau1 = tau_d
    tau2 = tau_r*tau_d / (tau_r + tau_d)
    return tau1, tau2


def get_syn_norm(tau_r, tau_d):
    tau1, tau2 = convert_tau_rdto12(tau_r, tau_d)
    # nsyns = len(np.array(tau1))
    dt = 0.01
    t = np.arange(0, 1000+dt/2, dt)
    val = np.exp(-t/tau1) - np.exp(-t/tau2)
    return 1 / max(val)


def print_nn_params(fdir, prefix,
                    cell_types, cnt_map, cell_type_params, syn_type_params, 
                    ids_target=None, extern_types=None, spk_times=None, extern_syn_type_params=None, d_noise=0, g_noise_syn=0, g_noise_ext=0, overwrite=False, ischeck=True):

    # check parameter
    n_type_cell_max = max(cell_types)
    n_type_syn_max = np.max(cnt_map)
    n_type_max = -100
    if n_type_cell_max != len(cell_type_params)-1:
        print("check the # of cell type parameters")
        return
    if n_type_syn_max != len(syn_type_params)-1:
        print("check the # of syn type parameters")
        return
    if ids_target is not None:
        for etype in extern_types:
            for j in etype:
                if  j > n_type_max:
                    n_type_max = j
        if n_type_max != len(extern_syn_type_params)-1:
            print("check the # of extern type parameters")
            return

    # check does the file exist
    fnames = [
        prefix+'_cell.csv', prefix+'_syn.csv', prefix+'_t_spike.csv'
    ]

    if ischeck:
        for fname in fnames:
            ch = check_f_exist(fdir, fname)
            if ch:
                print("%s exists"%(os.path.join(fdir, fname)))    
                if not overwrite:
                    return
            # else:
            #     os.remove(os.path.join(fdir, fname))
    
    # save ext spikes
    ids_save = []
    extern_saves = []
    if ids_target is not None:
    # save spike times
        with open(os.path.join(fdir, prefix+'_t_spike.csv'), 'w') as fid:
            n = 0
            fid.write(",len,spike_times\n")
            for i in range(len(ids_target)):
                if (len(spk_times[i]) > 0):
                    for _ in ids_target[i]:
                        fid.write("%d"%(n))
                        fid.write(",%d"%(len(spk_times[i])))
                        for t in spk_times[i]:
                            fid.write(",%5.3f"%(t))
                        fid.write("\n")

                    ids_save.append(ids_target[i])
                    extern_saves.append(extern_types[i])
                    n = n + 1

    # convert params
    ncells, df_cell = convert_cell_params(cell_types, cell_type_params)
    nsyns, df_syn1, id_pres, id_posts = convert_syn_params(cnt_map, syn_type_params, is_syn=True, d_noise=d_noise, g_noise=g_noise_syn)
    if ids_target is not None:
        nexts, df_syn2, _, _ = convert_syn_params(extern_saves, extern_syn_type_params, is_syn=False, g_noise=g_noise_ext)
        df_syn = pd.concat([df_syn1, df_syn2])
    else:
        nexts=0
        df_syn = df_syn1

    # save info file
    with open(os.path.join(fdir, prefix+'_info.csv'), 'w') as fid:
        fid.write(time.strftime("%Y-%m-%d-%H-%M-%S")+'\n')
        fid.write('ncells %d\n'%(ncells))
        fid.write('nsyns %d\n'%(nsyns))
        fid.write('nexts %d\n'%(nexts))
        # write connectivity, nsyns
        fid.write('pre/post synaptic neuron id, target neuronn id\n')
        for i in range(nsyns):
            fid.write('%d,%d,%d\n'%(i, id_pres[i], id_posts[i]))
        
        if ids_target is not None:
            i=0
            for target in ids_save:
                for n in target:
                    fid.write('%d,%d\n'%(i, n))
                    i += 1

    # save file
    df_cell.to_csv(os.path.join(fdir, prefix+'_cell.csv'), sep=',', mode='w', header=True, float_format="%.6f", index=True)
    df_syn.to_csv(os.path.join(fdir, prefix+'_syn.csv'), sep=',', mode='w', header=True, float_format="%.6f", index=True)
    


    # print("done\n");


def convert_cell_params(cell_types, cell_type_params):
    # initialize parameter dictionaries
    cell_params = dict((
        ('tau', []), ('r', []), ('e', []), ('vth', []), ('vahp', []), ('v0', []), ('vmax', []), ('t_refrac', [])
    ))
    ncells = len(cell_types)

    for i in range(ncells):
        ctype = cell_types[i]
        cell_params['tau'].append(cell_type_params[ctype]['tau'])
        cell_params['r'].append(cell_type_params[ctype]['r'])
        cell_params['e'].append(cell_type_params[ctype]['e'])
        cell_params['vth'].append(cell_type_params[ctype]['vth'])
        cell_params['v0'].append(cell_type_params[ctype]['v0'])
        cell_params['vahp'].append(cell_type_params[ctype]['vahp'])
        cell_params['vmax'].append(cell_type_params[ctype]['vmax'])
        cell_params['t_refrac'].append(cell_type_params[ctype]['t_refrac'])

    # print params    
    df_cell = pd.DataFrame(cell_params)
    df_cell = df_cell.astype(float)
    # df_cell.to_csv(os.path.join(fdir, prefix+'_cell.csv'), sep=',', mode='w', header=True, float_format="%.6f")
    return ncells, df_cell


def convert_syn_params(cnt_map, syn_type_params, is_syn=False, d_noise=0, g_noise=0):
    syn_params = dict((
        ('tau1', []), ('tau2', []), ('A', []), ('gmax', []), ('e', []), ('d', [])
    ))
    for syn in syn_type_params:
        get_syn_type_params(syn)

    nsyns = 0

    id_pres = []
    id_posts = []

    for i in range(len(cnt_map)):
        target = cnt_map[i]
        for j in range(len(target)):
            n = target[j]
            if n != -1:
                syn_params['tau1'].append(syn_type_params[n]['tau1'])
                syn_params['tau2'].append(syn_type_params[n]['tau2'])
                syn_params['A'].append(syn_type_params[n]['A'])
                syn_params['gmax'].append(syn_type_params[n]['gmax']+np.random.normal(loc=0, scale=g_noise))
                syn_params['e'].append(syn_type_params[n]['e'])
                syn_params['d'].append(abs(syn_type_params[n]['d']+np.random.normal(loc=0, scale=d_noise)))
                if is_syn:
                    id_pres.append(i)
                    id_posts.append(j)
                nsyns += 1

    df_syn = pd.DataFrame(syn_params)
    df_syn = df_syn.astype(float)

    return nsyns, df_syn, id_pres, id_posts


def check_f_exist(fdir, fname):
    fnames = os.listdir(fdir)
    for f in fnames:
        if fname in f:
            return True
    return False


def readOut(fname):
    times = []
    vals = []
    tspks = []
    with open(fname, "r") as fid:
        line = fid.readline()
        while line:
            tmp = line[:-1].split(',')
            times.append(float(tmp[0]))
            vals.append([float(v) for v in tmp[1:]])
            line = fid.readline()
    times = np.array(times)
    vals = np.array(vals)
    # for i in range(vals.shape[1]):
    #     tspks.append(times[vals[:, i] == 30])
    return times, vals#, tspks
