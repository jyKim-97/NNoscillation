
import json
import numpy as np


def CreateERnetwork(pre_types, post_types, probs):

    ntk = []
    n = len(post_types)
    for i, c in enumerate(pre_types):
        p = np.random.uniform(low=0, high=1, size=n)
        ntk.append([])
        for j in range(n):
            if i != j:
                if p[j] < probs[c][post_types[j]]:
                    ntk[i].append(j)
                    
    return ntk


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


def convert_syn_param(syn_params):
    """ tau_r, tau_d -> tau1, tau2, A (normalization constant) """
    tau1, tau2 = convert_tau_rdto12(syn_params['tau_r'], syn_params['tau_d'])
    A = get_syn_norm(syn_params['tau_r'], syn_params['tau_d'])
    syn_params['tau1'] = tau1
    syn_params['tau2'] = tau2
    syn_params['A'] = A


def printSimulInfo(fname, ntk_info, ext_info, d_noise_std=0):
    """
    ntk_info
    "cell_params", "syn_params", "cell_types", "syn_types", "adj_list", "cell_type_names", "syn_type_names"
    ext_info
    "input_params", "ext_types", "syn_types", "adj_list", "t_infos", "fs", "ext_type_names", "syn_type_names"
    """
    data = dict()

    """ Print Summary Info """
    data["NumofCells"] = len(ntk_info["cell_types"])
    data["NumofSyns"] = 0
    data["NumofPinputs"] = 0
    data["NumOfPinputSynapses"] = 0

    """ Print Cell Info """
    N = len(ntk_info["cell_types"])
    cell_infos = []
    for i in range(N):
        info = dict()

        c = ntk_info["cell_types"][i]
        params = ntk_info["cell_params"][c]
        
        info["type"]     = ntk_info["cell_type_names"][c]
        info["tau"]      = params["tau"]
        info["r"]        = params["r"]
        info["e"]        = params["e"]
        info["v"]        = np.random.normal(loc=-65, scale=1)
        info["vth"]      = params["vth"]
        info["vahp"]     = params["vahp"]
        info["vmax"]     = params["vmax"]
        info["t_refrac"] = params["t_refrac"]

        cell_infos.append(info)
    data["Cells"] = cell_infos

    """ Print Synapse Info """
    # convert params
    for param in ntk_info["syn_params"]:
        convert_syn_param(param)

    syn_infos = []
    for i in range(N):
        c_pre = ntk_info["cell_types"][i]
        for j in ntk_info["adj_list"][i]:
            c_post = ntk_info["cell_types"][j]

            tp = ntk_info["syn_types"][c_pre][c_post]
            params = ntk_info["syn_params"][tp]

            info = dict()
            
            info["type"]    = ntk_info["syn_type_names"][tp]
            info["tau1"]    = params["tau1"]
            info["tau2"]    = params["tau2"]
            info["A"]       = params["A"]
            info["gmax"]    = params["gmax"]
            info["e"]       = params["e"]
            info["d"]       = params["d"] + np.random.normal(loc=0, scale=d_noise_std)
            info["id_pre"]  = i
            info["id_post"] = j

            syn_infos.append(info)

    data["NumofSyns"] = len(syn_infos)
    data["Synapses"] = syn_infos

    """ Print Stim Info - Poisson input """
    Next = len(ext_info["ext_types"])
    input_infos = []
    for i in range(Next):

        info = dict()
        tp = ext_info["ext_types"][i]

        info["type"]   = ext_info["ext_type_names"][tp]
        info["tstart"] = t_infos[i][0]
        info["tend"]   = t_infos[i][1]
        info["f"]      = f[i]

        input_infos.append(info)

    data["NumofPinputs"] = len(input_infos)
    data["Pinputs"] = input_infos

    """ Print Stim Synapse info"""
    # convert params
    for param in ext_info["input_params"]:
        convert_syn_param(param)

    input_syn_infos = []
    for i in range(Next):
        c_pre = ext_info["ext_types"][i]
        for j in ext_info["adj_list"][i]:
            c_post = ntk_info["cell_types"][j]

            tp = ext_info["syn_types"][c_pre][c_post]
            params = ext_info["input_params"][tp]

            info = dict()

            info["type"]    = ext_info["syn_type_names"][tp]
            info["tau1"]    = params["tau1"]
            info["tau2"]    = params["tau2"]
            info["A"]       = params["A"]
            info["gmax"]    = params["gmax"]
            info["e"]       = params["e"]
            info["d"]       = params["d"]
            info["id_pre"]  = i
            info["id_post"] = j

            input_syn_infos.append(info)
    
    data["NumOfPinputSynapses"] = len(input_syn_infos)
    data["PinputSynapses"]      = input_syn_infos

    json_data = json.dumps(data)
    
    with open(fname, "w", encoding="utf-8") as fid:
        json.dump(data, fid, ensure_ascii=False, indent="\t")


def printNetworkInfo(fname, adj_list, cell_types, cell_type_names):
    n = len(adj_list)
    with open(fname, "w") as fid:
        fid.write("DL n=%d\n"%(n))
        fid.write("format = edgelist1\n")
        fid.write("labels:\n")
        for c in cell_types[:-1]:
            fid.write("%s, "%(cell_type_names[c]))
        fid.write("%s\n"%(cell_type_names[cell_types[-1]]))
        fid.write("data:\n")
        for i in range(n):
            for j in adj_list[i]:
                fid.write("%d %d\n"%(i+1, j+1))


if __name__ == '__main__':

    np.random.seed(1000)

    # set PN cell type & cells parameter
    params_pn = {'tau':20, 'r':100, 'e':-65, 'vth':-40, 'vahp':-80, 'v0':-65, 'vmax':30, 't_refrac':10}
    params_pv = {'tau':5, 'r':100, 'e':-65, 'vth':-50, 'vahp':-80, 'v0':-65, 'vmax':30, 't_refrac':0}

    # set synapse type & parameter
    syn_pn2pn = {'gmax':2e-3, 'tau_r':0.3, 'tau_d':1, 'e':0, 'd':1}
    syn_pn2pv = {'gmax':1e-3, 'tau_r':0.3, 'tau_d':1, 'e':0, 'd':1}
    syn_pv2pn = {'gmax':0, 'tau_r':1, 'tau_d':3, 'e':-80, 'd':1}
    syn_pv2pv = {'gmax':1e-3, 'tau_r':1, 'tau_d':3, 'e':-80, 'd':1}

    # connection probability
    probs = [[0.2, 0.2],
             [0.3, 0.3]]
    # PN, IN
    ncells = [10, 2]
    nall = ncells[0] + ncells[1]
    cell_types = []
    for i in range(nall):
        if i < ncells[0]:
            cell_types.append(0)
        else:
            cell_types.append(1)

    # create network
    syn_types = [[0, 1], [2, 3]]
    ntk = CreateERnetwork(cell_types, cell_types, probs)
    ntk_info = {"cell_params": [params_pn, params_pv],
                "syn_params": [syn_pn2pn, syn_pn2pv, syn_pv2pn, syn_pv2pv],
                "cell_types": cell_types,
                "syn_types": syn_types,
                "adj_list": ntk,
                "cell_type_names": ["e", "i"],
                "syn_type_names": ["syn_e2e", "syn_e2i", "syn_i2e", "syn_i2i"]
                }

    # create external inputs
    # ext_pn2pn = {'gmax':3e-4, 'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':0}
    # ext_pn2pv = {'gmax':3e-4, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':0}
    ext_pn2pn = {'gmax':0.5, 'tau_r':0.2, 'tau_d':4.6, 'e':0, 'd':0}
    ext_pn2pv = {'gmax':0.5, 'tau_r':0.1, 'tau_d':1.3, 'e':0, 'd':0}

    next = [5]
    nall = next[0]
    
    probs_ext = [[0.2, 0.1]]

    # create ext network
    ext_types = [0 for i in range(nall)]
    ntk_ext = CreateERnetwork(ext_types, cell_types, probs_ext)

    # set Poisson input frequencies
    f = [60 for i in range(nall)]

    # set t0 & t1
    t_infos = np.zeros([nall, 2])
    for i in range(nall):
        t_infos[i][0] = 1000
        t_infos[i][1] = 9000
    
    ext_info = {"input_params": [ext_pn2pn, ext_pn2pv],
                "ext_types": ext_types,
                "syn_types": [[0, 1]],
                "adj_list": ntk_ext,
                "t_infos": t_infos,
                "fs": f,
                "ext_type_names": ["ext_e"],
                "syn_type_names": ["syn_e2e", "syn_e2i"]
                }

    printSimulInfo("test.json", ntk_info, ext_info)

    printNetworkInfo("test.dl", ntk, cell_types, ["e", "i"])


    
    

    



