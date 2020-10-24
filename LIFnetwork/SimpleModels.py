import numpy as np
from tqdm import tqdm

global _dt, _tmax, _nitr, _times, _seed
_dt = 0.01
_tmax = 100
_nitr = int(_tmax/_dt)
_times = np.arange(0, _tmax+_dt/2, _dt)


def set_times(tmax=100, dt=0.01):
    global _dt, _tmax, _nitr, _times
    _dt = dt
    _tmax = tmax
    _nitr = int(_tmax/_dt)
    _times = np.arange(0, _tmax+_dt/2, _dt)


def set_seed(seed=None):
    global _seed
    _seed = seed
    np.random.seed(_seed)


class CellNetwork:
    ### equations
    ## LIF model
    # tau*dV/dt = -(V - E) + R*I 
    ## synapse
    # I=g*(V-E); g=gmax*(t-onset)/tau * exp(-(t-onset-tau)/tau)
    # LIF model updated by R*I -> I = Iext - Isyns
    ### params - LIF models
    # , gmax, tau, es, pre_id, post_id, Iext, target_id
    # cell props, (v0, vth, vmax, vahp, ev, tahp, tau, r), shape=(n_cells,)
    # syn props, (gbar_syn, es, tau_r, tau_d, delay, B, pre_id, post_id), shape=(n_syns,)
    # I ext props, (Iext, target_id), shape=(n_exts, nt), shape_id=(n_exts,?) target cell number
    # g ext props(Poisson input), (g_ext, e_ext, gbar_ext, target_id), shape=(n_erxts, nt), shape_(target_id, gbar_ext)=(n_exts, ?)
    ### units
    # tau (ms), v (mV), t (ms), e (mV), R (mOhm), I (uV)
    # g (mSiemens)
    ### update, v_cells -> i_syns -> ...
    def __init__(self, **params):
        # set additional vars
        if 'std' in params.keys():
            self.std = params['std']
        else:
            self.std = 0

        self._init_cell_objs(**params)
        if 'syn' not in params.keys():
            self._init_syn_objs(**params)
        else:
            self.update = self.update_w_no_syn            
        self._init_ext_objs(**params)
        self.count = 0
    
    def _init_cell_objs(self, **params):
        # read args - cell
        # (n,)
        self.v0 = np.array(params['v0']).reshape(-1)
        self.vth = np.array(params['vth']).reshape(-1)
        self.vmax = np.array(params['vmax']).reshape(-1)
        self.vahp = np.array(params['vahp']).reshape(-1)
        self.ev = np.array(params['ev']).reshape(-1)
        self.tahp = np.array(params['tahp']).reshape(-1)
        self.tau = np.array(params['tau']).reshape(-1)
        self.r = np.array(params['r']).reshape(-1)
        # allocate variables
        self.nn = len(self.v0)
        self.vcells = np.zeros([_nitr+1, self.nn])
        self.vcells[0] = self.v0*np.ones(self.nn)
        self.fire_bool = np.zeros(self.nn, dtype=bool)
        self.tspk = -100 * np.ones(self.nn)
        # initialize spike obj
        self.spks = []
        for i in range(self.nn):
            self.spks.append([])

        # temporal array
        self.ialls = np.zeros([_nitr, self.nn])
    
    def _init_syn_objs(self, **params):
        # use double exponential function
        # g = gbar*P*(v-e), P = B*(exp(-(t-onset-d)/tau1) - exp(-(-t-onset-d)/tau2))
        # read args, gbar_syn, es, tau_r, tau_d, delay, pre_id, post_id
        self.gbar_syn = np.array(params['gbar_syn']).reshape(-1)
        self.es = np.array(params['es']).reshape(-1)
        tau_r = np.array(params['tau_r']).reshape(-1)
        tau_d = np.array(params['tau_d']).reshape(-1)
        self.tau1, self.tau2 = convert_tau_rdto12(tau_r, tau_d)
        if 'delay' in params.keys():
            self.delay = np.array(params['delay']).reshape(-1)
        else:
            self.delay = np.zeros(self.tau1.shape)
        self.B = np.array(params['B']).reshape(-1)
        # shape=(nsyns,)
        self.pre_id = tuple(params['pre_id'])
        self.post_id = np.array(params['post_id'])
        # alloc variables
        self.nsyn = len(self.gbar_syn)
        self.isyns = np.zeros([_nitr+1, self.nsyn])
        self.onset = np.ones([self.nsyn]) * (-100)
        self.g = np.zeros([_nitr+1, self.nsyn])
    
    def _init_ext_objs(self, **params):
        # external current; C*dV/dt+\sum{I_syns} = \sum{I_exts}
        self.iext = None # current clamp
        self.ghat_ext = None # Poisson input
        target_id = params['target_id']
        if 'Iext' in params.keys():
            self.iext = self._arange_ext_objs(np.array(params['Iext']), target_id)
        elif 'g_ext' in params.keys():
            # I = -ghat * (g*v- gE), get ghat(_nitr, self.nn)
            self.ghat_ext = self._arange_ext_objs(np.array(params['g_ext']), target_id) # (_nitr, self.nn)
            # get g, gE (self.nn,)
            self.gbar_ext = np.zeros(self.nn)
            self.gE_ext = np.zeros(self.nn)
            e = np.array(params['e_ext']).reshape(-1)
            gbar = params['gbar_ext'] # have same size with target_id
            for n in range(len(target_id)): # n th external input
                ids = target_id[n]
                for i, icell in enumerate(ids): # icell th cell
                    self.gbar_ext[icell] += gbar[n][i]
                    self.gE_ext[icell] += gbar[n][i] * e[n]

    def _arange_ext_objs(self, vals, target_id):
        vals_ext = np.zeros([_nitr, self.nn])
        if len(vals.shape) == 1:
            vals = vals.reshape(-1, 1)
        for n in range(len(target_id)):
            ids = target_id[n]
            for i in ids:
                vals_ext[:, i] += vals[:, n]
        return vals_ext
    
    def run(self):
        for i in tqdm(range(_nitr), ncols=120):
            self.update()

    def apply_white_noise(self):
        self.vcells[self.count+1] += np.random.normal(loc=0, scale=self.std, size=self.nn)

    def update(self):
        self.update_cells()
        self.update_syns()
        self.apply_white_noise()
        self.count += 1

    def update_cells(self):
        # copy vs
        vs = self.vcells[self.count].copy()
        if any(self.fire_bool):
            vs[self.fire_bool] = self.vahp[self.fire_bool]
        # get iall = Iext-Isyn
        # iall = self.sum_ext_current(self.iexts[self.count], self.target_id)
        iall = self.sum_ext_current(vs)
        iall -= self.sum_syn_current(self.isyns[self.count], self.post_id)
        self.ialls[self.count] = iall
        # block input current from afterhyperpolarization period
        id_ahp = _times[self.count] - self.tspk < self.tahp
        iall[id_ahp] = 0
        # use RK4 method
        dv1 = self.f_lif(vs, iall) * _dt
        dv2 = self.f_lif(vs+dv1/2, iall) * _dt
        dv3 = self.f_lif(vs+dv2/2, iall) * _dt
        dv4 = self.f_lif(vs+dv3, iall) * _dt
        # update
        self.vcells[self.count+1] = vs + 1/6*(dv1+dv2+dv3+dv4)
        # check firing
        self.fire_bool = self.vcells[self.count+1] > self.vth
        if any(self.fire_bool):
            # find synapses on pre-
            id_fire = np.where(self.fire_bool)[0]
            self.tspk[id_fire] = _times[self.count]
            self.vcells[[self.count+1], id_fire] = self.vmax[id_fire]
            # add spike time
            for i in id_fire:
                self.spks[i].append(_times[self.count+1])
            # set onset
            id_post_syns = []
            for i in id_fire:
                id_post_syns.append(np.where(self.pre_id == i)[0])
            for i in id_post_syns:
                self.onset[i] = _times[self.count]

    def f_lif(self, v, I):
        return (-(v - self.ev) + self.r*I) / self.tau

    def update_syns(self):
        # use double exponential function
        syn_bool = _times[self.count] - self.onset < self.tau1*5
        if any(syn_bool):
            g = self.gbar_syn[syn_bool] * f_syn(_times[self.count], self.onset[syn_bool], self.delay[syn_bool], self.tau1[syn_bool], self.tau2[syn_bool], self.B[syn_bool])
            self.g[self.count+1][syn_bool] = g
            self.isyns[self.count+1][syn_bool] = g * (self.vcells[[self.count+1], self.post_id[syn_bool]] - self.es[syn_bool])

    ##### one cell
    def update_w_no_syn(self):
        self.update_cells_w_no_syns()
        self.apply_white_noise()
        self.count += 1

    def update_cells_w_no_syns(self):
        # copy vs
        vs = self.vcells[self.count].copy()
        if any(self.fire_bool):
            vs[self.fire_bool] = self.vahp[self.fire_bool]
        # get iall = Iext-Isyn
        iall = self.sum_ext_current(vs)
        # block input current from afterhyperpolarization period
        id_ahp = _times[self.count] - self.tspk < self.tahp
        iall[id_ahp] = 0
        # use RK4 method
        dv1 = self.f_lif(vs, iall) * _dt
        dv2 = self.f_lif(vs+dv1/2, iall) * _dt
        dv3 = self.f_lif(vs+dv2/2, iall) * _dt
        dv4 = self.f_lif(vs+dv3, iall) * _dt
        # update
        self.vcells[self.count+1] = vs + 1/6*(dv1+dv2+dv3+dv4)
        # check firing
        self.fire_bool = self.vcells[self.count+1] > self.vth
        if any(self.fire_bool):
            # find synapses on pre-
            id_fire = np.where(self.fire_bool)[0]
            self.tspk[id_fire] = _times[self.count]
            self.vcells[[self.count+1], id_fire] = self.vmax[id_fire]
            # add spike time
            for i in id_fire:
                self.spks[i].append(_times[self.count+1])
    ##### other functions
    def sum_syn_current(self, Isyn, ids):
        current = np.zeros(self.nn)
        for i, n in enumerate(ids):
            current[n] += Isyn[i]
        return current

    def sum_ext_current(self, v):
        if self.iext is not None:
            return self.iext[self.count]
        elif self.ghat_ext is not None:
            return -self.ghat_ext[self.count] * (self.gbar_ext*v - self.gE_ext)
        else:
            return np.zeros(self.nn)

    # def get_spike(self):
    #     ids = np.where(self.fire_bool)[0]
    #     for i in ids:
    #         self.spks[i].append(_times[self.count+1])


def get_params(params_cell, params_syn, cell_types, ntk):
    # cell_types, 0~N-1
    # params_cell(N, dict), tau, R, vth, v0, vmax, e 
    # params_syn(N, dict), gmax, tau, e
    # ntk(N, N), adjacency matrix with synapse type
    ####
    # parameter setting & init
    Ncells = len(cell_types)
    read_names_cell = ['tau', 'r', 'vth', 'v0', 'vmax', 'vahp', 'ev', 'tahp']
    write_names_cell = ['tau', 'r', 'vth', 'v0', 'vmax', 'vahp', 'ev', 'tahp']
    read_names_syn = ['gbar_syn', 'tau_r', 'tau_d', 'es']
    write_names_syn = ['gbar_syn', 'tau_r', 'tau_d', 'es']
    params_ntk = dict()
    pre_ids = []
    post_ids = []
    # cell params_ntk
    for s in write_names_cell:
        params_ntk[s] = []
    # syn parms
    for s in write_names_syn:
        params_ntk[s] = []
    params_ntk['B']= []
    ####
    # parameter allocation - cell
    for n in range(Ncells):
        ctype = cell_types[n]
        for sin, sout in zip(read_names_cell, write_names_cell):
            params_ntk[sout].append(params_cell[ctype][sin])
    # Network connection - synapse
    # get normalize constant B
    nsyntype = len(params_syn)
    Bs = np.zeros(nsyntype)
    for i in range(nsyntype):
        Bs[i] = get_syn_norm(params_syn[i]['tau_r'], params_syn[i]['tau_d'])
    for i in range(Ncells):
        for j in range(Ncells):
            stype = ntk[i][j]
            if stype != -1:
                for sin, sout in zip(read_names_syn, write_names_syn):
                    params_ntk[sout].append(params_syn[stype][sin])
                params_ntk['B'].append(Bs[stype])
                pre_ids.append(i)
                post_ids.append(j)
    params_ntk['pre_id'] = pre_ids
    params_ntk['post_id'] = post_ids
    return params_ntk


def current_clamp(t0, t1, amp):
    ii = np.zeros(_nitr)
    t = np.arange(0, _tmax, _dt)
    ii[(t > t0) & (t < t1)] = amp
    return ii


def gPoisson(p, tau_r, tau_d, delay=0, t0=None, t1=None):
    tau1, tau2 = convert_tau_rdto12(tau_r, tau_d)
    B = get_syn_norm(tau_r, tau_d)
    # with double exponential function
    p_bd = p * _dt
    rands = np.random.uniform(low=0, high=1, size=_nitr)
    # rands[(_times[1:]>=t0) & (_times[1:]<=t1)] = 1
    if (t0 is None) and (t1 is None):
        isonset = np.ones(_nitr, dtype=bool)
    elif t0 is None:
        isonset = _times[1:]<=t1
    elif t1 is None:
        isonset = _times[1:]>=t0
    else:
        isonset = (_times[1:]>=t0) & (_times[1:]<=t1)
    isonset = isonset & (rands < p_bd)
    # calculate alpha function
    g = np.zeros(_nitr)
    onset = -100
    for i in range(_nitr):
        if isonset[i]:
            onset = _times[i]
        if _times[i] - onset < tau_d*5:
            g[i] = f_syn(_times[i], onset, delay, tau1, tau2, B)
    # return 1 * 
    return g


def convert_tau_rdto12(tau_r, tau_d):
    tau_r = np.array(tau_r)
    tau_d = np.array(tau_d)
    tau1 = tau_d.copy()
    tau2 = tau_r*tau_d / (tau_r + tau_d)
    return tau1, tau2


def f_syn(t, onset, d, tau1, tau2, B):
    if isinstance(onset, (np.ndarray)):
        ind = t - onset > d
        g = np.zeros(onset.shape)
        g[ind] = B[ind] * (np.exp(-(t-onset[ind]-d[ind])/tau1[ind]) - np.exp(-(t-onset[ind]-d[ind])/tau2[ind]))
        return g
    elif isinstance(t, (np.ndarray)):
        ind = t - onset > d
        g = np.zeros(t.shape)
        g[ind] = B * (np.exp(-(t[ind]-onset-d)/tau1) - np.exp(-(t[ind]-onset-d)/tau2))
        return g
    else:
        if t-onset < d:
            return 0
        else:
            return B*(np.exp(-(t-onset-d)/tau1) - np.exp(-(t-onset-d)/tau2))


def get_syn_norm(tau_r, tau_d):
    tau1, tau2 = convert_tau_rdto12(tau_r, tau_d)
    # nsyns = len(np.array(tau1))
    shape = np.array(tau1).shape
    if len(shape) == 0:
        nsyns = 1
    else:
        nsyns = shape[0]
    # for i in range(nsyns):
    val = f_syn(_times, 0, 0, tau1, tau2, 1)
    B = 1/max(val)
    return B


class LIFmodel:
    ####### description ########
    # cm*dv/dt = -g*(v - el) + I
    # -> tau*dv/dt = -(v - el) + R*I 
    # tau (ms), v (mV), t (ms), el (mV), R(mOhm), I (uV)
    #############################
    def __init__(self, tau=20, R=10, vth=-50, v0=-65, vmax=30, el=-65):
        # add constants
        self.tau = tau
        self.R = R
        self.v0 = v0
        self.vth = vth
        self.vmax = vmax
        self.el = el
        self._init_obj()

    def _init_obj(self):
        self.v = self.v0
        self.flag_fire = False
        self.Iext = None
        self.n = 0 # count

    def update_all(self, Is):
        self.v = self.v0
        self.vs = np.zeros(_nitr+1)
        self.vs[0] = self.v
        for i in tqdm(range(_nitr)):
            self.update(Is[i])
            self.vs[i+1] = self.v

    def _f(self, v, I):
        return (-(v - self.el) + self.R*I) / self.tau

    def update(self, I):
        # use RK4 method
        if self.flag_fire:
            self.flag_fire = False
            self.v = self.v0
        dv1 = self._f(self.v, I)*_dt
        dv2 = self._f(self.v+dv1/2, I)*_dt
        dv3 = self._f(self.v+dv2/2, I)*_dt
        dv4 = self._f(self.v+dv3, I)*_dt
        # update
        self.v += 1/6*(dv1+dv2+dv3+dv4)
        if self.v > self.vth:
            # vnew = vmax
            self.v = self.vmax
            self.flag_fire = True
    
    def add_current(self, ii):
        # ii, (n, ) list
        if self.Iext is None:
            self.Iext = ii
        else:
            self.Iext = np.vstack((self.Iext, ii))

    def reset(self):
        self._init_obj()


if __name__=='__main__':

    import matplotlib.pyplot as plt
    # cell props, (v0, vth, vmax, ev, tau, R), shape=(n_cells,)
    # syn props, (gmax, tau, es, pre_id, post_id), shape=(n_syns,), shape_id=(n_cells, ?)
    # I ext props, (Iext, target_id), shape=(n_exts, nt), shape_id=(n_cells, ?)
    set_times(tmax=100, dt=0.05)
    # cell props
    v0 = [-65, -65]
    vth = [-50, -50]
    vmax = [30, 30]
    ev = [-65, -65]
    tau_cell = [20, 20]
    r = [120, 120]
    # syn props
    gmax = 0.01
    # tau_syn = 1
    tau_r = 0.2
    tau_d = 5
    es = 0
    B = get_syn_norm(tau_r, tau_d)
    pre_id = [0]
    post_id = [1]
    # I ext
    # Iext = np.zeros(int(tmax//dt))
    # Iext[int(20*dt):int(80*dt)] = 20
    # Iext = 2 * np.ones(int(tmax/dt)).reshape(-1, 1)
    Iext = current_clamp(20, 80, 1)
    target_id = [[0]]
    
    
    ntk = CellNetwork(v0=v0, vth=vth, vmax=vmax, ev=ev, tau=tau_cell, r=r,
                    gbar_syn=gmax, tau_r=tau_r, tau_d=tau_d, es=es, B=B, pre_id=pre_id, post_id=post_id,
                    Iext=Iext, target_id=target_id)
    ntk.run()

    # plot
    plt.figure(dpi=200)
    plt.subplot(3,1,1)
    plt.plot(_times, ntk.vcells[:, 0], 'k', lw=0.5)
    plt.ylabel('voltage')
    plt.subplot(3,1,2)
    # plt.plot(_times, ntk.isyns, 'r', lw=0.5)
    plt.plot(_times, ntk.g, 'r', lw=0.5)
    plt.ylabel('current')
    plt.subplot(3,1,3)
    plt.plot(_times, ntk.vcells[:, 1], 'k', lw=0.5)
    plt.ylabel('voltage')
    plt.xlabel('time')
    plt.tight_layout()
    plt.show()
