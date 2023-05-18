import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from tqdm import tqdm
import torch

def def_pars():

    #PVAR
    alpha_T  = 7.227552687697e-3
    alpha_C = 6.467024536808274e-6
    KG = 57.158526714327490
    hON_min = 2.020162670696560e-7
    hON_max = 2.001991006981458e-5
    KL = 1.985087802239813e3
    nL = 1.354848683543133
    nG = 1.555710381562409
    hC = 1.538847970212344
    gamma_mol2fluo = 0.254897149724905
    L0 = 196.3930019326237

    pvar = [alpha_T, alpha_C, KG, hON_min, hON_max, KL, nL, nG, hC, gamma_mol2fluo, L0]

    #PFIX
    Np = 10
    KD = 1300
    KC = 0.33333333333333333
    kappa = 90
    nr = 12221
    nT = 597
    nC = 219
    g0 = 0.098685
    PhiR_max = 0.5470
    PhiR0 = 0.0660
    rho_cell = 2e9
    nu = 0.1921
    AE = 10.5e3

    pfix = [Np, KD, KC, kappa, nr, nT, nC, g0, PhiR_max, PhiR0, rho_cell, nu, AE]

    return pvar, pfix

def get_init_cond(tmax = 200, L0=None):
    if L0 is None:
        pvar, _ = def_pars()
        L0 = pvar[-1]
    initial_sol = solve_ivp(ode_fun, [0, tmax*60], [100, 100, 100, 0.5], args=(L0,1.58369475/60), method='BDF', t_eval=[0, tmax*60], atol=1e-10, rtol=1e-10)

    return initial_sol.y[:, -1]

def lambda_p_fun(C, Tt):
    pvar, pfix = def_pars()
    _, _, _, _, _, _, _, _, hC, _, _ = pvar
    _, KD, KC, kappa, _, nT, nC, g0, PhiR_max, PhiR0, rho_cell, nu, AE = pfix

    PhiS = (nT * Tt + nC * C) / rho_cell
    lambda_p = (PhiR_max - PhiR0 - PhiS) * (nu * g0) / (g0 + nu * (1 + (AE/KD)/(1+((C/(kappa * KC)) ** hC))))

    return lambda_p

def compute_groups(L, Td, lambda_p):
    pvar, pfix = def_pars()
    _, _, KG, hON_min, hON_max, KL, nL, nG, _, _, _ = pvar
    Np, _, _, _, nr, _, _, g0, _, _, rho_cell, _, _ = pfix

    hON = hON_min + (hON_max - hON_min) * (L ** nL) / (L ** nL + KL ** nL)
    gON = Np * Td ** nG / (Td ** nG + KG ** nG)

    # Unbound ribosomes
    r_u = rho_cell * lambda_p / (nr*g0)

    return hON, gON, r_u

def ode_fun(t,x, L, lambda_c=0.0263949125):
    # Variables     
    Tt, Td, C, phi_p = x

    # Parameters
    pvar, pfix = def_pars()
    alpha_T, alpha_C, _, _, _, _, _, _, _, _, _ = pvar
    _, _, _, _, _, _, _, _, _, _, _, nu, _ = pfix

    # Positive variables
    Tt = np.maximum(Tt, np.zeros_like(Tt))
    Td = np.maximum(Td, np.zeros_like(Td))
    C = np.maximum(C, np.zeros_like(C))
    phi_p = np.maximum(phi_p, np.zeros_like(phi_p))

    # Growth Rate
    lambda_p = lambda_p_fun(C, Tt)

    # Useful functions
    hON, gON, r_u = compute_groups(L, Td, lambda_p)

    dTtdt = (alpha_T * lambda_p * r_u / nu - lambda_p * Tt).squeeze() # dTt/dt
    dTddt = (hON * (Tt - 2 * (Td + gON)) ** 2 - lambda_p * (Td+gON)).squeeze() # dTd/dt
    dCdt = (alpha_C * (r_u / nu) * gON - lambda_p * C).squeeze() # dC/dt
    dphidt = ((lambda_p - lambda_c) * (1 - phi_p) * phi_p).squeeze() # dphi_p/dt

    ddt = np.stack([dTtdt, dTddt, dCdt, dphidt], axis=0).squeeze()

    return ddt

def integrate_OL(phi_p0, L_arr, t_arr):
    initial_cond = get_init_cond()
    initial_cond[-1] = phi_p0

    if np.isscalar(L_arr):
        L_arr = np.array([L_arr])
    if np.isscalar(t_arr):
        t_arr = np.array([t_arr])

    if len(L_arr) != len(t_arr):
        raise ValueError('L_arr and t_arr must have the same length')

    output_y = []
    output_t = []
    t_init = 0
    for L, tmax in zip(L_arr, t_arr):
        sol = solve_ivp(ode_fun, [t_init, t_init+tmax*60], initial_cond, args=(L,), method='BDF', t_eval=np.linspace(t_init, t_init+tmax*60, 1000), atol=1e-10, rtol=1e-10)
        output_y.append(sol.y)
        output_t.append(sol.t)
        initial_cond = sol.y[:,-1]
        t_init = t_init+tmax*60
    
    output_y = np.hstack(output_y)
    output_t = np.hstack(output_t) / 60

    return output_t, output_y

def integrate_CL(phi_p0, L0, pid_par, sp_arr, t_arr, sampling_time=0.5, time_resolution=0.1):
    Kp, Ki, Kd, Kbc = pid_par

    initial_cond = get_init_cond(L0=L0)
    initial_cond[-1] = phi_p0

    if np.isscalar(sp_arr):
        sp_arr = np.array([sp_arr])
    if np.isscalar(t_arr):
        t_arr = np.array([t_arr])

    def sp_fun(t):
        if t >= np.sum(t_arr):
            sp = sp_arr[-1]
        else:
            sp = sp_arr[np.argmax(np.cumsum(t_arr) > t)]
        return sp
    
    if len(sp_arr) != len(t_arr):
        raise ValueError('sp_arr and t_arr must have the same length')

    output_y = []
    output_t = []
    output_L = []
    output_sp = []

    t_first = 0
    t_last = np.sum(t_arr)
    t_next = t_first + sampling_time

    # Initial Controller Action
    # Compute error
    sp = sp_fun(0)
    error = sp - phi_p0
    error_old = error
    error_accum = 0

    # Compute integral
    error_accum = error_accum + error * (sampling_time * 60)
    
    # Compute output
    u_compute = round(Kp * error + Kd * (error - error_old)/(sampling_time * 60) + Ki * error_accum)

    # Saturation
    u = np.clip(u_compute, 0, 800)
    L = u #* 800

    # Back-calculation (anti-windup)
    bc = u - u_compute
    error_accum = error_accum + bc * Kbc

    error_old = error
    output_y.append(initial_cond.reshape(-1,1))
    output_t.append(np.array([0]))
    output_L.append(np.array([L]))
    output_sp.append(np.array(sp_fun(0)))

    while t_next <= t_last:
        # print(t_next, L)
        t_eval = np.arange(t_first*60, t_next*60, time_resolution)
        sol = solve_ivp(ode_fun, [t_first*60, t_next*60], initial_cond, args=(L,), method='BDF', t_eval=t_eval, atol=1e-10, rtol=1e-10)
        output_y.append(sol.y[:,1:])
        output_t.append(sol.t[1:])
        output_L.append(np.ones(len(sol.t[1:])) * L)
        output_sp.append(np.ones(len(sol.t[1:])) * sp_fun(t_next))
        initial_cond = sol.y[:,-1]
        t_first = t_next
        t_next = t_next + sampling_time

        # Compute new L
        phi_p = sol.y[3,-1]
        sp = sp_fun(t_next)
        error = sp - phi_p

        # If setpoint has changed
        # if sp != sp_old:
        #     error_old = error
        #     error_accum = 0
        #     sp_old = sp
        #     u0 = u

        # Compute integral
        error_accum = error_accum + error * (sampling_time * 60)
        
        # Compute output
        u_compute = np.round(Kp * error + Kd * (error - error_old)/(sampling_time * 60) + Ki * error_accum)

        # Saturation
        u = np.clip(u_compute, 0, 800)
        L = u #* 800

        # Back-calculation (anti-windup)
        bc = u - u_compute
        error_accum = error_accum + bc * Kbc

        error_old = error

    output_y = np.hstack(output_y)
    output_t = np.hstack(output_t) / 60
    output_L = np.hstack(output_L)
    output_sp = np.hstack(output_sp)

    return output_t, output_y, output_L, output_sp


def optimize_gains(x0=None):
    def fmin_fun(x, info):
        Kp_exp, Ki_exp, Kd_exp, Kbc_exp = x
        Kp = 10**Kp_exp
        Ki = 10**Ki_exp
        Kd = 10**Kd_exp
        Kbc = 10**Kbc_exp

        u0 = 0.5
        t_arr = np.array([50, 50, 50])
        sp_arr = np.array([0.2, 0.7, 0.2])
        t, y, L, sp = integrate_CL(u0, 196.3930, (Kp, Ki, Kd, Kbc), sp_arr, t_arr, sampling_time=0.5)
        output = np.average((sp - y[3,:])**2)

        if info['Nfeval']%1 == 0:
            print('{0:4d}   {1: 3.6f}   {2: 3.6f}   {3: 3.6f}   {4: 3.6f}   {5: 3.6e}'.format(info['Nfeval'], Kp, Ki, Kd, Kbc, output))
        info['Nfeval'] += 1
        return output

    print('{0:4s}   {1:9s}   {2:9s}   {3:9s}   {4:9s}   {5:9s}'.format('Iter', ' Kp', ' Ki', ' Kd', 'Kbc', 'f(X)')  )
    if x0 is None:
        x0 = np.array([1, -0.3, 1, -1])
    
    bounds = [(-2, 4), (-2, 4), (-2, 4), (-2, 4)]
    res = minimize(fmin_fun, x0, args=({'Nfeval':0},), method='Nelder-Mead', options={'maxiter': 10000}, bounds=bounds)
    
    return res

def pid_pars(Kp=9.96358683,  Ki=0.3627385,  Kd=11.52560098,  Kbc=0.15180135):
    return Kp, Ki, Kd, Kbc

def datagen(n_traj, sp_per_traj, pid_parameters=None, tmax=50, time_resolution=0.1):
    t_arr = np.array([tmax] * sp_per_traj)
    y_out = []
    L_out = []
    sp_out = []
    t_out = []

    if pid_parameters is None:
        pid_parameters = pid_pars()

    for i in tqdm(range(n_traj)):
        sp_arr = np.random.rand(sp_per_traj) * 0.6 + 0.2
        L0 = np.random.rand() * 800
        u0 = np.random.rand() * 0.6 + 0.2
        t, y, L, sp = integrate_CL(u0, L0, pid_parameters, sp_arr, t_arr, sampling_time=0.5, time_resolution=time_resolution*60)
        y_out.append(y)
        L_out.append(L)
        sp_out.append(sp)
        t_out.append(t)

    return np.array(t_out), np.swapaxes(np.array(y_out),1,2), np.array(L_out), np.array(sp_out)

def time_embed(phi, L, n_embed, n_gap=1):
    phi_out_t0 = []
    phi_out_t1 = []
    L_out_t0 = []
    for j in range(phi.shape[0]):
        for i in range(n_embed*n_gap, phi.shape[1]-1):
            phi_out_t0.append(phi[j,i-n_embed*n_gap:i:n_gap])
            phi_out_t1.append(phi[j,i-n_embed*n_gap+1:i+1:n_gap])
            L_out_t0.append(L[j,i-n_embed:i])
    return np.array(phi_out_t0), np.array(L_out_t0), np.array(phi_out_t1)


class Dataset(torch.utils.data.Dataset):
    def __init__(self, phi, L, phi_t1):
        self.phi = phi
        self.L = L
        self.phi_t1 = phi_t1

    def __len__(self):
        return self.phi.shape[0]

    def __getitem__(self, idx):
        phi = torch.from_numpy(self.phi[idx]).float()
        L = torch.from_numpy(self.L[idx]).float()
        phi_t1 = torch.from_numpy(self.phi_t1[idx]).float()

        return phi, L, phi_t1