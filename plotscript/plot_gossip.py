import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import csv


def unicast_ode(t, p):
    dp = [0, 0, 0, 0, 0, 0]
    diffrate = 10.0
    passrate = 1.0
    AI = 0
    AU = 1
    AS = 2
    PI = 3
    PU = 4
    PS = 5
    # Ai
    dp[AI] = -diffrate * p[AI] + passrate * p[PI]
    # Au
    dp[AU] = -diffrate * p[AU] + passrate * p[PU]
    # As
    dp[AS] = -diffrate * p[AS] + passrate * p[PS]
    # Pi
    dp[PI] = +diffrate * p[AI] - passrate * p[PI]
    # Pu
    dp[PU] = +diffrate * p[AU] * (p[PU] + p[PS]) / (p[PI] + p[PU] + p[PS])-passrate * p[PU]-diffrate * p[AI] * p[PU] / (p[PI] + p[PU] + p[PS])
    # Pu
    dp[PS] = +diffrate * p[AI] * p[PU] / (p[PI] + p[PU] + p[PS]) +diffrate * (p[AS] + p[AU]) * (p[PI]) / (p[PI] + p[PU] + p[PS]) +diffrate * p[AS] * (p[PU] + p[PS]) / (p[PI] + p[PU] + p[PS]) -passrate * p[PS]
    return dp


def broadcast_ode(t, p):
    dp = [0, 0, 0, 0, 0, 0]
    diffrate = 10.0
    passrate = 1.0
    prob = 1.0
    k = 10
    AI = 0
    AU = 1
    AS = 2
    PI = 3
    PU = 4
    PS = 5

    # Ai
    dp[AI] = -diffrate * p[AI] + passrate * p[PI]

    # Au
    dp[AU] = -diffrate * p[AU] + passrate * p[PU]

    # Ax
    dp[AS] = -diffrate * p[AS] + passrate * p[PS]

    # Pi
    dp[PI] = +diffrate * p[AI] -passrate * p[PI] +diffrate * p[AI] * k * p[PU] / (p[PI] + p[PU] + p[PS]) +diffrate * p[AI] * k * p[PS] / (p[PI] + p[PU] + p[PS])\
             -diffrate * (p[AU] + p[AS]) * k * p[PI] / (p[PI] + p[PU] + p[PS]);

    # Pu
    dp[PS] = +diffrate * p[AS] -passrate * p[PS] -diffrate * p[AI] * k * p[PS] / (p[PI] + p[PU] + p[PS]) +diffrate * (p[AU] + p[AS]) * k * p[PI] / (p[PI] + p[PU] + p[PS])

    dp[PU] = +diffrate * p[AU] -passrate * p[PU] -diffrate * p[AI] * k * p[PU] / (p[PI] + p[PU] + p[PS])
    return dp


def load_data_file(file, scale=1):
    t = []
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t.append(float(row[0]))
            data.append(float(row[1])/(100*scale))
    return t, data


def solve_ode( ode_system ):
    return solve_ivp(ode_system, [0, 20],[0.00, 0.00, 0.00, 0.20, 0.80, 0.00], dense_output=True)


def setup_legend_and_fonts(title,file):
    plt.legend(fontsize=15,loc='best')
    plt.title(title,fontsize=20)
    plt.ylim(-0.05, 1.1)
    plt.xlim(-0.05, 20)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    plt.savefig(file)
    plt.show()


def load_simulation_data(source_dir, prefix, scale):
    t, ai_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__ai_.data',scale=scale)
    _, as_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__as_.data',scale=scale)
    _, pi_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__pi_.data',scale=scale)
    _, ps_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__ps_.data',scale=scale)
    _, au_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__au_.data',scale=scale)
    _, pu_data = load_data_file(source_dir+prefix+'_'+str(scale)+'__pu_.data',scale=scale)
    return t, ai_data, au_data, pi_data, pu_data, as_data, ps_data


def plot_all_simulation_data(source_dir, prefix, scale):
    time, ai_data, au_data, pi_data, pu_data, as_data, ps_data = load_simulation_data(source_dir, prefix, scale)
    plt.plot(time, ai_data, label='AI')
    plt.plot(time, au_data, label='AU')
    plt.plot(time, pi_data, label='PI')
    plt.plot(time, pu_data, label='PU')
    plt.plot(time, as_data, label='AS')
    plt.plot(time, ps_data, label='PS')
    setup_legend_and_fonts()


def plot_uninformed_simulation_data(source_dir, prefix, scale):
    time, _, au_data, _, pu_data, _, _ = load_simulation_data(source_dir, prefix, scale)
    data = [ au_data[i] + pu_data[i] for i in range(0,len(time))]
    plt.plot(time, data,label='AU+PU')
    setup_legend_and_fonts()


def plot_unicast_uninformed_simulation_data(source_dir, scale):
    plot_uninformed_simulation_data(source_dir,'u',scale)


def plot_broadcast_uninformed_simulation_data(source_dir, scale):
    plot_uninformed_simulation_data(source_dir,'bc',scale)


def plot_unicast_uninformed_simulation_data_with_ode(source_dir, scale):
    time, _, au_data, _, pu_data, _, _ = load_simulation_data(source_dir, 'u', scale)
    data = [ au_data[i] + pu_data[i] for i in range(0,len(time))]
    plt.plot(time, data,label='AU+PU')
    sol = solve_ivp( unicast_ode, [0,20], [0.00, 0.00, 0.00, 0.20, 0.80, 0.00], dense_output=True)
    t = np.linspace(0, 20, 100)
    z = sol.sol(t)
    fluiddata = [ z[1][i]+z[4][i] for i in range(0,len(t)) ]
    plt.plot(t, fluiddata, label='AU+PU ODE')
    setup_legend_and_fonts('Fraction of Informed Agents (N='+str(scale)+')', 'gossip_u_uninformed_'+str(scale)+".png")


def plot_broadcast_uninformed_simulation_data_with_ode(source_dir, scale):
    time, _, au_data, _, pu_data, _, _ = load_simulation_data(source_dir, 'bc', scale)
    data = [ au_data[i] + pu_data[i] for i in range(0,len(time))]
    plt.plot(time, data,label='AU+PU')
    sol = solve_ivp( broadcast_ode , [0,20], [0.00, 0.00, 0.00, 0.20, 0.80, 0.00], dense_output=True)
    t = np.linspace(0, 20, 100)
    z = sol.sol(t)
    fluiddata = [ z[1][i]+z[4][i] for i in range(0,len(t)) ]
    plt.plot(t, fluiddata, label='AU+PU ODE')
    setup_legend_and_fonts('Fraction of Informed Agents (N='+str(scale)+')', 'gossip_b_uninformed_'+str(scale)+".png")


def plot_unicast_informed_simulation_data_with_ode(source_dir, scale):
    time, ai_data, _, pi_data, _, ax_data, px_data = load_simulation_data(source_dir, 'u', scale)
    data = [ ai_data[i] + pi_data[i] + ax_data[i] + px_data[i] for i in range(0,len(time))]
    plt.plot(time, data,label='AI+PI+AS+PS')
    sol = solve_ivp( unicast_ode , [0,20], [0.00, 0.00, 0.00, 0.20, 0.80, 0.00], dense_output=True)
    t = np.linspace(0, 20, 100)
    z = sol.sol(t)
    fluiddata = [ z[0][i]+z[2][i]+z[3][i]+z[5][i] for i in range(0,len(t)) ]
    plt.plot(t, fluiddata, label='AI+PI+AS+PS ODE')
    setup_legend_and_fonts('Fraction of Informed Agents (N='+str(scale)+')', 'gossip_u_informed_'+str(scale)+".png")


def plot_broadcast_informed_simulation_data_with_ode(source_dir, scale):
    time, ai_data, _, pi_data, _, ax_data, px_data = load_simulation_data(source_dir, 'bc', scale)
    data = [ ai_data[i] + pi_data[i] + ax_data[i] + px_data[i] for i in range(0,len(time))]
    plt.plot(time, data,label='AI+PI+AS+PS')
    sol = solve_ivp( broadcast_ode , [0,20], [0.00, 0.00, 0.00, 0.20, 0.80, 0.00], dense_output=True)
    t = np.linspace(0, 20, 100)
    z = sol.sol(t)
    fluiddata = [ z[0][i]+z[2][i]+z[3][i]+z[5][i] for i in range(0,len(t)) ]
    plt.plot(t, fluiddata, label='AI+PI+AS+PS ODE')
    setup_legend_and_fonts('Fraction of Informed Agents (N='+str(scale)+')', 'gossip_b_informed_'+str(scale)+".png")


if __name__=='__main__':
    dir = '../data/'
    for scale in [100]:
        plot_broadcast_informed_simulation_data_with_ode(dir,scale)
        plot_unicast_informed_simulation_data_with_ode(dir,scale)
        plot_broadcast_uninformed_simulation_data_with_ode(dir,scale)
        plot_unicast_uninformed_simulation_data_with_ode(dir,scale)


