import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
import csv


def red_blue_ode(t, p):
    r, b, rt, bt = p
    dp = [0, 0, 0, 0]
    lambda_a = 1.0
    lambda_t = 1.0
    p_t = 0.5
    p_c = 0.5
    k = 10.0
    flag_r_0 = 1.0 if r > 0 else 0.0
    flag_b_0 = 1.0 if b > 0 else 0.0
    flag_bt_0 = 1.0 if bt > 0 else 0.0
    flag_rt_0 = 1.0 if rt > 0 else 0.0
    flag_r_rt_0 = 1.0 if r+rt > 0 else 0.0
    flag_b_bt_0 = 1.0 if b+bt > 0 else 0.0
    # R_INDEX = 0;
    dp[0] = -flag_r_0 * lambda_a * k * r / (r + bt) * p_t * r - flag_r_0 * lambda_a * k * r / (r + bt) * p_t * rt \
            + flag_r_rt_0 * lambda_t * rt - flag_r_0 * lambda_t * r / (r + rt) * rt \
            + flag_b_0 * lambda_t * (b / (b + bt)) * bt \
            + flag_bt_0 * lambda_t * (bt / (b + bt)) * bt \
            + flag_rt_0 * lambda_a * p_c * k * rt / (rt + b) * (b + bt)

    # B_INDEX = 1;
    dp[1] = -flag_b_0 * lambda_a * k * b / (b + rt) * p_t * b \
            - flag_b_0 * lambda_a * k * b / ((b + rt)) * p_t * bt \
            + flag_b_bt_0 * lambda_t * bt \
            + flag_r_0 * lambda_t * r / (r + rt) * rt \
            - flag_b_0 * lambda_t * b / (b + bt) * bt \
            + flag_rt_0 * lambda_t * (rt / (r + rt)) * rt \
            + flag_bt_0 * lambda_a * p_c * k * (bt / (r + bt)) * (r + rt)

    # RT_INEDX = 2;
    dp[2] = -flag_rt_0 * lambda_a * p_c * k * rt / (rt + b) * (b + bt) \
            - flag_r_rt_0 * lambda_t * rt \
            - flag_rt_0 * lambda_t * (rt / (r + rt)) * rt \
            + flag_r_0 * lambda_a * k * (r / (r + bt)) * p_t * r \
            + flag_r_0 * lambda_a * k * (r / (r + bt)) * p_t * rt

    # BT_INDEX = 3;
    dp[3] = -flag_bt_0 * lambda_a * k * (bt / (r + bt)) * p_c * (r + rt) \
            - flag_b_bt_0 * lambda_t * bt \
            - flag_bt_0 * lambda_t * (bt / (b + bt)) * bt \
            + flag_b_0 * lambda_a * k * b / (b + rt) * p_t * b \
            + flag_b_0 * lambda_a * k * b / (b + rt) * p_t * bt

    return dp

def solve_ode( ):
    return solve_ivp(red_blue_ode, [0, 10], [0.97, 0.01, 0.01, 0.01], dense_output=True)


def load_data_file(file, scale=1):
    t = []
    data = []
    with open(file, 'r') as csvfile:
        plots = csv.reader(csvfile, delimiter='\t')
        for row in plots:
            t.append(float(row[0]))
            data.append(float(row[1])/(100*scale))
    return t, data


def load_simulation_data(source_dir, scale):
    t, r_data = load_data_file(source_dir+'rb_'+str(scale)+'__r_.data',scale=scale)
    _, b_data = load_data_file(source_dir+'rb_'+str(scale)+'__b_.data',scale=scale)
    _, rt_data = load_data_file(source_dir+'rb_'+str(scale)+'__rt_.data',scale=scale)
    _, bt_data = load_data_file(source_dir+'rb_'+str(scale)+'__bt_.data',scale=scale)
    return t, r_data, b_data, rt_data, bt_data


def setup_legend_and_fonts(title,file):
    plt.legend(fontsize=15,loc='best')
    plt.title(title,fontsize=20)
    plt.ylim(-0.05, 1.1)
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)
    plt.xlabel('Time units',fontsize=15)
    plt.ylabel('% Population',fontsize=15)
    plt.savefig(file)
    plt.show()


def plot_all_simulation_data(source_dir, scale):
    time, r_data, b_data, rt_data, bt_data = load_simulation_data(source_dir, scale)
    plt.plot(time, r_data,label='R')
    plt.plot(time, b_data,label='B')
    plt.plot(time, rt_data,label='RT')
    plt.plot(time, bt_data,label='BT')
    setup_legend_and_fonts('Simulation (N='+str(scale)+")",'ac_sim_'+str(scale)+'.png')



def plot_red_blue_simulation_data(source_dir, scale):
    time, r_data, b_data, rt_data, bt_data = load_simulation_data(source_dir, scale)
    red = [ r_data[i]+rt_data[i] for i in range(0,len(time))]
    blue = [ b_data[i]+bt_data[i] for i in range(0,len(time))]
    plt.plot(time, red, label='R+RT')
    plt.plot(time, blue, label='B+BT')
    setup_legend_and_fonts('Simulation (N='+str(scale)+')', 'ac_sim_rb_'+str(scale)+'.png')


def plot_red_simulation_data_with_ode(source_dir, scale):
    time, r_data, _, rt_data, _ = load_simulation_data(source_dir, scale)
    sol = solve_ode();
    t = np.linspace(0, 10, 100)
    z = sol.sol(t)
    plt.plot(time, r_data,label='R')
    plt.plot(time, rt_data,label='RT')
    plt.plot(t, z[0],label='R ODE')
    plt.plot(t, z[2],label='RT ODE')
    setup_legend_and_fonts('Fluid approximation and Simulation (N='+str(scale)+')', 'ac_sim_ode_r_rt_'+str(scale)+'.png')


def plot_blue_simulation_data_with_ode(source_dir, scale):
    time, _, b_data, _, bt_data = load_simulation_data(source_dir, scale)
    sol = solve_ode();
    t = np.linspace(0, 10, 100)
    z = sol.sol(t)
    plt.plot(time, b_data,label='B')
    plt.plot(time, bt_data,label='BT')
    plt.plot(t, z[1],label='B ODE')
    plt.plot(t, z[3],label='BT ODE')
    setup_legend_and_fonts('Fluid approximation and Simulation (N='+str(scale)+')', 'ac_sim_ode_b_bt_'+str(scale)+'.png')


def plot_red_blue_simulation_data_with_ode(source_dir, scale):
    time, r_data, b_data, rt_data, bt_data = load_simulation_data(source_dir, scale)
    sol = solve_ode();
    t = np.linspace(0, 10, 100)
    z = sol.sol(t)
    red = [ r_data[i]+rt_data[i] for i in range(0,len(time))]
    blue = [ b_data[i]+bt_data[i] for i in range(0,len(time))]
    red_ode = [ z[0][i]+z[2][i] for i in range(0,len(t)) ]
    blue_ode = [ z[1][i]+z[3][i] for i in range(0,len(t)) ]
    plt.plot(time, red, label='R+RT')
    plt.plot(time, blue, label='B+BT')
    plt.plot(t, red_ode, label='R+RT ODE')
    plt.plot(t, blue_ode, label='B+BT ODE')
    setup_legend_and_fonts('Fluid approximation and Simulation (N='+str(scale)+')', 'ac_sim_ode_rb_'+str(scale)+'.png')


if __name__=='__main__':
    dir = '../data/'
    #dir = '/Users/loreti/Desktop/DATA/'
    plot_all_simulation_data(dir,10)
    plot_red_blue_simulation_data(dir,10)
    for scale in [1, 10, 100, 1000]:
        plot_blue_simulation_data_with_ode(dir, scale)
        plot_red_simulation_data_with_ode(dir, scale)
        plot_red_blue_simulation_data_with_ode(dir, scale)

