import sys

sys.path.append("../devices")
from power_mosfet import power_mosfet

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

import numpy as np

import matplotlib.pyplot as plt

def main():

    Vd_sweep = np.linspace(-5, 0.5, 1000)
    Id       = np.zeros(len(Vd_sweep))

    vgs_sweep  = np.arange(0.8, 5.0,  0.05)
    vds_sweep  = np.arange(0.0, 12.0, 0.01)

    i_ds_n1    = np.zeros([len(vgs_sweep), len(vds_sweep)])

    mosfet_dut = power_mosfet()

    mosfet_dut.set_params(KP=67.9211, VTO=2.08819, L_LAMBDA=0.0038193, L=100e-6, W=100e-6)
    mosfet_dut.set_ref(ref=0)
    mosfet_dut.set_vgs_ref(vgs_ref=1)
    mosfet_dut.i_d_ref = 0

    scb = Scoreboard(num_edges=2, num_sys_vars=0, degen_mtrx=0)

    fig, ax1 = plt.subplots()

    for vgs_idx, vgs in enumerate(vgs_sweep):
        for vds_idx, vds in enumerate(vds_sweep):

            scb.v[0] = vds
            scb.v[1] = vgs
            mosfet_dut.get_dependent_current(scb)
            i_ds_n1[vgs_idx][vds_idx] = scb.i[0]
            print(i_ds_n1[vgs_idx][vds_idx])

        l = "Vgs = {:.2f}".format(vgs)
        plt.plot(vds_sweep, i_ds_n1[vgs_idx], label=l)

    plt.show()

if __name__ == "__main__":

    main()
