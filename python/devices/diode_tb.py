import sys

sys.path.append("../devices")
from diode import diode_id

sys.path.append("../eqn")
from eqn_syn import Circuit, Scoreboard

import numpy as np

import matplotlib.pyplot as plt

def main():

    Vd_sweep = np.linspace(-5, 0.5, 1000)
    Id       = np.zeros(len(Vd_sweep))

    diode_dut = diode_id()

    diode_dut.set_params(IS=75e-12, N=1.0, EG=0.7, CJO=26e-12, M=0.5, IBV=5e-6, BV=400, TT=4.2e-6)
    diode_dut.set_ref(ref=0)
    diode_dut.i_d_ref = 0

    scb = Scoreboard(num_edges=1, num_sys_vars=0, degen_mtrx=0)

    for idx, Vd in enumerate(Vd_sweep):

        scb.v[0] = Vd
        diode_dut.get_dependent_current(scb)
        Id[idx] = scb.i[0]

    fig, ax1 = plt.subplots()

    plt.grid()
    ax1.plot(Vd_sweep, Id)

    plt.show()

if __name__ == "__main__":

    main()
