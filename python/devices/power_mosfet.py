from base_devices import ElementType, TwoPortElement

import math

import numpy as np

class power_mosfet(TwoPortElement):

    def set_params(self, KP=2e-5, VTO=1.0, L_LAMBDA=0.0, LD=0.0, L=1, W=1):
        self.KP       = KP
        self.VTO      = VTO
        self.L_LAMBDA = L_LAMBDA
        self.LD       = LD
        self.L        = L
        self.W        = W

        self.vth      = self.VTO

    def set_vgs_ref(self, vgs_ref):
        self.vgs_ref = vgs_ref

    def check_sat(self):

        # Vgs > Vth
        # Vds > Vgs - Vth
        # Vgs > 0
        # Vds > 0

        return (self.vgs > self.vth) and \
               (self.vds > (self.vgs - self.vth)) and \
               (self.vgs > 0) and (self.vds > 0)

    def check_tri(self):

        # Vgs > Vth
        # Vds < Vgs - Vth
        # Vgs > 0
        # Vds > 0

        return (self.vgs > self.vth) and \
               (self.vds <= (self.vgs - self.vth)) and \
               (self.vgs > 0) and \
               (self.vds > 0)

    def check_off(self):

        # Vgs < Vth

        return (self.vgs < self.vth)

    def sat(self):

        id = (self.KP / 2) * (self.W / (self.L - (2 * self.LD))) * \
            (self.vgs - self.vth)**2 * (1 + self.L_LAMBDA * self.vds)

        return id

    def tri(self):

        id = self.KP * (self.W / self.L - (2 * self.LD)) * \
            ( (self.vgs - self.vth - (self.vds / 2)) * (self.vds * (1 + self.L_LAMBDA * self.vds) ) )

        return id

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_degen_current(self, x, sys, t, bf, vm_sys_dy):
        return

    def get_current(self, scb):
        scb.i[self.i_x_ref] = scb.x[self.i_x_ref]

    def get_dependent_current(self, scb):

        self.vds = scb.v[self.ref]
        self.vgs = scb.v[self.vgs_ref]

        if self.check_tri():
            scb.i[self.i_d_ref] = self.tri()
        elif self.check_sat():
            scb.i[self.i_d_ref] = self.sat()
        else:
            scb.i[self.i_d_ref] = (1e-9 * self.vds)

    def get_region(self, x):

        self.vds = x[self.ref]
        self.vgs = x[self.vgs_ref]

        if self.check_tri():
            return 1
        elif self.check_sat():
            return 0
        else:
            return -1

    def get_weight(self):
        return 3

    def set_dependencies(self, ckt):

        self.i_x_ref = ckt.num_edges
        self.i_d_ref = ckt.num_edges + 1

        # Add additional row/column to Qf, Bf and degen
        ckt.qf = np.vstack([ckt.qf, np.zeros( ckt.qf.shape[0])])
        ckt.qf = np.hstack((ckt.qf, np.zeros((ckt.qf.shape[0], 1))))
        ckt.qf = np.vstack([ckt.qf, np.zeros( ckt.qf.shape[0])])
        ckt.qf = np.hstack((ckt.qf, np.zeros((ckt.qf.shape[0], 1))))

        ckt.bf = np.vstack([ckt.bf, np.zeros( ckt.bf.shape[0])])
        ckt.bf = np.hstack((ckt.bf, np.zeros((ckt.bf.shape[0], 1))))
        ckt.bf = np.vstack([ckt.bf, np.zeros( ckt.bf.shape[0])])
        ckt.bf = np.hstack((ckt.bf, np.zeros((ckt.bf.shape[0], 1))))

        ckt.qf[:, self.i_x_ref] = ckt.qf[:, self.ref]
        ckt.qf[:, self.ref]     = np.zeros(ckt.qf.shape[0])

        # Update Qf for degen
        ckt.qf[ckt.num_edges][self.i_x_ref]  = 1
        ckt.qf[ckt.num_edges][self.i_d_ref] = -1

        ckt.num_edges += 2

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        return
