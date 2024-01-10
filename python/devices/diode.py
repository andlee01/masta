from base_devices import ElementType, TwoPortElement

import math

import numpy as np

class diode_id(TwoPortElement):

    def __init__(self):
        TwoPortElement.__init__(self)

        self.type = ElementType.current_src

    def set_params(self,     \
                   IS=50e-15,\
                   RS=0.0,   \
                   N=1.0,    \
                   TT=0,     \
                   CJO=0,    \
                   VJ=1.0,   \
                   M=0.5,    \
                   EG=1.11,  \
                   XTI=3.0,  \
                   KF=0,     \
                   AF=1.0,   \
                   FC=0.5,   \
                   BV=100.0, \
                   IBV=1e-3, \
                   T=273):

        self.IS=IS
        self.RS=RS
        self.N=N
        self.TT=TT
        self.CJO=CJO
        self.VJ=VJ
        self.M=M
        self.EG=EG
        self.XTI=XTI
        self.KF=KF
        self.AF=AF
        self.FC=FC
        self.BV=BV
        self.IBV=IBV
        self.T=T

        k = 1.380649e-23
        q = 1.60217663e-19

        self.Vt = (k * self.T) / q

        self.GMIN = 1e-15

        assert self.IBV >= (self.IS * self.BV / self.Vt), "Non-convergence on regions c and d"

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_degen_current(self, x, sys, t, bf, vm_sys_dy):
        return

    def get_current(self, scb):
        scb.i[self.i_x_ref] = scb.x[self.i_x_ref]

    def get_dependent_current(self, scb):

        self.Vd = scb.v[self.ref]

        if self.Vd >= (-5 * self.N * self.Vt): # region (a)
            scb.i[self.i_d_ref] = self.IS * (math.exp(self.Vd / (self.N * self.Vt)) -1) + (self.GMIN * self.Vd)
        elif -self.BV < (self.Vd * -5 * self.N * self.Vt): # region (b)
            scb.i[self.i_d_ref] = -self.IS + self.GMIN * self.Vd
        elif self.Vd == -self.BV: # region (c)
            scb.i[self.i_d_ref] = -self.IBV
        else: # self.Vd < -self.BV region (d)
            scb.i[self.i_d_ref] = self.IS * (math.exp((-self.BV * self.Vd) / self.Vt) -1 + (self.BV / self.Vt))

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
