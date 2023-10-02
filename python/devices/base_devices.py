from enum import Enum
from abc import ABC, abstractmethod

import numpy as np

class ElementType(Enum):
    undefined   = 0
    resistor    = 1
    capacitor   = 2
    inductor    = 3
    voltage_src = 4
    current_src = 5

class TwoPortElement(ABC):

    def __init__(self):
        self.value       = 0.0
        self.type        = ElementType.undefined
        self.sys_var     = False
        self.sys_var_ref = -1
        self.degen_ref   = -1

    def set_value(self, value):
        self.value = value

    def set_type(self, type):
        self.type = type

    def get_type(self):
        return self.type

    def set_instance(self, instance):
        self.instance = instance

    def set_ref(self, ref):
        self.ref = ref

    def set_sys_var(self, ref):
        self.sys_var     = True
        self.sys_var_ref = ref

    def set_degen_ref(self, ref):
        self.degen_ref = ref

    def get_ref(self):
        return self.ref

    def get_degen_ref(self):
        return self.degen_ref

    @abstractmethod
    def get_voltage(self, scb):
        pass

    @abstractmethod
    def get_current(self, scb):
        pass

    @abstractmethod
    def get_degen_current(self, scb):
        pass

    @abstractmethod
    def get_dy(self, x, sys):
        pass

    @abstractmethod
    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):
        pass

    @abstractmethod
    def get_weight(self):
        pass

class resistor(TwoPortElement):

    # sys vector variable is current through resistor
    #  - V = i * R
    def get_voltage(self, scb):
        return scb.x[self.ref] * self.value

    def get_current(self, scb):
        return scb.x[self.ref]

    def get_degen_current(self, scb):
        return

    def get_weight(self):
        return 2

    def get_dy(self, x, sys):
        return sys

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):
        return


class capacitor(TwoPortElement):

    def get_voltage(self, scb):

        if self.sys_var_ref != -1:
            return scb.sys[self.sys_var_ref]

        # non sys-var cap
        #  - Voltage is dependent variable and not solved for
        return 0

    def get_current(self, scb):
        return scb.x[self.ref]

    def get_degen_current(self, scb):

        bf_loop = scb.degen_mtrx[self.ref,:].copy()
        bf_loop[self.ref] = 0

        # i = C *ddt(V)
        # Multiply by -1 to account for moving of ddt to RHS
        #  - i.e. ddt(Vc1) + ddt(Vc2) + ddt(V) = 0
        #                               ddt(V) = -1 * ( ddt(Vc1) + ddt(Vc2) )
        return -1 * np.dot(bf_loop, scb.dv) * self.value

    def get_weight(self):
        return 0

    def get_dy(self, x, sys):

        if self.sys_var_ref != -1:
            sys[self.sys_var_ref] = x[self.ref] / self.value

        return sys

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):
        return

class inductor(TwoPortElement):

    def get_voltage(self, scb):
        return scb.x[self.ref]

    def get_current(self, scb):
        return scb.sys[self.sys_var_ref]

    def get_degen_current(self, scb):
        return

    def get_weight(self):
        return 4

    def get_dy(self, x, sys):
        sys[self.sys_var_ref] = x[self.ref] / self.value

        return sys

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):
        return

class voltage_src(TwoPortElement):

    def get_voltage(self, scb):
        return self.value

    def get_current(self, scb):
        return scb.x[self.ref]

    def get_degen_current(self, scb):
        return

    def get_weight(self):
        return 0

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):

        linalg_b             += -1 * self.get_voltage(x=x, t=t) * linalg_bf[:,self.ref]
        linalg_bf[:,self.ref] = 0

        return [linalg_qf, linalg_bf, linalg_b]

    def get_dy(self, x, sys):
        return sys

class current_src(TwoPortElement):

    def get_voltage(self, scb):
        return scb.x[self.ref]

    def get_current(self, scb):
        return self.value

    def get_degen_current(self, scb):
        return

    def get_weight(self):
        return 3

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):

        linalg_b             += -1 * self.get_current(x=x, t=t) * linalg_qf[:,self.ref]
        linalg_qf[:,self.ref] = 0

        return [linalg_qf, linalg_bf, linalg_b]

    def get_dy(self, x, sys):
        return sys

class vccs_l1_mosfet(TwoPortElement):

    def set_params(self, KP, vth, l_lambda, L, W):
        self.KP       = KP
        self.vth      = vth
        self.l_lambda = l_lambda
        self.L        = L
        self.W        = W

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

        if self.vds < 0 or self.vgs < 0:
            return 0

        vdssat = self.vgs - self.vth

        id = (self.KP / 2) * (self.W / self.L) * \
            (self.vgs - self.vth)**2 * (1 + self.l_lambda * (self.vds - vdssat))

        return id

    def tri(self):

        if self.vds < 0 or self.vgs < 0:
            return 0

        id = self.KP * (self.W / self.L) * \
            ( ((self.vgs - self.vth) * self.vds) - (self.vds**2 / 2))

        return id

    def get_voltage(self, scb):
        return scb.x[self.ref]

    def get_degen_current(self, x, sys, t, bf, vm_sys_dy):
        return

    def get_current(self, scb):

        self.vds = scb.x[self.ref]
        self.vgs = scb.x[self.vgs_ref]

        if self.check_tri():
            return self.tri()
        elif self.check_sat():
            return self.sat()
        else:
            return (1e-9 * self.vds)
            return 0.0

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

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):

        v_step = 1e-3

        x_upper = x + v_step
        val_upper = self.get_current(x=x_upper, t=t)

        x_lower = x - v_step
        val_lower = self.get_current(x=x_lower, t=t)

        geq = (val_upper - val_lower) / (2 * v_step)

        f_idc = self.get_current(x=x, t=t)
        ieq   = f_idc - geq * x[self.ref]

        linalg_b             += -1 * ieq * linalg_qf[:,self.ref]
        linalg_qf[:,self.ref] = geq * linalg_qf[:,self.ref]

        return [linalg_qf, linalg_bf, linalg_b]

    def get_dy(self, x, sys):
        return sys
