from enum import Enum
from abc import ABC, abstractmethod

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
        self.dependant   = False

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

    def get_ref(self):
        return self.ref

    def get_dependant(self):
        return self.dependant

    @abstractmethod
    def get_voltage(self, x, t):
        pass

    @abstractmethod
    def get_current(self, x, t):
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
    def get_voltage(self, x, t):
        return x[self.ref] * self.value

    def get_current(self, x, t):
        return x[self.ref]

    def get_weight(self):
        return 2

class voltage_src(TwoPortElement):

    def get_voltage(self, x, t):
        return self.value

    def get_current(self, x, t):
        return x[self.ref]

    def get_weight(self):
        return 0

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):

        linalg_b             += -1 * self.get_voltage(x=x, t=t) * linalg_bf[:,self.ref]
        linalg_bf[:,self.ref] = 0

        return [linalg_qf, linalg_bf, linalg_b]

class current_src(TwoPortElement):

    def get_voltage(self, x, t):
        return x[self.ref]

    def get_current(self, x, t):
        return self.value

    def get_weight(self):
        return 3

    def upd_linalg_mtrx(self, x, sys, t, linalg_qf, linalg_bf, linalg_b):

        linalg_b             += -1 * self.get_current(x=x, t=t) * linalg_qf[:,self.ref]
        linalg_qf[:,self.ref] = 0

        return [linalg_qf, linalg_bf, linalg_b]

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

    def get_voltage(self, x, t):
        return x[self.ref]

    def get_current(self, x, t):

        self.vds = x[self.ref]
        self.vgs = x[self.vgs_ref]

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
