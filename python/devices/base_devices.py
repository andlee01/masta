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

class current_src(TwoPortElement):

    def get_voltage(self, x, t):
        return x[self.ref]

    def get_current(self, x, t):
        return self.value

    def get_weight(self):
        return 3

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
