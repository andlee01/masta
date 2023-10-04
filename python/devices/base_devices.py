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

        # Edge reference in circuit
        self.ref         = -1

        # System variable reference
        #  - set to -1 if edge is not a system variable
        self.sys_var_ref = -1

        # Degeneration edge reference
        #  - Additional variable inserted into the circuit to constrain the
        #    current of a dependent capacitor or the voltage of a dependent
        #    inductor
        self.degen_ref   = -1

        # Dependent sources
        # -----------------
        # An independent voltage source adds 1 variable into the system equations.
        # This variable corresponds to the edge current.
        # The edge voltage is a function of time only.
        #
        # An independent current source adds 1 variable into the system equations.
        # This variable corresponds to the edge voltage.
        # The edge current is a function of time only.
        #
        # A dependent current source adds 3 variables into the system equations.
        # The first variable is the edge voltage.
        # The second variable is the edge current.
        # The third variable is the edge current constraint and is calculated from
        # the guess edge voltages and currents.
        #
        # The solver flow is as follows:
        # 1. Each edge calculates its current using the current value of x and sys
        #    only.
        #
        #      - An independent current source calculates its current based on t
        #        and updates I[ref].
        #      - A dependent current source updates I[i_d_ref] = x[i_d_ref].
        #
        # 2. Each edge calculates its voltage using the current value of x and sys
        #    only.
        #
        #      - All current sources updates V[ref] = x[ref].
        #
        # 3. Dependent current sources calculate their dependent current values using
        #    the I and V vector from steps 1. and 2. only.
        #
        #      - Dependent current sources calculate the dependent current using V and I
        #        and updates I[i_x_ref] = f(V, I)
        #
        # 4. For every dependent source an equation is inserted that into the system of equations
        #    that asserts that I[i_d_ref] = I[i_x_ref]. This constrains the guessed value of a
        #    dependent current source to be equal to the true, dependent value.
        #
        # By following the above steps, potential problems of circular dependency are avoided.
        # For example, if I2 = f(V1) and V1 = f(I1) then I2 cannot calculate its value until V1 is
        # known, but V1 cannot calculate its value intil I1 is known. The additional, inserted
        # variables break this loop.

        # Dependent variable reference
        #  - Guessed current for a dependent current source or voltage for a
        #    dependent voltage source
        self.i_d_ref     = -1

        # An additional constraint variable that's inserted alongside i_d_ref
        #  - Calculates the dependent true value based on the guessed voltages
        #    and currents
        self.i_x_ref     = -1

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
    def set_dependencies(self, ckt):
        pass

    @abstractmethod
    def get_voltage(self, scb):
        pass

    @abstractmethod
    def get_current(self, scb):
        pass

    @abstractmethod
    def get_dependent_current(self, scb):
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

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):
        return

    def set_dependencies(self, ckt):
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
        bf_loop = scb.degen_mtrx[self.ref,:].copy()
        bf_loop[self.ref] = 0

        # i = C *ddt(V)
        # Multiply by -1 to account for moving of ddt to RHS
        #  - i.e. Vc1 + Vc2 + V = 0
        #                     V = -1 * ( Vc1 + Vc2 )
        return -1 * np.dot(bf_loop, scb.v)

    def get_current(self, scb):
        return scb.x[self.ref]

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):

        bf_loop = scb.degen_mtrx[self.ref,:].copy()
        bf_loop[self.ref] = 0

        # i = C *ddt(V)
        # Multiply by -1 to account for moving of ddt to RHS
        #  - i.e. ddt(Vc1) + ddt(Vc2) + ddt(V) = 0
        #                               ddt(V) = -1 * ( ddt(Vc1) + ddt(Vc2) )
        return -1 * np.dot(bf_loop, scb.dv) * self.value

    def set_dependencies(self, ckt):
        return

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

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):
        return

    def set_dependencies(self, ckt):
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

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):
        return

    def set_dependencies(self, ckt):
        return

    def set_dependencies(self, ckt):
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

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):
        return

    def set_dependencies(self, ckt):
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
        return scb.x[self.i_x_ref]

    def get_dependent_current(self, scb):

        self.vds = scb.v[self.ref]
        self.vgs = scb.v[self.vgs_ref]

        if self.check_tri():
            return self.tri()
        elif self.check_sat():
            return self.sat()
        else:
            return (1e-9 * self.vds)

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
