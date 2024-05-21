from enum import Enum
from abc import ABC, abstractmethod

import numpy as np

import math

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

        # For LTI Circuits only
        #  - Dependency reference for dependent source
        self.src_dep_ref = -1

        # For LTI Circuits only
        #  - Voltage and current sources can be either:
        #    1. Independent current/voltage sources
        #    2. Dependent current/voltage sources
        #    3. Constant valued current/voltage sources
        #
        #  - Only category 1. above is classed as an LTI system input
        self.is_input = False

        # Category 3
        self.is_const = False

        # For LTI Circuits only
        self.is_output  = False
        self.output_ref = -1

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

    def get_value(self):
        return self.value

    def get_sys_var_ref(self):
        return self.sys_var_ref

    def get_degen_ref(self):
        return self.degen_ref

    def set_src_dep_ref(self, ref):
        self.src_dep_ref = ref

    def get_src_dep_ref(self):
        return self.src_dep_ref

    def set_is_input(self):
        self.is_input = True

    def get_is_input(self):
        return self.is_input

    def set_is_output(self):
        self.is_output = True

    def get_is_output(self):
        return self.is_output

    def set_output_ref(self, ref):
        self.output_ref = ref

    def get_output_ref(self):
        return self.output_ref

    def set_is_const(self):
        self.is_const = True

    def get_is_const(self):
        return self.is_const

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
    def get_ddt(self, scb):
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
    def get_weight(self):
        pass

class resistor(TwoPortElement):

    def __init__(self, ramp=False, ramp_ddt=1e6):
        TwoPortElement.__init__(self)

        self.type = ElementType.resistor

    # sys vector variable is current through resistor
    #  - V = i * R
    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref] * self.value

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

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

    def get_ddt(self, scb):
        return


class capacitor(TwoPortElement):

    def __init__(self):
        TwoPortElement.__init__(self)

        self.type = ElementType.capacitor

    def get_voltage(self, scb):

        if self.sys_var_ref != -1:
            scb.v[self.ref] = scb.sys[self.sys_var_ref]

            return

        # non sys-var cap
        #  - Voltage is dependent variable and not solved for
        bf_loop = scb.degen_mtrx[self.ref,:].copy()
        bf_loop[self.ref] = 0

        # i = C *ddt(V)
        # Multiply by -1 to account for moving of ddt to RHS
        #  - i.e. Vc1 + Vc2 + V = 0
        #                     V = -1 * ( Vc1 + Vc2 )
        scb.v[self.ref] = -1 * np.dot(bf_loop, scb.v)

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):

        bf_loop = scb.degen_mtrx[self.ref,:].copy()
        bf_loop[self.ref] = 0

        # i = C *ddt(V)
        # Multiply by -1 to account for moving of ddt to RHS
        #  - i.e. ddt(Vc1) + ddt(Vc2) + ddt(V) = 0
        #                               ddt(V) = -1 * ( ddt(Vc1) + ddt(Vc2) )
        scb.i[self.degen_ref] = -1 * np.dot(bf_loop, scb.dv) * self.value

    def set_dependencies(self, ckt):
        return

    def get_weight(self):
        return 1

    def get_dy(self, x, sys):
        if self.sys_var_ref != -1:
            sys[self.sys_var_ref] = x[self.ref] / self.value

        return sys

    def get_ddt(self, scb):
        scb.dv[self.ref] = scb.x[self.ref] / self.value

class inductor(TwoPortElement):

    def __init__(self, ramp=False, ramp_ddt=1e6):
        TwoPortElement.__init__(self)

        self.type = ElementType.inductor

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_current(self, scb):
        scb.i[self.ref] = scb.sys[self.sys_var_ref]

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

    def get_ddt(self, scb):
        return

class voltage_src(TwoPortElement):

    def __init__(self, ramp=False, ramp_ddt=1e6, delay=0):
        TwoPortElement.__init__(self)

        self.type = ElementType.voltage_src

        self.ramp     = ramp
        self.ramp_ddt = ramp_ddt
        self.delay    = delay

    def get_voltage(self, scb):

        if self.ramp:

            if scb.t > self.delay:
                vin = self.ramp_ddt * (scb.t - self.delay)
                if vin > self.value:
                    vin = self.value
            else:
                vin = 0
        else:
            vin = self.value

        scb.v[self.ref] = vin

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

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

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):

        if scb.v[self.ref] < self.value:
            scb.dv[self.ref] = self.ramp_ddt
        else:
            scb.dv[self.ref] = 0

class sine_voltage_src(TwoPortElement):

    def __init__(self, omega=1e6, mag=0.6, phi=0, bias=0):
        TwoPortElement.__init__(self)

        self.type = ElementType.voltage_src

        self.omega    = omega
        self.mag      = mag
        self.phi      = phi
        self.bias     = bias

    def get_voltage(self, scb):
        scb.v[self.ref] = self.bias + (self.mag * math.sin(self.omega * scb.t + self.phi))

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

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

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        scb.dv[self.ref] = 0

class current_src(TwoPortElement):

    def __init__(self):
        TwoPortElement.__init__(self)

        self.type = ElementType.current_src

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_current(self, scb):
        scb.i[self.ref] = self.value

    def get_dependent_current(self, scb):
        return 0.0

    def get_degen_current(self, scb):
        return

    def set_dependencies(self, ckt):
        return

    def get_weight(self):
        return 3

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        return

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

    def print_stats(self):
        print (self.vgs)
        print (self.vds)

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

class vccs_l1_mosfet_small(TwoPortElement):

    def __init__(self):
        TwoPortElement.__init__(self)

        self.type = ElementType.current_src

    def ro(self):

        vdssat = self.vgs - self.vth

        idsat = (self.KP / 2) * (self.W / self.L) * vdssat**2

        ro = 1 / (self.l_lambda * idsat)
        return ro

    # equation 9.22
    def gm(self):
        beta = self.KP * (self.W / self.L)

        gm = math.sqrt(2 * beta * self.Ids)
        return gm

    def set_params(self, KP, vth, l_lambda, L, W):
        self.KP       = KP
        self.vth      = vth
        self.l_lambda = l_lambda
        self.L        = L
        self.W        = W

    def get_op(self, op, strt):

        #self.Ids = op.i[self.i_x_ref]
        #self.vgs = op.v[self.vgs_ref]


        self.Ids = op.i[strt]
        self.Ids = 20e-6
        self.vgs = op.v[strt+1]

        self.gm       = self.gm()
        self.ro       = self.ro()

        return [self.gm, self.ro]

    def get_op_t(self, op, i_x_ref, vgs_ref):

        self.Ids = op.i[i_x_ref]
        self.vgs = op.v[vgs_ref]

        self.gm       = self.gm()
        self.ro       = self.ro()

        return [self.gm, self.ro]

    def set_vgs_ref(self, vgs_ref):
        self.vgs_ref = vgs_ref

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_degen_current(self, x, sys, t, bf, vm_sys_dy):
        return

    def get_current(self, scb):
        scb.i[self.i_x_ref] = scb.x[self.i_x_ref]

    def get_dependent_current(self, scb):

        self.vgs = scb.v[self.vgs_ref]

        scb.i[self.i_d_ref] = self.gm * self.vgs

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

    def get_value(self):
        return self.gm
