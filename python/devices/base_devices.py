from enum import Enum
from abc import ABC, abstractmethod

import numpy as np

from scipy import signal

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
        self.instance    = "default"

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
        # A dependent current source has its normal edge-voltage unknown and
        # adds one branch-current unknown plus one constitutive equation.
        #
        # The solver flow is as follows:
        # 1. Each edge calculates its current using the current value of x and sys
        #    only.
        #
        #      - An independent current source calculates its current based on t
        #        and updates I[ref].
        #      - A dependent current source updates I[i_x_ref] = x[i_x_ref].
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
        #        and updates I[ref] = f(V, I).
        #
        # 4. For every dependent source an equation is inserted that into the system of equations
        #    that asserts that I[i_x_ref] = I[ref]. This constrains the guessed value of a
        #    dependent current source to be equal to the true, dependent value.
        #
        # By following the above steps, potential problems of circular dependency are avoided.
        # For example, if I2 = f(V1) and V1 = f(I1) then I2 cannot calculate its value until V1 is
        # known, but V1 cannot calculate its value intil I1 is known. The additional, inserted
        # variable and equation break this loop.

        # Legacy dependent-value reference used by dependent voltage sources.
        self.i_d_ref     = -1

        # Guessed branch-current reference for a dependent current source.
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
        #  - Only category 1 & 3. above are classed as an LTI system input
        self.is_input = False
        self.input_ref = -1

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

    def get_instance(self):
        return self.instance

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

    def set_input_ref(self, ref):
        self.input_ref = ref

    def get_input_ref(self):
        return self.input_ref

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
    def get_dependent_voltage(self, scb):
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

    def get_dependent_voltage(self, scb):
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

    def get_dependent_voltage(self, scb):
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

    def get_dependent_voltage(self, scb):
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

    def get_dependent_voltage(self, scb):
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

    def set_params(self, omega, mag, phi, bias):
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

    def get_dependent_voltage(self, scb):
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

    def get_time_series_input(self, tr):
        return self.bias + (self.mag * np.sin(self.omega * tr + self.phi))

class pulse_voltage_src(TwoPortElement):

    def __init__(self, T=1e-6, t_on=0.5e-6, t_r=0, t_f=0, t_del=0, mag=0.6, bias=0):
        TwoPortElement.__init__(self)

        self.type = ElementType.voltage_src

        self.t_on  = t_on
        self.T     = T
        self.t_r   = t_r
        self.t_f   = t_f
        self.t_del = t_del
        self.mag   = mag
        self.bias  = bias

    def get_voltage(self, scb):

        pulse_t = (scb.t + self.t_del) % self.T

        if(pulse_t < self.t_on):
            if(pulse_t < self.t_r):
                F = (1/self.t_r)*pulse_t
            elif(pulse_t >= self.t_on - self.t_f):
                F = 1 - (1/self.t_f)*(pulse_t - (self.t_on - self.t_r))
            else:
                F = 1
        else:
            F = 0

        scb.v[self.ref] = self.bias + (self.mag * F)
        #print (scb.v[self.ref])
        #print (scb.t)
        #scb.v[self.ref] = 0.0

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

    def get_dependent_current(self, scb):
        return 0.0

    def get_dependent_voltage(self, scb):
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

    def get_dependent_voltage(self, scb):
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

class sine_current_src(TwoPortElement):

    def __init__(self, omega=1e6, mag=0.6, phi=0, bias=0):
        TwoPortElement.__init__(self)

        self.type = ElementType.current_src

        self.omega    = omega
        self.mag      = mag
        self.phi      = phi
        self.bias     = bias

    def set_params(self, omega, mag, phi, bias):
        self.omega    = omega
        self.mag      = mag
        self.phi      = phi
        self.bias     = bias

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_current(self, scb):
        scb.i[self.ref] = self.bias + (self.mag * math.sin(self.omega * scb.t + self.phi))

    def get_dependent_current(self, scb):
        return 0.0

    def get_dependent_voltage(self, scb):
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

    def get_time_series_input(self, tr):
        return self.bias + (self.mag * np.sin(self.omega * tr + self.phi))

class vccs_l1_mosfet(TwoPortElement):

    def set_params(self, KP, vth, l_lambda, L, W):
        self.KP       = KP
        self.vth      = vth
        self.l_lambda = l_lambda
        self.L        = L
        self.W        = W

        self.diode_threshold = 0.7
        self.diode_is        = 1e-14
        self.vt_diode        = 0.026

        # Smoothing is deliberately expressed in volts.  A continuously
        # differentiable device law is much easier for Newton-like circuit
        # solvers to handle than hard region changes.
        self.vov_smoothing   = 25e-3
        self.vds_smoothing   = 1e-3
        self.gmin            = 1e-12

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

    def check_diode_threshold(self):

        return self.vds < self.diode_threshold

    def diode(self):

        if self.check_diode_threshold():
            self.diode_id = self.diode_is * (np.exp((-self.vds - self.diode_threshold) / self.vt_diode) - 1)
        else:
            self.diode_id = 0

    @staticmethod
    def _softplus(x, width):
        """Overflow-safe approximation of max(x, 0)."""
        return width * np.logaddexp(0.0, x / width)

    @staticmethod
    def _smooth_positive(x, width):
        return 0.5 * (x + np.sqrt(x * x + width * width))

    def sigmoid(self, x, k=50.0):
        # Equivalent to scipy.special.expit without another dependency here.
        z = np.clip(k * x, -60.0, 60.0)
        return 1.0 / (1.0 + np.exp(-z))

    def region_weights(self):
        vgs = self.vgs
        vds = self.vds

        vov = vgs - self.vth          # overdrive
        sat_boundary = vds - vov

        k = 50.0  # smoothing sharpness

        # turn-on (off → on)
        w_on = self.sigmoid(vov, k)

        # triode vs sat split
        w_sat = self.sigmoid(sat_boundary, k) * w_on
        w_tri = (1.0 - self.sigmoid(sat_boundary, k)) * w_on

        # off weight
        w_off = 1.0 - w_on

        return w_tri, w_sat, w_off

    def sat(self):

        #if self.vds < 0 or self.vgs < 0:
        #    return 0

        #vdssat = self.vgs - self.vth
        vdssat = np.log(1 + np.exp(self.vgs - self.vth))  # softplus

        id = (self.KP / 2) * (self.W / self.L) * \
            (self.vgs - self.vth)**2 * (1 + self.l_lambda * (self.vds - vdssat))

        return id

    def tri(self):

        #if self.vds < 0 or self.vgs < 0:
        #    return 0

        id = self.KP * (self.W / self.L) * \
            ( ((self.vgs - self.vth) * self.vds) - (self.vds**2 / 2))

        return id

    def _forward_current_and_derivatives(self, vds, vgs):
        """Smooth Level-1-like current for non-negative drain voltage.

        tanh limits the effective channel voltage to the overdrive, producing
        one differentiable expression for linear and saturation operation.
        """
        beta = self.KP * (self.W / self.L)
        threshold_arg = (vgs - self.vth) / self.vov_smoothing
        vov = self._softplus(vgs - self.vth, self.vov_smoothing)
        dvov_dvgs = self.sigmoid(threshold_arg, k=1.0)
        vds_pos = self._smooth_positive(vds, self.vds_smoothing)
        root = np.sqrt(vds * vds + self.vds_smoothing**2)
        dvds_pos_dvds = 0.5 * (1.0 + vds / root)
        # At extreme trial voltages softplus can underflow to zero.  tanh is
        # already fully saturated beyond this range, so clipping preserves the
        # device law while avoiding inf*0 in its derivative.
        ratio = np.clip(vds_pos / max(vov, 1e-15), 0.0, 40.0)
        tanh_ratio = np.tanh(ratio)
        sech2_ratio = 1.0 - tanh_ratio * tanh_ratio
        vds_eff = vov * tanh_ratio
        dveff_dvds_pos = sech2_ratio
        dveff_dvov = tanh_ratio - ratio * sech2_ratio
        channel = beta * (vov * vds_eff - 0.5 * vds_eff * vds_eff)
        modulation = 1.0 + self.l_lambda * vds_pos
        current = channel * modulation

        dchannel_dvov = beta * (
            vds_eff + (vov - vds_eff) * dveff_dvov
        )
        dchannel_dvds_pos = beta * (
            (vov - vds_eff) * dveff_dvds_pos
        )
        dcurrent_dvgs = dchannel_dvov * modulation * dvov_dvgs
        dcurrent_dvds = (
            dchannel_dvds_pos * modulation + channel * self.l_lambda
        ) * dvds_pos_dvds
        return current, dcurrent_dvds, dcurrent_dvgs

    def _forward_current(self, vds, vgs):
        return self._forward_current_and_derivatives(vds, vgs)[0]

    def current_and_derivatives(self, vds, vgs):
        """Return Ids, gds and gm from one analytic device evaluation."""
        forward, forward_gds, forward_gm = \
            self._forward_current_and_derivatives(vds, vgs)
        reverse, reverse_gds, reverse_gm = \
            self._forward_current_and_derivatives(-vds, vgs - vds)
        ids = forward - reverse + self.gmin * vds
        gds = forward_gds + reverse_gds + reverse_gm + self.gmin
        gm = forward_gm - reverse_gm
        return ids, gds, gm

    def drain_current(self, vds, vgs):
        """Source/drain-symmetric, continuously differentiable drain current."""
        return self.current_and_derivatives(vds, vgs)[0]

    def small_signal_params(self, vds, vgs):
        """Differentiate the same law used by the DC solver."""
        _, gds, gm = self.current_and_derivatives(vds, vgs)
        return gm, max(gds, self.gmin)

    def get_voltage(self, scb):
        scb.v[self.ref] = scb.x[self.ref]

    def get_degen_current(self, x, sys, t, bf, vm_sys_dy):
        return

    def get_current(self, scb):
        scb.i[self.i_x_ref] = scb.x[self.i_x_ref]

    def get_dependent_current(self, scb):
        self.vds = scb.v[self.ref]
        self.vgs = scb.v[self.vgs_ref]
        # Reuse the original branch's I-vector slot for the calculated model
        # current.  The same index in V/x remains the branch-voltage unknown.
        scb.i[self.ref] = self.drain_current(self.vds, self.vgs)

    def get_dependent_voltage(self, scb):
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

    def get_op(self, x):

        self.vds = x[self.ref]
        self.vgs = x[self.vgs_ref]

        return self.vds, self.vgs, x[self.i_x_ref]

    def get_op_t(self, x):

        vds = x[self.ref]
        vgs = x[self.vgs_ref]
        gm, gds = self.small_signal_params(vds, vgs)
        return gm, 1.0 / gds

    def get_weight(self):
        return 3

    def set_dependencies(self, ckt):

        self.i_x_ref = ckt.num_edges

        # Add one branch-current unknown and its constitutive equation.
        ckt.qf = np.vstack([ckt.qf, np.zeros( ckt.qf.shape[0])])
        ckt.qf = np.hstack((ckt.qf, np.zeros((ckt.qf.shape[0], 1))))

        ckt.bf = np.vstack([ckt.bf, np.zeros( ckt.bf.shape[0])])
        ckt.bf = np.hstack((ckt.bf, np.zeros((ckt.bf.shape[0], 1))))

        ckt.qf[:, self.i_x_ref] = ckt.qf[:, self.ref]
        ckt.qf[:, self.ref]     = np.zeros(ckt.qf.shape[0])

        # i_guess - i_model(Vds, Vgs) = 0.  I[self.ref] contains
        # the locally calculated model current.
        ckt.qf[self.i_x_ref][self.i_x_ref] = 1
        ckt.qf[self.i_x_ref][self.ref] = -1

        ckt.num_edges += 1

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        return

class vccs_l1_mosfet_small(TwoPortElement):

    def __init__(self):
        TwoPortElement.__init__(self)

        self.type = ElementType.current_src

    def get_ro(self):

        vdssat = self.vgs - self.vth

        idsat = (self.KP / 2) * (self.W / self.L) * vdssat**2

        ro = 1 / (self.l_lambda * idsat)
        return ro

    # equation 9.22
    def get_gm(self):
        beta = self.KP * (self.W / self.L)
        gm = math.sqrt(2 * beta * self.Ids)
        return gm

    def set_params(self, KP, vth, l_lambda, L, W):
        self.KP       = KP
        self.vth      = vth
        self.l_lambda = l_lambda
        self.L        = L
        self.W        = W

    def get_op(self):

        return [self.gm, self.ro]

    def get_op_t(self, op, i_x_ref, vgs_ref, vds_ref=None):

        self.Ids = op.i[i_x_ref]
        self.vds = op.v[self.ref if vds_ref is None else vds_ref]
        self.vgs = op.v[vgs_ref]

        # Linearize the exact large-signal law rather than a separate hard-
        # region square-law approximation.
        dc_model = vccs_l1_mosfet()
        dc_model.set_params(self.KP, self.vth, self.l_lambda, self.L, self.W)
        self.gm, gds = dc_model.small_signal_params(self.vds, self.vgs)
        self.ro = 1.0 / gds

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

        scb.i[self.ref] = self.gm * self.vgs

    def get_dependent_voltage(self, scb):
        return 0.0

    def get_weight(self):
        return 3

    def set_dependencies(self, ckt):

        self.i_x_ref = ckt.num_edges

        # Add one branch-current unknown and its constitutive equation.
        ckt.qf = np.vstack([ckt.qf, np.zeros( ckt.qf.shape[0])])
        ckt.qf = np.hstack((ckt.qf, np.zeros((ckt.qf.shape[0], 1))))

        ckt.bf = np.vstack([ckt.bf, np.zeros( ckt.bf.shape[0])])
        ckt.bf = np.hstack((ckt.bf, np.zeros((ckt.bf.shape[0], 1))))

        ckt.qf[:, self.i_x_ref] = ckt.qf[:, self.ref]
        ckt.qf[:, self.ref]     = np.zeros(ckt.qf.shape[0])

        ckt.qf[self.i_x_ref][self.i_x_ref] = 1
        ckt.qf[self.i_x_ref][self.ref] = -1

        ckt.num_edges += 1

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        return

    def get_value(self):
        return self.gm

class vcvs_src(TwoPortElement):

    def __init__(self, func=None):
        TwoPortElement.__init__(self)

        self.type = ElementType.voltage_src

        # Store callback
        self.func = func if func is not None else self.default_function

    def default_function(self, scb):
        return 0  # fallback behavior

    def get_voltage(self, scb):
        scb.v[self.i_x_ref] = scb.x[self.i_x_ref]

    def get_current(self, scb):
        scb.i[self.ref] = scb.x[self.ref]

    def get_dependent_current(self, scb):
        pass

    def get_dependent_voltage(self, scb):
        scb.v[self.i_d_ref] = self.func(scb)

    def get_degen_current(self, scb):
        return

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

        ckt.bf[:, self.i_x_ref] = ckt.bf[:, self.ref]
        ckt.bf[:, self.ref]     = np.zeros(ckt.bf.shape[0])

        # Update Bf for degen
        ckt.bf[ckt.num_edges][self.i_x_ref]  = 1
        ckt.bf[ckt.num_edges][self.i_d_ref] = -1

        ckt.num_edges += 2

    def get_weight(self):
        return 0

    def get_dy(self, x, sys):
        return sys

    def get_ddt(self, scb):
        pass
