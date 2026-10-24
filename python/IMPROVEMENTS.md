# Planned simulator improvements

## Implicit transient integration with companion models

Replace the current explicit ODE-based transient path with a SPICE-like
implicit timestep loop using capacitor and inductor companion models.

This does **not** require changing the circuit topology equations:

$$
Q_f I + B_f V = 0
$$

`Qf`, `Bf`, and the Kirchhoff residual remain unchanged. Only the transient
constitutive laws used to populate `I` and `V` change.

### Backward-Euler capacitor companion

Starting with:

$$
i_C = C\frac{dv_C}{dt}
$$

Backward Euler gives:

$$
v_C^{n+1} = v_C^n + \frac{\Delta t}{C}i_C^{n+1}
$$

In the existing branch-current formulation, treat the capacitor like a
timestep-dependent Thevenin companion:

```python
v[ref] = previous_voltage + (dt / capacitance) * x[ref]
```

The Jacobian stamp is:

```python
dcurrent[ref, ref] = 1.0
dvoltage[ref, ref] = dt / capacitance
```

### Backward-Euler inductor companion

Starting with:

$$
v_L = L\frac{di_L}{dt}
$$

Backward Euler gives:

$$
i_L^{n+1} = i_L^n + \frac{\Delta t}{L}v_L^{n+1}
$$

For an inductor whose branch variable is voltage:

```python
i[ref] = previous_current + (dt / inductance) * x[ref]
```

The Jacobian stamp is:

```python
dvoltage[ref, ref] = 1.0
dcurrent[ref, ref] = dt / inductance
```

### Proposed transient loop

1. Choose `t_next` and `dt`.
2. Keep capacitor-voltage and inductor-current history fixed at the last
   accepted timestep.
3. Evaluate independent sources at `t_next`; source derivatives are not
   required.
4. Stamp dynamic companions and nonlinear-device companions.
5. Solve the existing `Qf I + Bf V = 0` residual with damped Newton iteration.
6. If Newton converges and the local error is acceptable, commit the new
   capacitor voltages and inductor currents.
7. If it fails, restore the last accepted state, reduce `dt`, and retry.
8. Adapt the next timestep using a local truncation-error estimate and known
   source breakpoints.

History must never be updated during Newton iterations. It is committed only
after a timestep has been accepted.

### Follow-on work

- Start with backward Euler because it is robust and numerically dissipative.
- Add trapezoidal or Gear/BDF integration after the basic timestep controller
  is reliable.
- Add timestep rejection, rollback, minimum/maximum timestep limits, and
  source breakpoint scheduling.
- Define consistent initial-condition handling for capacitor-only loops,
  inductor cutsets, and ideal source discontinuities.
- Consider small ESR or other physical parasitics when an ideal voltage step
  would otherwise imply an impulsive capacitor current.

## Charge-based voltage-dependent capacitors

Implement nonlinear capacitors from stored charge rather than freezing a
voltage-dependent capacitance and treating it as linear:

$$
i_C = \frac{dQ}{dt}
$$

This avoids ambiguity in the meaning of `C(V)` and ensures that transient
current corresponds to the actual change in stored model charge.

### Clarify what `C(V)` means

For a linear capacitor, `Q(V)` is a straight line and the following two
definitions have the same value:

$$
C_{\mathrm{static}}(V)=\frac{Q(V)}{V}
$$

$$
C_{\mathrm{diff}}(V)=\frac{dQ}{dV}
$$

For a nonlinear capacitor, `Q(V)` is curved and they are different:

- `Q(V) / V` is the average slope from the origin to the operating point.
- `dQ / dV` is the local slope of the charge curve at the operating point.

The word "capacitance" is therefore ambiguous unless the model states which
quantity it supplies.

If `C(V)` is differential capacitance,

$$
C_{\mathrm{diff}}(V)=\frac{dQ}{dV}
$$

then:

$$
Q(V)=Q(V_0)+\int_{V_0}^{V}C_{\mathrm{diff}}(u)\,du
$$

and:

$$
i=C_{\mathrm{diff}}(V)\frac{dV}{dt}
$$

If a model instead defines charge as

$$
Q(V)=C(V)V
$$

then its differential capacitance is:

$$
\frac{dQ}{dV}=C(V)+V\frac{dC}{dV}
$$

Using only `C(V) * dV/dt` in the latter case omits the second term. Prefer an
explicit `Q(V)` model plus its derivative `dQ/dV`.

#### Worked example

Let:

$$
C(V)=C_0(1+\alpha V)
$$

If this function means differential capacitance,

$$
C(V)=\frac{dQ}{dV}
$$

then integrate it to obtain charge:

$$
Q(V)=\int_0^V C_0(1+\alpha u)\,du
=C_0\left(V+\frac{\alpha V^2}{2}\right)
$$

and:

$$
i=C_0(1+\alpha V)\frac{dV}{dt}
$$

If the model instead explicitly defines:

$$
Q(V)=C(V)V
$$

then:

$$
Q(V)=C_0(V+\alpha V^2)
$$

and:

$$
\frac{dQ}{dV}=C_0(1+2\alpha V)
$$

so:

$$
i=C_0(1+2\alpha V)\frac{dV}{dt}
$$

For example, with `C0 = 1 nF`, `alpha = 1/V`, and `V = 1 V`, the supplied
function evaluates to `2 nF`. If it is differential capacitance, the stored
charge is `1.5 nC` and the local capacitance is `2 nF`. If it is the factor in
`Q = C(V)V`, the stored charge is `2 nC` and the local differential
capacitance is `3 nF`.

Voltage-dependent capacitance curves such as varactor curves usually specify
differential or small-signal capacitance. In that case, integrate the curve to
obtain charge. Only use `Q = C(V)V` when that relationship is explicitly part
of the model definition.

Prefer an API that makes both quantities explicit:

```python
charge = model.get_charge(voltage)
differential_capacitance = model.get_dcharge_dvoltage(voltage)
```

### Backward-Euler charge equation

For a two-terminal nonlinear capacitor:

$$
i_C^{n+1}=\frac{Q(V_C^{n+1})-Q(V_C^n)}{\Delta t}
$$

Using an auxiliary branch-current unknown, add the constitutive constraint:

$$
F_C=i_C-\frac{Q(V_C)-Q_{\mathrm{previous}}}{\Delta t}=0
$$

The branch-current unknown continues to participate in the normal `Qf` KCL
equations. The capacitor voltage continues to participate in `Bf` KVL
equations.

At Newton iteration `k`, calculate:

$$
C_{\mathrm{diff}}^k=\left.\frac{dQ}{dV}\right|_{V^k}
$$

and stamp the constraint Jacobian as:

```python
jacobian[constraint_row, current_ref] = 1.0
jacobian[constraint_row, voltage_ref] = -(dQ_dV / dt)
```

The equivalent Newton companion is:

$$
i_C \approx G_{\mathrm{eq}}V_C+I_{\mathrm{eq}}
$$

where:

$$
G_{\mathrm{eq}}=\frac{C_{\mathrm{diff}}^k}{\Delta t}
$$

and:

$$
I_{\mathrm{eq}}=
\frac{Q(V^k)-Q_{\mathrm{previous}}}{\Delta t}
-G_{\mathrm{eq}}V^k
$$

At convergence, the numerical charge transfer satisfies:

$$
i_C^{n+1}\Delta t=Q(V_C^{n+1})-Q(V_C^n)
$$

### History and timestep acceptance

Keep `previous_charge` fixed throughout every Newton iteration for a timestep:

```python
candidate_charge = capacitor.get_charge(candidate_controls)
candidate_current = (candidate_charge - previous_charge) / dt
```

Commit the candidate charge only after the timestep is accepted. If Newton or
the local-error check rejects the timestep, discard the candidate state and
retain the last accepted charge.

Suggested interface:

```python
class NonlinearCapacitor:
    def get_charge(self, controls):
        ...

    def get_charge_derivatives(self, controls):
        # Return dQ/dV1, dQ/dV2, ...
        ...

    def commit(self, controls):
        self.previous_charge = self.get_charge(controls)

    def rollback(self):
        ...
```

### Multi-terminal charge models

For terminal charge

$$
Q_k=Q_k(V_1,V_2,\ldots,V_n)
$$

use:

$$
i_k^{n+1}=\frac{Q_k(V^{n+1})-Q_k(V^n)}{\Delta t}
$$

and stamp:

$$
\frac{\partial i_k}{\partial V_j}
=\frac{1}{\Delta t}\frac{\partial Q_k}{\partial V_j}
$$

The derivatives form the nonlinear capacitance matrix:

$$
C_{kj}=\frac{\partial Q_k}{\partial V_j}
$$

For a multi-terminal device, enforce terminal charge conservation:

$$
\sum_k Q_k=0
$$

which guarantees:

$$
\sum_k i_k=0
$$

Prefer calculating all but one terminal charge and defining the final charge
as the negative sum, for example:

```python
q_source = -(q_gate + q_drain + q_bulk)
```

For semiconductor devices, use conserved terminal-charge functions rather
than a collection of independently evaluated voltage-dependent capacitors.

### What charge conservation means

Charge cannot appear or disappear inside an ideal capacitor model. Any charge
entering one terminal must be balanced by charge on the other terminals.

For a two-terminal capacitor, define the charge on the positive terminal as
`Q(V)` and the charge on the negative terminal as `-Q(V)`:

$$
Q_+ + Q_- = Q(V)-Q(V)=0
$$

The corresponding terminal currents are:

$$
i_+ = \frac{dQ_+}{dt}
$$

$$
i_- = \frac{dQ_-}{dt}=-i_+
$$

Therefore the capacitor automatically satisfies KCL:

$$
i_+ + i_- = 0
$$

For a multi-terminal device such as a MOSFET, define terminal charges such as
`Qg`, `Qd`, `Qs`, and `Qb`. A charge-conserving model must satisfy:

$$
Q_g+Q_d+Q_s+Q_b=0
$$

Differentiating this identity gives:

$$
i_g+i_d+i_s+i_b=0
$$

Thus the dynamic currents generated by changing stored charge cannot violate
KCL. A robust implementation can calculate three terminal charges and define
the fourth as the negative sum:

```python
q_source = -(q_gate + q_drain + q_bulk)
```

This enforces conservation even in the presence of floating-point rounding.

#### Conservation across a numerical timestep

With backward Euler, calculate each terminal current from its actual charge
change:

$$
i_k^{n+1}=\frac{Q_k^{n+1}-Q_k^n}{\Delta t}
$$

Summing all terminal currents gives:

$$
\sum_k i_k^{n+1}
=\frac{\sum_k Q_k^{n+1}-\sum_k Q_k^n}{\Delta t}
=0
$$

provided the terminal charges sum to zero at both states. The discrete
transient model therefore conserves charge exactly, apart from numerical
roundoff.

This is stronger than freezing a voltage-dependent capacitance and using:

$$
i \approx C(V)\frac{V^{n+1}-V^n}{\Delta t}
$$

That approximation calculates a current but does not necessarily make the
integrated current equal the device model's actual charge change:

$$
Q(V^{n+1})-Q(V^n)
$$

Over many timesteps, those differences can accumulate into artificial charge.
They can also make independently modelled terminal capacitances produce
currents whose sum is not zero.

#### Implementation rule

Treat terminal charge as the persistent state and current as the consequence
of its change:

```python
candidate_charge = device.get_terminal_charges(candidate_voltages)
terminal_currents = (candidate_charge - accepted_charge) / dt
```

Keep `accepted_charge` fixed throughout Newton iteration. Commit
`candidate_charge` only after accepting the timestep; discard it when a
timestep is rejected. This preserves both charge conservation and correct
rollback behaviour.
