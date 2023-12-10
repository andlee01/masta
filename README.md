# Electrical Network Equation Generation

## Introduction
Masta is a tool that automates the generation of electrical networks circuit equations. It can be used to solve both linear and non-linear circuits. Linear circuits can also be translated into LTI systems.

## Example

![Example electrical network.](/assets/circuit1_rser_tex.svg)

Consider the electrical network shown above which is constructed in Masta as follows:

```python
# Create a new Masta circuit
ckt = Circuit()

Vs_val   = 1.0
C2_val   = 1e-9
L_val    = 1e-3
C1_val   = 2e-9
C3_val   = 3e-9
R2_val   = 1e3
Rser_val = 1e-3

# Add Vs
Vs = voltage_src()
Vs.set_value(Vs_val)
Vs.set_instance("Vs")
ckt.add_edge(4, 0, Vs)

# Add C2
C2 = capacitor()
C2.set_value(C2_val)
C2.set_instance("C2")
ckt.add_edge(1, 3, C2)

# Add L
L = inductor()
L.set_value(L_val)
L.set_instance("L")
ckt.add_edge(1, 2, L)

# Add C1
C1 = capacitor()
C1.set_value(C1_val)
C1.set_instance("C1")
ckt.add_edge(2, 0, C1)

# Add C3
C3 = capacitor()
C3.set_value(C3_val)
C3.set_instance("C3")
ckt.add_edge(2, 3, C3)

# Add R2
R2 = resistor()
R2.set_value(R2_val)
R2.set_instance("R2")
ckt.add_edge(3, 0, R2)

# Add Rser
Rser = resistor()
Rser.set_value(Rser_val)
Rser.set_instance("Rser")
ckt.add_edge(4, 1, Rser)
```

The circuit can be initilised by calling:

```
ckt.init_circuit()
```

This generates the network KCL and KVL equations. Masta uses the spanning-tree method of circuit analysis. It produces two matrices representing the electrical topology of the network.

* A matrix $`\mathbf{Q_f}`$ representing the KCL basic cut equations.
* A matrix $`\mathbf{B_f}`$ representing the KVL basic loop equations.

### Unknown Variables

Each network element contributes 1 unknown variable into the solution of the $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ equations:

| Element                    | Unknown Variable |
| -------------------------- | ---------------- |
| Capacitor                  |    $` i_C `$     |
| Inductor                   |    $` V_L `$     |
| Resistor                   |    $` i_R `$     |
| Independent Voltage Source |    $` i_{V_S} `$ |
| Independent Current Source |    $` V_{i_S} `$ |

The unknown variable vector of the electrical network shown above can therefore be written as:

$$ \vec{v} = \begin{bmatrix} i_{V_S} \\ i_{C2} \\ V_{L} \\ i_{C2} \\ i_{C3} \\ i_{R2} \\ i_{R_{SER}} \end{bmatrix} $$

### State Vector

Solutions of the $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ equations require capacitor voltages $`V_C`$ and inductor currents $`i_L`$. These quantities are selected as the state vector in an initial value problem (IVP) and can be treated as constants in solutions of the $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ equations.

For a complete description of the state vector, see [State Vector](doc/state_vector.md).

| Element                    | State Vector     |
| -------------------------- | ---------------- |
| Capacitor                  |    $` V_C `$     |
| Inductor                   |    $` i_L `$     |

The state vector of the electrical network shown above can therefore be written as:

$$ \vec{x} = \begin{bmatrix} V_{C1} \\ V_{C2} \\ V_{C3} \\ i_L \end{bmatrix} $$

### Independent Sources

Independent source values are functions of time only and therefore can be treated as constants in solutions of the $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ equations.

| Element                    | Constant         |
| -------------------------- | ---------------- |
| Independent Voltage Source |    $` V_{V_S} `$ |
| Independent Current Source |    $` i_{i_S} `$ |

Dependent sources require more consideration and are described in more detail in [Dependent Sources](doc/dependent_sources.md).

### $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ Solution

Formulating the vectors $`\vec{v}`$, $`\vec{x}`$ and, along with the current time $`t`$ allows for the definition of the branch voltage and current vectors ($`\vec{I}`$ and $`\vec{V}`$). The full set of network equations can now be written as:

$$ \begin{eqnarray}
\mathbf{Q_f} \cdot \vec{I} = 0  \\
\mathbf{B_f} \cdot \vec{V} = 0 \end{eqnarray} $$

The above equations can be solved in Masta using a Python equation solver:

```python
def circuit_eqn(v, x, ckt, t):

    # current Im
    I = ckt.get_im(x=v, sys=x, t=t)

    # voltage Vm
    V = ckt.get_vm(x=v, sys=x, t=t)

    # Copy matrices
    qf_num = ckt.qf.copy()
    bf_num = ckt.bf.copy()

    eqn_qf = np.dot(qf_num, I)
    eqn_bf = np.dot(bf_num, V)

    return eqn_qf + eqn_bf

def dypc(t, x, ckt):

    # Initial value (guess) for solution
    root_start = np.ones(ckt.num_edges)

    # Solve for I, V
    root = optimize.root(circuit_eqn, root_start, args=(x, ckt, t), tol=1e-7)
```

The call to the `Circuit` functions `get_im` and `get_vm` return the branch current and voltage vectors respectively. Each network `TwoPortElement` must implement a `get_current` and `get_voltage` method that calculates the branch current and voltage from the input vectors $`\vec{v}`$ and $`\vec{x}`$ and current time $`t`$.

### Formulation of IVP
Solving the $`\mathbf{Q_f}`$ and $`\mathbf{B_f}`$ equations will calculate capacitor currents and inductor voltages. Noting that:

$$ \begin{eqnarray}
i_C  = C \frac{dV_C}{dt}  \\
V_L  = L \frac{di_L}{dt} \end{eqnarray} $$

Then the definition of the derivative of the state vector follows:

$$
\dot{\vec{x}} = \begin{bmatrix} \frac{i_{C_0}}{C_0} \\ ...\\ \frac{i_{C_n}}{C_n} \\ \frac{V_{L_0}}{L_0} \\ ... \\ \frac{V_{L_n}}{L_n} \end{bmatrix}
$$

A python ODE solver can be used to simulate the transient behaviour of the electrical network.

```python
def dypc(t, x, ckt):

    # Initial value (guess) for solution at time t
    root_start = np.ones(ckt.num_edges)

    # Solve for I, V
    root = optimize.root(circuit_eqn, root_start, args=(x, ckt, t), tol=1e-7)

    # Calculate the state vector derivative
    return ckt.get_dy(x=root.x)

def ode_solve(ckt):

    num_sys_vars = ckt.get_num_sys_vars()

    # Simulate from t=0 to t=50e-6 using 10000 time steps
    tstep = 10000
    t     = np.linspace(0, 50e-6, tstep)

    # Initial value of state vector (x)
    x0    = np.zeros (num_sys_vars)

    r = ode(dypc).set_integrator('lsoda', method='bdf', atol=1e-9, rtol=1e-9)
    r.set_initial_value(x0, 0.0)
    r.set_f_params(ckt)

    y    = np.empty((tstep, num_sys_vars))
    y[0] = x0

    k = 1
    while r.successful() and k < tstep:
        r.integrate(t[k])

        y[k] = r.y

        k += 1

    return y
```

The results of the simulation can be interrogated and plotted for any branch:

```python
for idx in range(len(t)):
  VC3[idx] = ckt.scb.v[ckt.get_edge_info(ckt.get_ref_from_instance("C3")).ref]

  iR2[idx] = ckt.scb.i[ckt.get_edge_info(ckt.get_ref_from_instance("R2")).ref]
```

As an example, plotting the state vector transient response:

![Example electrical network.](/assets/litmus0_rser1.svg)