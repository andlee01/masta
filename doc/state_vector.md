# Independent State Variables
It is not always possible to include all capacitor voltages in the state vector. Such situations arise when the electrical network exhibits basic loops that contain only capacitors and independent voltage sources. This results in a capacitor voltage that is not an independent variable. In the example above $`V_{C3} = -V_{IN} + V_{C2} + V_{C1}`$. Dependence amongst the state variables prohibts the arbitary assignment of initial values. Dependent variables therefore can not form part of the state vector (in this instance $`V_{C3}`$). The handling of such variables is described in detail later.

The variable dependence can be avoided by breaking the offending loop with the insertion a small value series resistor. This is shown in the figure below.

![Example electrical network.](/assets/circuit1_rser_tex.svg)

Similar to capacitor voltages, it is not always possible to include all inductor currents in the state vector. Such situations arise when the circuit exhibits nodes that are connected only to inductors and independent current sources. The handling of such variables is described in detail later.