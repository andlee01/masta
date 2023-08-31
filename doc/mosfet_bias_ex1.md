# MOSFET BIAS EX1

The diode connected MOSFET (M1) conducts a drain current $i_{REF}$, assuming that $V_{CC}$ is large enough to ensure that M1 operates in the saturation region.

In this configuration $V_{DS} = V_{GS}$. The biasing current $i_{REF}$ establishes a bias voltage $V_{BIAS}$ between the gate-source of M1.

![title][def0]

[def0]: mosfet_bias_ex1.svg

The figure below shows the *i-v* characteristics of M1 (labelled $i_{DS}M_1$) as $i_{REF}$ is increased from 1 $\mu$ A to 150 $\mu$ A, superimposed against the *i-v* characteristics of a matched (equally sized) MOSFET. With reference to the figure below, it can be observed that when the drain currents of M1 and the matched MOSFET are equal, the gate-source voltages are also equal (since $V_{DS} = V_{GS} = V_{BIAS}$ of M1). 

![title][plt0]

[plt0]: plt_mosfet_bias_ex1.svg

$V_{BIAS}$ is often used to control the drain current of a MOSFET in a current mirror, such as that shown in the figure below. If channel length modulation is neglected, then the drain currents flowing in M1 and M2 are equal (assuming both M1 and M2 are operating in the satuation region).

![title][def1]

[def1]: mosfet_bias_ex1_mirror.svg