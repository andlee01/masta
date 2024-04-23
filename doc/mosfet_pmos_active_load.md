# MOSFET PMOS ACTIVE LOAD

A PMOS current mirror is shown below. The current sources generate a current according to:

$$ i_{REF+} = i_{BASE} + i_{REF} \\ i_{REF-} = i_{BASE} - i_{REF}$$

In the following analysis $i_{REF}$ is swept from $-0.1 \mu A$ to $0.1 \mu A$ and $i_{BASE} = 20 \mu A$.

![title][def0]

[def0]: pmos_active_load.svg

## Large Signal Behaviour
The plot below shows $V_{i_{REF+}}$ and $V_{i_{REF-}}$ as $i_{REF}$ is swept. It can be seen that:
* As M1 is in a diode connected configuration, $V_{i_{REF+}}$ remains almost constant over the range
* $V_{i_{REF-}}$ varies by approximately $ \pm 0.8V$.

![title][def1]

[def1]: mosfet_pmos_active_load_viref.svg

The above behaviour can be understood by considering the plot below. It shows the *i-v* characteristics of M2 at the end points of $i_{REF}$. The green trace corresponds to $i_{REF+} = 20.1 \mu A$ and the red trace corresponds to $i_{REF+} = 19.9 \mu A$. $V_{SG}M_2$ is set by M1 which is dependent on the current $i_{REF+}$. The horizontal lines are plotted at the corresponding end point values of $i_{REF-}$. The intersection of the horizontal line and the *i-v* curve can be used to determine $V_{SD}M_2$.

![title][def2]

[def2]: mosfet_pmos_active_load_iref.svg

## Small Signal Behaviour
Converting the above circuit to a small signal equivalent involves establishing a bias point with $i_{REF+} = i_{REF-} = 20 \mu$. As both M1 and M2 operate in their saturation regions, the equivalent circuit is as shown below, where:

$$ i_{REF+} = i_{REF} \\ i_{REF-} = -i_{REF}$$

![title][def3]

[def3]: pmos_active_load_small_signal.svg

For M1, $i_{REF+}$ is carried almost entirely by current source ${g_m}_{M1}$. This results in almost no current flowing through ${R_O}_{M1}$ and therefore $V_{i_{REF+}}$ is small in magnitude for any given $i_{REF+}$.
For M2, $V_{SG}M_2 = -V_{i_{REF+}}$ and therefore current source ${g_m}_{M2}$ produces a current of opposite polarity to $i_{REF-}$ which must flow through ${R_O}_{M2}$. Therefore, $V_{i_{REF-}}$ must increase in magnitude to counter ${g_m}_{M2}$ with $i_{REF-}$ carried entirely by ${R_O}_{M2}$.