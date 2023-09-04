# MOSFET PMOS ACTIVE LOAD

A PMOS current mirror is shown below. The current sources generate a current according to:

$$ i_{REF+} = i_{BASE} + i_{REF} \\ i_{REF-} = i_{BASE} - i_{REF}$$

In the following analysis $i_{REF}$ is swept from $-0.1 \mu A$ to $0.1 \mu A$ and $i_{BASE} = 20 \mu A$.

![title][def0]

[def0]: pmos_active_load.svg

The plot below shows $V_{i_{REF+}}$ and $V_{i_{REF-}}$ as $i_{REF}$ is swept. It can be seen that:

* When $i_{REF+} < i_{REF-}$ ($i_{REF} < 0$) then $V_{SD}M_2 > V_{SD}M_1$
* When $i_{REF+} > i_{REF-}$ ($i_{REF} > 0$) then $V_{SD}M_2 < V_{SD}M_1$

The change in $V_{SD}M_2$ is almost linear over the $i_{REF}$ region as both *M1* and *M2* operate in their saturation region.

![title][def1]

[def1]: mosfet_pmos_active_load_viref.svg

The above behaviour can be understood by considering the plot below. It shows the *i-v* characteristics of M2 at the end points of $i_{REF}$. The green trace corresponds to $i_{REF+} = 20.1 \mu A$ and the red trace corresponds to $i_{REF+} = 19.9 \mu A$. $V_{SG}M_1$ (and also $V_{SG}M_2$) increases as $i_{REF+}$ increases. The horizontal lines are plotted at the corresponding end point values of $i_{REF-}$. The intersection of the horizontal line and the *i-v* curve can be used to determine $V_{SD}M_2$. Notice that if $i_{REF-}$ were to decrease much further, then *M2* would enter its triode region and $V_{SD}M_2$ would no longer change in a linear fashion.

![title][def2]

[def2]: mosfet_pmos_active_load_iref.svg