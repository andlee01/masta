# COMMON SOURCE AMPLIFIER

A common source amplifier is shown below.

![title][def0]

[def0]: common_source_amp.svg

## Large Signal Behaviour
The plot below shows $V_{out}$ and $i_{DS}$ as $V_{in}$ is swept from $0V$ to $5V$ where $R_S = 0 \Omega$. 

![title][def1]

[def1]: common_source_amp_transfer.svg

The plot below shows $V_{out}$ and $i_{DS}$ as $V_{in}$ is swept from $0V$ to $5V$ where $R_S = 50k \Omega$. 

![title][def2]

[def2]: common_source_amp_degen_transfer.svg

It can be seen that the source degeneration amplifer ($R_S = 50k \Omega$) exhibits a more linear response than the amplifer without source degeneration, allowing a wider range of $V_{in}$ values before $M1$ triodes, at the expense of $V_{out}$ swing, and therefore gain.

This can be understood with the aid of the load line graph.

The plot below shows $R_S = 0 \Omega$.


![title][def3]

[def3]: common_source_amp_loadline.svg

The plot below shows $R_S = 50k \Omega$. The load line is plotted for $R_D$.

![title][def4]

[def4]: common_source_amp_degen_loadline.svg

Comparing the two load lines, it is clear that the degeneration resistor reduces variations in $V_{GS}$ for a given $V_{in}$ swing. With the degeneration resistor in place, any increase in $i_{DS}$, resulting from an increase in $V_{in}$, also results in an increase in voltage drop across $R_S$, which in turn limits the effective increase of $V_{GS}$. This self correcting mechanism works as negative feedback and serves to limit changes in $V_{GS}$ for a given change in $V_{in}$.

The negative feedback serves to improve amplifier linearity. As $i_{DS} \propto (V_{GS} - V_{th})^2$ in the saturation region, this quadratic relationship is inherently non-linear and produces a progressively less linear response as $V_{GS}$ increases in magnitide.

## Small Signal Behaviour

![title][def5]

[def5]: common_source_amp_small_signal.svg

![title][def6]

[def6]: common_source_amp_gm.svg

Note that even though $i_{DS}$ is nonlinear in the saturation region (particularly noticable on $g_{m_0}$) the transconductance appears linear. This is due to the square-root dependence of $g_m$ on $i_{DS}$, which tends to smooth out the nonlinearity of $i_{DS}$.