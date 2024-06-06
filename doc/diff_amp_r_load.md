# Differential Amplifier with Resistive Load

A differential amplifier with resistive load is shown below. The voltages sources generate a balanced differental signal.

![title][def0]

[def0]: diff_amp_r_load.svg

In the following analysis $i_{REF}$ is kept constant and equal to a nominal value of $200 \mu A$.

## Non-linear Analysis

### Harmonic Distortion

The applied differential voltage contains a single frequency at $f = 1 MHz$

$$ 
\begin{aligned}
V_{LO_{P}} &= A sin(\omega t) \\ 
V_{LO_{N}} &= A sin(\omega t + \pi)
\end{aligned}
$$

where $A = 0.3V$ and $\omega = 2 \pi f$.

![title][def1]

[def1]: diff_amp_r_load_harmonics_fft.svg

The output signal contains odd order harmonics of the input frequency.

![title][def2]

[def2]: diff_amp_r_load_harmonics_trans.svg

### Intermodulation Products

The applied differential voltage contains a two frequencies at $f_0 = 1 MHz$ and $f_1 = 0.8 MHz$

$$ 
\begin{aligned}
V_{LO_{P}} &= A sin(\omega_0 t) + A sin(\omega_1 t)\\ 
V_{LO_{N}} &= A sin(\omega_0 t + \pi) + A sin(\omega_1 t + \pi)
\end{aligned}
$$

where $A = 0.1V$ to $0.3V$, $\omega_0 = 2 \pi f_0$, $\omega_1 = 2 \pi f_1$.

![title][def3]

[def3]: diff_amp_r_load_intermodulation_fft.svg

The output signal contains third-order intermodulation products at $2 \omega_0 - \omega_1$ and $2 \omega_1 - \omega_0$.

![title][def4]

[def4]: diff_amp_r_load_intermodulation_trans.svg

### Non-linear Gain

The applied differential voltage contains a single frequency at $f = 1 MHz$

$$ 
\begin{aligned}
V_{LO_{P}} &= A sin(\omega t) \\ 
V_{LO_{N}} &= A sin(\omega t + \pi)
\end{aligned}
$$

where $A$ is swept from $-0.3V$ to $0.3V$ and $\omega = 2 \pi f$.

![title][def5]

[def5]: diff_amp_r_load_non_linear_gain.svg

The gain is a non-linear function of the input signal leading to the harmonic distoration and intermodulation products witnessed above.