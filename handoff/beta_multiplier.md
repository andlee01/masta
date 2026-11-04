# Beta multiplier handover

Updated: 2026-11-04. This directory is temporary and can be removed when the
handover is no longer needed. On another computer, open this repository in Codex
and ask it to read this file before continuing.

## Progress and committed work

The NMOS-only common-current sweep, Plotly graph, circuit drawing, and Quarto
documentation are implemented and pushed to `origin/master`.

- `85a2b07`: initial NMOS current sweep and Quarto documentation.
- `7ed3777`: plot both NMOS drain voltages and expand the explanation.
- `00ce8f9`: add persistent Python environment dependency list and ignore `.venv/`.

Relevant files:

- `python/circuits/beta_multiplier_nmos_tb.py`: sweep, operating-point validation,
  and interactive Plotly output.
- `python/networks/beta_multiplier_nmos_ntwrk.py`: two NMOS devices, source
  resistor, and two equal ideal current sources; uses the existing MOS models.
- `doc/beta_multiplier/beta_multiplier_nmos.qmd`: explanation and embedded graph.
- `doc/beta_multiplier/beta_multiplier_nmos.html`: standalone Plotly graph,
  committed so the documentation can embed it.
- `doc/beta_multiplier/beta_multiplier_nmos.tex` and `.svg`: circuit drawing.
- `doc/_quarto.yaml`: website navigation and graph resource configuration.
- `python/requirements-beta-multiplier.txt`: Python dependencies.
- `python/tests/test_beta_multiplier_numerics.py`: existing numerical regression
  tests for the models and full beta multiplier solver.

The full beta multiplier also has existing code in
`python/circuits/beta_multiplier_circuit_tb.py`,
`python/networks/beta_multiplier_ntwrk.py`, and `python/subckt/beta_multiplier.py`.
The new sweep does not implement the amplifier-assisted complete circuit.

## Sweep setup and interpretation

M1 is diode-connected with its source grounded. M2 shares M1's gate and has a
6.5 kohm source resistor. W/L is 10/2 for M1 and 40/2 for M2. The model parameters
include KP = 120 uA/V^2, Vth = 0.8 V, and lambda = 0.0125 V^-1.

Two ideal sources impose equal branch currents, swept from 0.1 to 30 uA over
600 points. The plot shows both drain voltages, the conventional M2 saturation
boundary VD2 = VG - Vth, and a 5 V drain-voltage budget. Green shading and
annotations identify the main acceptable current interval.

With ideal square-law devices and no channel-length modulation, the nonzero
equal-current solution is about 19.72 uA. With the existing smooth model, the
main acceptable interval is approximately 18.1521 to 20.4054 uA. A separate
near-threshold acceptable sample is a consequence of the smooth model and the
conventional region classification; it is not the main reference-current range.

This experiment shows what the lower network supports when equal currents are
imposed. It does not establish startup, remove the full circuit's zero-current
solution, or directly measure supply rejection of an amplifier-assisted circuit.

## Python setup on another computer

From the repository root on Linux/WSL, with Python venv support installed:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r python/requirements-beta-multiplier.txt
python python/circuits/beta_multiplier_nmos_tb.py
```

The last command updates the committed standalone graph. To validate without
changing it, use `--output /tmp/beta_multiplier_check.html`.

On this computer the standard venv command failed because system ensurepip was
unavailable. A temporary installation of `virtualenv` was used to create the
project's persistent `.venv` instead. The environment is now functional and
ignored by Git. Do not copy or commit `.venv`; recreate it on each computer.
Dependencies are listed without version pins.

Earlier sessions installed Plotly and control into `/tmp/masta-nmos-deps`.
That temporary directory disappeared. The project `.venv` replaces that setup;
no temporary PYTHONPATH setting is needed now.

## Verification completed

On 2026-11-04, using the new project `.venv`:

```bash
.venv/bin/python python/circuits/beta_multiplier_nmos_tb.py --output /tmp/beta_multiplier_venv_check.html
.venv/bin/python -m unittest discover -s python/tests -p test_beta_multiplier_numerics.py
```

All 600 operating points validated, the main interval was 18.1521 to
20.4054 uA, and all 8 numerical tests passed. Syntax and Git whitespace checks
also passed for the sweep changes. Quarto rendering was not rerun during this
environment setup.

## Where the circuit discussion stopped

The user wants a clear physical/small-signal explanation of why an amplifier-
assisted beta multiplier is less sensitive to VDD. Repeating that equal currents
imply supply independence was unhelpful: equal currents can both change together.

The discussion assumed an amplifier arrangement that equalises the upper PMOS
drain voltages, with matched PMOS devices sharing source and gate connections.
Verify the actual schematic before deriving a complete loop gain or supply
transfer function. No topology-specific amplifier small-signal derivation has
been completed here.

Useful established points:

- Finite PMOS output resistance creates a current disturbance when VSD changes.
  For matched PMOS devices with common sources and gates, the differential
  small-signal current is gop times the difference in their VSD variations.
  Equalising VSD suppresses that differential disturbance; it does not by itself
  prove that the common current is unchanged.
- In the ideal NMOS model, the lower pair and resistor constrain equal currents
  to a particular nonzero operating current. If k = I1/I2, the result is
  I2 = (2/R^2) * (sqrt(k/beta1) - 1/sqrt(beta2))^2, provided the unsquared
  voltage equation admits a positive-current solution. I2 labels the resistor
  branch current, not universally a particular mirror transistor.
- For beta2 = 4*beta1 near k = 1, a small fractional ratio change produces
  approximately twice that fractional change in I2. A fixed mirror-ratio error
  shifts the current; variation of that error with VDD creates line sensitivity.
- Common-current changes require analysis of the entire feedback circuit.
  Equal PMOS VSD does not change the devices' intrinsic output resistance.
  Avoid asserting a specific feedback resistance multiplier without deriving it
  from the actual amplifier topology.

The most recent question focused on the diode-connected NMOS. With its source
fixed, gate tied to drain, and body contribution absent, its small-signal
resistance is 1/(gm + gds), approximately 1/gm. A saturated MOS with fixed VGS has
output resistance 1/gds. Thus, for the same imposed current change, the
diode-connected device's VDS changes much less, by roughly a factor gm*ro.
The diode-connected device is itself normally in saturation; the distinction
is that its gate follows its drain. This is a promising starting point for the
next explanation, tied to the actual circuit and the plotted drain voltages.

## Repository and conversation housekeeping

There are many unrelated untracked scripts, figures, backup files, and Python
caches in the working tree. They were intentionally excluded from these commits;
do not stage the entire working tree indiscriminately.

Original beta multiplier work conversation ID:
`01a16e1a-f97c-7c60-91ea-7031c94fffb4`.
This follow-up conversation ID:
`01a17eaf-6a59-7243-a692-80945e4fdb9c`.

Conversation history is separate from Git. These IDs alone do not transfer the
local sessions to another computer. This document provides the context needed
to start a fresh conversation with a clone of the repository.
