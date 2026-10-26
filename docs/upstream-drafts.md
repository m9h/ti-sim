# Upstream drafts

**Status:** drafts, 2026-10-25. **Nothing has been posted.** Morgan posts, or decides not to.

## 1. Grossman et al. 2017: the phantom amplitude definition (question to the authors)

**To:** Nir Grossman (Imperial) / Ed Boyden (MIT), or as a question on the Cell supplement.

We reproduced the figure 2 saline-phantom table from the stated geometry (50 mm × 10 mm
cylinder, 1 mm wires on the rim at the trapezoid corners, 1 kHz / 1.02 kHz, 1 mA per
pair, σ = 0.333 S/m as in figure S2) with a 2-D finite element model checked against the
analytic line-source field. Every peak position and half-width in the table is
reproduced (layout C: E_AM,y at (0, 0) with HWHM (0.42, 0.50) R vs (0.42, 0.49); D and E
steer to x = -0.23 and -0.34 R vs -0.2 and -0.35 R). The amplitudes, however, come out a
constant factor ~7 above the table's simulated values in layouts C-E (1.51 vs 0.21 V/m
in C). Could you tell us (a) whether the tabulated |E_AM| is the peak-to-trough
modulation (equation 1) or half of it, (b) whether "1 mA" is amplitude or peak-to-peak,
and (c) the conductivity used for the simulated phantom (0.333 S/m as in S2, or that of
the physical saline)? Code: github.com/m9h/ti-sim, `examples/grossman_phantom.py`.

## 2. SimNIBS: `get_maxTI` docstring and a tetrahedral TI example

Small: `TI_utils.get_maxTI` returns the peak-to-trough modulation amplitude (2|E₂| for
parallel equal fields), which is what Grossman's equation 1 gives, but the docstring says
"maximal modulation amplitude", which some readers take as half of that; a sentence
would settle it. Also an example computing TI on tetrahedra in deep ROIs (both shipped
examples visualise on the GM surface), since the deep targets are the point of TI.
Offered as a PR if wanted; our `examples/ernie_ti.py` does it.

## 3. TI-Toolbox: cross-check offer

Offer to cross-check their exhaustive search against `ti_sim.montage` on a shared charm
head and the same ROI, and to share the MIDA-through-create_mesh route if they want a
MIDA option (they cannot ship MIDA either).
