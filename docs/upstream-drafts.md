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

## 4. SimNIBS: electrode current normalisation in `fem.tdcs` (issue draft, not posted)

Observation (SimNIBS 4.6.0, Fedora build): `fem.tdcs` solves the Dirichlet problem with 0/1 V
on the electrode plug surfaces and scales the solution to the requested current using
`_calc_flux_electrodes`, the flux of the element-wise J through the skin–electrode interface
triangles, averaged over the two electrodes. On meshes with thin single-layer pads this
estimate is biased low, so the returned fields are too high: on a 262k-tet gmsh 4-layer sphere
(10 mm, 2 mm saline pads) the exact Galerkin current of the returned solution is 1.061 mA for a
requested 1 mA; on ernie F3–F4 (19.6 mm, 2 mm pads) 2.27 mA for 2 mA (+13.5 %), with 4 + 2 mm
gel/rubber 2.19 mA (+9.3 %); the SimNIBS example 40 mm electrodes +4.3 / +9.4 %. The logged
"current calibration error" (|flux_a − flux_b| / mean) is 1–3 % in these cases, so it does not
reveal the bias.

Check used: with fields "veEJ", ∫ J·E dV over all tetrahedra equals ΔV · I_h exactly for the
P1 solution, where ΔV is the potential difference between the two plug surfaces (tags 2100+)
and I_h = a(v, χ_electrode) is the consistent discrete flux. An independent DOLFINx solve on
the same mesh and the analytic 4-shell sphere agree with the fields rescaled by I_nominal / I_h
(0.0262 vs 0.0259 V/m per pair at the centre; 279 vs 280 Ω).

`tdcs_leadfield` is not affected: it prescribes the current through `assemble_rhs` and only
logs the flux estimate; re-solving all 75 electrodes of an ernie 10-10 lead field with
Dirichlet conditions and the energy-identity current reproduces SimNIBS's columns to 0.1 %
(median; 2.2 % worst electrode).

Suggestion: in `fem.tdcs`, compute the current as a(v, χ) = (K v)·χ with χ the indicator of
the electrode's Dirichlet nodes (one sparse matrix–vector product with the already assembled
stiffness matrix), or equivalently as the dissipated power divided by ΔV, or use the
right-hand-side formulation the lead field already uses; keep the two-electrode flux
comparison as the mesh-quality warning. Reproducer: ti-sim `scripts/gmsh-sphere-mesh.py` +
`ti_sim.head.discrete_current`.

