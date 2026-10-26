#!/usr/bin/env bats
# End-to-end pipeline for ti-sim: toolchain -> venv -> unit tests -> Grossman phantom ->
# layered sphere (label volume -> SimNIBS) -> ernie head (two montages) -> optional
# montage search and MIDA -> summary. Numbers are asserted from the examples' JSON.
#
#   bats --timing --print-output-on-failure tests/pipeline.bats   (or scripts/run-pipeline.sh)
#
# Environment (all optional):
#   TI_SIM_DATA        where scripts/fetch-ernie.sh puts the head (default: ./data)
#   TI_SIM_M2M         SimNIBS m2m_ernie directory (default: $TI_SIM_DATA/m2m_ernie)
#   TI_SIM_FETCH=1     allow the ernie download when missing (413 MB)
#   TI_SIM_WORK        SimNIBS/mesh cache (default: $TI_SIM_DATA/ti-sim-work)
#   TI_SIM_RESULTS     JSON/PNG output directory (default: ./results)
#   TI_SIM_LEADFIELD=1 run the lead-field montage search (20-40 min)
#   TI_SIM_MIDA        your MIDA directory (scripts/check-mida.sh); enables the MIDA stage
#   OMP_NUM_THREADS    CPU threads for SimNIBS (default 4)

bats_require_minimum_version 1.5.0

setup_file() {
    export ROOT
    ROOT="$(cd "$BATS_TEST_DIRNAME/.." && pwd)"
    export TI_SIM_DATA="${TI_SIM_DATA:-$ROOT/data}"
    export TI_SIM_M2M="${TI_SIM_M2M:-$TI_SIM_DATA/m2m_ernie}"
    export TI_SIM_WORK="${TI_SIM_WORK:-$TI_SIM_DATA/ti-sim-work}"
    export TI_SIM_RESULTS="${TI_SIM_RESULTS:-$ROOT/results}"
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-4}"
    export MPLBACKEND=Agg
    export UV="${UV:-uv}"
    mkdir -p "$TI_SIM_RESULTS"
    rm -f "$TI_SIM_RESULTS"/.stage-*
    cd "$ROOT"
}

setup() { cd "$ROOT"; }

py() { "$UV" run --no-sync python "$@"; }

json_get() {
    python3 -I - "$1" "$2" <<'PY'
import json, sys
v = json.load(open(sys.argv[1]))
for k in sys.argv[2].split("."):
    v = v[int(k)] if isinstance(v, list) else v[k]
print(v)
PY
}

assert_range() {
    local v
    v="$(json_get "$1" "$2")"
    python3 -I -c "import sys; v, lo, hi = map(float, sys.argv[1:]); sys.exit(0 if lo <= v <= hi else 1)" \
        "$v" "$3" "$4" || { echo "$2 = $v not in [$3, $4] ($1)"; return 1; }
    echo "$2 = $v in [$3, $4]"
}

need_stage() { [ -f "$TI_SIM_RESULTS/.stage-$1" ] || skip "stage $1 did not pass"; }
mark_stage() { : > "$TI_SIM_RESULTS/.stage-$1"; }

@test "00 toolchain: uv, system python >= 3.12, bats" {
    run -0 "$UV" --version
    echo "$output"
    run -0 /usr/bin/python3 -c 'import sys; assert sys.version_info >= (3, 12), sys.version; print(sys.version)'
    echo "$output"
    bats --version
}

@test "01 uv sync" {
    if [ "${TI_SIM_NO_SYNC:-0}" = 1 ]; then
        echo "sync skipped: TI_SIM_NO_SYNC=1 (container image, read-only venv)"
    else
        run -0 "$UV" sync
        echo "$output" | tail -3
    fi
    mark_stage sync
}

@test "02 system SimNIBS visible in the venv" {
    need_stage sync
    run -0 py -c 'import simnibs, skfem; print("simnibs", simnibs.__version__, "skfem", skfem.__version__)'
    echo "$output"
    mark_stage simnibs
}

@test "10 unit tests (envelope, membrane, phantom, tissues, MIDA mapping, montage search, report)" {
    need_stage sync
    run -0 "$UV" run --no-sync pytest -q -p no:cacheprovider -m "not slow" --deselect tests/test_head.py \
        --junitxml "$TI_SIM_RESULTS/unit.xml"
    echo "$output" | tail -3
    mark_stage unit
}

@test "20 ernie example head present (scripts/fetch-ernie.sh when TI_SIM_FETCH=1)" {
    if [ ! -f "$TI_SIM_M2M/ernie.msh" ]; then
        [ "${TI_SIM_FETCH:-0}" = 1 ] || skip "no $TI_SIM_M2M and TI_SIM_FETCH != 1"
        run -0 scripts/fetch-ernie.sh "$TI_SIM_DATA"
        echo "$output" | tail -2
    fi
    [ -f "$TI_SIM_M2M/ernie.msh" ]
    mark_stage ernie
}

@test "30 Grossman 2017 phantom: geometry of the envelope peaks matches the paper" {
    need_stage unit
    local out="$TI_SIM_RESULTS/grossman_phantom.json"
    TI_SIM_OUT="$out" run -0 py examples/grossman_phantom.py
    echo "$output"
    # layout C, E_AM,y: centred, HWHM (0.42, 0.49) R; E_AM,x at (0, -0.59) R
    assert_range "$out" layouts.C.y.pos_R.0 -0.05 0.05
    assert_range "$out" layouts.C.y.pos_R.1 -0.1 0.1
    assert_range "$out" layouts.C.y.hwhm_R.0 0.35 0.5
    assert_range "$out" layouts.C.y.hwhm_R.1 0.4 0.6
    python3 -I - "$out" <<'PY'
import json, sys
c = json.load(open(sys.argv[1]))["layouts"]["C"]["x"]
# layout C is top/bottom symmetric, so E_AM,x has twin peaks at (0, +-0.59) R; the paper lists one
assert abs(c["pos_R"][0]) < 0.05 and 0.5 < abs(c["pos_R"][1]) < 0.7, c["pos_R"]
PY
    # steering: 1:2.5 and 1:4 move the y-peak by -0.2 R and -0.35 R
    assert_range "$out" layouts.D.y.pos_R.0 -0.3 -0.12
    assert_range "$out" layouts.E.y.pos_R.0 -0.45 -0.25
    # amplitude: a constant factor above the paper's simulated values (open question, see docs)
    assert_range "$out" layouts.C.y.ratio_to_sim 5 9
    assert_range "$out" layouts.E.y.ratio_to_sim 5 9
    mark_stage phantom
}

@test "31 layered sphere: label volume -> SimNIBS mesh -> TI, envelope deeper than |E1| and steered by the ratio" {
    need_stage simnibs
    local out="$TI_SIM_RESULTS/layered_sphere.json"
    TI_SIM_OUT="$out" run -0 py examples/layered_sphere.py
    echo "$output"
    python3 -I - "$out" <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))["ratios"]
for k, row in r.items():  # the envelope peaks deeper than the unmodulated field under the electrodes
    assert row["am_max_depth_mm"] > row["e1_max_depth_mm"] + 10, (k, row)
# the envelope is bounded by the weaker pair's field, so a 1:4 split puts the peak on pair 1's
# side (x < 0) and lowers the envelope at the centre relative to 1:1
assert r["1:4"]["am_peak_x_over_R"] < 0, r["1:4"]
assert r["1:4"]["am_at_centre"] < r["1:1"]["am_at_centre"], (r["1:4"], r["1:1"])
PY
    mark_stage sphere
}

@test "40 ernie: SimNIBS example and Wessel 2023 montages, deep ROIs in the literature range" {
    need_stage ernie
    local out="$TI_SIM_RESULTS/ernie_ti.json"
    TI_SIM_OUT="$out" run -0 py examples/ernie_ti.py
    echo "$output"
    assert_range "$out" "montages.wessel2023_striatum.rois.L putamen.am.median" 0.08 0.6
    assert_range "$out" "montages.wessel2023_striatum.rois.R putamen.am.median" 0.08 0.6
    assert_range "$out" "montages.simnibs_example.rois.GM.am.max" 0.1 2.0
    python3 -I - "$out" <<'PY'
import json, sys
m = json.load(open(sys.argv[1]))["montages"]
for name, row in m.items():
    assert row["focality_0.2"]["am_mm3"] <= row["focality_0.2"]["tacs_mm3"], name
PY
    mark_stage head
}

@test "41 ernie montage search (TI_SIM_LEADFIELD=1): lead field + exhaustive two-pair search" {
    need_stage ernie
    [ "${TI_SIM_LEADFIELD:-0}" = 1 ] || skip "set TI_SIM_LEADFIELD=1 (20-40 min lead field)"
    local out="$TI_SIM_RESULTS/ernie_montage_search.json"
    TI_SIM_OUT="$out" run -0 py examples/ernie_montage_search.py
    echo "$output"
    assert_range "$out" best.0.roi_median 0.05 1.0
    mark_stage search
}

@test "42 ernie placement sensitivity (TI_SIM_SENSITIVITY=1): +-5 mm per electrode, Yamada 2026 ranges" {
    need_stage ernie
    [ "${TI_SIM_SENSITIVITY:-0}" = 1 ] || skip "set TI_SIM_SENSITIVITY=1 (17 SimNIBS runs, ~20 min)"
    local out="$TI_SIM_RESULTS/ernie_placement_sensitivity.json"
    TI_SIM_OUT="$out" run -0 py examples/ernie_placement_sensitivity.py
    echo "$output" | tail -6
    python3 -I - "$out" <<'PY'
import json, sys
s = json.load(open(sys.argv[1]))["summary"]
for e, row in s.items():  # Yamada: 0.4-15 % intensity change and 0.6-6.9 mm hotspot shift for +-10 mm grids
    for roi, pct in row["intensity_sensitivity_pct"].items():
        assert pct < 40, (e, roi, pct)
    assert row["hotspot_shift_max_mm"] < 30, (e, row)
PY
    mark_stage sensitivity
}

@test "50 MIDA (TI_SIM_MIDA): mapping, mesh and, with electrode positions, the 2023 montages" {
    need_stage simnibs
    [ -n "${TI_SIM_MIDA:-}" ] || skip "set TI_SIM_MIDA (scripts/check-mida.sh)"
    run -0 scripts/check-mida.sh "$TI_SIM_MIDA"
    echo "$output"
    local out="$TI_SIM_RESULTS/mida_ti.json"
    TI_SIM_OUT="$out" run -0 py examples/mida_ti.py
    echo "$output"
    assert_range "$out" n_tets 100000 50000000
    mark_stage mida
}

@test "90 summary" {
    for s in sync simnibs unit ernie phantom sphere head search sensitivity mida; do
        if [ -f "$TI_SIM_RESULTS/.stage-$s" ]; then echo "  $s: ok"; else echo "  $s: skipped/failed"; fi
    done
    ls -la "$TI_SIM_RESULTS"/*.json 2>/dev/null || true
}
