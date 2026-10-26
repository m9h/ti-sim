import pytest

from ti_sim import tissues


def test_tables_cover_the_simnibs_tags():
    for name, t in tissues.TABLES.items():
        assert set(t) == set(range(1, 13)), name
        assert all(v > 0 for v in t.values())
    assert tissues.SIMNIBS[2] == 0.275 and tissues.RAMPERSAD[2] == 0.276


def test_table_from_env(monkeypatch):
    monkeypatch.delenv("TI_SIM_COND", raising=False)
    assert tissues.table_from_env() == ("simnibs", tissues.SIMNIBS)
    monkeypatch.setenv("TI_SIM_COND", "rampersad,4=0.0085")
    label, t = tissues.table_from_env()
    assert label == "rampersad,4=0.0085" and t[4] == 0.0085 and t[2] == 0.276
    monkeypatch.setenv("TI_SIM_COND", "nope")
    with pytest.raises(KeyError):
        tissues.table_from_env()


def test_apply_to_tdcslist():
    sim_struct = pytest.importorskip("simnibs.sim_struct")
    tdcs = sim_struct.TDCSLIST()
    tissues.apply_to_tdcslist(tdcs, tissues.RAMPERSAD)
    assert tdcs.cond[1].value == 0.276 and tdcs.cond[1].name == "GM"
