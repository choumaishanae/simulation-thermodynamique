import numpy as np
import pytest
from engine import simulate, binary_diagram
from reference import load_reference

@pytest.mark.parametrize('filename', ['hydrocarbures.dwxml','eau_methanol.dwxml'])
def test_saved_cases_conserve_matter_energy(filename):
    r = load_reference(filename)['streams']['1']['0']
    p = r['properties']
    s = simulate(filename,p['temperature'],p['temperature'],p['pressure'],p['molarflow'],r['composition'])
    assert abs(s['Q']) < 1e-8
    assert max(abs(s['component_error'])) < 1e-6
    assert abs(s['energy_error']) < 1e-6
    assert abs(s['mass_in']-s['mass_liq']-s['mass_vap']) < 1e-6
    assert s['fugacity_error'] < 1e-5

def test_heating_and_flow_scaling():
    a = simulate('eau_methanol.dwxml',298.15,353.15,101325,50/3600,[.6,.4])
    b = simulate('eau_methanol.dwxml',298.15,353.15,101325,100/3600,[.6,.4])
    assert a['Q'] > 0
    assert b['Q'] == pytest.approx(2*a['Q'])
    assert b['beta'] == pytest.approx(a['beta'])
    assert b['V'] == pytest.approx(2*a['V'])

def test_single_vapor():
    r = simulate('hydrocarbures.dwxml',450,450,101325,1,[.1,.2,.3,.4])
    assert r['beta'] == pytest.approx(1)
    assert r['L'] == 0

def test_invalid_composition():
    with pytest.raises(ValueError):
        simulate('eau_methanol.dwxml',298,353,101325,1,[.6,.5])

def test_txy():
    rows, failures = binary_diagram(101325,11)
    assert not failures
    assert len(rows) == 11
    assert rows[0]['bulle_K'] == pytest.approx(337.7,abs=4)
    assert rows[-1]['bulle_K'] == pytest.approx(373.15,abs=4)
    assert all(r['rosee_K'] >= r['bulle_K']-1e-4 for r in rows)
