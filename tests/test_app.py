from pathlib import Path
import pytest
from streamlit.testing.v1 import AppTest

APP=str(Path(__file__).parents[1]/'app.py')

@pytest.mark.parametrize('case',[0,1,2,3])
def test_presets(case):
    app=AppTest.from_file(APP,default_timeout=60).run()
    app.sidebar.selectbox[0].select_index(case).run()
    app.sidebar.button[0].click().run()
    assert not app.exception
    assert not app.error
    assert len(app.metric)==4

def test_flow_unit_conversion():
    app=AppTest.from_file(APP,default_timeout=60).run()
    app.sidebar.selectbox[0].select_index(3).run()
    app.sidebar.selectbox[1].select('mol/h').run()
    assert app.sidebar.number_input[3].value==pytest.approx(50.)
    app.sidebar.button[0].click().run()
    q=float(app.metric[3].value)
    app.sidebar.selectbox[1].select('kg/s').run()
    app.sidebar.button[0].click().run()
    assert not app.exception and not app.error
    assert float(app.metric[3].value)==pytest.approx(q)
    app.sidebar.selectbox[1].select('kmol/h').run()
    assert app.sidebar.number_input[3].value==pytest.approx(.05)

def test_invalid_composition_is_reported():
    app=AppTest.from_file(APP,default_timeout=60).run()
    app.sidebar.selectbox[0].select_index(2).run()
    app.sidebar.number_input[4].set_value(.7)
    app.sidebar.button[0].click().run()
    assert not app.exception
    assert len(app.error)==1
    assert len(app.metric)==0

def test_diagram_button():
    app=AppTest.from_file(APP,default_timeout=60).run()
    app.button[0].click().run()
    assert not app.exception and not app.error
    assert len(app.get('plotly_chart'))==1
    assert len(app.dataframe[0].value)==51
