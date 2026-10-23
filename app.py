import json
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from engine import simulate, binary_diagram, make_model, NAMES
from reference import load_reference

st.set_page_config(page_title='Simulations thermodynamiques', page_icon='⚗️', layout='wide')
st.title('Simulations thermodynamiques')
st.caption('Mélanges · Chauffage · Séparation liquide–vapeur')
CASES = {
    'Hydrocarbures — référence DWSIM': ('hydrocarbures.dwxml', None, None),
    'Hydrocarbures — exercice 100 kmol/h, 100 psia, 200 °F': ('hydrocarbures.dwxml', 200*5/9-32*5/9, 100*6894.757293168),
    'Eau / méthanol — référence DWSIM': ('eau_methanol.dwxml', None, None),
    'Eau / méthanol — exercice chauffage à 80 °C': ('eau_methanol.dwxml', 80., None),
}
case = st.sidebar.selectbox('Simulation', list(CASES))
filename, target, exercise_p = CASES[case]
ref = load_reference(filename)
names = [c['Name'] for c in ref['compounds']]
saved = ref['streams']['1']['0']
p = saved['properties']
with st.sidebar:
    st.subheader('Données d’entrée')
    flow_unit = st.selectbox('Unité du débit', ['kmol/h', 'mol/h', 'mol/s', 'kg/s'],key='unit_'+case)
    flow_key, previous_key = 'flow_'+case, 'previous_unit_'+case
    default_flow = 100. if '100 kmol' in case else (0.05 if '80 °C' in case else p['molarflow']*3.6)
    if flow_key not in st.session_state:
        st.session_state[flow_key] = default_flow
        st.session_state[previous_key] = 'kmol/h'
    if flow_unit != st.session_state[previous_key]:
        old = st.session_state[previous_key]
        mw = sum(float(c['Molar_Weight'])*st.session_state.get('z_'+case+'_'+n,z)
                 for c,n,z in zip(ref['compounds'],names,saved['composition']))
        if mw<=0:
            mw=sum(float(c['Molar_Weight'])*z for c,z in zip(ref['compounds'],saved['composition']))
        to_mols = {'kmol/h':1/3.6, 'mol/h':1/3600, 'mol/s':1., 'kg/s':1000/mw}
        st.session_state[flow_key] *= to_mols[old]/to_mols[flow_unit]
        st.session_state[previous_key] = flow_unit
    with st.form('inputs_'+case):
        tin = st.number_input('Température entrée (°C)', min_value=-200., max_value=500., value=p['temperature']-273.15)
        tout = st.number_input('Température après chauffage / séparation (°C)', min_value=-200., max_value=500., value=target if target is not None else p['temperature']-273.15)
        pressure = st.number_input('Pression absolue (bar)', min_value=0.001, max_value=500., value=(exercise_p or p['pressure'])/1e5, format='%.6f')
        flow = st.number_input('Débit ('+flow_unit+')', min_value=1e-9, key=flow_key, format='%.9f')
        st.write('Composition molaire')
        zs = [st.number_input(NAMES[n], min_value=0., max_value=1., value=z, step=0.01, format='%.6f',key='z_'+case+'_'+n)
              for n,z in zip(names, saved['composition'])]
        st.caption(f'Somme des fractions : {sum(zs):.6f}')
        submitted = st.form_submit_button('Calculer', type='primary', width='stretch')

st.caption('Moteur Python PRSV2-VL, vérifié contre DWSIM 10.2.11 : équilibre, enthalpies et diagramme binaire.')
st.caption('Entrée → Chauffage à température imposée → Séparateur à la même température et pression. '
           'Pas de perte de pression ; le séparateur est adiabatique après le chauffage.')
tab_calc, tab_plot, tab_ref = st.tabs(['Résultats et bilans', 'Diagramme T-x-y', 'Références et validation'])

if submitted:
    st.session_state['calculation'] = None
    try:
        _, constants, _ = make_model(filename)
        F = flow/3.6 if flow_unit == 'kmol/h' else flow/3600 if flow_unit == 'mol/h' else flow
        if flow_unit == 'kg/s':
            F = flow*1000/sum(z*m for z,m in zip(zs, constants.MWs))
        with st.spinner('Calcul de l’équilibre liquide–vapeur…'):
            result = simulate(filename, tin+273.15, tout+273.15, pressure*1e5, F, zs)
        st.session_state['calculation'] = (case, result, tin, tout, pressure, F, zs)
    except Exception as exc:
        st.error(f'Calcul impossible : {exc}')

with tab_calc:
    calc = st.session_state.get('calculation')
    if calc and calc[0] == case:
        _, r, ti, to, pr, F, composition = calc
        st.caption(f'Calcul effectué : entrée {ti:g} °C → sortie {to:g} °C ; {pr:g} bar absolus.')
        a,b,c = st.columns(3)
        a.metric('Fraction vapeur molaire', f'{r["beta"]:.6f}')
        b.metric('Liquide (kmol/h)', f'{r["L"]*3.6:.6g}')
        c.metric('Vapeur (kmol/h)', f'{r["V"]*3.6:.6g}')
        rows = []
        for label, f, mass, h in [('Entrée', F, r['mass_in'],r['hin']),
                                  ('Liquide',r['L'],r['mass_liq'],r['hl']),
                                  ('Vapeur',r['V'],r['mass_vap'],r['hv'])]:
            rows.append({'Flux':label,'T (°C)':ti if label=='Entrée' else to,'P (bar abs)':pr,
                         'Débit (kmol/h)':f*3.6,'Débit (kg/s)':mass,'Flux enthalpique (kW)':h})
        st.dataframe(pd.DataFrame(rows), hide_index=True, width='stretch')
        comps = pd.DataFrame({'Composé':[NAMES[n] for n in names], 'z entrée':composition,
                              'x liquide':r['x'] if r['L']>0 else [None]*len(names),
                              'y vapeur':r['y'] if r['V']>0 else [None]*len(names)})
        st.dataframe(comps, hide_index=True, width='stretch')
        st.caption('Une phase absente n’a pas de composition définie (cellules vides).')
        st.metric('Puissance de chauffage Q (kW)', f'{r["Q"]:.8g}')
        st.write('Convention : Q > 0 pour la chaleur reçue. H entrée + Q = H liquide + H vapeur.')
        st.write(f'Résidu énergétique : {r["energy_error"]:.3e} kW · '
                 f'Résidu matière maximal : {max(abs(r["component_error"])):.3e} mol/s · '
                 f'Résidu ln(f liquide / f vapeur) : {r["fugacity_error"]:.3e}')
        st.caption('Enthalpies : corrélations Cp du fichier original + départ PRSV2-VL de DWSIM 10.2.11. '
                   'Référence gaz parfait à 298,15 K. Les valeurs non arrondies servent aux bilans.')
    else:
        st.write('Renseignez les données à gauche puis cliquez sur Calculer.')

@st.cache_data(show_spinner=False,max_entries=16)
def diagram_cached(P):
    return binary_diagram(P)

with tab_plot:
    st.write('Eau / méthanol à pression constante')
    dp = st.number_input('Pression du diagramme (bar absolus)', min_value=0.01, max_value=10., value=1.01325, format='%.5f')
    if st.button('Tracer le diagramme T-x-y'):
        with st.spinner('Calcul des points de bulle et de rosée…'):
            rows, failures = diagram_cached(dp*1e5)
        if rows:
            df = pd.DataFrame(rows)
            fig = go.Figure()
            fig.add_scatter(x=df.fraction_eau, y=df.bulle_K, name='Bulle (x eau)', mode='lines')
            fig.add_scatter(x=df.fraction_eau, y=df.rosee_K, name='Rosée (y eau)', mode='lines')
            fig.update_layout(xaxis_title='Fraction molaire eau : x (liquide) / y (vapeur)', yaxis_title='Température (K)', template='plotly_white')
            st.plotly_chart(fig, width='stretch')
            st.dataframe(df, hide_index=True)
        if failures:
            st.warning(f'{len(failures)} points non convergés ont été exclus. Le diagramme est incomplet.')
        st.caption('Équilibre vapeur–liquide uniquement ; pas de calcul LLE/SLE. Extrémités calculées par saturation des composés purs.')

with tab_ref:
    st.write(f'Fichier original : {filename} — DWSIM {ref["version"]}')
    st.write('Références indépendantes calculées avec le moteur officiel DWSIM 10.2.11. '
             'Les résultats de l’application sont recalculés en Python ; ils ne sont pas lus dans ces références.')
    if ref['missing_package_ids']:
        st.caption('Le fichier original contient des liens vers un paquet thermodynamique absent. '
                   'La validation utilise le seul paquet PRSV2-VL effectivement enregistré. '
                   'Les anciennes enthalpies sauvegardées peuvent donc différer du recalcul.')
    st.caption('L’énoncé eau/méthanol demande 50 mol/h ; le fichier contient 0,01389 mol/s, '
               'soit 50,004 mol/h. Les compositions 0,6 / 0,4 sont molaires dans les fichiers. '
               'Les pourcentages manuscrits de l’exercice hydrocarbures sont traités comme molaires.')
    if st.button('Comparer avec DWSIM 10.2.11'):
        comparisons = []
        records=json.loads((Path(__file__).parent/'data'/'dwsim_validation.json').read_text(encoding='utf-8'))['records']
        for r in records:
            if r.get('kind')!='flash_tight':
                continue
            try:
                result = simulate(r['file'],r['T'],r['T'],r['P'],1.,r['z'])
                comparisons.append({'Cas':r['file'],'T (°C)':r['T']-273.15,'P (bar)':r['P']/1e5,
                    'β DWSIM':r['beta'],'β Python':result['beta'],'Écart β':result['beta']-r['beta'],
                    'Écart enthalpie (J/mol)':result['hin']*1000-r['H']})
            except Exception as exc:
                st.error(f'{r["file"]} : {exc}')
        st.dataframe(pd.DataFrame(comparisons),hide_index=True)
        st.caption('Référence DWSIM : tolérances PT interne/externe à 10⁻¹². '
                   'Les réglages d’origine, moins stricts, expliquent de petits écarts avec les captures.')
    with st.expander('Modèle et limites de comparaison'):
        st.write('PRSV2-VL de DWSIM 10.2.11 : règle de mélange Van Laar, interactions eau/méthanol '
                 'directionnelles −0,0968 / −0,0839, constantes et corrélations Cp des fichiers originaux. '
                 'Le calcul énergétique reprend H_PR_MIX, effectivement utilisé par ce paquet, malgré '
                 'le libellé LeeKesler hérité dans les réglages XML. '
                 'La validation couvre les deux mélanges des exercices : 14 états PT, 16 propriétés de phase '
                 'et 16 températures de bulle/rosée. Une seule phase liquide et une phase vapeur ; '
                 'pas de séparation liquide–liquide ni de phase solide. Sources et licence GPL dans le projet.')


