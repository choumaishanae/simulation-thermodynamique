"""Read saved DWSIM values as references, never as calculated results."""
from pathlib import Path
import xml.etree.ElementTree as ET

DATA = Path(__file__).parent / 'data'

def load_reference(name):
    root = ET.parse(DATA / name).getroot()
    tags = {g.findtext('Name'): g.findtext('Tag') for g in root.find('GraphicObjects')}
    streams = {}
    for obj in root.find('SimulationObjects'):
        if not obj.findtext('Type', '').endswith('MaterialStream'):
            continue
        phases = {}
        for phase in obj.findall('Phases/Phase'):
            props = {}
            for node in phase.findall('Properties'):
                for item in node:
                    try:
                        props[item.tag] = float(item.text)
                    except (ValueError, TypeError):
                        pass
            phases[phase.findtext('ID')] = dict(properties=props, composition=[
                float(c.findtext('MoleFraction', '0')) for c in phase.findall('Compounds/Compound')])
        streams[tags[obj.findtext('Name')]] = phases
    compounds = [{n.tag: n.text for n in c if len(n) == 0} for c in root.find('Compounds')]
    package_ids={p.findtext('ID') for p in root.find('PropertyPackages')}
    assigned={p.text for p in root.findall('SimulationObjects/SimulationObject/PropertyPackage') if p.text}
    return dict(compounds=compounds, streams=streams, version=root.findtext('GeneralInfo/BuildVersion'),
                missing_package_ids=sorted(assigned-package_ids))
