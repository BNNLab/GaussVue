import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
import streamlit as st
import streamlit.components.v1 as components
from .gaussian_io import xyz_frame


def render(cube_file, fchk_file, numbers, coords, style, bg):
        st.markdown('**Molecular orbital surfaces** (positive phase: blue; negative phase: pink)')
        st.caption('For accurate Gaussian orbital surfaces, use the Gaussian cubegen utility on the original .fchk. No orbital coefficients are guessed.')
        cube_bytes=cube_file.getvalue() if cube_file else None
        if fchk_file:
            if shutil.which('cubegen'):
                st.success('Gaussian cubegen is available on this machine.')
                orbital=st.text_input('Orbital selection for cubegen','Homo',help='Examples: Homo, Lumo, MO=15, AMO=15, BMO=15')
                if st.button('Generate orbital cube from .fchk'):
                    import re
                    valid=bool(re.fullmatch(r'(Homo|Lumo|MO=\d+|AMO=\d+|BMO=\d+)',orbital,re.IGNORECASE))
                    if not valid:st.error('Enter Homo, Lumo, MO=n, AMO=n or BMO=n.')
                    else:
                        with st.spinner('Running cubegen…'):
                            try:
                                with tempfile.TemporaryDirectory() as d:
                                    fchk=Path(d)/'input.fchk'; out=Path(d)/'orbital.cube'
                                    fchk.write_bytes(fchk_file.getvalue())
                                    proc=subprocess.run(['cubegen','0',orbital,str(fchk),str(out),'45','h'],
                                                        capture_output=True,text=True,timeout=180,cwd=d)
                                    if proc.returncode!=0 or not out.exists():
                                        raise RuntimeError((proc.stderr or proc.stdout)[-1200:])
                                    st.session_state['generated_cube']=out.read_bytes()
                                    st.session_state['generated_cube_name']=orbital
                            except Exception as exc:st.error(f'cubegen failed: {exc}')
            else:
                st.warning('Gaussian cubegen was not found on PATH. Generate a .cube file on a machine with Gaussian installed and upload it here.')
                st.code('cubegen 0 Homo calculation.fchk homo.cube 45 h\ncubegen 0 Lumo calculation.fchk lumo.cube 45 h',language='bash')
        if cube_bytes is None and 'generated_cube' in st.session_state:
            cube_bytes=st.session_state['generated_cube']
        if cube_bytes:
            try:
                cube_text=cube_bytes.decode('utf-8')
                iso=st.slider('Orbital isosurface threshold',0.005,0.15,0.035,0.005)
                if len(coords):
                    orbital_xyz=xyz_frame(numbers,coords[-1],'Final Gaussian geometry')
                else:
                    # 3Dmol can display atomic geometry from cube atom records by itself.
                    orbital_xyz=None
                import json
                payload=json.dumps({'cube':cube_text,'xyz':orbital_xyz,'iso':iso,'bg':bg,'style':style}).replace('</','<\\/')
                html='''<!doctype html><html><head><script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.5.5/3Dmol-min.js"></script><style>html,body,#v{margin:0;width:100%;height:100%;background:#AABDE1}</style></head><body><div id="v"></div><script>
    const p=PAYLOAD;const v=$3Dmol.createViewer('v',{backgroundColor:p.bg,antialias:true});
    if(p.xyz){v.addModel(p.xyz,'xyz')}else{v.addModel(p.cube,'cube')}
    v.setStyle({},{stick:{radius:0.14},sphere:{scale:0.27}});
    v.addVolumetricData(p.cube,'cube',{isoval:p.iso,color:'#38bdf8',opacity:0.78});
    v.addVolumetricData(p.cube,'cube',{isoval:-p.iso,color:'#fb7185',opacity:0.78});
    v.zoomTo();v.render();</script></body></html>'''.replace('PAYLOAD',payload)
                components.html(html,height=600)
                st.download_button('Download cube',cube_bytes,file_name='orbital.cube')
                if len(coords):st.caption('The molecular geometry is from the log file; ensure it corresponds to the uploaded .fchk/.cube.')
            except UnicodeDecodeError:
                st.error('This does not appear to be a text-format Gaussian cube file.')
        else:
            st.info('Upload a .cube file, or upload .fchk and use cubegen to generate one.')
