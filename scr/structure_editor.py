import numpy as np
import streamlit as st
from periodictable import elements
from .geometry import inferred_bonds, edit_geometry, smiles_to_3d
from .gaussian_io import xyz_frame, pdb_block, gaussian_input
from .editor_component import PICKER
from .fragments import FRAGMENTS, attach_fragment


def render(numbers, coords, data, bg):
    st.subheader('Visual structure editor')
    st.caption('Build a 3D molecule from SMILES or import a Gaussian geometry. Click atoms in the viewer to select them; selected atoms are highlighted in yellow and numbered in selection order.')
    a,b = st.columns([3,2])
    with a:
        smiles = st.text_input('SMILES', value='CCO', help='Paste a SMILES string; explicit hydrogens are generated automatically.')
        if st.button('Generate 3D from SMILES',type='primary'):
            try:
                en,ex,eb,bt = smiles_to_3d(smiles.strip())
                st.session_state['editor_nums'],st.session_state['editor_xyz'],st.session_state['editor_bonds'] = en,ex,eb
                st.session_state['editor_bond_types'] = bt
                st.session_state.pop('editor_live_signature', None)
                st.session_state['editor_history'] = []
                st.session_state['editor_revision'] = st.session_state.get('editor_revision',0)+1
                st.session_state['editor_origin'] = 'smiles'
                st.success('3D structure generated with RDKit ETKDG and MMFF/UFF where supported.')
            except Exception as err:
                st.error(f'SMILES generation failed: {err}')
    with b:
        if len(coords):
            source_step = st.number_input('Gaussian geometry step',min_value=1,max_value=len(coords),value=len(coords))
            if st.button('Import Gaussian geometry'):
                st.session_state['editor_nums'] = list(map(int,numbers))
                st.session_state['editor_xyz'] = np.array(coords[source_step-1],float).copy()
                st.session_state['editor_bonds'] = inferred_bonds(st.session_state['editor_nums'],st.session_state['editor_xyz'])
                st.session_state['editor_bond_types'] = {tuple(sorted(e)):'Single' for e in st.session_state['editor_bonds']}
                st.session_state.pop('editor_live_signature', None)
                st.session_state['editor_history'] = []
                st.session_state['editor_revision'] = st.session_state.get('editor_revision',0)+1
                st.session_state['editor_origin'] = 'gaussian'
        else:
            st.caption('Upload a Gaussian log to import an optimisation geometry, or generate a molecule from SMILES.')

    if 'editor_xyz' not in st.session_state:
        st.info('Generate a structure from SMILES to begin.')
    else:
        en = st.session_state['editor_nums']
        ex = st.session_state['editor_xyz']
        eb = st.session_state['editor_bonds']
        st.info(f'{len(en)} atoms · {len(eb)} bonds · Click to select atoms, click again to deselect.')
        left,right = st.columns([3,2])
        with right:
            operation = st.selectbox('Editing operation', ['Bond length','Bond angle','Dihedral angle','Add atom','Add fragment','Remove atom','Add bond','Bond type'])
            required = {'Bond length':2,'Bond angle':3,'Dihedral angle':4,'Add atom':1,'Add fragment':1,'Remove atom':1,'Add bond':2,'Bond type':2}[operation]
            st.caption(f'Select {required} atom(s) in order: ' + {
                'Bond length':'A → B (move B side)',
                'Bond angle':'A → B (vertex) → C (move C side)',
                'Dihedral angle':'A → B → C → D (rotate D side)',
                'Add atom':'attachment atom', 'Add fragment':'attachment atom', 'Remove atom':'atom to delete', 'Bond type':'two atoms to connect or disconnect', 'Add bond':'two atoms to connect'}[operation])
            symbol=None; value=None
            specs = {
                'Add atom': ('Bond length (Å)', 0.3, 5.0, 1.090, 0.001),
                'Bond length': ('Target length (Å)', 0.3, 10.0, 1.500, 0.001),
                'Bond angle': ('Target angle (°)', 0.1, 179.9, 109.5, 0.1),
                'Dihedral angle': ('Target dihedral (°)', -180.0, 180.0, 180.0, 0.1),
            }
            bond_type = st.selectbox('Bond type', ['Single','Double','Aromatic','Triple','H-bond','No bond']) if operation in ('Bond type','Add bond') else None
            if operation == 'Add fragment':
                st.markdown('**Fragment library**')
                category = st.selectbox('Category', list(FRAGMENTS), key='fragment_category')
                fragment_name = st.selectbox('Fragment', list(FRAGMENTS[category]), key='fragment_name')
                fragment_smiles, _, default_distance = FRAGMENTS[category][fragment_name]
                st.caption(f'Fragment SMILES: `{fragment_smiles}`')
                fragment_distance = st.number_input('Attachment distance (Å)', min_value=0.5,
                    max_value=4.0, value=float(default_distance), step=0.001, format='%.3f')
                st.caption('Select one atom in the 3D viewer. An attached hydrogen is replaced if available.')
            if operation == 'Add atom':
                symbol=st.text_input('Element symbol','H').strip().capitalize()
            if operation in specs:
                label,low,high,initial,increment = specs[operation]
                # Streamlit deletes widget keys when a widget is not rendered.
                # Persist values in non-widget keys so switching operation never
                # resets a previously entered target to its default.
                stored_key=f'editor_saved_value_{operation}'
                number_key=f'editor_number_{operation}'
                slider_key=f'editor_slider_{operation}'
                if stored_key not in st.session_state:
                    st.session_state[stored_key]=initial
                if number_key not in st.session_state:
                    st.session_state[number_key]=st.session_state[stored_key]
                if slider_key not in st.session_state:
                    st.session_state[slider_key]=st.session_state[stored_key]
                def from_number(n=number_key,sl=slider_key,store=stored_key):
                    st.session_state[store]=st.session_state[n]
                    st.session_state[sl]=st.session_state[n]
                def from_slider(n=number_key,sl=slider_key,store=stored_key):
                    st.session_state[store]=st.session_state[sl]
                    st.session_state[n]=st.session_state[sl]
                display_format='%.3f' if operation in ('Bond length','Add atom') else '%.1f'
                value=st.number_input(label,min_value=low,max_value=high,step=increment,
                                      format=display_format,key=number_key,on_change=from_number)
                st.slider('Adjust '+label,min_value=low,max_value=high,step=increment,
                          format=display_format,key=slider_key,on_change=from_slider)
                value=st.session_state[stored_key]
        with left:
            revision=st.session_state.get('editor_revision',0)
            atoms=[{'elem':elements[int(z)].symbol,'x':float(x[0]),'y':float(x[1]),'z':float(x[2])}
                   for z,x in zip(en,ex)]
            selection=PICKER(
                data={'atoms':atoms,'bonds':eb,'bg':bg,'maxPick':required,
                      'operation':operation,'value':value,'cameraKey':'gaussian_editor_camera',
                      'bondTypes':{'-'.join(map(str,k)):v for k,v in st.session_state.get('editor_bond_types',{}).items()}},
                default={'selected':[]},
                on_selected_change=lambda: None,
                key=f'picker_{revision}_{operation}',
                height=580,
            )
        picked = list(selection.selected or [])
        with right:
            st.write('**Selected:**', ' → '.join(f'{elements[int(en[i])].symbol}{i+1}' for i in picked) if picked else 'None')
            if st.button('Clear selection'):
                st.session_state['editor_revision'] = st.session_state.get('editor_revision',0)+1
                st.session_state.pop('editor_live_signature',None)
                st.rerun()

            # Changing a numeric target applies immediately. Each selection uses
            # a fixed baseline to avoid cumulative error while moving a slider.
            numeric_ops = ('Bond length','Bond angle','Dihedral angle')
            signature = (operation, tuple(picked), st.session_state.get('editor_revision',0))
            if st.session_state.get('editor_live_signature') != signature:
                st.session_state['editor_live_signature'] = signature
                st.session_state['editor_live_base'] = np.asarray(ex,float).copy()
                st.session_state['editor_live_last_value'] = value
                st.session_state['editor_live_undo_saved'] = False
            elif operation in numeric_ops and len(picked)==required and value != st.session_state.get('editor_live_last_value'):
                try:
                    baseline = st.session_state['editor_live_base']
                    covalent_edges=[edge for edge in eb if st.session_state.get('editor_bond_types',{}).get(tuple(sorted(edge)))!='H-bond']
                    new_nums,new_xyz,new_bonds = edit_geometry(en,baseline,covalent_edges,operation,picked,value)
                    if not st.session_state.get('editor_live_undo_saved'):
                        st.session_state.setdefault('editor_history',[]).append((list(en),baseline.copy(),list(eb),dict(st.session_state.get('editor_bond_types',{}))))
                        st.session_state['editor_live_undo_saved']=True
                    st.session_state['editor_xyz']=new_xyz
                    st.session_state['editor_live_last_value']=value
                    st.rerun()
                except Exception as err:
                    st.error(f'Live edit could not be applied: {err}')
                    st.session_state['editor_live_last_value']=value

            if operation == 'Add fragment':
                if st.button('Attach fragment', type='primary', disabled=len(picked)!=1):
                    try:
                        result = attach_fragment(en, ex, eb,
                            st.session_state.get('editor_bond_types', {}),
                            picked[0], fragment_name, fragment_distance)
                        st.session_state.setdefault('editor_history', []).append(
                            (list(en), ex.copy(), list(eb), dict(st.session_state.get('editor_bond_types', {}))))
                        (st.session_state['editor_nums'], st.session_state['editor_xyz'],
                         st.session_state['editor_bonds'], st.session_state['editor_bond_types']) = result
                        st.session_state['editor_revision'] = st.session_state.get('editor_revision', 0) + 1
                        st.session_state.pop('editor_live_signature', None)
                        st.rerun()
                    except Exception as exc:
                        st.error(f'Fragment attachment failed: {exc}')

            if operation in ('Bond type','Add bond'):
                if len(picked)==2:
                    pair=tuple(sorted(picked))
                    current=st.session_state.get('editor_bond_types',{}).get(pair,'No bond')
                    st.caption(f'Current bond: {current}')
                if st.button('Add bond' if operation=='Add bond' else 'Change bond type',type='primary',disabled=len(picked)!=2):
                    pair=tuple(sorted(picked))
                    if operation=='Add bond' and pair in [tuple(sorted(e)) for e in eb]:
                        st.warning('A bond already exists between these atoms. Use Bond type to modify it.')
                    elif operation=='Add bond' and bond_type=='No bond':
                        st.warning('Select a bond type other than No bond.')
                    else:
                        st.session_state.setdefault('editor_history',[]).append((list(en),ex.copy(),list(eb),dict(st.session_state.get('editor_bond_types',{}))))
                        if bond_type=='No bond':
                            st.session_state['editor_bonds']=[e for e in eb if tuple(sorted(e))!=pair]
                            st.session_state['editor_bond_types'].pop(pair,None)
                        else:
                            if pair not in [tuple(sorted(e)) for e in eb]:
                                st.session_state['editor_bonds'].append(pair)
                            st.session_state['editor_bond_types'][pair]=bond_type
                        st.session_state['editor_revision']=st.session_state.get('editor_revision',0)+1
                        st.rerun()
            elif operation in ('Add atom','Remove atom'):
                if st.button('Add atom' if operation=='Add atom' else 'Remove atom',type='primary',disabled=len(picked)!=required):
                    try:
                        if operation=='Remove atom' and len(en)<=1:
                            raise ValueError('Cannot remove the last atom.')
                        result=edit_geometry(en,ex,eb,operation,picked,value,symbol)
                        st.session_state.setdefault('editor_history',[]).append((list(en),ex.copy(),list(eb),dict(st.session_state.get('editor_bond_types',{}))))
                        st.session_state['editor_nums'],st.session_state['editor_xyz'],st.session_state['editor_bonds']=result
                        if operation=='Remove atom':
                            deleted=picked[0]
                            st.session_state['editor_bond_types']={(i-(i>deleted),j-(j>deleted)):kind for (i,j),kind in st.session_state['editor_bond_types'].items() if i!=deleted and j!=deleted}
                        else:
                            st.session_state['editor_bond_types'][tuple(sorted((picked[0],len(en))))]='Single'
                        st.session_state['editor_revision']=st.session_state.get('editor_revision',0)+1
                        st.rerun()
                    except Exception as err:
                        st.error(f'Edit not applied: {err}')
            else:
                st.caption('Numeric edits are applied automatically when the value changes. No Apply button is needed.')
            if st.button('Undo last edit',disabled=not st.session_state.get('editor_history')):
                old=st.session_state['editor_history'].pop()
                st.session_state['editor_nums'],st.session_state['editor_xyz'],st.session_state['editor_bonds']=old[:3]
                if len(old)>3: st.session_state['editor_bond_types']=old[3]
                st.session_state['editor_revision']=st.session_state.get('editor_revision',0)+1
                st.session_state.pop('editor_live_signature',None)
                st.rerun()
        with st.expander('Atom numbers and coordinates'):
            import pandas as pd
            st.dataframe(pd.DataFrame([{'Index':i+1,'Element':elements[int(z)].symbol,
                'X (Å)':float(x[0]),'Y (Å)':float(x[1]),'Z (Å)':float(x[2])}
                for i,(z,x) in enumerate(zip(en,ex))]),hide_index=True,use_container_width=True)
        st.divider()
        st.subheader('Export structure')
        en,ex,eb=st.session_state['editor_nums'],st.session_state['editor_xyz'],st.session_state['editor_bonds']
        c1,c2=st.columns(2)
        with c1:
            command=st.text_input('Gaussian route commands','Opt Freq')
            functional=st.text_input('Functional / method','B3LYP')
            basis=st.text_input('Basis set','6-31G(d)')
            solvation=st.selectbox('Solvation model',['None','SMD','PCM','CPCM'])
            solvent=st.text_input('Solvent','Water',disabled=solvation=='None')
            default_charge=int(getattr(data,'charge',0)) if st.session_state.get('editor_origin')=='gaussian' and data is not None else 0
            default_mult=int(getattr(data,'mult',1)) if st.session_state.get('editor_origin')=='gaussian' and data is not None else 1
            charge=st.number_input('Charge',value=default_charge,step=1)
            multiplicity=st.number_input('Multiplicity',min_value=1,value=default_mult,step=1)
        with c2:
            title=st.text_input('Title','Gaussian Explorer edited structure')
            chk=st.text_input('Checkpoint filename','edited.chk')
            memory=st.text_input('Memory','8GB')
            nproc=st.number_input('CPU cores',min_value=1,max_value=256,value=8)
            st.warning('Verify chemical validity, formal charge, spin multiplicity and generated geometry before Gaussian calculations. Geometry edits do not automatically minimise the structure.')
        try:
            gjf=gaussian_input(en,ex,command,functional,basis,solvation,solvent,charge,multiplicity,title,nproc,memory,chk)
            with st.expander('Gaussian input preview'):
                st.code(gjf,language='text')
            x,y,z=st.columns(3)
            x.download_button('Download XYZ',xyz_frame(en,ex,'Edited structure'),'edited.xyz')
            y.download_button('Download PDB',pdb_block(en,ex,eb,st.session_state.get('editor_bond_types',{})),'edited.pdb')
            z.download_button('Download Gaussian GJF',gjf,'edited.gjf')
        except Exception as err:
            st.error(f'Export error: {err}')
