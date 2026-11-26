import numpy as np
import streamlit as st
from .gaussian_io import xyz_frame
from .viewer import viewer


def render(numbers, coords, energies, log_file, style, bg, spin):
    if not len(coords):
        st.info('No data to visualise yet!')
    else:
        left,right=st.columns([3,2])
        with right:
            log_identity = (log_file.name, hash(log_file.getvalue()))
            if st.session_state.get('active_log_identity') != log_identity:
                st.session_state['active_log_identity'] = log_identity
                st.session_state['geom_step'] = len(coords)
                st.session_state['last_plot_selection'] = None
            if 'pending_geom_step' in st.session_state:
                st.session_state['geom_step'] = st.session_state.pop('pending_geom_step')
            step=st.slider('Geometry step',1,len(coords),len(coords),key='geom_step')
            playback=st.toggle('Animate optimisation trajectory',False,disabled=len(coords)<2)
            st.caption('Geometry coordinates are taken directly from Gaussian output (Å).')
            if len(energies):
                import plotly.graph_objects as go
                st.subheader('SCF energy history')
                e = energies / 27.211386245988
                finite = e[np.isfinite(e)]
                if len(finite):
                    y_min, y_max = float(np.min(finite)), float(np.max(finite))
                    padding = max((y_max - y_min) * 0.05, 1e-6)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(
                        x=np.arange(1, len(e) + 1),
                        y=e,
                        mode='lines+markers',
                        line=dict(color='#38bdf8', width=2),
                        marker=dict(size=5),
                        name='SCF Energy'
                    ))
                    fig.update_layout(
                        xaxis_title='SCF Step',
                        yaxis_title='Energy (Hartree)',
                        yaxis=dict(
                            range=[y_min - padding, y_max + padding],
                            tickformat='.6f'
                        ),
                        height=300,
                        margin=dict(l=20, r=20, t=20, b=20),
                        showlegend=False
                    )
                    # Only map SCF points to geometry when the counts match.
                    aligned = len(energies) == len(coords)
                    fig.update_layout(clickmode='event+select')
                    selection = st.plotly_chart(
                        fig, use_container_width=True, key='scf_energy_chart',
                        on_select='rerun' if aligned else 'ignore',
                        selection_mode='points',
                    )
                    if aligned:
                        points = selection.selection.points
                        if points:
                            clicked = int(points[-1]['point_index']) + 1
                            if st.session_state.get('last_plot_selection') != clicked:
                                st.session_state['last_plot_selection'] = clicked
                                st.session_state['pending_geom_step'] = clicked
                                st.rerun()
                        st.caption('Click an SCF energy marker to display its geometry.')
                    else:
                        st.caption('Click-to-geometry navigation is disabled because SCF energies and saved geometries do not have a reliable one-to-one mapping.')
                else:
                    st.warning('No finite SCF energies are available to plot.')
                if len(energies)!=len(coords):
                    st.caption('SCF energy entries may not correspond one-to-one to saved geometry frames.')
            xyz=xyz_frame(numbers,coords[step-1],f'Geometry {step}')
            st.download_button('Download selected geometry (.xyz)',xyz,file_name=f'geometry_{step}.xyz')
        with left:
            if playback:
                frame_strings=[xyz_frame(numbers,c,f'Step {i+1}') for i,c in enumerate(coords)]
                viewer(frames=frame_strings,style=style,bg=bg,spin=spin,animate=True)
            else:
                viewer(xyz=xyz,style=style,bg=bg,spin=spin)
