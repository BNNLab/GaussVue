import numpy as np
import streamlit as st
from .gaussian_io import xyz_frame
from .viewer import viewer


def render(numbers, coords, frequencies, displacements, data, style, bg, spin):
    if not len(frequencies) or not len(displacements) or not len(coords):
        st.info('Vibrational modes were not found. Upload a Gaussian frequency calculation with printed normal-mode displacement vectors.')
    else:
        left,right=st.columns([3,2])
        with right:
            idx=st.selectbox('Normal mode',list(range(len(frequencies))),
                             format_func=lambda i:f'Mode {i+1} | {frequencies[i]:.2f} cm⁻¹'+(' (imaginary)' if frequencies[i]<0 else ''))
            amp=st.slider('Displacement amplitude (visual scale)',0.1,2.5,0.7,0.1)
            nframes=st.slider('Animation frames',8,40,20,4)
            st.metric('Frequency',f'{frequencies[idx]:.2f} cm⁻¹')
            if hasattr(data,'vibirs') and idx<len(data.vibirs):
                st.metric('IR intensity',f'{data.vibirs[idx]:.2f} km/mol')
            st.caption('Animation is a visualisation of the normal-mode eigenvector, not a time-resolved molecular dynamics simulation.')
        with left:
            base=coords[-1]
            disp=displacements[idx]
            if disp.shape!=base.shape:
                st.error('Normal-mode displacement shape does not match the final geometry.')
            else:
                frames=[xyz_frame(numbers,base+amp*np.sin(2*np.pi*k/nframes)*disp,f'Mode {idx+1}') for k in range(nframes)]
                viewer(frames=frames,style=style,bg=bg,spin=spin,animate=True)
        import pandas as pd
        with st.expander('Frequency table'):
            tbl={'Mode':np.arange(1,len(frequencies)+1),'Frequency (cm-1)':frequencies}
            if hasattr(data,'vibirs') and len(data.vibirs)==len(frequencies):tbl['IR intensity (km/mol)']=data.vibirs
            st.dataframe(pd.DataFrame(tbl),use_container_width=True,hide_index=True)
