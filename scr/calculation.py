import numpy as np
import streamlit as st


def render(data):
    if data is None:
        st.info('No data to visualise yet!')
    else:
        info={
            'Method':str(getattr(data,'methods','Not parsed')),
            'Basis':str(getattr(data,'metadata',{}).get('basis_set','Not parsed')),
            'Charge':str(getattr(data,'charge','Not parsed')),
            'Multiplicity':str(getattr(data,'mult','Not parsed')),
            'Optimization status':str(getattr(data,'optdone','Not parsed')),
        }
        for k,v in info.items():st.write(f'**{k}:** {v}')
        if hasattr(data,'moenergies'):
            import pandas as pd
            for s,levels in enumerate(data.moenergies):
                st.subheader(f'MO energies: spin channel {s+1}')
                homo=int(data.homos[s]) if hasattr(data,'homos') and s<len(data.homos) else -1
                st.dataframe(pd.DataFrame({'MO':np.arange(1,len(levels)+1),'Energy (eV)':levels,
                                           'Occupancy label':['HOMO' if i==homo else 'LUMO' if i==homo+1 else '' for i in range(len(levels))]}),
                             use_container_width=True,hide_index=True)
        if hasattr(data,'freeenergy'):st.write(f'**Gibbs free energy:** {data.freeenergy:.9f} Eh')
        if hasattr(data,'enthalpy'):st.write(f'**Enthalpy:** {data.enthalpy:.9f} Eh')
