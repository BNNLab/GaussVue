import os
import tempfile
import streamlit as st
from cclib.io import ccread
from periodictable import elements

@st.cache_data(show_spinner='Parsing Gaussian output…')
def parse_log(raw):
    # cclib supports file-like objects for Gaussian output.
    with tempfile.NamedTemporaryFile(suffix='.log', delete=False) as tmp:
        tmp.write(raw)
        filename = tmp.name
    try:
        data = ccread(filename)
        if data is None:
            raise ValueError('cclib could not identify this file as a supported calculation.')
        return data
    finally:
        os.unlink(filename)


def xyz_frame(numbers, coords, label='Molecule'):
    lines = [str(len(numbers)), label]
    for z, xyz in zip(numbers, coords):
        symbol = elements[int(z)].symbol
        lines.append(f'{symbol:<3} {xyz[0]: .8f} {xyz[1]: .8f} {xyz[2]: .8f}')
    return '\n'.join(lines) + '\n'


def pdb_block(nums, xyz, bonds, bond_types=None):
    lines = ['HEADER    GAUSSIAN EXPLORER EDITED STRUCTURE']
    for i,(z,pos) in enumerate(zip(nums,xyz), 1):
        sym = elements[int(z)].symbol
        lines.append(f'HETATM{i:5d} {sym:>4s} MOL A   1    '
                     f'{pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}'
                     f'{1.:6.2f}{0.:6.2f}          {sym:>2s}  ')
    for i,j in bonds:
        kind=(bond_types or {}).get(tuple(sorted((i,j))),'Single')
        if kind=='H-bond': continue
        repeat={'Double':2,'Triple':3}.get(kind,1)
        lines.append(f'CONECT{i+1:5d}'+f'{j+1:5d}'*repeat)
    return '\n'.join(lines + ['END',''])


def gaussian_input(nums, xyz, command, functional, basis, solvation, solvent,
                   charge, multiplicity, title, nproc, memory, chk):
    # The route accepts user-supplied Gaussian keywords without executing them.
    route = ' '.join(v.strip() for v in [command, f'{functional.strip()}/{basis.strip()}'
                                       if functional.strip() and basis.strip() else '',
                                       f'SCRF=(SMD,Solvent={solvent.strip()})' if solvation == 'SMD' else
                                       f'SCRF=(PCM,Solvent={solvent.strip()})' if solvation == 'PCM' else
                                       f'SCRF=(CPCM,Solvent={solvent.strip()})' if solvation == 'CPCM' else ''] if v.strip())
    if not route: raise ValueError('Specify at least one Gaussian route keyword.')
    if any('\n' in x or '\r' in x for x in [memory,chk,title,route]):
        raise ValueError('Gaussian header fields must be single-line.')
    if int(multiplicity) < 1: raise ValueError('Multiplicity must be at least 1.')
    lines = [f'%chk={chk}', f'%mem={memory}', f'%nprocshared={int(nproc)}',
             f'#p {route}', '', title, '', f'{int(charge)} {int(multiplicity)}']
    for z,pos in zip(nums,xyz):
        lines.append(f'{elements[int(z)].symbol:<3} {pos[0]: .8f} {pos[1]: .8f} {pos[2]: .8f}')
    return '\n'.join(lines + ['',''])

