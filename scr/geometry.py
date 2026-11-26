import numpy as np
from periodictable import elements

# --------------------------- Geometry editor helpers ---------------------------
# Atom indices shown in the UI are 1-based; internal arrays are 0-based.
# Connectivity is inferred from covalent radii and should be reviewed, especially
# for metals, weak complexes and unusually long/short bonds.
COVALENT_RADII = {
    1: .31, 5: .84, 6: .76, 7: .71, 8: .66, 9: .57, 14: 1.11,
    15: 1.07, 16: 1.05, 17: 1.02, 35: 1.20, 53: 1.39,
}


def inferred_bonds(nums, xyz):
    result = []
    for i in range(len(nums)):
        for j in range(i + 1, len(nums)):
            ri = COVALENT_RADII.get(int(nums[i]), 1.15)
            rj = COVALENT_RADII.get(int(nums[j]), 1.15)
            distance = float(np.linalg.norm(xyz[i] - xyz[j]))
            if .35 < distance < 1.22 * (ri + rj):
                result.append((i, j))
    return result


def graph_component(n, bonds, start, cut=None):
    graph = [[] for _ in range(n)]
    for i, j in bonds:
        if cut is not None and {i, j} == set(cut):
            continue
        graph[i].append(j)
        graph[j].append(i)
    seen = {start}
    todo = [start]
    while todo:
        for nxt in graph[todo.pop()]:
            if nxt not in seen:
                seen.add(nxt)
                todo.append(nxt)
    return seen


def unit(vector):
    norm = float(np.linalg.norm(vector))
    if norm < 1e-10:
        raise ValueError('The selected atoms define a degenerate direction.')
    return vector / norm


def rotate_group(xyz, indices, origin, axis, degrees):
    result = xyz.copy()
    k = unit(axis)
    theta = np.deg2rad(degrees)
    for i in indices:
        v = result[i] - origin
        result[i] = origin + (v * np.cos(theta) + np.cross(k, v) * np.sin(theta)
                              + k * np.dot(k, v) * (1 - np.cos(theta)))
    return result


def angle_deg(a, b, c):
    v1, v2 = unit(a-b), unit(c-b)
    return float(np.rad2deg(np.arccos(np.clip(np.dot(v1, v2), -1, 1))))


def dihedral_deg(a, b, c, d):
    axis = unit(c-b)
    v = (a-b) - np.dot(a-b, axis)*axis
    w = (d-c) - np.dot(d-c, axis)*axis
    v, w = unit(v), unit(w)
    return float(np.rad2deg(np.arctan2(np.dot(np.cross(v, w), axis), np.dot(v, w))))


def edit_geometry(nums, xyz, bonds, operation, ids, value=None, symbol=None):
    nums, xyz, bonds = list(nums), np.asarray(xyz, float).copy(), list(bonds)
    n = len(nums)
    if any(i < 0 or i >= n for i in ids):
        raise ValueError('Atom index is outside the current structure.')
    if len(set(ids)) != len(ids):
        raise ValueError('Select distinct atoms.')
    if operation == 'Add atom':
        anchor = ids[0]
        z = elements.symbol(symbol).number
        outward = xyz[anchor] - xyz.mean(axis=0)
        if np.linalg.norm(outward) < 1e-8:
            outward = np.array([1., 0., 0.])
        new_xyz = xyz[anchor] + unit(outward) * float(value)
        xyz = np.vstack([xyz, new_xyz])
        nums.append(z)
        bonds.append((anchor, n))
    elif operation == 'Remove atom':
        atom = ids[0]
        nums.pop(atom)
        xyz = np.delete(xyz, atom, axis=0)
        bonds = [(i - (i > atom), j - (j > atom)) for i, j in bonds if i != atom and j != atom]
    elif operation == 'Bond length':
        a, b = ids
        if {a,b} not in [set(edge) for edge in bonds]:
            raise ValueError('These atoms are not connected by an inferred bond.')
        moving = graph_component(n, bonds, b, (a,b))
        if a in moving:
            raise ValueError('This bond is part of a ring. A fragment-preserving bond edit is ambiguous.')
        shift = (float(value) - np.linalg.norm(xyz[b]-xyz[a])) * unit(xyz[b]-xyz[a])
        xyz[list(moving)] += shift
    elif operation == 'Bond angle':
        a, b, c = ids
        if {b,c} not in [set(edge) for edge in bonds]:
            raise ValueError('The central and last atoms must be bonded.')
        moving = graph_component(n, bonds, c, (b,c))
        if b in moving or a in moving:
            raise ValueError('Cannot rotate this fragment without moving other selected atoms (ring/bridge).')
        v1, v2 = unit(xyz[a]-xyz[b]), unit(xyz[c]-xyz[b])
        axis = np.cross(v2, v1)
        current = angle_deg(xyz[a], xyz[b], xyz[c])
        if np.linalg.norm(axis) < 1e-9:
            axis = np.cross(v2, [1.,0.,0.])
            if np.linalg.norm(axis) < 1e-9: axis = np.cross(v2, [0.,1.,0.])
        xyz = rotate_group(xyz, moving, xyz[b], axis, current-float(value))
    elif operation == 'Dihedral angle':
        a,b,c,d = ids
        if {b,c} not in [set(edge) for edge in bonds]:
            raise ValueError('The central atoms must be bonded.')
        moving = graph_component(n, bonds, d, (b,c))
        if b in moving or a in moving:
            raise ValueError('Cannot rotate the selected side of a ring or bridged bond.')
        current = dihedral_deg(xyz[a],xyz[b],xyz[c],xyz[d])
        xyz = rotate_group(xyz, moving, xyz[c], xyz[c]-xyz[b], float(value)-current)
    return nums, xyz, bonds



def smiles_to_3d(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid SMILES string.')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 2026
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError('RDKit could not embed this molecule in 3D.')
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=300)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=300)
    nums = [a.GetAtomicNum() for a in mol.GetAtoms()]
    xyz = np.array(mol.GetConformer().GetPositions(), dtype=float)
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    types = {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))): ('Aromatic' if b.GetIsAromatic() else {'SINGLE':'Single','DOUBLE':'Double','TRIPLE':'Triple'}.get(str(b.GetBondType()), 'Single')) for b in mol.GetBonds()}
    return nums, xyz, bonds, types


