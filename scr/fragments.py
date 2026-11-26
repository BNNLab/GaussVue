"""Small library of prebuilt substituents and rigid attachment geometry."""
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# The attachment atom is the atom at the specified zero-based index.
FRAGMENTS = {
    'Alkyl': {
        'Methyl (Me)': ('C', 0, 1.52),
        'Ethyl (Et)': ('CC', 0, 1.52),
        'n-Propyl': ('CCC', 0, 1.52),
        'Isopropyl (i-Pr)': ('C(C)C', 0, 1.52),
        'tert-Butyl (t-Bu)': ('C(C)(C)C', 0, 1.54),
        'Trifluoromethyl (CF3)': ('C(F)(F)F', 0, 1.53),
    },
    'Oxygen and sulfur': {
        'Hydroxyl (OH)': ('O', 0, 1.43),
        'Methoxy (OMe)': ('OC', 0, 1.43),
        'Ethoxy (OEt)': ('OCC', 0, 1.43),
        'Thiol (SH)': ('S', 0, 1.82),
        'Methylthio (SMe)': ('SC', 0, 1.82),
        'Acetoxy (OAc)': ('OC(C)=O', 0, 1.43),
    },
    'Nitrogen': {
        'Amino (NH2)': ('N', 0, 1.47),
        'Methylamino (NHMe)': ('NC', 0, 1.47),
        'Dimethylamino (NMe2)': ('N(C)C', 0, 1.47),
        'Nitro (NO2)': ('[N+](=O)[O-]', 0, 1.48),
        'Azido (N3)': ('N=[N+]=[N-]', 0, 1.47),
    },
    'Carbonyl and unsaturated': {
        'Formyl (CHO)': ('C=O', 0, 1.50),
        'Acetyl (COMe)': ('C(=O)C', 0, 1.51),
        'Carboxyl (COOH)': ('C(=O)O', 0, 1.51),
        'Carboxamide (CONH2)': ('C(=O)N', 0, 1.51),
        'Cyano (CN)': ('C#N', 0, 1.46),
        'Vinyl': ('C=C', 0, 1.50),
        'Ethynyl': ('C#C', 0, 1.46),
    },
    'Aromatic and cyclic': {
        'Phenyl (Ph)': ('c1ccccc1', 0, 1.49),
        'Benzyl (Bn)': ('Cc1ccccc1', 0, 1.52),
        'Pyridyl (2-pyridyl)': ('c1ccccn1', 0, 1.49),
        'Cyclopropyl': ('C1CC1', 0, 1.51),
        'Cyclohexyl': ('C1CCCCC1', 0, 1.52),
    },
    'Halogens': {
        'Fluoro (F)': ('F', 0, 1.36),
        'Chloro (Cl)': ('Cl', 0, 1.77),
        'Bromo (Br)': ('Br', 0, 1.94),
        'Iodo (I)': ('I', 0, 2.14),
    },
}


def _unit(v):
    n = np.linalg.norm(v)
    if n < 1e-9:
        raise ValueError('Degenerate attachment direction')
    return v / n


def _rotation(a, b):
    """Rotation mapping unit vector a onto b."""
    a, b = _unit(a), _unit(b)
    v = np.cross(a, b)
    c = float(np.dot(a, b))
    if c < -0.999999:
        basis = np.array([1., 0., 0.]) if abs(a[0]) < 0.9 else np.array([0., 1., 0.])
        axis = _unit(np.cross(a, basis))
        return 2 * np.outer(axis, axis) - np.eye(3)
    if c > 0.999999:
        return np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K / (1 + c)


def attach_fragment(numbers, xyz, bonds, bond_types, anchor, name, bond_length=None):
    """Attach a fragment at a selected host atom, replacing one explicit H if present.

    Returns (numbers, coords, bonds, bond_types). Existing non-H host positions
    are preserved; the new fragment is rigidly oriented, not force-field relaxed.
    """
    entry = next((items[name] for items in FRAGMENTS.values() if name in items), None)
    if entry is None:
        raise ValueError('Unknown fragment')
    smiles, attachment, default_length = entry
    length = float(default_length if bond_length is None else bond_length)
    if not 0.5 <= length <= 4.0:
        raise ValueError('Attachment distance must be between 0.5 and 4.0 Angstrom')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = 2026
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError('Could not generate fragment coordinates')
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    elif AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    frag = np.array(mol.GetConformer().GetPositions())
    numbers = list(map(int, numbers))
    xyz = np.array(xyz, float).copy()
    bonds = [tuple(map(int, e)) for e in bonds]
    bond_types = dict(bond_types)
    if not 0 <= anchor < len(numbers):
        raise ValueError('Invalid host attachment atom')
    # Remove one explicitly bonded H from the host, when available.
    attached_h = next((j if i == anchor else i for i, j in bonds
                       if (i == anchor and numbers[j] == 1) or (j == anchor and numbers[i] == 1)), None)
    if attached_h is not None:
        outward = xyz[attached_h] - xyz[anchor]
    else:
        neighbors = [j if i == anchor else i for i, j in bonds if i == anchor or j == anchor]
        outward = xyz[anchor] - np.mean(xyz[neighbors], axis=0) if neighbors else np.array([1., 0., 0.])
    outward = _unit(outward)
    # Remove one hydrogen from fragment attachment atom to make the new bond.
    frag_h = next((n.GetIdx() for n in mol.GetAtomWithIdx(attachment).GetNeighbors()
                   if n.GetAtomicNum() == 1), None)
    if frag_h is None:
        neighbor_ids = [n.GetIdx() for n in mol.GetAtomWithIdx(attachment).GetNeighbors()]
        local_out = (frag[attachment] - np.mean(frag[neighbor_ids], axis=0)) if neighbor_ids else np.array([1., 0., 0.])
    else:
        local_out = frag[frag_h] - frag[attachment]
    rot = _rotation(local_out, -outward)
    shifted = (frag - frag[attachment]) @ rot.T + xyz[anchor] + length * outward
    remove_host = {attached_h} if attached_h is not None else set()
    keep_host = [i for i in range(len(numbers)) if i not in remove_host]
    host_map = {old: new for new, old in enumerate(keep_host)}
    keep_frag = [i for i in range(mol.GetNumAtoms()) if i != frag_h]
    frag_map = {old: len(keep_host) + k for k, old in enumerate(keep_frag)}
    new_nums = [numbers[i] for i in keep_host] + [mol.GetAtomWithIdx(i).GetAtomicNum() for i in keep_frag]
    new_xyz = np.vstack([xyz[keep_host], shifted[keep_frag]])
    new_bonds = []
    new_types = {}
    for i, j in bonds:
        if i in host_map and j in host_map:
            pair = tuple(sorted((host_map[i], host_map[j])))
            new_bonds.append(pair)
            new_types[pair] = bond_types.get(tuple(sorted((i, j))), 'Single')
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in frag_map and j in frag_map:
            pair = tuple(sorted((frag_map[i], frag_map[j])))
            new_bonds.append(pair)
            new_types[pair] = ('Aromatic' if b.GetIsAromatic() else
                               {'DOUBLE':'Double', 'TRIPLE':'Triple'}.get(str(b.GetBondType()), 'Single'))
    pair = tuple(sorted((host_map[anchor], frag_map[attachment])))
    new_bonds.append(pair)
    new_types[pair] = 'Single'
    return new_nums, new_xyz, new_bonds, new_types
