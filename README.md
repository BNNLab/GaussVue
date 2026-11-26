# GaussVue ![GaussVue](images/GaussVue_logo_small.png)

Run from the project root:

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

Put your existing `BNNLab_v3.png` in `images/`. The app will run without it.

- `app.py`: configuration, sidebar, file uploads, tab orchestration
- `style.css`: styling
- `scr/gaussian_io.py`: Gaussian parsing and XYZ/PDB/GJF exports
- `scr/geometry.py`: geometry operations and SMILES generation
- `scr/viewer.py`: 3Dmol viewer for optimisation and vibrations
- `scr/editor_component.py`: bidirectional 3Dmol atom-picker
- `scr/optimisation.py`, `vibrations.py`, `orbitals.py`, `calculation.py`, `structure_editor.py`: tab UIs

The existing calculation/editor behaviour is retained; no chemistry algorithms were intentionally changed. Bond order is not encoded in Gaussian XYZ-style input.

## Premade R-group fragments

In **Structure editor**, choose **Add fragment**, select a category and substituent,
click an attachment atom in the 3D viewer, adjust the attachment distance and
press **Attach fragment**. The selected host atom loses one explicit hydrogen if
available. The new fragment is positioned rigidly and can be undone.

The library is in `scr/fragments.py`. Add new entries to `FRAGMENTS` as
`'Display name': ('SMILES', attachment_atom_index, default_bond_length_A)`.
The geometry is an initial guess only. Check valence, steric clashes and
formal charge, then optimise before production Gaussian calculations.
