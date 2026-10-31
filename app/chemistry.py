"""RDKit analysis and SVG rendering for the molecule explorer."""

import re

from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

from .condensed import MAX_ATOMS, MAX_INPUT, normalize_condensed, parse_condensed

# SMARTS distinguish adjacent functional classes so, for example, an ester is
# not also reported as an ether or a ketone.
_GROUPS = [
    ("carboxylic-acid", "Carboxylic acid", "#ec6b57", "[CX3](=[OX1])[OX2H1]"),
    ("ester", "Ester", "#e59b36", "[CX3](=[OX1])[OX2H0][#6;!$(C=O)]"),
    ("amide", "Amide", "#8b76d1", "[CX3](=[OX1])[NX3]"),
    ("aldehyde", "Aldehyde", "#e07b9a", "[CX3;H2,$([CX3H1]-[#6])](=[OX1])"),
    ("ketone", "Ketone", "#da8d45", "[#6][CX3](=[OX1])[#6]"),
    ("alcohol", "Alcohol", "#35a98a", "[CX4][OX2H1]"),
    ("phenol", "Phenol", "#4c9b67", "[c][OX2H1]"),
    ("amine", "Amine", "#6588cf", "[NX3;$([N]-[#6]);!$(N-[C,S,P]=[O,S,N]);!$(N=O);!$(N-O)]"),
    ("ether", "Ether", "#58a9bc", "[#6;!$(C=O)][OX2H0][#6;!$(C=O)]"),
    ("alkene", "Alkene", "#c0a131", "[CX3]=[CX3]"),
    ("alkyne", "Alkyne", "#a3a34d", "[CX2]#[CX2]"),
    ("nitrile", "Nitrile", "#9872bd", "[CX2]#[NX1]"),
    ("nitro", "Nitro", "#cc6f9c", "[NX3+](=[OX1])[O-]"),
    ("organic-halide", "Organic halide", "#6ea76c", "[#6][F,Cl,Br,I]"),
    ("aromatic-ring", "Aromatic ring", "#729dbb", None),
]
_QUERIES = [(key, name, color, Chem.MolFromSmarts(smarts) if smarts else None)
            for key, name, color, smarts in _GROUPS]
_HYBRIDIZATION = {Chem.HybridizationType.SP: "sp", Chem.HybridizationType.SP2: "sp2",
                  Chem.HybridizationType.SP3: "sp3", Chem.HybridizationType.S: "s"}


def _parse_input(text: str, input_format: str) -> tuple[Chem.Mol, str]:
    if not isinstance(text, str) or not text.strip():
        raise ValueError("Enter a molecule, for example CCO or CH3CH2OH.")
    if len(text) > MAX_INPUT:
        raise ValueError(f"Please keep input to {MAX_INPUT} characters or fewer.")
    if input_format not in ("auto", "smiles", "condensed"):
        raise ValueError("Choose auto, smiles, or condensed as the input format.")
    text = text.strip()
    if input_format == "auto":
        normalized = normalize_condensed(text)
        # Bracket notation always belongs to SMILES. Explicit hydrogens outside
        # brackets mark conventional condensed formulas, including molecular-only
        # formulas that must receive the connectivity error rather than a guess.
        input_format = "condensed" if "[" not in normalized and re.search(r"H(?:\d|(?=[A-Z()]|$))", normalized) else "smiles"
    if input_format == "condensed":
        return parse_condensed(text), input_format
    if any(char.isspace() for char in text):
        raise ValueError("Use only the SMILES structure, without names or spaces.")
    # Input errors are returned to the UI instead of emitted as server log noise.
    params = Chem.SmilesParserParams()
    params.parseName = False
    params.allowCXSMILES = False
    with rdBase.BlockLogs():
        mol = Chem.MolFromSmiles(text, params)
    if mol is None:
        raise ValueError("This SMILES could not be read. Check atom symbols, ring numbers, "
                         "brackets, and bonds; try CCO for ethanol.")
    if mol.GetNumHeavyAtoms() == 0 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("Enter a concrete molecule containing at least one non-hydrogen atom.")
    if mol.GetNumAtoms() > MAX_ATOMS:
        raise ValueError(f"Please use a structure with at most {MAX_ATOMS} atoms.")
    # Explicit standalone H atoms do not have skeletal drawing or hybridization
    # value. Preserve isotope H where RDKit requires it.
    return Chem.RemoveHs(mol), input_format


def _color(hex_color: str) -> tuple[float, float, float]:
    return tuple(int(hex_color[i:i + 2], 16) / 255 for i in (1, 3, 5))


def analyze_molecule(text: str, input_format: str = "auto", show_hybridization: bool = False,
                     highlight_groups: list[str] | None = None) -> dict:
    mol, resolved_format = _parse_input(text, input_format)
    atoms = [{"index": a.GetIdx(), "symbol": a.GetSymbol(),
              "hybridization": _HYBRIDIZATION.get(a.GetHybridization(), "other"),
              "formal_charge": a.GetFormalCharge(), "hydrogen_count": a.GetTotalNumHs(includeNeighbors=True),
              "aromatic": a.GetIsAromatic()} for a in mol.GetAtoms()]
    groups = []
    atom_colors, bond_colors = {}, {}
    selected = set(highlight_groups or [])
    for key, name, color, query in _QUERIES:
        matches = mol.GetSubstructMatches(query) if query is not None else tuple(
            tuple(ring) for ring in mol.GetRingInfo().AtomRings()
            if all(mol.GetAtomWithIdx(index).GetIsAromatic() for index in ring)
        )
        if not matches:
            continue
        indices = sorted({index for match in matches for index in match})
        groups.append({"id": key, "name": name, "color": color,
                       "matches": [list(match) for match in matches], "atom_indices": indices})
        if key in selected:
            for index in indices:
                atom_colors[index] = _color(color)
            for match in matches:
                # Query bond mapping avoids highlighting incidental bonds merely
                # because both endpoints belong to some match of the same class.
                endpoints = [(match[b.GetBeginAtomIdx()], match[b.GetEndAtomIdx()])
                             for b in query.GetBonds()] if query is not None else list(zip(match, match[1:] + match[:1]))
                for start, end in endpoints:
                    actual = mol.GetBondBetweenAtoms(start, end)
                    bond_colors[actual.GetIdx()] = _color(color)
    rdDepictor.Compute2DCoords(mol)
    if show_hybridization:
        for atom, info in zip(mol.GetAtoms(), atoms):
            atom.SetProp("atomNote", f"{info['index'] + 1}: {info['hybridization']}")
    drawer = rdMolDraw2D.MolDraw2DSVG(900, 560)
    options = drawer.drawOptions()
    options.clearBackground = False
    options.padding = 0.12
    options.annotationFontScale = 0.7
    options.addAtomIndices = False
    options.bondLineWidth = 2.2
    options.useBWAtomPalette()
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol, highlightAtoms=list(atom_colors),
                                      highlightBonds=list(bond_colors),
                                      highlightAtomColors=atom_colors,
                                      highlightBondColors=bond_colors)
    drawer.TagAtoms(mol)
    drawer.FinishDrawing()
    warnings = []
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        warnings.append("This structure contains radical atoms; review their intended valence.")
    if Chem.GetMolFrags(mol) and len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("This input contains multiple disconnected fragments.")
    if any(info["hybridization"] == "other" for info in atoms):
        warnings.append("Some atoms use a hybridization outside sp, sp2, or sp3; inspect those atoms individually.")
    return {"canonical_smiles": Chem.MolToSmiles(mol),
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": round(Descriptors.MolWt(mol), 3),
            "atom_count": mol.GetNumHeavyAtoms(), "bond_count": mol.GetNumBonds(),
            "input_format": resolved_format, "svg": drawer.GetDrawingText(),
            "atoms": atoms, "functional_groups": groups, "warnings": warnings}
