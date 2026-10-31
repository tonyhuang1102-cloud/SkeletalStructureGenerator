import pytest
from rdkit import Chem

from app.chemistry import analyze_molecule


@pytest.mark.parametrize(("formula", "smiles"), [
    ("CH3CH2OH", "CCO"), ("CH3COOH", "CC(=O)O"),
    ("CH3COCH3", "CC(=O)C"), ("CH3CHO", "CC=O"),
    ("CH3CH(CH3)OH", "CC(C)O"), ("CH3(CH2)4CH3", "CCCCCC"),
    ("(CH3)3COH", "CC(C)(C)O"), ("CH3COOCH2CH3", "CC(=O)OCC"),
    ("CH3CH=CH2", "CC=C"), ("CH3C≡N", "CC#N"),
    ("CH3CH2NH2", "CCN"), ("CH3OCH3", "COC"),
    ("CH3CH(CH3)2", "CC(C)C"), ("CH3(CH2CH2)3CH3", "CCCCCCCC"),
    ("CH3C(=O)OH", "CC(=O)O"), ("CH3Cl", "CCl"),
    ("CH₃CH₂OH", "CCO"), ("CH4", "C"), ("OH2", "O"),
])
def test_condensed_connectivity(formula, smiles):
    result = analyze_molecule(formula, "condensed")
    assert result["canonical_smiles"] == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    assert result["input_format"] == "condensed"


@pytest.mark.parametrize("formula", [
    "C2H6O", "C6H6", "CH3CH3O", "CH3CH2", "CH5", "CH3(CH3)CH3",
    "CH3(CH2)0CH3", "CH3(CH2)9999CH3", "CH3(CH2)149CH3",
    "CH3CH(CH3", "CH3CH2OH)", "CH3()OH", "CH3C=", "CH3C==O",
    "CH3Ph", "CH3..OH", "(CH3)3", "=CH2", "CH3C(=O)(=O)OH",
])
def test_invalid_condensed_formula_is_rejected(formula):
    with pytest.raises(ValueError):
        analyze_molecule(formula, "condensed")


@pytest.mark.parametrize(("text", "expected"), [
    ("CCO", "smiles"), ("CH3CH2OH", "condensed"),
    ("CH3COOH", "condensed"), ("[NH4+]", "smiles"),
    ("c1cc[nH]c1", "smiles"),
])
def test_auto_detection(text, expected):
    assert analyze_molecule(text)["input_format"] == expected


@pytest.mark.parametrize(("smiles", "present", "absent"), [
    ("CCO", {"alcohol"}, {"ether", "phenol"}),
    ("CC(=O)O", {"carboxylic-acid"}, {"alcohol", "ketone", "ether"}),
    ("CC(=O)OC", {"ester"}, {"carboxylic-acid", "ketone", "ether"}),
    ("CC(=O)OC(=O)C", set(), {"ester", "ether", "carboxylic-acid"}),
    ("CC(=O)N", {"amide"}, {"amine", "ketone"}),
    ("CC(C)=O", {"ketone"}, {"aldehyde", "ester"}),
    ("CC=O", {"aldehyde"}, {"ketone"}),
    ("C=O", {"aldehyde"}, {"ketone"}),
    ("O=CO", {"carboxylic-acid"}, {"aldehyde"}),
    ("O=COC", {"ester"}, {"aldehyde"}),
    ("NC=O", {"amide"}, {"aldehyde", "amine"}),
    ("c1ccccc1O", {"phenol", "aromatic-ring"}, {"alcohol", "alkene"}),
    ("c1ccoc1", {"aromatic-ring"}, {"alkene", "ether"}),
    ("CCN", {"amine"}, {"amide"}),
    ("N", set(), {"amine"}),
    ("NN", set(), {"amine"}),
    ("COC", {"ether"}, {"alcohol"}),
    ("C=C", {"alkene"}, {"alkyne"}),
    ("C#C", {"alkyne"}, {"alkene", "nitrile"}),
    ("CC#N", {"nitrile"}, {"alkyne", "amine"}),
    ("C[N+](=O)[O-]", {"nitro"}, {"amine"}),
    ("CCBr", {"organic-halide"}, set()),
])
def test_functional_group_specificity(smiles, present, absent):
    actual = {group["id"] for group in analyze_molecule(smiles)["functional_groups"]}
    assert present <= actual
    assert not absent & actual


def test_atom_data_and_molecular_properties():
    result = analyze_molecule("CC(=O)O")
    assert result["molecular_formula"] == "C2H4O2"
    assert result["molecular_weight"] == pytest.approx(60.052, abs=0.01)
    assert result["atom_count"] == 4
    assert result["bond_count"] == 3
    assert [atom["hybridization"] for atom in result["atoms"]] == ["sp3", "sp2", "sp2", "sp2"]
    assert [atom["index"] for atom in result["atoms"]] == [0, 1, 2, 3]
    assert result["atoms"][0]["hydrogen_count"] == 3
    assert analyze_molecule("C#N")["atoms"][0]["hybridization"] == "sp"


def test_svg_highlights_and_picking():
    result = analyze_molecule("CCO", highlight_groups=["alcohol"], show_hybridization=True)
    assert "<svg" in result["svg"]
    assert "atom-selector atom-0" in result["svg"]
    assert "note" in result["svg"]
    assert "#35A98A" in result["svg"]
    assert result["functional_groups"][0]["atom_indices"] == [1, 2]
    assert result["functional_groups"][0]["matches"] == [[1, 2]]


def test_multiple_matches_have_stable_indices():
    result = analyze_molecule("OCCO")
    group = result["functional_groups"][0]
    assert group["id"] == "alcohol"
    assert len(group["matches"]) == 2
    assert group["atom_indices"] == [0, 1, 2, 3]


@pytest.mark.parametrize("text", ["", " ", "this-is-not-smiles", "C1CC", "CCO ethanol", "CCO<script>", "*", "[H][H]", "C" * 151, "C" * 501])
def test_invalid_or_oversized_smiles(text):
    with pytest.raises(ValueError):
        analyze_molecule(text, "smiles")


def test_ambiguous_formula_auto_does_not_guess():
    with pytest.raises(ValueError, match="connectivity"):
        analyze_molecule("C2H6O")


def test_fragment_and_radical_warnings():
    assert "disconnected" in analyze_molecule("[Na+].[Cl-]")["warnings"][0]
    assert "radical" in analyze_molecule("[CH3]")["warnings"][0]


def test_long_valid_chain_is_bounded_but_allowed():
    assert analyze_molecule("CH3(CH2)148CH3")["atom_count"] == 150


def test_attached_hydrogens_include_preserved_isotopes():
    result = analyze_molecule("[2H]C", "smiles")
    carbon = next(atom for atom in result["atoms"] if atom["symbol"] == "C")
    assert carbon["hydrogen_count"] == 4
