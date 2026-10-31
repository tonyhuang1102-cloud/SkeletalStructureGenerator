"""A deliberately bounded parser for explicit, acyclic condensed formulas.

Hydrogen counts are constraints, not guesses. Parenthesized chain units with two
open valences repeat in sequence; terminal substituents repeat as branches.
"""

from dataclasses import dataclass, field
import re

from rdkit import Chem

MAX_ATOMS = 150
MAX_INPUT = 500
_SUBSCRIPTS = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_ATOM = re.compile(r"(Cl|Br|[BCNOPSFI])(?:H(\d*))?")
_VALENCE = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "B": 3,
            "F": 1, "Cl": 1, "Br": 1, "I": 1}


@dataclass
class Fragment:
    atoms: list[tuple[str, int]] = field(default_factory=list)
    bonds: list[tuple[int, int, int]] = field(default_factory=list)
    first: int = 0
    tail: int = 0
    attachment_order: int = 1

    def append(self, other: "Fragment", anchor: int | None = None,
               order: int | None = None) -> tuple[int, int]:
        offset = len(self.atoms)
        if offset + len(other.atoms) > MAX_ATOMS:
            raise ValueError(f"Please use a structure with at most {MAX_ATOMS} atoms.")
        self.atoms.extend(other.atoms)
        self.bonds.extend((a + offset, b + offset, n) for a, b, n in other.bonds)
        first, tail = other.first + offset, other.tail + offset
        if anchor is not None:
            self.bonds.append((anchor, first, order or other.attachment_order))
        return first, tail

    def open_valence(self, index: int) -> int:
        symbol, hydrogens = self.atoms[index]
        used = sum(n for a, b, n in self.bonds if index in (a, b))
        return _VALENCE[symbol] - hydrogens - used


def normalize_condensed(text: str) -> str:
    return re.sub(r"\s+", "", text.translate(_SUBSCRIPTS)).replace("≡", "#").replace("−", "-")


class _Parser:
    def __init__(self, text: str):
        self.text = text
        self.pos = 0

    def atom(self) -> tuple[str, int]:
        match = _ATOM.match(self.text, self.pos)
        if not match:
            raise ValueError("Use explicit groups such as CH3CH2OH or CH3CH(CH3)OH. "
                             "Ring abbreviations and molecular formulas need SMILES instead.")
        self.pos = match.end()
        count = match.group(2)
        hydrogens = 0 if count is None else int(count or "1")
        if hydrogens > 4:
            raise ValueError("A group has too many hydrogens. Check its condensed formula.")
        return match.group(1), hydrogens

    def multiplicity(self) -> int:
        match = re.match(r"\d+", self.text[self.pos:])
        if not match:
            return 1
        self.pos += len(match.group())
        count = int(match.group())
        if not 1 <= count <= MAX_ATOMS:
            raise ValueError(f"Repeat counts must be between 1 and {MAX_ATOMS}.")
        return count

    def sequence(self, depth: int = 0) -> Fragment:
        if depth > 20:
            raise ValueError("This formula has too many nested branches.")
        result = Fragment()
        cursor = None
        pending_order = 1
        leading: list[tuple[Fragment, int]] = []
        expect_atom = True
        while self.pos < len(self.text) and self.text[self.pos] != ")":
            char = self.text[self.pos]
            if char in "-=#":
                if pending_order != 1 or (expect_atom and result.atoms):
                    raise ValueError("Check the bond symbols in your condensed formula.")
                pending_order = {"-": 1, "=": 2, "#": 3}[char]
                if not result.atoms:
                    result.attachment_order = pending_order
                    pending_order = 1
                self.pos += 1
                expect_atom = True
                if self.pos == len(self.text) or self.text[self.pos] in ")-=#":
                    raise ValueError("A bond must be followed by an atom or group.")
                continue
            if char == "(":
                self.pos += 1
                part = self.sequence(depth + 1)
                if self.pos >= len(self.text) or self.text[self.pos] != ")":
                    raise ValueError("Close every opening parenthesis in the formula.")
                self.pos += 1
                repeated = self.pos < len(self.text) and self.text[self.pos].isdigit()
                count = self.multiplicity()
                if not part.atoms:
                    raise ValueError("Parentheses must contain an atom or group.")
                if cursor is None:
                    leading.append((part, count))
                    continue
                # A repeated unit that retains a second open valence extends the
                # chain. Saturated-end fragments instead repeat as substituents.
                chain = repeated and part.attachment_order == 1 and (
                    part.open_valence(part.tail) >= (2 if part.first == part.tail else 1)
                )
                for _ in range(count):
                    _, tail = result.append(part, cursor, pending_order if expect_atom else None)
                    if chain:
                        cursor = tail
                pending_order = 1
                expect_atom = False
                continue
            if char.isdigit():
                raise ValueError("A molecular formula does not specify connectivity. "
                                 "Use CH3CH2OH rather than C2H6O, or switch to SMILES.")
            atom = self.atom()
            part = Fragment(atoms=[atom])
            # Conventional CO / CHO shorthand means a carbonyl. Its oxygen is
            # a side atom; the following group remains attached to the carbon.
            if atom[0] == "C" and atom[1] <= 1:
                oxygen = _ATOM.match(self.text, self.pos)
                if oxygen and oxygen.group(1) == "O" and oxygen.group(2) is None:
                    self.pos = oxygen.end()
                    part.atoms.append(("O", 0))
                    part.bonds.append((0, 1, 2))
            first, tail = result.append(part, cursor, pending_order)
            if cursor is None:
                result.first = first
                for branch, count in leading:
                    for _ in range(count):
                        result.append(branch, first)
                leading.clear()
            cursor = tail
            pending_order = 1
            expect_atom = False
        if leading:
            # A single parenthesized formula is allowed, but an unattached bundle
            # of repeated substituents has no specified central atom.
            if not result.atoms and len(leading) == 1 and leading[0][1] == 1:
                return leading[0][0]
            raise ValueError("Repeated leading groups need a central atom, for example (CH3)3COH.")
        if cursor is not None:
            result.tail = cursor
        return result


def parse_condensed(text: str) -> Chem.Mol:
    if not text or len(text) > MAX_INPUT:
        raise ValueError(f"Enter a condensed formula of at most {MAX_INPUT} characters.")
    normalized = normalize_condensed(text)
    parser = _Parser(normalized)
    fragment = parser.sequence()
    if fragment.attachment_order != 1:
        raise ValueError("A formula cannot begin with a double or triple bond.")
    if parser.pos != len(normalized):
        raise ValueError("There is an extra closing parenthesis in the formula.")
    if not fragment.atoms:
        raise ValueError("Enter a condensed formula such as CH3CH2OH.")
    # This catches both overfilled atoms and incomplete groups before RDKit can
    # infer missing hydrogen or produce a radical from an underfilled atom.
    for index, (symbol, _) in enumerate(fragment.atoms):
        if fragment.open_valence(index) != 0:
            raise ValueError(f"The hydrogen count or bonds at {symbol} (atom {index + 1}) "
                             "do not complete its usual valence. Use explicit hydrogen "
                             "counts, = / # bonds, or SMILES for charged or unusual structures.")
    mol = Chem.RWMol()
    for symbol, hydrogens in fragment.atoms:
        atom = Chem.Atom(symbol)
        atom.SetNumExplicitHs(hydrogens)
        atom.SetNoImplicit(True)
        mol.AddAtom(atom)
    bond_types = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE, 3: Chem.BondType.TRIPLE}
    for a, b, order in fragment.bonds:
        mol.AddBond(a, b, bond_types[order])
    result = mol.GetMol()
    try:
        Chem.SanitizeMol(result)
    except Exception as exc:
        raise ValueError("This condensed formula gives an invalid valence. Check hydrogen counts and bonds.") from exc
    return result
