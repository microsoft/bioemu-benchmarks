"""Taken from https://www.cup.uni-muenchen.de/ch/compchem/tink/as.html"""

RESTYPE_1TO3: dict[str, str] = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "U": "SEC",  # Selenocysteine
    "O": "PYL",  # Pyrrolysine
    "X": "UNK",  # Unknown
    "B": "NLE",  # Norleucine, Found in DEShaw WW_DOMAIN, isomer of Leucine found in some bacterial strains.
    "Z": "GLX",  # Glutamine or glutamic acid. Used for when we can't distinguish between GLU and GLN.
}

RESTYPE_3TO1: dict[str, str] = {v: k for k, v in RESTYPE_1TO3.items()}


def get_aa1code_from_aa3code(aa3code: str) -> str:
    if aa3code in RESTYPE_3TO1:
        return RESTYPE_3TO1[aa3code]
    else:
        return "X"
