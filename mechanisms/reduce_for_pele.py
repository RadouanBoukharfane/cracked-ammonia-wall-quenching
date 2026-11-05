"""Removes inert/carbon species not needed for NH3/H2/N2/air (AR, HE, OH*, C species) and the
reactions involving them (including explicit colliders such as 'H2 + AR'), drops their third-body
efficiencies and the reaction notes, and writes a PelePhysics-ready YAML. Works on the YAML tree
(Cantera's writer chokes on some notes converted from Chemkin files).

  python reduce_for_pele.py IN.yaml OUT.yaml
"""
import re
import sys

import cantera as ct
from ruamel.yaml import YAML

DROP = {"AR", "HE", "OH*", "CO", "CO2", "HCO", "CH4", "C2H6"}


def species_in(eq):
    toks = re.split(r"\s*(?:<=>|=>|\+|\(\+\s*M\s*\)|\(|\))\s*", eq)
    return {re.sub(r"^\d+\s*", "", t).strip() for t in toks if t.strip()}


def main(src, dst):
    y = YAML()
    y.width = 4096
    doc = y.load(open(src))
    names = [s["name"] for s in doc["species"]]
    drop = DROP & set(names)
    phase = doc["phases"][0]
    phase["species"] = [s for s in phase["species"] if s not in drop]
    phase["elements"] = [e for e in phase["elements"] if e not in ("C", "Ar", "He")]
    doc["species"] = [s for s in doc["species"] if s["name"] not in drop]
    for s in doc["species"]:
        s.pop("note", None)
    rx = []
    for r in doc["reactions"]:
        if species_in(r["equation"]) & drop:
            continue
        r.pop("note", None)
        if "efficiencies" in r:
            for k in list(r["efficiencies"]):
                if k in drop:
                    del r["efficiencies"][k]
        rx.append(r)
    doc["reactions"] = rx
    with open(dst, "w") as f:
        y.dump(doc, f)
    full, g = ct.Solution(src), ct.Solution(dst)
    print(f"{src}: {full.n_species} sp / {full.n_reactions} r -> {dst}: {g.n_species} sp / {g.n_reactions} r "
          f"(removed {sorted(drop)})")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
