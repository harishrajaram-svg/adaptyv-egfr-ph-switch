"""Fetch structure files into analysis/02-tnf/structures/ if absent, and return the path.

Why this exists: contacts.py and per_partner.py both read `<ID>.cif` from the CURRENT WORKING
DIRECTORY, and the 13 CIFs they need were downloaded into a session scratchpad that no longer
exists. `.gitignore` excludes `*.pdb` outside targets/, so a clean clone has neither. That made
`outbox/02-pk-three-site-plan.txt`'s claim -- "artifacts, all reproducible from a clean clone" --
false for everything downstream of per_partner.json, which is most of the day-1 census.
playbook s26: "it reproduces" means in a clean clone, not on your machine.
"""

import os
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
STRUCTURES = os.path.join(HERE, "structures")
RCSB = "https://files.rcsb.org/download/{}"


def _get(name):
    os.makedirs(STRUCTURES, exist_ok=True)
    path = os.path.join(STRUCTURES, name)
    if not os.path.exists(path):
        url = RCSB.format(name)
        print(f"fetching {url}")
        urllib.request.urlretrieve(url, path)
    if os.path.getsize(path) == 0:
        os.remove(path)
        raise SystemExit(f"REFUSING: {url} returned an empty file; removed the stub rather than "
                         f"letting a later read see zero contacts")
    return path


def cif(pdb_id):
    return _get(f"{pdb_id.upper()}.cif")


def pdb(pdb_id):
    return _get(f"{pdb_id.upper()}.pdb")
