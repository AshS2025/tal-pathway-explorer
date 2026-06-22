# -*- coding: utf-8 -*-
"""
rmg_thermo_server.py
====================

WHAT THIS IS
------------
A long-running "server" script. It loads RMG's thermodynamic database
into memory ONCE (slow — about 10 seconds), then sits waiting for
SMILES strings to be sent in on its standard input. For each SMILES it
receives, it computes the standard enthalpy of formation (Hf, "delta-H
of formation") and prints the value back on its standard output.

WHO RUNS IT
-----------
This script runs INSIDE the `rmg_env` conda environment (Python 2.7).
You will NOT run this script directly. Instead, the Python 3 client
(`src/rmg_thermo.py`) spawns this script as a subprocess and sends it
work via stdin/stdout pipes.

WHY THE TWO-PROCESS DESIGN
--------------------------
Your main project uses Python 3.13. RMG only works on Python 2.7. The
two interpreters cannot share imports, so we use the operating system's
stdin/stdout pipes as the "phone line" between them. The Python 3
side asks questions, this Python 2 side answers them.

INPUT / OUTPUT PROTOCOL
-----------------------
   stdin  : one SMILES string per line, e.g.  "CCO\n"
   stdout : one value per line, either:
              - a number in kJ/mol, e.g.  "-235.30\n"
              - the literal string "NO_THERMO\n" if RMG can't compute it
              - a startup message "READY" or "DB_NOT_FOUND" first

The Python 3 client uses these markers as a handshake to know when
the server is alive and ready to accept queries.

UNITS
-----
RMG returns enthalpies in Joules per mole (J/mol) internally.
This server converts to kJ/mol before printing — please make sure
DORAnet's `max_rxn_thermo_change` threshold is in kJ/mol too when
plugging in this calculator.
"""

import os
import sys

# ----------------------------------------------------------------------
# STEP 1: try to import RMG. If this fails the env is wrong.
# ----------------------------------------------------------------------
try:
    from rmgpy.molecule import Molecule
    from rmgpy.species import Species
    from rmgpy.data.thermo import ThermoDatabase
except ImportError as exc:
    # Print to stderr so the client can see it even if stdout is being
    # read for the protocol.
    sys.stderr.write("IMPORT_FAILED: " + str(exc) + "\n")
    sys.stderr.write(
        "This script must be run inside the rmg_env conda environment.\n"
    )
    sys.exit(1)


# ----------------------------------------------------------------------
# STEP 2: find RMG's thermo database directory.
# RMG's database is a separate folder of YAML/data files. It is usually
# pointed to by the RMG_DATABASE environment variable, OR sits next to
# the rmgpy install. We try the common locations.
# ----------------------------------------------------------------------
def find_rmg_database():
    """
    Returns the path to RMG's `database/thermo` folder, or None if we
    can't find one.

    Beginner note: "database" here means a folder of files RMG reads
    from disk at startup — NOT a SQL database or anything like that.
    """
    # Highest priority: the RMG_DATABASE environment variable.
    env = os.environ.get("RMG_DATABASE")
    if env and os.path.isdir(os.path.join(env, "thermo")):
        return os.path.join(env, "thermo")

    # Fallback: common conda install layout. The conda package
    # sometimes drops the database under the env's share/ tree.
    import rmgpy
    rmgpy_dir = os.path.dirname(rmgpy.__file__)
    # ...rmg_env/Lib/site-packages/rmgpy  → go up two levels to env root
    env_root = os.path.normpath(os.path.join(rmgpy_dir, "..", "..", ".."))
    candidates = [
        # The conda package ships the database at share/rmgdatabase/
        os.path.join(env_root, "share", "rmgdatabase", "thermo"),
        os.path.join(env_root, "share", "rmgpy", "database", "thermo"),
        os.path.join(env_root, "rmg-database", "input", "thermo"),
        # If user cloned rmg-database in their home directory:
        os.path.expanduser(os.path.join("~", "rmg-database", "input", "thermo")),
    ]
    for c in candidates:
        if os.path.isdir(c):
            return c
    return None


thermo_path = find_rmg_database()
if thermo_path is None:
    # No database found — tell the client and exit. The client's error
    # handling will surface this in a meaningful way.
    sys.stdout.write("DB_NOT_FOUND\n")
    sys.stdout.flush()
    sys.stderr.write(
        "Could not locate RMG's thermo database. Either:\n"
        "  1) Set the RMG_DATABASE environment variable to the path of\n"
        "     your cloned rmg-database repo, OR\n"
        "  2) Clone https://github.com/ReactionMechanismGenerator/RMG-database\n"
        "     into ~\\rmg-database\n"
    )
    sys.exit(2)


# ----------------------------------------------------------------------
# STEP 3: load the thermo database into memory.
# This is the slow step (about 10 seconds). It happens once, when the
# server starts up. After that, individual SMILES queries are fast.
# ----------------------------------------------------------------------
sys.stderr.write("Loading RMG ThermoDatabase from " + thermo_path + "...\n")
db = ThermoDatabase()
db.load(
    thermo_path,
    # Empty list = load all libraries the database knows about. For a
    # smaller / faster startup you could pass specific library names.
    libraries=[],
    depository=False,
)
sys.stderr.write("ThermoDatabase loaded.\n")


# ----------------------------------------------------------------------
# STEP 4: handshake — tell the client we're ready to accept SMILES.
# The client is reading our stdout and waiting for the literal word
# "READY" before it starts sending queries.
# ----------------------------------------------------------------------
sys.stdout.write("READY\n")
sys.stdout.flush()


# ----------------------------------------------------------------------
# STEP 5: the main loop. Read one SMILES per line, compute Hf, write
# the result. Loop until stdin closes (which happens when the Python 3
# client disconnects or exits).
# ----------------------------------------------------------------------
def compute_hf_kjmol(smiles_str):
    """
    Convert a SMILES string into an Hf value in kJ/mol, or return None
    if RMG can't handle this molecule.

    Beginner note: this function is the actual chemistry part. Steps:
      1. Parse the SMILES into an RMG Molecule object.
      2. Wrap that Molecule in a Species (RMG's thermo API works on
         Species, not bare Molecules — because Species can have multiple
         resonance forms, isotopes, etc.).
      3. Ask the ThermoDatabase for the thermo data of this Species.
         RMG tries: exact database match first, then group-additivity
         (Joback-style) estimation if no exact match.
      4. Pull out the enthalpy at 298 K (room temperature) — the
         standard way to report Hf.
      5. Convert from J/mol to kJ/mol.
    """
    try:
        molecule = Molecule().fromSMILES(smiles_str)
    except Exception:
        return None  # bad SMILES, RMG won't parse it

    try:
        species = Species(molecule=[molecule])
        # getThermoData returns a NASAPolynomial or Wilhoit object that
        # contains the thermodynamic functions. We just need Hf(298).
        thermo = db.getThermoData(species)
        if thermo is None:
            return None
        hf_joules_per_mol = thermo.getEnthalpy(298.0)
        return hf_joules_per_mol / 1000.0     # kJ/mol
    except Exception:
        return None


while True:
    # Read one line from the client. If the client closes the pipe,
    # readline() returns an empty string and we exit the loop.
    line = sys.stdin.readline()
    if not line:
        break

    smiles_input = line.strip()
    if not smiles_input:
        # Blank line — ignore and wait for the next one.
        continue

    value = compute_hf_kjmol(smiles_input)

    if value is None:
        sys.stdout.write("NO_THERMO\n")
    else:
        # Six decimal places is plenty for kJ/mol. Format with sign.
        sys.stdout.write("{0:.6f}\n".format(value))

    # CRITICAL: flush after every reply. Without this, the Python 3
    # client would block forever waiting for output that's stuck in
    # this process's stdout buffer.
    sys.stdout.flush()

# Reaching here means the client disconnected. Clean exit.
sys.stderr.write("Client disconnected, shutting down.\n")
