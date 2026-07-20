"""
scripts/dora_xgb_server.py — runs INSIDE the `dora_xgb` conda env.

DORA-XGB's dependencies (xgboost 1.6.2, mordred, old numpy) can't live in
the main Python 3.13 env, so — exactly like the RMG thermo bridge — it
runs as a subprocess in its own env and talks over stdin/stdout.

Protocol (one line each):
    stdout: "READY"                once the model has loaded
            "LOAD_FAILED"          if the model couldn't be built
    stdin:  "<reactants>>products>"  a reaction SMILES to score
    stdout: "0.975"                its feasibility score in [0, 1]
            "NO_SCORE"             if that reaction couldn't be scored

Any library chatter (warnings, load messages) is forced to stderr so it
never corrupts the one-line-per-message protocol on stdout.
"""
import os
import sys

# Keep stdout pristine for the protocol: send everything the libraries
# print to stderr instead, and write protocol lines to the real stdout.
_real_stdout = sys.stdout
sys.stdout = sys.stderr


def emit(msg):
    _real_stdout.write(msg + "\n")
    _real_stdout.flush()


try:
    from DORA_XGB import DORA_XGB
    import DORA_XGB as _pkg

    _cof = os.path.join(
        os.path.dirname(_pkg.__file__),
        "cofactors",
        "expanded_cofactors_no_stereochem.tsv",
    )
    _model = DORA_XGB.feasibility_classifier(
        cofactor_positioning="by_descending_MW",
        max_species=4,
        fp_type="ecfp4",
        nBits=2048,
        model_type="main",
        cofactors_filepath=_cof,
    )
except Exception as e:  # noqa: BLE001
    sys.stderr.write("DORA-XGB load failed: %r\n" % (e,))
    sys.stderr.flush()
    emit("LOAD_FAILED")
    sys.exit(1)

emit("READY")

for _line in sys.stdin:
    _rxn = _line.strip()
    if not _rxn:
        continue
    try:
        _score = float(_model.predict_proba(_rxn))
        emit("%.6f" % _score)
    except Exception as e:  # noqa: BLE001
        sys.stderr.write("score failed for %s: %r\n" % (_rxn, e))
        sys.stderr.flush()
        emit("NO_SCORE")
