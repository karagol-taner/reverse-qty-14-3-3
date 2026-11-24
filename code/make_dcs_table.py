"""Join the conservation column onto the per-residue fluctuation and score it."""

import json
import os
import sys

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
    "HSD": "H", "HSE": "H", "HSP": "H", "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C", "ASH": "D", "GLH": "E", "LYN": "K",
}
REFERENCE = "P27348"


def fatal(msg):
    sys.exit("STOP: " + msg)


def read_tsv(path):
    """-> list of dicts, keyed by the header row."""
    if not os.path.isfile(path):
        fatal("no %s" % path)
    with open(path) as fh:
        lines = [L.rstrip("\n") for L in fh if L.strip()]
    head = lines[0].split("\t")
    out = []
    for n, L in enumerate(lines[1:], start=2):
        parts = L.split("\t")
        if len(parts) != len(head):
            fatal("line %d of %s holds %d fields against %d in the header"
                  % (n, os.path.basename(path), len(parts), len(head)))
        out.append(dict(zip(head, parts)))
    return out


def rescale(raw):
    """0 to 10 across everything given, as the published code rescales it."""
    lo, hi = min(raw), max(raw)
    span = hi - lo
    return [0.0 if span == 0 else (v - lo) / span * 10.0 for v in raw]


def main(project):
    results = os.path.join(project, "md", "results")
    rmsf = read_tsv(os.path.join(results, "rmsf_table.tsv"))
    cons_path = os.path.join(project, "data",
                             "conservation_%s.tsv" % REFERENCE)
    cons = {int(r["i"]): r for r in read_tsv(cons_path)}
    print("%d residues measured, %d positions with conservation"
          % (len(rmsf), len(cons)))

    missing = sorted({int(r["residue"]) for r in rmsf} - set(cons))
    if missing:
        fatal("no conservation for residue %s. The alignment covers %d to %d of "
              "%s, so the build and the reference are not on one numbering."
              % (missing[:5], min(cons), max(cons), REFERENCE))

    bad = []
    for r in rmsf:
        want = cons[int(r["residue"])]["A_i"]
        got = THREE_TO_ONE.get(r["resname_native"].upper())
        if got is None:
            fatal("residue %s of chain %s is %s in the native build and that "
                  "name is not in the table of residue names, so it cannot be "
                  "checked against the alignment"
                  % (r["residue"], r["chain"], r["resname_native"]))
        if got != want:
            bad.append("%s%s is %s in the build and %s in %s"
                       % (r["chain"], r["residue"], got, want, REFERENCE))
    if bad:
        fatal("the build and the alignment disagree at %d of %d positions, so "
              "the numbering has slipped and joining them would score the wrong "
              "residues. First few: %s" % (len(bad), len(rmsf), "; ".join(bad[:5])))
    print("every residue name in the native build matches %s at its number"
          % REFERENCE)

    c = [float(cons[int(r["residue"])]["conservation"]) for r in rmsf]
    out, scores = [], {}
    for b in ("native", "converted"):
        d = [float(r["norm_dynamics_" + b]) for r in rmsf]
        scores["dcs_" + b] = rescale([ci * di for ci, di in zip(c, d)])
        scores["rcs_" + b] = rescale([ci * (1.0 - di) for ci, di in zip(c, d)])

    cols = ["chain", "residue", "A_i", "resname_native", "resname_converted",
            "conservation", "coverage", "norm_dynamics_native",
            "norm_dynamics_converted", "dcs_native_0_10", "dcs_converted_0_10",
            "rcs_native_0_10", "rcs_converted_0_10", "converted_position",
            "interface", "terminal"]
    p = os.path.join(results, "dcs_table.tsv")
    with open(p, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for k, r in enumerate(rmsf):
            cr = cons[int(r["residue"])]
            fh.write("\t".join([
                r["chain"], r["residue"], cr["A_i"], r["resname_native"],
                r["resname_converted"], "%.6f" % c[k], cr["coverage"],
                r["norm_dynamics_native"], r["norm_dynamics_converted"],
                "%.4f" % scores["dcs_native"][k],
                "%.4f" % scores["dcs_converted"][k],
                "%.4f" % scores["rcs_native"][k],
                "%.4f" % scores["rcs_converted"][k],
                r["converted_position"], r["interface"], r["terminal"]]) + "\n")
    print("wrote", os.path.relpath(p, project))

    meta = os.path.join(project, "data", "conservation_%s.json" % REFERENCE)
    depth = json.load(open(meta))["reweighting"]["effective_depth_per_column"]
    conv = [k for k, r in enumerate(rmsf) if r["converted_position"] == "1"]
    for b in ("native", "converted"):
        s = scores["dcs_" + b]
        rank = sorted(range(len(s)), key=lambda k: -s[k])
        print("\n%s build, highest dynamic conserved score" % b)
        for k in rank[:10]:
            r = rmsf[k]
            tag = [t for t, on in (("converted", r["converted_position"] == "1"),
                                   ("interface", r["interface"] == "1"),
                                   ("terminal", r["terminal"] == "1")) if on]
            print("  %s%-4s %-3s  %5.2f   %s"
                  % (r["chain"], r["residue"], cons[int(r["residue"])]["A_i"],
                     s[k], ",".join(tag) or "-"))
        print("  the %d converted positions average %.2f against %.2f over all "
              "%d residues"
              % (len(conv), sum(s[k] for k in conv) / len(conv),
                 sum(s) / len(s), len(s)))
    print("\nalignment depth behind the conservation column: %.1f effective "
          "sequences per position" % depth)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1
         else os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
