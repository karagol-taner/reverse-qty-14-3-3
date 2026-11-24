"""Cartoon renders of the deposited dimer for the graphical abstract."""

import json
import os
import sys

import numpy as np
import pymol2

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
PDB = os.path.join(PROJECT, "md", "inputs", "5IQP_dimer_native.pdb")
FIGS = os.path.join(PROJECT, "figures")

NATIVE = "#117A65"
CONVERTED = "#7D3C98"
MARK = "#C2185B"
DEEP = "#7A0E36"
WIDTH, HEIGHT = 2400, 1700


def fatal(msg):
    raise SystemExit("STOP: " + msg)


def rgb(h):
    h = h.lstrip("#")
    return [int(h[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]


def unit(v):
    n = float(np.linalg.norm(v))
    if n < 1e-9:
        fatal("a frame vector came out as zero, so the orientation is not "
              "defined")
    return v / n


def frame(cmd, chains, groove):
    """Right is the subunit axis and up is the groove face."""
    ca = {c: cmd.get_coords("d and chain %s and name CA" % c) for c in chains}
    for c, v in ca.items():
        if v is None or len(v) < 20:
            fatal("chain %s carries too few CA atoms to orient on" % c)
    cen = np.vstack(list(ca.values())).mean(0)
    x = unit(ca[chains[1]].mean(0) - ca[chains[0]].mean(0))
    sel = "d and name CA and resi " + "+".join(str(r) for r in groove)
    gv = cmd.get_coords(sel)
    if gv is None or len(gv) < 2 * len(groove) - 2:
        fatal("the groove residues did not select in both chains, so up is "
              "not defined")
    up = unit(gv.mean(0) - cen)
    up = unit(up - x * float(up @ x))
    z = unit(np.cross(x, up))
    y = unit(np.cross(z, x))
    R = np.vstack([x, y, z])
    if float((gv.mean(0) - cen) @ y) <= 0:
        fatal("the groove face did not resolve above the centre, so the "
              "orientation would put the bilayer on the wrong side")
    return R, cen


def main():
    if not os.path.isfile(PDB):
        fatal("no %s" % PDB)
    V = json.load(open(os.path.join(PROJECT, "data",
                                    "derived_values.json")))
    conv = sorted(V["design"]["lead_positions"])
    iface = sorted(V["interface"]["lead_iface"]["conv_in_iface"])
    if not conv or not iface:
        fatal("no converted or interface positions in the derived file")

    with pymol2.PyMOL() as p:
        cmd = p.cmd
        cmd.load(PDB, "d")
        chains = sorted(cmd.get_chains("d"))
        if len(chains) != 2:
            fatal("the build carries %d chains, the panel shows a dimer"
                  % len(chains))

        for name, h in (("c_native", NATIVE), ("c_conv", CONVERTED),
                        ("c_mark", MARK), ("c_deep", DEEP)):
            cmd.set_color(name, rgb(h))

        cmd.set("auto_zoom", 0)
        cmd.hide("everything")
        cmd.show("cartoon", "d")
        cmd.color("c_native", "d and chain %s" % chains[0])
        cmd.color("c_conv", "d and chain %s" % chains[1])

        cmd.set("ray_opaque_background", 0)
        cmd.set("antialias", 2)
        cmd.set("ray_trace_mode", 1)
        cmd.set("ray_trace_gain", 0.10)
        cmd.set("ray_trace_color", "black")
        cmd.set("cartoon_fancy_helices", 1)
        cmd.set("cartoon_transparency", 0)
        cmd.set("specular", 0.12)
        cmd.set("ambient", 0.38)
        cmd.set("direct", 0.52)
        cmd.set("depth_cue", 0)
        cmd.set("ray_shadows", 0)
        cmd.set("orthoscopic", 1)
        cmd.set("sphere_quality", 3)

        groove = sorted(int(k) for k in V["conservation"]["groove"])
        if len(groove) < 3:
            fatal("the derived file names %d groove residues, too few to "
                  "orient on" % len(groove))
        R, cen = frame(cmd, chains, groove)
        view = list(R.flatten()) + [0.0, 0.0, -180.0] + list(cen) \
            + [100.0, 260.0, -20.0]
        cmd.set_view(view)
        cmd.zoom("d", 3.5)

        out = os.path.join(FIGS, "ga_native.png")
        cmd.png(out, width=WIDTH, height=HEIGHT, dpi=400, ray=1)
        if not os.path.isfile(out):
            fatal("PyMOL did not write %s" % out)
        print("wrote figures/ga_native.png")
        before = list(cmd.get_view())

        sel = "d and name CA and resi " + "+".join(str(r) for r in conv)
        if cmd.count_atoms(sel) == 0:
            fatal("none of the converted positions selected in the build")
        isel = "d and name CA and resi " + "+".join(str(r) for r in iface)
        for name, where, colour, scale in (("marks", sel, "c_mark", 0.60),
                                           ("imarks", isel, "c_deep", 1.35)):
            cmd.create(name, where)
            if cmd.count_atoms(name) == 0:
                fatal("the %s object came out empty" % name)
            cmd.hide("everything", name)
            cmd.show("spheres", name)
            cmd.color(colour, name)
            cmd.set("sphere_scale", scale, name)

        after = list(cmd.get_view())
        if max(abs(a - b) for a, b in zip(before, after)) > 1e-4:
            fatal("the camera moved between the two renders, so the panels "
                  "would not superimpose")
        out = os.path.join(FIGS, "ga_converted.png")
        cmd.png(out, width=WIDTH, height=HEIGHT, dpi=400, ray=1)
        if not os.path.isfile(out):
            fatal("PyMOL did not write %s" % out)
        print("wrote figures/ga_converted.png")
    return 0


if __name__ == "__main__":
    sys.exit(main())
