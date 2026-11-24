"""Draw every figure from the derived files, an svg beside each png."""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgba_array

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
DATA = os.path.join(PROJECT, "data")
RESULTS = os.path.join(PROJECT, "md", "results")
SERIES = RESULTS
INPUTS = os.path.join(PROJECT, "md", "inputs")
FIGS = os.path.join(PROJECT, "figures")

NATIVE = "#117A65"
CONVERTED = "#7D3C98"
CLASS = {"Q": "#117A65", "T": "#C2185B", "Y": "#7D3C98"}
GREY = "#868E96"
PALE = "#CED4DA"
INK = "#212529"
CMAP = LinearSegmentedColormap.from_list("rqty", ["#DEE2E6", "#117A65",
                                                  "#C2185B"])
DIMER_VIEW = (40, -90)
DPI = 400

plt.rcParams.update({
    "font.size": 8, "axes.labelsize": 8, "axes.titlesize": 8.5,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.edgecolor": INK, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": INK, "ytick.color": INK, "savefig.facecolor": "white",
    "figure.facecolor": "white", "axes.facecolor": "white",
    "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none",
})
MISSING = []


def fatal(msg):
    sys.exit("STOP: " + msg)


def load(name, root=DATA):
    p = os.path.join(root, name)
    if not os.path.isfile(p):
        fatal("no %s" % p)
    with open(p) as fh:
        return json.load(fh)


def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    for ext in ("png", "svg"):
        fig.savefig(os.path.join(FIGS, "%s.%s" % (name, ext)),
                    dpi=DPI if ext == "png" else None, facecolor="white",
                    bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("  wrote figures/%s.png and .svg" % name)


def panel_label(ax, letter, dx=-0.055, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9.5,
            fontweight="bold", va="top", ha="left")


def read_tsv(path):
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(head, L.rstrip("\n").split("\t")))
                for L in fh if L.strip()]


def read_xvg(path):
    """-> (x, y) arrays. Lengths stay in the nanometres GROMACS writes."""
    x, y = [], []
    for L in open(path, errors="replace"):
        s = L.strip()
        if not s or s[0] in "#@":
            continue
        p = s.split()
        if len(p) >= 2:
            x.append(float(p[0]))
            y.append(float(p[1]))
    return np.array(x), np.array(y)


def ca_trace(pdb):
    """-> {chain: (residue numbers, N by 3 coordinates)} in file order."""
    out = {}
    for L in open(pdb, errors="replace"):
        if not L.startswith("ATOM") or L[12:16].strip() != "CA":
            continue
        ch = L[21]
        out.setdefault(ch, ([], []))
        out[ch][0].append(int(L[22:26]))
        out[ch][1].append((float(L[30:38]), float(L[38:46]), float(L[46:54])))
    return {c: (np.array(r), np.array(v)) for c, (r, v) in out.items()}


def atoms_of(pdb, chain, resid):
    """-> {atom name: coordinate} for one residue."""
    out = {}
    for L in open(pdb, errors="replace"):
        if not L.startswith("ATOM") or L[21] != chain:
            continue
        if int(L[22:26]) != resid:
            continue
        out[L[12:16].strip()] = np.array([float(L[30:38]), float(L[38:46]),
                                          float(L[46:54])])
    return out


def smooth(points, per=6):
    """Catmull-Rom through the points, so the trace reads as a tube."""
    p = np.asarray(points, float)
    if len(p) < 3:
        return p
    ext = np.vstack([p[0], p, p[-1]])
    t = np.linspace(0, 1, per, endpoint=False)[:, None]
    out = []
    for i in range(len(p) - 1):
        p0, p1, p2, p3 = ext[i], ext[i + 1], ext[i + 2], ext[i + 3]
        out.append(0.5 * (2 * p1 + (-p0 + p2) * t
                          + (2 * p0 - 5 * p1 + 4 * p2 - p3) * t ** 2
                          + (-p0 + 3 * p1 - 3 * p2 + p3) * t ** 3))
    out.append(p[-1][None, :])
    return np.vstack(out)


def camera(elev, azim):
    """An orthographic camera. The rows are screen right, screen up and depth."""
    e, a = np.radians(elev), np.radians(azim)
    eye = np.array([np.cos(e) * np.cos(a), np.cos(e) * np.sin(a), np.sin(e)])
    right = np.cross([0.0, 0.0, 1.0], eye)
    if np.linalg.norm(right) < 1e-9:
        right = np.array([1.0, 0.0, 0.0])
    right = right / np.linalg.norm(right)
    up = np.cross(eye, right)
    return np.stack([right, up / np.linalg.norm(up), eye])


def project(pts, cam):
    """-> (n, 3): screen x, screen y and depth, which grows toward the eye."""
    return np.atleast_2d(np.asarray(pts, float)) @ cam.T


def fit(ax, pts, pad=0.0):
    """Frame a projected drawing so that it fills the panel at equal scale."""
    pos = ax.get_position()
    fw, fh = ax.figure.get_size_inches()
    shape = (pos.width * fw) / (pos.height * fh)
    xy = np.asarray(pts, float)[:, :2]
    lo, hi = xy.min(0) - pad, xy.max(0) + pad
    mid = (lo + hi) / 2.0
    w, h = hi - lo
    if w / h < shape:
        w = h * shape
    else:
        h = w / shape
    ax.set_xlim(mid[0] - w / 2, mid[0] + w / 2)
    ax.set_ylim(mid[1] - h / 2, mid[1] + h / 2)
    ax.set_axis_off()


class Trace:
    """Backbone segments gathered across chains and painted back to front."""

    def __init__(self):
        self.segs, self.cols, self.depth = [], [], []

    def add(self, pts, colours):
        seg = np.stack([pts[:-1, :2], pts[1:, :2]], axis=1)
        col = to_rgba_array(colours)
        if len(col) == 1:
            col = np.repeat(col, len(seg), axis=0)
        self.segs.append(seg)
        self.cols.append(col)
        self.depth.append((pts[:-1, 2] + pts[1:, 2]) / 2.0)

    def draw(self, ax, lw=2.4, zorder=2):
        order = np.argsort(np.concatenate(self.depth))
        segs = np.concatenate(self.segs)[order]
        cols = np.concatenate(self.cols)[order]
        ax.add_collection(LineCollection(segs, colors=cols, linewidths=lw,
                                         capstyle="round", zorder=zorder))


def sticks(ax, atoms, names, colour, cam, lw=1.6):
    """Bonds between the named atoms of one residue, by distance."""
    keys = [n for n in names if n in atoms]
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            if np.linalg.norm(atoms[a] - atoms[b]) < 1.75:
                p = project([atoms[a], atoms[b]], cam)
                ax.plot(p[:, 0], p[:, 1], color=colour, lw=lw,
                        solid_capstyle="round", zorder=4)


def figure_design(V, uni):
    lead = V["design"]["lead"]
    rec = uni["uniprot"][lead]
    seq, helices = rec["seq"], rec["helix"]
    positions = V["design"]["lead_positions"]
    in_iface = V["interface"]["lead_iface"]["conv_in_iface"]
    last = V["interface"]["lead_iface"]["last"]

    fig = plt.figure(figsize=(7.1, 4.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1.0], width_ratios=[1, 1],
                          hspace=0.72, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    for a, b in helices:
        ax.add_patch(plt.Rectangle((a - 0.5, -0.32), b - a + 1, 0.64,
                                   facecolor="#DEE2E6", edgecolor="none",
                                   zorder=1))
    ax.plot([1, len(seq)], [0, 0], color=GREY, lw=0.7, zorder=0)
    seen = set()
    for p in positions:
        a = seq[p - 1]
        ax.plot([p, p], [0.34, 0.95], lw=0.9, color=CLASS[a], zorder=3,
                label=a if a not in seen else None)
        seen.add(a)
    for k, p in enumerate(sorted(in_iface)):
        ax.plot([p], [1.16], marker="v", ms=4.2, color=INK, zorder=4)
        ax.annotate("%s%d" % (seq[p - 1], p), xy=(p, 1.26 + 0.34 * (k % 2)),
                    ha="center", va="bottom", fontsize=7)
    ax.axvline(last + 0.5, color=GREY, lw=0.7, ls=(0, (3, 2)), zorder=2)
    ax.annotate("crystal ends at %d" % last, xy=(last - 2, -0.70), ha="right",
                va="bottom", fontsize=6.8, color=GREY)
    ax.set_xlim(0, len(seq) + 1)
    ax.set_ylim(-0.80, 1.95)
    ax.set_yticks([])
    ax.set_xlabel("residue of %s" % V["design"]["lead_acc"])
    ax.spines["left"].set_visible(False)
    handles = [plt.Line2D([], [], color=CLASS[a], lw=1.4,
                          label={"Q": "Q to L", "T": "T to V",
                                 "Y": "Y to F"}[a]) for a in "QTY"]
    handles.append(plt.Line2D([], [], color="#DEE2E6", lw=5,
                              label="annotated helix"))
    ax.legend(handles=handles, frameon=False, ncol=4, loc="lower center",
              bbox_to_anchor=(0.5, 1.02), handlelength=1.4, columnspacing=1.2)
    panel_label(ax, "a", dx=-0.035, dy=1.22)

    rows = V["design"]["rows"]
    genes = [r["gene"] for r in rows]
    y = np.arange(len(rows))[::-1]
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.barh(y, [r["conv_n"] for r in rows], color=PALE, height=0.62,
             label="outside the interface")
    ax2.barh(y, [len(V["interface"]["projection"][g]["conv_in_iface"])
                 for g in genes], color=CONVERTED, height=0.62,
             label="in the interface")
    for yy, r in zip(y, rows):
        ax2.annotate("%d" % r["conv_n"], xy=(r["conv_n"] + 0.7, yy),
                     va="center", fontsize=7)
    ax2.set_yticks(y, genes)
    ax2.set_xlim(0, 34)
    ax2.set_xlabel("converted positions")
    ax2.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0),
               ncol=2, handlelength=1.1, handletextpad=0.4, columnspacing=1.0)
    panel_label(ax2, "b", dx=-0.26, dy=1.20)

    ax3 = fig.add_subplot(gs[1, 1])
    for yy, r in zip(y, rows):
        ax3.plot([r["gravy_native"], r["gravy_conv"]], [yy, yy], color=GREY,
                 lw=0.9, zorder=1)
    ax3.scatter([r["gravy_native"] for r in rows], y, s=16, color=NATIVE,
                zorder=2, label="native")
    ax3.scatter([r["gravy_conv"] for r in rows], y, s=16, color=CONVERTED,
                zorder=2, label="converted")
    ax3.axvline(0, color=INK, lw=0.7, ls=(0, (3, 2)))
    ax3.set_yticks(y, ["" for _ in genes])
    ax3.set_xlabel("GRAVY, Kyte and Doolittle")
    ax3.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0),
               ncol=2, handletextpad=0.3, columnspacing=1.0)
    panel_label(ax3, "c", dx=-0.10, dy=1.20)
    save(fig, "fig_design")


def draw_dimer(ax, V, pdb, tr):
    """The dimer as a traced fold, with the converted positions marked."""
    cam = camera(*DIMER_VIEW)
    chains = sorted(tr)
    converted = set(V["design"]["lead_positions"])
    iface = sorted(V["interface"]["lead_iface"]["conv_in_iface"])
    fold, drawn, marks, rings = Trace(), [], [], []
    for c, col in zip(chains, (NATIVE, CONVERTED)):
        resid, xyz = tr[c]
        line = project(smooth(xyz), cam)
        fold.add(line, col)
        drawn.append(line)
        q = project(xyz, cam)
        marks.append(q[[i for i, r in enumerate(resid) if r in converted]])
        if c == chains[0]:
            rings.append(q[[list(resid).index(r) for r in iface]])
    fold.draw(ax, lw=2.8)
    m, r = np.vstack(marks), np.vstack(rings)
    back = np.array([np.any(np.hypot(*(np.vstack(drawn)[:, :2] - p[:2]).T)
                            < 2.2) and np.any(
                        (np.hypot(*(np.vstack(drawn)[:, :2] - p[:2]).T) < 2.2)
                        & (np.vstack(drawn)[:, 2] > p[2] + 3.0)) for p in m])
    for sel, z in ((~back, 5), (back, 1)):
        ax.scatter(m[sel, 0], m[sel, 1], s=20, color=CLASS["T"],
                   edgecolor="white", linewidth=0.5, zorder=z)
    ax.scatter(r[:, 0], r[:, 1], s=74, facecolor="none", edgecolor=INK,
               linewidth=1.1, zorder=6)
    fit(ax, project(np.vstack([tr[c][1] for c in chains]), cam), pad=1.4)
    handles = [plt.Line2D([], [], color=NATIVE, lw=2, label="subunit A"),
               plt.Line2D([], [], color=CONVERTED, lw=2, label="subunit B"),
               plt.Line2D([], [], marker="o", ls="", mfc=CLASS["T"],
                          mec="white", label="converted position"),
               plt.Line2D([], [], marker="o", ls="", mfc="none", mec=INK,
                          label="ringed: positions "
                                + ", ".join("%d" % p for p in iface))]
    ax.legend(handles=handles, frameon=False, loc="lower left",
              bbox_to_anchor=(-0.02, 1.02), ncol=2, handletextpad=0.3,
              columnspacing=1.0, labelspacing=0.35, handlelength=1.3)


def draw_closeup(ax2, V, pdb, tr):
    """The interface tyrosine and the two partners its hydroxyl reaches."""
    chains = sorted(tr)
    ytyr = max(int(k) for k in V["interface"]["lead_iface"]["polar"])
    pairs = V["interface"]["lead_iface"]["polar"][str(ytyr)]
    best = None
    for a in chains:
        for b in chains:
            if a == b:
                continue
            ty = atoms_of(pdb, a, ytyr)
            if "OH" not in ty:
                continue
            d = []
            for c in pairs:
                num = int("".join(ch for ch in c["partner"] if ch.isdigit()))
                pa = atoms_of(pdb, b, num)
                if c["partner_atom"] not in pa:
                    d = None
                    break
                d.append(np.linalg.norm(ty["OH"] - pa[c["partner_atom"]]))
            if d and (best is None or max(d) < best[0]):
                best = (max(d), a, b)
    if best is None:
        MISSING.append("the tyrosine contact was not found in the structure")
        return None
    _, ca, cb = best
    tint = dict(zip(chains, (NATIVE, CONVERTED)))
    ty = atoms_of(pdb, ca, ytyr)
    ring = ["CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"]
    shown = [ty[a] for a in ring if a in ty]
    partners = []
    for c in pairs:
        num = int("".join(ch for ch in c["partner"] if ch.isdigit()))
        pa = atoms_of(pdb, cb, num)
        partners.append((c, num, pa))
        shown += [pa[n] for n in pa if n not in ("N", "C", "O")]

    nrm = np.cross(ty["CD1"] - ty["CG"], ty["CD2"] - ty["CG"])
    nrm = nrm / np.linalg.norm(nrm)
    cam = camera(float(np.degrees(np.arcsin(np.clip(nrm[2], -1, 1)))),
                 float(np.degrees(np.arctan2(nrm[1], nrm[0]))))
    extent = [project(np.vstack(shown), cam)]

    sticks(ax2, ty, ring, tint[ca], cam, lw=2.4)
    oh = project(ty["OH"], cam)[0]
    ax2.scatter([oh[0]], [oh[1]], s=34, color=CLASS["T"], edgecolor="white",
                linewidth=0.5, zorder=6)
    stem = ty["CG"] - ty["CZ"]
    stem = stem / max(np.linalg.norm(stem), 1e-6)
    at = project(ty["CB"] + 1.9 * stem, cam)[0]
    ax2.text(at[0], at[1], "Tyr%d" % ytyr, fontsize=7.5, color=tint[ca],
             ha="center", va="center")
    extent.append(at[None, :])

    for k, (c, num, pa) in enumerate(partners):
        names = [n for n in pa if n not in ("N", "C", "O")]
        sticks(ax2, pa, names, tint[cb], cam, lw=1.8)
        q = project(pa[c["partner_atom"]], cam)[0]
        ax2.plot([oh[0], q[0]], [oh[1], q[1]], color=INK, lw=1.1,
                 ls=(0, (2, 2)), zorder=5)
        d = q[:2] - oh[:2]
        away = np.array([-d[1], d[0]]) / max(np.linalg.norm(d), 1e-9)
        mid = (oh[:2] + (0.44 + 0.13 * k) * d
               + (0.62 if k == 0 else -0.62) * away)
        ax2.text(mid[0], mid[1], "%.2f" % c["distance"], fontsize=7, color=INK,
                 ha="center", va="center")
        extent.append(np.append(mid, 0.0)[None, :])
        far = max((pa[t] for t in names),
                  key=lambda v: np.linalg.norm(v - ty["OH"]))
        out = far - pa[c["partner_atom"]]
        at = project(far + 2.4 * out / max(np.linalg.norm(out), 1e-6), cam)[0]
        ax2.text(at[0], at[1], "%s%d" % (c["partner"][:3].title(), num),
                 fontsize=7.5, color=tint[cb], ha="center", va="center")
        extent.append(at[None, :])

    fit(ax2, np.vstack(extent), pad=0.9)
    return ca, cb


def figure_interface(V):
    pdb = os.path.join(INPUTS, "5IQP_dimer_native.pdb")
    if not os.path.isfile(pdb):
        MISSING.append("fig_interface: no %s" % pdb)
        return
    tr = ca_trace(pdb)
    fig = plt.figure(figsize=(7.1, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.34,
                          wspace=0.30, left=0.09, right=0.99, top=0.885,
                          bottom=0.075)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[0, 1])
    draw_dimer(axa, V, pdb, tr)
    pair = draw_closeup(axb, V, pdb, tr)
    if pair is None:
        plt.close(fig)
        return
    axb.annotate("subunit %s to subunit %s, angstrom" % pair, xy=(0.0, 1.02),
                 xycoords="axes fraction", va="bottom", fontsize=7, color=GREY)
    panel_label(axa, "a", dx=-0.04, dy=1.26)
    panel_label(axb, "b", dx=-0.04, dy=1.26)

    rows = V["interface"]["rows"]
    labels = ["%s\n%s" % (r["pdb"], r["gene"]) for r in rows]
    x = np.arange(len(rows))
    ax = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])
    ax.bar(x, [r["buried"] for r in rows], color=PALE, width=0.62,
           label="buried per subunit")
    ax.bar(x, [r["area_on_conv"] for r in rows], color=CONVERTED, width=0.62,
           label="on the converted positions")
    for xx, r in zip(x, rows):
        ax.annotate("%.0f%%" % r["share_pct"], xy=(xx, r["area_on_conv"] + 28),
                    ha="center", fontsize=7, color=CONVERTED)
        ax.annotate("%d residues" % r["n_iface"], xy=(xx, r["buried"] + 28),
                    ha="center", fontsize=6.8, color=GREY)
    ax.set_xticks(x, labels)
    ax.set_ylabel("buried surface, $\\mathrm{\\AA}^2$")
    ax.set_ylim(0, 1250)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.02),
              ncol=1, handlelength=1.2, labelspacing=0.35)
    panel_label(ax, "c", dx=-0.19, dy=1.26)

    for xx, r in zip(x, rows):
        pol = r["polar"]
        y = max(int(k) for k in pol)
        for c in pol[str(y)]:
            if c["atom"] != "OH":
                continue
            mark = "o" if c["partner"].startswith("ASP") else "s"
            ax2.scatter([xx], [c["distance"]], s=26, marker=mark,
                        facecolor="white", edgecolor=NATIVE, zorder=3,
                        linewidth=1.1)
        ax2.annotate("Tyr%d" % y, xy=(xx, 3.33), ha="center", fontsize=6.8,
                     color=GREY)
    ax2.axhline(3.2, color=GREY, lw=0.7, ls=(0, (3, 2)))
    ax2.annotate("hydrogen bond cutoff", xy=(len(rows) - 0.55, 3.21),
                 ha="right", va="bottom", fontsize=6.5, color=GREY)
    ax2.set_xticks(x, labels)
    ax2.set_xlim(-0.6, len(rows) - 0.4)
    ax2.set_ylim(2.5, 3.45)
    ax2.set_ylabel("hydroxyl to partner, $\\mathrm{\\AA}$")
    handles = [plt.Line2D([], [], marker="o", ls="", mfc="white", mec=NATIVE,
                          label="to aspartate"),
               plt.Line2D([], [], marker="s", ls="", mfc="white", mec=NATIVE,
                          label="to arginine")]
    ax2.legend(handles=handles, frameon=False, loc="lower left",
               bbox_to_anchor=(0.0, 1.02), ncol=2, handletextpad=0.2,
               columnspacing=1.0)
    panel_label(ax2, "d", dx=-0.19, dy=1.26)
    save(fig, "fig_interface")


def figure_criterion(V):
    need = ["sasa_dimer.xvg", "sasa_sub1.xvg", "sasa_sub2.xvg",
            "mindist_dimer.xvg", "numcont_dimer.xvg"]
    have = {}
    for b in ("native", "converted"):
        d = os.path.join(SERIES, b)
        gone = [f for f in need if not os.path.isfile(os.path.join(d, f))]
        if gone:
            MISSING.append("fig_criterion: md/results/%s needs %s"
                           % (b, ", ".join(gone)))
            return
        have[b] = d
    ns = V["md"]["ns"]["native"]
    fig, ax = plt.subplots(2, 2, figsize=(7.1, 4.3),
                           gridspec_kw={"wspace": 0.34, "hspace": 0.52})
    axes = [ax[0, 0], ax[0, 1], ax[1, 0], ax[1, 1]]
    for b, col in (("native", NATIVE), ("converted", CONVERTED)):
        t, ab = read_xvg(os.path.join(have[b], "sasa_dimer.xvg"))
        _, a1 = read_xvg(os.path.join(have[b], "sasa_sub1.xvg"))
        _, a2 = read_xvg(os.path.join(have[b], "sasa_sub2.xvg"))
        n = min(len(t), len(ab), len(a1), len(a2))
        buried = (a1[:n] + a2[:n] - ab[:n]) * 100.0 / 2.0
        axes[0].plot(t[:n] / 1000.0, buried, lw=0.6, color=col, label=b)
        tc, nc = read_xvg(os.path.join(have[b], "numcont_dimer.xvg"))
        axes[1].plot(tc / 1000.0, nc, lw=0.6, color=col)
        td, md = read_xvg(os.path.join(have[b], "mindist_dimer.xvg"))
        axes[2].plot(td / 1000.0, md * 10.0, lw=0.6, color=col)
    axes[1].set_ylabel("inter-subunit contacts")
    hold = V["md"]["thresholds"]["hold"] * V["md"]["buried"]["native"]
    axes[0].axhline(hold, color=GREY, lw=0.8, ls=(0, (3, 2)))
    axes[0].annotate("declared hold threshold", xy=(2, hold + 20),
                     fontsize=6.5, color=GREY)
    axes[0].set_ylabel("buried per subunit, $\\mathrm{\\AA}^2$")
    axes[0].set_ylim(0, None)
    axes[2].set_ylabel("shortest subunit distance, $\\mathrm{\\AA}$")

    pos = sorted(V["interface"]["lead_iface"]["conv_in_iface"])
    width = 0.34
    xs = np.arange(len(pos))
    for k, (b, col) in enumerate((("native", NATIVE), ("converted", CONVERTED))):
        vals = [np.mean([V["md"]["per_position"]["%d_sub%d" % (p, s)][b + "_A"]
                         for s in (1, 2)]) for p in pos]
        axes[3].bar(xs + (k - 0.5) * width, vals, width=width, color=col,
                    label=b)
    axes[3].set_xticks(xs, ["%d" % p for p in pos])
    axes[3].set_xlabel("converted interface position")
    axes[3].set_ylabel("distance to partner, $\\mathrm{\\AA}$")
    axes[3].set_ylim(0, 3.0)
    for a in (axes[0], axes[1], axes[2]):
        a.set_xlabel("time, ns")
        a.set_xlim(0, ns)
    for a, letter in zip(axes, "abcd"):
        panel_label(a, letter, dx=-0.26)
    fig.legend(*axes[0].get_legend_handles_labels(), frameon=False,
               loc="center", bbox_to_anchor=(0.5, 0.955), ncol=2,
               handlelength=1.4, columnspacing=1.4)
    save(fig, "fig_criterion")


def figure_fluctuation(V):
    p = os.path.join(RESULTS, "rmsf_table.tsv")
    if not os.path.isfile(p):
        MISSING.append("fig_fluctuation: no %s" % p)
        return
    rows = read_tsv(p)
    chains = sorted({r["chain"] for r in rows})
    fig, axes = plt.subplots(len(chains), 1, figsize=(7.1, 1.7 * len(chains)),
                             sharex=True, gridspec_kw={"hspace": 0.35})
    axes = np.atleast_1d(axes)
    for ax, ch in zip(axes, chains):
        sel = [r for r in rows if r["chain"] == ch]
        x = [int(r["residue"]) for r in sel]
        for key, col, lab in (("rmsf_native_A", NATIVE, "native"),
                              ("rmsf_converted_A", CONVERTED, "converted")):
            ax.plot(x, [float(r[key]) for r in sel], lw=0.8, color=col,
                    label=lab)
        top = max(float(r["rmsf_converted_A"]) for r in sel)
        for r in sel:
            if r["converted_position"] == "1":
                ax.plot([int(r["residue"])], [top * 1.06], marker="v", ms=3.4,
                        color=GREY, clip_on=False)
        ax.set_ylabel("C$\\alpha$ fluctuation, $\\mathrm{\\AA}$")
        ax.set_title("subunit %s" % ch, fontsize=8, loc="left")
    axes[0].legend(frameon=False, ncol=2, loc="upper center",
                   bbox_to_anchor=(0.5, 1.52), handlelength=1.4)
    axes[0].annotate("markers are the converted positions", xy=(0.0, 1.34),
                     xycoords="axes fraction", fontsize=7, color=GREY)
    axes[-1].set_xlabel("residue of %s" % V["design"]["lead_acc"])
    save(fig, "fig_fluctuation")


def figure_fluctuation_map(V):
    p = os.path.join(RESULTS, "rmsf_table.tsv")
    pdb = os.path.join(INPUTS, "5IQP_dimer_native.pdb")
    if not (os.path.isfile(p) and os.path.isfile(pdb)):
        MISSING.append("fig_fluctuation_map: needs rmsf_table.tsv and the dimer")
        return
    rows = read_tsv(p)
    val = {b: {(r["chain"], int(r["residue"])): float(r["rmsf_%s_A" % b])
               for r in rows} for b in ("native", "converted")}
    tr = ca_trace(pdb)
    lo = min(min(v.values()) for v in val.values())
    hi = max(max(v.values()) for v in val.values())
    iface = V["interface"]["lead_iface"]["conv_in_iface"]

    cam = camera(*DIMER_VIEW)
    frame = project(np.vstack([tr[c][1] for c in tr]), cam)
    pad = 1.4
    left, right, top, bottom, wspace = 0.012, 0.988, 0.935, 0.20, 0.04
    pw = (right - left) * 7.1 / (2.0 + wspace)
    shape = (np.ptp(frame[:, 1]) + 2 * pad) / (np.ptp(frame[:, 0]) + 2 * pad)
    fig = plt.figure(figsize=(7.1, pw * shape / (top - bottom)))
    gs = fig.add_gridspec(1, 2, left=left, right=right, top=top,
                          bottom=bottom, wspace=wspace)
    for k, b in enumerate(("native", "converted")):
        ax = fig.add_subplot(gs[0, k])
        fold = Trace()
        for c in sorted(tr):
            resid, xyz = tr[c]
            v = np.array([val[b].get((c, int(r)), np.nan) for r in resid])
            keep = ~np.isnan(v)
            pts = smooth(xyz[keep], per=6)
            w = np.interp(np.linspace(0, 1, len(pts)),
                          np.linspace(0, 1, int(keep.sum())), v[keep])
            fold.add(project(pts, cam), CMAP((np.clip(w, lo, hi) - lo)
                                             / (hi - lo))[:-1])
        fold.draw(ax, lw=2.4)
        resid, xyz = tr[sorted(tr)[0]]
        q = project(xyz, cam)[[list(resid).index(r) for r in sorted(iface)]]
        ax.scatter(q[:, 0], q[:, 1], s=34, facecolor="none", edgecolor=INK,
                   linewidth=0.9, zorder=6)
        fit(ax, frame, pad=pad)
        ax.set_title(b, fontsize=8.5, pad=3)
    cax = fig.add_axes([0.37, 0.095, 0.26, 0.026])
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(vmin=lo, vmax=hi))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("C$\\alpha$ fluctuation, $\\mathrm{\\AA}$", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    save(fig, "fig_fluctuation_map")


def figure_conservation(V):
    path = os.path.join(DATA, "conservation_P27348.tsv")
    if not os.path.isfile(path):
        fatal("no %s" % path)
    rows = read_tsv(path)
    pos = np.array([int(r["i"]) for r in rows])
    cons = np.array([float(r["conservation"]) for r in rows])
    converted = set(V["design"]["lead_positions"])
    last = V["interface"]["lead_iface"]["last"]
    in_iface = V["interface"]["lead_iface"]["conv_in_iface"]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(7.1, 2.7),
                                  gridspec_kw={"wspace": 0.26,
                                               "width_ratios": [1.35, 1]})
    ax.plot(pos, cons, lw=0.8, color=GREY)
    sel = [i for i, p in enumerate(pos) if p in converted]
    ax.scatter(pos[sel], cons[sel], s=13, color=CONVERTED, zorder=3,
               label="converted position")
    for p in in_iface:
        ax.annotate("%d" % p, xy=(p, cons[p - 1] + 0.035), ha="center",
                    fontsize=6.8, color=INK)
    ax.axhline(V["conservation"]["mean_all"], color=NATIVE, lw=0.8,
               ls=(0, (3, 2)), label="mean over the chain")
    ax.set_xlabel("residue of %s" % V["design"]["lead_acc"])
    ax.set_ylabel("conservation")
    ax.set_xlim(0, len(pos) + 1)
    ax.set_ylim(0, 0.82)
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.0, 1.0),
              ncol=2, handletextpad=0.4)
    panel_label(ax, "a", dx=-0.07, dy=1.20)

    a = [c for p, c in zip(pos, cons) if p <= last and p in converted]
    b = [c for p, c in zip(pos, cons) if p <= last and p not in converted]
    parts = ax2.violinplot([b, a], positions=[0, 1], widths=0.7,
                           showextrema=False, showmedians=False)
    for body, col in zip(parts["bodies"], [GREY, CONVERTED]):
        body.set_facecolor(col)
        body.set_alpha(0.25)
        body.set_edgecolor("none")
    for k, (vals, col) in enumerate([(b, GREY), (a, CONVERTED)]):
        ax2.scatter(np.full(len(vals), k) + np.linspace(-0.11, 0.11, len(vals)),
                    vals, s=3.5, color=col, alpha=0.75, zorder=3)
        ax2.plot([k - 0.26, k + 0.26], [np.mean(vals)] * 2, color=INK, lw=1.3,
                 zorder=4)
    ax2.set_xticks([0, 1], ["not converted\nn = %d" % len(b),
                            "converted\nn = %d" % len(a)])
    ax2.set_ylabel("conservation")
    ax2.set_ylim(0, 0.82)
    ax2.annotate("difference of means %+.3f, permutation p = %.2f"
                 % (V["conservation"]["difference"], V["conservation"]["p"]),
                 xy=(0.0, 1.02), va="bottom", fontsize=6.5,
                 xycoords="axes fraction")
    panel_label(ax2, "b", dx=-0.17, dy=1.20)
    save(fig, "fig_conservation")


def main():
    V = load("derived_values.json")
    uni = load("paralogue_uniprot.json")
    print("building the figures")
    figure_design(V, uni)
    figure_interface(V)
    figure_criterion(V)
    figure_fluctuation(V)
    figure_fluctuation_map(V)
    figure_conservation(V)
    for m in MISSING:
        print("  skipped, " + m)


if __name__ == "__main__":
    main()
