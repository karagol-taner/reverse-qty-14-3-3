"""Dimer-interface analysis for the 14-3-3 reverse-QTY series."""
import collections, hashlib, json, os, re, subprocess, sys
import numpy as np
import mdtraj as md
from Bio import Align
from Bio.Align import substitution_matrices

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
CACHE = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "structures", "pdb_cache")
RAW = os.path.join(ROOT, "data", "paralogue_uniprot.json")
RQTY = {"Q": "L", "T": "V", "Y": "F"}
HELICAL = ("H", "G", "I")
CUTOFFS = (3.5, 4.0, 4.5, 5.0, 5.5, 6.0)
DSASA = (1.0, 5.0, 10.0, 20.0, 30.0)
GENES = ("YWHAB", "YWHAQ", "YWHAG", "YWHAH", "YWHAZ", "SFN", "YWHAE")
SYSTEMS = (("2BQ0", "A", "B", "YWHAB", "2"),
           ("8EQ8", "A", "B", "YWHAB", None),
           ("5IQP", "A", "B", "YWHAQ", None),
           ("6F08", "A", "B", "YWHAZ", None))


def fetch(pid, ext):
    path = os.path.join(CACHE, f"{pid}.{ext}")
    if not os.path.exists(path):
        os.makedirs(CACHE, exist_ok=True)
        url = f"https://files.rcsb.org/download/{pid}.{ext}"
        r = subprocess.run(["curl", "-sS", "--fail", "--max-time", "120", "-o", path, url])
        if r.returncode != 0 or not os.path.exists(path):
            raise SystemExit(f"fetch failed: {url}")
    size = os.path.getsize(path)
    if size < 10000:
        raise SystemExit(f"{path} is {size} bytes, which is too small to be a structure")
    return path, hashlib.sha256(open(path, "rb").read()).hexdigest()[:16]


def oper(cif_path, oid):
    txt = open(cif_path).read()
    i = txt.find("_pdbx_struct_oper_list.id")
    if i < 0:
        raise SystemExit(f"no _pdbx_struct_oper_list in {cif_path}")
    blk = txt[i:txt.find("\n#", i)]
    for row in re.split(r"\n(?=\d+ )", blk):
        if not re.match(r"^%s " % re.escape(oid), row):
            continue
        nums = [float(x) for x in re.findall(r"-?\d+\.\d+", row)]
        if len(nums) < 12:
            raise SystemExit(f"operator {oid} in {cif_path} carries {len(nums)} numbers")
        m = np.array(nums[:12]).reshape(3, 4)
        return m[:, :3], m[:, 3]
    raise SystemExit(f"operator {oid} not found in {cif_path}")


def resolution(cif_path):
    txt = open(cif_path).read()
    for key in ("_refine.ls_d_res_high", "_reflns.d_resolution_high",
                "_em_3d_reconstruction.resolution"):
        m = re.search(re.escape(key) + r"\s+(\d+\.?\d*)", txt)
        if m:
            return float(m.group(1))
    raise SystemExit(f"no resolution record in {cif_path}")


def composition(pdb_path):
    het = collections.Counter()
    chains = collections.Counter()
    for ln in open(pdb_path):
        if ln[:6] == "HETATM" and ln[17:20].strip() not in ("HOH", "MSE", "NIY"):
            het[ln[17:20].strip()] += 1
        elif ln[:6] == "ATOM  ":
            chains[ln[21]] += 1
    return dict(het_components=dict(het),
                chains_large=sorted(c for c, n in chains.items() if n > 1000),
                chains_small=sorted(c for c, n in chains.items() if n <= 1000))


def read_pdb(path, chain):
    rows = []
    for ln in open(path):
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        if ln[21] != chain or ln[16] not in " A" or ln[76:78].strip() == "H":
            continue
        rows.append(ln)
    if not rows:
        raise SystemExit(f"chain {chain} carries no atoms in {path}")
    return rows


def rewrite(lines, coords, chain):
    return [l[:21] + chain + l[22:30] + "%8.3f%8.3f%8.3f" % tuple(c) + l[54:]
            for l, c in zip(lines, coords)]


def coords(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])


def dimer_pdb(pid, ch_a, ch_b, op_id, out):
    pdb, sha_p = fetch(pid, "pdb")
    cif, sha_c = fetch(pid, "cif")
    a = read_pdb(pdb, ch_a)
    b = read_pdb(pdb, ch_b)
    xb = coords(b)
    if op_id is not None:
        r, t = oper(cif, op_id)
        xb = xb @ r.T + t
    with open(out, "w") as fh:
        fh.writelines(rewrite(a, coords(a), "A"))
        fh.writelines(rewrite(b, xb, "C"))
        fh.write("END\n")
    return sha_p, sha_c


def series_table(raw):
    out = {}
    for g in GENES:
        rec = raw["uniprot"][g]
        seq, hel = rec["seq"], rec["helix"]
        pos = [p for p, c in enumerate(seq, 1)
               if c in RQTY and any(a <= p <= b for a, b in hel)]
        conv = list(seq)
        for p in pos:
            conv[p - 1] = RQTY[conv[p - 1]]
        kd = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
              "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
              "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
        gravy = lambda s: sum(kd[c] for c in s) / len(s)
        out[g] = dict(acc=rec["acc"], length=len(seq),
                      helix_pct=100.0 * sum(b - a + 1 for a, b in hel) / len(seq),
                      qty=sum(1 for c in seq if c in RQTY), conv=pos,
                      conv_pct=100.0 * len(pos) / len(seq),
                      gravy_native=gravy(seq), gravy_conv=gravy("".join(conv)))
    return out


def analyse(pid, ch_a, ch_b, gene, op_id, raw, tmp, tmp_au=None):
    tmp_au = tmp_au or tmp
    seq = raw["uniprot"][gene]["seq"]
    out_pdb = os.path.join(tmp, f"{pid}_dimer.pdb")
    sha_p, sha_c = dimer_pdb(pid, ch_a, ch_b, op_id, out_pdb)
    t = md.load(out_pdb)
    t = t.atom_slice(t.top.select("protein"))
    cid = {c.chain_id: c.index for c in t.top.chains}
    if set(cid) != {"A", "C"}:
        raise SystemExit(f"{pid}: built dimer carries chains {sorted(cid)}")
    res = list(t.top.chain(cid["A"]).residues)
    bad = [(r.resSeq, r.code) for r in res
           if not (1 <= r.resSeq <= len(seq)) or seq[r.resSeq - 1] != r.code]
    if bad:
        raise SystemExit(f"{pid} chain {ch_a} disagrees with {gene} at {bad[:5]}")
    dssp = dict(zip([r.resSeq for r in res],
                    md.compute_dssp(t.atom_slice(t.top.select(f"chainid {cid['A']}")),
                                    simplified=False)[0]))
    obs = {r.resSeq: r.code for r in res}
    qty = sorted(p for p, a in obs.items() if a in RQTY)
    conv = sorted(p for p in qty if dssp[p] in HELICAL)

    heavy = lambda ch: [a.index for a in t.top.atoms
                        if a.residue.chain.index == cid[ch] and a.element.symbol != "H"]
    ia, ic = heavy("A"), heavy("C")
    xa, xc = t.xyz[0][ia] * 10.0, t.xyz[0][ic] * 10.0
    ra = [t.top.atom(i).residue.resSeq for i in ia]
    d = np.sqrt(((xa[:, None, :] - xc[None, :, :]) ** 2).sum(-1))
    mind = collections.defaultdict(lambda: 1e9)
    for i, u in enumerate(ra):
        m = float(d[i].min())
        if m < mind[u]:
            mind[u] = m
    sd = md.shrake_rupley(t, mode="residue")[0] * 100.0
    ta = t.atom_slice(t.top.select(f"chainid {cid['A']}"))
    sa = md.shrake_rupley(ta, mode="residue")[0] * 100.0
    idx = [r.index for r in t.top.chain(cid["A"]).residues]
    ds = {r.resSeq: float(sa[i]) - float(sd[idx[i]]) for i, r in enumerate(res)}
    buried = sum(v for v in ds.values() if v > 0)
    if buried < 300.0:
        raise SystemExit(f"{pid}: {buried:.0f} A^2 buried, which is a crystal contact "
                         f"and not a dimer interface")
    contact = {c: sorted(u for u, m in mind.items() if m <= c) for c in CUTOFFS}
    sasa_if = {th: sorted(u for u, v in ds.items() if v > th) for th in DSASA}
    iface = sorted(set(contact[5.0]) | set(sasa_if[1.0]))

    polar = {}
    sidechain = {"N", "O"}
    for p in [q for q in conv if q in set(iface)]:
        src = [a for a in t.top.atoms
               if a.residue.chain.index == cid["A"] and a.residue.resSeq == p
               and a.element.symbol in sidechain and a.name not in ("N", "O", "OXT")]
        found = []
        for a in src:
            for b in t.top.atoms:
                if b.residue.chain.index != cid["C"] or b.element.symbol not in sidechain:
                    continue
                dd = float(np.linalg.norm(t.xyz[0][a.index] - t.xyz[0][b.index]) * 10.0)
                if dd <= 4.0:
                    found.append(dict(atom=a.name, partner=f"{b.residue.name}{b.residue.resSeq}",
                                      partner_atom=b.name, distance=round(dd, 2)))
        polar[str(p)] = sorted(found, key=lambda r: r["distance"])

    au_buried = None
    if op_id is not None:
        au = os.path.join(tmp_au, f"{pid}_au.pdb")
        dimer_pdb(pid, ch_a, ch_b, None, au)
        tau = md.load(au)
        tau = tau.atom_slice(tau.top.select("protein"))
        caid = {c.chain_id: c.index for c in tau.top.chains}
        s_pair = md.shrake_rupley(tau, mode="residue")[0] * 100.0
        s_one = md.shrake_rupley(tau.atom_slice(tau.top.select(f"chainid {caid['A']}")),
                                 mode="residue")[0] * 100.0
        idx_au = [r.index for r in tau.top.chain(caid["A"]).residues]
        au_buried = float(sum(max(float(s_one[i]) - float(s_pair[idx_au[i]]), 0.0)
                              for i in range(len(idx_au))))
    pdb_path = os.path.join(CACHE, f"{pid}.pdb")
    cif_path = os.path.join(CACHE, f"{pid}.cif")
    return dict(pdb=pid, gene=gene, sha256_pdb=sha_p, sha256_cif=sha_c,
                resolution=resolution(cif_path), composition=composition(pdb_path),
                observed=len(res), first=res[0].resSeq, last=res[-1].resSeq,
                helix_pct=100.0 * sum(1 for v in dssp.values() if v in HELICAL) / len(dssp),
                beta_pct=100.0 * sum(1 for v in dssp.values() if v in ("E", "B")) / len(dssp),
                qty_observed=len(qty), conv_dssp=conv,
                buried_per_subunit=buried,
                contact={str(k): v for k, v in contact.items()},
                sasa={str(k): v for k, v in sasa_if.items()},
                iface=iface, conv_in_iface=[p for p in conv if p in set(iface)],
                polar_contacts=polar, asymmetric_unit_pair_buried=au_buried,
                area_on_conv=sum(ds[p] for p in conv if p in set(iface)),
                dsasa={str(k): round(v, 2) for k, v in ds.items()},
                mind={str(k): round(v, 2) for k, v in mind.items()},
                dssp={str(k): str(v) for k, v in dssp.items()})


def project(raw, ref_gene, iface, table):
    al = Align.PairwiseAligner()
    al.substitution_matrix = substitution_matrices.load("BLOSUM62")
    al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, "global"
    ref = raw["uniprot"][ref_gene]["seq"]
    out = {}
    for g in GENES:
        s = raw["uniprot"][g]["seq"]
        if g == ref_gene:
            m, ident = {p: p for p in range(1, len(s) + 1)}, 100.0
        else:
            a = al.align(ref, s)[0]
            m, same = {}, 0
            for (i0, i1), (j0, j1) in zip(*a.aligned):
                for k in range(i1 - i0):
                    m[i0 + k + 1] = j0 + k + 1
                    same += ref[i0 + k] == s[j0 + k]
            ident = 100.0 * same / len(ref)
        eq = {m[p] for p in iface if p in m}
        hit = [p for p in table[g]["conv"] if p in eq]
        out[g] = dict(identity_to_ref=ident, iface_mapped=len(eq), conv_in_iface=hit,
                      residues=[f"{s[p-1]}{p}" for p in hit])
    return out


def main():
    raw = json.load(open(RAW))
    tmp = os.path.join(CACHE, "built")
    os.makedirs(tmp, exist_ok=True)
    table = series_table(raw)
    per = [analyse(pid, a, b, g, op, raw, tmp) for pid, a, b, g, op in SYSTEMS]
    ref = next(r for r in per if r["pdb"] == "2BQ0")
    proj = project(raw, "YWHAB", ref["iface"], table)
    checks = {}
    for r in per:
        own = sorted(r["conv_in_iface"])
        checks[r["pdb"]] = dict(
            gene=r["gene"], own=own, projected=sorted(proj[r["gene"]]["conv_in_iface"]),
            agree=own == sorted(proj[r["gene"]]["conv_in_iface"]))
    for pid, c in checks.items():
        if not c["agree"]:
            raise SystemExit(f"{pid}: interface set from coordinates {c['own']} does not match "
                             f"the projection {c['projected']}")
    res = dict(series_table={g: dict(v, conv_n=len(v["conv"])) for g, v in table.items()},
               structures=per, projection=proj, cross_structure_checks=checks)
    path = os.path.join(ROOT, "data", "iface_1433.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    print(f"wrote {path}")
    for r in per:
        print(f"{r['pdb']} {r['gene']} {r['resolution']} A, het "
              f"{r['composition']['het_components']}, "
              f"{len(r['composition']['chains_large'])} large and "
              f"{len(r['composition']['chains_small'])} small chains")
        print(f"{r['pdb']} {r['gene']}: {r['observed']} observed residues "
              f"{r['first']}..{r['last']}, {r['helix_pct']:.1f} per cent helix, "
              f"{r['beta_pct']:.1f} per cent beta, {r['buried_per_subunit']:.0f} A^2 buried, "
              f"{len(r['iface'])} interface residues, "
              f"{len(r['conv_in_iface'])} of {len(r['conv_dssp'])} convertible in the interface "
              f"{r['conv_in_iface']}")
    for g, p in proj.items():
        print(f"{g}: {len(p['conv_in_iface'])} of {len(table[g]['conv'])} convertible in the "
              f"interface {p['residues']}, alignment identity to P31946 {p['identity_to_ref']:.1f} per cent")


if __name__ == "__main__":
    main()
