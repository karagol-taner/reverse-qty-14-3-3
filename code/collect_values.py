"""Collect every derived value of the project into data/derived_values.json."""

import json
import os
import random
import statistics as st
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
DATA = os.path.join(PROJECT, "data")
RESULTS = os.path.join(PROJECT, "md", "results")
PLACEHOLDER = "PLACEHOLDER, not yet measured"


def fatal(msg):
    sys.exit("STOP: " + msg)


def load(path, required=True):
    if not os.path.isfile(path):
        if required:
            fatal("no %s" % path)
        return None
    with open(path) as fh:
        return json.load(fh)


def rng(values, fmt="%.1f"):
    """-> 'low to high' over the values, for a range stated in prose."""
    return "%s to %s" % (fmt % min(values), fmt % max(values))


def main():
    iface = load(os.path.join(DATA, "iface_1433.json"))
    uni = load(os.path.join(DATA, "paralogue_uniprot.json"))
    lit = load(os.path.join(DATA, "interface_literature.json"))
    survey = load(os.path.join(DATA, "pdb_survey_1433.json"))
    spec = load(os.path.join(DATA, "prep_5IQP_build_spec.json"))
    crit = load(os.path.join(RESULTS, "kill_criterion.json"))
    rep = {b: load(os.path.join(RESULTS, "%s_report.json" % b))
           for b in ("native", "converted")}
    cons_meta = load(os.path.join(DATA, "conservation_P27348.json"))
    rmsf = load(os.path.join(RESULTS, "rmsf.json"), required=False)

    V = {}

    S = iface["series_table"]
    order = sorted(S, key=lambda g: -S[g]["conv_n"])
    V["design"] = dict(
        paralogues=len(S),
        rows=[dict(gene=g, acc=S[g]["acc"], length=S[g]["length"],
                   helix_pct=S[g]["helix_pct"], qty=S[g]["qty"],
                   conv_n=S[g]["conv_n"], conv_pct=S[g]["conv_pct"],
                   gravy_native=S[g]["gravy_native"],
                   gravy_conv=S[g]["gravy_conv"]) for g in order],
        length_range=rng([S[g]["length"] for g in S], "%d"),
        helix_range=rng([S[g]["helix_pct"] for g in S]),
        conv_n_range=rng([S[g]["conv_n"] for g in S], "%d"),
        conv_pct_range=rng([S[g]["conv_pct"] for g in S]),
        beta_note="essentially no beta content in any deposited structure",
        crosses_zero=[g for g in S if S[g]["gravy_conv"] > 0],
        lead="YWHAQ", lead_acc=S["YWHAQ"]["acc"],
        lead_conv_n=S["YWHAQ"]["conv_n"],
        lead_conv_pct=S["YWHAQ"]["conv_pct"],
        lead_qty=S["YWHAQ"]["qty"],
        lead_length=S["YWHAQ"]["length"],
        lead_gravy_native=S["YWHAQ"]["gravy_native"],
        lead_gravy_conv=S["YWHAQ"]["gravy_conv"],
        lead_positions=S["YWHAQ"]["conv"],
        lead_unconverted=S["YWHAQ"]["qty"] - S["YWHAQ"]["conv_n"])

    seq = uni["uniprot"]["YWHAQ"]["seq"]
    if len(seq) != S["YWHAQ"]["length"]:
        fatal("the YWHAQ sequence is %d residues and the series table says %d"
              % (len(seq), S["YWHAQ"]["length"]))
    by_class = {}
    for p in S["YWHAQ"]["conv"]:
        by_class[seq[p - 1]] = by_class.get(seq[p - 1], 0) + 1
    if sum(by_class.values()) != S["YWHAQ"]["conv_n"]:
        fatal("the converted positions do not all carry Q, T or Y")
    three = {"Q": "Gln", "T": "Thr", "Y": "Tyr"}
    lead5 = [r for r in iface["structures"] if r["pdb"] == "5IQP"]
    if len(lead5) != 1:
        fatal("%d records for 5IQP in iface_1433.json" % len(lead5))
    V["design"]["lead_iface_names"] = [
        three[seq[p - 1]] + str(p)
        for p in sorted(lead5[0]["conv_in_iface"])]
    V["design"]["lead_by_class"] = by_class
    V["design"]["lead_class_sentence"] = ", ".join(
        "%d %s to %s" % (by_class[a], {"Q": "glutamine", "T": "threonine",
                                       "Y": "tyrosine"}[a],
                         {"Q": "leucine", "T": "valine",
                          "Y": "phenylalanine"}[a])
        for a in "QTY" if a in by_class)

    rows = []
    for s in iface["structures"]:
        share = 100.0 * s["area_on_conv"] / s["buried_per_subunit"]
        rows.append(dict(pdb=s["pdb"], gene=s["gene"],
                         resolution=s["resolution"], first=s["first"],
                         last=s["last"], observed=s["observed"],
                         helix_pct=s["helix_pct"], beta_pct=s["beta_pct"],
                         buried=s["buried_per_subunit"],
                         n_iface=len(s["iface"]),
                         conv_in_iface=s["conv_in_iface"],
                         area_on_conv=s["area_on_conv"], share_pct=share,
                         polar=s["polar_contacts"]))
    tyr = []
    for s in iface["structures"]:
        y = max(int(k) for k in s["polar_contacts"])
        for c in s["polar_contacts"][str(y)]:
            tyr.append(c["distance"])
    V["interface"] = dict(
        rows=rows, structures=len(rows),
        genes=sorted({r["gene"] for r in rows}),
        n_iface_range=rng([r["n_iface"] for r in rows], "%d"),
        buried_range=rng([r["buried"] for r in rows], "%.0f"),
        conv_in_iface_counts=sorted({len(r["conv_in_iface"]) for r in rows}),
        share_range=rng([r["share_pct"] for r in rows]),
        tyr_bond_range=rng(tyr, "%.2f"), tyr_bonds=len(tyr),
        projection=iface["projection"],
        cross_checks=iface["cross_structure_checks"],
        all_agree=all(v["agree"] for v in
                      iface["cross_structure_checks"].values()),
        lead_iface=[r for r in rows if r["pdb"] == "5IQP"][0])

    byg = iface["projection"]
    V["interface"]["per_paralogue"] = [
        dict(gene=g, conv_n=S[g]["conv_n"], in_iface=len(byg[g]["conv_in_iface"]),
             residues=byg[g]["residues"], identity=byg[g]["identity_to_ref"])
        for g in order]
    V["interface"]["in_iface_counts"] = sorted(
        {len(byg[g]["conv_in_iface"]) for g in byg})

    V["build"] = dict(spec={k: spec[k] for k in spec
                            if not isinstance(spec[k], (list, dict))},
                      raw=spec)

    ifc = crit["interface"]
    V["md"] = dict(
        verdict=crit["verdict"], ns=crit["ns"], frames=crit["frames"],
        thresholds=crit["thresholds"],
        burn_ns=crit["thresholds"]["burn_fraction"] * crit["ns"]["native"],
        buried=ifc["buried"], contacts=ifc["contacts"],
        mindist_tail=crit["converted_mindist_tail_A"],
        per_position=crit["per_position_mindist"],
        rmsd=crit["rmsd_A"], polar=crit["native_polar_pairs"],
        atoms={b: rep[b]["atoms"] for b in rep},
        written=crit["written"])
    V["md"]["position_82_change"] = [
        crit["per_position_mindist"]["82_sub%d" % k]["change_A"]
        for k in (1, 2)]
    V["md"]["position_8_15_change"] = [
        crit["per_position_mindist"]["%d_sub%d" % (p, k)]["change_A"]
        for p in (8, 15) for k in (1, 2)]
    V["md"]["asp21_occupancy"] = [crit["native_polar_pairs"]
                                  ["%s82OH to ASP21 OD2" % c]["occupancy_3_5"]
                                  for c in "AB"]
    V["md"]["arg18_occupancy"] = [crit["native_polar_pairs"]
                                  ["%s82OH to ARG18 NE" % c]["occupancy_3_5"]
                                  for c in "AB"]
    V["md"]["lys75_occupancy"] = [crit["native_polar_pairs"]
                                  ["%s8NE2 to LYS75 O" % c]["occupancy_3_5"]
                                  for c in "AB"]
    V["md"]["lys75_distance"] = [crit["native_polar_pairs"]
                                 ["%s8NE2 to LYS75 O" % c]["mean"]
                                 for c in "AB"]
    V["md"]["lys75_crystal_A"] = [c["distance"] for c in
                                  V["interface"]["lead_iface"]["polar"]["8"]]

    path = os.path.join(DATA, "conservation_P27348.tsv")
    if not os.path.isfile(path):
        fatal("no %s" % path)
    cons = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(head, line.rstrip("\n").split("\t")))
            cons[int(r["i"])] = (r["A_i"], float(r["conservation"]),
                                 float(r["coverage"]))
    modelled = [i for i in cons if i <= V["interface"]["lead_iface"]["last"]]
    conv = set(S["YWHAQ"]["conv"])
    a = [cons[i][1] for i in modelled if i in conv]
    b = [cons[i][1] for i in modelled if i not in conv]
    obs = st.mean(a) - st.mean(b)
    pool, hits, N = a + b, 0, 20000
    rnd = random.Random(0)
    for _ in range(N):
        rnd.shuffle(pool)
        if abs(st.mean(pool[:len(a)]) - st.mean(pool[len(a):])) >= abs(obs):
            hits += 1
    rank = sorted(modelled, key=lambda k: -cons[k][1])
    groove = [49, 56, 127, 128]
    V["conservation"] = dict(
        depth_per_position=cons_meta["reweighting"]["effective_depth_per_column"],
        effective_depth=cons_meta["reweighting"]["effective_depth"],
        sequences_scored=cons_meta["alignment"]["sequences_scored"],
        sequences_returned=cons_meta["alignment"]["sequences_returned"],
        release=cons_meta["alignment"]["uniprot_release"],
        pfam=cons_meta["alignment"]["query"],
        reweight=cons_meta["reweighting"]["identity"],
        mean_all=st.mean(cons[i][1] for i in cons),
        modelled=len(modelled),
        mean_converted=st.mean(a), mean_rest=st.mean(b),
        n_converted=len(a), n_rest=len(b),
        difference=obs, permutations=N, p=(hits + 1) / (N + 1),
        groove={str(i): dict(residue=cons[i][0], conservation=cons[i][1],
                             rank=rank.index(i) + 1) for i in groove},
        interface_conversions={
            str(i): dict(residue=cons[i][0], conservation=cons[i][1],
                         rank=rank.index(i) + 1)
            for i in V["interface"]["lead_iface"]["conv_in_iface"]},
        ranked_of=len(modelled))

    path = os.path.join(RESULTS, "dcs_table.tsv")
    if not os.path.isfile(path):
        V["dcs"] = dict(available=False, note=PLACEHOLDER)
    else:
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            tab = [dict(zip(head, ln.rstrip("\n").split("\t"))) for ln in fh]
        if not tab:
            fatal("%s carries a header and no rows" % path)
        bypos = {}
        for r in tab:
            bypos.setdefault(int(r["residue"]), []).append(r)
        keys = ["conservation", "dcs_native_0_10", "dcs_converted_0_10",
                "rcs_native_0_10", "rcs_converted_0_10",
                "norm_dynamics_native", "norm_dynamics_converted"]
        for k in keys + ["converted_position"]:
            if k not in head:
                fatal("%s carries no %s column" % (path, k))
        rec = {}
        for p, rs in bypos.items():
            rec[p] = {k: st.mean(float(x[k]) for x in rs) for k in keys}
            rec[p]["converted"] = rs[0]["converted_position"] == "1"
            rec[p]["interface"] = rs[0]["interface"] == "1"
            rec[p]["resname"] = (rs[0]["resname_native"],
                                 rs[0]["resname_converted"])
            rec[p]["chains"] = len(rs)
        inconv = [p for p in rec if rec[p]["converted"]]
        other = [p for p in rec if not rec[p]["converted"]]
        if len(inconv) != len(S["YWHAQ"]["conv"]):
            fatal("%s marks %d converted positions, the rule converts %d"
                  % (path, len(inconv), len(S["YWHAQ"]["conv"])))
        scores = {}
        for k in keys:
            x = [rec[p][k] for p in inconv]
            y = [rec[p][k] for p in other]
            d = st.mean(x) - st.mean(y)
            pool, hits, N = x + y, 0, 20000
            rnd = random.Random(0)
            for _ in range(N):
                rnd.shuffle(pool)
                if abs(st.mean(pool[:len(x)]) - st.mean(pool[len(x):])) >= abs(d):
                    hits += 1
            scores[k] = dict(converted=st.mean(x), rest=st.mean(y),
                             difference=d, p=(hits + 1) / (N + 1))
        iface = V["interface"]["lead_iface"]["conv_in_iface"]
        V["dcs"] = dict(
            available=True, positions=len(rec), n_converted=len(inconv),
            n_rest=len(other), permutations=20000, scores=scores,
            interface={str(i): {k: rec[i][k] for k in keys}
                       for i in iface if i in rec},
            converted_positions=[
                dict(position=p, resname_native=rec[p]["resname"][0],
                     resname_converted=rec[p]["resname"][1],
                     interface=rec[p]["interface"], chains=rec[p]["chains"],
                     **{k: rec[p][k] for k in keys})
                for p in sorted(inconv)],
            table_rows=len(tab))

    if rmsf is None:
        V["rmsf"] = dict(available=False, note=PLACEHOLDER)
    else:
        bc = rmsf["by_class"]
        V["rmsf"] = dict(
            available=True, selection=rmsf["selection"],
            residues=rmsf["residues"], mean=rmsf["mean_A"],
            by_class={c: {b: bc[b][c] for b in bc} for c in bc["native"]},
            largest=rmsf["largest_increases"],
            measured_with=rmsf["measured_with"],
            units=rmsf.get("units", ""))

    V["generated_by"] = "code/collect_values.py"
    out = os.path.join(DATA, "derived_values.json")
    with open(out, "w") as fh:
        json.dump(V, fh, indent=1)
    print("wrote", os.path.relpath(out, PROJECT))
    print("  %d paralogues, %d structures, verdict %s, conservation depth %.1f"
          % (V["design"]["paralogues"], V["interface"]["structures"],
             V["md"]["verdict"], V["conservation"]["depth_per_position"]))
    print("  per-residue fluctuation: %s"
          % ("present" if V["rmsf"]["available"] else PLACEHOLDER))


if __name__ == "__main__":
    main()
