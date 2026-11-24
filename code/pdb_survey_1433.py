"""Experimental structures per 14-3-3 paralogue, classified by what else they carry."""
import json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "data", "pdb_survey_1433.json")
SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GQL = "https://data.rcsb.org/graphql"
ACC = {"YWHAQ": "P27348", "YWHAB": "P31946", "YWHAG": "P61981", "YWHAH": "Q04917",
       "YWHAZ": "P63104", "SFN": "P31947", "YWHAE": "P62258"}
FUSION_SLACK = 10

Q = """query($ids:[String!]!){entries(entry_ids:$ids){rcsb_id
 struct{title} exptl{method}
 rcsb_entry_info{resolution_combined polymer_entity_count nonpolymer_entity_count
                 deposited_polymer_monomer_count}
 rcsb_primary_citation{pdbx_database_id_DOI journal_abbrev year}
 polymer_entities{entity_poly{rcsb_sample_sequence_length}
   rcsb_polymer_entity_align{reference_database_accession
     aligned_regions{length}}}}}"""


def call(url, payload=None):
    cmd = ["curl", "-sS", "--max-time", "90"]
    if payload is not None:
        cmd += ["-H", "Content-Type: application/json", "-X", "POST",
                "-d", json.dumps(payload)]
    cmd.append(url)
    for _ in range(4):
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip().startswith("{"):
            return json.loads(r.stdout)
        time.sleep(3)
    raise SystemExit("fatal: request failed, nothing recorded: %s" % url)


def entries_for(acc):
    payload = {"query": {"type": "terminal", "service": "text", "parameters": {
        "attribute": "rcsb_polymer_entity_container_identifiers."
                     "reference_sequence_identifiers.database_accession",
        "operator": "exact_match", "value": acc}},
        "request_options": {"paginate": {"start": 0, "rows": 500},
                            "results_verbosity": "compact"},
        "return_type": "entry"}
    d = call(SEARCH, payload)
    ids = d.get("result_set", []) or []
    return ids, d.get("total_count", len(ids))


def aligned(ent, acc):
    best = 0
    for pe in ent.get("polymer_entities") or []:
        for al in pe.get("rcsb_polymer_entity_align") or []:
            if al.get("reference_database_accession") != acc:
                continue
            n = sum(r["length"] for r in al.get("aligned_regions") or [])
            best = max(best, n)
    return best


def seq_len(ent):
    return max([(pe.get("entity_poly") or {}).get("rcsb_sample_sequence_length") or 0
                for pe in ent.get("polymer_entities") or []] or [0])


def main():
    out = {}
    for gene, acc in ACC.items():
        ids, total = entries_for(acc)
        if not ids:
            raise SystemExit("fatal: no entries returned for %s (%s)" % (gene, acc))
        ents = []
        for i in range(0, len(ids), 100):
            ents += call(GQL, {"query": Q, "variables": {"ids": ids[i:i + 100]}})["data"]["entries"]
        rows = []
        for e in ents:
            info = e["rcsb_entry_info"]
            n_al, n_seq = aligned(e, acc), seq_len(e)
            rows.append(dict(
                id=e["rcsb_id"], method=(e["exptl"] or [{}])[0].get("method", ""),
                resolution=(info.get("resolution_combined") or [None])[0],
                polymer_entities=info.get("polymer_entity_count"),
                nonpolymer_entities=info.get("nonpolymer_entity_count"),
                monomers=info.get("deposited_polymer_monomer_count"),
                aligned=n_al, sequence_length=n_seq,
                fused=bool(n_al and n_seq - n_al > FUSION_SLACK),
                title=(e["struct"] or {}).get("title", ""),
                doi=(e.get("rcsb_primary_citation") or {}).get("pdbx_database_id_DOI"),
                journal=(e.get("rcsb_primary_citation") or {}).get("journal_abbrev"),
                year=(e.get("rcsb_primary_citation") or {}).get("year")))
        alone = [r for r in rows if r["polymer_entities"] == 1 and not r["fused"]]
        clean = [r for r in alone if r["nonpolymer_entities"] == 0]
        out[gene] = dict(acc=acc, n_entries=total, n_examined=len(rows),
                         truncated=total > len(rows), entries=rows,
                         paralogue_alone=[r["id"] for r in alone],
                         paralogue_alone_no_heteroatom=[r["id"] for r in clean],
                         with_paper=[r["id"] for r in alone if r["doi"]],
                         best_resolution=min([r["resolution"] for r in rows
                                              if r["resolution"]] or [None]))
        print("%-6s %s  %3d entries%s, %2d of the paralogue alone, %d of those "
              "with no heteroatom, %d with a paper, best resolution %s"
              % (gene, acc, total, " (%d examined)" % len(rows) if total > len(rows) else "",
                 len(alone), len(clean),
                 len([r for r in alone if r["doi"]]),
                 out[gene]["best_resolution"]))
        for r in sorted(alone, key=lambda r: (r["resolution"] is None, r["resolution"])):
            print("    %s %s %s het %s mono %s %s"
                  % (r["id"], r["method"][:18], r["resolution"],
                     r["nonpolymer_entities"], r["monomers"], (r["doi"] or "no paper")))
    json.dump(out, open(OUT, "w"), indent=1)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
