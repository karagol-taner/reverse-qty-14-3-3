"""Fetch the seven human 14-3-3 paralogues from UniProt."""
import datetime, json, os, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), "data", "paralogue_uniprot.json")
ACC = {"YWHAQ": "P27348", "YWHAB": "P31946", "YWHAG": "P61981", "YWHAH": "Q04917",
       "YWHAZ": "P63104", "SFN": "P31947", "YWHAE": "P62258"}
FIELDS = "accession,sequence,ft_helix,ft_strand,ft_transmem,ft_chain,ft_signal"


def get(acc):
    url = ("https://rest.uniprot.org/uniprotkb/%s?format=json&fields=%s" % (acc, FIELDS))
    for _ in range(4):
        r = subprocess.run(["curl", "-sS", "--max-time", "60", url],
                           capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip().startswith("{"):
            return json.loads(r.stdout)
        time.sleep(3)
    raise SystemExit("fatal: could not fetch %s, nothing written" % acc)


def ranges(rec, kind):
    return sorted([f["location"]["start"]["value"], f["location"]["end"]["value"]]
                  for f in rec.get("features", []) if f["type"] == kind)


def main():
    out = {}
    for gene, acc in ACC.items():
        d = get(acc)
        if d.get("primaryAccession") != acc:
            raise SystemExit("fatal: %s returned %r" % (acc, d.get("primaryAccession")))
        seq = d["sequence"]["value"]
        if len(seq) != d["sequence"]["length"]:
            raise SystemExit("fatal: %s sequence is %d residues against the %d reported"
                             % (acc, len(seq), d["sequence"]["length"]))
        hel = ranges(d, "Helix")
        if not hel:
            raise SystemExit("fatal: %s carries no helix annotation" % acc)
        out[gene] = dict(acc=acc, seq=seq, helix=hel, strand=ranges(d, "Beta strand"),
                         transmem=ranges(d, "Transmembrane"),
                         mature=([1, len(seq)] if not ranges(d, "Chain")
                                 else ranges(d, "Chain")[0]),
                         signal=ranges(d, "Signal"))
        print("%-6s %s  %3d residues, %2d helices, %d strands"
              % (gene, acc, len(seq), len(hel), len(out[gene]["strand"])))
    if len(out) != len(ACC):
        raise SystemExit("fatal: fetched %d of %d paralogues" % (len(out), len(ACC)))
    json.dump(dict(uniprot=out,
                   provenance="rest.uniprot.org/uniprotkb, fields %s" % FIELDS,
                   retrieved=datetime.date.today().isoformat()),
              open(OUT, "w"), indent=1)
    print("wrote %s" % OUT)


if __name__ == "__main__":
    main()
