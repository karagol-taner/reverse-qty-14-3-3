"""Per-position conservation for 14-3-3 theta, for the dynamic conserved score."""

import json
import os
import sys
import urllib.parse
import urllib.request
import datetime

import numpy as np
from pyfamsa import Aligner, Sequence

REFERENCE = "P27348"
PFAM = "PF00244"
MIN_LENGTH_FRACTION = 0.6
MAX_LENGTH_FRACTION = 1.4
MIN_SEQUENCE_COVERAGE = 0.5
REWEIGHT_IDENTITY = 0.8
THREADS = max(1, os.cpu_count() or 1)

BACKGROUND = {
    "A": 0.0787, "C": 0.0157, "D": 0.0531, "E": 0.0636, "F": 0.0394,
    "G": 0.0691, "H": 0.0229, "I": 0.0592, "K": 0.0595, "L": 0.0958,
    "M": 0.0239, "N": 0.0427, "P": 0.0469, "Q": 0.0395, "R": 0.0530,
    "S": 0.0694, "T": 0.0585, "V": 0.0671, "W": 0.0137, "Y": 0.0323,
}
AA = "".join(sorted(BACKGROUND))
GAP = "-"
SYMBOLS = AA + GAP


def fetch(url):
    req = urllib.request.Request(url, headers={"User-Agent": "conservation_1433"})
    with urllib.request.urlopen(req, timeout=600) as r:
        return r.read().decode(), dict(r.headers)


def read_fasta(text):
    """-> list of (header, sequence), in file order."""
    out, name, buf = [], None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                out.append((name, "".join(buf)))
            name, buf = line[1:].strip(), []
        elif name is not None:
            buf.append(line.strip())
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def main(out_dir):
    stamp = datetime.datetime.now(datetime.timezone.utc)

    ref_text, _ = fetch("https://rest.uniprot.org/uniprotkb/%s.fasta" % REFERENCE)
    ref = read_fasta(ref_text)
    if len(ref) != 1:
        sys.exit("STOP: %d records returned for %s" % (len(ref), REFERENCE))
    ref_name, ref_seq = ref[0]
    print("reference %s, %d residues" % (REFERENCE, len(ref_seq)))

    query = "(xref:pfam-%s)" % PFAM
    url = ("https://rest.uniprot.org/uniprotkb/stream?query=%s&format=fasta"
           "&compressed=false" % urllib.parse.quote(query))
    text, headers = fetch(url)
    release = headers.get("X-UniProt-Release", "unrecorded")
    release_date = headers.get("X-UniProt-Release-Date", "unrecorded")
    records = read_fasta(text)
    print("%d sequences carrying %s, UniProt release %s"
          % (len(records), PFAM, release))

    lo = MIN_LENGTH_FRACTION * len(ref_seq)
    hi = MAX_LENGTH_FRACTION * len(ref_seq)
    kept, seen = [], set()
    for name, seq in records:
        s = seq.upper().replace("U", "C").replace("O", "K")
        if not (lo <= len(s) <= hi):
            continue
        if s in seen:
            continue
        seen.add(s)
        kept.append((name, s))
    if ref_seq not in seen:
        kept.insert(0, (ref_name, ref_seq))
        print("the reference was not in the family set and was added")
    print("%d after the length filter and exact duplicates" % len(kept))

    aligner = Aligner(guide_tree="upgma", threads=THREADS)
    msa = aligner.align([Sequence(n.split()[0].encode(), s.encode())
                         for n, s in kept])
    aligned = {r.id.decode(): r.sequence.decode() for r in msa}
    print("aligned, %d columns" % len(next(iter(aligned.values()))))

    ref_key = [k for k in aligned if REFERENCE in k]
    if len(ref_key) != 1:
        sys.exit("STOP: the reference matched %d rows of the alignment"
                 % len(ref_key))
    ref_row = aligned[ref_key[0]]
    cols = [j for j, c in enumerate(ref_row) if c != GAP]
    if "".join(ref_row[j] for j in cols) != ref_seq:
        sys.exit("STOP: the alignment's reference row does not read back as the "
                 "reference sequence, so the column mapping would be wrong")
    print("%d reference columns" % len(cols))

    code = {c: i for i, c in enumerate(SYMBOLS)}
    gap_i = code[GAP]
    rows = []
    for name, seq in aligned.items():
        r = [code.get(seq[j], gap_i) for j in cols]
        if sum(1 for v in r if v != gap_i) >= MIN_SEQUENCE_COVERAGE * len(cols):
            rows.append(r)
    X = np.array(rows, dtype=np.int8)
    print("%d sequences cover at least %.0f%% of the reference columns"
          % (len(X), 100 * MIN_SEQUENCE_COVERAGE))

    N, L = X.shape
    onehot = np.zeros((N, L * len(SYMBOLS)), dtype=np.float32)
    onehot[np.arange(N)[:, None].repeat(L, 1).ravel(),
           (np.arange(L) * len(SYMBOLS))[None, :].repeat(N, 0).ravel()
           + X.ravel()] = 1.0
    neighbours = np.zeros(N, dtype=np.int32)
    step = max(1, int(5e7 // max(N, 1)))
    for a in range(0, N, step):
        block = onehot[a:a + step] @ onehot.T
        neighbours[a:a + step] = (block >= REWEIGHT_IDENTITY * L).sum(1)
        print("  weights %d/%d" % (min(a + step, N), N), flush=True)
    w = 1.0 / neighbours
    n_eff = float(w.sum())
    print("effective depth %.1f sequences from %d" % (n_eff, N))

    freq = np.zeros((L, len(AA)), dtype=np.float64)
    for a in range(len(AA)):
        freq[:, a] = ((X == a) * w[:, None]).sum(0)
    occupancy = freq.sum(1)
    coverage = occupancy / n_eff
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(occupancy[:, None] > 0, freq / occupancy[:, None], 0.0)
    b = np.array([BACKGROUND[a] for a in AA])
    with np.errstate(divide="ignore", invalid="ignore"):
        terms = np.where(f > 0, f * np.log(f / b), 0.0)
    rel_entropy = terms.sum(1)
    ceiling = -np.log(b.min())
    conservation = np.clip(rel_entropy / ceiling, 0.0, 1.0)

    os.makedirs(out_dir, exist_ok=True)
    tsv = os.path.join(out_dir, "conservation_%s.tsv" % REFERENCE)
    with open(tsv, "w") as fh:
        fh.write("i\tA_i\tconservation\trelative_entropy_nats\tcoverage\n")
        for k in range(L):
            fh.write("%d\t%s\t%.6f\t%.6f\t%.4f\n"
                     % (k + 1, ref_seq[k], conservation[k], rel_entropy[k],
                        coverage[k]))

    meta = dict(
        written=stamp.isoformat(timespec="seconds"),
        reference=REFERENCE, reference_length=len(ref_seq),
        quantity=("relative entropy of each column's weighted amino-acid "
                  "frequencies against the background below, in nats, divided "
                  "by -log(min background) so that it lands in 0 to 1"),
        not_evcouplings=("computed here from the alignment described below. It "
                         "is the same quantity as the conservation column of an "
                         "EVcouplings frequencies file but not the same "
                         "alignment, so values are not interchangeable with an "
                         "EVcouplings run"),
        alignment=dict(
            source="UniProtKB", query=query, url=url,
            uniprot_release=release, uniprot_release_date=release_date,
            sequences_returned=len(records),
            after_length_filter_and_duplicates=len(kept),
            length_filter=[MIN_LENGTH_FRACTION, MAX_LENGTH_FRACTION],
            aligner="FAMSA via pyfamsa, UPGMA guide tree",
            columns_in_alignment=len(ref_row),
            reference_columns=len(cols),
            minimum_sequence_coverage=MIN_SEQUENCE_COVERAGE,
            sequences_scored=int(N)),
        reweighting=dict(identity=REWEIGHT_IDENTITY,
                         rule="one over the number of sequences at or above "
                              "that identity, gap against gap counted as a "
                              "match",
                         effective_depth=n_eff,
                         effective_depth_per_column=n_eff / len(cols)),
        background=BACKGROUND,
        background_source="Robinson and Robinson 1991",
        note=("relative entropy is measured against the background, so a column "
              "held at a rare residue scores above a column held just as "
              "firmly at a common one. The dynamics-aware pipeline's "
              "conservation column behaves the same way"),
        coverage=dict(mean=float(coverage.mean()), min=float(coverage.min()),
                      min_at_position=int(coverage.argmin()) + 1))
    with open(os.path.join(out_dir, "conservation_%s.json" % REFERENCE),
              "w") as fh:
        json.dump(meta, fh, indent=1)
    print("wrote", tsv)
    print("effective depth per column %.1f" % (n_eff / len(cols)))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data")
