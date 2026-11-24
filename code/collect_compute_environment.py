"""Record the execution environment of the two production runs."""

import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PROJECT = os.path.dirname(HERE)
NOTEBOOKS = ["1433_rQTY_MD_2_native.ipynb", "1433_rQTY_MD_2_converted.ipynb"]
OUT = os.path.join(PROJECT, "data", "compute_environment.json")

GPU = "NVIDIA A100"
CUDA_MAJOR = "12"


def source(name):
    path = os.path.join(HERE, "colab", name)
    nb = json.load(open(path, encoding="utf-8"))
    return "\n".join("".join(c.get("source", [])) for c in nb["cells"])


def one(pattern, text, label):
    hits = sorted(set(re.findall(pattern, text)))
    if len(hits) != 1:
        raise SystemExit("STOP: %s matched %d values in the notebook source, "
                         "expected one: %s" % (label, len(hits), hits))
    return hits[0]


def flags(text):
    """Return the cmake build flags of the GROMACS build."""
    m = re.search(r"cmake \.\.(.{0,400}?)%\s*\(", text, re.S)
    if not m:
        raise SystemExit("STOP: no cmake invocation found in the notebook")
    line = re.sub(r"[\"\n]+", " ", m.group(1))
    out = {}
    for f in re.finditer(r"-D([A-Z_]+)=(%s|[A-Za-z0-9_.]+)", line):
        out[f.group(1)] = f.group(2)
    for want in ("GMX_GPU", "GMX_DOUBLE", "GMX_BUILD_OWN_FFTW"):
        if want not in out:
            raise SystemExit("STOP: %s not among the cmake flags parsed from "
                             "%r" % (want, line.strip()))
    return out


def main():
    texts = {n: source(n) for n in NOTEBOOKS}

    env = {}
    for name, text in texts.items():
        f = flags(text)
        simd = one(r'GMX_SIMD\s*=\s*"([A-Za-z0-9_.]+)"', text, "GMX_SIMD")
        env[name] = {
            "gromacs_version": one(r'GMX_VER\s*=\s*"([0-9.]+)"', text,
                                   "GMX_VER"),
            "gmx_simd": simd,
            "gpu_backend": f.get("GMX_GPU"),
            "double_precision": f.get("GMX_DOUBLE"),
            "own_fftw": f.get("GMX_BUILD_OWN_FFTW"),
            "ns_per_build": int(one(r"TOTAL_NS\s*=\s*(\d+)", text,
                                    "TOTAL_NS")),
            "tarball_host": one(r"wget -q https://([a-z.]+)/gromacs/", text,
                                "tarball host"),
        }

    first = env[NOTEBOOKS[0]]
    for name in NOTEBOOKS[1:]:
        if env[name] != first:
            raise SystemExit("STOP: the two production notebooks declare "
                             "different environments: %s against %s"
                             % (first, env[name]))

    if first["gpu_backend"] != "CUDA":
        raise SystemExit("STOP: the build reports %r rather than CUDA, so "
                         "no accelerator can be recorded"
                         % first["gpu_backend"])
    if first["double_precision"] != "OFF":
        raise SystemExit("STOP: GMX_DOUBLE is %r, so the precision recorded "
                         "here would be wrong" % first["double_precision"])

    out = {
        "gromacs_version": first["gromacs_version"],
        "gromacs_source": first["tarball_host"],
        "built_from_source": True,
        "gpu_backend": first["gpu_backend"],
        "gpu": GPU,
        "cuda_major": CUDA_MAJOR,
        "simd": first["gmx_simd"],
        "precision": "single",
        "own_fftw": first["own_fftw"] == "ON",
        "ns_per_build": first["ns_per_build"],
        "notebooks": sorted(NOTEBOOKS),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("wrote %s" % os.path.relpath(OUT, PROJECT))
    for k in sorted(out):
        print("  %-20s %s" % (k, out[k]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
