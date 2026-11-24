# The obligate 14-3-3 dimer interface holds reverse-QTY substitutions for membrane compatibility

**Taner Karagöl<sup>1,\*,‡</sup>, Alper Karagöl<sup>2,3,\*,‡</sup>**

<sup>1</sup>McGovern Institute for Brain Research, Massachusetts Institute of Technology, Cambridge, MA, USA  
<sup>2</sup>Broad Institute of MIT and Harvard, Cambridge, MA, USA  
<sup>3</sup>Department of Mechanical Engineering, Massachusetts Institute of Technology, Cambridge, MA, USA  

## Abstract

The QTY code and its reverse exchange polar and apolar residues on helices by a rule written on sequence and secondary structure alone. The rule therefore cannot see a contact that exists only in an assembled oligomer, and every target converted so far folds as a monomer. We applied the reverse code, glutamine to leucine, threonine to valine and tyrosine to phenylalanine at helical positions, to the seven human 14-3-3 paralogues, which function as obligate dimers. The rule reaches 22 to 30 positions per paralogue and crosses zero hydropathicity in 14-3-3 theta alone. Across four deposited structures and three paralogues, three converted positions fall inside the dimer interface in every one and account for 18.9 to 25.2 per cent of the surface buried on dimer formation. The deepest is a tyrosine whose hydroxyl donates two cross-interface hydrogen bonds. The converted dimer was run against its native control in water for 100 nanoseconds under thresholds declared beforehand. It keeps 1.00 of the native buried interfacial area and 1.01 of its inter-subunit contacts, against a hold threshold of 0.50. Deleting the hydroxyl produces a local retreat of 0.38 to 0.40 angstrom, and a contact the crystal places at 2.64 angstrom does not form in the native run. Reverse-QTY conversion therefore reaches an obligate dimer interface and leaves it intact over the window tested, with the burden quantified per structure before any construct is built.

**Keywords:** Reverse-QTY code; 14-3-3 proteins; protein design; dimer interface; molecular dynamics; membrane compatibility

## Repository layout

```
GitHub/
├── README.md                               This file
├── LICENSE                                 MIT licence for the code
├── Supplementary_Information.pdf           Supplementary tables and figure
├── code/                                   Analysis scripts and the simulation notebooks
│   ├── 00_fetch_sequences.py               UniProt entries and helix ranges for the seven paralogues
│   ├── pdb_survey_1433.py                  Deposited structures per paralogue, classified
│   ├── iface_1433.py                       Convertible set and dimer interface from the coordinates
│   ├── conservation_1433.py                Per-position conservation from the family alignment
│   ├── make_dcs_table.py                   Dynamic and rigid conserved scores per residue
│   ├── collect_values.py                   Every derived value gathered into one file
│   ├── collect_compute_environment.py      GROMACS build and accelerator read off the production notebooks
│   ├── build_figures.py                    Figures 1 to 5 and Supplementary Figure S1
│   ├── render_ga_structures.py             Cartoon renders of the dimer for the graphical abstract
│   ├── build_graphical_abstract.py         The graphical abstract
│   └── colab/                              GROMACS runs on Colab, in order
│       ├── 1433_rQTY_MD_1_prepare.ipynb    Fetch 5IQP and construct the native and converted dimers
│       ├── 1433_rQTY_MD_2_native.ipynb     Native dimer in water, GROMACS on a GPU
│       ├── 1433_rQTY_MD_2_converted.ipynb  Converted dimer in water, same protocol
│       ├── 1433_rQTY_MD_3_compare.ipynb    The one-sided criterion applied to the two runs
│       └── 1433_rQTY_MD_4_rmsf.ipynb       Per-residue fluctuation of both systems
├── data/                                   Derived data files
│   ├── paralogue_uniprot.json              Sequences and annotated helices, seven accessions
│   ├── pdb_survey_1433.json                Every PDB entry per paralogue with its classification
│   ├── iface_1433.json                     Convertible positions and the measured interface per structure
│   ├── interface_literature.json           Interface positions reported in the literature
│   ├── conservation_P27348.json            Alignment provenance and parameters
│   ├── conservation_P27348.tsv             Conservation per position of the lead paralogue
│   ├── prep_5IQP_build_spec.json           Substitutions and chains of the two simulated systems
│   ├── derived_values.json                 Every value the analysis produced, in one file
│   └── compute_environment.json            GROMACS version, build flags and accelerator of the two runs
├── figures/                                Figures of the paper, an svg beside each png
│   ├── fig_conservation.png                Figure 5, conservation of the converted positions
│   ├── fig_criterion.png                   Figure 3, the interface through both trajectories
│   ├── fig_design.png                      Figure 1, the conversion across the family
│   ├── fig_fluctuation.png                 Supplementary Figure S1, fluctuation per residue
│   ├── fig_fluctuation_map.png             Figure 4, the fluctuation on the fold
│   ├── fig_interface.png                   Figure 2, the dimer interface and the tyrosine
│   ├── ga_converted.png                    The same cartoon with the converted positions marked
│   ├── ga_native.png                       Cartoon of the dimer, for the graphical abstract
│   ├── graphical_abstract.png              The graphical abstract
│   └── graphical_abstract_toc.png          The graphical abstract at the table of contents size
└── md/                                     Simulation inputs and results
    ├── inputs/                             What the two systems were built from
    │   ├── build_spec.json                 Build specification written by the prepare notebook
    │   ├── substitutions_YWHAQ.csv         The 30 substitutions of the lead paralogue
    │   ├── CHARMM-GUI-checklist.txt        Settings used at CHARMM-GUI for both systems
    │   └── mdp/                            GROMACS parameters per system, minimisation, equilibration, production
    └── results/                            What the two runs produced
        ├── converted_report.json           Run record of the converted system
        ├── dcs_table.tsv                   Conservation joined to fluctuation, with both scores
        ├── kill_criterion.json             Verdict, ratios and the thresholds declared beforehand
        ├── native_report.json              Run record of the native system
        ├── rmsf.json                       Per-residue fluctuation, both systems
        ├── rmsf_converted_A.tsv            Fluctuation per residue, converted system, chain A
        ├── rmsf_converted_B.tsv            Fluctuation per residue, converted system, chain B
        ├── rmsf_native_A.tsv               Fluctuation per residue, native system, chain A
        ├── rmsf_native_B.tsv               Fluctuation per residue, native system, chain B
        ├── rmsf_table.tsv                  Fluctuation per residue, both systems side by side
        ├── native/                         Trajectory series of the native system
        └── converted/                      Trajectory series of the converted system
```

## Dataset summary

| File | Rows | Description |
| --- | --- | --- |
| `data/paralogue_uniprot.json` | 7 | Sequence, length and annotated helix ranges for each human 14-3-3 paralogue, from UniProt. |
| `data/pdb_survey_1433.json` | 645 | Every deposited entry per paralogue, marked for whether it is a dimer of the paralogue alone. |
| `data/iface_1433.json` | 4 | Convertible Q/T/Y set per paralogue and the measured dimer interface of each structure, over both interface definitions and their cutoff sweeps. |
| `data/conservation_P27348.tsv` | 245 | Conservation per position of 14-3-3 theta, with the coverage behind each column. |
| `md/results/rmsf_table.tsv` | 459 | C-alpha fluctuation per residue and chain for both systems, with the converted and interface flags. |
| `md/results/dcs_table.tsv` | 459 | Conservation joined to the fluctuation, with the dynamic and rigid conserved scores of each system. |
| `data/derived_values.json` | - | Every value the analysis produced, collected from the files that produced it. |

The 30 substitutions of the lead paralogue are listed in `md/inputs/substitutions_YWHAQ.csv` and recorded with the chains of both simulated systems in `data/prep_5IQP_build_spec.json`.

## Reproducing the analysis

The scripts take no arguments and resolve their own paths, so they run from anywhere in the tree.

```
python3 code/00_fetch_sequences.py
python3 code/pdb_survey_1433.py
python3 code/iface_1433.py
python3 code/conservation_1433.py
python3 code/make_dcs_table.py
python3 code/collect_values.py
python3 code/build_figures.py
```

Coordinates are not redistributed here. `code/colab/1433_rQTY_MD_1_prepare.ipynb` fetches the deposited assembly from RCSB, checks it, derives the converted positions and writes the two dimers into `md/inputs/`, which `code/iface_1433.py` and `code/build_figures.py` then read. The notebooks run in order: part 1 prepares the two systems, the two part 2 notebooks run them in water on a GPU, part 3 applies the criterion and part 4 measures the per-residue fluctuation. The GROMACS parameter files of both runs are under `md/inputs/mdp/`.

## Software dependencies

| Domain | Packages |
| --- | --- |
| Data wrangling | `numpy` |
| Plotting | `matplotlib` |
| Structural | `mdtraj`, `biopython` |
| Alignment | `pyfamsa`, or `famsa` on the path |
| External APIs | `requests`, `curl` (system) |
| Simulation | GROMACS 2024 with CUDA, CHARMM36m, run on Colab |

Tested on Python 3.10 and 3.11.

## Data availability

* Sequences and helix annotations are from UniProtKB (`https://rest.uniprot.org/`).
* Structures are from the RCSB Protein Data Bank (`https://www.rcsb.org/`) and are fetched and checked rather than redistributed here.
* The family alignment is built from the Pfam PF00244 members of UniProt, at the release recorded in `data/conservation_P27348.json`.
* Every derived file needed to redraw the figures is included under `data/`, `md/results/` and `figures/`.

## Citation

Karagöl T., Karagöl A. *The obligate 14-3-3 dimer interface holds reverse-QTY substitutions for membrane compatibility.* (2026).

## Contact

Taner Karagöl, taner.karagol@gmail.com, karagolt@mit.edu  
Alper Karagöl, alper.karagol@gmail.com, akaragol@mit.edu  

## License

The code in this repository is released under the MIT License. UniProt, Pfam and RCSB data are subject to the licenses of those providers.
