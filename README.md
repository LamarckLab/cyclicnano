# cyclicnano

De novo design of cyclic protein oligomers for symmetric nanoparticle assembly:
RFdiffusion backbones, objective geometric measurement, ProteinMPNN sequences and
AlphaFold2 self-consistency validation, in one unattended command.

Built to produce the C5 and C3 components of I53 icosahedral nanoparticles, but
nothing in the pipeline is specific to that target: the symmetry order is a
parameter, so the same run designs any cyclic homo-oligomer from C2 to C12.

The name is the scope. *Cyclic*, because a Cn point group is defined by a single
n-fold rotation axis, which is what distinguishes it from the dihedral and
polyhedral groups. *Nano*, for the particles those components assemble into.

[![CI](https://github.com/LamarckLab/cyclicnano/actions/workflows/ci.yml/badge.svg)](https://github.com/LamarckLab/cyclicnano/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green)](LICENSE)

---

## What this is

An I53 icosahedral nanoparticle is built from two cyclic building blocks: a C5 pentamer
on the five-fold axes and a C3 trimer on the three-fold axes. Designing one therefore
starts by designing and validating those two oligomers.

That front half of the protocol is a fixed sequence of four tool invocations with a
human making judgement calls in between. This package removes the human from the loop
by turning every one of those judgement calls into a threshold in a config file.

Nothing in stages 00-04 is specific to I53, or even to C5 and C3. The symmetry order is
a parameter, so the same pipeline designs any cyclic homo-oligomer from C2 to C12. I53
is the application; cyclic oligomer design is the tool.

## Pipeline

```
00  generate    RFdiffusion, symmetric mode          ->  N backbones
01  measure     29 geometric metrics per backbone    ->  metrics table (no filtering by default)
02  design      ProteinMPNN, tied positions          ->  M sequences per backbone
03  fold        AlphaFold2 monomer, single-sequence  ->  predicted structures
04  validate    CA RMSD vs design, success rate      ->  ranked backbones
```

Stages 05-08 (fusion search, icosahedral assembly, junction redesign, reporting) are
not implemented yet. See [Status](#status).

## Design principles

Three rules make the pipeline reproducible and the thresholds tunable after the fact.

**Compute everything.** Every metric is computed for every candidate regardless of
whether any rule references it. Metrics that only the future fusion stage will need,
such as terminal helix lengths, are already recorded.

**Filter declaratively, and only when asked.** No rules ship enabled: every backbone
reaches sequence design, and stage 01 records rather than rejects. Thresholds, when
wanted, live in YAML as expressions over metric names, never in code:

```yaml
backbone_filter:
  rules:
    - "n_clash == 0"
    - "helix_frac >= 0.35"
    - "rg_ratio <= 1.30"
```

Because the metrics are recorded either way, a filter can be applied after a run
finishes by adding rules and replaying stage 01 alone.

**Never delete.** A rejected candidate keeps all of its metrics and gains a pass flag
plus the rules it failed. Retuning a threshold replays one CPU-bound stage instead of
a GPU-day, and the rejected population stays available for analysis.

Together these mean the thresholds are not buried assumptions. `filter_summary.json`
reports how many candidates each individual rule rejected, so it is always visible
which threshold is doing the work.

## Installation

```bash
git clone https://github.com/LamarckLab/cyclicnano.git
cd cyclicnano
python -m pip install -e ".[parquet]"
```

The package itself depends only on numpy, pandas and PyYAML. The external tools are
called as subprocesses and are not Python dependencies:

| Stage | Tool | Notes |
|-------|------|-------|
| 00 | [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) | symmetric oligomer mode |
| 02 | [ProteinMPNN](https://github.com/dauparas/ProteinMPNN) | tied positions required |
| 03 | [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) | AlphaFold2 weights |

Point the config at them via `paths.*` and name their environments under
`backend.envs`. AlphaFold2 parameters are downloaded by localcolabfold and are subject
to DeepMind's own licence terms; check those before any non-academic use.

## Usage

A run is described by two config files: `base.yaml` holds everything the pentamer and
the trimer must share, and a second file supplies the symmetry and the run name. The
split is deliberate, since settings copied between two self-contained files drift
apart and make the two sets of results incomparable.

```bash
# validate the configuration without running anything
cyclicnano check --config configs/runs/base.yaml --config configs/runs/c5.yaml                  --profile configs/profiles/amax.yaml

# print the constructed tool commands without executing them
cyclicnano run --config configs/runs/base.yaml --config configs/runs/c5.yaml                --profile configs/profiles/amax.yaml --dry-run

# full run; swap c5.yaml for c3.yaml to design the trimer instead
cyclicnano run --config configs/runs/base.yaml --config configs/runs/c5.yaml                --profile configs/profiles/amax.yaml

# re-run one stage after changing a threshold, repeating no GPU work
cyclicnano run --config ... --stages 01 --force

# inspect one backbone without the pipeline
cyclicnano metrics path/to/backbone.pdb --symmetry C5
```

Every run writes `logs/resolved_config.yaml`, the fully merged configuration including
any `--set` overrides, so a result directory records the settings that produced it.

## Output

```
runs/<name>/
├── tables/
│   ├── backbones.csv        one row per backbone, all metrics, pass flags, rankings
│   └── sequences.csv        one row per designed sequence, RMSD and pLDDT
├── 00_backbones/            RFdiffusion output, one PDB per backbone
├── 01_backbone_filter/
│   ├── filter_summary.json  the rules, and the rejection count of each
│   └── passed/              the backbones that survived; stage 02 reads this
├── 02_sequences/            one ProteinMPNN FASTA per backbone
├── 03_predictions/          two files per sequence: AlphaFold2 PDB and scores
├── 04_delivery/
│   ├── delivery.csv         one row per backbone: success rate, best RMSD, verdict
│   └── <backbone_id>/       backbone, every sequence, every prediction
└── logs/commands.jsonl      every subprocess, its exit code and its runtime
```

Self-consistency scoring has no directory. It writes columns on the backbone table,
and the ranked view of those columns is `delivery.csv`; a second ranking file would
be the same numbers under another name, free to disagree with the first.

Delivery covers every backbone that passed the **shape** filter, flagged with its
stage-04 verdict rather than filtered by it. The RMSD cutoffs are uncalibrated and a
near miss is worth looking at, so the information is kept and the judgement left to
the reader.

`backbones.csv` is the artifact to keep. It holds every metric for every candidate,
including the rejected ones, which is what makes a filter ablation possible after the
run rather than requiring a second one.

## Two traps this package handles for you

**RFdiffusion contig length is the total, not per chain.** A C5 of 60-residue subunits
needs `contigmap.contigs=[300-300]`. Passing 60 silently produces a much smaller
oligomer. The conversion lives in `Config.total_length`.

**ProteinMPNN must tie symmetric positions.** Without `--homooligomer 1` each chain of
a Cn backbone is designed independently and the sequences cannot form a homo-oligomer.
`design.tie_chains: false` is rejected at startup.

A third is a matter of interpretation rather than a bug: AlphaFold2 is run in
single-sequence mode. A de novo design has no homologues, so an MSA search returns
either nothing or unrelated hits, and the resulting confidence is not meaningful.
`fold.msa_mode` is flagged at startup if it is set to anything else.

## Testing

```bash
pytest                     # 169 checks
python tests/run_all.py    # same suite without pytest
```

The suite runs against synthetic Cn oligomers built from ideal helices, so the
symmetry order, chain RMSD, helix fraction and clash count are all known in advance.
It needs no GPU, no model weights and no external tools, and it covers stages 01 and
04 end to end by synthesising the outputs of the GPU stages. Three real RFdiffusion
C5 backbones are shipped as fixtures, one per stage-01 verdict: an alpha/beta one,
a two-helix hairpin and a three-helix bundle.

## What the numbers rest on

The two filters are not equally well founded, and it is worth knowing which is which.

**Stage 01, shape.** Not calibrated here, and not meant to be. The shipped preset,
`configs/presets/pure_three_helix.yaml`, requires no beta at all and exactly three
helices, on two experimental grounds: beta-rich designs express poorly, and the beta
that appears in the generated set is mostly sheet meeting sheet edge-on across the
subunit interface, which does not hold an oligomer together. Neither objection is
visible to any metric computed downstream. An earlier preset derived four thresholds
from 28 hand-labelled backbones; the rule here overrides those labels deliberately,
and the file says so and says what it cost.

The preset is strict and the yield reflects it: about 2.5 percent of C5 backbones
survive it, measured twice at 200 backbones each. That is why `n_backbones` is 500.
C3 runs richer, around 9 percent.

Relaxing the strand rule was tried and rejected. Admitting up to 0.10 strand doubles
the count, and the backbones it admits carry noticeably more loop; the rule stays at
zero on that inspection.

**Stage 04, self-consistency.** `rmsd_cutoff: 1.0`, `plddt_cutoff: 80`,
`success_rate >= 0.3` is **not calibrated against anything**.
They are reasonable and deliberately strict, since the literature usually calls a
design self-consistent at 2.0 angstrom, but no measurement or experimental outcome
sets them. Treat them as a starting point.

**Only the monomer is validated.** AlphaFold2 folds a single chain, so a passing
design is one whose sequence encodes the designed subunit fold. Nothing here tests
that copies of it assemble into the ring; that question is left to experiment.

## Status

Implemented and tested: stages 00-04, the filter engine, the run manifest, and the
conda backend.

Not yet implemented:

- **Stages 05-08** - fusion search, icosahedral assembly with clash and closure
  checks, junction redesign, HTML reporting.
- **Container images** - the backend abstraction already covers Docker and Apptainer;
  the images themselves are not built yet.
- **Threshold calibration** - the shipped thresholds are conservative starting points
  derived from general design practice, not values validated against experimental
  outcomes. Treat them as a starting point and retune against your own run.

## Citation

See [CITATION.cff](CITATION.cff).

## License

MIT. See [LICENSE](LICENSE).
