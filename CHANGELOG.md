# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added (generation)
- `generate.noise_scale`, setting RFdiffusion's `denoiser.noise_scale_ca` and
  `noise_scale_frame` together as the tool's own documentation advises. Default 0.5,
  chosen on measurement: 40 C5 backbones per condition with identical seeds gave
  0, 5 and 4 backbones of three or more helices at noise 1.0, 0.5 and 0.0, and 5, 11
  and 26 with no beta at all.
- 0.0 was rejected despite producing the most helix. It takes the two-helix count to
  34 of 40 and the median longest helix from 19 to 28 residues, which is the
  long-hairpin fold the protocol rules out: more helix is not the same as more of the
  wanted topology. The diversity cost at 0.5 is small, median pairwise CA RMSD
  unchanged at 10.6 A and near-duplicate pairs 1 -> 3 out of 780.

### Fixed (dry run)
- `--dry-run` now prints each constructed command, which is what the README has
  always claimed it does. It only wrote them to `logs/commands.jsonl`.
- `--dry-run` no longer fails on an empty run directory. Stages 02 and 03 raised
  because nothing upstream existed, so the flag broke in the one situation it is for,
  before the first real run, and the ProteinMPNN and colabfold commands could never
  be inspected at all. Both build their commands from directory paths alone, so they
  are now shown without the data.
- The CPU stages honour `--dry-run`. Taking no Runner, they previously did their real
  work: a dry run rewrote the folder stage 02 reads and laid out a delivery tree from
  whatever the tables happened to hold.

### Changed (stage-01 rule set)
- `target.n_backbones` 200 -> 500. A volume decision, not a clever one: the stage-01
  pass rate is about 2.5 percent for C5 and nothing tried has moved it.
- `configs/presets/low_strand_compact.yaml` removed. It was kept to reproduce the C3
  results generated with it; those have been discarded, so nothing depends on it.
- Relaxing the strand rule to 0.10 was measured and rejected. It takes 5 backbones in
  200 to 11, and four of the six admitted have strands too far from the neighbouring
  subunit to form an edge-on sheet, but they carry visibly more loop on inspection.
- The shipped preset is now `configs/presets/pure_three_helix.yaml`: no beta at all
  and exactly three helices. It replaces `low_strand_compact.yaml`, which capped the
  strand fraction at 0.15 and the distance from the oligomer centroid to its furthest
  CA at 28 A.
- Both earlier rules went for experimental reasons. The beta that appears in the
  generated set is mostly sheet meeting sheet edge-on across the subunit interface,
  which does not hold an oligomer together, and beta-rich designs express poorly; a
  fraction of 0.15 still admits that. The two-helix hairpin that dominates the C5
  output has helices too long to design against.
- The radius ceiling was dropped rather than retuned. At the thresholds under
  consideration (24 to 27 A) it would have rejected the only qualifying C5 backbone,
  whose radius is 27.9 A, taking the yield to zero. The metric is still computed and
  written to the backbone table.
- `n_backbones` 100 -> 200 and `n_seq_per_backbone` 10 -> 20. Measured on 100
  backbones of each symmetry, the new rule passes 1 C5 and 9 C3, so fewer backbones
  survive and each is sampled harder.

### Changed (output layout)
- `05_delivery/` is now `04_delivery/`, and holds one folder per backbone that passed
  **stage 01** rather than one PDB per backbone that passed stage 04. Each folder
  carries the backbone, every sequence written for it, and the prediction and score
  file for each of those sequences; the stage-04 verdict is a column in
  `delivery.csv` instead of a filter. The cutoffs are uncalibrated and a near miss is
  worth looking at.
- Self-consistency scoring no longer writes a directory. `04_validate/
  ranked_backbones.csv` is gone; its content is `04_delivery/delivery.csv`, so the
  two cannot disagree. The CLI exposes five stage tokens, `00` to `04`, with
  validation running as the first half of `04`.
- Stage 01 writes `01_backbone_filter/passed/`, a copy of the surviving backbones,
  and stage 02 reads that folder instead of re-deriving the selection from the table.
  It is rebuilt on every run, so a tightened rule cannot leave a rejected backbone
  behind for the GPU.
- `02_sequences/` holds the FASTA files directly; the `seqs/` level ProteinMPNN
  insists on is flattened away and the jsonl intermediates go to `logs/`.
- `03_predictions/` holds exactly two files per sequence. colabfold's own config and
  citation list move to `logs/colabfold/`.
- `filter_summary.json` now restates the rules it applied and lists the backbones
  that passed, so the output tree is readable without the config that produced it.
- Test fixtures renamed for the fold they contain rather than the verdict of a filter
  that no longer exists, and a third added: `rfdiffusion_c5_60aa_alpha_beta.pdb`,
  `_two_helix.pdb` and `_three_helix.pdb`, one per stage-01 verdict.

### Added
- Stage 00: symmetric backbone generation with RFdiffusion, with the per-chain to
  total contig-length conversion handled internally.
- Stage 01: objective backbone filtering. Twenty-nine geometric metrics are computed
  for every backbone and the pass/fail decision is expressed as config rules.
- Stage 02: ProteinMPNN sequence design with tied positions enforced.
- Stage 03: AlphaFold2 monomer prediction via localcolabfold in single-sequence mode.
- Stage 04: self-consistency RMSD and per-backbone success-rate rollup.
- Declarative filter engine with an AST whitelist, per-rule rejection accounting and
  startup syntax validation.
- Two-grain run manifest (backbones, sequences) with partial-upsert semantics and
  per-stage resume markers.
- Backend abstraction covering conda, Docker and Apptainer.
- Test suite of 169 checks against synthetic oligomers with known ground truth,
  plus three real RFdiffusion backbones as fixtures.

### Validated against real tools
- Stages 00 and 02 were executed against RFdiffusion and ProteinMPNN on a real GPU
  host. Symmetric output parses as five chains of 60 residues with continuous
  residue numbering, and tied-position design returns identical sequences on every
  chain. Both are now pinned by fixtures in `tests/data/`.
- Filter thresholds were recalibrated against 24 real backbones. The shipped
  defaults pass 71% of them and `configs/presets/compact_subunit.yaml` passes 33%.

### Changed
- Run directories hold only results. A 50-backbone run went from 467 MB to 42 MB by
  turning off RFdiffusion trajectory output, dropping .trb metadata and discarding
  the colabfold products nothing reads; both cleanups are config switches.
- Backbone filenames are continuous across batches: design_0 .. design_N-1 instead
  of a restart inside every batch.
- `run.seed` defaults to 1. ProteinMPNN reads `--seed 0` as a request for a random
  seed, so the previous default of 0 silently made sequence design irreproducible.

### Added (stage 05)
- `05_delivery/` collects every design that passed stage 04 as a pentamer PDB, a
  FASTA of its best sequence, and a row in `delivery.csv`.
- Backbone filtering is off by default. Stage 01 measures every backbone and passes
  all of them to sequence design; rules are opt-in. Metrics are recorded either way,
  so a filter can be applied afterwards by replaying stage 01 alone.
- Default run size is 50 backbones with 10 sequences each.
- A dry run no longer writes stage completion markers, which previously made the
  following real run skip the stages it had only simulated.
- The emitted RFdiffusion command now reproduces a verified hand-run invocation
  exactly: lowercase symmetry token, explicit `hydra.run.dir`, and no speculative
  keys. `inference.seed` does not exist in either shipped config and would have
  aborted the run; it is opt-in through `generate.seed_key`.
- Guiding potentials are off by default, matching that baseline invocation.
- GPU pinning is optional: `compute.gpus: null` leaves `CUDA_VISIBLE_DEVICES` unset.

### Added
- Gyration-tensor shape descriptors (`shape_anisotropy`, `asphericity`,
  `acylindricity`, `axis_ratio`) for the subunit, and the same set over the whole
  assembly, since the design target is a globular ring rather than a globular chain.
- `n_helices`, separating helix count from total secondary-structure element count.
- Optional `r_0` and `d_0` parameters on the oligomer-contact potential.

### Removed
- `configs/presets/globular_3helix.yaml`. It passed 1 backbone in 24 on real output
  and demanded a three-helix all-alpha fold that RFdiffusion rarely produces at these
  chain lengths. Replaced by `compact_subunit.yaml`.

### Not yet implemented
- Stage 03 execution: AlphaFold2 is not deployed yet, so only its command
  construction is covered.
- Fusion search and icosahedral assembly (stages 05-08).
- Container images.
