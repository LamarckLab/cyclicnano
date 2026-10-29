# Pipeline stages

Each stage reads the run manifest, does its work, writes results back, and drops a
marker in `.done/`. A stage with a marker is skipped on the next invocation unless
`--force` is passed, so an interrupted run resumes where it stopped.

```
00_backbones ──> 01_backbone_filter ──> 02_sequences ──> 03_predictions ──> 04_delivery
   GPU                   CPU                 GPU              GPU              CPU
```

Five directories for six steps. Self-consistency scoring sits between prediction and
delivery but writes no directory of its own: everything it produces is a column on
the backbone table, and the ranked view of those columns is the delivery table. It
runs as the first half of stage 04.

## 00 - Backbone generation

RFdiffusion in symmetric mode, invoked in batches of `generate.batch_size` with a seed
derived from `run.seed` plus the batch index.

The only subtlety is the contig length. The symmetric config of RFdiffusion expects
the **total** residue count across all chains, so a C5 of 60-residue subunits is
`[300-300]`. `Config.total_length` performs `sym_order x monomer_length`; the stage
never sees the per-chain number.

Output: `00_backbones/batch*_*.pdb`, one row per backbone in `tables/backbones.csv`.

## 01 - Backbone measurement and optional filtering

Computes 29 metrics per backbone (see [metrics.md](metrics.md)) and evaluates the
`backbone_filter.rules` expressions against them. Sets `bb_pass` and records
`bb_fail_reason` listing every rule that failed.

**The rule list is empty by default**, so nothing is rejected and every backbone
reaches sequence design. The metrics are recorded regardless, which is what makes
filtering a later decision rather than a commitment made before the run.

A backbone whose PDB cannot be parsed is marked failed rather than crashing the run.

Output:

- metrics on every backbone row of `tables/backbones.csv`;
- `01_backbone_filter/filter_summary.json`, restating the rules and giving the
  rejection count attributable to each one individually;
- `01_backbone_filter/passed/*.pdb`, a copy of the backbones that survived.

The copy is the handover to stage 02, which reads the folder rather than the table:
ProteinMPNN takes a directory, and one source of truth cannot fall out of step with
itself. It is rebuilt from scratch on every run, so tightening a rule removes the
backbones it now rejects instead of leaving them for the next stage to pick up.

This is the only stage worth re-running on its own. It is CPU-bound and takes seconds,
so changing a threshold and running `cyclicnano run --stages 01 --force` is cheap.

## 02 - Sequence design

ProteinMPNN over `01_backbone_filter/passed/`. Three calls:
`parse_multiple_chains.py`, then `make_tied_positions_dict.py --homooligomer 1`, then
`protein_mpnn_run.py`.

Tied positions are mandatory. Without them each chain of a Cn backbone is designed
independently and the resulting sequences cannot assemble; `design.tie_chains: false`
is rejected during `cyclicnano check`.

The output FASTA carries the chains joined by `/`. The parser drops the first entry
(the input sequence), splits on `/`, and records whether the chains came back
identical as a sanity check that tying actually took effect.

**A backbone's sequences depend on which other backbones were in the run.**
ProteinMPNN is invoked once over the whole folder with a single seed, so the position
of a backbone in that set determines where its sampling starts in the random stream.
Observed directly: the same backbone, byte-identical structure, had 15 of 20 sequences
fold back when it was one of 5 in the folder and 17 of 20 when it was one of 16. The
run is reproducible against another run with the same input set; `success_rate` is
therefore a property of the backbone *and* the batch it was designed in, not of the
backbone alone, and small differences in it should not be read as real.

Output: `02_sequences/design_N_mpnn_seqs.fa`, one per backbone, and one row per designed
sequence in `tables/sequences.csv`. ProteinMPNN writes into a `seqs/` subdirectory;
the stage lifts the files out of it and removes it, and the jsonl intermediates go to
`logs/`, so the stage directory holds FASTA files and nothing else.

## 03 - Structure prediction

localcolabfold in `--msa-mode single_sequence`, monomer only, templates off.

An MSA is not merely unnecessary for a de novo design, it is actively misleading: the
sequence has no homologues, so the search returns nothing useful and any confidence
derived from it is not interpretable. `cyclicnano check` flags any other `fold.msa_mode`.

Only the monomer is predicted. This validates that the designed sequence encodes the
designed fold; it does not validate that the subunits assemble into the ring. That
limitation is inherent to the protocol, not to this implementation.

Output: `03_predictions/` holding exactly two files per sequence, the rank-1 PDB and
its score JSON; `pred_path`, `mean_plddt`, `min_plddt` and `ptm` written back to the
sequence table. colabfold writes eight files per sequence and two run-level files of
its own; the plots, the single-sequence a3m and the PAE dumps are deleted once their
contents have been read, and its config and citation list move to `logs/colabfold/`.
Set `fold.keep_all_outputs: true` to keep everything.

## 04 - Self-consistency

Per sequence: CA RMSD between the predicted monomer and chain A of the design
backbone, via Kabsch superposition. A sequence passes when the RMSD is at or below
`validate.rmsd_cutoff` and its mean pLDDT is at or above `validate.plddt_cutoff`.

Per backbone: `success_rate` is the fraction of its sequences that pass, and is the
quantity the ranking is built on. `best_rmsd` and `median_rmsd` are recorded alongside.

Additional cutoffs listed in `validate.also_record` are stored as separate columns, so
a stricter or looser definition of success can be applied later without re-predicting
anything.

Output: columns on the backbone and sequence tables. No directory.

## 04 - Delivery

One folder per backbone that passed **stage 01**, not per backbone that passed the
RMSD cutoff. The cutoffs are uncalibrated, a strict shape filter leaves few survivors,
and a backbone whose twenty sequences all land at 1.4 A is a different kind of near
miss from one whose best is 6 A. Discarding both as failures throws that away, so the
verdict is recorded as a column and the ranking puts the ones that passed on top.

```
04_delivery/
├── delivery.csv              one row per backbone: success rate, best RMSD, verdict
├── delivery.json             the same table, for readers without pandas
└── <backbone_id>/
    ├── <backbone_id>.pdb     the design oligomer, as RFdiffusion wrote it
    ├── sequences.fasta       every sequence, annotated, the best one tagged BEST
    ├── summary.csv           one row per sequence: RMSD, pLDDT, pass
    └── predictions/          <sequence_id>.pdb and <sequence_id>_scores.json
```

Nothing here is new. It is the flat output of stages 02 and 03 regrouped by backbone
under short names, so that one candidate can be looked at without cross-referencing
three directories. The structure files are poly-glycine backbones as RFdiffusion wrote
them; the designed sequence is not threaded onto them.

## Not yet implemented

| Stage | Purpose |
|---|---|
| 05 | Fusion search: scan building-block pairs and helical fusion registers for geometries whose C5 and C3 axes meet at the icosahedral 37.3774 degree angle. |
| 06 | Assembly: expand the fused subunit by the 60 icosahedral rotations, then check for clashes and confirm the shell actually closes. |
| 07 | Junction redesign: redesign the fusion helix and the new inter-component contacts with the building-block bodies held fixed. |
| 08 | Reporting. |
