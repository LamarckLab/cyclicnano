# Changelog

All notable changes to this project are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added (web browser: the ring)
- A sequence row now opens on the whole designed oligomer - the pentamer or the
  trimer - instead of on the single predicted subunit. A sequence exists to make a
  ring, and a subunit on its own does not show whether it does.
- It is assembled when the browser asks for it, from the two files it needs: the
  monomer prediction laid onto every chain of the design backbone, one Kabsch fit per
  subunit. No file is stored, because the result is derivable in milliseconds and one
  per sequence would be tens of megabytes of the same thing. `GET
  /api/runs/{run}/assembly?path=&design=`, both paths resolved exactly as any other
  file, so neither can leave the delivery folder.
- The placement moved to `cyclicnano.assembly`, where stage 05 and the page share it.
  They have to: the energy and the picture then describe one geometry, and a test
  compares the two callers' output byte for byte to keep it that way.
- The note under the viewer says what the ring is not. Nothing in this pipeline
  predicts the assembly - ipTM went with the multimer prediction - so the ring is the
  geometry the design asks for carrying the side chains the sequence implies, the fit
  is rigid-body with no relax, and neighbouring subunits may clash. A picture of a
  pentamer invites the reading that the sequence makes one, which is exactly the
  claim this pipeline no longer tests.
- The monomer stays one click away in the same select, because RMSD and pLDDT were
  measured on it and nothing else.
- On the real runs the clash is not a quibble. Measured on design_27 seq011, adjacent
  subunits come within 0.51 A of each other heavy atom to heavy atom, with 585 pairs
  inside 5 A: they interpenetrate. The relaxed dimer stage 05 scores is built by the
  same placement and is credible only because a restrained relax follows it, which
  this does not have.
- Clicking a ddG still switches to the relaxed dimer the energy came from. It does not
  draw the interface side chains over it; the residue numbers below the viewer say the
  same thing in text, which is the level the question is asked at.

### Changed (web browser: the distributions are curves)
- The three distribution columns are Gaussian kernel density curves with a rug of the
  real values beneath, and each row's own minimum and maximum printed at the ends,
  replacing the bar histograms. Asked for, and better: the shape is readable at a
  glance and the spread is a number rather than an inference.
- The track stays shared across rows and only the labels are per row. An axis rescaled
  per row would draw the same curve for a backbone spanning 0.4-0.6 and one spanning
  0.4-4.0, and a column exists to be read downwards. `PER_ROW_AXIS` flips it.
- The bandwidth is a fixed fraction of the shared axis rather than Silverman's rule.
  These distributions are bimodal and a per-row rule from the standard deviation
  smooths that into one hump: on design_27 it gives 0.64 A on a 4 A axis, wide enough
  to merge two lobes into a shape that is not in the data. Verified by rendering the
  curves as text from the real runs, where design_230 shows its two RMSD lobes and its
  two pLDDT lobes separately.
- The rug is not decoration. A kernel estimate over twenty points can imply a
  smoothness the data does not have, so every measured value is drawn and the curve
  can be checked against it in the same cell.
- Axis bounds now come from the data rather than from a guess. Measured over 1220
  sequences and 371 energies: RMSD 0.30 to 24.08, pLDDT 51.51 to 97.87, energy -77.01
  to +0.85. The pLDDT axis was 40-100 and is now 50-100; the energy axis was -100 to
  +10 and is now -85 to +5; RMSD stays at 0-4 because the decision is at 1.0 and
  resolution near the cutoff beats room for a long tail.

### Changed (stage 05: no ipTM, a ddG verdict)
- The AlphaFold2-multimer prediction is gone, and with it ipTM. The group stopped
  trusting those predictions after finding clashes in the trimers and pentamers. The
  removal was clean, because the dimer the energy is measured on was always built from
  the stage-03 monomer: nothing about dG depended on the oligomer prediction. Stage 05
  is now **CPU-only**.
- What it costs is recorded in the stage docstring rather than left implicit: ipTM was
  the only number in this pipeline that asked whether a sequence assembles into the
  ring at all, and dG asks the narrower question of whether the designed interface
  scores well on a geometry that is assumed rather than predicted. A clashing multimer
  prediction is also not unambiguous evidence against the predictor - it can be
  evidence against the design - and the two cannot be told apart from the prediction
  alone.
- `interface.ddg_cutoff`, -25 REU, from the group. A sequence below it gets
  `interface_pass`; a backbone gets `n_interface_pass` and an `interface_success_rate`
  over the sequences that were *scored*, not over all of them, since the rest were
  already rejected on the monomer and counting them would report the monomer filter
  twice under another name.
- `final_pass` is unchanged: it still rests on the monomer success rate alone. The two
  rates sit beside each other so a dropped backbone can be read as having failed on
  the fold or on the interface, rather than merely as having failed.
- `interface_pass` is None when there is no energy, never False: "not measured" and
  "measured and above the cutoff" are different facts.
- Stage 05 no longer re-scores a sequence that already has an energy. The relax is
  about a minute and a half a structure, so re-running the stage to apply a changed
  cutoff would have spent hours reproducing numbers already on disk - 371 sequences on
  the measured runs. `--stages 05 --force` now re-reads the cutoff and sends only the
  unscored sequences to Rosetta.
- Fixed: stage 05 opened its Rosetta manifest inside `logs/`, which until now was
  created as a side effect of writing the multimer FASTA. Removing that prediction
  left the stage writing into a directory nothing had made.
- `ipTM` and the multimer confidences stay in `tables/sequences.csv` on runs that
  already hold them; nothing rewrites a table backwards.

### Changed (web browser: the table, the viewer, two languages)
- The backbone table shows the whole distribution of RMSD, pLDDT and ddG for each
  backbone's sequences, as sparkline histograms with the cutoff marked, in place of
  the single best RMSD and best pLDDT, which said almost nothing about twenty
  sequences. Bin edges are fixed in the script, not derived per row: per-row bins
  would give every row its own x-axis, so two histograms of the same shape would mean
  different things.
- The raw values are sent and binned in the browser, because the browser holds the
  cutoffs; sending bins would put the bin edges in one place and the thresholds they
  are read against in another. They come from the same `sequences.json` the detail
  page shows, so the two cannot disagree.
- Mol*'s own panels are switched off - no left panel, sequence track or selection mode
  - leaving the viewport, its expand button and two selects.
- The page is English by default, with a Chinese toggle. The Chinese side translates
  prose only and leaves RMSD, pLDDT, ddG, REU and the tool names alone. Strings live
  in `static/i18n.js` with English as the per-key fallback, so an incomplete
  translation shows English rather than a key name. `store.py` stopped returning
  display text for this and now returns a structure's kind: an API that returns
  labels can only speak one language.
- The theme toggle is a sun or a moon instead of the word "theme", showing the mode a
  click moves to.
- The energy column header is `ddG`, the name the group uses, while the stored field
  and the exported column stay `dG_interface`. The column tooltip states what the
  number actually is.
- Fixed: the ddG column's tooltip key did not exist, so it rendered as the literal
  string `col.ddg.title`. Caught by a check that every key the page asks for is
  present in both languages.

### Added (web browser: self-diagnosis)
- `POST /api/clientlog`. This service is developed on a machine that cannot open the
  page, so a viewer that fails to render left no console anyone could read. The page
  now reports whether Mol* loaded, whether `Viewer.create` exists, whether the browser
  has WebGL and which renderer backs it, the host element's measured size, and the
  caught exception - on screen and into the server log. It truncates its payload,
  because an unauthenticated route that appends to a file is a way to fill a disk.

### Added (web browser: submit and export)
- Sequence export from the run page, as FASTA or as CSV, in three scopes:
  `delivered` (passing sequences of passing backbones, 43 of them on the measured C5
  run), `passed` (plus the ones that folded back on a backbone whose overall rate did
  not) and `all` (320, failures included). FASTA carries its metrics in the
  description line because the file travels on its own; the CSV is written with a
  byte-order mark so Excel opens it as UTF-8 on a double click. Both name the energy
  column `dG_interface`, the stored field, so a file that outlives the page stays
  traceable.
- A submit form, off unless the service is started with `--allow-submit`. That flag is
  the whole of the access control: there is no authentication, so whoever can reach
  the port could otherwise spend GPU hours on a shared machine. Without it the submit
  routes answer 403 and the page draws no form.
- Every submitted field is declared with a type and a range, and the form the browser
  renders is built from the same declaration the server enforces, so the two cannot
  drift. The command is built as a list and spawned without a shell: `../evil` as a
  run name, 99999 backbones and `/etc/passwd` as a filter preset are each refused with
  a 400 naming the field.
- One job at a time, since two pipelines on this machine would contend for the same
  cards; a submit while one runs gets a 409 naming the job that holds the slot.
- A job runs in its own session and cancelling signals the process group. Measured on
  a real run the group held four processes - the CLI, conda, a shell and RFdiffusion
  itself - and signalling only the parent would have left the last of those holding a
  card. After a cancel the group is empty and the cards are back to 1 MiB. The record
  is marked cancelled before the signal, so the reaper does not relabel a deliberate
  `-15` as a failure.
- A stage selector, because stages 00 and 01 alone take minutes and say how many
  backbones survive the filter, while finding that out after ProteinMPNN and two
  rounds of AlphaFold costs hours. Resubmitting under an existing name resumes it.
- The log is followed by byte offset rather than by line: the pipeline writes to the
  file while the page reads it, and a line-based read would lose or duplicate the
  partial line at the boundary.

### Changed (web browser: the viewer)
- The viewer is Mol*, the one RCSB serves, with its own control panel and sequence
  track in place of a hand-rolled toolbar: that panel already does representations,
  colour themes, selections and measurements better. The structure list and a
  colour-theme shortcut stay on this side, because Mol* cannot know about the 55
  structures one backbone has.
- Its pLDDT theme is applicable here, verified rather than assumed: it needs defined
  B-factors on a model not marked experimental, and the colabfold predictions have the
  first and none of the records that would trigger the second. The four thresholds and
  colours were read out of the bundle, so the legend under the viewport is the colours
  the viewport is showing. The theme is withheld for the design backbone, whose
  B-factors are 0.00 throughout and which would otherwise render entirely as the
  "very low confidence" orange.
- Structure overlay removed; one structure at a time. The coordinate-frame caveat it
  existed to show is stated in the note under the viewer instead.
- A run is now recognised by `tables/backbones.csv` as well as by the delivery
  manifest. Requiring the manifest made a run stopped after stage 01 invisible on the
  page that would have shown its funnel, which is the only thing such a run has to
  show.

### Fixed (web browser: downloads)
- Downloading a structure saved a `.txt` file. The route served every request as
  `text/plain` and the URL carries no filename, so the browser invented one from the
  media type. `?download=1` now sends `application/octet-stream` with the real
  filename, while the viewer keeps `text/plain`.
- `_alive` treated only `ProcessLookupError` as "no such process". A pid that cannot
  exist raises a plain `OSError` on some platforms, which propagated out of a liveness
  probe and broke the job listing.

### Changed (web browser)
- The interface energy column is labelled `dG_separated`, Rosetta's own field name,
  rather than dG or ddG. It is the complex score minus the two separated and repacked
  chains, in REU: a binding energy, not a ddG, which would need a reference state to
  subtract and a de novo design has none. REU is also not kcal/mol, so neither
  spelling is a free energy. Naming the column after the field it holds asserts
  nothing false; `DG_LABEL` in `static/app.js` relabels it, and the stored field stays
  `dG_interface`.
- Structure overlay removed. One structure at a time; the coordinate-frame caveat it
  existed to make visible is now stated in the note under the viewer.
- Viewer rebuilt: five representations (cartoon, cartoon with side-chain sticks,
  backbone tube, ball and stick, van der Waals surface), three colourings, a toolbar
  over the viewport for camera reset, spin, expand and PNG export, and a soft outline
  view style.
- Chains are coloured with the blue ordinal ramp, not one hue per chain. The chains of
  a cyclic homo-oligomer are identical copies in ring order, so the colour encodes a
  position rather than an identity, and the ramp clears the ordinal gate on the dark
  viewport while five categorical hues do not separate under colour-vision deficiency
  (measured: the default first five slots fail the all-pairs normal-vision floor at
  dE 12.9).
- pLDDT colouring uses AlphaFold's four bands instead of a `roygb` rainbow, which was
  a multi-hue ramp for a magnitude. Secondary-structure colouring is applied with an
  explicit colour function rather than 3Dmol's `ssJmol`, which paints helix magenta
  and would have left the legend describing colours the viewport was not showing.
- Typography, spacing, surfaces and dark mode reworked against a validated palette.
  Every magnitude mark is one series and therefore one colour; the funnel is a shared
  linear scale with counts printed outside each bar; pass and fail carry a glyph so no
  state is colour alone. Dark mode is stepped for its own surface and stamped before
  first paint.

### Fixed (web browser)
- Downloading a structure saved a `.txt` file. The route served every request as
  `text/plain` and the URL carries no filename, so the browser invented one from the
  media type. `?download=1` now sends `application/octet-stream` with the real
  filename, while the viewer keeps `text/plain` because it parses the bytes itself.
- `_archive` and the other hidden directories answered 403 instead of 404, and the
  `startswith("_")` check in `run_dir` was unreachable because the name pattern
  rejected a leading underscore first. Hidden now means hidden: the pattern admits the
  name and the explicit check answers "no such run", so a browser is not told that
  `_archive` exists.

### Added (web browser)
- `cyclicnano.webui`, a read-only web page over a finished run: one card per run
  showing its funnel, the sortable delivery table, the per-sequence metrics, and a
  3D viewer. Optional extra, `pip install -e ".[webui]"`; the pipeline itself still
  imports no web framework. Run it with `python -m cyclicnano.webui`.
- It reads `04_delivery/` only, which stage 04 already made self-contained, so the
  browser never reaches into the stage directories. `webui/store.py` is standard
  library over that tree and is tested without a server.
- The funnel is read from `tables/` rather than counted from the delivery folder.
  The delivery holds only what survived stage 01, so counting it would report 16
  generated C5 backbones instead of 500 and make the filter look free. A stage that
  has not run renders as an em dash, not a zero.
- Four viewable structures per sequence, taken from `sequences.json` rather than
  globbed, and each flagged with whether it is in the design coordinate frame. Only
  the relaxed dimer is both a prediction and in that frame, so it is the only one the
  viewer will overlay on the backbone; the AlphaFold2 outputs come out in their own
  frame and overlaying them raw would show two objects apart for a reason unrelated
  to design quality.
- pLDDT colouring reads the B-factor column, which colabfold fills with per-residue
  confidence. Verified against `seqNNN_af2_scores.json`: both span 79.12 to 98.38 on
  `design_27__seq011`.
- Delivery folder names carry a rank prefix, so `10_design_169` sorts before
  `2_design_230` as a string. Rank is sorted as a number, and a backbone is located
  by id rather than by prefix, so a bookmarked URL survives a re-rank.
- The file route resolves paths after following symlinks and serves only `.pdb`,
  `.json` and `.csv`, so a query string cannot walk out of the delivery folder into
  the logs and configs beside it.

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
