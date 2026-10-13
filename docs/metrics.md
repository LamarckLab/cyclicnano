# Backbone metrics

Every metric below is computed for every generated backbone in stage 01, whether or
not a filter rule references it, and stored in `tables/backbones.csv`. All of them are
deterministic functions of backbone coordinates and depend only on numpy.

Any metric name can be used in a `backbone_filter.rules` expression.

## Composition and shape

| Metric | Unit | Definition |
|---|---|---|
| `n_chains` | count | Chains parsed from the oligomer PDB. |
| `n_res_monomer` | count | Residues in chain A. |
| `ss_string` | - | Per-residue H/E/L assignment for chain A. |
| `helix_frac` | fraction | Fraction of chain A assigned H. |
| `strand_frac` | fraction | Fraction assigned E. |
| `loop_frac` | fraction | Fraction assigned L. |
| `n_sse` | count | Number of helix and strand runs. |
| `n_helices` | count | Number of helix runs alone, i.e. the subunit topology. |
| `loop_max_len` | residues | Longest contiguous loop run. |

Secondary structure is assigned from the CA trace using the P-SEA criteria (Labesse
et al., 1997): inter-CA distances at separations 2, 3 and 4, plus the virtual bond
angle and virtual torsion. Runs shorter than four residues (helix) or three (strand)
are demoted to loop.

This is an approximation. It is used in preference to DSSP so that the metric layer
needs no external binary and gives identical results on every machine, which matters
more here than the last few percent of assignment accuracy. A DSSP backend can be
added if exact agreement with published assignments is ever required.

**Why they filter.** Loop-dominated backbones are poorly designable: ProteinMPNN will
produce sequences for them, but those sequences rarely fold back. A long loop run is
a flexible segment that AlphaFold2 will place somewhere other than where the design
put it, inflating RMSD for reasons that have nothing to do with the fold.

## Compactness

| Metric | Unit | Definition |
|---|---|---|
| `rg` | angstrom | Radius of gyration of chain A CA atoms. |
| `rg_ratio` | ratio | `rg` divided by 2.2 x N^0.38, the compact-globular reference. |
| `contact_order` | fraction | Relative contact order over CA pairs within 8 A separated by at least 3 residues. |

`rg_ratio` near 1.0 is compact. Above roughly 1.3 the subunit is elongated: a two-helix
hairpin lands around 1.2 to 1.4 and a single straight helix above 1.5.

**Why they filter.** An elongated subunit has less buried core per residue and is
harder to design. High relative contact order correlates with slow and unreliable
folding, and in practice with low self-consistency in stage 04.

## Shape

| Metric | Unit | Definition |
|---|---|---|
| `asphericity` | ratio | Gyration-tensor asphericity, normalised by Rg squared. |
| `acylindricity` | ratio | Gyration-tensor acylindricity, normalised by Rg squared. |
| `shape_anisotropy` | 0 to 1 | Relative shape anisotropy (kappa squared): 0 for a sphere, 1 for a rod. |
| `axis_ratio` | ratio | Longest principal axis divided by the shortest. |

Computed from the eigenvalues of the gyration tensor of the chain-A CA coordinates.

`rg_ratio` measures overall size against a reference and cannot tell a compact ball
from an equally sized disc. `shape_anisotropy` measures the shape itself, which is
what the phrase "the subunit should be roughly spherical" actually means.

Reference values from the test fixtures:

| Shape | `shape_anisotropy` | `axis_ratio` |
|---|---|---|
| Ideal sphere | 0.005 | 1.1 |
| Three-helix bundle | 0.24 | 2.0 |
| Two-helix hairpin | 0.46 | 8.5 |
| Straight helix | 0.98 | 16.1 |

**Why they filter.** A compact three-helix subunit is the topology that works best as
an icosahedral building block, and these two metrics separate it from an elongated
hairpin far more cleanly than `rg_ratio` does.

Note the calibration: a helix bundle is inherently somewhat prolate, so a threshold
tight enough to look "spherical" in the everyday sense would reject exactly the folds
that are wanted. A cutoff of 0.30 accepts the bundle at 0.24 and rejects the hairpin
at 0.46. These reference values come from idealised synthetic fixtures; real
RFdiffusion output is less regular, so expect to loosen them against the per-rule
counts in `filter_summary.json`.

Neither metric is used by the shipped preset, which selects on secondary structure
instead. They remain computed and available as rules.

## Symmetry

| Metric | Unit | Definition |
|---|---|---|
| `sym_angle_deg` | degrees | Rotation angle relating chain A to chain B. |
| `sym_order_detected` | - | 360 divided by `sym_angle_deg`. |
| `sym_order_ok` | bool | Whether the detected order matches the requested Cn within 0.1. |
| `sym_rmsd` | angstrom | Largest CA RMSD between chain A and any other chain. |

The rotation is recovered by Kabsch superposition of chain A onto chain B. The
rotation axis is the eigenvector of that rotation matrix with eigenvalue 1, and the
angle follows from its trace.

**Why they filter.** This is a quality-control check on the generator, not on the
design. RFdiffusion enforces symmetry during denoising, but a run can still drift.
`sym_order_ok` catches an oligomer that came out with the wrong stoichiometry, and
`sym_rmsd` above a few tenths of an angstrom means the chains are not genuinely
equivalent, which makes tied-position sequence design in stage 02 incoherent.

## Packing and contacts

| Metric | Unit | Definition |
|---|---|---|
| `n_clash_intra` | count | Backbone heavy-atom pairs within one chain closer than 2.5 A, excluding neighbouring residues. |
| `n_clash_inter` | count | Backbone heavy-atom pairs between different chains closer than 2.5 A. |
| `n_clash` | count | Sum of the two. |
| `min_interchain_dist` | angstrom | Closest approach between backbone atoms of different chains. |
| `iface_contacts_per_chain` | count | Mean inter-chain CA-CA pairs within 8 A, per chain. |

Both clash terms use the same criterion: two backbone heavy atoms closer than 2.5 A
overlap. A CA-CA proxy is deliberately not used for the inter-chain case, because
ordinary interface contacts sit at 4 to 6 A and would all register as clashes.

**Why they filter.** `n_clash == 0` is the one rule that should never be relaxed; a
backbone with overlapping atoms is not a structure. `iface_contacts_per_chain` is the
opposite failure: subunits arranged on a ring so wide that they never touch. Both
failure modes occur in RFdiffusion symmetric output and neither is visible from the
sequence.

## Ring geometry

| Metric | Unit | Definition |
|---|---|---|
| `pore_radius` | angstrom | Smallest perpendicular distance from the symmetry axis to any CA. |
| `max_radius` | angstrom | Largest such distance. |
| `height` | angstrom | Extent along the symmetry axis. |

Distances are measured from the axis recovered in the symmetry block, through the
centroid of all CA atoms.

**Why they matter.** `pore_radius` is the central channel of the ring. It is not
filtered by default because the right value depends on the application, but it is the
metric to constrain if the cage is intended to encapsulate cargo or to exclude it.

## Recorded for the fusion stage

| Metric | Unit | Definition |
|---|---|---|
| `nterm_helix_len` | residues | Length of the leading helix run, 0 if the N terminus is not helical. |
| `cterm_helix_len` | residues | Length of the trailing helix run, 0 if the C terminus is not helical. |

These are not filtered in the current pipeline and exist only because the fusion stage
will need them.

Joining a C5 subunit to a C3 subunit into a single chain requires a rigid helical
fusion, not a flexible linker: with no designed interface between the two components,
the icosahedral geometry is held entirely by the fusion, and a flexible linker leaves
the relative orientation of the two oligomers undetermined. A rigid fusion needs a
terminal helix on each partner, long enough to overlap at several registers.

Recording these now means the fusion stage can be calibrated against backbones that
have already been generated and validated, instead of requiring a regeneration.

## Adding a metric

1. Compute it in `cyclicnano.geometry.backbone_metrics` and return it in the dict.
2. Add it to `test_metric_set_is_complete` in `tests/test_geometry.py`.
3. Add a test showing it separates a good fixture from a bad one.
4. Document it here, including the failure mode it is meant to catch.
5. Only then consider adding a rule for it to `configs/default.yaml`.

Step 5 is deliberately last. A metric is useful as soon as it is recorded; making it a
filter is a separate decision that should be justified by looking at its distribution
across a real run.
