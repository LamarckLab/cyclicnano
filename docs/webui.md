# Web browser

A read-only web page over a finished run: the funnel of each run, the delivery table
of each backbone, the metrics of each sequence, and the structures in a viewer that
can be rotated.

It submits nothing. No route starts a pipeline stage, touches a GPU or writes into
the output tree, so the page cannot disturb a run in progress and a careless click
cannot cost compute. Reading and running are separate concerns and this is the
reading half.

## Running it

```bash
conda activate lmk_Cyclicnano
pip install -e ".[webui]"                     # fastapi + uvicorn
python -m cyclicnano.webui.fetch_vendor       # once: the Mol* bundle

# read-only
python -m cyclicnano.webui --root /data/lmk/cyclicnano_outputs --port 8080

# read-only plus the submit form
python -m cyclicnano.webui --port 8080 \
    --allow-submit --profile configs/profiles/amax.yaml
```

`--root` defaults to `$CYCLICNANO_OUTPUTS`, then to `/data/lmk/cyclicnano_outputs`.
The default bind is `0.0.0.0`; pass `--host 127.0.0.1` to keep it on the machine.

**There is no authentication.** Anyone who can reach the port can read every run,
and with `--allow-submit` can start one. That is deliberate for a service on a lab network and wrong for anything else: do not
put it on a public address.

`fetch_vendor` downloads Mol*, pinned by version and checked against a hash, into
`static/vendor/`, which git ignores. Skipping it is not fatal — the page falls back to
a CDN — but then the viewer needs the browser to reach the open internet.

## What a run has to look like

A directory under the root counts as a run when it holds either
`tables/backbones.csv` or `04_delivery/delivery.json`. The table alone is enough on
purpose: running stages 00 and 01 by themselves is the cheap way to see whether the
filter leaves anything before spending GPU hours on the rest, and such a run has a
funnel and no delivery. Requiring the manifest made exactly that run invisible on the
page that would have shown its only result. Its card is marked 未到交付阶段 and its
backbone table says so rather than rendering empty, which would read as "everything
failed".

Directories whose name starts with an underscore are skipped. That is where pilots,
archived attempts and the job records live; some of them hold a delivery manifest from
an earlier attempt, and listing those beside the real runs would be the most
misleading thing the page could do.

Everything the page shows comes from two places:

| Shown | Read from |
| --- | --- |
| the funnel on the run card | `tables/backbones.csv`, `tables/sequences.csv` |
| the backbone table | `04_delivery/delivery.json` |
| the sequence metrics and the structure paths | `04_delivery/<folder>/sequences.json` |

The funnel comes from the tables rather than from counting the delivery folder,
because the tables still hold the backbones and sequences that were dropped, and the
attrition is the whole interest of a funnel. A stage that has not run shows an em
dash, not a zero: "nothing passed" and "not attempted" are different facts.

## The three levels

1. **Runs.** One card per run, each a funnel: generated, passed the shape filter,
   sequences designed, predictions made, sequences that folded back, backbones
   delivered.
2. **One run.** The delivery table, passing backbones first, every column sortable.
   Rank sorts as a number, so `10_design_169` does not land ahead of `2_design_230`.
3. **One backbone.** Its sequences with RMSD, pLDDT and ddG, beside a viewer.

The backbone table carries the whole distribution of each metric rather than a single
best value, which said almost nothing about twenty sequences. Each cell is a Gaussian
kernel density curve with the real values as a rug of ticks beneath it, and this
backbone's own minimum and maximum printed at the ends.

The track is shared by every row and only the end labels are per row. An axis rescaled
to each row's range would draw the same curve for a backbone spanning 0.4-0.6 and one
spanning 0.4-4.0, and the point of a column is reading down it; so the curve says where
on a common scale the sequences sit while the labels say how far they actually spread.
`PER_ROW_AXIS` in `static/app.js` flips that if the other reading is wanted.

The bandwidth is a fixed fraction of the shared axis, not Silverman's rule.
Deliberately: these distributions are bimodal, since a sequence either folds back near
0.5 angstrom or misses by three, and a per-row rule computed from the standard
deviation smooths exactly that into one hump - measured on `design_27` it gives 0.64
angstrom on a 4 angstrom axis, wide enough to merge the two lobes into a shape that is
not in the data. A fixed narrow kernel keeps the lobes apart and smooths every row the
same way.

The rug is there because a kernel estimate over twenty points can imply a smoothness
the data does not have. The curve is the shape; the ticks are the evidence, and one can
be checked against the other without leaving the cell. Height is normalised per row, so
a taller peak means "more concentrated here" within that backbone and not "more
sequences than that backbone"; the count is in the tooltip.

The axis bounds come from the measured runs - 1220 sequences and 371 energies - not
from a guess: RMSD spans 0.30 to 24.08 with a 99th percentile of 16.09, pLDDT 51.51 to
97.87, energy -77.01 to +0.85. pLDDT and energy are given axes that cover their data.
RMSD stops at 4 rather than 24 on purpose, because the decision is made at 1.0 and
resolution near the cutoff is worth more than room for a tail of designs that missed by
an order of magnitude. A value past either end is drawn at the edge and its true size is
printed as that row's maximum, so nothing is silently clipped.

Two success rates sit beside each other on level 2 and level 3: the monomer one, from
RMSD and pLDDT, and the interface one, the share of scored sequences below the ddG
cutoff. Only the first decides delivery.

## The viewer

Four kinds of structure. The three with files behind them are taken from
`sequences.json` rather than by globbing the folder, because that file is the
pipeline's statement of which structure belongs to which sequence:

| | What it is | Frame |
| --- | --- | --- |
| `design_<id>_bb.pdb` | the filtered RFdiffusion backbone, main chain only | design |
| `seqNNN_af2_C<n>.pdb` | the whole ring, assembled per request, no file on disk | design |
| `seqNNN_af2.pdb` | the no-MSA monomer prediction; B-factors are pLDDT | AlphaFold2 |
| `seqNNN_dimer_relaxed.pdb` | the relaxed dimer that `dG_separated` comes from | design |

The ring is what a sequence row opens on, because a sequence exists to make one and
a single subunit on its own does not show whether it does. It is the monomer
prediction placed on every chain of the design backbone - one Kabsch fit per
subunit, carried out when the browser asks for it, since it is derivable in
milliseconds and storing one per sequence would be tens of megabytes of the same
thing. Nothing in this pipeline predicts the assembly: the ring is the geometry the
design asks for, carrying the side chains the sequence implies, and it is not
evidence that the sequence assembles into it. The fit is rigid-body with no relax,
so neighbouring subunits may clash. The note under the viewer says all of this, in
both languages, because a picture of a pentamer invites exactly the reading it does
not support.

One structure at a time. The frame column is not a display option but a caveat the
note under the viewer repeats: the AlphaFold2 outputs come out in their own
coordinate frame, so their RMSD to the design is computed after superposition.
Measured on `design_27__seq011`, the monomer sits 19.8 angstrom from design chain A
as written, against 2.8 for the dimer.

The dimer is the one prediction in the designed frame, because stage 05 builds it by
superposing the predicted monomer onto chains A and B of the design. That
superposition is exact to a median 0.55 angstrom over the 43 evaluated C5 sequences;
the restrained relax then displaces the result by a median 2.7 angstrom as a
near-rigid body, changing its internal geometry by only 0.8 and the distance between
the two subunit centroids by between -0.6 and +0.9. So the interface `dG_separated`
describes is the designed one — what drifted is the whole object.

Mol*'s own panels are switched off: no left panel, no sequence track, no selection
mode. What is left is the viewport, its expand button, and two selects of this page's
own. The panels were more than this page needs, and the structure list has to live
here regardless, because Mol* cannot know about the fifty-odd files one backbone has.

Colouring, each with the scale legend under the viewport:

- **by chain**, using the blue ordinal ramp rather than one hue per chain. The chains
  of a cyclic homo-oligomer are identical copies in ring order, so their colour
  encodes a position, not an identity; light to dark is the way round the ring. The
  ramp also clears the ordinal gate on this viewport, which five distinct hues do
  not: past three, categorical hues stop separating under colour-vision deficiency.
- **by pLDDT**, in AlphaFold's own four bands (>= 90, 70-90, 50-70, < 50). This is the
  one multi-hue scale on the page, kept because every structural biologist reads
  those colours without a legend and a single-hue ramp here would be less legible,
  not more. It reads the B-factor column, which colabfold fills with per-residue
  confidence: verified against `seqNNN_af2_scores.json`, both span 79.12 to 98.38 on
  `design_27__seq011`. Disabled for the design backbone, which has no pLDDT.
- **by secondary structure**, coloured here rather than with 3Dmol's built-in
  `ssJmol` or `ssPyMOL` schemes, which paint helix magenta or red and would leave the
  legend describing colours the viewport is not showing.

The toolbar over the viewport resets the camera, toggles a slow spin, expands the
viewport to most of the window height, and saves a PNG of the current view.

A path that `sequences.json` names and disk lacks is shown as missing rather than
dropped from the list, so a truncated copy is visible instead of silent.

## Exporting sequences

The run page hands over the sequences as a file, because the next step after reading
the table is ordering something.

| Scope | What it holds |
| --- | --- |
| `delivered` | passing sequences from backbones whose success rate passed |
| `passed` | plus the ones that folded back on a backbone whose rate did not |
| `all` | everything in the delivery, failures included |

`delivered` is the default and the narrowest, because it is the one people act on. On
the measured C5 run that is 43 sequences across 4 backbones; `all` is 320 across 16.

FASTA carries its metrics in the description line, since the file travels on its own
and the numbers are what tell two records apart:

```
>design_27__seq000 C5 sc_rmsd=0.52 mean_plddt=94.2 dG_interface=-31.4 len=60
IEKLLKLYYYVELANTALILAEVTGSEKFLEIAERLAERIDEVDPELAELLRELAREIKE
```

The CSV is written with a byte-order mark so Excel opens it as UTF-8 on a double
click, which is the whole reason a spreadsheet was asked for. Both name the energy
column `dG_interface`, the field the pipeline stores, rather than the label the table
shows: a file that outlives this page has to stay traceable to where its numbers came
from.

## Submitting a run

Off unless the service was started with `--allow-submit`, and that flag is the whole
of the access control. There is no authentication, so whoever can reach the port can
spend GPU hours on a machine other people are using; making it a deployment decision
is the point. Without the flag the submit routes answer 403 and the page draws no form
at all rather than offering a button that fails.

Three rules hold the rest together.

**Every field is declared with a type and a range**, and nothing reaches the command
line that has not been through `validate`. The form the browser renders is built from
the same declaration the server enforces, so the two cannot drift. The command is
built as a list and spawned without a shell, so a value cannot become an argument:
`../evil` as a run name, `99999` backbones and `/etc/passwd` as a filter preset are
each refused with a 400 that names the field and says why.

**One job at a time.** The design server is shared, and two concurrent pipelines would
contend for the same cards. A submit while one is running gets a 409 naming the job
that holds the slot.

**A job gets its own session, so cancelling kills the whole tree.** Measured on a real
run: the process group held four processes, the CLI, conda, a shell and RFdiffusion
itself, and signalling only the parent would have left the last of those holding a
card. After a cancel the group is empty and the cards are back to 1 MiB. The record is
marked cancelled *before* the signal goes out, so the reaper, which only sees an exit
code, does not relabel a deliberate `-15` as a failure.

The stage selector is worth using. Stages 00 and 01 alone take minutes and say how
many backbones survive the filter; finding that out after ProteinMPNN and two rounds
of AlphaFold have run costs hours. A run submitted under a name that already exists
resumes it, skipping the stages that are done, which is also how a `00-01` run gets
extended to a full delivery without repeating the GPU work.

Job records live in `<root>/_jobs/<id>/`, one `job.json` and one `log.txt` each. The
record carries the exact invocation, so a run can be repeated from a terminal without
reconstructing it from the form.

### When the service restarts under a running job

The job keeps running, since it has its own session, but its exit code is lost with
the process that would have reaped it. Such a job is reported as state unknown rather
than guessed either way, and the note points at the pipeline's own `.done` marker,
which is the only evidence left. Calling it finished because the process is no longer
there would report a crash as a success.

## Two languages

English by default, with a toggle in the header. The Chinese side translates the prose
only: RMSD, pLDDT, ddG, REU, RFdiffusion, ProteinMPNN, AlphaFold2, Rosetta, FASTA, PDB
and the rest keep their usual spelling, because that is what they are called in the lab
and in the literature, and translating them would make the page harder to read.

Every string lives in `static/i18n.js`, keyed, with English as the fallback for a key
the other table is missing - so a half-finished translation shows English rather than a
key name. The structure labels and the notes under the viewer were moved out of the API
for this: `store.py` now returns a structure's *kind*, and the page decides what to
call it. An API that returns display text can only ever speak one language.

The theme toggle beside it is a sun or a moon, showing the mode a click would move to.
Both choices are remembered in `localStorage`, and the theme is stamped on the document
before first paint so a dark-mode reader never sees a white flash.

## When the viewer does not draw

This service is developed on a machine that cannot open the page, so a viewer that
fails to render leaves no console anyone can read. The page therefore reports on
itself. On a failure - and on demand, from the link under the viewer - it gathers
whether Mol* loaded, whether its `Viewer.create` exists, whether the browser has a
WebGL context and which renderer backs it, the host element's measured size and
whether it is in the document, and the caught exception with a short stack. That goes
on screen and to `POST /api/clientlog`, which writes it to the server log.

So the next failure can be diagnosed from the log rather than from a description.
Access-log evidence is worth reading first: if the browser never requested
`/api/runs/<run>/file`, the viewer failed before loading a structure, which points at
the plugin or at WebGL rather than at the data.

The endpoint writes to the log and nowhere else, and truncates its payload, because an
unauthenticated route that appends to a file is a way to fill a disk.

## On the name `ddG`

The column header and the number under it deliberately disagree, and this is the note
that keeps that honest.

The number is Rosetta's `dG_separated`: `InterfaceAnalyzerMover` takes the score of
the complex and subtracts the sum of the two separated and repacked chains, in REU.
That is a binding energy. It is **not** a ddG, which is the *difference between two*
binding energies and needs a reference state to subtract: a wild type against a
mutant, a parent against a redesign. These are de novo designs, so there is nothing to
subtract. And since REU is not kcal/mol, neither spelling is a free energy in the
thermodynamic sense; strictly it is a score difference.

The header says `ddG` anyway, because that is the name the group uses for it. The
tooltip on the column states what the number is regardless, so the label is a courtesy
and not a claim, and the exported files keep the stored field name `dG_interface` so
nothing that leaves this page carries the courtesy with it.

`DG_LABEL` at the top of `static/app.js` is the only line to change if the name is
ever settled. Nothing downstream reads it.

## Routes

| Route | Returns |
| --- | --- |
| `GET /api/health` | whether the output root is readable |
| `GET /api/runs` | every run with its funnel |
| `GET /api/runs/{run}/backbones` | the delivery table |
| `GET /api/runs/{run}/backbones/{id}` | one backbone, its sequences, its structures |
| `GET /api/runs/{run}/file?path=` | one file under `04_delivery/`, as `text/plain` |
| `GET /api/runs/{run}/file?path=&download=1` | the same bytes as a named attachment |
| `GET /api/runs/{run}/assembly?path=&design=` | the ring, built from those two files |
| `GET /api/runs/{run}/interface?path=&cutoff=` | the interface residues of a dimer |
| `GET /api/runs/{run}/export?format=&scope=` | the sequences as FASTA or CSV |
| `GET /api/jobs` | the slot, the history, and the form declaration |
| `POST /api/jobs` | start a run (403 without `--allow-submit`) |
| `GET /api/jobs/{id}` | one job, its state checked against the system |
| `GET /api/jobs/{id}/log?offset=` | the log from a byte offset |
| `POST /api/jobs/{id}/cancel` | signal the job's process group |
| `POST /api/clientlog` | what the browser could not say out loud |

`/api/docs` serves the generated OpenAPI page.

The log is followed by byte offset rather than by line, because the pipeline is
writing to the file while the page reads it and a line-based read would lose or
duplicate the partial line at the boundary.

The two forms of the file route exist because the viewer and a save want different
things from the same bytes. The viewer parses them itself, which wants `text/plain`;
a save wants a filename and a media type the browser will not rename. Serving a PDB
as `text/plain` to a download is what produced `.txt` files — the URL carries no
filename, so the browser invents one from the media type.

The file route is the one real attack surface, since its path arrives in a query
string. It is checked after resolution rather than by inspecting the string, so a
symlink pointing out of the tree is caught too, and the suffix must be `.pdb`,
`.json` or `.csv`. The delivery tree sits beside logs and configs on the same
machine, and that allowlist is what keeps a read-only viewer from becoming a file
server.

## Colour

Marks come from a validated palette rather than from taste. Every magnitude mark on
the page is one series, so it is one colour; the funnel uses a shared linear scale
across its six stages, which makes the last bars very short — that attrition is the
finding, and every count is printed outside its bar so no value is readable only as a
length. Pass and fail use the fixed status palette and always carry a glyph, so no
state is signalled by colour alone.

Dark mode is chosen rather than flipped: its own steps from the same ramps, declared
under both the OS media query and an explicit theme attribute so the header toggle
wins either way. The choice is remembered in `localStorage` and stamped before first
paint, so a dark-mode reader never sees a white flash.

## Testing

`tests/test_webui.py` covers the reading layer and the export with the standard
library alone, so it runs on a machine with no web stack. `tests/test_webui_jobs.py`
covers validation, command construction and the job state machine, spawning only a
harmless sleeper to check that a cancel really kills its process group.
`tests/test_webui_http.py` covers the routes and needs FastAPI and `httpx`; if FastAPI
is absent it stands down **and says so**, because `tests/run_all.py` counts any test
that does not raise as a pass and a silent skip would be indistinguishable from
coverage.

That warning is not hypothetical. The health route gained a field, and the test
asserting its exact shape kept passing on the development machine, where FastAPI is
not installed, then failed the moment the suite ran on the server. Run the suite where
the service runs.

## Operational notes

Nothing is cached. A run directory is a few hundred kilobytes of CSV and JSON, and a
stage that has just finished should appear on the next page load rather than when a
cache decides to expire.

`sudo` on the design server needs a password, so this cannot be a systemd unit. It
runs as a plain background process and **does not come back after a reboot**. To
survive one, add a user crontab entry:

```cron
@reboot . $HOME/mambaforge/etc/profile.d/conda.sh && conda activate lmk_Cyclicnano && python -u -m cyclicnano.webui --port 8080 >> /data/lmk/cyclicnano_webui/server.log 2>&1
```
