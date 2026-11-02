/* cyclicnano browser: hash-routed.
 *
 *   #/                 the runs, each as a funnel
 *   #/<run>            the delivery table of one run
 *   #/<run>/<id>       one backbone: its sequences and its structures
 *   #/_jobs            submit a run, and the recent jobs
 *   #/_jobs/<id>       one job and its log
 *
 * The reading side is read-only and so is this. Every number shown is one the
 * pipeline already wrote, so a figure on screen can always be traced back to a
 * file under 04_delivery/.
 *
 * All display text goes through `t()` in i18n.js. Nothing here is
 * language-specific except the keys.
 */
'use strict';

/* The column header, and the field it actually holds.
 *
 * These deliberately disagree. The number is Rosetta's `dG_separated`: the score
 * of the complex minus the two separated and repacked chains, in REU. That is a
 * binding energy, not a ddG - a ddG is the difference between two of them and
 * needs a reference state to subtract, which a de novo design does not have.
 *
 * The header says ddG anyway, because that is the name the group uses. The
 * tooltip states what the number is regardless, so the label is a courtesy and
 * not a claim, and the exported files keep the stored field name. If the name is
 * ever settled, DG_LABEL is the only line to change. */
const DG_LABEL = 'ddG';

/* ------------------------------------------------------------------ helpers */

function el(tag, attrs, ...kids) {
  const node = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k === 'class') node.className = v;
    else if (k === 'text') node.textContent = v;
    else if (k.startsWith('on')) node.addEventListener(k.slice(2), v);
    else node.setAttribute(k, v === true ? '' : v);
  }
  for (const kid of kids.flat()) {
    if (kid === null || kid === undefined || kid === false) continue;
    node.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  }
  return node;
}

function icon(paths) {
  const ns = 'http://www.w3.org/2000/svg';
  const root = document.createElementNS(ns, 'svg');
  root.setAttribute('viewBox', '0 0 16 16');
  root.setAttribute('aria-hidden', 'true');
  for (const d of paths) {
    const node = document.createElementNS(ns, 'path');
    node.setAttribute('d', d);
    root.append(node);
  }
  return root;
}

const SUN = ['M8 5.1a2.9 2.9 0 1 0 0 5.8 2.9 2.9 0 0 0 0-5.8z',
             'M8 .9v2M8 13.1v2M.9 8h2M13.1 8h2',
             'M2.98 2.98l1.4 1.4M11.62 11.62l1.4 1.4M2.98 13.02l1.4-1.4M11.62 4.38l1.4-1.4'];
const MOON = ['M13.4 10.6A5.9 5.9 0 0 1 5.4 2.6a5.7 5.7 0 1 0 8 8z'];

const isNil = (v) => v === null || v === undefined || v === '' ||
                     (typeof v === 'number' && !Number.isFinite(v));

const fixed = (v, d) => (isNil(v) ? '—' : Number(v).toFixed(d));

/** A number cell. Missing is an em dash, never a zero: "not measured" and
 *  "measured as zero" are different facts and these tables hold both. */
function numCell(value, digits) {
  const td = el('td', {class: 'num'});
  if (isNil(value)) td.append(el('span', {class: 'nil', text: '—'}));
  else td.textContent = Number(value).toFixed(digits === undefined ? 2 : digits);
  return td;
}

/** A pass/fail cell. The glyph and the column name carry the state as well as
 *  the colour, so it is never colour alone. */
function flagCell(value) {
  if (isNil(value)) return el('td', {}, el('span', {class: 'nil', text: '—'}));
  return el('td', {}, el('span', {class: 'flag ' + (value ? 'yes' : 'no'),
                                  text: value ? '✓' : '✗'}));
}

/** A rate, as a thin bar beside a count — never behind it.
 *
 *  The count rather than the ratio, because the ratio alone leaves out the
 *  denominator and the denominator moves: a monomer rate is over all twenty
 *  sequences, while an interface rate is over however many were scored, which
 *  on the measured runs ranges from 6 to 20. "5/20" and "17/17" are both 1.00
 *  and 0.25 away from each other, and reading that off a bar length is not the
 *  same as being told.
 */
function meterCell(rate, passed, total) {
  const td = el('td', {class: 'num'});
  if (isNil(rate) && isNil(passed)) {
    td.append(el('span', {class: 'nil', text: '—'}));
    return td;
  }
  const width = isNil(rate)
    ? (total ? (passed / total) * 100 : 0)
    : Math.max(0, Math.min(1, rate)) * 100;
  const count = (isNil(passed) || isNil(total))
    ? (isNil(rate) ? '—' : Number(rate).toFixed(2))
    : `${passed}/${total}`;
  td.append(el('span', {class: 'meter', title: isNil(rate) ? null : Number(rate).toFixed(3)},
    el('span', {class: 'track'}, el('span', {class: 'fill', style: `width:${width}%`})),
    el('span', {class: 'count', text: count})));
  return td;
}

async function getJSON(url) {
  const r = await fetch(url, {cache: 'no-store'});
  if (!r.ok) {
    let detail = r.statusText;
    try { detail = (await r.json()).detail || detail; } catch (e) { /* not JSON */ }
    throw new Error(`${r.status} ${detail}`);
  }
  return r.json();
}

async function postJSON(url, body) {
  const r = await fetch(url, {
    method: 'POST',
    headers: {'Content-Type': 'application/json'},
    body: JSON.stringify(body || {}),
  });
  let payload = null;
  try { payload = await r.json(); } catch (e) { /* not JSON */ }
  if (!r.ok) throw new Error((payload && payload.detail) || `${r.status} ${r.statusText}`);
  return payload;
}

const fileURL = (run, path, download) =>
  `/api/runs/${encodeURIComponent(run)}/file?path=${encodeURIComponent(path)}` +
  (download ? '&download=1' : '');

const short = (id) => (id || '').split('__').pop() || t('none');

/* ------------------------------------------------------- mini distributions */

/* The axis is shared by every row; only the end labels are per row.
 *
 * The point of a column is reading down it, and an axis rescaled to each row's
 * own range would draw the same curve for a backbone spanning 0.4-0.6 and one
 * spanning 0.4-4.0, which is worse than no chart. So the track is fixed and the
 * numbers at the ends are this backbone's own minimum and maximum: the curve
 * says where on a common scale the sequences sit, the labels say how far they
 * actually spread. The fixed span is in the column header's tooltip.
 *
 * To make the axis follow each row instead, set PER_ROW_AXIS.
 */
const PER_ROW_AXIS = false;

/* Measured across both delivered runs, 1220 sequences and 371 energies:
 *   sc_rmsd        0.30 … 24.08, 99th percentile 16.09
 *   mean_plddt    51.51 … 97.87, 1st percentile 59.41
 *   dG_interface  -77.01 … +0.85
 * The pLDDT and energy axes therefore cover their data. The RMSD axis stops at
 * 4 on purpose rather than at 24: the decision is made at 1.0, so resolution
 * near the cutoff is worth more than room for a tail of designs that missed by
 * an order of magnitude. A value past either end is drawn at the edge and its
 * true size is printed as the row's maximum, so nothing is silently clipped. */
const DISTRIBUTIONS = {
  sc_rmsd: {lo: 0, hi: 4, better: 'low', digits: 2, unit: 'Å'},
  mean_plddt: {lo: 50, hi: 100, better: 'high', digits: 1, unit: ''},
  dG_interface: {lo: -85, hi: 5, better: 'low', digits: 1, unit: 'REU'},
};

/* One bandwidth for every row, a fixed fraction of the shared axis.
 *
 * Not Silverman's rule, deliberately. These distributions are bimodal - a
 * sequence either folds back near 0.5 A or misses by 3 - and a per-row rule
 * computed from the standard deviation smooths exactly that into one blob:
 * measured on design_27 it gives 0.64 A on a 4 A axis, wide enough to merge the
 * two lobes into a single hump that is not in the data. A fixed, narrow kernel
 * keeps the lobes apart, and keeps every row smoothed the same way so the
 * shapes can be compared.
 */
const KDE_BANDWIDTH = 0.04;          // of the axis span
const KDE_SAMPLES = 64;

/** A density curve over the shared axis, with the real values as a rug beneath.
 *
 *  The curve is the shape; the rug is the evidence. A kernel estimate over
 *  twenty points can imply a smoothness the data does not have, so every value
 *  is also drawn as a tick: what the curve claims can be checked against what
 *  was measured without leaving the cell.
 *
 *  Height is normalised per row, so a taller peak means "more concentrated
 *  here" within this backbone and not "more sequences than that backbone". The
 *  count is in the tooltip, and the rate columns beside it carry the counts.
 */
function sparkDensity(field, values, cutoff) {
  const spec = DISTRIBUTIONS[field];
  const all = values || [];
  const numbers = all.filter((v) => !isNil(v)).map(Number).sort((a, b) => a - b);
  const missing = all.length - numbers.length;

  if (!numbers.length) {
    const why = missing ? `n=0, ${missing} not measured` : 'n=0';
    return el('span', {class: 'nil', title: why, text: '—'});
  }

  const least = numbers[0];
  const most = numbers[numbers.length - 1];

  // A per-row axis still needs a width when every value is identical.
  let lo = spec.lo;
  let hi = spec.hi;
  if (PER_ROW_AXIS) {
    const pad = Math.max((most - least) * 0.08, Math.abs(most || 1) * 0.02);
    lo = least - pad;
    hi = most + pad;
  }
  const span = (hi - lo) || 1;

  const W = 100;
  const H = 24;
  const at = (v) => ((v - lo) / span) * W;
  const clampX = (x) => Math.max(0, Math.min(W, x));

  // The kernel is summed over the true values, not over clamped ones, so a
  // value past the end of the axis still pulls the curve towards that edge
  // rather than piling up exactly on it.
  const h = span * KDE_BANDWIDTH;
  const ys = [];
  for (let i = 0; i <= KDE_SAMPLES; i++) {
    const x = lo + (span * i) / KDE_SAMPLES;
    let density = 0;
    for (const v of numbers) {
      const z = (x - v) / h;
      density += Math.exp(-0.5 * z * z);
    }
    ys.push(density);
  }
  const peak = Math.max(...ys) || 1;

  const points = ys.map((d, i) => {
    const x = (W * i) / KDE_SAMPLES;
    const y = H - (d / peak) * (H - 1.5);
    return `${x.toFixed(2)},${y.toFixed(2)}`;
  });

  const ns = 'http://www.w3.org/2000/svg';
  const svg = document.createElementNS(ns, 'svg');
  svg.setAttribute('viewBox', `0 0 ${W} ${H + 5}`);
  svg.setAttribute('preserveAspectRatio', 'none');
  svg.setAttribute('class', 'spark-svg');
  svg.setAttribute('aria-hidden', 'true');

  const add = (tag, attrs) => {
    const node = document.createElementNS(ns, tag);
    for (const [k, v] of Object.entries(attrs)) node.setAttribute(k, v);
    svg.append(node);
    return node;
  };

  // The side of the axis that passes, washed in the curve's own hue: it
  // annotates the scale rather than encoding a second quantity.
  if (!isNil(cutoff)) {
    const edge = clampX(at(Number(cutoff)));
    const x0 = spec.better === 'low' ? 0 : edge;
    const width = spec.better === 'low' ? edge : W - edge;
    if (width > 0.5) {
      add('rect', {x: x0, y: 0, width, height: H, class: 'spark-pass'});
    }
  }

  add('path', {d: `M0,${H} L${points.join(' L')} L${W},${H} Z`, class: 'spark-fill'});
  add('path', {d: `M${points.join(' L')}`, class: 'spark-line'});

  // Every measured value, so the curve can be checked against the data.
  for (const v of numbers) {
    const x = clampX(at(v)).toFixed(2);
    add('line', {x1: x, x2: x, y1: H + 1, y2: H + 4.5, class: 'spark-tick'});
  }

  if (!isNil(cutoff)) {
    const x = clampX(at(Number(cutoff))).toFixed(2);
    add('line', {x1: x, x2: x, y1: 0, y2: H + 4.5, class: 'spark-cut'});
  }

  const unit = spec.unit ? ` ${spec.unit}` : '';
  const summary = `n=${numbers.length}`
    + (missing ? `, ${missing} not measured` : '')
    + `, ${least.toFixed(spec.digits)} … ${most.toFixed(spec.digits)}${unit}`
    + `, median ${numbers[Math.floor(numbers.length / 2)].toFixed(spec.digits)}`
    + (isNil(cutoff) ? '' : `, cutoff ${Number(cutoff).toFixed(spec.digits)}${unit}`)
    + (PER_ROW_AXIS ? '' : `, axis ${spec.lo} … ${spec.hi}${unit}`);

  return el('div', {class: 'spark', title: summary},
    svg,
    el('div', {class: 'spark-axis'},
      el('span', {text: least.toFixed(spec.digits)}),
      el('span', {text: most.toFixed(spec.digits)})));
}

/* -------------------------------------------------------------------- state */

const state = {
  run: null,
  backboneId: null,
  detail: null,
  sort: {},              // one {key, dir} per table, kept across repaints
  seqId: null,
  structPath: null,
  colorMode: 'chain',
  viewer: null,
  cutoffs: {},           // of the run being looked at
};

/* Nodes of the detail view that are repainted without rebuilding the page. */
const nodes = {};

/* The element every page paints into. `stale` is how a render that has been
 * overtaken stops before it replaces a newer page's contents. */
const view = () => document.getElementById('view');

function crumbs(...parts) {
  const nav = document.getElementById('crumbs');
  nav.replaceChildren();
  parts.forEach((p, i) => {
    if (i) nav.append(el('span', {class: 'sep', text: '/'}));
    nav.append(p.href ? el('a', {href: p.href, text: p.text}) : el('span', {text: p.text}));
  });
}

/** Show a failure, and say so where it can be read.
 *
 *  The page is developed on a machine that cannot open it, so an error that
 *  only ever appears in a div is an error nobody here learns about. This one
 *  cost a round trip: `controls()` referenced a `current` that an earlier
 *  rewrite had removed, and the page said so on screen while the server log
 *  said nothing at all.
 */
function fail(err) {
  // No token guard: this is the one paint that must always happen. A render
  // that failed has nothing newer to be overtaken by, and swallowing the
  // message would leave a blank page.
  view().replaceChildren(el('div', {class: 'error', text: String(err.message || err)}));
  viewerError = err instanceof Error ? err : new Error(String(err));
  report('render-failed');
}

/* --------------------------------------------------------- sortable tables */

/** A table whose headers sort it.
 *
 *  Sorting happens in the browser over rows already loaded, which is why the
 *  server hands over a whole delivery table in one response: a run has tens of
 *  backbones, not thousands, and a sort that needs a round trip is a sort
 *  nobody uses. A column with `value: false` is not sortable, which is how the
 *  histogram columns opt out: a distribution has no one number to order by.
 */
function sortableTable(id, columns, rows, opts) {
  opts = opts || {};
  const sort = state.sort[id] || (state.sort[id] = {key: null, dir: 1});

  if (sort.key) {
    const col = columns.find((c) => c.key === sort.key);
    if (col && col.value !== false) {
      const valueOf = (r) => (col.value ? col.value(r) : r[col.key]);
      rows = rows.slice().sort((a, b) => {
        const x = valueOf(a), y = valueOf(b);
        if (isNil(x) && isNil(y)) return 0;
        if (isNil(x)) return 1;                        // missing values sink
        if (isNil(y)) return -1;
        if (typeof x === 'number' && typeof y === 'number') return (x - y) * sort.dir;
        if (typeof x === 'boolean') return (x === y ? 0 : x ? -1 : 1) * sort.dir;
        return String(x).localeCompare(String(y), 'en') * sort.dir;
      });
    }
  }

  const head = el('tr');
  for (const col of columns) {
    const active = sort.key === col.key;
    const label = t(col.label);
    const sortable = col.value !== false;
    head.append(el('th', {
      class: (col.left ? 'left' : 'num') + (sortable ? '' : ' nosort'),
      title: col.title ? t(col.title) : label,
      onclick: sortable ? () => {
        if (active) sort.dir = -sort.dir;
        else { sort.key = col.key; sort.dir = col.desc ? -1 : 1; }
        (opts.onsort || render)();
      } : null,
    }, label, active ? el('span', {class: 'caret', text: sort.dir > 0 ? ' ▲' : ' ▼'}) : null));
  }

  const body = el('tbody');
  for (const row of rows) {
    const tr = el('tr', {
      class: opts.selected && opts.selected(row) ? 'selected' : null,
      onclick: opts.onclick ? () => opts.onclick(row) : null,
    });
    for (const col of columns) tr.append(col.cell(row));
    body.append(tr);
  }

  return el('div', {class: 'tablewrap'}, el('table', {}, el('thead', {}, head), body));
}

/* ------------------------------------------------------------ level 1: runs */

const FUNNEL = [
  {key: 'n_generated', label: 'funnel.generated'},
  {key: 'n_passed_filter', label: 'funnel.filtered'},
  {key: 'n_sequences', label: 'funnel.sequences'},
  {key: 'n_predictions', label: 'funnel.predictions'},
  {key: 'n_sequences_pass', label: 'funnel.folded'},
  {key: 'n_final_pass', label: 'funnel.delivered', final: true},
];

/** The funnel of one run.
 *
 *  One series, so one colour: the bar length carries the count and the row
 *  position carries the stage. The scale is linear and shared across all six
 *  rows, which makes the last bars very short — that attrition is the finding,
 *  not a rendering problem, and the count is printed outside every bar so no
 *  value is readable only as a length.
 */
function funnelChart(run) {
  const max = Math.max(...FUNNEL.map((s) => run[s.key] || 0), 1);
  const chart = el('div', {class: 'funnel'});
  let previous = null;

  for (const stage of FUNNEL) {
    const n = run[stage.key];
    const kept = (!isNil(n) && !isNil(previous) && previous)
      ? t('funnel.kept', (100 * n / previous).toFixed(1)) : '';
    const tip = t(stage.label + '.title') + kept;
    const mark = stage.final ? ' row-final' : '';
    chart.append(
      el('div', {class: 'k' + mark, title: tip, text: t(stage.label)}),
      el('div', {class: 'track' + mark, title: tip},
         el('div', {class: 'fill', style: `width:${isNil(n) ? 0 : (100 * n / max)}%`})),
      el('div', {class: 'n' + mark + (isNil(n) ? ' nil' : ''), title: tip,
                 text: isNil(n) ? '—' : String(n)}));
    if (!isNil(n)) previous = n;
  }
  return chart;
}

/* ------------------------------------------------- level 2: one run's table */

/** One project card. The funnel where there is one, and what to do where there
 *  is not: a project with a name and no results is a real state, and an empty
 *  funnel of six em dashes says less than a sentence does. */
function projectCard(run) {
  const badges = [
    run.symmetry ? el('span', {class: 'badge sym', text: run.symmetry}) : null,
    el('span', {class: 'badge ' + (run.complete ? 'done' : ''),
                text: run.complete ? t('runs.complete') : t('runs.incomplete')}),
    run.started && !run.delivered
      ? el('span', {class: 'badge warn', title: t('runs.undelivered.title'),
                    text: t('runs.undelivered')}) : null,
  ];
  // A project whose own tables measured a different symmetry than it declared is
  // worth seeing rather than quietly resolving in favour of one of them.
  if (run.symmetry_declared && run.symmetry && run.symmetry_declared !== run.symmetry) {
    badges.push(el('span', {class: 'badge bad',
      text: t('projects.mismatch', run.symmetry_declared, run.symmetry)}));
  }

  const bits = [];
  if (!isNil(run.n_generated) && !isNil(run.n_final_pass) && run.n_generated) {
    bits.push(t('runs.rate', (100 * run.n_final_pass / run.n_generated).toFixed(1)));
  }
  if (!isNil(run.n_interface)) bits.push(t('runs.interface', run.n_interface));
  if (run.created) {
    bits.push(t('projects.created', new Date(run.created * 1000).toLocaleDateString()));
  }
  bits.push(t('runs.updated', new Date(run.updated * 1000).toLocaleString()));

  return el('a', {class: 'card', href: `#/${encodeURIComponent(run.run)}`},
    el('div', {class: 'card-head'},
      el('span', {class: 'name', text: run.title || run.run}),
      ...badges.filter(Boolean)),
    run.title && run.title !== run.run
      ? el('div', {class: 'cardpath', text: run.run}) : null,
    run.started
      ? funnelChart(run)
      : el('div', {class: 'emptyproject'},
          el('div', {class: 'emptytag', text: t('projects.empty')}),
          el('div', {class: 'hint', text: t('projects.empty.hint')})),
    run.note ? el('div', {class: 'cardnote', text: run.note}) : null,
    el('div', {class: 'aside', text: bits.join(' · ')}));
}

async function renderProjects() {
  const token = renderToken;
  crumbs({text: t('projects.title')});
  // The job service is asked about too, because whether a project can be
  // created at all is the same switch as whether a job can be started.
  const [runs, info] = await Promise.all([
    getJSON('/api/runs'),
    getJSON('/api/jobs').catch(() => ({enabled: false})),
  ]);

  const head = el('div', {class: 'pagehead'},
    el('div', {},
      el('h1', {text: t('projects.title')}),
      el('p', {class: 'sub', text: t('projects.sub')})),
    info.enabled
      ? el('a', {class: 'btnlink', href: '#/_new', text: t('projects.new')}) : null);

  if (!runs.length) {
    if (!isCurrentRender(token)) return;
  view().replaceChildren(head,
      el('div', {class: 'panel notice'},
        el('p', {text: t('projects.none')}),
        el('p', {class: 'hint', text: t('projects.none.detail')})));
    return;
  }

  if (!isCurrentRender(token)) return;
  view().replaceChildren(head,
    el('div', {class: 'runs'}, runs.map((run) => projectCard(run))));
}

async function renderNewProject() {
  const token = renderToken;
  crumbs({text: t('projects.title'), href: '#/'}, {text: t('new.title')});
  const info = await getJSON('/api/jobs').catch(() => ({enabled: false}));

  if (!info.enabled) {
    if (!isCurrentRender(token)) return;
  view().replaceChildren(
      el('h1', {text: t('new.title')}),
      el('div', {class: 'panel notice'},
        el('p', {text: t('new.disabled')}),
        el('p', {class: 'hint', text: t('jobs.disabled.why')})));
    return;
  }

  // The symmetry choices come from the server's own field spec, so the form
  // cannot offer a value the server would refuse.
  const symmetryField = (info.spec.fields || []).find((f) => f.name === 'symmetry');
  const choices = (symmetryField && symmetryField.choices) || ['C3', 'C5'];

  const draft = {name: '', title: '', symmetry: 'C5', monomer_length: 60, note: ''};
  const problem = el('div', {});

  function row(key, labelKey, input) {
    const id = `np_${key}`;
    input.id = id;
    const bind = () => { draft[key] = input.value; };
    input.addEventListener('input', bind);
    input.addEventListener('change', bind);
    return el('div', {class: 'formrow'},
      el('label', {for: id, text: t(labelKey)}),
      input,
      el('div', {class: 'fieldnote', text: t(`${labelKey}.help`)}));
  }

  const create = el('button', {class: 'primary', type: 'submit', text: t('new.create')});

  const form = el('form', {class: 'panel jobform', onsubmit: async (e) => {
    e.preventDefault();
    create.disabled = true;
    problem.replaceChildren();
    try {
      const made = await postJSON('/api/runs', draft);
      location.hash = `#/${encodeURIComponent(made.name)}`;
    } catch (err) {
      problem.replaceChildren(el('div', {class: 'error', text: err.message}));
      create.disabled = false;
    }
  }},
    el('div', {class: 'formgrid'},
      row('name', 'new.name',
          el('input', {type: 'text', autocomplete: 'off', required: true})),
      row('title', 'new.label', el('input', {type: 'text', autocomplete: 'off'})),
      row('symmetry', 'new.symmetry',
          el('select', {}, choices.map((c) => el('option', {
            value: c, selected: c === draft.symmetry, text: c})))),
      row('monomer_length', 'new.monomer',
          el('input', {type: 'number', value: 60, min: 20, max: 400, step: 1})),
      row('note', 'new.note', el('input', {type: 'text', autocomplete: 'off'}))),
    problem,
    el('div', {class: 'controls'}, create,
       el('a', {class: 'btnlink', href: '#/', text: t('new.cancel')})));

  if (!isCurrentRender(token)) return;
  view().replaceChildren(
    el('h1', {text: t('new.title')}),
    el('p', {class: 'sub', text: t('new.sub')}),
    form);
}

/** The submit panel a project carries, with the run name and symmetry fixed.
 *
 *  One project is one run directory, so those two are the project rather than
 *  parameters of the job: letting the form change them would let a submission
 *  made from inside a project write somewhere else.
 */
function submitPanel(meta, info, open) {
  if (!info.enabled) return null;
  const fixed = {name: meta.run};
  if (meta.symmetry) fixed.symmetry = meta.symmetry;
  return el('details', {class: 'panel submitpanel', open: open ? true : null},
    el('summary', {text: t('jobs.submit.here')}),
    el('p', {class: 'hint', text: t('jobs.submit.locked')}),
    jobFormPanel(info, fixed));
}

async function renderRun(run) {
  const token = renderToken;
  crumbs({text: t('projects.title'), href: '#/'}, {text: run});
  const [meta, rows, info] = await Promise.all([
    getJSON(`/api/runs/${encodeURIComponent(run)}`),
    getJSON(`/api/runs/${encodeURIComponent(run)}/backbones`),
    getJSON('/api/jobs').catch(() => ({enabled: false})),
  ]);
  if (info.spec) formSpec = info.spec;

  const subtitle = [meta.run,
                    meta.symmetry || null,
                    meta.monomer_length_declared
                      ? t('run.sub.monomer', meta.monomer_length_declared).trim()
                      : null].filter(Boolean).join(' · ');

  // Nothing has run yet, so the funnel and the table have nothing to say and the
  // form is the whole page.
  if (!meta.started) {
    if (!isCurrentRender(token)) return;
  view().replaceChildren(
      el('h1', {text: meta.title || run}),
      el('p', {class: 'sub', text: subtitle}),
      meta.note ? el('p', {class: 'sub', text: meta.note}) : null,
      el('div', {class: 'panel notice'},
        el('p', {text: t('projects.empty')}),
        el('p', {class: 'hint', text: t('projects.empty.hint')})),
      submitPanel(meta, info, true) || submitDisabledPanel());
    return;
  }

  if (!rows.length) {
    if (!isCurrentRender(token)) return;
  view().replaceChildren(
      el('h1', {text: t('run.nodelivery.title', meta.title || run)}),
      el('p', {class: 'sub', text: t('run.nodelivery.sub')}),
      el('div', {class: 'panel notice'},
        el('p', {class: 'hint', text: t('run.nodelivery.hint')})),
      submitPanel(meta, info, true) || submitDisabledPanel());
    return;
  }

  // The thresholds the curves are read against live in each backbone's
  // sequences.json, which the list response does not carry. They are a property
  // of the run, so one is fetched and used for every row.
  try {
    const sample = await getJSON(`/api/runs/${encodeURIComponent(run)}/backbones/`
                                 + encodeURIComponent(rows[0].backbone_id));
    state.cutoffs = sample.criteria || {};
  } catch (err) {
    state.cutoffs = {};                                // the curves still draw
  }

  const passed = rows.filter((r) => r.group === 'passed').length;
  const first = rows[0] || {};

  if (!isCurrentRender(token)) return;
  view().replaceChildren(
    el('h1', {text: meta.title || run}),
    el('p', {class: 'sub', text: subtitle}),
    el('p', {class: 'sub', text: t('run.sub', rows.length, passed,
      first.symmetry ? ` · ${first.symmetry}` : '',
      first.n_res_monomer ? t('run.sub.monomer', first.n_res_monomer) : '')}),
    submitPanel(meta, info, false),
    exportBar(run),
    sortableTable(`backbones:${run}`, BACKBONE_COLUMNS, rows, {
      onclick: (r) => {
        location.hash = `#/${encodeURIComponent(run)}/${encodeURIComponent(r.backbone_id)}`;
      },
    }));
}


const rankOf = (row) => {
  const head = String(row.folder || '').split('/').pop().split('_')[0];
  return /^\d+$/.test(head) ? Number(head) : null;
};

const dist = (row, field) => ((row.distribution || {})[field] || []);

/* The counts behind the monomer rate, taken from the same arrays the curves are
 * drawn from, so the fraction in the header and the shape beside it cannot
 * disagree.
 *
 * The denominator is every sequence, not every sequence that has an RMSD. That
 * is what the rate the pipeline stored is over - `success_rate` is n_pass over
 * len(sequences) - and a sequence whose prediction is missing still counts
 * against the backbone. Dividing by the measured ones instead would quietly
 * flatter a backbone with a failed prediction. */
const countPass = (row) => dist(row, 'sc_pass').filter(Boolean).length;
const countAll = (row) => dist(row, 'sc_pass').length;
const cutoff = (key) => (isNil(state.cutoffs[key]) ? null : state.cutoffs[key]);

const BACKBONE_COLUMNS = [
  {key: 'rank', label: 'col.rank', title: 'col.rank.title', left: true, value: rankOf,
   cell: (r) => el('td', {class: 'left num'},
                   isNil(rankOf(r)) ? el('span', {class: 'nil', text: '—'}) : String(rankOf(r)))},
  {key: 'backbone_id', label: 'col.backbone', left: true,
   cell: (r) => el('td', {class: 'id left', text: r.backbone_id})},
  {key: 'final_pass', label: 'col.delivered', title: 'col.delivered.title', desc: true,
   cell: (r) => flagCell(r.final_pass)},
  {key: 'success_rate', label: 'col.success', title: 'col.success.title', desc: true,
   cell: (r) => meterCell(r.success_rate, countPass(r), countAll(r))},
  // The three distributions replace the single best RMSD and best pLDDT, which
  // said almost nothing about twenty sequences.
  {key: 'rmsd_dist', label: 'col.rmsd', title: 'col.rmsd.dist.title', value: false,
   cell: (r) => el('td', {class: 'spark-cell'},
                   sparkDensity('sc_rmsd', dist(r, 'sc_rmsd'), cutoff('rmsd_cutoff')))},
  {key: 'plddt_dist', label: 'col.plddt', title: 'col.plddt.dist.title', value: false,
   cell: (r) => el('td', {class: 'spark-cell'},
                   sparkDensity('mean_plddt', dist(r, 'mean_plddt'),
                                  cutoff('plddt_cutoff')))},
  {key: 'ddg_dist', label: 'col.ddg', title: 'col.ddg.dist.title', value: false,
   cell: (r) => el('td', {class: 'spark-cell'},
                   sparkDensity('dG_interface', dist(r, 'dG_interface'),
                                  cutoff('ddg_cutoff')))},
  {key: 'interface_success_rate', label: 'col.iface.rate', title: 'col.iface.rate.title',
   desc: true,
   cell: (r) => meterCell(r.interface_success_rate, r.n_interface_pass,
                          r.n_interface_scored)},
  {key: 'n_seqs', label: 'col.seqs', title: 'col.seqs.title',
   cell: (r) => numCell(r.n_seqs, 0)},
  {key: 'n_predictions', label: 'col.preds', title: 'col.preds.title',
   cell: (r) => numCell(r.n_predictions, 0)},
  {key: 'final_fail_reason', label: 'col.failreason', left: true,
   cell: (r) => el('td', {class: 'left reason', text: r.final_fail_reason || ''})},
];

/* --------------------------------------------------------------- the export */

const EXPORT_SCOPES = ['delivered', 'passed', 'all'];

/** Two formats, because they are read by different things. FASTA goes to whoever
 *  orders the synthesis and carries its metrics in the description line, since it
 *  travels on its own. The CSV is written with a byte-order mark so Excel opens
 *  it as UTF-8 on a double click. */
function exportBar(run) {
  const base = `/api/runs/${encodeURIComponent(run)}/export`;
  const scope = el('select', {onchange: () => update()},
    EXPORT_SCOPES.map((s) => el('option', {
      value: s, text: t(`export.scope.${s}`), title: t(`export.scope.${s}.help`)})));
  const fasta = el('a', {class: 'btnlink', text: t('export.fasta')});
  const csv = el('a', {class: 'btnlink', text: t('export.csv')});
  const help = el('span', {class: 'hint'});

  function update() {
    fasta.href = `${base}?format=fasta&scope=${scope.value}`;
    csv.href = `${base}?format=csv&scope=${scope.value}`;
    help.textContent = t(`export.scope.${scope.value}.help`);
  }
  update();

  return el('div', {class: 'panel exportbar'},
    el('div', {class: 'controls'}, el('label', {}, t('export.scope'), scope), fasta, csv),
    help,
    el('p', {class: 'hint', text: t('export.note', DG_LABEL)}));
}

/* ------------------------------------------- level 3: one backbone + viewer */

const SEQUENCE_COLUMNS = [
  {key: 'seq_index', label: 'col.index', left: true,
   cell: (r) => el('td', {class: 'left num'}, String(r.seq_index), r.best ? ' ' : null,
                   r.best ? el('span', {class: 'badge best', text: t('detail.best')}) : null)},
  {key: 'sc_pass', label: 'col.pass', title: 'col.pass.title', desc: true,
   cell: (r) => flagCell(r.sc_pass)},
  {key: 'sc_rmsd', label: 'col.rmsd', title: 'col.rmsd.title',
   cell: (r) => numCell(r.sc_rmsd, 2)},
  {key: 'mean_plddt', label: 'col.plddt', title: 'col.plddt.title', desc: true,
   cell: (r) => numCell(r.mean_plddt, 1)},
  // Clicking the energy asks to see where it came from, so this cell does more
  // than select the row: it switches the viewer to the relaxed dimer the number
  // was measured on and opens the interface beneath it. The row click alone
  // leaves the current kind of structure in place, which is right everywhere
  // else and wrong here.
  {key: 'dG_interface', label: 'col.ddg', title: 'col.ddg.title',
   cell(r) {
     const td = numCell(r.dG_interface, 1);
     if (!isNil(r.dG_interface)) {
       td.classList.add('clickable');
       td.title = t('col.ddg.open');
       td.addEventListener('click', (e) => {
         e.stopPropagation();                     // not the row's plain select
         focusInterface(r.sequence_id);
       });
     }
     return td;
   }},
  {key: 'interface_pass', label: 'col.iface', title: 'col.iface.title', desc: true,
   cell: (r) => flagCell(r.interface_pass)},
  {key: 'sequence', label: 'col.sequence', left: true,
   cell: (r) => el('td', {class: 'left seq', title: r.sequence, text: r.sequence})},
];

async function renderBackbone(run, backboneId) {
  const token = renderToken;
  crumbs({text: t('projects.title'), href: '#/'},
         {text: run, href: `#/${encodeURIComponent(run)}`},
         {text: backboneId});

  const detail = await getJSON(`/api/runs/${encodeURIComponent(run)}/backbones/`
                               + encodeURIComponent(backboneId));
  state.detail = detail;
  state.cutoffs = detail.criteria || {};

  if (!detail.sequences.some((s) => s.sequence_id === state.seqId)) {
    state.seqId = detail.best_sequence_id ||
      (detail.sequences.length ? detail.sequences[0].sequence_id : null);
    state.structPath = null;
  }
  if (!detail.structures.some((s) => s.path === state.structPath && !s.missing)) {
    state.structPath = defaultStructure();
  }

  const row = detail.row || {};
  const crit = detail.criteria || {};

  nodes.tableHost = el('div', {});
  nodes.controlHost = el('div', {class: 'controls top'});
  nodes.scaleHost = el('div', {class: 'scale'});
  nodes.noteHost = el('div', {class: 'structnote'});
  nodes.ifaceHost = el('div', {});
  nodes.linkHost = el('div', {class: 'controls bottom'});
  nodes.msg = el('div', {class: 'viewer-msg', hidden: true});
  nodes.stage = STAGE;                                 // the one Mol* instance

  const criteria = isNil(crit.ddg_cutoff)
    ? t('detail.criteria.noddg', fixed(crit.rmsd_cutoff, 1),
        fixed(crit.plddt_cutoff, 0), DG_LABEL)
    : t('detail.criteria', fixed(crit.rmsd_cutoff, 1), fixed(crit.plddt_cutoff, 0),
        DG_LABEL, fixed(crit.ddg_cutoff, 0));

  const left = el('div', {class: 'panel'},
    el('h2', {text: t('detail.sequences')}),
    el('div', {class: 'criteria', text: criteria}),
    el('div', {class: 'facts'},
      fact(t('fact.success'), fixed(detail.success_rate, 2)),
      fact(t('fact.passed'), `${detail.n_pass ?? '—'}`, `/ ${detail.n_sequences ?? '—'}`),
      fact(t('fact.bestrmsd'), fixed(row.best_rmsd, 2), 'Å'),
      fact(t('fact.bestplddt'), fixed(row.best_plddt, 1)),
      fact(t('fact.iface'), fixed(row.interface_success_rate, 2),
           t('fact.ifacecount', row.n_interface_pass ?? '—',
             row.n_interface_scored ?? '—')),
      fact(t('fact.delivered'), detail.group === 'passed' ? t('yes') : t('no'))),
    nodes.tableHost,
    el('p', {class: 'hint', text: t('detail.pickrow')}));

  const right = el('div', {class: 'panel viewer-panel'},
    el('h2', {text: t('detail.structure')}),
    nodes.controlHost,
    nodes.msg,
    nodes.stage,
    nodes.scaleHost,
    nodes.linkHost,
    nodes.noteHost,
    nodes.ifaceHost);

  if (!isCurrentRender(token)) return;
  view().replaceChildren(
    el('h1', {text: backboneId}),
    el('p', {class: 'sub',
             text: `${run} · ${detail.folder}`
                 + (row.symmetry ? ` · ${row.symmetry}` : '')
                 + (row.n_res_monomer ? t('run.sub.monomer', row.n_res_monomer) : '')}),
    el('div', {class: 'detail'}, left, right));

  paintDetail();
  resizeViewer();                                      // the host just moved
  await showStructure();
}

function fact(key, value, unit) {
  return el('div', {class: 'fact'},
    el('div', {class: 'k', text: key}),
    el('div', {class: 'v'}, value, unit ? ' ' : null,
       unit ? el('span', {class: 'u', text: unit}) : null));
}

/** Repaint the parts of the detail view that depend on the selection, leaving
 *  the viewport — and therefore the camera — alone. */
function paintDetail() {
  nodes.tableHost.replaceChildren(
    sortableTable(`sequences:${state.run}:${state.backboneId}`,
                  SEQUENCE_COLUMNS, state.detail.sequences, {
                    selected: (r) => r.sequence_id === state.seqId,
                    onclick: (r) => selectSequence(r.sequence_id),
                    onsort: paintDetail,
                  }));
  nodes.controlHost.replaceChildren(...controls());
  nodes.linkHost.replaceChildren(...links());
  nodes.noteHost.replaceChildren(...structureNote());
  paintScale();
  paintInterface();
}

function selectSequence(sequenceId) {
  if (sequenceId === state.seqId) return;
  state.seqId = sequenceId;
  const mine = state.detail.structures.filter(
    (s) => s.sequence_id === sequenceId && !s.missing);
  const current = currentStructure();
  // Keep looking at the same kind of structure for the newly chosen sequence.
  const same = current ? mine.find((s) => s.key === current.key) : null;
  const pick = same || mine.find((s) => s.key === 'monomer') || mine[0];
  state.structPath = pick ? pick.path : state.structPath;
  paintDetail();
  showStructure();
}

function defaultStructure() {
  const structures = (state.detail.structures || []).filter((s) => !s.missing);
  const mine = structures.filter((s) => s.sequence_id === state.seqId);
  const pick = mine.find((s) => s.key === 'monomer') || mine[0] || structures[0];
  return pick ? pick.path : null;
}

/* ------------------------------------------------------------- the structure
 *
 * Mol*, the viewer RCSB serves, with its own panels switched off: the structure
 * list here knows about files Mol* cannot, and the rest of its UI was more than
 * this page needs. What is left is the viewport, its expand button, and two
 * selects.
 */

const MOLSTAR_VIEWER_OPTIONS = {
  layoutIsExpanded: false,
  layoutShowControls: false,          // no left or right panel
  layoutShowSequence: false,
  layoutShowLeftPanel: false,
  layoutShowRemoteState: false,
  layoutShowLog: false,
  viewportShowExpand: true,           // the one control worth keeping
  viewportShowSelectionMode: false,
  viewportShowAnimation: false,
};

/** The viewport background, read from the stylesheet so the theme owns it.
 *
 *  One source of truth: `--stage` is defined once per theme in style.css, the
 *  element behind the canvas uses it, and Mol* is told the same value. Nothing
 *  here decides what light or dark looks like.
 */
function stageColor() {
  const value = getComputedStyle(document.documentElement)
    .getPropertyValue('--stage').trim();
  return /^#[0-9a-f]{6}$/i.test(value) ? value : '#fcfcfb';
}

/** Repaint the live viewport after the theme changes.
 *
 *  `viewportBackgroundColor` is read only while the plugin is being built, and
 *  the plugin is built once per page load, so a toggle has to reach the canvas
 *  directly. This is the same call Mol* makes for that option itself, and it
 *  takes a number where the option takes a string.
 */
function applyStageBackground() {
  if (!state.viewer) return;
  try {
    const colour = parseInt(stageColor().slice(1), 16);
    state.viewer.plugin.canvas3d.setProps({renderer: {backgroundColor: colour}});
  } catch (e) {
    /* an older build, or no canvas yet; the next page load picks it up */
  }
}

/* Mol* takes over the element it is given, so the host holds nothing else and is
 * built once: creating the plugin costs a five-megabyte bundle, and losing it
 * would lose the camera, so it has to survive a sort, a row click and a
 * navigation. It is moved into each detail view rather than rebuilt. */
const STAGE = el('div', {class: 'stage', id: 'molstar-host'});

const COLOR_THEMES = [
  {key: 'chain', theme: 'chain-id', label: 'colour.chain'},
  {key: 'plddt', theme: 'plddt-confidence', label: 'colour.plddt'},
  {key: 'ss', theme: 'secondary-structure', label: 'colour.ss'},
];

/* The four thresholds Mol* writes into its own pLDDT theme, read out of the
 * bundle rather than assumed, so this legend is the colours the viewport is
 * actually showing. They are AlphaFold's own bands. */
const PLDDT_BANDS = [
  {color: '#0053D6', label: 'legend.plddt.vhigh'},
  {color: '#65CBF3', label: 'legend.plddt.high'},
  {color: '#FFDB13', label: 'legend.plddt.low'},
  {color: '#FF7D45', label: 'legend.plddt.vlow'},
];

const interfaceURL = (run, path) =>
  `/api/runs/${encodeURIComponent(run)}/interface`
  + `?path=${encodeURIComponent(path)}`;

let viewerStarting = null;
let viewerError = null;

const currentStructure = () =>
  (state.detail.structures || []).find((s) => s.path === state.structPath) || null;

const structLabel = (s) => (s.seq_short ? `${s.seq_short} · ` : '') + t('struct.' + s.key);

/** Whether a colouring can say anything about this structure.
 *
 *  The design backbone carries B-factor 0.00 on every atom, so Mol*'s pLDDT
 *  theme would paint the whole thing the "very low confidence" orange. That is
 *  not a low-confidence model, it is a file with no confidence in it, so the
 *  option is withheld rather than offered as a wrong answer.
 */
function themeApplies(key) {
  const current = currentStructure();
  if (key === 'plddt') return !!current && current.key !== 'backbone';
  return true;
}

function controls() {
  const structures = state.detail.structures || [];
  const current = currentStructure();

  const select = el('select', {
    onchange: (e) => {
      state.structPath = e.target.value;
      // Picking something out of "other sequences" moves the table's selection
      // with it, so the metrics on the left always belong to what is on screen.
      // Without this the two halves of the page could describe different
      // sequences and nothing would say so.
      const chosen = (state.detail.structures || [])
        .find((x) => x.path === e.target.value);
      if (chosen && chosen.sequence_id) state.seqId = chosen.sequence_id;
      paintDetail();
      showStructure();
    },
  });
  const groups = [
    [t('viewer.group.backbone'), structures.filter((s) => s.key === 'backbone')],
    [t('viewer.group.current', short(state.seqId)),
     structures.filter((s) => s.sequence_id === state.seqId)],
    [t('viewer.group.other'),
     structures.filter((s) => s.sequence_id && s.sequence_id !== state.seqId)],
  ];
  for (const [label, items] of groups) {
    if (!items.length) continue;
    const group = el('optgroup', {label});
    for (const s of items) {
      group.append(el('option', {
        value: s.path, selected: s.path === state.structPath, disabled: s.missing,
        text: structLabel(s) + (s.missing ? t('viewer.missing') : ''),
      }));
    }
    select.append(group);
  }

  const colour = el('select', {
    onchange: (e) => { state.colorMode = e.target.value; paintDetail(); showStructure(); },
  }, COLOR_THEMES.map((th) => el('option', {
    value: th.key, selected: state.colorMode === th.key, disabled: !themeApplies(th.key),
    text: t(th.label) + (themeApplies(th.key) ? '' : t('colour.na')),
  })));

  return [
    el('label', {class: 'grow'}, t('viewer.structure'), select),
    el('label', {}, t('viewer.colour'), colour),
  ];
}

function links() {
  const current = currentStructure();
  const seq = (state.detail.sequences || []).find((s) => s.sequence_id === state.seqId);
  return [
    current ? el('a', {href: fileURL(state.run, current.path, true),
                       text: t('viewer.download', current.path.split('/').pop())}) : null,
    seq ? el('a', {href: '#', onclick(e) {
      e.preventDefault();
      navigator.clipboard.writeText(seq.sequence).then(
        () => { e.target.textContent = t('viewer.copied'); },
        () => { e.target.textContent = t('viewer.copyfail'); });
    }, text: t('viewer.copyseq', short(seq.sequence_id))}) : null,
    el('a', {href: '#', class: 'quiet', onclick(e) {
      e.preventDefault();
      runSelfCheck();
    }, text: t('viewer.diagnose')}),
  ].filter(Boolean);
}

/** The scale legend, for the one colouring whose colours are known here.
 *
 *  Mol*'s chain and secondary-structure themes draw from its own palettes, which
 *  this page does not control, so it names the theme instead of drawing swatches
 *  that might not match. A legend that disagrees with the viewport is worse than
 *  no legend.
 */
function paintScale() {
  const items = [];
  if (state.colorMode === 'plddt') {
    items.push(el('span', {class: 'item', text: t('legend.plddt')}));
    for (const band of PLDDT_BANDS) {
      items.push(el('span', {class: 'item'},
        el('span', {class: 'sw', style: `background:${band.color}`}), t(band.label)));
    }
    items.push(el('span', {class: 'item', text: t('legend.plddt.source')}));
  } else {
    const chosen = COLOR_THEMES.find((th) => th.key === state.colorMode);
    items.push(el('span', {class: 'item',
      text: t('legend.molstar', chosen ? chosen.theme : 'chain-id')}));
  }
  nodes.scaleHost.replaceChildren(...items);
}

function structureNote() {
  const current = currentStructure();
  if (!current) return [];
  return [document.createTextNode(t('note.' + current.key, DG_LABEL))];
}

function message(text) {
  if (!nodes.msg) return;
  nodes.msg.textContent = text || '';
  nodes.msg.hidden = !text;
}

/* ---------------------------------------------------- the interface residues
 *
 * Stage 05 records how good the interface is and not which residues make it, so
 * they are computed on request from the structure the energy was measured on.
 * Only the relaxed dimer: an AlphaFold monomer has no second chain, and the
 * design backbone is poly-glycine with only its main chain, where the closest
 * atoms across two subunits sit 5.42 angstrom apart - just outside the 5
 * angstrom contact cutoff, so it would honestly report nothing and look broken.
 */

let ifaceToken = 0;

function switchSequence(sequenceId) {
  if (!sequenceId || sequenceId === state.seqId) return false;
  state.seqId = sequenceId;
  return true;
}

/** Show the relaxed dimer of one sequence, and its interface. */
function focusInterface(sequenceId) {
  const dimer = (state.detail.structures || []).find(
    (x) => x.key === 'dimer' && x.sequence_id === sequenceId && !x.missing);
  switchSequence(sequenceId);
  if (dimer) state.structPath = dimer.path;
  paintDetail();
  showStructure();
}

function paintInterface() {
  const host = nodes.ifaceHost;
  if (!host) return;
  const current = currentStructure();
  if (!current || current.key !== 'dimer') {
    host.replaceChildren();
    return;
  }

  const token = ++ifaceToken;
  host.replaceChildren(el('div', {class: 'hint', text: t('iface.loading')}));

  getJSON(interfaceURL(state.run, current.path))
    .then((data) => {
      if (token !== ifaceToken) return;            // a faster click won
      const chains = Object.entries(data.chains || {});
      if (!chains.length) {
        host.replaceChildren(el('div', {class: 'hint', text: t('iface.none')}));
        return;
      }
      host.replaceChildren(
        el('div', {class: 'ifacebox'},
          el('div', {class: 'ifacehead'},
            el('strong', {text: t('iface.title')}),
            el('span', {class: 'hint',
                        text: t('iface.sub', data.n_residues, data.cutoff)})),
          ...chains.map(([chain, residues]) => el('div', {class: 'ifacerow'},
            el('span', {class: 'ifacechain', text: t('iface.chain', chain)}),
            el('span', {class: 'ifaceres', text: residues.join(' ')})))));
    })
    .catch((err) => {
      if (token !== ifaceToken) return;
      host.replaceChildren(el('div', {class: 'hint',
                                      text: t('iface.fail', err.message)}));
    });
}

/* ------------------------------------------------------------- the self-check
 *
 * This page cannot be opened from where it is developed, so when the viewer does
 * not draw there is no console to read. This gathers what would have been in it
 * and posts it to the server log, which can be read without asking anyone to
 * copy text out of a browser.
 */

function webglReport() {
  try {
    const canvas = document.createElement('canvas');
    for (const kind of ['webgl2', 'webgl', 'experimental-webgl']) {
      const gl = canvas.getContext(kind);
      if (!gl) continue;
      const dbg = gl.getExtension('WEBGL_debug_renderer_info');
      return {
        context: kind,
        vendor: dbg ? gl.getParameter(dbg.UNMASKED_VENDOR_WEBGL) : null,
        renderer: dbg ? gl.getParameter(dbg.UNMASKED_RENDERER_WEBGL) : null,
        version: gl.getParameter(gl.VERSION),
      };
    }
    return {context: null, note: 'no WebGL context of any kind'};
  } catch (err) {
    return {context: null, error: String((err && err.message) || err)};
  }
}

function selfCheck() {
  const ms = window.molstar;
  const host = nodes.stage;
  return {
    at: new Date().toISOString(),
    url: location.href,
    route: location.hash || '#/',
    molstar: typeof ms,
    molstar_viewer: ms ? typeof ms.Viewer : null,
    molstar_create: (ms && ms.Viewer) ? typeof ms.Viewer.create : null,
    molstar_keys: ms ? Object.keys(ms).slice(0, 24) : null,
    viewer_built: !!state.viewer,
    stage: host ? {
      in_document: document.body.contains(host),
      width: host.clientWidth,
      height: host.clientHeight,
      children: host.childElementCount,
    } : null,
    webgl: webglReport(),
    structure: state.structPath,
    error: viewerError ? {
      name: viewerError.name || null,
      message: String(viewerError.message || viewerError),
      stack: String(viewerError.stack || '').split('\n').slice(0, 6).join('\n'),
    } : null,
    user_agent: navigator.userAgent,
  };
}

async function report(reason) {
  const data = selfCheck();
  try {
    await postJSON('/api/clientlog', {reason, report: data});
  } catch (err) { /* the page is more useful than the report */ }
  return data;
}

async function runSelfCheck() {
  const data = await report('manual');
  nodes.noteHost.replaceChildren(
    el('details', {class: 'cmdbox', open: true},
      el('summary', {text: t('viewer.diagnose')}),
      el('pre', {class: 'cmd', text: JSON.stringify(data, null, 2)}),
      el('p', {class: 'hint', text: t('viewer.diagnose.hint')})));
}

/* ---------------------------------------------------------------- the viewer */

/** The Mol* plugin, built at most once per page load. */
async function ensureViewer() {
  if (state.viewer) return state.viewer;
  if (!viewerStarting) {
    // A CSS colour string, not a number. Mol* normalises this option through a
    // parser that calls .toLowerCase() on it, so passing 0x1a1a19 threw a
    // TypeError out of Viewer.create and nothing was ever drawn - which is what
    // the viewer's first version did, and why its panel stayed empty.
    viewerStarting = window.molstar.Viewer.create(
      STAGE, {...MOLSTAR_VIEWER_OPTIONS, viewportBackgroundColor: stageColor()});
  }
  state.viewer = await viewerStarting;
  return state.viewer;
}

/** Tell Mol* its element changed size.
 *
 *  Needed after the host is moved into a freshly built detail view, and after
 *  the window changes. `handleResize` is the direct call; the resize event is
 *  the fallback, because Mol* listens for it itself.
 */
function resizeViewer() {
  if (!state.viewer) return;
  try {
    state.viewer.handleResize();
  } catch (e) {
    window.dispatchEvent(new Event('resize'));
  }
}

async function showStructure() {
  const current = currentStructure();
  if (!nodes.stage) return;
  if (!current) { message(t('viewer.nostructure')); return; }
  if (!window.molstar) {
    viewerError = new Error('window.molstar is undefined');
    message(t('viewer.noengine'));
    report('no-molstar');
    return;
  }

  message(t('viewer.loading', structLabel(current)));
  let viewer;
  try {
    viewer = await ensureViewer();
  } catch (err) {
    viewerError = err;
    viewerStarting = null;                             // let a retry try again
    message(t('viewer.initfail', (err && err.message) || err));
    report('create-failed');
    return;
  }
  if (currentStructure() !== current) return;          // a faster click won

  const url = new URL(fileURL(state.run, current.path), location.href).href;
  const label = current.path.split('/').pop();
  const theme = (COLOR_THEMES.find((th) => th.key === state.colorMode)
                 || COLOR_THEMES[0]).theme;

  try { await viewer.plugin.clear(); } catch (e) { /* nothing loaded yet */ }

  try {
    await viewer.loadStructureFromUrl(url, 'pdb', false, {
      representationParams: {theme: {globalName: theme}}, label,
    });
  } catch (themeError) {
    // A theme that does not apply must not cost us the structure.
    try {
      await viewer.loadStructureFromUrl(url, 'pdb', false, {label});
      message(t('viewer.themefail'));
      paintScale();
      resizeViewer();
      return;
    } catch (err) {
      viewerError = err;
      message(t('viewer.loadfail', (err && err.message) || err));
      report('load-failed');
      return;
    }
  }

  viewerError = null;
  paintScale();
  resizeViewer();
  message('');
}

/* ------------------------------------------------------------- level 0: jobs */

let logTimer = null;
let logOffset = 0;

function stopFollowing() {
  if (logTimer) { clearInterval(logTimer); logTimer = null; }
}

const stamp = (x) => (x ? new Date(x * 1000).toLocaleString() : '—');

function duration(job) {
  const end = job.finished || Math.floor(Date.now() / 1000);
  const secs = Math.max(0, end - (job.started || end));
  if (secs < 90) return t('time.seconds', secs);
  if (secs < 5400) return t('time.minutes', (secs / 60).toFixed(0));
  return t('time.hours', (secs / 3600).toFixed(1));
}

const JOB_BADGE = {running: 'running', finished: 'done', failed: 'bad',
                   cancelled: '', unknown: 'warn'};

let formSpec = null;
const formState = {};
let warnHost = null;

function field(spec, fixed) {
  const id = `f_${spec.name}`;
  if (fixed && Object.prototype.hasOwnProperty.call(fixed, spec.name)) {
    // Decided by the project, so it is reported rather than offered. With one
    // project per run directory, letting the form change the run name would let
    // a job submitted from inside a project write somewhere else.
    formState[spec.name] = fixed[spec.name];
    return el('div', {class: 'formrow'},
      el('label', {text: t(spec.label)}),
      el('div', {class: 'fixedvalue', text: String(fixed[spec.name])}),
      el('div', {class: 'fieldnote'}, el('code', {text: spec.key})));
  }
  let input;
  if (spec.kind === 'choice') {
    input = el('select', {id}, spec.choices.map((c) => el('option', {
      value: c, selected: String(formState[spec.name]) === String(c), text: c})));
  } else if (spec.kind === 'name') {
    input = el('input', {id, type: 'text', value: formState[spec.name] ?? '',
                         autocomplete: 'off'});
  } else {
    input = el('input', {
      id, type: 'number', value: formState[spec.name] ?? spec.default,
      step: spec.kind === 'float' ? '0.05' : '1',
      min: isNil(spec.minimum) ? null : spec.minimum,
      max: isNil(spec.maximum) ? null : spec.maximum,
    });
  }
  const take = () => { formState[spec.name] = input.value; paintWarn(); };
  input.addEventListener('input', take);
  input.addEventListener('change', take);

  const range = (!isNil(spec.minimum) && !isNil(spec.maximum))
    ? ` · ${spec.minimum} – ${spec.maximum}` : '';
  // `spec.label` is a key: the server sends keys so the form can be rendered in
  // either language. Help is optional, so its key is asked for rather than
  // looked up blind, which would print `jobs.field.seed.help` at the reader.
  const help = `${spec.label}.help`;
  return el('div', {class: 'formrow'},
    el('label', {for: id, text: t(spec.label)}),
    input,
    el('div', {class: 'fieldnote'},
      el('code', {text: spec.key}), range,
      has(help) ? el('div', {text: t(help)}) : null));
}

/** The one number worth shouting about before a click costs hours. */
function paintWarn() {
  if (!warnHost || !formSpec) return;
  const notes = [];
  for (const spec of formSpec.fields) {
    if (isNil(spec.warn_above)) continue;
    const value = Number(formState[spec.name]);
    if (Number.isFinite(value) && value > spec.warn_above) {
      notes.push(t('jobs.warn.above', spec.label, value, spec.warn_above));
    }
  }
  if (formState.stages === 'all' && Number(formState.n_backbones) >= 200) {
    notes.push(t('jobs.warn.full'));
  }
  warnHost.replaceChildren(...notes.map((text) => el('div', {class: 'warnline', text})));
}

/** The parameter form. `fixed` names the fields the caller has already decided,
 *  which is how the same form serves the global page and a project's own page. */
function jobFormPanel(info, fixed) {
  const problem = el('div', {});
  warnHost = el('div', {class: 'warnbox'});

  if (!Object.keys(formState).length) {
    for (const spec of formSpec.fields) formState[spec.name] = spec.default;
    formState.preset = formSpec.presets[0];
    formState.stages = 'all';
  }

  const preset = el('select', {id: 'f_preset',
    onchange: (e) => { formState.preset = e.target.value; }},
    formSpec.presets.map((p) => el('option', {
      value: p, selected: formState.preset === p,
      text: p === 'none' ? t('jobs.preset.none') : p})));

  const stages = el('select', {id: 'f_stages',
    onchange: (e) => { formState.stages = e.target.value; paintWarn(); }},
    Object.keys(formSpec.stage_sets).map((k) => el('option', {
      value: k, selected: formState.stages === k, text: t(`jobs.stages.${k}`)})));

  const submit = el('button', {class: 'primary', type: 'submit', text: t('jobs.start')});

  const form = el('form', {class: 'jobform', onsubmit: async (e) => {
    e.preventDefault();
    submit.disabled = true;
    problem.replaceChildren();
    try {
      const job = await postJSON('/api/jobs', formState);
      location.hash = `#/_jobs/${encodeURIComponent(job.id)}`;
    } catch (err) {
      problem.replaceChildren(el('div', {class: 'error', text: err.message}));
      submit.disabled = false;
    }
  }},
    el('div', {class: 'formgrid'},
      formSpec.fields.map((spec) => field(spec, fixed)),
      el('div', {class: 'formrow'},
        el('label', {for: 'f_preset', text: t('jobs.preset')}), preset,
        el('div', {class: 'fieldnote', text: t('jobs.preset.help')})),
      el('div', {class: 'formrow'},
        el('label', {for: 'f_stages', text: t('jobs.stages')}), stages,
        el('div', {class: 'fieldnote', text: t('jobs.stages.help')}))),
    warnHost,
    problem,
    el('div', {class: 'controls'}, submit,
       el('span', {class: 'hint', text: info.running ? t('jobs.busy') : t('jobs.free')})));

  paintWarn();
  return form;
}

function submitDisabledPanel() {
  return el('div', {class: 'panel notice'},
    el('p', {text: t('jobs.disabled')}),
    el('p', {class: 'hint', text: t('jobs.disabled.why')}),
    el('pre', {class: 'cmd',
               text: 'python -u -m cyclicnano.webui --port 8080 \\\n'
                   + '    --allow-submit --profile configs/profiles/amax.yaml'}));
}

async function renderJobs() {
  const token = renderToken;
  crumbs({text: t('projects.title'), href: '#/'}, {text: t('nav.submit')});
  stopFollowing();

  const info = await getJSON('/api/jobs');
  formSpec = info.spec;

  if (!info.enabled) {
    if (!isCurrentRender(token)) return;
  view().replaceChildren(el('h1', {text: t('jobs.title')}), submitDisabledPanel());
    return;
  }

  const parts = [el('h1', {text: t('jobs.title')}),
                 el('p', {class: 'sub', text: t('jobs.sub')})];
  if (info.running) parts.push(jobCard(info.running, true));
  parts.push(el('div', {class: 'panel'},
    el('h2', {text: t('jobs.params')}),
    jobFormPanel(info, null)));
  if (info.jobs.length) {
    parts.push(el('div', {class: 'panel'},
      el('h2', {text: t('jobs.history')}),
      el('div', {class: 'joblist'}, info.jobs.map((j) => jobRow(j)))));
  }
  if (!isCurrentRender(token)) return;
  view().replaceChildren(...parts);
  paintWarn();
}

function stateBadge(job) {
  return el('span', {class: `badge ${JOB_BADGE[job.state] || ''}`.trim(),
                     text: t(`job.state.${job.state}`)});
}

function jobRow(job) {
  return el('a', {class: 'jobrow', href: `#/_jobs/${encodeURIComponent(job.id)}`},
    stateBadge(job),
    el('span', {class: 'jid', text: job.run}),
    el('span', {class: 'hint', text: `${job.stages} · ${job.preset}`}),
    el('span', {class: 'spacer'}),
    el('span', {class: 'hint', text: `${stamp(job.started)} · ${duration(job)}`}));
}

function jobCard(job, brief) {
  const bits = [
    el('div', {class: 'card-head'},
      el('span', {class: 'name', text: job.run}),
      stateBadge(job),
      job.resumed ? el('span', {class: 'badge warn', text: t('job.resumed')}) : null),
    el('div', {class: 'facts'},
      fact(t('job.id'), job.id),
      fact(t('job.started'), stamp(job.started)),
      fact(t('job.elapsed'), duration(job)),
      fact(t('job.stages'), String(job.stages)),
      isNil(job.returncode) ? null : fact(t('job.returncode'), String(job.returncode))),
  ];
  if (job.note) bits.push(el('div', {class: 'warnline', text: job.note}));
  if (brief) {
    bits.push(el('p', {class: 'hint'},
      el('a', {href: `#/_jobs/${encodeURIComponent(job.id)}`, text: t('job.showlog')})));
  }
  return el('div', {class: 'panel'}, ...bits.filter(Boolean));
}

async function renderJob(jobId) {
  const token = renderToken;
  crumbs({text: t('projects.title'), href: '#/'},
         {text: t('nav.submit'), href: '#/_jobs'},
         {text: jobId});
  stopFollowing();
  logOffset = 0;

  const job = await getJSON(`/api/jobs/${encodeURIComponent(jobId)}`);
  const logBox = el('pre', {class: 'log', text: ''});
  const head = el('div', {});
  const problem = el('div', {});

  const cancel = el('button', {text: t('job.cancel'), onclick: async () => {
    cancel.disabled = true;
    try {
      await postJSON(`/api/jobs/${encodeURIComponent(jobId)}/cancel`, {});
      await refresh();
    } catch (err) {
      problem.replaceChildren(el('div', {class: 'error', text: err.message}));
      cancel.disabled = false;
    }
  }});

  function paintHead(current) {
    head.replaceChildren(jobCard(current, false));
    cancel.hidden = current.state !== 'running';
  }

  async function pullLog() {
    try {
      const chunk = await getJSON(
        `/api/jobs/${encodeURIComponent(jobId)}/log?offset=${logOffset}`);
      if (chunk.text) {
        // Stick to the bottom only if the reader is already there, so scrolling
        // back through a long log is not yanked away on the next poll.
        const atBottom = logBox.scrollHeight - logBox.scrollTop - logBox.clientHeight < 40;
        logBox.append(document.createTextNode(chunk.text));
        logOffset = chunk.offset;
        if (atBottom) logBox.scrollTop = logBox.scrollHeight;
      }
    } catch (err) { /* a poll that fails is retried by the next tick */ }
  }

  async function refresh() {
    const current = await getJSON(`/api/jobs/${encodeURIComponent(jobId)}`);
    paintHead(current);
    await pullLog();
    if (current.state !== 'running') stopFollowing();
    return current;
  }

  if (!isCurrentRender(token)) return;
  view().replaceChildren(
    el('h1', {text: t('job.title', job.run)}),
    head,
    el('div', {class: 'panel'},
      el('h2', {text: t('job.log')}),
      el('div', {class: 'controls'}, cancel,
         el('span', {class: 'hint', text: t('job.follow')})),
      problem,
      logBox,
      el('details', {class: 'cmdbox'},
        el('summary', {text: t('job.command')}),
        el('pre', {class: 'cmd', text: job.command || ''}))));

  paintHead(job);
  await pullLog();
  if (job.state === 'running') logTimer = setInterval(refresh, 3000);
}

/* ------------------------------------------------------------------ routing */

function parseHash() {
  const parts = (location.hash || '#/').replace(/^#\/?/, '').split('/')
    .filter((p) => p.length).map(decodeURIComponent);
  return {run: parts[0] || null, backboneId: parts[1] || null};
}

/* Renders are async and a reader can navigate faster than one finishes, so the
 * slower of two overlapping renders would otherwise paint over the newer page.
 * Each render takes a number and the view-building helpers check it. */
let renderToken = 0;
const isCurrentRender = (token) => token === renderToken;

async function render() {
  const token = ++renderToken;
  const {run, backboneId} = parseHash();
  state.run = run;
  state.backboneId = backboneId;
  try {
    if (run === '_jobs') {
      state.detail = null;
      if (backboneId) await renderJob(backboneId);
      else await renderJobs();
    } else if (run === '_new') { state.detail = null; await renderNewProject(); }
    else if (!run) { state.detail = null; await renderProjects(); }
    else if (!backboneId) { state.detail = null; await renderRun(run); }
    else await renderBackbone(run, backboneId);
  } catch (err) {
    stopFollowing();
    fail(err);
  }
}

window.addEventListener('hashchange', () => {
  const {run, backboneId} = parseHash();
  if (run !== state.run || backboneId !== state.backboneId) {
    state.seqId = null;                                // a different backbone
    state.structPath = null;
  }
  stopFollowing();                                     // leaving a job page
  render();
});

// Mol* watches the window itself; this only covers a layout change that does
// not raise one, such as the detail grid switching to a single column.
window.addEventListener('resize', resizeViewer);

/* ------------------------------------------------------- theme and language */

function isDark() {
  const root = document.documentElement;
  return root.dataset.theme
    ? root.dataset.theme === 'dark'
    : window.matchMedia('(prefers-color-scheme: dark)').matches;
}

function paintThemeButton() {
  const button = document.getElementById('theme');
  // The icon shows the mode a click would move to, which is the convention a
  // reader already knows from every other toggle of this kind.
  button.replaceChildren(icon(isDark() ? SUN : MOON));
  button.title = t('theme.toggle');
}

function paintLangButton() {
  const button = document.getElementById('lang');
  button.textContent = t('lang.name');
  button.title = t('lang.toggle');
}

document.getElementById('theme').addEventListener('click', () => {
  document.documentElement.dataset.theme = isDark() ? 'light' : 'dark';
  try {
    localStorage.setItem('cyclicnano.theme', document.documentElement.dataset.theme);
  } catch (e) { /* blocked */ }
  paintThemeButton();
  applyStageBackground();
});

document.getElementById('lang').addEventListener('click', () => {
  setLang(otherLang());
  paintLangButton();
  paintThemeButton();
  render();                      // the viewer survives; only the text is rebuilt
});

window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', () => {
  if (document.documentElement.dataset.theme) return;   // an explicit choice wins
  paintThemeButton();
  applyStageBackground();
});

/* Anything that escapes a handler. `fail` covers the render chain, which is
 * where errors usually land, but an exception inside a click listener or a
 * promise nobody awaited would otherwise be invisible from here. */
window.addEventListener('error', (e) => {
  viewerError = e.error instanceof Error ? e.error : new Error(String(e.message));
  report('window-error');
});
window.addEventListener('unhandledrejection', (e) => {
  viewerError = e.reason instanceof Error ? e.reason : new Error(String(e.reason));
  report('unhandled-rejection');
});

getJSON('/api/health').then((h) => {
  document.getElementById('rootnote').textContent =
    h.root + (h.ok ? '' : t('root.unreadable'));
}).catch(() => {});

paintThemeButton();
paintLangButton();
render();
