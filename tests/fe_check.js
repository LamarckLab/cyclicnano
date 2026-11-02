/* Render the page against a fake DOM and assert on what comes out.
 *
 * NOT part of the Python suite. `pytest` and `tests/run_all.py` both collect
 * only `test_*.py`, so this file is invisible to them; it is run by hand.
 *
 *     npm install jsdom                       # once, anywhere
 *     python -m cyclicnano.webui --port 8080 &
 *     node tests/fe_check.js
 *
 * Why it exists. The page is developed on a machine that cannot open it, so a
 * JavaScript error that only ever appears in a div on screen is an error nobody
 * here learns about. One such error - `controls()` reading a `current` that an
 * earlier rewrite had deleted - took the whole structure page down, and was
 * found by a reader rather than by anything in this repository. `node --check`
 * cannot see it, because the file parses; a reference check over call sites
 * cannot see it, because `current` is read as a variable and never called.
 * Running the page is what sees it.
 *
 * Mol* is deliberately left out: it needs WebGL, which jsdom has none of. The
 * viewer then reports that it did not load and the rest of the page renders,
 * which is both the correct behaviour and the half worth checking here.
 *
 * Two things this file learned the hard way, so they are not undone.
 *
 * The route goes in the document's URL rather than being assigned to
 * `location.hash` after the scripts load. Assigning it queues a hashchange that
 * jsdom delivers once app.js has installed its listener, so the page renders
 * twice, which a browser opening the same URL never does.
 *
 * Nothing waits on a fixed timeout. Every wait polls for the condition it
 * actually needs. The first version of this file passed when run alone and
 * failed when run after the other routes, purely from being given less of the
 * event loop, and reported a bug in the page that was a bug in itself.
 */
'use strict';
const fs = require('fs');
const path = require('path');

const ORIGIN = process.env.CYCLICNANO_URL || 'http://127.0.0.1:8080';
const STATIC = path.join(__dirname, '..', 'src', 'cyclicnano', 'webui', 'static');

let JSDOM;
try {
  ({JSDOM} = require('jsdom'));
} catch (err) {
  // Say which failure it was. "Not installed" was printed for an installed
  // jsdom too: version 27 is ESM-only and `require` of it throws
  // ERR_REQUIRE_ESM, which looked like a missing package for long enough to
  // waste a deploy. `npm install jsdom@26` is the version this file wants.
  const missing = err && err.code === 'MODULE_NOT_FOUND';
  console.error(missing
    ? 'jsdom is not installed. npm install jsdom@26, then run this again.'
    : `jsdom is installed but cannot be required (${err && err.code}): `
      + `${(err && err.message) || err}
`
      + 'This file is CommonJS, so it needs a CommonJS jsdom: npm install jsdom@26.');
  process.exit(2);
}

const failures = [];
const check = (ok, what) => {
  console.log(`    ${ok ? 'ok  ' : 'FAIL'}  ${what}`);
  if (!ok) failures.push(what);
};

const clean = (p) => ((p && p.name ? `${p.name}: ` : '') + String((p && p.message) || p));

/** Poll until `fn()` is truthy, or give up and record it. */
async function waitFor(fn, what, timeout = 10000, step = 50) {
  const until = Date.now() + timeout;
  for (;;) {
    const value = fn();
    if (value) return value;
    if (Date.now() > until) {
      check(false, `timed out waiting for ${what}`);
      return value;
    }
    await new Promise((r) => setTimeout(r, step));
  }
}

/** A page opened at one route, rendered, and handed back for inspection. */
async function page(route) {
  const html = fs.readFileSync(path.join(STATIC, 'index.html'), 'utf8')
    // The vendored bundle and the two page scripts are injected below instead,
    // so their order is controlled and Mol* is left out.
    .replace(/<script\b[\s\S]*?<\/script>/g, '');

  // `dangerously`, because the files have to run as real <script> elements: a
  // classic script's top-level const and function declarations go into the
  // global lexical environment where the next script can see them, while the
  // same code through eval in strict mode keeps them to itself.
  const dom = new JSDOM(html, {url: `${ORIGIN}/${route}`, runScripts: 'dangerously',
                               pretendToBeVisual: true});
  const {window} = dom;

  // Two things jsdom lacks that every browser has.
  window.matchMedia = (query) => ({
    matches: false, media: query,
    addEventListener() {}, removeEventListener() {},
    addListener() {}, removeListener() {},
  });
  window.fetch = (url, opts) => globalThis.fetch(new URL(String(url), ORIGIN).href, opts);

  const problems = [];
  window.addEventListener('error', (e) => problems.push(e.error || e.message));
  window.addEventListener('unhandledrejection', (e) => problems.push(e.reason));

  for (const name of ['i18n.js', 'app.js']) {
    const tag = window.document.createElement('script');
    tag.textContent = fs.readFileSync(path.join(STATIC, name), 'utf8');
    window.document.body.appendChild(tag);
  }

  const view = window.document.getElementById('view');
  // Rendered means a heading or an error is on the page, not that some number
  // of milliseconds has passed.
  await waitFor(() => view.querySelector('h1, .error'), `${route} to render`);

  return {
    window, view, problems,
    error: () => view.querySelector('.error'),
    headings: () => [...view.querySelectorAll('h1, h2')].map((h) => h.textContent),
    rows: () => view.querySelectorAll('tbody tr').length,
    renders: () => window.eval('renderToken'),
    state: (expr) => window.eval(`state.${expr}`),
  };
}

/* -------------------------------------------------------------------- routes */

const ROUTES = [
  {route: '#/', expect: 'Projects'},
  {route: '#/c5', expect: 'c5', rows: true},
  {route: '#/c5/design_27', expect: 'design_27', rows: true},
  {route: '#/_jobs', expect: 'Submit a run'},
  {route: '#/_new', expect: 'New project'},
];

async function everyRouteRenders() {
  console.log('\n  every route renders, and none of them reports an error');
  for (const {route, expect, rows} of ROUTES) {
    const out = await page(route);
    const error = out.error();
    const detail = error ? ` (.error: ${error.textContent})`
      : (out.problems.length ? ` (${out.problems.map(clean).join('; ')})` : '');
    check(!error && !out.problems.length, `${route} renders cleanly${detail}`);
    check(out.headings().some((h) => h.includes(expect)),
          `${route} is headed ${JSON.stringify(expect)}, `
          + `got ${JSON.stringify(out.headings())}`);
    // Opening a URL paints once. More than that means something is routing
    // twice, and the second paint would overwrite work done between them.
    check(out.renders() === 1,
          `${route} rendered ${out.renders()} time(s), expected 1`);
    if (rows) check(out.rows() > 0, `${route} has table rows, got ${out.rows()}`);
  }
}

/* ----------------------------------------------------- the structure page */

async function theStructurePageWorks() {
  console.log('\n  the structure page: its columns, and what a row and a ddG cell show');
  const out = await page('#/c5/design_27');
  const {view, window} = out;
  const chosen = () => {
    const select = view.querySelector('.controls.top select');
    return select && select.selectedIndex >= 0
      ? select.options[select.selectedIndex].textContent : '';
  };

  const headers = [...view.querySelectorAll('thead th')].map((h) => h.textContent.trim());
  check(headers.join('|') === '#|Pass|RMSD|pLDDT|ddG|Interface|Sequence',
        `sequence columns are named plainly, got ${JSON.stringify(headers)}`);
  check(view.querySelectorAll('tbody tr.selected').length === 1,
        'exactly one sequence is selected');

  // What a sequence is for. The monomer is a click away in the same select; it
  // is not what the page opens on.
  check(/pentamer|trimer|-mer/.test(chosen()),
        `the page opens on the assembly, got ${JSON.stringify(chosen())}`);

  const ddg = [...view.querySelectorAll('td.clickable')];
  check(ddg.length > 0, `the measured ddG cells are clickable, got ${ddg.length}`);

  // Clicking an energy asks to see where it came from.
  const before = out.state('structPath');
  check(ddg[0].isConnected, 'the cell being clicked is still in the document');
  ddg[0].dispatchEvent(new window.MouseEvent('click', {bubbles: true}));

  await waitFor(() => /relaxed dimer/.test(chosen()), 'the dimer to be shown');
  check(/relaxed dimer/.test(chosen()), `clicking ddG shows the dimer, got ${chosen()}`);
  check(out.state('structPath') !== before,
        'and it is a different file than was shown before');
  check(view.querySelectorAll('tbody tr.selected').length === 1,
        'and the table still marks exactly one sequence');

  await waitFor(() => view.querySelector('.ifacebox'), 'the interface residues');
  const panel = view.querySelector('.ifacebox');
  check(!!panel && /chain A/.test(panel.textContent),
        'and the interface residues are listed per chain');

  // Then a row of another sequence. The dimer was an explicit request about one
  // sequence and does not become the page's new habit.
  const rows = [...view.querySelectorAll('tbody tr:not(.selected)')];
  if (!rows.length) {
    check(false, 'there is another sequence to click');
    return;
  }
  const dimerPath = out.state('structPath');
  rows[0].dispatchEvent(new window.MouseEvent('click', {bubbles: true}));
  await waitFor(() => out.state('structPath') !== dimerPath, 'the row to be selected');
  check(/pentamer|trimer|-mer|relaxed dimer/.test(chosen()),
        `clicking a row shows an assembly, got ${JSON.stringify(chosen())}`);
  check(view.querySelectorAll('tbody tr.selected').length === 1,
        'and exactly one sequence is marked again');
}

/* ----------------------------------------------------------------- the run */

(async () => {
  console.log(`cyclicnano front-end check against ${ORIGIN}`);
  try {
    const health = await globalThis.fetch(`${ORIGIN}/api/health`);
    if (!health.ok) throw new Error(String(health.status));
  } catch (err) {
    console.error(`\n  the service is not answering at ${ORIGIN}: ${clean(err)}`);
    process.exit(2);
  }

  await everyRouteRenders();
  await theStructurePageWorks();

  console.log(`\n  ${failures.length ? `${failures.length} failed` : 'all checks passed'}`);
  process.exit(failures.length ? 1 : 0);
})();
