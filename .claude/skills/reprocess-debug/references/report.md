# Publishing the post-mortem

Distilled from `docs/agent_debug.md` §12. Reference examples kept at
`data/failure-postmortem.html`, `data/failure-postmortem-batch5.html`,
`data/failure-postmortem-run3.html`.

A failure post-mortem is a deliverable with an audience, so it goes out as an Artifact — a
private page on claude.ai the user can choose to share — not as terminal scrollback.

## §procedure

1. **Load the `artifact-design` skill before writing a single line of the file.** Required,
   Markdown reports included. Format is part of the design decision.
2. Run `bin/triage.py --json` first and build the page from that JSON, not from numbers retyped
   out of the terminal. Retyping is where reconciliation errors come from.
3. Write the HTML with `Write`, keeping the source in the repo as
   `data/failure-postmortem-batch<N>.html` so a later session can edit and redeploy rather than
   rebuild. (The existing files are inconsistently named — `-batch5`, `-run3`,
   `run-673887-`; use `-batch<N>` for new ones.)
4. Call `Artifact` with `file_path`, a one-sentence `description` (it becomes the gallery card
   subtitle), and a `favicon`. `🧬` was used for the run reports — keep it, so the family is
   recognisable.
5. To update: **edit the same file and call `Artifact` with the same `file_path`** — it
   redeploys to the same URL. A different path claims a new URL. From a *different*
   conversation you must pass the artifact's `url` explicitly, or you create a duplicate.
6. **Record the URL in `references/run-index.md`**, in the run's row. This is the step that gets
   skipped, because once the artifact is published the work feels finished — and then the next
   session has to list every artifact and guess from titles. Look the URL up there rather than
   with `action: "list"`.

## §mechanics

* **Write page content only.** No `<!DOCTYPE>`, `<html>`, `<head>` or `<body>` — the file is
  wrapped in that skeleton at publish time. Start with `<title>`, then `<style>`, then markup.
* **`<title>` must be in the first 8 KB** and should be a *name* of two to four words, not a
  summary and not a name-plus-explainer after a dash. The explanation goes in `description`.
  Keep title and favicon stable across redeploys.
* **Self-contained.** A strict CSP blocks every external host except Google Fonts. Inline all
  CSS and JS. No CDN scripts, no remote images, no `fetch`. Downloads the page starts itself
  are inert in the viewer sandbox, so never offer a file via a link.
* **Wide content scrolls inside its own box.** Every table and code block goes in a wrapper with
  `overflow-x: auto`, or the body scrolls sideways on a phone. The failure tables need this.
* **Mermaid renders natively** via `<pre class="mermaid">`; no library needed.

## §theme

The one bug that makes a report unreadable rather than merely ugly. The viewer has **three**
states: an explicit choice stamps `data-theme="dark"`/`"light"` on the root, and the default
"system" setting **stamps nothing** — there, only `prefers-color-scheme` distinguishes them.

```css
:root { /* complete light palette, every token defined here */ }
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) { /* redefine only the tokens */ }
}
:root[data-theme="dark"] { /* redefine them again so the toggle wins both ways */ }
```

* Never give a colour its only definition inside a media or `[data-theme]` block — in the
  un-stamped state it never applies, and you get one theme's text on the other's ground.
* Style components through tokens only.
* `body` needs an explicit `background` from a token. A transparent body borrows the host page's
  ground and breaks in one theme.
* Before publishing, scan the stylesheet for any colour declared only behind a media or
  `[data-theme]` selector.

## §house-style

Reuse this so successive post-mortems read as one family. A cool, green-biased neutral set —
chosen for the subject — with semantic colours kept separate from the accent.

```css
/* light */
--paper:#f6f8f5; --surface:#fff;    --surface-2:#eef2ee;
--ink:#16211d;   --ink-2:#3d4b45;   --muted:#6c7a73;
--rule:#d3dcd6;  --rule-soft:#e4eae5;
--accent:#1c6b58; --accent-dim:#2f8570;
--critical:#a1382a; --warn:#97701f; --ok:#2f6f4a;
/* dark */
--paper:#0f1512; --surface:#161e1a; --surface-2:#1c2621;
--ink:#e3eae5;   --ink-2:#b3c0b9;   --muted:#86948c;
--rule:#2c3a33;  --rule-soft:#232f29;
--accent:#58c0a5; --accent-dim:#3f9b83;
--critical:#e08272; --warn:#d6ae5c; --ok:#74c193;
```

Type: **monospace carries the structure** — headings, eyebrows, labels, every number — with a
serif for prose. That inversion suits a report made of accessions, exit codes and percentages.
No webfonts; system stacks only:

```css
--mono: ui-monospace,"SF Mono",SFMono-Regular,"Cascadia Mono","Roboto Mono",Menlo,Consolas,monospace;
--serif: Charter,"Bitstream Charter","Sitka Text",Cambria,"Noto Serif",Georgia,serif;
```

Devices that carried real information:

* Single column, `max-width: 60rem`; prose capped near `44rem`.
* `font-variant-numeric: tabular-nums` everywhere digits line up.
* A **scoreboard** of headline counts at the top.
* A **stage strip** of failures per pipeline stage, top-border colour-coded by severity — this
  is where "SRA2FASTQ shows 0 failures but caused 141" becomes visible.
* A **stacked verdict bar** splitting correct-rejection / too-strict / pipeline-bug, with `flex`
  weights set to the raw counts so the widths *are* the data. `triage.py`'s `by_verdict` gives
  these directly.
* Findings as cards with a severity-coloured left border, a big tabular count, a verdict chip.
* Numbered markers **only** where order is real — a priority list, not decoration.
* Fixed-width evidence blocks holding verbatim log lines and measured numbers.
* Lay out sibling groups with flex/grid `gap`, not per-element margins.

## §before-publishing

* **Make the numbers reconcile.** Category counts must sum to the total; a stacked bar's
  segments must sum to the whole. `triage.py` asserts this for its own output — carry the same
  discipline into the page, and check it in code, not by eye.
* **Separate permanent from recovered.** Every headline number should be the permanent count,
  and the page should say so. Batch 3's raw 111 is 92 once self-healed tasks are excluded.
* **State the method and its limits.** Say how categories were assigned, that every task is
  accounted for, and which claims rest on a spot check rather than the full set. Mark estimates
  as estimates.
* Structural check for unclosed tags:

```bash
python3 -c "
import html.parser
class P(html.parser.HTMLParser):
    def __init__(s): super().__init__(); s.st=[]
    def handle_starttag(s,t,a):
        if t not in ('br','img','hr','meta','link','input'): s.st.append(t)
    def handle_endtag(s,t):
        if s.st and s.st[-1]==t: s.st.pop()
        else: print('MISMATCH',t,'near',s.st[-3:])
p=P(); p.feed(open('data/failure-postmortem-batch6.html').read()); print('unclosed:',p.st)"
```

* **Relay the findings in the chat reply too.** The artifact's contents are not shown to the
  user automatically — hand over the link *and* the headline conclusions.

Never publish a page that impersonates a real person or organisation, or presents fabricated
records as genuine. Reports on your own analysis of this pipeline are fine.
