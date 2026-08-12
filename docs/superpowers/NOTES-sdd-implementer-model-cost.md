# Note: is Opus-low a cheaper SDD implementer than Sonnet-high?

Date: 2026-07-28. Input for a plan; nothing built yet.

## Verdict

Plausible, not established. The price sheet does **not** settle it, and the usual
argument for it (low effort = fewer thinking tokens) aims at the smallest part of
the bill.

## Break-even is one number

Opus 5 $5/$25 per MTok; Sonnet 5 $3/$15 ($2/$10 intro **through 2026-08-31**).
The Opus:Sonnet ratio is identical on input and output (5/3 = 25/15), so the
in/out mix cancels and the break-even is just the price ratio:

| Window | Opus:Sonnet per token | Opus-low wins iff it burns |
|---|---|---|
| now → Aug 31 | 2.50× | **< 40%** of Sonnet-high's total tokens |
| after Aug 31 | 1.67× | **< 60%** of Sonnet-high's total tokens |

## Where the tokens actually go (measured, not assumed)

From `~/.claude/projects/*/*.jsonl`, 13,829 Opus 4.8 turns + 1,036 Sonnet 4.6
turns (structural proxy for the 5-series):

| model | turns | $ cache read | $ cache write | $ output | output share |
|---|---:|---:|---:|---:|---:|
| opus-4-8 | 13,829 | 1266 | 1233 | 536 | **18%** |
| sonnet-4-6 | 1,036 | 29 | 25 | 11 | **17%** |

~1,500 output tokens/turn against ~191,000 input tokens/turn. **~82% of spend is
the input side** — context re-sent every turn. Halving thinking tokens saves ~9%,
nowhere near the 40–60% needed.

So the real lever is **turn count**, not thinking depth. Low effort helps only if
it consolidates tool calls enough to cut turns; each turn re-bills the whole context.

## The variable that probably decides it

Our escalation rule (CLAUDE.md) sends BLOCKED / NEEDS_CONTEXT / DONE_WITH_CONCERNS
to `sdd-implementer-escalated` (Opus **high**), which pays for the failed low run too.
Roughly: with `r` = Opus-low's token ratio and `p` = escalation rate, need `r + p < 0.6`.
A 20% escalation rate means low effort must run under 40% of tokens just to break even —
and low effort is what raises `p`.

## Two facts that block naive measurement

1. **`sdd-implementer` has never actually been dispatched.** 325 real subagent
   dispatches across all projects: 264 `general-purpose`, 9 `claude`, 52 default.
   Zero `sdd-implementer`. Those SDD runs inherited the parent's model/effort, so
   the pinned Opus-low config is untested.
2. **Subagent turns are not logged in the transcripts** — 0 sidechain records — so
   per-subagent usage cannot be read from the JSONL.

Measurement channel is `~/.claude/stats-cache.json`: `modelUsage` (cumulative
per-model) and `dailyModelTokens` (per-day per-model). Opus and Sonnet are distinct
model IDs, so a mixed run separates cleanly by model. `/cost` also works per session.

## To build

- `~/.claude/agents/sdd-implementer-sonnet.md` — sibling of `sdd-implementer`,
  identical body, frontmatter `model: sonnet` / `effort: high`. Does not exist yet.
- Snapshot/diff script over `stats-cache.json` → per-task cost at $5/$25 vs $3/$15.
- A/B: one plan, ~10 tasks alternating between the two agents in one session.
  Record escalation count per arm — that is the real variable, not tokens.

---

## Appendix: plan status (from the `other_libs/climt` checkout)

Checkbox state in `docs/superpowers/plans/` is **stale and unreliable** — the
land/ocean/ice plans read "not started" but shipped in PR #219. Verified against
the repo instead:

**Open / outstanding**

- `2026-06-16-cork-table-optimizer.md` — not started (1/22 target files exist).
- `2026-05-18-climt-docs-feature.md` — not started (3/17).
- `2026-07-19-in-browser-nongrey-rce-demo.md` — in progress; current branch
  `feature/pyodide-cork-prep` is this work (25/41).
- `2026-07-22-boundary-layer-and-jit-tridiagonal.md` — Tasks 1–3 done
  (`_core/tridiagonal.py` + `simple_boundary_layer/` exist), **Tasks 4–5 not done**:
  `scipy` is still imported in `_core/snow_ice_column.py` and
  `_components/second_best/processes/subsurface.py`.
- `2026-05-16-cork-co2-band-refinement.md` — 4/5 tasks.
- `2026-04-20-cork-radiation-phase4.md` — 9/26 tasks.

**Cleanup**

Three `picket-fence`/`pf-` plan files are pre-rename copies of their `cork` twins
(phase-3 pair differs by 32 lines out of 1435: `PicketFence`→`Cork`,
`CorkLongwave`→`CorkLongwaveRadiation`). They will mislead any agent that reads
them. Candidates for deletion:

- `2026-04-14-picket-fence-radiation-phase3.md`
- `2026-04-20-picket-fence-radiation-phase4.md`
- `2026-06-04-pf-co2-adjustable-hifi-earth-lw.md`
