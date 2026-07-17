# DEPRECATED as a standalone repository

As of **2026-07-17**, `GRiDCodeGenerator` is no longer developed as a separate
repository. It has been folded directly into
[**A2R-Lab/GRiD**](https://github.com/A2R-Lab/GRiD) under the
`GRiDCodeGenerator/` directory, via a history-preserving `git subtree` merge
(the full commit history is retained inside GRiD).

## Where to develop now

- **Repo:** https://github.com/A2R-Lab/GRiD
- **Path:** `GRiDCodeGenerator/`

Codegen edits and their GRiD consumers are now a single atomic commit — no more
submodule-pointer bump.

## Status of this repository

This standalone repository is kept only as a **frozen mirror** of the pre-fold
history. Do not open PRs or push new work here — send all changes to
`A2R-Lab/GRiD`.
