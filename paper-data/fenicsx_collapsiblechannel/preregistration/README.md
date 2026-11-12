# Preregistration files (written before the runs they cover; timestamps self-asserted, not externally verifiable)

No local filesystem path occurs in any of these files, so the shipped copies are **byte-identical** to the working-tree copies in the solver repository (nothing was removed or edited). The hashes recorded at freeze time (in `.sha256` files and result notes) refer to the file as it was when frozen; where amendments were appended afterwards, the shipped file is the later one and the frozen hash cannot be checked from this package alone (git history of the private solver repository holds the frozen text). Both hashes are given.

| file | commit(s) | sha256 recorded at freeze | sha256 of the shipped copy |
|---|---|---|---|
| `PREREG_COLLAPSIBLE_PRODUCTION.md` (production matrix; the title still says "UNHASHED DRAFT" because the hash fixes that text) | text 7037ba6 (21:53:09), hash file 664de8c (21:53:40), first production submission 21:56:46 (git author times, UTC) | `3b92a32f1882414a876833f4e42a48e70cc56b05f640e0fcb254867e43d30ff6` | same |
| `PREREG_COLLAPSIBLE_PRODUCTION_AMEND1.md` (stage-3 options; written after the stage-1 S1 and S2 cells existed, before every stage-3 job: AMEND1 commit 4f0f6eb 22:21:54, first stage-3 submission 22:22:04) | 4f0f6eb | `b89475b4fd7a9898679f3faf90490d594e2757aff859b6b3771bb8febadb8987` | same |
| `PREREG_COLLAPSIBLE_PRODUCTION_AMEND2.md` (MF1 open-loop replay, before the MF1 runs: commit 6eb35c6 00:36:45, first MF1 job 00:37:10) | 6eb35c6 | `8818abfcbd5938227ec7aeac29e7e113068b8d3f844e7740c2f4fddf29a930b4` | same |
| `PREREG_COL1.md` (structural preflight; amendments A1-A5 appended afterwards) | ebd9f2e (prereg), f065a27 (results with amendments) | `754e17045dbf61918acdd0c370483fb674a07eaf1640b3b597a2af2ac8e8f981` (frozen text, not the shipped file) | `13c781384f51001d87125c43b435fcf54ff59ee35a20389d019d596c693794e9` (with A1-A5) |
| `PREREG_COL2.md` (fluid preflight, beam probe) | 08c052f | `c445f002ee25bc452e584affe39681774a3f131fda29ff2a4bc67fd9cf1665d4` | same |
| `PREREG_COL2_AMEND1.md` / `_AMEND2.md` / `_AMEND3.md` | 15071d3 / a383337 / 0d6ffb9 (A3 written while job 16778 was still running, before the jobs it covers) | `cc30a7ace191b3f1e4de46509e92bbc69cbd696bae714e10f93b480ce66f14af` / `a4f829e47620a3d01adfae18dafb4e63bbc7ffbb164f4539f0bcf2974d547748` / `ad6a0cdd11d92de699ec2cd56fb10c047aa9ca3663d9dca0ccf4f5bb5eda5392` | same |
| `PREREG_COL3_S1.md` (P2-triangle wall gates, with A1 and A2 appended) | 09c0822 (S1 prereg, 21:12:25), A1 4c70255 (before the S2 runs), A2 50b56db | frozen S1 text `519dc03e7d93edafd4a1bc4599ac0b76cdebc1e4d9ba3abea7e235df76706cce`; after A1 `6676fe8c92c40e296cf95a8fc7c278400f9dbe85fd529d15facadf614c3dca83`; after A2 `dccc5b111e3c1d6375d4dc5cdbf73cb1efc24002998f66c78b43196071c450b8` | `dccc5b111e3c1d6375d4dc5cdbf73cb1efc24002998f66c78b43196071c450b8` (after A2) |

Caveats (also in the supplement): the timestamps are self-asserted (git author time and Slurm submit time, no external timestamping); the production runs report `git_dirty` true (untracked run output); AMEND1 was written with knowledge of the stage-1 S1 and S2 results, so its choices (dropping FM-W, WL3 and the first MF1) were not blind.

## Added 2026-11-11 (review request of PR #71): MF1 beam refinement to 320 elements
| file | commit | sha256 (identical to the shipped copy; no paths in the files) |
|---|---|---|
| `PREREG_MF1_BEAM320.md` (written before the run; solver repository c7350ea, this branch cc392c1) | before job 18902 (19:59 UTC) | `c871c6e5d1b65db17fc62c96c3ae460368ebf40790926f7fbed58e18ccea549c` |
| `PREREG_MF1_BEAM320_AMEND1.md` (after the first 320-element run failed its convergence flag, before the rerun; solver repository 9023fe7, this branch ef5652a) | before jobs 18911 and 18912 (20:13 UTC) | `aa7623c6ea843d3a4e761342eb85221eec917ddc57a4ecfc03dbf884ca7b68cd` |
Timestamps self-asserted. The amendment changes only the Newton stopping criterion (residual < 1e-4 and update < 1e-10 m), chosen from the measured residual floor.
