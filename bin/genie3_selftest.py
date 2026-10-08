#!/usr/bin/env python3
"""Self-test for modal_genie3.py -- no GPU, no Modal account, ~0.1s.

WHAT THIS EXISTS TO CATCH, in order of how much it would cost:

1. THE WRONG TARGET. Genie 3's own BinderBench TNF-alpha problem uses 1TNF, which carries LEU at
   mature 143 where the assay construct carries ASP, and BinderBench's `common` interface set
   CONTAINS that residue. A generator pointed at it designs against a different protein at the
   epitope core. bin/build_genie3_problem.py guards the laptop side; the wrapper must also guard
   the container side, on the bytes actually being folded. Both are asserted here.
2. COLABFOLD CREEPING BACK IN. The whole reason this arm is affordable in the remaining time is
   that `genie3 generate` runs without the evaluation stage, so no ColabFold and no MSA-server
   query. An `evaluation:` block in the emitted experiment yaml would silently reintroduce both.
3. AN UNPINNED REVISION. Mosaic ran unpinned for an entire campaign and made METHODS 13's
   reproducibility claim wider than the evidence. A branch name here would repeat it.
4. LOST WORK. Same lesson as esmfold2: output that exists only as a return value dies with a
   timeout, a crash or the local client, and --detach does not help. The Volume commit has to sit
   inside the per-file loop.

usage: python3 bin/genie3_selftest.py       (run from the repo root)
"""
import re
import sys
from pathlib import Path

MOD = Path("biomodals/modal_genie3.py")
BUILDER = Path("bin/build_genie3_problem.py")


def main():
    if not MOD.exists():
        # biomodals/ IS GITIGNORED. The wrapper is tracked only as patches/modal_genie3.patch,
        # so in a fresh clone this file has nothing to check. The old message said "run from
        # the repo root", which sent a reader looking for the wrong problem (found 2026-10-08
        # running the suite in a clean checkout).
        patch = Path("patches/modal_genie3.patch")
        print(f"SKIP  {MOD} is absent. biomodals/ is gitignored by design; the wrapper is",
              file=sys.stderr)
        print(f"      tracked as {patch}" + (" (present)" if patch.exists() else " (MISSING)"),
              file=sys.stderr)
        print(f"      Recreate it with:  mkdir -p biomodals && "
              f"git apply --directory=biomodals {patch}", file=sys.stderr)
        print(f"      This is a SKIP, not a pass: nothing was verified.", file=sys.stderr)
        return 0 if patch.exists() else 1
    src = MOD.read_text()
    lines = src.splitlines()

    # --- 1. the target guard, on both sides ---
    assert 'int(line[22:26]) == 143' in src, "the container has no residue-143 check"
    assert 'aa != "ASP"' in src, "the container does not require ASP143"
    assert "REFUSING" in src, "the container check does not refuse, it only warns"
    b = BUILDER.read_text()
    assert 'SENTINEL = (143, "D")' in b, "the builder's Asp143 sentinel is gone"
    assert "MAX_TARGET_CHAINS = 3" in b, "the target chain count is no longer pinned"
    print("target       Asp143 required in the builder AND re-checked in the container on the")
    print("             bytes being folded; target chains pinned at 3 (max_n_chain 4 - binder)")

    # --- 2. generation only ---
    assert '"generate"' in src, "the wrapper does not call the generate subcommand"
    assert '"run"' not in src, "the wrapper may be calling `genie3 run` (generate + evaluate)"
    # the emitted yaml is a single f-string; it must not declare an evaluation stage
    yaml_blocks = re.findall(r'cfg\.write_text\(f?"""(.*?)"""\)', src, re.S)
    assert len(yaml_blocks) == 1, f"expected exactly one emitted experiment yaml, got {len(yaml_blocks)}"
    y = yaml_blocks[0]
    assert "evaluation:" not in y, (
        "the emitted experiment yaml declares an evaluation stage -- that pulls ColabFold and "
        "queries the MSA server, which is exactly what this arm avoids")
    for key in ("paths:", "rootdir:", "dataset:", "source: target", "cond_strategy:", "n_sample:"):
        assert key in y, f"the emitted yaml is missing {key!r}"
    # ColabFold may be DISCUSSED in prose -- the docstring explains why we avoid it. What must
    # not happen is installing it or naming it in config, so strip comments and docstrings first.
    code = "\n".join(l for l in lines if not l.lstrip().startswith("#"))
    code = re.sub(r'"""(.*?)"""', "", code, flags=re.S)
    assert "colabfold" not in code.lower(), (
        "ColabFold appears in executable code or config, not just in prose")
    install = re.search(r"uv_pip_install\((.*?)\)", src, re.S)
    assert install and "colabfold" not in install.group(1).lower(), (
        "ColabFold is in the image's install list")
    print("scope        `genie3 generate` only; emitted yaml has no evaluation stage, so no")
    print("             ColabFold and no MSA-server query (generation reads no MSA)")

    # --- 3. pinned ---
    ref = re.search(r'GENIE3_GIT_REF = "([0-9a-f]+)"', src)
    assert ref and len(ref.group(1)) == 40, (
        "GENIE3_GIT_REF is not a full 40-character commit sha -- a branch or tag can move")
    assert "@{GENIE3_GIT_REF}" in src, "the image does not install at the pinned ref"
    assert "/pretrained/v1/checkpoints/step=600000.ckpt" in src, "the checkpoint path is not set"
    assert src.count("WEIGHTS}/pretrained") >= 2, (
        "checkpoint and config must both be given as absolute paths under WEIGHTS, or genie3 "
        "resolves its defaults relative to the working directory")
    print(f"pins         genie3 at {ref.group(1)[:12]}..., weights baked at build time,")
    print("             checkpoint and config both absolute under /weights")

    # --- 4. durability ---
    assert "runs_vol.commit()" in src, "nothing ever commits -- the Volume would stay empty"
    assert "modal volume get genie3-runs" in src, "the recovery command is not documented"
    ci = next(i for i, l in enumerate(lines) if "runs_vol.commit()" in l)
    fi = max(i for i, l in enumerate(lines[:ci]) if re.match(r"\s*for .* in files:", l))
    commit_indent = len(lines[ci]) - len(lines[ci].lstrip())
    loop_indent = len(lines[fi]) - len(lines[fi].lstrip())
    assert commit_indent > loop_indent, (
        f"runs_vol.commit() is at indent {commit_indent}, the file loop at {loop_indent} -- "
        "the commit is OUTSIDE the loop, so partial work is still lost")
    print(f"durability   commit at indent {commit_indent} inside the file loop "
          f"(indent {loop_indent}) -- per-sample, not per-run")

    # --- 5. phase 0 must actually measure the two numbers it exists for ---
    for token in ("per_sample_s", "peak_device_gb", "PHASE0.txt", "nvidia-smi"):
        assert token in src, f"the phase-0 measurement is missing {token!r}"
    # 2026-10-08: the first version read torch.cuda.max_memory_allocated() in the PARENT, while
    # genie3 runs in a child process, so it reported 0.00 GB as if it were a measurement.
    assert "max_memory_allocated" not in code, (     # prose may explain the old bug
        "memory is read from this process's torch allocator; genie3 runs in a subprocess, so "
        "that reads 0.00 GB and reports an instrumentation failure as a measurement")
    assert "UNMEASURED" in src, "a failed memory poll must say unmeasured, not report zero"
    print("phase 0      seconds/sample measured, device memory polled via nvidia-smi (NOT the")
    print("             parent allocator), written to PHASE0.txt so it survives the client")

    print("\nselftest OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
