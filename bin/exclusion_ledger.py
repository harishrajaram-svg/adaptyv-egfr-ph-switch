#!/usr/bin/env python3
"""The exclusion ledger, keyed on binder SEQUENCE and classified by reason.

WHY THIS REPLACES METHODS SECTION 10'S PROSE LEDGER. In a note of 2026-10-05 the reviewer
asked that the 38-molecule ledger be opened up again: de-duplicate its aliases on exact
binder sequence, work up the further high-ranking molecules that had been left out, and
measure how far the lists overlap before any counts are added together. They also asked
that three failure modes be kept strictly apart -- ineligibility, assessment that was
never adequate, and rejections resting on the pH gate alone -- and that the last of those
be returned to consideration, with one documented assessment applied uniformly across
that reopened pool and across the designs currently submitted.

Section 10 held two lists and never reconciled them:

  LIST A  45 run names check_discards warns on: measured, outranking a submitted design,
          rejected anyway. Section 10 says these are 38 distinct molecules -- the gate is
          NAME-keyed on the discard side, so alias pairs double-count. The ledger built to
          catch this project's keystone trap was itself keyed on names.
  LIST B  molecules that outrank the weakest submitted design on the pooled
          sequence-keyed table and were never candidates. Called 28, then 59.

Adding 38 and 59 into "97 excluded" is wrong if the lists intersect, and nobody had
checked. This establishes the overlap first, then classifies the union into exactly one of
three reasons, which are not equivalent:

  ELIGIBILITY   fails a rule that disqualifies it whatever it scored -- novelty Level < 3
                or an expression-QC flag. Not reopenable.
  ASSESSMENT    never adequately measured: < MIN_N poses, or no affinity on either
                species. Excluded for want of evidence, not on evidence.
  GATE-ONLY     eligible and adequately measured, excluded only by a pH-gate threshold or
                a judgement call. REOPENED: per the review, pH-gate-dependent rejections
                are "unsupported by a validated selection rule", so a gate-only exclusion
                is not a finding.

--gate then applies the SAME two-partner histidine-only gate that ranks the finalists to
every reopened molecule over its own human-leg poses, closing the gap Section 10 admitted:
"all 59 are ranked on the SUPERSEDED target-only basis -- the all-site product was
computed only for the twelve submitted designs."

    exclusion_ledger.py           # reconcile + classify (fast)
    exclusion_ledger.py --gate    # also run the uniform gate over the reopened pool (slow)
"""
import csv, glob, json, os, re, statistics as st, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pose_seq_index
import check_discards as cd

MR = 'analysis/01-egfr/master_rank.json'
CSV = 'submissions/01-egfr.csv'
NOV = 'analysis/01-egfr/novelty_LEVELS.tsv'
QC = 'analysis/01-egfr/express_qc.tsv'
OUT = 'analysis/01-egfr/exclusion_ledger.json'
# Gate results cached BY SEQUENCE so reclassifying is free. The gate is ~1.2s a
# pose over ~400 poses; without this, every change to the classification rules
# costs another 8 minutes and the temptation is to stop re-deriving it.
GATE_CACHE = 'analysis/01-egfr/exclusion_gate_cache.json'
MIN_N = 5


def eligibility_by_name():
    """design name -> (novelty verdict, qc flags). Both source tables are NAME-keyed.

    Reads EVERY analysis/01-egfr/novelty_*.tsv, not just novelty_LEVELS.tsv. The first
    version of this function read that one file, which holds 16 names, and the resulting
    ledger reported "ELIGIBILITY: 0" over a 75-molecule union at 0% join coverage -- a
    vacuous zero presented as a clean result, which is the exact failure this ledger
    exists to catch. A `level` column is present in only two of the twelve files; where it
    is absent the boolean passes/clears_gate columns are used and the verdict is recorded
    as derived rather than read.
    """
    nov, flags = {}, {}
    for f in sorted(glob.glob('analysis/01-egfr/novelty_*.tsv')):
        try:
            rows = list(csv.DictReader(open(f), delimiter='\t'))
        except OSError:
            continue
        for r in rows:
            d = r.get('design') or r.get('name') or ''
            nm = os.path.basename(d).rsplit('.', 1)[0]
            if not nm:
                continue
            rec = nov.setdefault(nm, {'level': None, 'passes': None, 'src': []})
            rec['src'].append(os.path.basename(f))
            if r.get('level') not in (None, ''):
                try:
                    lv = int(r['level'])
                    rec['level'] = lv if rec['level'] is None else max(rec['level'], lv)
                except (ValueError, TypeError):
                    pass
            for k in ('passes', 'clears_gate'):
                if r.get(k) not in (None, ''):
                    v = str(r[k]).strip().lower() in ('true', '1', 'yes')
                    rec['passes'] = v if rec['passes'] is None else (rec['passes'] or v)
    if os.path.exists(QC):
        for r in csv.DictReader(open(QC), delimiter='\t'):
            flags[r['name']] = (r.get('flags') or '').strip()
    return nov, flags


def classify(row, nov, flags, poses):
    """Exactly one reason, in precedence order: eligibility, then assessment, then gate.

    Also returns an eligibility_status string, because "assessed and passed" and "never
    assessed" are different states and collapsing them lets an unchecked molecule read as
    eligible. Expression QC was only ever run on the 20 submission candidates, so for most
    of the excluded pool that axis is genuinely unassessed and says so.
    """
    names = row.get('names') or []
    recs = [nov[n] for n in names if n in nov]
    levels = [r['level'] for r in recs if r['level'] is not None]
    passes = [r['passes'] for r in recs if r['passes'] is not None]
    qc = [flags[n] for n in names if n in flags]
    bad_qc = [n + ':' + flags[n] for n in names
              if n in flags and flags[n] and flags[n] != 'clean']

    if levels:
        elig = 'level >= 3' if max(levels) >= 3 else 'level ' + str(max(levels))
    elif passes:
        elig = 'derived from passes=' + str(any(passes))
    else:
        elig = 'NOT ASSESSED (no novelty record)'
    if not qc:
        elig += '; QC NOT RUN'
    elif bad_qc:
        elig += '; QC ' + bad_qc[0]

    if levels and max(levels) < 3:
        return 'ELIGIBILITY', 'novelty Level ' + str(max(levels)) + ' < 3', elig
    if passes and not any(passes) and not levels:
        return 'ELIGIBILITY', 'fails the novelty gate (derived from passes)', elig
    if bad_qc:
        return 'ELIGIBILITY', 'expression QC: ' + bad_qc[0], elig
    if (row.get('ratio_n') or 0) < MIN_N:
        return 'ASSESSMENT', 'only ' + str(row.get('ratio_n') or 0) + ' pH poses (< ' + str(MIN_N) + ')', elig
    if not poses:
        return 'ASSESSMENT', 'no human-leg pose on disk; cannot apply the uniform gate', elig
    if row.get('hu_med') is None and row.get('mo_med') is None:
        return 'ASSESSMENT', 'no affinity measurement on either species', elig
    return 'GATE-ONLY', 'eligible and measured; excluded on the pH gate or a judgement call', elig


def main():
    t0 = time.time()
    sub = {r['sequence'].strip().upper(): r for r in csv.DictReader(open(CSV))}
    mr = json.load(open(MR))
    by_seq = {r['seq']: r for r in mr}
    nov, flags = eligibility_by_name()
    idx = pose_seq_index.load()

    def human_poses(seq):
        out = []
        for d in idx.get(seq, []):
            if os.path.basename(d.rstrip('/')).endswith(('_hu', '_human')):
                out += sorted(glob.glob(os.path.join(d, '*.cif')))
        return sorted(set(out))

    # ---- LIST A: the 45 warn names, resolved to sequences -------------------
    name2seq = cd.design_seqs()
    prim = cd.primary_values()
    measured = cd.measured_on_instrument()
    worst_sub_tgt = min(float(r['ph_ratio_target_only_SUPERSEDED']) for r in sub.values())
    listA_names, unresolved = [], []
    for d, v in sorted(prim.items(), key=lambda kv: -kv[1]):
        if v < worst_sub_tgt:
            continue
        q = name2seq.get(d) or name2seq.get(cd.stub(d))
        if q is None:
            unresolved.append((d, v)); continue
        if q in sub:
            continue                       # a submitted design is not an exclusion
        if q in measured:
            listA_names.append((d, v, q))
    listA_seqs = {q for _, _, q in listA_names}

    # ---- LIST B: outrank the weakest finalist on the pooled table -----------
    listB_seqs = {r['seq'] for r in mr
                  if r['seq'] not in sub and (r.get('ratio') or 0) >= worst_sub_tgt
                  and (r.get('ratio_n') or 0) >= MIN_N}

    both = listA_seqs & listB_seqs
    union = listA_seqs | listB_seqs

    print(f"LIST A  45 warn run names -> {len(listA_names)} names over "
          f"{len(listA_seqs)} DISTINCT SEQUENCES "
          f"({len(listA_names) - len(listA_seqs)} alias collapses)")
    if unresolved:
        print(f"        {len(unresolved)} warn name(s) could not be resolved to a sequence: "
              f"{', '.join(n for n, _ in unresolved[:4])}")
    print(f"LIST B  {len(listB_seqs)} sequences outranking the weakest finalist "
          f"({worst_sub_tgt:.3f}x target-only) with n >= {MIN_N}")
    print(f"OVERLAP {len(both)} sequence(s) in both lists")
    note = ("the lists are DISJOINT, so the union is their sum"
            if not both else
            f"NOT {len(listA_seqs)} + {len(listB_seqs)} = "
            f"{len(listA_seqs)+len(listB_seqs)}; {len(both)} are counted in both")
    print(f"UNION   {len(union)} distinct excluded molecules -- {note}")

    rows = {}
    for q in sorted(union):
        r = by_seq.get(q) or {'seq': q, 'names': []}
        poses = human_poses(q)
        cat, why, elig = classify(r, nov, flags, poses)
        rows[q] = dict(seq=q, names=(r.get('names') or [])[:4],
                       in_list_a=q in listA_seqs, in_list_b=q in listB_seqs,
                       ratio_target_only=r.get('ratio'), ratio_n=r.get('ratio_n') or 0,
                       hu_med=r.get('hu_med'), mo_med=r.get('mo_med'),
                       n_human_poses=len(poses), category=cat, reason=why,
                       eligibility_status=elig)
    from collections import Counter
    c = Counter(v['category'] for v in rows.values())
    print("\nCLASSIFICATION of the union:")
    for k in ('ELIGIBILITY', 'ASSESSMENT', 'GATE-ONLY'):
        print(f"  {k:<12} {c.get(k,0):>3}")
    reopened = [v for v in rows.values() if v['category'] == 'GATE-ONLY']
    print(f"\nREOPENED (gate-only): {len(reopened)} molecules, "
          f"{sum(v['n_human_poses'] for v in reopened)} human-leg poses to gate")

    cache = {}
    if os.path.exists(GATE_CACHE):
        cache = json.load(open(GATE_CACHE))
    for v in rows.values():                       # attach whatever is already cached
        if v['seq'] in cache:
            v['gate'] = cache[v['seq']]

    if '--gate' in sys.argv:
        from ph_gate_multisite import score_pose
        todo = [v for v in reopened if v['seq'] not in cache]
        print(f"\napplying the uniform histidine-only gate to the reopened pool "
              f"({len(cache)} cached, {len(todo)} to compute)...", flush=True)
        for i, v in enumerate(sorted(todo, key=lambda x: -(x['ratio_target_only'] or 0)), 1):
            per = []
            for cif in human_poses(v['seq']):
                r = score_pose(cif)
                if r and 'error' not in r:
                    per.append(r)
            if not per:
                v['gate'] = {'error': 'every pose failed to score'}
                cache[v['seq']] = v['gate']; continue
            his = [p['product_his_only'] for p in per]
            med = st.median(his)
            v['gate'] = dict(n_poses=len(per), his_only_median=round(med, 4),
                             his_only_min=round(min(his), 4), his_only_max=round(max(his), 4),
                             spread_over_median=round((max(his)-min(his))/med, 3) if med else None)
            cache[v['seq']] = v['gate']
            json.dump(cache, open(GATE_CACHE, 'w'))     # checkpoint every molecule
            print(f"  [{i}/{len(todo)}] {(v['names'] or ['?'])[0][:44]:<46} "
                  f"tgt {v['ratio_target_only']:>6.3f} -> his-only {med:>7.3f} "
                  f"(n={len(per)}, spread {v['gate']['spread_over_median']})", flush=True)

    json.dump(dict(generated='2026-10-05', min_n=MIN_N,
                   worst_finalist_target_only=worst_sub_tgt,
                   list_a_names=len(listA_names), list_a_sequences=len(listA_seqs),
                   list_b_sequences=len(listB_seqs), overlap=len(both), union=len(union),
                   unresolved_warn_names=[n for n, _ in unresolved],
                   molecules=rows), open(OUT, 'w'), indent=1)
    print(f"\nwrote {OUT} in {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
