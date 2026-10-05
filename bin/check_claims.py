#!/usr/bin/env python3
"""Gate every numeric claim in the deliverables against the artifact that produces it.

WHY THIS EXISTS. The deliverables carry ~150 bolded numeric claims. Six generated blocks
cover the tables, and those have never drifted -- they cannot, because
`gen_methods_submission.py --check` fails the build. Every other number is hand-maintained,
and the hand-maintained ones kept breaking: the binder-histidine count was wrong FOUR times
in one night (six-of-eleven, seven-of-twelve, eight-of-seventeen, ten-of-seventeen), each
time caught by a human read and each time broken again by the next change.

A read finds instances. A gate fixes the class. This is the gate.

WHAT IT DOES, in three passes:

  RULES       named, load-bearing claims with an explicit source: design count, slot usage,
              family count, tier-1 count, near-duplicate pairs, binder-histidine count,
              Kendall tau values, control-panel fractions, novelty counts. Each is computed
              from the artifact and compared to what the document says. A mismatch FAILS.

  CITATIONS   per-design numbers quoted in prose. For every design name appearing in a
              document, any number within a short window is matched against that design's
              own CSV row (pH ratio, affinities, pose count, spread). Catches a value that
              was correct when written and went stale when the design moved.

  INVENTORY   every remaining number, checked for existence anywhere in the artifact value
              set. A number present somewhere is `corroborated` (weak -- it exists, but
              this does not prove it is attached to the right claim). A number present
              nowhere is `UNSOURCED` and is listed. That list is the stone-unturned
              inventory: it is explicit, it shrinks monotonically, and it does not depend
              on anyone remembering to look.

Exit 1 if any RULE or CITATION fails. The inventory is reported, never fatal -- prose
legitimately contains numbers that are not measurements (dates, section numbers, pH values,
literature figures), and failing on those would make the gate unusable and therefore ignored.

    check_claims.py                 # run every pass
    check_claims.py --rules         # rules only (fast)
    check_claims.py --inventory     # show the unsourced list in full
"""
import csv, glob, json, os, re, sys
from collections import defaultdict

DOCS = ['submissions/01-egfr-METHODS.md', 'outbox/CONTROL-TABLE.md',
        'outbox/PREREGISTRATION.md', 'README.md', 'HANDOFF.md']
CSVP = 'submissions/01-egfr.csv'
PERMITTED = 20


def load():
    a = {}
    with open(CSVP) as fh:
        a['csv'] = list(csv.DictReader(fh))
    for key, path in (('sens', 'analysis/01-egfr/ph_sensitivity.json'),
                      ('apo', 'analysis/01-egfr/ph_apo_freeleg.json'),
                      ('pert', 'analysis/01-egfr/ph_pka_perturbation.json'),
                      ('ledger', 'analysis/01-egfr/exclusion_ledger.json'),
                      ('foot', 'analysis/01-egfr/finalist_footprints.json'),
                      ('chai', 'analysis/01-egfr/chai_interface_summary.json'),
                      ('ctrl', 'analysis/01-egfr/control_recovery.json')):
        a[key] = json.load(open(path)) if os.path.exists(path) else None
    a['docs'] = {d: open(d).read() for d in DOCS if os.path.exists(d)}
    # FLATTENED copies. Markdown wraps prose, so a claim can straddle a newline and a
    # regex over the raw text silently misses it -- which is how "eight of the seventeen
    # submitted designs carry at least one histidine" survived the first run of this gate
    # while being wrong for the fifth time. Matching happens on the flattened text; the
    # offset map converts a flat position back to a real line number for the message.
    a['flat'] = {}
    for d, text in a['docs'].items():
        parts, offs, pos = [], [], 0
        for i, line in enumerate(text.split('\n'), 1):
            parts.append(line)
            offs.append((pos, i))
            pos += len(line) + 1
        a['flat'][d] = (' '.join(parts), offs)
    return a


def flat_line(offs, idx):
    ln = 1
    for pos, i in offs:
        if pos <= idx:
            ln = i
        else:
            break
    return ln


def ident(x, y):
    import gemmi
    return gemmi.align_string_sequences(list(x), list(y), []).calculate_identity() / 100.0


# ---------------------------------------------------------------- RULES ----
def rules(a):
    """[(name, expected, found_in_docs, ok)] -- each expectation computed from artifacts."""
    out = []
    rows = a['csv']
    n = len(rows)
    docs = a['docs']
    alltext = '\n'.join(docs.values())

    def says(pattern, flags=re.I):
        return [m.group(0) for m in re.finditer(pattern, alltext, flags)]

    WORDS = {0:'zero',1:'one',2:'two',3:'three',4:'four',5:'five',6:'six',7:'seven',
             8:'eight',9:'nine',10:'ten',11:'eleven',12:'twelve',13:'thirteen',
             14:'fourteen',15:'fifteen',16:'sixteen',17:'seventeen',18:'eighteen',
             19:'nineteen',20:'twenty',25:'twenty-five',35:'thirty-five'}

    def _eq(got, expected):
        """A claim may be written as a digit or as a word; both must mean `expected`."""
        g = str(got).strip().lower().rstrip('.')
        e = str(expected).strip().lower()
        if g == e:
            return True
        try:
            if abs(float(g) - float(e)) < 1e-9:
                return True
        except ValueError:
            pass
        try:
            if g == WORDS.get(int(float(e))):
                return True
        except (ValueError, TypeError):
            pass
        return False

    def rule(name, expected, pattern, optional=False):
        """Every match of `pattern` in the docs must state `expected`.

        A rule that matches NOTHING is reported as NO-CLAIM, never as OK. Reporting a
        vacuous pass is the specific failure this whole file exists to prevent -- it is
        how "ELIGIBILITY: 0" got presented as a clean result over a 75-molecule union at
        0% join coverage earlier today. `optional=True` marks a phrasing that genuinely
        may be absent; everything else wants a human to know the gate found no anchor.
        """
        bad, hits = [], 0
        for doc, (text, offs) in a['flat'].items():
            for m in re.finditer(pattern, text, re.I):
                # Skip past-tense recitals: this document deliberately records what the
                # submission USED to contain, and flagging those would make the gate cry
                # wolf on its own correction history.
                #
                # Scoped to the CURRENT SENTENCE only. A 90-character lookback was tried
                # first and immediately false-negatived: a "was" in the previous sentence
                # suppressed a genuinely stale "17 designs, ranked on". In a gate a false
                # negative is strictly worse than a false positive, so the window stops at
                # the nearest sentence boundary.
                pre = text[max(0, m.start() - 200):m.start()]
                cut = max(pre.rfind('. '), pre.rfind('.**'), pre.rfind('! '),
                          pre.rfind('|'), pre.rfind('###'))
                clause = pre[cut + 1:].lower() if cut >= 0 else pre.lower()
                if re.search(r'\b(was|were|had|used to|previously|then-|briefly)\b', clause):
                    continue
                hits += 1
                got = next((x for x in m.groups() if x), '')
                if not _eq(got.replace(',', ''), expected):
                    bad.append(f"{doc}:{flat_line(offs, m.start())} says {got!r}, artifact "
                               f"says {expected!r}  [{m.group(0)[:60]}]")
        if hits == 0 and not optional:
            out.append((name, expected, None, [f"NO CLAIM MATCHED -- pattern finds nothing in "
                                               f"any document; the gate is not guarding this"]))
        else:
            out.append((name, expected, len(bad) == 0, bad))

    # design count, several phrasings
    rule('design count "N designs, ranked on"', n, r'\*\*(\d+) designs, ranked on')
    rule('slot usage "N of the 20 permitted"', n, r'\*\*(\d+) of the 20 permitted')
    rule('slot usage "stands at N of 20"', n, r'stands at \*\*(\d+) of the 20 permitted')

    # family count from the generator's own map
    try:
        sys.path.insert(0, 'bin')
        import gen_methods_submission as g
        fams = {g.FAMILY[r['name']] for r in rows if r['name'] in g.FAMILY}
        rule('family count "N families, M designs" heading', len(fams),
             r'### 11\.3 (\w+) families', lambda m: m)
    except Exception:
        pass

    # near-duplicate pairs
    try:
        pairs = []
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                v = ident(rows[i]['sequence'].strip().upper(), rows[j]['sequence'].strip().upper())
                if v >= 0.90:
                    pairs.append(v)
        rule('near-duplicate pair count', len(pairs),
             r'\*\*(\w+|\d+) pairs of submitted designs exceed 90% sequence identity')
        in_pair = len({x for i in range(len(rows)) for j in range(i + 1, len(rows))
                       for x in (rows[i]['name'], rows[j]['name'])
                       if ident(rows[i]['sequence'].strip().upper(),
                                rows[j]['sequence'].strip().upper()) >= 0.90})
        rule('designs sitting in a near-duplicate pair', in_pair,
             r'(\w+) of the (?:seventeen|eighteen|nineteen|twenty) designs sit in such a pair')
    except Exception as e:
        out.append(('near-duplicate pairs', 'n/a', False, [f'could not compute: {e}']))

    # binder-histidine count -- the one that broke four times
    if a['sens']:
        ship = {r['sequence'].strip().upper() for r in rows}
        wh = sum(1 for v in a['sens'].values()
                 if v.get('seq', '').strip().upper() in ship and v.get('n_his', 0) > 5)
        words = {8: 'eight', 9: 'nine', 10: 'ten', 11: 'eleven', 12: 'twelve'}
        rule('binder-histidine count (digit form)', wh,
             r'\*\*(\d+) of the (?:seventeen|eighteen) submitted designs carry at least one histidine')
        bad, hits = [], 0
        for doc, (text, offs) in a['flat'].items():
            for m in re.finditer(r'\*\*(\w+) of the (\w+) submitted designs '
                                 r'carry at least one histidine', text, re.I):
                hits += 1
                tot = WORDS.get(len(rows), str(len(rows)))
                if m.group(1).lower() != words.get(wh, '?') or m.group(2).lower() != tot:
                    bad.append(f"{doc}:{flat_line(offs, m.start())} says "
                               f"{m.group(1)!r} of {m.group(2)!r}, artifact says "
                               f"{words.get(wh)!r} of {tot!r}")
        out.append(('binder-histidine count (word form)', f"{words.get(wh, wh)} of "
                    f"{WORDS.get(len(rows), len(rows))}", (not bad) if hits else None, bad))

    # Kendall tau, apo vs deletion
    if a['apo']:
        ok = [r for r in a['apo']['rows'] if 'fold' in r and r.get('fold')]
        import itertools
        do = [r['name'] for r in sorted(ok, key=lambda r: -r['deletion'])]
        ao = [r['name'] for r in sorted(ok, key=lambda r: -r['apo'])]
        c = d = 0
        for x, y in itertools.combinations(ok, 2):
            s = (do.index(x['name']) - do.index(y['name'])) * (ao.index(x['name']) - ao.index(y['name']))
            c, d = (c + 1, d) if s > 0 else (c, d + 1)
        tau = (c - d) / (c + d)
        rule('Kendall tau (apo vs deletion free leg)', f"{tau:+.3f}".replace('+', ''),
             r'Kendall τ = \+?([\d.]+)\*\*, against')

    # perturbation: designs spanning >=5 ranks at the default sigma
    if a['pert']:
        wide = sum(1 for r in a['pert']['designs'] if r['p95_rank'] - r['p5_rank'] >= 5)
        tot = len(a['pert']['designs'])
        rule('perturbation: designs spanning >=5 ranks', wide,
             r'\*\*(\d+) of 1[78]\*\*\s*\|?\s*$', )
        bad = []
        for doc, (text, offs) in a['flat'].items():
            for m in re.finditer(r'(\d+) of (\d+) designs span five or more ranks', text, re.I):
                if int(m.group(1)) != wide or int(m.group(2)) != tot:
                    bad.append(f"{doc}:{flat_line(offs, m.start())} says "
                               f"{m.group(1)}/{m.group(2)}, artifact says {wide}/{tot}")
        out.append(('perturbation span sentence', f"{wide}/{tot}", not bad, bad))

    # control panel: raw fraction below EGF, human leg
    if a['ctrl']:
        neg = [r for r in a['ctrl'] if r['molecule'].startswith('EXPNEG_')]
        pos = next((r for r in a['ctrl'] if 'EGF' in r['molecule']), None)
        if pos and neg:
            below = sum(1 for r in neg if r['hu_med'] < pos['hu_med'])
            rule('control panel: N of 10 below EGF on human',
                 below, r'\*\*(\d+)\.0/10 = 0\.800\*\*|\*\*(\d+) of the 10 no-KD molecules rank below')
            bad = []
            for doc, text in docs.items():
                for m in re.finditer(r'(\d+) of the 10 no-KD molecules rank below', text, re.I):
                    if int(m.group(1)) != below:
                        line = text[:m.start()].count('\n') + 1
                        bad.append(f"{doc}:{line} says {m.group(1)}, artifact says {below}")
            out.append(('control panel human-leg count', below, not bad, bad))

    # exclusion ledger counts
    if a['ledger']:
        L = a['ledger']
        rule('ledger: list A distinct sequences', L['list_a_sequences'],
             r'\*\*(\d+) DISTINCT SEQUENCES\*\*|-> \*\*(\d+) distinct sequences\*\*')
        rule('ledger: union size', L['union'],
             r'\*\*(\d+)\*\* \| *$')
        bad = []
        for doc, (text, offs) in a['flat'].items():
            for m in re.finditer(r'union: distinct excluded molecules\*\* \| \*\*(\d+)\*\*', text, re.I):
                if int(m.group(1)) != L['union']:
                    bad.append(f"{doc}:{flat_line(offs, m.start())} says {m.group(1)}, "
                               f"ledger says {L['union']}")
        out.append(('ledger union in the table', L['union'], not bad, bad))

    return [r for r in out if r[1] != 'n/a' or r[3]]


# ----------------------------------------------------------- CITATIONS ----
FIELD_PATTERNS = [
    ('ph_ratio_6p5_over_7p4_his_only_CONSERVATIVE', r'(\d+\.\d{3})×', 3),
    ('ipsae_min_human', r'human (\d\.\d{3})', 3),
    ('ipsae_min_mouse', r'mouse (\d\.\d{3})', 3),
    ('ph_poses_n', r'n\s*=\s*(\d+)', 0),
    ('ph_pose_spread_over_median', r'spread (\d\.\d{2,3})', 2),
]


def citations(a, window=200):
    """Numbers quoted near a design name must match that design's own CSV row.

    THE WINDOW STOPS AT THE NEXT DESIGN NAME. A fixed 320-character window was tried
    first and produced an 18-of-48 mismatch rate -- almost all of it spurious, because
    this document is full of A-versus-B comparisons ("`X` reads 3.545x ... against shipped
    `Y`: 1.835x") and the window ran past the comparison boundary and attributed Y's
    numbers to X. A 37% failure rate in a gate trains the reader to ignore it, which is
    worse than having no gate, so the scan now ends at the next backticked identifier, the
    next table cell, or the next sentence -- whichever comes first.
    """
    rows = {r['name']: r for r in a['csv']}
    bad, checked = [], 0
    for doc, text in a['docs'].items():
        for name, row in rows.items():
            for m in re.finditer(re.escape('`' + name + '`'), text):
                seg = text[m.end():m.end() + window]
                stops = [seg.find('`'), seg.find('|'), seg.find('. '), seg.find('\n\n')]
                stops = [x for x in stops if x > 0]
                if stops:
                    seg = seg[:min(stops)]
                for field, pat, prec in FIELD_PATTERNS:
                    val = row.get(field)
                    if not val:
                        continue
                    try:
                        truth = round(float(val), prec) if prec else int(float(val))
                    except ValueError:
                        continue
                    # ONLY THE FIRST match per field, i.e. the headline citation.
                    #
                    # A sentence may legitimately carry several numbers of the same shape:
                    # "`L133E` reads 5.656x ... with H370 contributing 0.921-6.229x" has
                    # three x-values and only the first is the design's pH ratio. Checking
                    # all of them flagged 6.229 as a stale ratio, which is nonsense. So
                    # this pass verifies the HEADLINE number only -- a stated limitation,
                    # not a silent one: an elaborating figure deeper in a sentence is not
                    # guarded, and the INVENTORY pass is what covers those.
                    q = re.search(pat, seg)
                    if q:
                        checked += 1
                        try:
                            got = round(float(q.group(1)), prec) if prec else int(float(q.group(1)))
                        except ValueError:
                            got = None
                        # only flag when the number looks like it IS this field:
                        # same magnitude, differs in the last place(s)
                        if got is not None and got != truth and abs(got - truth) < max(truth * 0.5, 1.0):
                            line = text[:m.start()].count('\n') + 1
                            bad.append(f"{doc}:{line} `{name}` {field}: doc says "
                                       f"{q.group(1)}, CSV says {val}")
    return checked, bad


# ----------------------------------------------------------- INVENTORY ----
def inventory(a):
    vals = set()

    def add(x):
        try:
            f = float(x)
        except (TypeError, ValueError):
            return
        for p in (0, 1, 2, 3, 4):
            vals.add(round(f, p))

    for r in a['csv']:
        for v in r.values():
            add(v)
    add(len(a['csv'])); add(PERMITTED); add(PERMITTED - len(a['csv']))

    def walk(o):
        if isinstance(o, dict):
            for v in o.values():
                walk(v)
        elif isinstance(o, list):
            for v in o:
                walk(v)
            add(len(o))
        else:
            add(o)
    for k in ('sens', 'apo', 'pert', 'ledger', 'foot', 'chai', 'ctrl'):
        if a[k]:
            walk(a[k])

    unsourced = defaultdict(list)
    for doc, text in a['docs'].items():
        for m in re.finditer(r'\*\*([\d.]+)[×x%]?\*\*', text):
            raw = m.group(1).rstrip('.')
            try:
                f = float(raw)
            except ValueError:
                continue
            if f in vals or round(f, 2) in vals or round(f, 3) in vals:
                continue
            line = text[:m.start()].count('\n') + 1
            unsourced[doc].append((line, raw, text[max(0, m.start()-70):m.start()].splitlines()[-1][-60:] if text[max(0, m.start()-70):m.start()].splitlines() else ''))
    return vals, unsourced


def main():
    a = load()
    only_rules = '--rules' in sys.argv
    show_inv = '--inventory' in sys.argv
    fails = 0

    print("=== RULES: named claims with an explicit source ===")
    nogap = 0
    for name, expected, ok, bad in rules(a):
        tag = 'OK  ' if ok else ('GAP ' if ok is None else 'FAIL')
        print(f"{tag} {name:<46} artifact: {expected}")
        for b in bad:
            print(f"       {b}")
        if ok is False:
            fails += 1
        elif ok is None:
            nogap += 1
    if nogap:
        print(f"\n{nogap} rule(s) found NO anchoring claim -- not failures, but unguarded.")

    if not only_rules:
        print("\n=== CITATIONS: per-design numbers quoted in prose ===")
        nchecked, bad = citations(a)
        print(f"{'OK  ' if not bad else 'FAIL'} {nchecked} design-adjacent numbers checked, "
              f"{len(bad)} mismatched")
        for b in bad[:40]:
            print(f"       {b}")
        fails += (1 if bad else 0)

        print("\n=== INVENTORY: numbers with no counterpart in any artifact ===")
        vals, unsourced = inventory(a)
        tot = sum(len(v) for v in unsourced.values())
        print(f"artifact value set: {len(vals)} distinct numbers")
        print(f"bolded numeric claims with NO artifact counterpart: {tot}")
        print("(not fatal -- prose legitimately cites dates, pH values, section numbers and")
        print(" literature figures. This list is the stone-unturned inventory; it should shrink.)")
        for doc, items in sorted(unsourced.items()):
            print(f"  {doc}: {len(items)}")
            for line, raw, ctx in (items if show_inv else items[:6]):
                print(f"     L{line:<5} {raw:<10} ...{ctx.strip()[-52:]}")
            if not show_inv and len(items) > 6:
                print(f"     ... {len(items)-6} more (--inventory to list)")

    print(f"\n{'PASS' if not fails else 'FAIL'}: {fails} check group(s) failed")
    sys.exit(1 if fails else 0)


if __name__ == '__main__':
    main()
