#!/usr/bin/env python3
"""What does "no detectable binding at 7.4" actually require? Reproduces METHODS section 1.1.

Run it: python3 bin/linkage_limits.py


Wyman linkage: d(ln K)/d(pH) = -dNu(H+), so over a window of dpH the maximum achievable
enhancement is 10^(n*dpH) where n is the number of protons taken up on binding. That bound
is thermodynamic, not a property of any model -- it is the same relation behind the Bohr
effect. Coupling between sites cannot beat it; it can only help you approach it.

So the question "can a single site do it" has two parts:
  (a) what ratio can one site give?            -> bounded by 10^0.9
  (b) what ratio does "no detectable" need?    -> depends on KD(6.5) and the assay ceiling
"""
DPH = 7.4 - 6.5


def ratio(pka_free, pka_bound, lo=6.5, hi=7.4):
    K = lambda ph: (1 + 10 ** (pka_bound - ph)) / (1 + 10 ** (pka_free - ph))
    return K(lo) / K(hi)


print(f"ONE-PROTON BOUND over {DPH} pH units: {10 ** DPH:.3f}x")
for n in (1, 2, 3, 4):
    print(f"  n = {n} protons -> ceiling {10 ** (n * DPH):>9.1f}x")

print("\n(a) WHAT ONE SITE CAN GIVE, as a function of its free pKa.")
print("    Best case is pKa_bound -> +inf, i.e. the site is fully protonated when bound.")
print(f"    {'pKa_free':>9}{'best ratio':>12}   note")
for pf in (4.0, 5.0, 5.5, 6.0, 6.22, 6.5, 7.0, 7.4, 8.0):
    best = ratio(pf, 99.0)
    note = "  <- H433, the site we used" if abs(pf - 6.22) < 0.01 else ""
    print(f"    {pf:>9.2f}{best:>12.3f}{note}")
print("    The ceiling approaches 7.94x only as pKa_free falls well below 6.5: a site that")
print("    is already deprotonated at BOTH pHs when free gives the full proton's worth.")

print("\n(b) WHAT 'NO DETECTABLE BINDING AT 7.4' REQUIRES.")
print("    Detection fails when KD(7.4) exceeds the assay's upper quantifiable limit.")
CEIL = 10_000.0  # nM, a typical SPR/BLI upper limit
print(f"    Taking the limit as {CEIL:.0f} nM:")
print(f"    {'KD(6.5)':>10}{'needed ratio':>14}{'protons needed':>16}")
import math
for kd in (1.0, 10.0, 50.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0):
    need = CEIL / kd
    n = math.log10(need) / DPH
    print(f"    {kd:>8.0f} nM{need:>13.0f}x{n:>15.2f}")

print("\n(c) THE CROSSOVER: the weakest binder that one proton can switch 'off'.")
one = 10 ** DPH
print(f"    One site at its 7.94x ceiling takes KD(6.5) = {CEIL / one:,.0f} nM to {CEIL:,.0f} nM.")
print(f"    So a single site satisfies the stated criterion only if the binder is WEAKER")
print(f"    than about {CEIL / one:,.0f} nM at pH 6.5 -- and at H433's actual free pKa of 6.22,")
print(f"    the ceiling is {ratio(6.22, 99.0):.3f}x, so weaker than {CEIL / ratio(6.22, 99.0):,.0f} nM.")

print("\n(d) WHAT G532 IMPLIES. Its published SPR ratio is 13.26x on human EGFR.")
g = 13.26
print(f"    13.26x exceeds the one-proton bound of {one:.2f}x, so G532 must exchange")
print(f"    more than one proton: n >= {math.log10(g) / DPH:.2f}. It is multi-site by arithmetic,")
print("    which is the direct evidence that the two-site route is reachable.")
