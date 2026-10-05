# Problem 2 day-1 free arithmetic. Playbook Part I section 1 + 2, re-derived for a REVERSED switch.
# Problem 1: bind at pH 6.5, off at 7.4  -> protons TAKEN UP on binding.
# Problem 2: bind at pH 7.4, off at 6.0  -> protons RELEASED on binding. Opposite sign.
HI, LO = 7.4, 6.0
DPH = HI - LO

def K_factor(pka_free, pka_bound, ph):
    return (1 + 10**(pka_bound - ph)) / (1 + 10**(pka_free - ph))

def ratio(pka_free, pka_bound):
    "K(7.4)/K(6.0) -- how much binding is LOST going acidic. >1 means the switch works."
    return K_factor(pka_free, pka_bound, HI) / K_factor(pka_free, pka_bound, LO)

print(f"WINDOW: {HI} -> {LO} = {DPH:.1f} pH units (problem 1 was 0.9)")
print("ONE-PROTON BOUND (thermodynamic, no model can beat it):")
for n in (1,2,3,4):
    print(f"  n={n} protons -> ceiling {10**(n*DPH):>12,.1f}x     (problem 1 n={n}: {10**(n*0.9):>10,.1f}x)")

print("\n(a) WHAT ONE SITE CAN GIVE. Best case pKa_bound -> -inf (binding forces it neutral).")
print(f"    {'pKa_free':>9}{'best ratio':>12}  note")
for pf in (6.0, 6.5, 6.7, 7.0, 7.4, 8.0, 9.0, 11.0):
    note = "  <- textbook free His" if pf==6.7 else ("  <- free Lys/Tyr territory" if pf>=9.0 else "")
    print(f"    {pf:>9.2f}{ratio(pf,-99):>12.2f}{note}")

print("\n(a2) REALISTIC shifts: a site whose pKa drops by dpKa on binding.")
print(f"    {'pKa_free':>9}" + "".join(f"{f'-{d}':>9}" for d in (1.0,1.5,2.0,2.5,3.0)))
for pf in (6.5, 6.7, 7.0, 7.5, 8.0):
    print(f"    {pf:>9.2f}" + "".join(f"{ratio(pf, pf-d):>9.2f}" for d in (1.0,1.5,2.0,2.5,3.0)))

print("\n(b) WHAT 'NO DETECTABLE BINDING AT pH 6.0' REQUIRES, assay ceiling 10,000 nM.")
CEIL = 10_000.0
print(f"    {'KD(7.4)':>10}{'ratio needed':>14}{'protons':>10}{'  (protons needed in problem 1 window)':>10}")
for kd in (1,10,50,100,500,1000,2000):
    r = CEIL/kd
    import math
    print(f"    {kd:>7} nM{r:>13,.0f}x{math.log10(r)/DPH:>10.2f}{math.log10(r)/0.9:>12.2f}")
print(f"\n(c) CROSSOVER: one proton at the full {10**DPH:.1f}x ceiling switches off anything weaker than {CEIL/10**DPH:,.0f} nM at 7.4.")
print(f"    At a realistic single-His ceiling of ~7x, the crossover is {CEIL/7:,.0f} nM.")
