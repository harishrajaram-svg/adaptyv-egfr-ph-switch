# Orientation: necessary, not sufficient -- and a correction

All 8 Mosaic seeds (tune01 + tune02, L76, 180 steps each):

| design | acid_O | acid_CA | O<CA | pass |
|---|---|---|---|---|
| tune01_s0 | 5.20 | 3.46 | no | no |
| tune01_s1 | 4.96 | 7.20 | yes | no |
| tune01_s2 | 3.17 | 5.52 | yes | YES |
| tune01_s3 | 9.91 | 7.55 | no | no |
| tune02_s0 | 15.00 | 15.58 | yes | no |
| tune02_s1 | 3.77 | 6.20 | yes | YES |
| tune02_s2 | 12.49 | 11.60 | no | no |
| tune02_s3 | 2.46 | 4.41 | yes | YES |

| | passed | failed |
|---|---|---|
| oriented (O<CA) | **3** | **2** |
| flipped (CA<O) | **0** | 3 |

## Correction

I said in conversation that this held "eight for eight -- every pass has O closer
than CA, every failure has CA closer." **The second half is false.** `tune01_s1`
(4.96 A) and `tune02_s0` (15.00 A) are both correctly oriented and both failed.
The honest claim is one-directional:

- All 3 passes are oriented and 0 flipped designs ever passed, so correct
  orientation looks **necessary**.
- 2 of 5 oriented designs still failed, so it is **not sufficient**.
  `tune02_s0` is oriented and 15 A away -- orientation with no proximity at all.

## What this means for the directional term

It still justifies building it: nothing flipped has ever passed, and flipped
orientation accounts for 3 of 5 failures. But it will not be enough alone,
because distance fails independently. The directional reward must be **added
alongside** the distance term, not substituted for it -- which is how it was
specified, and this is the evidence for why.

n=8, one length, two configurations. A design hint, not a measured law.
