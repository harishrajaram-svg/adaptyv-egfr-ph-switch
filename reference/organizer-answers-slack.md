# Organizer answers from the Proteinbase Slack — read 2026-10-04 11:25 AM

## 1. THE DEADLINE MOVED. Oct 6 23:59 AoE, not Oct 4.
Tudor-Stefan Cotet (Adaptyv), Oct 2, posted to both #general and
#anthropic_adaptyv_competition:
  "We are extending the deadline for the first challenge until **Oct 6 23:59 AoE**. This is
   due to some teams receiving their Modal credits a bit later. We will update all submission
   deadlines mentioned on the first challenge's page to match this."
AoE = UTC-12, so Oct 6 23:59 AoE = Oct 7 11:59 UTC = **Oct 7, 7:59 AM EDT**.
From 11:25 AM Oct 4 that is **~68.5 hours**, not the ~21 we were working to.
(An earlier 1-day extension was announced Sep 29; this supersedes it.)

## 2. NOVELTY IS A HARD AUTOMATIC GATE AT UPLOAD. We pass; two new candidates do not.
Submission rules as quoted on the site: "Designs need a novelty score 3/4 or higher
(checked automatically after upload)." Simon (Adaptyv): "you need to clear novelty 3 ...
you can use an existing scaffold but **the sequence should be <30% similar to anything
existing**."
CHECKED, 16 of the 20 submitted designs (the 4 Mosaic ones have no generator pose here):
best sequence identity ranges **0.105 to 0.197** -- all comfortably inside the 30% bar.
They fail only the STRUCTURAL half (TM >= 0.5), which is Level 4, not the Level 3 gate.
**The current submission is safe.**
BUT the two gc candidates found this morning are NOT:
    `h370_2site__boltzgen_egfr_2site_009`  sequence identity **0.745**  -> auto-rejected
    `h370_only__boltzgen_egfr_h370_020`    sequence identity **0.413**  -> auto-rejected
Neither can be submitted. `h370_020` was the one carrying the best mouse number in the
project (0.7072). Caught before promotion, by one gate.

## 3. MOUSE IS PROBABLY MEASURED ONLY AT pH 7.4 -- which puts the two criteria in conflict.
Amir Shanehsazzadeh (Anthropic), Oct 3: "I'm pretty sure it's **human EGFR at both pH and
mouse EGFR at neutral (7.4)** but Tudor can confirm."
Olga Lavinda (NYU) made the consequence explicit on Oct 4 at 10:22 AM, still unanswered:
  "cross-reactivity is defined as a KD ratio mouse/human of around 1 ... If mouse is measured
   only at 7.4, **a design that genuinely switches off at 7.4 cannot show a mouse/human ratio
   near 1 — the two headline criteria work against each other for any pH-switch design.**"
CONSEQUENCE FOR US: our tier-2 ordering sorts by mouse affinity descending. If mouse is read
only at 7.4, a high mouse number means the design still binds at 7.4, i.e. it switches LESS.
Our ranking may be rewarding the opposite of what we want. This needs a decision.

## 4. There is no minimum RATIO. The minimum AFFINITY question is still unanswered.
Amir, Oct 3, answering "is there a minimum ratio": "I'd say the ratio of affinities is itself
a criteria here. **I don't know that we have a minimum ratio.**"
ingmar (Provolut) asked the sharper question on Oct 2 and it has NOT been answered:
"Is there a minimum affinity at pH 6.5 for a design to count as pH-selective? ... with at most
~10x per proton over this pH interval and a 10 uM quantifiable ceiling, a design needs two
switching histidines to look on/off at micromolar affinity and three at 100 nM."
That is the same ceiling arithmetic we derived independently, from another team.

## 5. The assay, which we have not been designing against
  * target: **HEK293-expressed, fully human-glycosylated**, **tethered (inactive)** EGFR ECD.
    "we'll screen against the glycosylated, tethered EGFR, so you need to account for that
    when designing (e.g., exposed epitopes, that's why we suggested Domain III)."
  * **MES replaces HEPES** for the acidic condition. Buffer otherwise 10 mM HEPES / 150 mM
    NaCl / 0.2% Tween-20 / 3 mM EDTA, as in the Germinal paper.
  * BLI and SPR, in vitro. **Binders expressed by E. coli cell-free synthesis.**
  * Target is His-tagged. Amir: "The tag or really any component of the system that is there
    for screening should definitely not be targeted. Good call out given it's Histidine. Any
    such binder might look pH-selective but it would bind to anything with a His tag."
    -> our designs engage H433 on domain III, so this is fine, but it is worth stating in the
    methods that we checked.
  * RISK WE HAVE NOT ADDRESSED: full human glycans. We designed against an unglycosylated
    crystal structure. Domain III carries sequons. A design whose epitope overlaps a glycan
    may simply not bind.

## 6. Track 3 is judged mostly on METHODS, and in-silico scoring counts for little
Tudor: "**In silico scoring will only be a smaller component in the overall selection
criteria.** ... worse confidence scores on the rigid structure or on obstructed epitopes
won't be penalized much (if at all)."
So the methods document is not a formality -- it is the main deliverable.

## 7. Mechanics
  * **Max 20 sequences per team/Proteinbase account in a single collection.**
  * **You can submit once every 24 hours.** So an early safe submission can be improved later.
  * Whether a later submission replaces or adds to an earlier one was asked Oct 2 and is
    not answered in the channel.

## 8. The organizers themselves signposted Proton-PottsMPNN
Amir, Oct 3, @channel: "I'd like to make sure everyone is aware of **Proton-PottsMPNN** ...
which came out just this past week and is highly relevant for pH-sensitive binder design.
**The authors expedited the release to be available in time for the competition!**"
Code: https://github.com/christian-creator/ProtonPottsMPNN
With 68 hours rather than 21, this is now actually runnable, and it is the method built for
exactly this problem, endorsed in-channel by the organizer.
