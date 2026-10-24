# ORCHID precision check, supplement: reviewer instructions

**Why.** ORCHID releases data read from papers and handbooks by machine. To report how often what we release is right,
we ask people to check a random sample by hand. Your answers are the measurement, so answer what you see, not what you
expect. "No" and "can't tell" are as useful to us as "yes".

**What you get.** One page (your packet) with 34 items: 30 of your own (10 phase-registry entries R..., 10 ESPEI-ready
values E..., 10 assessed phase diagrams A...) and 4 that another reviewer also checks, so that we can measure agreement
(some of your own items are in another packet too). Please do not discuss items with other reviewers until everyone has
handed in.

**How to work.** Open the page in Chrome or Firefox (double-click the file, or the localhost link you were given).
Type your name at the top. Answers save in this browser as you go. Use **Next unanswered** to move on. When you finish,
press **Download results** and email the JSON file to the coordinator. Download a backup now and then too.
The page contains images of copyrighted papers and of the ASM handbook: do not upload, post or forward it.

## Answers

- **yes**: you found it in the source and it is right.
- **no**: it is wrong. Say what is wrong in the note (one line is enough, e.g. "page says 1053 K, not 1035 K").
- **can't tell**: you cannot decide (page unreadable, paper not accessible, value not on the page shown).
- **n/a**: the question does not apply (no structure asserted, no conversion, no printed uncertainty).

## R items: phase registry

- **Name and aliases**: do all names listed (TDB names, paper labels) denote this one phase of this system?
- **Structure**: only when a structure is asserted. Is the prototype, Strukturbericht, Pearson symbol and space group
  right for this phase? Use the source rows shown, or a standard reference (Pearson's, ICSD, the paper). Otherwise n/a.
- **TDB definition**: open the TDB file (link) at the line given. Are the PHASE and CONSTITUENT statements copied exactly?
- **Example papers**: open one or two DOIs. Do they discuss this phase (under any of its names)?

## E items: ESPEI-ready values

The page image is the page of the paper where the value was read (table, row and caption are given; click to enlarge).

- **Value**: is the ESPEI value the printed value, after the conversion stated on the item (for example kJ to J,
  °C to K, eV/atom to J/mol-atom, or a change of reference phase such as bcc to hcp Zr)?
- **Unit / basis**: per mole of atoms or per mole of formula, J vs kJ, activity vs activity coefficient.
- **Conditions**: temperature, pressure (1 atm is assumed when none is printed) and composition (mole fraction of the
  element named, e.g. X_ZR).
- **Phases** and **reference state**: are the phases and the reference states the ones the paper states?
- **Conversion**: if the item lists a conversion, is it done right? Otherwise n/a.
- **Uncertainty**: if the paper prints an uncertainty (±), is `printed` that number, unchanged? Are the other
  components plausible (reading = last printed digit; method = a stated prior)? n/a when nothing is printed and the
  components are only priors.

## A items: assessed phase diagrams (blind)

As in the earlier expert sample, you do not see ORCHID's answer. Name the phases present in the **orange-tinted field**
(the black ring is the released point inside it), from the drawing alone. Build the name with the phase buttons or type
it, e.g. `L+(Cu)`. Use `L#2` for a second liquid and `G` for gas. Choose **cannot tell** if the drawing does not let
you decide, **not one field** if the tint covers more than one field, and **undetermined in the source** if the
handbook itself leaves it open (a "?" or a line that stops). Spelling and order of phases do not matter. Set your
confidence.

**Time.** A packet takes about 2 to 3 hours: an R item about 5 minutes, an E item 5 to 8 minutes (finding the row on
the page), an A item 1 to 3 minutes. You can stop and come back: answers stay in the browser you used. Please hand in by
**27 October 2026**. Questions: ask the coordinator.
