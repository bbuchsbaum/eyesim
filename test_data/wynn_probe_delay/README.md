# Wynn probe-delay GazeWeave inputs

These are frozen, phase-resolved fixation inputs for the provisional
GazeWeave probe-delay persistence court. They were supplied by the repository
owner on 2026-08-20 and moved here from the repository root without changing
their contents.

The files come from a new recognition-memory study with a four-level newness
rating, not the experiment reported in Wynn, Ryan, and Buchsbaum (2020). The
legacy `wynn_probe_delay` directory name is retained to avoid breaking frozen
manifests and validation paths; it must not be interpreted as study provenance.

## Files

- `study_fix_input_new.csv`: 101,991 fixation rows from 46 participants.
  Each studied image was intended to appear four times. The table retains
  participant, trial, presentation, future retrieval condition, image,
  saliency, position, duration, and fixation timing fields.
- `testdelay_fix_input_matched.csv`: 29,699 fixation rows from the same 46
  participants. It retains the 1--4 newness rating and binary accuracy plus
  relative probe onset, inferred epoch source, probe/delay classification, and
  boundary-straddling status. The retrieval trial contains a 500-ms probe
  followed by a 2,500-ms delay.
- `manifest.csv`: immutable acquisition and integrity metadata.

The experimental image rectangle is x = 112--912 and y = 84--684 on a
1024-by-768 display. Court importers translate this rectangle to an
800-by-600 analysis coordinate system and retain the unmodified source
coordinates in these CSV files.

## Handling contract

- The source CSVs are read-only, explicitly ignored by Git, and excluded from
  the installed R package. They are local analysis inputs and must not be
  committed, pushed, or otherwise published.
- No raw participant rows are copied into validation results. Published
  artifacts may contain aggregate metrics, frozen pseudonymous fold IDs, and
  checksums only.
- While this provisional court is under development, its entire results
  directory is also Git-ignored; publishing even derived trial scores requires
  a later explicit decision by the repository owner.
- No separate data-only licence has been inferred. The repository's existing
  code licence does not by itself establish a redistribution licence for the
  participant-level tables.
- Any regenerated or corrected input receives a new filename and checksum;
  these files are never overwritten in place.

The prospective analysis contract is
`inst/validation/GAZEWEAVE-PROBE-DELAY-CONTROLS.md`.
