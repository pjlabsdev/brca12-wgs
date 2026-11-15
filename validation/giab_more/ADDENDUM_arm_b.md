# Addendum to the preregistration: code version of arm B

Frozen on 2026-11-15, before any whole-genome run of HG003 to HG007 (arm B) was started.

Arm A was run on `main` ce825de and missed the SNV threshold in two samples. One cause, poly-G
read tails in soft clips, was fixed afterwards (`main` a705662). Arm B therefore runs on a705662,
not on ce825de. Samples, truth sets (v4.2.1 for all five), comparison script and thresholds stay
exactly as in `PREREGISTRATION.md`.

What arm B can and cannot show: its reads (NovaSeq PCR-free 30x, Google) are different reads of
the same five people whose GIAB HiSeq reads were used to choose the fix. They were never run
through the pipeline, so read-level artefacts are independent; the true variants are the same.
A pass in arm B is therefore reported as "same people, new reads", not as new people.

Each sample is reported whether it passes or not; nothing is changed after the first arm B result.
