# Auxiliary scripts and programs
Various scripts and programs dedicated to astronomical observations and their post-processing.

`compare_registration/` — run every `ost_photometry.reduce` image-shift method on two frames (or raw datasets) and plot the differences.

`astro_imaging/` — reduce, select the best frames (FWHM / roundness / best X %), register all filters onto one grid and write weighted linear stacks per filter; re-stack with another selection; standalone frame-quality CLI. Replaces the Siril register → filter → weighted-stack loop.

`archive_pipeline/` — fetch an object or observation run from the OST data archive, classify frames, group calibration data by mount session (camera orientation from plate solving) with flat probabilities, reduce per group and stack every target.
