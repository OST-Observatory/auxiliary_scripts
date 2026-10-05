#! /usr/bin/python
# -*- coding: utf-8 -*-

"""
    Archive pipeline, step 1 of 3: fetch the frames of an object or an
    observation run from the OST data archive.

    Downloaded:
      * the science frames of the request (all targets of a run, or all
        frames of an object over all its runs),
      * every bias / dark / flat candidate of these runs and of the
        neighbouring runs within ``calib_window_days``.
    Recorded as metadata only (not downloaded): other light frames in that
    window; they show when another camera was mounted at the telescope.

    Files land in a content-addressed cache (``cache_dir``) and are never
    downloaded twice. The result is ``<workspace>/manifest.ecsv``, the input
    of ``2_classify_and_group.py``.

    Credentials: set OST_ARCHIVE_USER / OST_ARCHIVE_PASSWORD in the
    environment, or answer the prompt. Anonymous access works for public
    runs but is limited to 60 requests per minute.

    Instead of the archive, a local directory tree can be used
    (``local_directory``); the archive is then not contacted.

    Needs ost_photometry >= 0.6.
"""

############################################################################
#                                 Request                                  #
############################################################################
#   Either an object (all its runs) ...
object_name: str | None = None
# object_name: str | None = "M57"

#   ... or an observation run (all targets of that night).
run_name: str | None = "?"
# run_name: str | None = "2022-03-08"

#   Run mode only: restrict the science frames to these targets
#   (archive names, case and spaces ignored). ``None`` = all targets.
targets: list[str] | None = None
# targets: list[str] | None = ["M57", "M104"]

#   Alternative: a local directory tree with FITS files (no archive access).
local_directory: str | None = None

############################################################################
#                Additional options: only edit if necessary                #
############################################################################
#   Days around the run(s) searched for calibration frames.
calib_window_days: float = 7.0

#   Anonymous access (public runs only, slow rate limit).
anonymous: bool = False

#   Science exposures without matching darks (same camera, binning, gain,
#   offset, readout mode, exposure time, temperature +-2 K) or setups
#   without bias in the window: ask the archive's dark / bias finder and
#   download the frames of the closest run. Needs a login (not anonymous).
use_dark_finder: bool = True
finder_kinds: list[str] = ["dark", "bias"]

#   Darks match lights within max(dark_exptime_tolerance s,
#   fraction x exposure time) -- as in 2_classify_and_group.py.
dark_exptime_tolerance: float = 0.5
dark_exptime_tolerance_fraction: float = 0.05

#   Download also the context frames (other lights in the window). Only
#   their metadata is needed for the grouping.
download_context: bool = False

#   Archive URL.
archive_url: str = "https://polaris.astro.physik.uni-potsdam.de/data_archive"

#   File cache shared by all requests, and the working directory of this
#   request.
cache_dir: str = "~/.cache/ost_archive"
workspace: str = "archive_workspace/"

############################################################################
#                               Libraries                                  #
############################################################################

import time
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

from astropy import log

log.setLevel("ERROR")

from ost_photometry import style
from ost_photometry.archive import (
    ArchiveClient,
    fetch_dataset,
    manifest_from_directory,
    write_manifest,
)

############################################################################
#                                  Main                                    #
############################################################################


def _progress(i: int, n: int, name: str) -> None:
    if i == 1 or i == n or i % 25 == 0:
        print(f"   [{i}/{n}] {name}", flush=True)


if __name__ == "__main__":
    start_time = time.time()
    work = Path(workspace).expanduser()
    work.mkdir(parents=True, exist_ok=True)

    if local_directory:
        manifest = manifest_from_directory(Path(local_directory).expanduser())
        report_lines = [f"Local directory {local_directory}: {len(manifest)} FITS frames"]
    else:
        if (run_name in (None, "?")) == (object_name in (None, "?")):
            raise SystemExit(
                f"{style.Bcolors.FAIL}Set exactly one of object_name or run_name "
                f"(or local_directory).{style.Bcolors.ENDC}"
            )
        client = ArchiveClient(archive_url)
        if not anonymous:
            client.login()
        manifest, report = fetch_dataset(
            client,
            object_name=object_name if object_name not in (None, "?") else None,
            run_name=run_name if run_name not in (None, "?") else None,
            targets=targets,
            calib_window_days=calib_window_days,
            cache_dir=Path(cache_dir).expanduser(),
            download_context=download_context,
            use_dark_finder=use_dark_finder,
            finder_kinds=finder_kinds,
            exptime_tolerance=dark_exptime_tolerance,
            exptime_tolerance_fraction=dark_exptime_tolerance_fraction,
            progress=_progress,
        )
        report_lines = report.lines()
        client.logout()

    path = write_manifest(manifest, work / "manifest.ecsv")
    (work / "fetch_report.txt").write_text("\n".join(report_lines) + "\n")
    for line in report_lines:
        print("   " + line)
    print(style.Bcolors.OKGREEN + "   Done" + style.Bcolors.ENDC)
    print(f"   Manifest: {path}")
    print("   Next: 2_classify_and_group.py")
    print("--- %s minutes ---" % ((time.time() - start_time) / 60.0))
