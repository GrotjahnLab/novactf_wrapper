#!/usr/bin/env python3
"""
novaCTF wrapper -- implements the full pipeline documented in readme.md:

  1. Parse taSolution.log to find which views survived alignment, and use
     them as the -secs range for newstack (aligned + binned 2D stack).
  2. Copy the Warp .xml file in and parse it directly (Dose, Angles,
     GridCTF) to generate the .tlt, .defocus, and dose .txt files needed
     for the next stages -- restricted to the same views selected in step 1.
  3. Validate that the stack size matches the tlt/defocus files.
  4. novaCTF -Algorithm defocus            -> generate shifted defocus files.
  5. novaCTF -Algorithm ctfCorrection       -> CTF-correct every defocus-
     shifted stack, in parallel across --cores.
  6. clip flipyz + novaCTF -Algorithm filterProjections, also in parallel.
  7. novaCTF -Algorithm 3dctf               -> WBP reconstruction with 3D-CTF
     correction.
  8. trimvol -yz                            -> reorient back to XYZ.
  9. binvol -binning <final/stack>          -> produce the final binned
     tomogram.
 10. optional cleanup of intermediate files.

Every step logs its command (with a timestamp) to novaCTF_process.log before
running it, and every step is resumable: if its expected output already
exists, it is skipped (use --force to redo it anyway). This must be run
under WSL/Linux -- novaCTF and IMOD both require it.

Example:
    python3 novactf_pipeline.py \\
        --name MIM019_2_lam1_ts_002 \\
        --ta-solution-log taSolution.log \\
        --xml-dir ../.. \\
        --tomo-size 4096,5760,2400 --bin-factor 2 \\
        --pixel-size 1.66 --defocus-step 20 \\
        --final-bin 6 --cores 4
"""

import concurrent.futures
import glob
import os
import re
import subprocess
import sys
import threading
import xml.etree.ElementTree as ET
from datetime import datetime
from pathlib import Path

import click
import numpy as np

_log_lock = threading.Lock()


class NovaCTFError(RuntimeError):
    """Raised whenever a step fails, with a message meant to be read by a human."""


# --------------------------------------------------------------------------
# logging + subprocess helpers
# --------------------------------------------------------------------------

def log_line(log_path: Path, text: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with _log_lock:
        with open(log_path, "a") as fh:
            fh.write(f"[{timestamp}] {text}\n")


def run_cmd(cmd, log_path: Path, dry_run: bool = False):
    cmd_str = " ".join(str(c) for c in cmd)
    log_line(log_path, cmd_str)
    click.echo(f"+ {cmd_str}")

    if dry_run:
        return None

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        log_line(log_path, f"FAILED (exit {result.returncode}): {cmd_str}")
        raise NovaCTFError(
            f"Command failed (exit {result.returncode}): {cmd_str}\n"
            f"--- stderr ---\n{result.stderr.strip()}\n"
            f"--- stdout (tail) ---\n{result.stdout.strip()[-2000:]}"
        )
    return result


def skip_if_exists(path: Path, force: bool) -> bool:
    """Returns True if `path` already exists and we should skip regenerating it."""
    if path.exists() and not force:
        click.echo(f"-> {path} already exists, skipping (use --force to redo).")
        return True
    return False


# --------------------------------------------------------------------------
# taSolution.log parsing
# --------------------------------------------------------------------------

_HEADER_RE = re.compile(r"^\s*view\s+rotation\s+tilt\s+deltilt\s+mag", re.IGNORECASE)


def compact_ranges(views) -> str:
    """[4,5,6,...,36] -> '4-36'; [4,5,7,8] -> '4-5,7-8'"""
    views = sorted(set(views))
    ranges = []
    start = prev = views[0]
    for v in views[1:]:
        if v == prev + 1:
            prev = v
            continue
        ranges.append((start, prev))
        start = prev = v
    ranges.append((start, prev))
    return ",".join(str(a) if a == b else f"{a}-{b}" for a, b in ranges)


def expand_secs(secs: str) -> list:
    """Inverse of compact_ranges: '4-9,11-36' -> [4,5,...,9,11,...,36] (1-based, ordered)."""
    indices = []
    for part in secs.split(","):
        if "-" in part:
            a, b = part.split("-")
            indices.extend(range(int(a), int(b) + 1))
        else:
            indices.append(int(part))
    return indices


def parse_ta_solution_sections(log_path: Path) -> str:
    """
    Parse taSolution.log for the 'view rotation tilt deltilt mag ...' solution
    table and return the surviving view numbers as a newstack -secs string.
    Uses the LAST such table in the file (the final alignment solution).
    """
    lines = log_path.read_text().splitlines()

    last_header_index = None
    for i, line in enumerate(lines):
        if _HEADER_RE.match(line):
            last_header_index = i

    if last_header_index is None:
        raise NovaCTFError(
            f"Could not find a 'view rotation tilt deltilt mag ...' table in {log_path}"
        )

    views = []
    for line in lines[last_header_index + 1:]:
        tokens = line.split()
        if not tokens:
            break
        try:
            views.append(int(tokens[0]))
        except ValueError:
            break

    if not views:
        raise NovaCTFError(
            f"Found the solution table header in {log_path} but no view rows followed it"
        )

    secs = compact_ranges(views)
    click.echo(f"-> taSolution.log: {len(views)} views survived alignment ({secs})")
    return secs


# --------------------------------------------------------------------------
# IMOD header parsing (for validation)
# --------------------------------------------------------------------------

_DIMS_RE = re.compile(
    r"Number of columns, rows, sections\s*\.*\s*(\d+)\s+(\d+)\s+(\d+)"
)


def get_stack_dimensions(mrc_path: Path):
    result = subprocess.run(["header", str(mrc_path)], capture_output=True, text=True)
    if result.returncode != 0:
        raise NovaCTFError(f"Could not read header of {mrc_path}:\n{result.stderr}")
    m = _DIMS_RE.search(result.stdout)
    if not m:
        raise NovaCTFError(f"Could not parse dimensions out of header of {mrc_path}")
    return tuple(int(x) for x in m.groups())  # (nx, ny, nz)


# --------------------------------------------------------------------------
# pipeline steps
# --------------------------------------------------------------------------

def step_newstack(cfg, log_path):
    if skip_if_exists(cfg.nova_stack, cfg.force):
        return

    cmd = [
        "newstack",
        "-secs", cfg.secs,
        "-fromone",
        "-InputFile", f"{cfg.name}.mrc.st",
        "-OutputFile", str(cfg.nova_stack),
        "-TransformFile", f"{cfg.name}.mrc.xf",
        "-TaperAtFill", "1,0",
        "-TaperAtFill", "1,0",
        "-SizeToOutputInXandY", cfg.size_xy,
        "-OffsetsInXandY", "0.0,0.0",
        "-ImagesAreBinned", "1.0",
        "-BinByFactor", str(cfg.bin_factor),
        "-AntialiasFilter", "-1",
    ]
    run_cmd(cmd, log_path, cfg.dry_run)


def step_xml_and_tlt(cfg, log_path):
    """
    Copy the Warp .xml for this tilt-series in and parse it directly --
    pulling Dose, Angles (tilt), and GridCTF (defocus) values out -- to write
    the .tlt, .defocus, and dose .txt files novaCTF needs. Restricted to the
    same views selected by cfg.secs (taSolution.log / --secs), renumbered
    1..N in stack order, so line N always corresponds to projection N in the
    aligned stack.
    """
    if skip_if_exists(cfg.nova_tlt, cfg.force) and skip_if_exists(cfg.nova_defocus_base, cfg.force):
        return

    xml_src = cfg.xml_dir / f"{cfg.name}.mrc.xml"
    xml_dst = Path(f"{cfg.name}.mrc.xml")
    if xml_src.resolve() != xml_dst.resolve():
        run_cmd(["cp", str(xml_src), str(xml_dst)], log_path, cfg.dry_run)

    log_line(log_path, f"parse {xml_dst} -> {cfg.nova_tlt}, {cfg.nova_defocus_base}, {cfg.nova_dose}")
    click.echo(f"+ parsing {xml_dst} for tlt/defocus/dose (views {cfg.secs})")

    if cfg.dry_run:
        return

    selected = expand_secs(cfg.secs)

    tree = ET.parse(xml_dst)
    root = tree.getroot()

    dose_node = root.find("Dose")
    angles_node = root.find("Angles")
    defocus_node = root.find("GridCTF")
    if dose_node is None or angles_node is None or defocus_node is None:
        raise NovaCTFError(
            f"{xml_dst} is missing one of the expected <Dose>/<Angles>/<GridCTF> elements"
        )

    dose_all = [float(v) for v in dose_node.text.split("\n") if v.strip()]
    tlt_all = [float(v) for v in angles_node.text.split("\n") if v.strip()]
    defocus_all = [int(round(float(node.get("Value")) * 10000)) for node in defocus_node.findall("Node")]

    for label, values in (("Dose", dose_all), ("Angles", tlt_all), ("GridCTF", defocus_all)):
        if not values:
            raise NovaCTFError(f"Found no {label} values in {xml_dst}")

    max_index = max(selected)
    for label, values in (("Dose", dose_all), ("Angles", tlt_all), ("GridCTF", defocus_all)):
        if max_index > len(values):
            raise NovaCTFError(
                f"--secs/taSolution.log selects view {max_index}, but {xml_dst}'s {label} "
                f"array only has {len(values)} entries."
            )

    dose_sel = [dose_all[i - 1] for i in selected]
    tlt_sel = [tlt_all[i - 1] for i in selected]
    defocus_sel = [defocus_all[i - 1] for i in selected]

    np.savetxt(cfg.nova_tlt, tlt_sel, fmt="%f")
    np.savetxt(cfg.nova_dose, dose_sel, fmt="%f")

    lines = []
    for i, (tlt, defocus) in enumerate(zip(tlt_sel, defocus_sel)):
        line = f"{i + 1}\t{i + 1}\t{tlt}\t{tlt}\t{defocus}"
        if i == 0:
            line += "\t2"
        lines.append(line)
    Path(cfg.nova_defocus_base).write_text("\n".join(lines))

    click.echo(
        f"-> Wrote {len(selected)} lines to {cfg.nova_tlt}, {cfg.nova_defocus_base}, {cfg.nova_dose}"
    )


def step_validate(cfg, log_path):
    if cfg.dry_run:
        click.echo("(dry-run) skipping validation")
        return

    nx, ny, nz = get_stack_dimensions(cfg.nova_stack)
    tlt_lines = [l for l in cfg.nova_tlt.read_text().splitlines() if l.strip()]

    if len(tlt_lines) != nz:
        raise NovaCTFError(
            f"Mismatch: {cfg.nova_stack} has {nz} sections but {cfg.nova_tlt} has "
            f"{len(tlt_lines)} tilt angles. Check the -secs range against taSolution.log "
            f"and make sure the tlt file corresponds to the same views."
        )

    click.echo(f"-> Validated: {nz} sections in stack == {len(tlt_lines)} lines in tlt file.")
    log_line(log_path, f"Validated stack/tlt match: {nz} sections")


def step_defocus(cfg, log_path):
    marker = Path(f"{cfg.nova_defocus_base}_0")
    if skip_if_exists(marker, cfg.force):
        return

    cmd = [
        "novaCTF", "-Algorithm", "defocus",
        "-InputProjections", str(cfg.nova_stack),
        "-FULLIMAGE", cfg.size_xy,
        "-THICKNESS", str(cfg.thickness),
        "-TILTFILE", str(cfg.nova_tlt),
        "-SHIFT", "0.0,0.0",
        "-CorrectionType", cfg.correction_type,
        "-DefocusFileFormat", cfg.defocus_format,
        "-CorrectAstigmatism", "1" if cfg.correct_astigmatism else "0",
        "-DefocusFile", cfg.nova_defocus_base,
        "-PixelSize", str(cfg.pixel_size),
        "-DefocusStep", str(cfg.defocus_step),
    ]
    run_cmd(cmd, log_path, cfg.dry_run)


def count_defocus_stacks(cfg) -> int:
    n = len(glob.glob(f"{cfg.nova_defocus_base}_*"))
    if n == 0:
        if cfg.dry_run:
            click.echo("(dry-run) no defocus files yet -- assuming 1 for planning purposes")
            return 1
        raise NovaCTFError(
            f"No files matching {cfg.nova_defocus_base}_* -- did the defocus step run?"
        )
    return n


def _ctf_correction_one(cfg, log_path, i):
    corrected = Path(f"{cfg.nova_base}_corrected.mrc_{i}")
    if skip_if_exists(corrected, cfg.force):
        return
    cmd = [
        "novaCTF", "-Algorithm", "ctfCorrection",
        "-InputProjections", str(cfg.nova_stack),
        "-OutputFile", str(corrected),
        "-DefocusFile", f"{cfg.nova_defocus_base}_{i}",
        "-TILTFILE", str(cfg.nova_tlt),
        "-CorrectionType", cfg.correction_type,
        "-DefocusFileFormat", cfg.defocus_format,
        "-CorrectAstigmatism", "1" if cfg.correct_astigmatism else "0",
        "-PixelSize", str(cfg.pixel_size),
        "-AmplitudeContrast", str(cfg.amplitude_contrast),
        "-Cs", str(cfg.cs),
        "-Volt", str(cfg.volt),
    ]
    run_cmd(cmd, log_path, cfg.dry_run)


def step_ctf_correction(cfg, log_path, n):
    run_in_parallel(cfg, "CTF correction", n, lambda i: _ctf_correction_one(cfg, log_path, i))


def _flip_filter_one(cfg, log_path, i):
    corrected = Path(f"{cfg.nova_base}_corrected.mrc_{i}")
    flipped = Path(f"{cfg.nova_base}_corrected_flipped.mrc_{i}")
    filtered = Path(f"{cfg.nova_base}_filtered.mrc_{i}")

    if skip_if_exists(filtered, cfg.force):
        return

    if not skip_if_exists(flipped, cfg.force):
        run_cmd(["clip", "flipyz", str(corrected), str(flipped)], log_path, cfg.dry_run)

    cmd = [
        "novaCTF", "-Algorithm", "filterProjections",
        "-InputProjections", str(flipped),
        "-OutputFile", str(filtered),
        "-TILTFILE", str(cfg.nova_tlt),
        "-StackOrientation", "xz",
    ]
    run_cmd(cmd, log_path, cfg.dry_run)


def step_flip_filter(cfg, log_path, n):
    run_in_parallel(cfg, "flip + filter", n, lambda i: _flip_filter_one(cfg, log_path, i))


def run_in_parallel(cfg, label, n, task_fn):
    click.echo(f"Running {label} for {n} stacks across {cfg.cores} core(s)...")
    errors = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=cfg.cores) as pool:
        futures = {pool.submit(task_fn, i): i for i in range(n)}
        for future in concurrent.futures.as_completed(futures):
            i = futures[future]
            try:
                future.result()
            except Exception as exc:  # noqa: BLE001 - want to report and continue
                errors.append((i, exc))

    if errors:
        details = "\n".join(f"  stack {i}: {exc}" for i, exc in errors)
        raise NovaCTFError(f"{label} failed for {len(errors)} stack(s):\n{details}")


def step_reconstruction(cfg, log_path):
    if skip_if_exists(cfg.nova_rec, cfg.force):
        return
    cmd = [
        "novaCTF", "-Algorithm", "3dctf",
        "-InputProjections", f"{cfg.nova_base}_filtered.mrc",
        "-OutputFile", str(cfg.nova_rec),
        "-TILTFILE", str(cfg.nova_tlt),
        "-THICKNESS", str(cfg.thickness),
        "-FULLIMAGE", cfg.size_xy,
        "-SHIFT", cfg.shift,
        "-PixelSize", str(cfg.pixel_size),
        "-DefocusStep", str(cfg.defocus_step),
    ]
    run_cmd(cmd, log_path, cfg.dry_run)


def step_trim(cfg, log_path):
    if skip_if_exists(cfg.nova_trim, cfg.force):
        return
    run_cmd(["trimvol", "-yz", str(cfg.nova_rec), str(cfg.nova_trim)], log_path, cfg.dry_run)


def step_bin(cfg, log_path):
    if skip_if_exists(cfg.nova_bin, cfg.force):
        return

    if cfg.final_bin % cfg.bin_factor != 0:
        raise NovaCTFError(
            f"--final-bin ({cfg.final_bin}) must be an exact multiple of --bin-factor "
            f"({cfg.bin_factor}); got a non-integer ratio."
        )
    ratio = cfg.final_bin // cfg.bin_factor
    run_cmd(
        ["binvol", "-binning", str(ratio), str(cfg.nova_trim), str(cfg.nova_bin)],
        log_path,
        cfg.dry_run,
    )


def step_cleanup(cfg, log_path, n):
    if not cfg.cleanup:
        return

    click.echo("Cleaning up intermediate files...")
    patterns = []
    for i in range(n):
        patterns.append(f"{cfg.nova_base}_corrected.mrc_{i}")
        patterns.append(f"{cfg.nova_base}_corrected_flipped.mrc_{i}")
        patterns.append(f"{cfg.nova_base}_filtered.mrc_{i}")
    patterns.append(str(cfg.nova_rec))

    for p in patterns:
        path = Path(p)
        if path.exists():
            log_line(log_path, f"rm {path}")
            click.echo(f"- rm {path}")
            if not cfg.dry_run:
                path.unlink()


# --------------------------------------------------------------------------
# config object
# --------------------------------------------------------------------------

class Config:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

        self.name = kwargs["name"]
        self.nova_base = f"{self.name}.mrc_nova"
        self.nova_stack = Path(f"{self.nova_base}.mrc")
        self.nova_tlt = Path(f"{self.nova_base}.tlt")
        self.nova_defocus_base = f"{self.nova_base}.defocus"
        self.nova_dose = Path(f"{self.nova_base}.dose.txt")
        self.nova_rec = Path(f"{self.nova_base}.rec")
        self.nova_trim = Path(f"{self.nova_base}_trim.rec")
        self.nova_bin = Path(f"{self.nova_base}_bin{self.final_bin}.rec")


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

@click.command(context_settings={"help_option_names": ["-h", "--help"]})
@click.option("--name", required=True,
              help="Tilt-series base name, e.g. MIM019_2_lam1_ts_002 "
                   "(expects <name>.mrc.st, <name>.mrc.xf, <name>.mrc.xml to exist).")
@click.option("--workdir", default=".", show_default=True,
              help="Directory containing the input files; the script chdir's here first.")
@click.option("--secs", default=None,
              help="Manual newstack -secs range (e.g. '4-36'). If omitted, "
                   "parsed automatically from --ta-solution-log.")
@click.option("--ta-solution-log", type=click.Path(path_type=Path), default=Path("taSolution.log"),
              show_default=True, help="Path to taSolution.log, used to auto-detect --secs.")
@click.option("--xml-dir", type=click.Path(path_type=Path), default=None,
              help="Directory containing <name>.mrc.xml (the Warp xml file).")
@click.option("--tomo-size", required=True,
              help="UNBINNED tomogram X,Y,thickness in px, e.g. '4096,5760,2400'. Divided by "
                   "--bin-factor to get SizeToOutputInXandY/FULLIMAGE and THICKNESS.")
@click.option("--bin-factor", type=int, required=True, help="BinByFactor used in newstack.")
@click.option("--shift", default="0.0,0.0", show_default=True, help="SHIFT for reconstruction.")
@click.option("--correction-type", default="phaseflip", show_default=True,
              type=click.Choice(["phaseflip", "multiplication"]))
@click.option("--defocus-format", default="imod", show_default=True,
              type=click.Choice(["imod", "ctffind4", "gctf"]))
@click.option("--correct-astigmatism/--no-correct-astigmatism", default=True, show_default=True)
@click.option("--pixel-size", type=float, required=True,
              help="UNBINNED pixel size in Angstrom, e.g. 1.66. Converted to binned nm "
                   "internally (pixel_size * bin_factor / 10) for the novaCTF -PixelSize flag.")
@click.option("--defocus-step", type=float, required=True, help="Defocus step in nm.")
@click.option("--amplitude-contrast", type=float, default=0.07, show_default=True)
@click.option("--cs", type=float, default=2.7, show_default=True, help="Spherical aberration, mm.")
@click.option("--volt", type=float, default=300, show_default=True, help="Accelerating voltage, kV.")
@click.option("--cores", type=int, default=4, show_default=True,
              help="Parallel workers for CTF correction / flip+filter steps.")
@click.option("--final-bin", type=int, required=True,
              help="Final binning factor of the output tomogram (must be a multiple of --bin-factor).")
@click.option("--cleanup/--keep-intermediates", default=False, show_default=True,
              help="Delete corrected/flipped/filtered stacks and the untrimmed .rec once done.")
@click.option("--force", is_flag=True, default=False,
              help="Redo every step even if its output already exists (disables resuming).")
@click.option("--dry-run", is_flag=True, default=False, help="Print commands without running them.")
@click.option("--log-file", default="novaCTF_process.log", show_default=True,
              help="Log file name, created in --workdir.")
def main(**kwargs):
    """Run the full novaCTF correction pipeline end-to-end (resumable)."""
    os.chdir(kwargs.pop("workdir"))
    log_path = Path(kwargs.pop("log_file")).resolve()

    if kwargs["secs"] is None and not kwargs["ta_solution_log"].exists() and not kwargs["dry_run"]:
        raise click.UsageError(
            "Neither --secs nor a readable --ta-solution-log was provided; "
            "pass one of them so the newstack -secs range can be determined."
        )

    if kwargs.get("xml_dir") is None:
        kwargs["xml_dir"] = Path(".")

    bin_factor = kwargs["bin_factor"]

    tomo_size = kwargs.pop("tomo_size")
    parts = tomo_size.split(",")
    if len(parts) != 3:
        raise click.UsageError(
            f"--tomo-size must be 'X,Y,thickness' (three comma-separated integers), got {tomo_size!r}"
        )
    try:
        unbinned_x, unbinned_y, unbinned_thickness = (int(p) for p in parts)
    except ValueError:
        raise click.UsageError(
            f"--tomo-size must be three integers separated by commas, got {tomo_size!r}"
        )

    def to_binned(value, label):
        if value % bin_factor != 0:
            raise click.UsageError(
                f"--tomo-size {label} ({value}) is not evenly divisible by --bin-factor ({bin_factor})"
            )
        return value // bin_factor

    x = to_binned(unbinned_x, "X")
    y = to_binned(unbinned_y, "Y")
    thickness = to_binned(unbinned_thickness, "thickness")
    kwargs["size_xy"] = f"{x},{y}"
    kwargs["thickness"] = thickness
    click.echo(
        f"-> Binned tomogram size: {x},{y},{thickness} px "
        f"(unbinned {unbinned_x},{unbinned_y},{unbinned_thickness} / bin {bin_factor})"
    )

    pixel_size_angstrom = kwargs.pop("pixel_size")
    binned_pixel_size_nm = round((pixel_size_angstrom * bin_factor) / 10.0, 6)
    kwargs["pixel_size"] = binned_pixel_size_nm
    click.echo(
        f"-> Binned pixel size: {binned_pixel_size_nm} nm "
        f"(unbinned {pixel_size_angstrom} Å * bin {bin_factor} / 10)"
    )

    cfg = Config(**kwargs)

    log_line(log_path, f"=== novaCTF pipeline started for {cfg.name} ===")

    try:
        if cfg.secs is None:
            if cfg.ta_solution_log.exists():
                cfg.secs = parse_ta_solution_sections(cfg.ta_solution_log)
            elif cfg.dry_run:
                cfg.secs = "1-1"
                click.echo("(dry-run) no taSolution.log -- using placeholder --secs 1-1")
            else:
                raise NovaCTFError(f"{cfg.ta_solution_log} does not exist and no --secs was given")
        else:
            click.echo(f"-> Using manual --secs {cfg.secs}")

        step_newstack(cfg, log_path)
        step_xml_and_tlt(cfg, log_path)
        step_validate(cfg, log_path)
        step_defocus(cfg, log_path)

        n = count_defocus_stacks(cfg)
        click.echo(f"-> {n} defocus-shifted stacks to process")

        step_ctf_correction(cfg, log_path, n)
        step_flip_filter(cfg, log_path, n)
        step_reconstruction(cfg, log_path)
        step_trim(cfg, log_path)
        step_bin(cfg, log_path)
        step_cleanup(cfg, log_path, n)
    except NovaCTFError as exc:
        log_line(log_path, f"PIPELINE FAILED: {exc}")
        click.echo(f"\nERROR: {exc}\n", err=True)
        click.echo(
            "The pipeline stopped here. Fix the issue above and rerun the same "
            "command -- completed steps will be skipped automatically.",
            err=True,
        )
        sys.exit(1)

    log_line(log_path, f"=== novaCTF pipeline finished successfully for {cfg.name} ===")
    click.echo(f"\nDone. Final tomogram: {cfg.nova_bin}")


if __name__ == "__main__":
    main()
