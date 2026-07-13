# novaCTF_wrapper

A Python/Click wrapper (and tutorial) for producing novaCTF 3D-CTF-corrected
tomograms that match Warp-generated tomograms exactly, so they can be used
directly for template matching.

novaCTF ([Turoňová et al. 2017](https://doi.org/10.1016/j.jsb.2017.07.007))
does depth-dependent CTF correction during weighted back-projection (WBP)
reconstruction, but it's a multi-step, multi-command pipeline with a lot of
places to get a parameter wrong. `novactf_pipeline.py` runs the whole thing
end-to-end from one command, is resumable if it crashes partway through, and
logs every command it runs with a timestamp.

## Requirements

- **WSL or native Linux.** novaCTF, IMOD, and `xml2pytom.py` all need to run
  under Linux -- running `xml2pytom.py` on Windows produces `.tlt`/`.defocus`
  files with the wrong line endings/formatting.
- [IMOD](https://bio3d.colorado.edu/imod/) (`newstack`, `clip`, `header`,
  `trimvol`, `binvol`) on your `PATH`.
- [novaCTF](https://github.com/turonova/novaCTF) compiled and on your `PATH`.
- Python 3 with `click` installed (`pip install click`).
- A working `xml2pytom.py` (included) that turns a Warp `.xml` per-tilt-series
  file into the `.tlt`/`.defocus`/dose `.txt` files novaCTF needs. **Edit the
  constants at the top of `xml2pytom.py`** (`tomo_dir`, `ref_name`,
  `mask_name`, box/search sizes, GPU list, pixel size, etc.) for your dataset
  before running it -- it's a per-project script, not a CLI tool.

## Installation

```bash
git clone <this-repo-url>
cd novaCTF_wrapper
pip install -r requirements.txt
```

## Pipeline overview

`novactf_pipeline.py` runs these stages in order. Every stage checks whether
its expected output already exists and skips it if so, so if the pipeline
dies partway through (crash, disk full, wrong parameter caught mid-run), just
fix the problem and rerun the same command -- completed stages are skipped
automatically (`--force` to redo everything anyway).

1. **`taSolution.log` parsing** -- finds which tilt views survived alignment
   (the "view / rotation / tilt / deltilt / mag / ..." solution table) and
   turns them into a `newstack -secs` range, e.g. `4-36`. You can bypass this
   with `--secs` if you already know the range.
2. **`newstack`** -- builds the aligned, binned 2D stack from the raw `.st` +
   `.xf` files, using the sections from step 1.
3. **xml copy + `xml2pytom.py`** -- copies the Warp `.xml` for this
   tilt-series in and runs `xml2pytom.py` to produce the `.tlt` and
   `.defocus` files.
4. **Validation** -- compares the number of sections in the aligned stack
   against the number of lines in the `.tlt` file and fails immediately (with
   a clear message) if they don't match, instead of failing confusingly deep
   into the pipeline.
5. **`novaCTF -Algorithm defocus`** -- generates N defocus-shifted files,
   where N depends on `THICKNESS`, `PixelSize`, and `DefocusStep`.
6. **`novaCTF -Algorithm ctfCorrection`** -- CTF-corrects the stack once per
   defocus file, in parallel across `--cores`.
7. **`clip flipyz` + `novaCTF -Algorithm filterProjections`** -- reorients
   each corrected stack into the XZY layout novaCTF's reconstruction expects,
   then filters it. Also parallelized across `--cores`.
8. **`novaCTF -Algorithm 3dctf`** -- WBP reconstruction with 3D-CTF
   correction, combining all N filtered stacks into one tomogram.
9. **`trimvol -yz`** -- rotates the reconstruction back to standard XYZ
   orientation (novaCTF's output is left in the flipped XZY layout from
   step 7).
10. **`binvol`** -- produces the final tomogram at `--final-bin`, which must
    be an exact multiple of `--bin-factor`.
11. **Optional cleanup** -- deletes the corrected/flipped/filtered
    intermediate stacks and the untrimmed `.rec` once the final tomogram
    exists (`--cleanup`; off by default).

Reconstruction itself (step 8) is a single serial process -- novaCTF has no
internal multithreading/GPU support, so only steps 6 and 7 benefit from
`--cores`.

## Usage

```bash
python3 novactf_pipeline.py \
    --name MIM019_2_lam1_ts_002 \
    --ta-solution-log taSolution.log \
    --xml-dir /path/to/warp/xml/files \
    --tomo-size 4096,5760,2400 \
    --bin-factor 2 \
    --pixel-size 1.66 \
    --defocus-step 20 \
    --final-bin 6 \
    --cores 4
```

Run `python3 novactf_pipeline.py --help` for the full option list. The two
parameters worth calling out:

- **`--tomo-size X,Y,thickness`** is given in **unbinned** pixels (e.g. the
  full-frame camera size and your desired sample thickness before binning).
  It's divided by `--bin-factor` to get novaCTF's `FULLIMAGE`/`THICKNESS`.
- **`--pixel-size`** is the **unbinned** pixel size in **Angstrom**. It's
  converted to binned nanometers internally
  (`pixel_size_Å * bin_factor / 10`) for novaCTF's `-PixelSize` flag, since
  novaCTF wants nm at the binned sampling, not Å at full resolution.

A log of every command run (with timestamps) is written to
`novaCTF_process.log` in the working directory.

## Parameter reference

| CLI option | novaCTF/IMOD flag(s) | Meaning |
|---|---|---|
| `--tomo-size X,Y,thickness` (unbinned) | `-SizeToOutputInXandY`, `-FULLIMAGE`, `-THICKNESS` | Divided by `--bin-factor` first. `THICKNESS` is the full Z box height of the reconstruction, not a half-thickness. |
| `--bin-factor` | `-BinByFactor` | Binning applied in `newstack`; also the divisor for `--tomo-size` and multiplier for `--pixel-size`. |
| `--pixel-size` (unbinned, Å) | `-PixelSize` (binned, nm) | Converted as `pixel_size * bin_factor / 10`. |
| `--shift` | `-SHIFT` | Must be `0.0,0.0` for the defocus-generation step; only meaningful as a real Z-offset during reconstruction. |
| `--correction-type` | `-CorrectionType` | `phaseflip` (recommended) or `multiplication`. |
| `--defocus-format` | `-DefocusFileFormat` | `imod`, `ctffind4`, or `gctf`. |
| `--correct-astigmatism` | `-CorrectAstigmatism` | `1` or `0`. |
| `--defocus-step` | `-DefocusStep` | In nm; smaller steps are more accurate at the cost of more sub-stacks/compute. |
| `--amplitude-contrast`, `--cs`, `--volt` | `-AmplitudeContrast`, `-Cs`, `-Volt` | Microscope parameters, used only in the `ctfCorrection` step. |
| `--final-bin` | -- | Must be an exact multiple of `--bin-factor`; `binvol -binning` gets `final_bin / bin_factor`. |

## Notes from getting this working

- **`THICKNESS` must be identical everywhere it's used** (defocus generation
  *and* reconstruction). It determines how many defocus-shifted stacks get
  generated; changing it only for the reconstruction step will make novaCTF
  expect a different number of stacks than actually exist.
- **`PixelSize` must be the same in every novaCTF call**, and must be the
  *binned* value in nm, not the unbinned value or the value in Å. Getting
  this wrong silently changes how many defocus stacks get generated, which
  shows up as `ExceptionFileOpen` errors deep into a later step.
- A **low-contrast 3D-CTF-corrected tomogram is expected**, not a bug --
  that's the tradeoff of proper CTF correction with plain WBP. Don't add
  `RADIAL`/`FakeSIRTiterations` into this pipeline to "fix" it, since that
  would make the tomogram diverge from Warp's rather than match it. If you
  want something nicer to look at for picking, reconstruct a **separate**
  tomogram with plain IMOD `tilt` (no CTF correction) using
  `-FakeSIRTiterations 15`, from the plain aligned stack -- not from any of
  the CTF-corrected/flipped stacks in this pipeline.
- **Windows line endings (`\r\n`) in a `.tlt` file break `tilt`/novaCTF**
  with confusing "End of file before all values gotten" errors even when the
  line count is correct. Run `dos2unix` on any `.tlt` file that started life
  on Windows.
- **novaCTF has no built-in parallelization** (checked the source -- no
  OpenMP/pthreads/MPI/GPU anywhere). The only real speedup available is
  running the independent per-stack `ctfCorrection`/`filterProjections`
  calls concurrently, which is what `--cores` does here.
- **`novaCTF`'s reconstruction output is in XZY orientation**, not standard
  XYZ, because of the `clip flipyz` step before filtering. `trimvol -yz`
  undoes this; if `trimvol` errors reading the file, the fallback is
  `clip flipyz` followed by `clip flipz`.

## Files

- `novactf_pipeline.py` -- the full pipeline, run this.
- `xml2pytom.py` -- converts a Warp `.xml` into `.tlt`/`.defocus`/dose files;
  also emits a `pytom_match_template.py` submission script. Edit the
  constants at the top per dataset.
- `ctf_correction.sh`, `flip_filter.sh` -- earlier, single-purpose shell
  scripts covering steps 6 and 7 individually. Kept for reference/debugging;
  `novactf_pipeline.py` supersedes them.

## Citing novaCTF

Turoňová, B., Schur, F.K.M., Wan, W. and Briggs, J.A.G. *Efficient 3D-CTF
correction for cryo-electron tomography using NovaCTF improves subtomogram
averaging resolution to 3.4 Å.* J Struct Biol. 2017.
[doi:10.1016/j.jsb.2017.07.007](https://doi.org/10.1016/j.jsb.2017.07.007)
