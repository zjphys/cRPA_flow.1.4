# Workflow Quick Operation Guide

This guide covers every public command in `workflow.sh` and all calculation
stages. Run commands from the directory containing `workflow.sh`.

## 1. Setup

Required user input:

```text
POSCAR
```

Use a VASP 5-style POSCAR with element symbols on line 6. Provide an existing
`POTCAR`, or configure VASPKIT so the workflow can generate one. Review
`workflow.conf` before running, especially the VASP command, modules, Slurm
resources, pseudopotential setup, and INCAR templates.

Show the built-in help:

```bash
./workflow.sh --help
```

## 2. Stages

| Stage | Purpose | Main prerequisite |
|---|---|---|
| `00_relax` | Structural relaxation | `POSCAR`, `POTCAR` |
| `01_scf` | Self-consistent calculation | Relaxed `CONTCAR`, or input `POSCAR` |
| `02_dos` | DOS calculation | Completed `01_scf/CHGCAR` |
| `03_band` | Band structure | Completed `01_scf/CHGCAR`, generated `KPATH.in` |
| `04_wann` | Wannier construction | Completed SCF, DOS, and band data |
| `05_crpa` | cRPA calculation | Completed `04_wann` and selected Wannier states |

`02_dos` and `03_band` are independent children of `01_scf`. Wannier and cRPA
are deliberately prepared and launched separately.

## 3. Prepare and Run Stages 00–03

Generate the relaxation, SCF, DOS, and band inputs:

```bash
./workflow.sh prepare
```

Skip relaxation when the supplied POSCAR is already the desired structure:

```bash
./workflow.sh prepare --no-relax
```

Refresh workflow-owned stage directories after changing configuration:

```bash
./workflow.sh prepare --force
```

Run all prepared base stages sequentially in the current shell:

```bash
./workflow.sh run
```

Submit the base Slurm pipeline. SCF depends on relaxation when enabled; DOS and
bands both depend on SCF:

```bash
./workflow.sh submit
# Optional shared prefix: Pu-relax, Pu-scf, Pu-dos, Pu-band
./workflow.sh submit --job-name Pu
```

Both `--job-name Pu` and `--job-name=Pu` are accepted. Without this option,
the job names embedded in the generated `job.sh` files remain in effect.

Every generated `job.sh` is self-contained and may instead be submitted from
its own stage directory, for example:

```bash
cd 01_scf
sbatch job.sh
```

Submit downstream stages separately only after their prerequisite outputs
exist, or add the appropriate Slurm dependency yourself. The job synchronizes
required sibling-stage inputs, changes into its own directory, and runs the
command and setup captured when it was generated. Regenerate stage inputs and
jobs after changing `workflow.conf`.

Run exactly one prepared stage directly, mainly for debugging:

```bash
./workflow.sh execute 00_relax
./workflow.sh execute 01_scf
./workflow.sh execute 02_dos
./workflow.sh execute 03_band
./workflow.sh execute 04_wann
./workflow.sh execute 05_crpa
```

Check whether every registered stage is prepared, started, or finished:

```bash
./workflow.sh status
```

## 4. Plot Bands and DOS

After `02_dos` and `03_band` finish, plot all available element projections:

```bash
./workflow.sh postprocess --emin -3 --emax 4 --title "My material"
```

Plot selected elements and write PNG, PDF, and SVG files:

```bash
./workflow.sh postprocess \
  --elements Mn Sb \
  --format png --format pdf --format svg \
  --output mn_sb_band_dos
```

Plot orbital components for one element:

```bash
./workflow.sh postprocess \
  --orbital-element Mn \
  --orbitals dxy dyz dz2 dxz dx2-y2 \
  --emin -5 --emax 5
```

Reuse existing `PBAND_*.dat` and `PDOS_*.dat` without running VASPKIT again.
Spin-polarized `_UP.dat`/`_DW.dat` pairs are detected automatically:

```bash
./workflow.sh postprocess --reuse-data --dpi 300 --dos-max 20
```

After `04_wann` finishes, optionally overlay its interpolated bands on the
projected DFT band panel while keeping the same DOS panel:

```bash
./workflow.sh postprocess --wannier-bands --emin -3 --emax 4
```

The overlay reads `04_wann/INCAR` to determine `ISPIN`. It uses
`wannier90_band.dat` for `ISPIN=1`, or both `wannier90.1_band.dat` and
`wannier90.2_band.dat` for `ISPIN=2`. Raw Wannier energies are shifted by the
final `E-fermi` value in `04_wann/OUTCAR`; spin up is solid and spin down is
dashed. Use `--reuse-data` as well when the existing PBAND/PDOS files should
not be regenerated.

Other useful plotting controls are `--marker-scale`, `--vaspkit`, and repeated
`--format` options. The defaults are PNG and PDF output.

## 5. Prepare and Run Stage 04: Wannier

Choose positionally paired element/orbital projections and the number of
Wannier functions. UP and DW target bands are ranked independently. When a
channel's selected indices are noncontiguous, its largest contiguous run
drives the energy windows; tied runs prefer the one containing the
highest-ranked band. All selected bands still count toward `NUM_WANN` and
remain in `wannier_band_ranking.csv`.

```bash
./workflow.sh prepare-wannier \
  --elements Mn Sb \
  --orbitals d p
```

Without `--num-bands`, `NUM_WANN` is inferred from `01_scf/POSCAR`: each
paired element count is multiplied by `s=1`, `p=3`, `d=5`, or `f=7` and the
results are summed. Wannier ranking accepts only these aggregate shells because
the SCF stage uses `LORBIT=10`; component names such as `dxy` are rejected. Use
`--num-bands N` to override the inferred count.

Ranking is read directly from `01_scf/PROCAR`, not from the symmetry-line band
path. For every band and spin channel, the requested ion/shell projections are
summed over the complete irreducible SCF mesh using the `weight =` value written
by VASP, then normalized by the total k-point weight. `02_dos/EIGENVAL` still
provides the energy extrema and Wannier windows. A missing or empty SCF PROCAR
is a hard error; rerun or restart the SCF projection output before preparation.

The Gamma-centered Wannier mesh uses `KPR_WANN` from `workflow.conf`
(default `0.04`). Pass `--kpr VALUE` to override it for one preparation.

Example with every optional preparation control:

```bash
./workflow.sh prepare-wannier \
  --elements Mn Sb \
  --orbitals d p \
  --num-bands 22 \
  --kpr 0.04 \
  --frozen-margin 0.1 \
  --vaspkit vaspkit \
  --force
```

Run or submit only the prepared Wannier stage:

```bash
./workflow.sh run-wannier
./workflow.sh submit-wannier --job-name Pu  # Pu-wann
```

Inspect `04_wann/OUTCAR`, the Wannier90 output, interpolated bands, and
`wannier_band_ranking.csv` before choosing cRPA target states.

For collinear spin-polarized calculations, `prepare-wannier` reads the UP and
DW blocks of SCF `PROCAR` and ranks them independently. It selects the largest
contiguous run in each channel, then builds one common frozen window from the
intersection of the chosen-run and neighboring-band guard ranges.

## 6. Prepare and Run Stage 05: cRPA

After Wannier finishes, omit `--target-states` to select every state from `1`
through `NUM_WANN`:

```bash
./workflow.sh prepare-crpa
```

To select a subset, one-based inclusive ranges and individual indices may be
mixed:

```bash
./workflow.sh prepare-crpa --target-states 1-5 8 10-12
```

Override all optional cRPA preparation values when required:

```bash
./workflow.sh prepare-crpa \
  --target-states 1-10 \
  --nbandsgw 160 \
  --encutgw 350 \
  --kpar 4 \
  --force
```

`NBANDSGW` defaults to the effective completed-Wannier `NBANDS`, `ENCUTGW`
defaults to two-thirds of `ENCUT`, and `KPAR` comes from `CRPA_KPAR` in
`workflow.conf`. For spin-polarized calculations, add `ISPIN` and `MAGMOM`
manually to `INCAR_CRPA_TEMPLATE`; preparation checks the effective `ISPIN`
against the spin-channel count stored in `WANPROJ`.

Run or submit only the prepared cRPA stage:

```bash
./workflow.sh run-crpa
./workflow.sh submit-crpa --job-name Pu  # Pu-crpa
```

## 7. Complete Example

Local sequential workflow:

```bash
./workflow.sh prepare --no-relax
./workflow.sh run
./workflow.sh postprocess --elements Mn Sb --emin -4 --emax 4
./workflow.sh prepare-wannier \
  --elements Mn Sb --orbitals d p
./workflow.sh run-wannier
./workflow.sh prepare-crpa
./workflow.sh run-crpa
./workflow.sh status
```

Slurm workflow:

```bash
./workflow.sh prepare --no-relax
./workflow.sh submit --job-name Pu
# Wait for SCF, DOS, and band jobs to finish.
./workflow.sh prepare-wannier \
  --elements Mn Sb --orbitals d p
./workflow.sh submit-wannier --job-name Pu
# Wait for the Wannier job and inspect its results.
./workflow.sh prepare-crpa
./workflow.sh submit-crpa --job-name Pu
./workflow.sh status
```

`submit-wannier` and `submit-crpa` do not automatically add dependencies to
earlier jobs; invoke them only after their prerequisite calculations finish.
Their generated jobs may also be submitted directly with `cd 04_wann && sbatch
job.sh` or `cd 05_crpa && sbatch job.sh`.
