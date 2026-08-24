#!/usr/bin/env bash
set -euo pipefail

# POSCAR-first VASP workflow: relax -> SCF -> DOS/bands, then optional Wannier/cRPA.
# Optional site-specific values can be placed in workflow.conf.

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${WORKFLOW_CONFIG:-$ROOT_DIR/workflow.conf}"
[[ -f "$CONFIG_FILE" ]] && source "$CONFIG_FILE"

VASPKIT_BIN="${VASPKIT_BIN:-vaspkit}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
VASP_COMMAND="${VASP_COMMAND:-vasp_std}"
SUBMIT_COMMAND="${SUBMIT_COMMAND:-sbatch}"
KPR_RELAX="${KPR_RELAX:-0.05}"
KPR_SCF="${KPR_SCF:-0.04}"
KPR_DOS="${KPR_DOS:-0.03}"
KPR_WANN="${KPR_WANN:-0.04}"
EDIFF="${EDIFF:-1e-6}"
EDIFFG="${EDIFFG:--0.02}"
ISIF="${ISIF:-3}"
NSW="${NSW:-200}"
NEDOS="${NEDOS:-2000}"
ENCUT_FACTOR="${ENCUT_FACTOR:-1.30}"
ENCUT="${ENCUT:-auto}"
GENERATE_BANDS="${GENERATE_BANDS:-yes}"
RUN_RELAX="${RUN_RELAX:-yes}"
SBATCH_PARTITION="${SBATCH_PARTITION:-}"
SBATCH_NODES="${SBATCH_NODES:-1}"
SBATCH_NTASKS_PER_NODE="${SBATCH_NTASKS_PER_NODE:-1}"
SBATCH_CPUS_PER_TASK="${SBATCH_CPUS_PER_TASK:-1}"
SBATCH_TIME="${SBATCH_TIME:-24:00:00}"
SBATCH_EXTRA="${SBATCH_EXTRA:-}"
CRPA_SBATCH_PARTITION="${CRPA_SBATCH_PARTITION-$SBATCH_PARTITION}"
CRPA_SBATCH_NODES="${CRPA_SBATCH_NODES-$SBATCH_NODES}"
CRPA_SBATCH_NTASKS_PER_NODE="${CRPA_SBATCH_NTASKS_PER_NODE-$SBATCH_NTASKS_PER_NODE}"
CRPA_SBATCH_CPUS_PER_TASK="${CRPA_SBATCH_CPUS_PER_TASK-$SBATCH_CPUS_PER_TASK}"
CRPA_SBATCH_TIME="${CRPA_SBATCH_TIME-$SBATCH_TIME}"
CRPA_SBATCH_EXTRA="${CRPA_SBATCH_EXTRA-$SBATCH_EXTRA}"
CRPA_KPAR="${CRPA_KPAR:-$CRPA_SBATCH_NODES}"
EXECUTION_SETUP="${EXECUTION_SETUP:-}"

# Complete templates and stage commands can be replaced in workflow.conf.
if [[ -z "${INCAR_RELAX_TEMPLATE+x}" ]]; then
  INCAR_RELAX_TEMPLATE=$(cat <<'EOF'
SYSTEM = {{SYSTEM}}
ENCUT = {{ENCUT}}
PREC = Accurate
EDIFF = {{EDIFF}}
EDIFFG = {{EDIFFG}}
IBRION = 2
ISIF = {{ISIF}}
NSW = {{NSW}}
ISMEAR = 0
SIGMA = 0.05
LREAL = Auto
LASPH = .TRUE.
LWAVE = .FALSE.
LCHARG = .FALSE.
EOF
)
fi

if [[ -z "${INCAR_SCF_TEMPLATE+x}" ]]; then
  INCAR_SCF_TEMPLATE=$(cat <<'EOF'
SYSTEM = {{SYSTEM}} SCF
ENCUT = {{ENCUT}}
PREC = Accurate
EDIFF = {{EDIFF}}
ISMEAR = 0
SIGMA = 0.05
LREAL = Auto
LASPH = .TRUE.
LORBIT = 10
LWAVE = .TRUE.
LCHARG = .TRUE.
EOF
)
fi

if [[ -z "${INCAR_DOS_TEMPLATE+x}" ]]; then
  INCAR_DOS_TEMPLATE=$(cat <<'EOF'
SYSTEM = {{SYSTEM}} DOS
ENCUT = {{ENCUT}}
PREC = Accurate
EDIFF = {{EDIFF}}
ICHARG = 11
ISMEAR = -5
LREAL = Auto
LASPH = .TRUE.
LORBIT = 10
NEDOS = {{NEDOS}}
LWAVE = .FALSE.
LCHARG = .FALSE.
EOF
)
fi

if [[ -z "${INCAR_BAND_TEMPLATE+x}" ]]; then
  INCAR_BAND_TEMPLATE=$(cat <<'EOF'
SYSTEM = {{SYSTEM}} bands
ENCUT = {{ENCUT}}
PREC = Accurate
EDIFF = {{EDIFF}}
ICHARG = 11
ISMEAR = 0
SIGMA = 0.05
LREAL = Auto
LASPH = .TRUE.
LORBIT = 10
LWAVE = .FALSE.
LCHARG = .FALSE.
EOF
)
fi

if [[ -z "${INCAR_CRPA_TEMPLATE+x}" ]]; then
  INCAR_CRPA_TEMPLATE=$(cat <<'EOF'
SYSTEM = {{SYSTEM}} cRPA
ENCUT = {{ENCUT}}
ISMEAR = 0; SIGMA = 0.1
EDIFF = {{EDIFF}}
NBANDS = {{NBANDS}}
NBANDSGW = {{NBANDSGW}}

# SOC can be enabled in an overridden INCAR_CRPA_TEMPLATE when required.
# LSORBIT = .TRUE.
ISYM = -1

LMAXMIX = 6

ALGO = CRPA
NTARGET_STATES = {{TARGET_STATES}}
ENCUTGW = {{ENCUTGW}}

KPAR = {{CRPA_KPAR}}
EOF
)
fi

if [[ -z "${SBATCH_TEMPLATE+x}" ]]; then
  SBATCH_TEMPLATE=$(cat <<'EOF'
#SBATCH --job-name={{JOB_NAME}}
{{SBATCH_PARTITION_LINE}}
#SBATCH --nodes={{SBATCH_NODES}}
#SBATCH --ntasks-per-node={{SBATCH_NTASKS_PER_NODE}}
#SBATCH --cpus-per-task={{SBATCH_CPUS_PER_TASK}}
#SBATCH --time={{SBATCH_TIME}}
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
{{SBATCH_EXTRA}}
EOF
)
fi

if ! declare -p STAGE_COMMANDS >/dev/null 2>&1; then
  declare -A STAGE_COMMANDS=([default]="$VASP_COMMAND")
elif [[ -z "${STAGE_COMMANDS[default]+x}" ]]; then
  STAGE_COMMANDS[default]="$VASP_COMMAND"
fi

STAGES=(00_relax 01_scf 02_dos)
MARKER=".generated-by-poscar-workflow"
STAGE_FILE="$ROOT_DIR/.workflow-stages"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '==> %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }

usage() {
  cat <<'EOF'
Usage:
  ./workflow.sh prepare [--force] [--no-relax]
                                      Generate all inputs from POSCAR
  ./workflow.sh run                Run stages sequentially in this shell
  ./workflow.sh submit [--job-name PREFIX]
                                    Submit an afterok-linked Slurm pipeline
  ./workflow.sh execute STAGE      Run one prepared stage directly
  ./workflow.sh status             Show prepared/output state
  ./workflow.sh postprocess [ARGS] Plot element- or orbital-projected bands and DOS
  ./workflow.sh prepare-wannier --elements E... --orbitals O... [--num-bands N]
                                    Generate 04_wann after DOS/bands finish
  ./workflow.sh run-wannier        Run only the prepared 04_wann stage
  ./workflow.sh submit-wannier [--job-name PREFIX]
                                    Submit only the prepared 04_wann stage
  ./workflow.sh prepare-crpa [--target-states I...]
                                    Generate 05_crpa after 04_wann finishes
  ./workflow.sh run-crpa           Run only the prepared 05_crpa stage
  ./workflow.sh submit-crpa [--job-name PREFIX]
                                    Submit only the prepared 05_crpa stage

Only POSCAR is required as user-supplied scientific input. VASPKIT must be
configured with a licensed pseudopotential library to generate POTCAR. It also
generates Gamma-centered relaxation/SCF/DOS meshes and the symmetry-aware band
path. Copy workflow.conf.example to
workflow.conf only when site or calculation defaults need changing.
That file can also replace each complete INCAR template, the Slurm header,
shared runtime setup, and individual stage commands.

Use --no-relax when POSCAR is already the structure that should be used for
SCF, DOS, and band calculations.

Each generated stage/job.sh is self-contained. It may be submitted separately
with "cd STAGE && sbatch job.sh" after any required predecessor outputs exist.
Regenerate jobs after changing runtime commands, setup, resources, or headers.

Submission commands accept "--job-name PREFIX" or "--job-name=PREFIX". The
stage name is appended automatically, for example "--job-name Pu" submits
Pu-relax, Pu-scf, Pu-dos, and Pu-band.

After the DOS and band stages finish, run "workflow.sh postprocess". VASPKIT
tasks 213 and 113 generate element-projected band and DOS data. Plot
options such as "--emin -3 --emax 4 --title Material" are forwarded to
postprocess.py. Python, NumPy, Matplotlib, and VASPKIT are required.
Use "--orbital-element Am --orbitals s p d f" for aggregate shell projections,
or select components such as "--orbital-element Ni --orbitals dz2 dx2-y2".
The DOS panel always includes the selected element total(s).

To prepare a Wannier calculation after SCF and DOS finish (after inspecting the
optional band stage when available), use for example:
  workflow.sh prepare-wannier --elements Mn Sb --orbitals d p
Elements and orbitals are paired by position. NUM_WANN defaults to the sum of
the matching POSCAR atom counts times shell multiplicities (s=1, p=3, d=5,
f=7); --num-bands N overrides it. The command accepts --kpr VALUE
(default KPR_WANN, 0.04), --frozen-margin EV (default 0.1), and --force. It
ranks aggregate s/p/d/f weights directly from the full-zone 01_scf/PROCAR,
using VASP k-point integration weights, then generates 04_wann and records the
selected-band energy frontier. A missing SCF PROCAR is an error.

After 04_wann finishes successfully, prepare cRPA. By default all Wannier
states from 1 through NUM_WANN are excluded from screening:
  workflow.sh prepare-crpa
Use --target-states to select a subset, for example:
  workflow.sh prepare-crpa --target-states 1-5 8 10-12
Inclusive ranges and individual indices can be mixed, for example 1-5 8 10-12.
The command accepts --nbandsgw N, --encutgw EV, and --force. By default it
uses the completed Wannier NBANDS for NBANDSGW and two-thirds of ENCUT for
ENCUTGW. It copies the required Wannier restart inputs into 05_crpa.
EOF
}

is_enabled() {
  local wanted="$1"
  [[ -s "$STAGE_FILE" ]] && grep -Fxq "$wanted" "$STAGE_FILE"
}

validate_poscar() {
  [[ -s "$ROOT_DIR/POSCAR" ]] || die "POSCAR is missing or empty in $ROOT_DIR"
  local symbols counts
  symbols="$(awk 'NR==6 {print}' "$ROOT_DIR/POSCAR")"
  counts="$(awk 'NR==7 {print}' "$ROOT_DIR/POSCAR")"
  [[ "$symbols" =~ [A-Za-z] ]] ||
    die "VASP 4-style POSCAR detected. Add element symbols on line 6."
  [[ "$counts" =~ [0-9] ]] || die "Could not read atom counts on POSCAR line 7."
}

system_name() {
  local title
  title="$(head -n 1 "$ROOT_DIR/POSCAR" | tr -d '\r')"
  [[ -n "$title" ]] && printf '%s' "$title" || printf 'VASP calculation'
}

generate_potcar() {
  if [[ -s "$ROOT_DIR/POTCAR" ]]; then
    info "Reusing existing POTCAR"
    return
  fi
  command -v "$VASPKIT_BIN" >/dev/null 2>&1 ||
    die "POTCAR is absent and '$VASPKIT_BIN' is not available. Install/configure VASPKIT or provide POTCAR."

  local tmp
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/poscar-workflow.XXXXXX")"
  cp "$ROOT_DIR/POSCAR" "$tmp/POSCAR"
  info "Generating POTCAR with VASPKIT task 103"
  if ! (cd "$tmp" && printf '103\n' | "$VASPKIT_BIN"); then
    rm -rf "$tmp"
    die "VASPKIT failed to generate POTCAR. Check its pseudopotential path."
  fi
  [[ -s "$tmp/POTCAR" ]] || { rm -rf "$tmp"; die "VASPKIT completed but did not create POTCAR."; }
  cp "$tmp/POTCAR" "$ROOT_DIR/POTCAR"
  rm -rf "$tmp"
}

resolve_encut() {
  if [[ "$ENCUT" != "auto" ]]; then
    printf '%s' "$ENCUT"
    return
  fi
  awk -v factor="$ENCUT_FACTOR" '
    /ENMAX/ {
      line=$0
      sub(/^.*ENMAX[[:space:]]*=[[:space:]]*/, "", line)
      sub(/[;[:space:]].*$/, "", line)
      if ((line+0) > max) max=line+0
    }
    END {
      if (max <= 0) exit 1
      value=max*factor
      rounded=int((value+4.999999)/5)*5
      print rounded
    }' "$ROOT_DIR/POTCAR" ||
    die "Could not determine ENMAX from POTCAR; set ENCUT in workflow.conf."
}

generate_kpath() {
  [[ "$GENERATE_BANDS" == "yes" ]] || return 1
  if [[ -s "$ROOT_DIR/KPATH.in" ]]; then
    info "Reusing existing KPATH.in"
    return 0
  fi
  if ! command -v "$VASPKIT_BIN" >/dev/null 2>&1; then
    warn "VASPKIT unavailable; skipping the band stage."
    return 1
  fi

  local tmp
  tmp="$(mktemp -d "${TMPDIR:-/tmp}/poscar-kpath.XXXXXX")"
  cp "$ROOT_DIR/POSCAR" "$tmp/POSCAR"
  info "Generating symmetry-aware band path with VASPKIT task 303"
  if (cd "$tmp" && printf '303\n' | "$VASPKIT_BIN") && [[ -s "$tmp/KPATH.in" ]]; then
    cp "$tmp/KPATH.in" "$ROOT_DIR/KPATH.in"
    rm -rf "$tmp"
    return 0
  fi
  rm -rf "$tmp"
  warn "Band-path generation failed; DOS workflow remains available."
  return 1
}

generate_gamma_kpoints() {
  local stage="$1" kpr="$2" tmp
  command -v "$VASPKIT_BIN" >/dev/null 2>&1 ||
    die "'$VASPKIT_BIN' is required to generate $stage/KPOINTS."
  awk -v value="$kpr" 'BEGIN { exit !(value + 0 > 0) }' ||
    die "The KPR for $stage must be positive; received '$kpr'."

  tmp="$(mktemp -d "${TMPDIR:-/tmp}/poscar-kpoints.XXXXXX")"
  cp "$ROOT_DIR/POSCAR" "$tmp/POSCAR"
  info "Generating Gamma-centered $stage/KPOINTS with VASPKIT task 102 (KPR=$kpr)"
  if ! (cd "$tmp" && printf '102\n2\n%s\n' "$kpr" | "$VASPKIT_BIN"); then
    rm -rf "$tmp"
    die "VASPKIT failed to generate $stage/KPOINTS with KPR=$kpr."
  fi
  [[ -s "$tmp/KPOINTS" ]] || {
    rm -rf "$tmp"
    die "VASPKIT completed but did not create $stage/KPOINTS."
  }
  cp "$tmp/KPOINTS" "$ROOT_DIR/$stage/KPOINTS"
  rm -rf "$tmp"
}

render_template() {
  local output="$1" token placeholder
  shift
  while (( $# >= 2 )); do
    token="$1"
    placeholder="{{$token}}"
    output="${output//$placeholder/$2}"
    shift 2
  done
  if [[ "$output" =~ \{\{[A-Z0-9_]+\}\} ]]; then
    die "Unrecognized template token '${BASH_REMATCH[0]}'."
  fi
  printf '%s\n' "$output"
}

write_incar() {
  local stage="$1" name="$2" encut="$3" template
  case "$stage" in
    00_relax) template="$INCAR_RELAX_TEMPLATE" ;;
    01_scf) template="$INCAR_SCF_TEMPLATE" ;;
    02_dos) template="$INCAR_DOS_TEMPLATE" ;;
    03_band) template="$INCAR_BAND_TEMPLATE" ;;
    *) die "No INCAR template is defined for stage '$stage'." ;;
  esac

  render_template "$template" \
    SYSTEM "$name" \
    ENCUT "$encut" \
    EDIFF "$EDIFF" \
    EDIFFG "$EDIFFG" \
    ISIF "$ISIF" \
    NSW "$NSW" \
    KPR_RELAX "$KPR_RELAX" \
    KPR_SCF "$KPR_SCF" \
    KPR_DOS "$KPR_DOS" \
    NEDOS "$NEDOS" \
    SBATCH_NODES "$SBATCH_NODES" > "$ROOT_DIR/$stage/INCAR"
}

write_job() {
  local stage="$1" relax_enabled="${2:-no}" dir="$ROOT_DIR/$1" partition_line=""
  local job_partition="$SBATCH_PARTITION"
  local job_nodes="$SBATCH_NODES"
  local job_ntasks_per_node="$SBATCH_NTASKS_PER_NODE"
  local job_cpus_per_task="$SBATCH_CPUS_PER_TASK"
  local job_time="$SBATCH_TIME"
  local job_extra="$SBATCH_EXTRA"
  local job_command="${STAGE_COMMANDS[$stage]:-${STAGE_COMMANDS[default]}}"
  if [[ -n "$EXECUTION_SETUP" ]]; then
    job_command="$EXECUTION_SETUP"$'\n'"$job_command"
  fi
  if [[ "$stage" == 05_crpa ]]; then
    job_partition="$CRPA_SBATCH_PARTITION"
    job_nodes="$CRPA_SBATCH_NODES"
    job_ntasks_per_node="$CRPA_SBATCH_NTASKS_PER_NODE"
    job_cpus_per_task="$CRPA_SBATCH_CPUS_PER_TASK"
    job_time="$CRPA_SBATCH_TIME"
    job_extra="$CRPA_SBATCH_EXTRA"
  fi
  [[ -n "$job_partition" ]] &&
    partition_line="#SBATCH --partition=$job_partition"
  {
    printf '%s\n' '#!/usr/bin/env bash'
    render_template "$SBATCH_TEMPLATE" \
      JOB_NAME "${stage#*_}" \
      SBATCH_PARTITION_LINE "$partition_line" \
      SBATCH_NODES "$job_nodes" \
      SBATCH_NTASKS_PER_NODE "$job_ntasks_per_node" \
      SBATCH_CPUS_PER_TASK "$job_cpus_per_task" \
      SBATCH_TIME "$job_time" \
      SBATCH_EXTRA "$job_extra"
    printf '\nset -euo pipefail\n'
    printf 'export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"\n'
    printf '%s\n' \
      'job_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"' \
      'job_error() { printf "ERROR: %s\\n" "$*" >&2; exit 1; }' \
      'require_input() { [[ -s "$1" ]] || job_error "required predecessor input is missing or empty: $1"; }'
    printf '[[ -f "$job_dir/%s" ]] || job_error "stage marker is missing: $job_dir/%s"\n' \
      "$MARKER" "$MARKER"
    case "$stage" in
      01_scf)
        if [[ "$relax_enabled" == yes ]]; then
          printf '%s\n' \
            'require_input "$job_dir/../00_relax/CONTCAR"' \
            'cp -- "$job_dir/../00_relax/CONTCAR" "$job_dir/POSCAR"'
        fi
        ;;
      02_dos|03_band)
        printf '%s\n' \
          'require_input "$job_dir/../01_scf/POSCAR"' \
          'require_input "$job_dir/../01_scf/CHGCAR"' \
          'cp -- "$job_dir/../01_scf/POSCAR" "$job_dir/POSCAR"' \
          'cp -- "$job_dir/../01_scf/CHGCAR" "$job_dir/CHGCAR"'
        ;;
      04_wann)
        printf '%s\n' \
          'require_input "$job_dir/../01_scf/POSCAR"' \
          'require_input "$job_dir/../01_scf/CHGCAR"' \
          'require_input "$job_dir/../01_scf/POTCAR"' \
          'cp -- "$job_dir/../01_scf/POSCAR" "$job_dir/POSCAR"' \
          'cp -- "$job_dir/../01_scf/CHGCAR" "$job_dir/CHGCAR"' \
          'cp -- "$job_dir/../01_scf/POTCAR" "$job_dir/POTCAR"'
        ;;
    esac
    printf '%s\n' 'cd "$job_dir"'
    printf 'exec bash -lc %q\n' "$job_command"
  } > "$dir/job.sh"
  chmod +x "$dir/job.sh"
}

register_stage() {
  local stage="$1"
  if ! is_enabled "$stage"; then
    printf '%s\n' "$stage" >> "$STAGE_FILE"
  fi
}

prepare() {
  local force=no run_relax="$RUN_RELAX" option
  for option in "$@"; do
    case "$option" in
      --force) force=yes ;;
      --no-relax) run_relax=no ;;
      *) die "Unknown prepare option: $option" ;;
    esac
  done
  validate_poscar
  generate_potcar
  local name encut stage
  name="$(system_name)"
  encut="$(resolve_encut)"
  info "Using ENCUT=$encut eV (max POTCAR ENMAX x $ENCUT_FACTOR, rounded up)"

  local have_band=no
  generate_kpath && have_band=yes || true
  STAGES=(01_scf 02_dos)
  [[ "$run_relax" == yes ]] && STAGES=(00_relax "${STAGES[@]}")
  [[ "$have_band" == yes ]] && STAGES+=(03_band)

  for stage in "${STAGES[@]}"; do
    if [[ -e "$ROOT_DIR/$stage" && ! -f "$ROOT_DIR/$stage/$MARKER" ]]; then
      die "$stage exists and was not generated by this workflow; refusing to overwrite it."
    fi
    if [[ -e "$ROOT_DIR/$stage" && "$force" != yes ]]; then
      die "$stage already exists. Use prepare --force to refresh generated inputs."
    fi
    mkdir -p "$ROOT_DIR/$stage"
    : > "$ROOT_DIR/$stage/$MARKER"
    cp "$ROOT_DIR/POSCAR" "$ROOT_DIR/$stage/POSCAR"
    cp "$ROOT_DIR/POTCAR" "$ROOT_DIR/$stage/POTCAR"
  done

  if [[ "$run_relax" == yes ]]; then
    write_incar 00_relax "$name" "$encut"
    generate_gamma_kpoints 00_relax "$KPR_RELAX"
  fi
  write_incar 01_scf "$name" "$encut"
  generate_gamma_kpoints 01_scf "$KPR_SCF"
  write_incar 02_dos "$name" "$encut"
  generate_gamma_kpoints 02_dos "$KPR_DOS"
  if [[ "$have_band" == yes ]]; then
    write_incar 03_band "$name" "$encut"
    cp "$ROOT_DIR/KPATH.in" "$ROOT_DIR/03_band/KPOINTS"
  fi
  for stage in "${STAGES[@]}"; do write_job "$stage" "$run_relax"; done
  printf '%s\n' "${STAGES[@]}" > "$STAGE_FILE"
  info "Prepared stages: ${STAGES[*]}"
}

prepare_wannier_stage() {
  local argument
  local -a kpr_arguments=(--kpr "$KPR_WANN")
  [[ -f "$ROOT_DIR/prepare_wannier.py" ]] ||
    die "prepare_wannier.py is missing from $ROOT_DIR"
  command -v "$PYTHON_BIN" >/dev/null 2>&1 ||
    die "Python command '$PYTHON_BIN' is unavailable."
  for argument in "$@"; do
    if [[ "$argument" == --kpr || "$argument" == --kpr=* ]]; then
      kpr_arguments=()
      break
    fi
  done
  "$PYTHON_BIN" "$ROOT_DIR/prepare_wannier.py" \
    --root "$ROOT_DIR" --vaspkit "$VASPKIT_BIN" \
    "${kpr_arguments[@]}" "$@"
  write_job 04_wann
  register_stage 04_wann
  info "Prepared 04_wann. Use run-wannier or submit-wannier."
}

prepare_crpa_stage() {
  [[ -f "$ROOT_DIR/prepare_crpa.py" ]] ||
    die "prepare_crpa.py is missing from $ROOT_DIR"
  command -v "$PYTHON_BIN" >/dev/null 2>&1 ||
    die "Python command '$PYTHON_BIN' is unavailable."
  INCAR_CRPA_TEMPLATE="$INCAR_CRPA_TEMPLATE" \
    "$PYTHON_BIN" "$ROOT_DIR/prepare_crpa.py" \
      --root "$ROOT_DIR" --kpar "$CRPA_KPAR" "$@"
  write_job 05_crpa
  register_stage 05_crpa
  info "Prepared 05_crpa. Use run-crpa or submit-crpa."
}

sync_predecessor() {
  local stage="$1"
  case "$stage" in
    00_relax) ;;
    01_scf)
      if is_enabled 00_relax; then
        [[ -s "$ROOT_DIR/00_relax/CONTCAR" ]] || die "00_relax/CONTCAR is missing."
        cp "$ROOT_DIR/00_relax/CONTCAR" "$ROOT_DIR/01_scf/POSCAR"
      fi
      ;;
    02_dos|03_band|04_wann)
      [[ -s "$ROOT_DIR/01_scf/CHGCAR" ]] || die "01_scf/CHGCAR is missing."
      cp "$ROOT_DIR/01_scf/POSCAR" "$ROOT_DIR/$stage/POSCAR"
      cp "$ROOT_DIR/01_scf/CHGCAR" "$ROOT_DIR/$stage/CHGCAR"
      if [[ "$stage" == 04_wann ]]; then
        [[ -s "$ROOT_DIR/01_scf/POTCAR" ]] || die "01_scf/POTCAR is missing."
        cp "$ROOT_DIR/01_scf/POTCAR" "$ROOT_DIR/$stage/POTCAR"
      fi
      ;;
    05_crpa) ;;
    *) die "Unknown stage: $stage" ;;
  esac
}

execute_stage() {
  local stage="${1:-}" command
  [[ -n "$stage" && -f "$ROOT_DIR/$stage/$MARKER" ]] ||
    die "Stage '$stage' is not prepared."
  sync_predecessor "$stage"
  command="${STAGE_COMMANDS[$stage]:-${STAGE_COMMANDS[default]}}"
  info "Running $stage with: $command"
  if [[ -n "$EXECUTION_SETUP" ]]; then
    command="$EXECUTION_SETUP"$'\n'"$command"
  fi
  (cd "$ROOT_DIR/$stage" && bash -lc "$command")
}

run_all() {
  local stage
  for stage in 00_relax 01_scf 02_dos 03_band; do
    is_enabled "$stage" || continue
    execute_stage "$stage"
  done
}

run_wannier() {
  is_enabled 04_wann || die "04_wann is not prepared. Run prepare-wannier first."
  execute_stage 04_wann
}

run_crpa() {
  is_enabled 05_crpa || die "05_crpa is not prepared. Run prepare-crpa first."
  execute_stage 05_crpa
}

submit_stage_job() {
  local stage="$1"
  local submit_command="$SUBMIT_COMMAND" submit_directory
  local -a submit_options=(--parsable)
  shift
  if [[ -n "$SUBMIT_JOB_PREFIX" ]]; then
    submit_options+=(--job-name="${SUBMIT_JOB_PREFIX}-${stage#*_}")
  fi
  if [[ "$submit_command" == */* ]]; then
    submit_directory="$(cd "$(dirname "$submit_command")" && pwd)"
    submit_command="$submit_directory/$(basename "$submit_command")"
  fi
  (cd "$ROOT_DIR/$stage" && "$submit_command" "${submit_options[@]}" "$@" job.sh)
}

parse_submit_options() {
  SUBMIT_JOB_PREFIX=""
  while (( $# )); do
    case "$1" in
      --job-name)
        (( $# >= 2 )) && [[ -n "$2" && "$2" != --* ]] ||
          die "--job-name requires a non-empty prefix."
        SUBMIT_JOB_PREFIX="$2"
        shift 2
        ;;
      --job-name=*)
        SUBMIT_JOB_PREFIX="${1#*=}"
        [[ -n "$SUBMIT_JOB_PREFIX" ]] ||
          die "--job-name requires a non-empty prefix."
        shift
        ;;
      *)
        die "Unknown submit option '$1'. Use --job-name PREFIX."
        ;;
    esac
  done
}

submit_all() {
  parse_submit_options "$@"
  command -v "$SUBMIT_COMMAND" >/dev/null 2>&1 ||
    die "Submission command '$SUBMIT_COMMAND' is unavailable."
  is_enabled 01_scf || die "SCF stage is not prepared. Run prepare first."

  local output relax_id scf_id stage job_id
  if is_enabled 00_relax; then
    output="$(submit_stage_job 00_relax)"
    relax_id="${output%%;*}"
    [[ "$relax_id" =~ ^[0-9]+$ ]] || die "Could not parse job ID from: $output"
    info "Submitted 00_relax as job $relax_id"
    output="$(submit_stage_job 01_scf --dependency="afterok:$relax_id")"
  else
    output="$(submit_stage_job 01_scf)"
  fi
  scf_id="${output%%;*}"
  [[ "$scf_id" =~ ^[0-9]+$ ]] || die "Could not parse job ID from: $output"
  info "Submitted 01_scf as job $scf_id"

  # DOS and band calculations are independent children of the SCF calculation.
  for stage in 02_dos 03_band; do
    is_enabled "$stage" || continue
    output="$(submit_stage_job "$stage" --dependency="afterok:$scf_id")"
    job_id="${output%%;*}"
    [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse job ID from: $output"
    info "Submitted $stage as job $job_id"
  done
}

submit_wannier() {
  parse_submit_options "$@"
  command -v "$SUBMIT_COMMAND" >/dev/null 2>&1 ||
    die "Submission command '$SUBMIT_COMMAND' is unavailable."
  is_enabled 04_wann || die "04_wann is not prepared. Run prepare-wannier first."
  [[ -s "$ROOT_DIR/04_wann/job.sh" ]] ||
    die "04_wann/job.sh is missing. Run prepare-wannier again."

  local output job_id
  output="$(submit_stage_job 04_wann)"
  job_id="${output%%;*}"
  [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse job ID from: $output"
  info "Submitted 04_wann as job $job_id"
}

submit_crpa() {
  parse_submit_options "$@"
  command -v "$SUBMIT_COMMAND" >/dev/null 2>&1 ||
    die "Submission command '$SUBMIT_COMMAND' is unavailable."
  is_enabled 05_crpa || die "05_crpa is not prepared. Run prepare-crpa first."
  [[ -s "$ROOT_DIR/05_crpa/job.sh" ]] ||
    die "05_crpa/job.sh is missing. Run prepare-crpa again."

  local output job_id
  output="$(submit_stage_job 05_crpa)"
  job_id="${output%%;*}"
  [[ "$job_id" =~ ^[0-9]+$ ]] || die "Could not parse job ID from: $output"
  info "Submitted 05_crpa as job $job_id"
}

status() {
  local stage state
  for stage in 00_relax 01_scf 02_dos 03_band 04_wann 05_crpa; do
    is_enabled "$stage" || continue
    state="prepared"
    [[ -s "$ROOT_DIR/$stage/OUTCAR" ]] && state="started"
    [[ -s "$ROOT_DIR/$stage/OUTCAR" ]] &&
      grep -q 'General timing and accounting informations for this job' "$ROOT_DIR/$stage/OUTCAR" &&
      state="finished"
    printf '%-10s %s\n' "$stage" "$state"
  done
}

case "${1:-}" in
  prepare) shift; prepare "$@" ;;
  run) run_all ;;
  submit) shift; submit_all "$@" ;;
  execute) execute_stage "${2:-}" ;;
  status) status ;;
  postprocess) shift; exec "$PYTHON_BIN" "$ROOT_DIR/postprocess.py" --root "$ROOT_DIR" "$@" ;;
  prepare-wannier) shift; prepare_wannier_stage "$@" ;;
  run-wannier) run_wannier ;;
  submit-wannier) shift; submit_wannier "$@" ;;
  prepare-crpa) shift; prepare_crpa_stage "$@" ;;
  run-crpa) run_crpa ;;
  submit-crpa) shift; submit_crpa "$@" ;;
  -h|--help|help|"") usage ;;
  *) die "Unknown command '$1'. Run ./workflow.sh --help." ;;
esac
