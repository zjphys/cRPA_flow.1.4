#!/usr/bin/env bash
set -uo pipefail

# Prepare and optionally run/submit one isolated workflow for every POSCAR.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MODE="submit"
FORCE="no"
NO_RELAX="no"
DRY_RUN="no"
STRUCTURE_DIR=""
CALCULATION_DIR=""
MARKER=".batch-workflow-case"

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
info() { printf '==> %s\n' "$*"; }
warn() { printf 'WARNING: %s\n' "$*" >&2; }

usage() {
  cat <<'EOF'
Usage:
  ./batch_workflow.sh [OPTIONS] STRUCTURE_DIR [CALCULATION_DIR]

Create one calculation directory per input POSCAR, copy the version1 workflow
into it, run "workflow.sh prepare", and optionally run or submit the workflow.

Input layouts:
  STRUCTURE_DIR/POSCAR_Fe.vasp       Flat files named POSCAR*, *.vasp, or *.poscar
  STRUCTURE_DIR/Fe/POSCAR            POSCAR files in nested structure directories

Options:
  --mode prepare|submit|run          Default: submit
  --no-relax                         Pass --no-relax to workflow.sh prepare
  --force                            Refresh existing batch-managed calculations
  --dry-run                          Print the discovered mapping without changes
  -h, --help                         Show this help

CALCULATION_DIR defaults to a "calculations" directory next to STRUCTURE_DIR.
Existing calculation directories are skipped unless --force is supplied.
EOF
}

while (( $# > 0 )); do
  case "$1" in
    --mode)
      (( $# >= 2 )) || die "--mode requires a value."
      MODE="$2"
      shift 2
      ;;
    --no-relax)
      NO_RELAX="yes"
      shift
      ;;
    --force)
      FORCE="yes"
      shift
      ;;
    --dry-run)
      DRY_RUN="yes"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while (( $# > 0 )); do
        if [[ -z "$STRUCTURE_DIR" ]]; then
          STRUCTURE_DIR="$1"
        elif [[ -z "$CALCULATION_DIR" ]]; then
          CALCULATION_DIR="$1"
        else
          die "Too many positional arguments."
        fi
        shift
      done
      ;;
    -*)
      die "Unknown option: $1"
      ;;
    *)
      if [[ -z "$STRUCTURE_DIR" ]]; then
        STRUCTURE_DIR="$1"
      elif [[ -z "$CALCULATION_DIR" ]]; then
        CALCULATION_DIR="$1"
      else
        die "Too many positional arguments."
      fi
      shift
      ;;
  esac
done

[[ "$MODE" == "prepare" || "$MODE" == "submit" || "$MODE" == "run" ]] ||
  die "--mode must be prepare, submit, or run."
[[ -n "$STRUCTURE_DIR" ]] || { usage >&2; exit 2; }
[[ -d "$STRUCTURE_DIR" ]] || die "Structure directory does not exist: $STRUCTURE_DIR"
[[ -s "$SCRIPT_DIR/workflow.sh" ]] || die "workflow.sh is missing from $SCRIPT_DIR"

STRUCTURE_DIR="$(cd "$STRUCTURE_DIR" && pwd)"
if [[ -z "$CALCULATION_DIR" ]]; then
  CALCULATION_DIR="$(dirname "$STRUCTURE_DIR")/calculations"
elif [[ "$CALCULATION_DIR" != /* ]]; then
  CALCULATION_DIR="$(pwd)/$CALCULATION_DIR"
fi

declare -a POSCARS=()
while IFS= read -r -d '' poscar; do
  POSCARS+=("$poscar")
done < <(
  find "$STRUCTURE_DIR" \
    \( -path "$CALCULATION_DIR" -prune \) -o \
    \( -type f \
      \( -name 'POSCAR' -o -name 'POSCAR*' -o -iname '*.vasp' -o -iname '*.poscar' \) \
      -print0 \) | sort -z
)
(( ${#POSCARS[@]} > 0 )) ||
  die "No POSCAR, POSCAR*, *.vasp, or *.poscar files found in $STRUCTURE_DIR"

case_name() {
  local input="$1" relative stem parent name
  relative="${input#"$STRUCTURE_DIR"/}"
  stem="$(basename "$relative")"
  parent="$(dirname "$relative")"

  if [[ "$stem" == "POSCAR" && "$parent" != "." ]]; then
    name="${parent//\//__}"
  else
    name="$stem"
    name="${name%.*}"
    name="${name#POSCAR}"
    name="${name#_}"
    name="${name#-}"
    [[ -n "$name" ]] || name="POSCAR"
    if [[ "$parent" != "." ]]; then
      name="${parent//\//__}__${name}"
    fi
  fi

  name="$(printf '%s' "$name" | sed 's/[^A-Za-z0-9._-]/_/g')"
  printf '%s' "$name"
}

declare -A SEEN_NAMES=()
declare -a CASE_NAMES=()
for poscar in "${POSCARS[@]}"; do
  name="$(case_name "$poscar")"
  if [[ -n "${SEEN_NAMES[$name]+x}" ]]; then
    die "Inputs map to the same case name '$name': ${SEEN_NAMES[$name]} and $poscar"
  fi
  SEEN_NAMES["$name"]="$poscar"
  CASE_NAMES+=("$name")
done

info "Found ${#POSCARS[@]} structure(s); mode=$MODE"
for i in "${!POSCARS[@]}"; do
  printf '  %s -> %s\n' "${POSCARS[$i]}" "$CALCULATION_DIR/${CASE_NAMES[$i]}"
done
[[ "$DRY_RUN" == "yes" ]] && exit 0

mkdir -p "$CALCULATION_DIR" || die "Cannot create $CALCULATION_DIR"

copy_workflow_files() {
  local destination="$1" source base
  while IFS= read -r -d '' source; do
    base="$(basename "$source")"
    case "$base" in
      batch_workflow.sh|POSCAR|POTCAR|KPATH.in|.workflow-stages) continue ;;
    esac
    cp "$source" "$destination/$base" || return 1
  done < <(find "$SCRIPT_DIR" -maxdepth 1 -type f -print0)
  chmod +x "$destination/workflow.sh" || return 1
}

prepared=0
completed=0
skipped=0
failed=0

for i in "${!POSCARS[@]}"; do
  poscar="${POSCARS[$i]}"
  name="${CASE_NAMES[$i]}"
  destination="$CALCULATION_DIR/$name"

  if [[ -e "$destination" && ! -d "$destination" ]]; then
    warn "$name: destination exists and is not a directory; skipped."
    ((failed+=1))
    continue
  fi
  if [[ -d "$destination" && ! -f "$destination/$MARKER" ]]; then
    warn "$name: destination was not created by this batch script; skipped."
    ((failed+=1))
    continue
  fi
  if [[ -d "$destination" && "$FORCE" != "yes" ]]; then
    warn "$name: calculation already exists; skipped (use --force to refresh)."
    ((skipped+=1))
    continue
  fi
  if [[ -s "$destination/POSCAR" ]] && ! cmp -s "$poscar" "$destination/POSCAR"; then
    warn "$name: input differs from the existing calculation POSCAR; use a new calculation directory."
    ((failed+=1))
    continue
  fi

  mkdir -p "$destination" || { warn "$name: could not create destination."; ((failed+=1)); continue; }
  : > "$destination/$MARKER"
  if ! copy_workflow_files "$destination"; then
    warn "$name: failed to copy workflow files."
    ((failed+=1))
    continue
  fi
  if ! cp "$poscar" "$destination/POSCAR"; then
    warn "$name: failed to copy POSCAR."
    ((failed+=1))
    continue
  fi

  prepare_args=(prepare)
  [[ "$FORCE" == "yes" ]] && prepare_args+=(--force)
  [[ "$NO_RELAX" == "yes" ]] && prepare_args+=(--no-relax)

  info "$name: preparing"
  if ! (cd "$destination" && bash ./workflow.sh "${prepare_args[@]}"); then
    warn "$name: preparation failed; continuing with the next structure."
    ((failed+=1))
    continue
  fi
  ((prepared+=1))

  if [[ "$MODE" == "prepare" ]]; then
    ((completed+=1))
    continue
  fi

  info "$name: $MODE"
  if (cd "$destination" && bash ./workflow.sh "$MODE"); then
    ((completed+=1))
  else
    warn "$name: $MODE failed; continuing with the next structure."
    ((failed+=1))
  fi
done

printf '\nBatch summary: prepared=%d completed=%d skipped=%d failed=%d\n' \
  "$prepared" "$completed" "$skipped" "$failed"
(( failed == 0 )) || exit 1
