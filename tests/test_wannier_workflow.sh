#!/usr/bin/env bash
set -euo pipefail

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/version1-wannier-workflow-test.XXXXXX")"
trap 'case "$TEST_DIR" in */version1-wannier-workflow-test.*) rm -rf -- "$TEST_DIR" ;; esac' EXIT

cp "$SOURCE_DIR/workflow.sh" "$TEST_DIR/workflow.sh"
chmod +x "$TEST_DIR/workflow.sh"
mkdir -p "$TEST_DIR/01_scf"
for name in POSCAR POTCAR CHGCAR; do
  printf '%s data\n' "$name" > "$TEST_DIR/01_scf/$name"
done

cat > "$TEST_DIR/prepare_wannier.py" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$@" > "$(dirname "$0")/prepare-arguments"
root=""
while (( $# )); do
  if [[ "$1" == "--root" ]]; then
    root="$2"
    shift 2
  else
    shift
  fi
done
[[ -n "$root" ]]
mkdir -p "$root/04_wann"
: > "$root/04_wann/.generated-by-poscar-workflow"
printf 'mock INCAR\n' > "$root/04_wann/INCAR"
EOF
chmod +x "$TEST_DIR/prepare_wannier.py"

cat > "$TEST_DIR/mock_sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" > "$(dirname "$0")/submit-arguments"
pwd > "$(dirname "$0")/submit-working-directory"
printf '12345\n'
EOF
chmod +x "$TEST_DIR/mock_sbatch"

cat > "$TEST_DIR/workflow.conf" <<EOF
PYTHON_BIN="bash"
VASPKIT_BIN="mock-vaspkit"
SUBMIT_COMMAND="$TEST_DIR/mock_sbatch"
KPR_WANN=0.06
EXECUTION_SETUP='printf "tasks=%s\\n" "\${SLURM_NTASKS:-missing}" > setup-ran'
declare -A STAGE_COMMANDS=(
  [default]='printf "ran 04_wann\\n" > executed'
)
EOF

(
  cd "$TEST_DIR"
  ./workflow.sh prepare-wannier \
    --elements Mn Sb --orbitals d p
)
grep -Fx -- '--kpr' "$TEST_DIR/prepare-arguments"
grep -Fx -- '0.06' "$TEST_DIR/prepare-arguments"
! grep -Fx -- '--num-bands' "$TEST_DIR/prepare-arguments"

(
  cd "$TEST_DIR"
  ./workflow.sh prepare-wannier \
    --elements Mn Sb --orbitals d p --num-bands 22 --kpr 0.05
)

test -s "$TEST_DIR/04_wann/job.sh"
! grep -F 'workflow.sh' "$TEST_DIR/04_wann/job.sh"
! grep -F "$TEST_DIR" "$TEST_DIR/04_wann/job.sh"
grep -F '../01_scf/CHGCAR' "$TEST_DIR/04_wann/job.sh"
grep -F '../01_scf/POTCAR' "$TEST_DIR/04_wann/job.sh"
grep -Fx '04_wann' "$TEST_DIR/.workflow-stages"
grep -Fx -- '--vaspkit' "$TEST_DIR/prepare-arguments"
grep -Fx -- 'mock-vaspkit' "$TEST_DIR/prepare-arguments"
grep -Fx -- '--kpr' "$TEST_DIR/prepare-arguments"
grep -Fx -- '0.05' "$TEST_DIR/prepare-arguments"
! grep -Fx -- '0.06' "$TEST_DIR/prepare-arguments"

status="$("$TEST_DIR/workflow.sh" status)"
grep -Eq '^04_wann[[:space:]]+prepared$' <<< "$status"

"$TEST_DIR/workflow.sh" run
test ! -e "$TEST_DIR/04_wann/executed"

"$TEST_DIR/workflow.sh" run-wannier
grep -Fx 'ran 04_wann' "$TEST_DIR/04_wann/executed"
cmp "$TEST_DIR/01_scf/POSCAR" "$TEST_DIR/04_wann/POSCAR"
cmp "$TEST_DIR/01_scf/POTCAR" "$TEST_DIR/04_wann/POTCAR"
cmp "$TEST_DIR/01_scf/CHGCAR" "$TEST_DIR/04_wann/CHGCAR"

"$TEST_DIR/workflow.sh" submit-wannier
grep -Fx -- '--parsable' "$TEST_DIR/submit-arguments"
grep -Fx 'job.sh' "$TEST_DIR/submit-arguments"
test "$(cat "$TEST_DIR/submit-working-directory")" = "$TEST_DIR/04_wann"

"$TEST_DIR/workflow.sh" submit-wannier --job-name=Pu
grep -Fx -- '--job-name=Pu-wann' "$TEST_DIR/submit-arguments"
grep -Fx 'job.sh' "$TEST_DIR/submit-arguments"

printf 'updated POSCAR\n' > "$TEST_DIR/01_scf/POSCAR"
printf 'updated POTCAR\n' > "$TEST_DIR/01_scf/POTCAR"
printf 'updated CHGCAR\n' > "$TEST_DIR/01_scf/CHGCAR"
printf 'stale\n' > "$TEST_DIR/04_wann/POSCAR"
printf 'stale\n' > "$TEST_DIR/04_wann/POTCAR"
printf 'stale\n' > "$TEST_DIR/04_wann/CHGCAR"
rm -f "$TEST_DIR/04_wann/executed" "$TEST_DIR/04_wann/setup-ran"
mv "$TEST_DIR/workflow.sh" "$TEST_DIR/workflow.sh.saved"
mv "$TEST_DIR/workflow.conf" "$TEST_DIR/workflow.conf.saved"
(
  cd "$TEST_DIR/01_scf"
  SLURM_NTASKS=12 bash ../04_wann/job.sh
)
grep -Fx 'ran 04_wann' "$TEST_DIR/04_wann/executed"
grep -Fx 'tasks=12' "$TEST_DIR/04_wann/setup-ran"
cmp "$TEST_DIR/01_scf/POSCAR" "$TEST_DIR/04_wann/POSCAR"
cmp "$TEST_DIR/01_scf/POTCAR" "$TEST_DIR/04_wann/POTCAR"
cmp "$TEST_DIR/01_scf/CHGCAR" "$TEST_DIR/04_wann/CHGCAR"

printf 'wannier workflow functional test passed\n'
