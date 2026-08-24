#!/usr/bin/env bash
set -euo pipefail

SOURCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/version1-crpa-workflow-test.XXXXXX")"
trap 'case "$TEST_DIR" in */version1-crpa-workflow-test.*) rm -rf -- "$TEST_DIR" ;; esac' EXIT

cp "$SOURCE_DIR/workflow.sh" "$TEST_DIR/workflow.sh"
cp "$SOURCE_DIR/prepare_crpa.py" "$TEST_DIR/prepare_crpa.py"
chmod +x "$TEST_DIR/workflow.sh" "$TEST_DIR/prepare_crpa.py"

mkdir -p "$TEST_DIR/04_wann"
: > "$TEST_DIR/04_wann/.generated-by-poscar-workflow"
cat > "$TEST_DIR/04_wann/POSCAR" <<'EOF'
Shell cRPA material
1.0
1 0 0
0 1 0
0 0 1
Mn Se
1 2
Direct
0 0 0
0.25 0.25 0.25
0.75 0.75 0.75
EOF
cat > "$TEST_DIR/04_wann/INCAR" <<'EOF'
SYSTEM = Shell cRPA material Wannier
ENCUT = 500
EDIFF = 1e-6
NBANDS = 72
NUM_WANN = 10
EOF
cat > "$TEST_DIR/04_wann/OUTCAR" <<'EOF'
effective NBANDS = 80
General timing and accounting informations for this job:
EOF
for name in POTCAR KPOINTS CHGCAR WAVECAR; do
  printf '%s restart data\n' "$name" > "$TEST_DIR/04_wann/$name"
done
cat > "$TEST_DIR/04_wann/WANPROJ" <<'EOF'
# ISPIN NKPTS NB_TOT NW
1 1 80 10
EOF

cat > "$TEST_DIR/mock_sbatch" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" > "$(dirname "$0")/submit-arguments"
pwd > "$(dirname "$0")/submit-working-directory"
printf '24680\n'
EOF
chmod +x "$TEST_DIR/mock_sbatch"

cat > "$TEST_DIR/workflow.conf" <<EOF
PYTHON_BIN="python3"
SUBMIT_COMMAND="$TEST_DIR/mock_sbatch"
declare -A STAGE_COMMANDS=(
  [default]='printf "ran 05_crpa\\n" > executed'
)
SBATCH_PARTITION="regular"
SBATCH_NODES=2
SBATCH_NTASKS_PER_NODE=4
SBATCH_CPUS_PER_TASK=1
SBATCH_TIME="01:00:00"
SBATCH_EXTRA=""
CRPA_SBATCH_PARTITION="crpa"
CRPA_SBATCH_NODES=8
CRPA_SBATCH_NTASKS_PER_NODE=12
CRPA_SBATCH_CPUS_PER_TASK=2
CRPA_SBATCH_TIME="36:00:00"
CRPA_SBATCH_EXTRA="#SBATCH --exclusive"
CRPA_KPAR=7
EOF

(
  cd "$TEST_DIR"
  ./workflow.sh prepare-crpa --nbandsgw 70 --encutgw 300
)

test -s "$TEST_DIR/05_crpa/job.sh"
! grep -F 'workflow.sh' "$TEST_DIR/05_crpa/job.sh"
! grep -F "$TEST_DIR" "$TEST_DIR/05_crpa/job.sh"
! grep -F 'require_input "$job_dir/../' "$TEST_DIR/05_crpa/job.sh"
test -s "$TEST_DIR/05_crpa/WANPROJ"
test ! -e "$TEST_DIR/05_crpa/OUTCAR"
grep -Fx '05_crpa' "$TEST_DIR/.workflow-stages"
grep -Fx 'NBANDS = 80' "$TEST_DIR/05_crpa/INCAR"
grep -Fx 'NBANDSGW = 70' "$TEST_DIR/05_crpa/INCAR"
! grep -Eq '^[[:space:]]*ISPIN[[:space:]]*=' "$TEST_DIR/05_crpa/INCAR"
grep -Fx 'ENCUTGW = 300' "$TEST_DIR/05_crpa/INCAR"
grep -Fx 'NTARGET_STATES = 1 2 3 4 5 6 7 8 9 10' "$TEST_DIR/05_crpa/INCAR"
grep -Fx 'KPAR = 7' "$TEST_DIR/05_crpa/INCAR"
grep -Fx '#SBATCH --partition=crpa' "$TEST_DIR/05_crpa/job.sh"
grep -Fx '#SBATCH --nodes=8' "$TEST_DIR/05_crpa/job.sh"
grep -Fx '#SBATCH --ntasks-per-node=12' "$TEST_DIR/05_crpa/job.sh"
grep -Fx '#SBATCH --cpus-per-task=2' "$TEST_DIR/05_crpa/job.sh"
grep -Fx '#SBATCH --time=36:00:00' "$TEST_DIR/05_crpa/job.sh"
grep -Fx '#SBATCH --exclusive' "$TEST_DIR/05_crpa/job.sh"

status="$("$TEST_DIR/workflow.sh" status)"
grep -Eq '^05_crpa[[:space:]]+prepared$' <<< "$status"

"$TEST_DIR/workflow.sh" run
test ! -e "$TEST_DIR/05_crpa/executed"

"$TEST_DIR/workflow.sh" run-crpa
grep -Fx 'ran 05_crpa' "$TEST_DIR/05_crpa/executed"

"$TEST_DIR/workflow.sh" submit-crpa
grep -Fx -- '--parsable' "$TEST_DIR/submit-arguments"
grep -Fx 'job.sh' "$TEST_DIR/submit-arguments"
test "$(cat "$TEST_DIR/submit-working-directory")" = "$TEST_DIR/05_crpa"
! grep -q -- '--dependency' "$TEST_DIR/submit-arguments"

"$TEST_DIR/workflow.sh" submit-crpa --job-name Pu
grep -Fx -- '--job-name=Pu-crpa' "$TEST_DIR/submit-arguments"
grep -Fx 'job.sh' "$TEST_DIR/submit-arguments"
! grep -q -- '--dependency' "$TEST_DIR/submit-arguments"

rm -f "$TEST_DIR/05_crpa/executed"
mv "$TEST_DIR/workflow.sh" "$TEST_DIR/workflow.sh.saved"
mv "$TEST_DIR/workflow.conf" "$TEST_DIR/workflow.conf.saved"
(
  cd "$TEST_DIR/04_wann"
  bash ../05_crpa/job.sh
)
grep -Fx 'ran 05_crpa' "$TEST_DIR/05_crpa/executed"

printf 'cRPA workflow functional test passed\n'
