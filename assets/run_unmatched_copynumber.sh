#!/bin/bash
#BSUB -q oversubscribed
#BSUB -G team113-grp
#BSUB -R "select[mem>8000] rusage[mem=8000] span[hosts=1]"
#BSUB -M 8000

################################################################################
# run_unmatched_copynumber.sh - submit the Dermatlas unmatched copy-number nextflow pipeline
#
# Sections ([edit] = user-editable, [skip] = owned by the script):
#   [skip] FUNCTIONS
#   [skip] FAILURE REPORTING ............. exit trap for a failed launch
#   [skip] RESERVED VARIABLES ............ the env-var classes checked below
#   [edit] ENVIRONMENT SETUP ............. SOURCE_ME: a source_me.sh, or "none"
#   [edit] OPT-IN REPORTING .............. website / Slack / work-dir toggles
#                                          (defaults; the environment overrides)
#   [edit] MANUAL ENVIRONMENT OVERRIDES .. commented exports, one per variable
#   [skip] ENVIRONMENT VALIDATION
#   [edit] RUN CONFIGURATION ............. CONFIG, REVISION, LABEL
#   [skip] FILE SYSTEM SETUP ............. RUN_ID, the .lock, log/trace/stats/clone/cache paths
#   [skip] EXECUTION OF THE PIPELINE
#
# Usage 1 - managed (projectify) runs. The dermanager-generated source_me.sh
#   exports everything. Edit nothing; optionally flip the OPT-IN REPORTING
#   toggles.
#     cd <project_dir> && bsub ... < commands/<pipeline>/run_unmatched_copynumber.sh
#
# Usage 2 - git-clone runs (no dermanager source_me.sh). Provide the
#   environment either way:
#     a) set SOURCE_ME="none", then uncomment and fill in MANUAL ENVIRONMENT
#        OVERRIDES; or
#     b) write your own source_me.sh - copy the MANUAL ENVIRONMENT OVERRIDES
#        block into a file, uncomment, fill in - and point SOURCE_ME at it.
#   Then point CONFIG at your own unmatched_copynumber.config (typically with
#   your own subcohorts), and set the OPT-IN REPORTING toggles to "false" unless
#   you know your cohort's slug, sample-list version and API endpoint.
################################################################################

# -E: the ERR trap must be inherited by functions.
set -Eeuo pipefail

###################
#### FUNCTIONS ####
###################

function require_env() {
  local remedy="$1"
  shift
  local var missing=()
  for var in "$@"; do
    [[ -n "${!var:-}" ]] || missing+=("${var}")
  done
  if (( ${#missing[@]} )); then
    printf 'ERROR: the environment did not provide: %s\n' "${missing[*]}" >&2
    printf '%s\n' "${remedy}" >&2
    exit 1
  fi
}

function normalize_bool() {
  # An unrecognised value is fatal: a typo must not silently pick a side.
  local var="$1" val="${!1:-}"
  case "${val,,}" in
    true|t|yes|y|on|1)   printf -v "${var}" 'true' ;;
    false|f|no|n|off|0)  printf -v "${var}" 'false' ;;
    *)
      printf 'ERROR: %s must be a true/false value (got: "%s").\n' "${var}" "${val}" >&2
      printf 'Accepted: true/false, yes/no, on/off, 1/0 (any case). Set it in the shell, in\n' >&2
      printf 'source_me.sh, or edit the OPT-IN REPORTING section of this script.\n' >&2
      exit 1
      ;;
  esac
}

function check_for_source_me() {
  local source_me="${SOURCE_ME:-./source_me.sh}"
  if [[ ! -f "${source_me}" ]]; then
    printf 'ERROR: %s not found (working directory: %s).\n' "${source_me}" "${PWD}" >&2
    printf 'Submit from the project directory, e.g.\n' >&2
    printf '  cd <project_dir> && bsub ... < commands/run_unmatched_copynumber.sh\n' >&2
    printf '  (or go to the https://team113.sanger.ac.uk/dermatlas/cohorts/ page and re-create the run command for this pipeline for this cohort)\n' >&2
    exit 1
  fi
}

function sanitize() {
  local input="${1:-}"
  # printf, not echo: echo's trailing newline would become a trailing hyphen.
  printf '%s' "${input}" | tr '[:upper:]' '[:lower:]' | tr '[:space:]' '-' | tr -cd '[:alnum:]-_'
}

function truncate_string() {
  # Named for the coreutils `truncate` binary it must not shadow.
  local input="${1:-}"
  local max_length="${2:-}"
  if (( ${#input} > max_length )); then
    echo "${input:0:max_length}"
  else
    echo "${input}"
  fi
}

function set_run_id_from_label() {
  local label="${1:-}"
  local timestamp
  if [[ -z "${label}" ]]; then
    printf 'ERROR: set_run_id_from_label requires a label argument\n' >&2
    exit 1
  fi
  timestamp=$(date +'%Y%m%dT%H%M%S')
  local sanitized_label
  sanitized_label=$(sanitize "${label}")
  echo "${sanitized_label}_${timestamp}"
}

function set_run_id() {
  local study="${1:-}"
  local project="${2:-}"
  local cohort="${3:-}"
  local timestamp
  timestamp=$(date +'%Y%m%dT%H%M%S')
  if [[ -n "${study}" ]]; then
    study="${study}_"
  fi
  if [[ -n "${project}" ]]; then
    project="${project}_"
  fi
  if [[ -n "${cohort}" ]]; then
    local sanitized_cohort
    sanitized_cohort=$(sanitize "${cohort}")
    cohort="$(truncate_string "${sanitized_cohort}" 40)_"
  fi
  echo "${study}${project}${cohort}${timestamp}"
}

function json_escape() {
  local s="${1:-}"
  s="${s//\\/\\\\}"
  s="${s//\"/\\\"}"
  s="${s//$'\r'/}"
  s="${s//$'\t'/\\t}"
  s="${s//$'\n'/\\n}"
  printf '%s' "${s}"
}

function launcher_failure_ref() {
  # Mirrors the onComplete handler's run reference. Never angle brackets: this
  # reaches Slack, which parses <...> as a link.
  local ref="${COHORT_SLUG:-}"
  [[ -n "${ref}" ]] || ref="${LABEL:-}"
  [[ -n "${ref}" ]] || ref="unknown-cohort"
  printf '%s' "${ref}"
}

function launcher_failure_details() {
  # The trap can fire at any point, so each value is printed only once knowable.
  local status="${1:-}"
  printf 'Cohort: %s\n' "$(launcher_failure_ref)"
  printf 'Study: %s\n' "${STUDY:-unset}"
  printf 'Project: %s\n' "${PROJECT:-unset}"
  printf 'Pipeline: %s\n' "${PIPELINE_SLUG:-unset}${REVISION:+ (${REVISION})}"
  printf 'Exit status: %s\n' "${status}"
  if [[ -n "${_LAUNCHER_EXTRA_NOTE:-}" ]]; then
    printf '%s\n' "${_LAUNCHER_EXTRA_NOTE}"
  fi
  if [[ -n "${_LAST_ERR_CMD:-}" ]]; then
    printf 'Failed command (line %s): %s\n' "${_LAST_ERR_LINE:-?}" "${_LAST_ERR_CMD}"
  fi
  if [[ -n "${RUN_ID:-}" ]]; then
    printf 'Run id: %s\n' "${RUN_ID}"
  fi
  if [[ -n "${NXF_PULL_LOG_FILE:-}" && -f "${NXF_PULL_LOG_FILE}" ]]; then
    printf 'Nextflow pull log: %s\n' "${NXF_PULL_LOG_FILE}"
  fi
  # Often the only identifying values when source_me.sh never loaded.
  if [[ -n "${LSB_JOBID:-}" ]]; then
    printf 'LSF job: %s\n' "${LSB_JOBID}${LSB_JOBNAME:+ (${LSB_JOBNAME})}"
  fi
  if [[ -n "${LS_SUBCWD:-}" ]]; then
    printf 'Submitted from: %s\n' "${LS_SUBCWD}"
  fi
  printf 'Host: %s\n' "$(hostname -s 2>/dev/null || echo unknown)"
}

function launcher_failure_slack_message() {
  # Single-quoted formats: the backticks are Slack markup, not substitution.
  local status="${1:-}"
  printf ':octagonal_sign: *%s* failed before the pipeline was submitted - `%s`\n' \
         "${_LAUNCHER_LABEL}" "$(launcher_failure_ref)"
  printf '_Study: %s | Project: %s | Pipeline: %s_\n' \
         "${STUDY:-unset}" "${PROJECT:-unset}" "${PIPELINE_SLUG:-unset}"
  printf '```\n%s\n```' "$(launcher_failure_details "${status}")"
}

function report_launcher_failure() {
  # stderr unconditionally - it lands in the LSF job output, so a batch of
  # submissions never fails silently. Slack only if opted in and possible.
  local status="${1:-}" msg
  printf '\n' >&2
  printf 'ERROR: %s failed - the pipeline was never submitted.\n' "${_LAUNCHER_LABEL}" >&2
  launcher_failure_details "${status}" >&2

  if (( ${_TRAP_CAN_SLACK:-0} != 1 )); then
    printf 'NOTE: failed before the reporting opt-ins were read; no Slack message sent.\n' >&2
    return 0
  fi
  if [[ "${DERMATLAS_SLACK_NOTIFICATIONS:-false}" != "true" ]]; then
    printf 'NOTE: DERMATLAS_SLACK_NOTIFICATIONS is not "true"; no Slack message sent.\n' >&2
    return 0
  fi
  if [[ -z "${SLACK_WEBHOOK_URL:-}" ]]; then
    printf 'NOTE: SLACK_WEBHOOK_URL is unset; no Slack message sent.\n' >&2
    return 0
  fi
  if ! command -v curl >/dev/null 2>&1; then
    printf 'NOTE: curl was not found on PATH; no Slack message sent.\n' >&2
    return 0
  fi
  msg="$(launcher_failure_slack_message "${status}")"
  if curl -sS -X POST --max-time 15 \
          -H 'Content-Type: application/json; charset=utf-8' \
          --data "{\"text\":\"$(json_escape "${msg}")\"}" \
          "${SLACK_WEBHOOK_URL}" >/dev/null 2>&1; then
    printf 'Slack failure notification sent.\n' >&2
  else
    printf 'NOTE: the Slack failure notification could not be sent.\n' >&2
  fi
  return 0
}

function pipeline_dir_flock_scope() {
  # Would a flock here be honoured by other nodes? Never fatal: the process that
  # deletes data must fail closed on this, not the launcher.
  local dir="${1:-}" line fstype opts
  command -v findmnt >/dev/null 2>&1 || { printf 'unknown'; return 0; }
  line="$(findmnt -T "${dir}" -no FSTYPE,OPTIONS 2>/dev/null)" || { printf 'unknown'; return 0; }
  [[ -n "${line}" ]] || { printf 'unknown'; return 0; }
  read -r fstype opts <<< "${line}"
  # These keep flock inside one node, so two farm nodes would each think they
  # own this directory.
  case ",${opts}," in
    *,localflock,*|*,noflock,*|*,local_lock=all,*|*,local_lock=flock,*)
      printf 'node-local'; return 0 ;;
  esac
  case "${fstype}" in
    lustre|nfs|nfs4|gpfs|ceph|cephfs|beegfs|glusterfs) printf 'cluster' ;;
    *)                                                 printf 'node-local' ;;
  esac
  return 0
}

function lock_holder_description() {
  # A blank file is legitimate: the holder rewrites .lock in the instant after
  # acquiring. Angle brackets become parens for Slack.
  local desc=""
  if [[ -f "${LOCK_FILE:-}" ]]; then
    desc="$(head -c 4096 -- "${LOCK_FILE}" 2>/dev/null \
            | tr -c '[:print:]' ' ' | tr '<>' '()' | tr -s ' ')"
  fi
  if [[ "${desc}" =~ [^[:space:]] ]]; then
    printf '%s' "${desc}"
  else
    printf 'unknown (the holder had not yet recorded its identity)'
  fi
}

function acquire_pipeline_lock() {
  # The lock lives on the open file description, so the kernel drops it when this
  # process dies, however it dies - there is no release, and no stale lock. Nothing may
  # unlink ${LOCK_FILE}: the next run would create a fresh inode and lock that instead,
  # excluding nobody. ">>" not ">", which truncates on open and would destroy a live
  # holder's identity line before flock reports the loss. The fd is inherited by
  # `nextflow run` on purpose, so an orphaned nextflow still holds the directory.
  local scope rc=0 held_ino path_ino

  if [[ -L "${LOCK_FILE}" ]]; then
    printf 'ERROR: %s is a symlink; refusing to lock or write through it.\n' "${LOCK_FILE}" >&2
    exit 1
  fi
  if [[ -e "${LOCK_FILE}" && ! -f "${LOCK_FILE}" ]]; then
    printf 'ERROR: %s exists and is not a regular file.\n' "${LOCK_FILE}" >&2
    exit 1
  fi

  scope="$(pipeline_dir_flock_scope "${PIPELINE_DIR}")"
  if [[ "${scope}" == "node-local" ]]; then
    printf 'WARNING: flock on %s looks node-local (%s).\n' "${PIPELINE_DIR}" \
           "$(findmnt -T "${PIPELINE_DIR}" -no FSTYPE,OPTIONS 2>/dev/null || echo 'unknown mount')" >&2
    printf '         Runs on other nodes will NOT exclude each other, and a work-dir\n' >&2
    printf '         pruner running elsewhere will not see this lock. Continuing anyway.\n' >&2
  fi

  # Unchecked: a failed `exec` redirection exits the shell outright, and
  # _HOLDS_PIPELINE_LOCK is still 0, so no sentinel is written.
  exec {_LOCK_FD}>>"${LOCK_FILE}"
  chmod 0644 "${LOCK_FILE}" 2>/dev/null || true

  flock -n -E "${_LOCK_CONFLICT_RC}" -x "${_LOCK_FD}" || rc=$?
  if (( rc == _LOCK_CONFLICT_RC )); then
    _LAUNCHER_EXTRA_NOTE="Lock holder: $(lock_holder_description)"
    printf 'ERROR: another run already owns %s.\n' "${PIPELINE_DIR}" >&2
    printf '       %s\n' "${_LAUNCHER_EXTRA_NOTE}" >&2
    printf '       Lock file: %s\n' "${LOCK_FILE}" >&2
    printf '       Wait for that run to finish, or kill it, before submitting again.\n' >&2
    exit "${_LOCK_CONFLICT_RC}"
  elif (( rc != 0 )); then
    # flock's own sysexits codes: reporting these as contention would send someone
    # hunting for a run that does not exist.
    printf 'ERROR: flock failed on %s (exit %s). This is NOT contention - flock may be\n' "${LOCK_FILE}" "${rc}" >&2
    printf '       unsupported on this filesystem, or the descriptor is unusable.\n' >&2
    exit 1
  fi

  # If .lock was replaced between our open() and our flock() we hold an orphaned
  # inode and exclude nobody.
  held_ino="$(stat -Lc '%d:%i' "/proc/self/fd/${_LOCK_FD}" 2>/dev/null || true)"
  path_ino="$(stat -c '%d:%i' "${LOCK_FILE}" 2>/dev/null || true)"
  if [[ -n "${held_ino}" && "${held_ino}" != "${path_ino}" ]]; then
    printf 'ERROR: %s was replaced or removed while it was being locked.\n' "${LOCK_FILE}" >&2
    printf '       Nothing may unlink this file: doing so breaks mutual exclusion.\n' >&2
    exit 1
  fi

  # Only now: every exit above leaves the lock unheld.
  _HOLDS_PIPELINE_LOCK=1

  # A separate, truncating open: the fd above is O_APPEND and cannot be rewritten
  # through. Under the lock, so two runs cannot interleave here.
  printf 'run_id=%s pid=%s host=%s lsf_job=%s revision=%s flock_scope=%s started=%s\n' \
         "${RUN_ID:-unset}" "$$" "$(hostname -s 2>/dev/null || echo unknown)" \
         "${LSB_JOBID:-none}${LSB_JOBNAME:+ (${LSB_JOBNAME})}" "${REVISION:-unset}" \
         "${scope}" "$(date +'%Y-%m-%dT%H:%M:%S%z')" \
         > "${LOCK_FILE}" \
    || printf 'NOTE: could not record the lock holder in %s.\n' "${LOCK_FILE}" >&2

  # From here until the exit trap neither sentinel exists, which a pruner must read
  # as "in progress, or died without a verdict - do not touch".
  rm -f -- "${PIPELINE_DIR}/.completed_successfully" \
           "${PIPELINE_DIR}/.completed_with_error"

  printf 'Pipeline directory locked: %s\n' "${LOCK_FILE}"
  return 0
}

function write_completion_sentinel() {
  # The external work-dir pruner's input. The filename is the state; the contents are
  # an audit line for humans that no pruning decision should parse. Neither file exists
  # while a run is in progress, or after one that died without reaching a trap.
  #
  # The lock guard is load-bearing: a run that lost the race to flock reaches its trap
  # too, and writing .completed_with_error there would declare the live run failed and
  # invite the pruner to delete its work directory.
  local status="${1:-0}" name path outcome
  (( ${_HOLDS_PIPELINE_LOCK:-0} == 1 )) || return 0
  [[ -n "${PIPELINE_DIR:-}" && -d "${PIPELINE_DIR:-}" ]] || return 0

  if (( status == 0 )); then
    name=".completed_successfully"; outcome="success"
  elif (( status >= 128 )); then
    name=".completed_with_error";   outcome="killed"   # 128+n: bkill, MEMLIMIT, RUNLIMIT
  else
    name=".completed_with_error";   outcome="failed"
  fi
  path="${PIPELINE_DIR}/${name}"

  # Both: exactly one verdict may exist, and rm -f drops a symlink planted where the
  # sentinel goes rather than writing through it.
  rm -f -- "${PIPELINE_DIR}/.completed_successfully" \
           "${PIPELINE_DIR}/.completed_with_error"

  # One write, so a reader sees the file absent or whole. Still under the lock, since
  # the fd closes only when this process exits.
  if printf 'outcome=%s\nexit_status=%s\nrun_id=%s\nrevision=%s\nfinished=%s\nlsf_job_id=%s\nlsf_job_name=%s\nhost=%s\nwork_dir=%s\nwork_dir_disposition=%s\n' \
       "${outcome}" "${status}" "${RUN_ID:-unset}" "${REVISION:-unset}" \
       "$(date +'%Y-%m-%dT%H:%M:%S%z' 2>/dev/null || echo unknown)" \
       "${LSB_JOBID:-none}" "${LSB_JOBNAME:-none}" \
       "$(hostname -s 2>/dev/null || echo unknown)" \
       "${NXF_WORK:-unset}" "${_WORK_DIR_DISPOSITION:-kept}" \
       > "${path}"
  then
    chmod 0644 "${path}" 2>/dev/null || true
    printf 'Wrote completion sentinel: %s\n' "${path}"
  else
    printf 'NOTE: could not write the completion sentinel %s.\n' "${path}" >&2
  fi
  return 0
}

function on_launcher_exit() {
  # Re-emits ${status} unchanged so the LSF job still fails.
  local status="${1:-0}"
  trap - EXIT ERR INT TERM HUP     # never re-enter
  set +e                           # a failure here must not mask ${status}
  set +u                           # under `set -u` an unset ref kills the shell
                                   # even with `set +e`
  # Sentinel first: reporting can hang on curl or be cut short by a second signal.
  write_completion_sentinel "${status}"
  if (( status != 0 )); then
    report_launcher_failure "${status}"
  fi
  exit "${status}"
}

function work_dir_is_ours() {
  # An inherited NXF_WORK belongs to someone else.
  [[ -n "${PIPELINE_DIR:-}" && "${NXF_WORK:-}" == "${PIPELINE_DIR}/work" ]]
}

function cleanup_work_dir() {
  # Safe only because publishDir has already copied out the keepers (which requires
  # params.publish_dir_mode = 'copy'). Callers must do this on success only.
  local work="${NXF_WORK:-}"
  if ! work_dir_is_ours; then
    printf 'NOTE: refusing to delete "%s" - it is not the work directory this script created.\n' "${work}" >&2
    _WORK_DIR_DISPOSITION="not-removed-refused"
    return 0
  fi
  if [[ ! -d "${work}" ]]; then
    _WORK_DIR_DISPOSITION="absent"
    return 0
  fi
  printf 'Cleaning up work directory: %s\n' "${work}"
  if rm -rf "${work}"; then
    printf 'Work directory removed.\n'
    _WORK_DIR_DISPOSITION="removed"
  else
    printf 'NOTE: the work directory was not fully removed: %s\n' "${work}" >&2
    _WORK_DIR_DISPOSITION="partially-removed"
  fi
  return 0
}

function measure_work_dir() {
  # One walk, not `du -s` plus `du -s --inodes`: on Lustre the traversal is the whole
  # cost. Metadata only, so there is no block I/O for ionice to throttle - and Lustre
  # RPCs bypass the local I/O scheduler regardless.
  #
  # %b is st_blocks: allocated bytes, not apparent size. Symlinks are not followed, so
  # a stage-in symlink costs its own inode and not the BAM behind it. Hardlinks are
  # counted per link - an over-estimate, close enough for costing.
  local work="${NXF_WORK:-}" out rc=0
  work_dir_is_ours || return 0
  [[ -d "${work}" ]] || return 0

  printf 'Measuring work directory: %s\n' "${work}"
  # pipefail inside the substitution: PIPESTATUS outside would describe the
  # assignment, and a find that died mid-walk leaves awk printing a partial total
  # that reads exactly like a real one.
  out="$(set -o pipefail
         find "${work}" -xdev -printf '%b\n' 2>/dev/null \
           | awk '{ blocks += $1; files++ } END { printf "%.0f %d\n", blocks * 512, files }')" \
    || rc=$?
  if (( rc != 0 )) || [[ -z "${out}" ]]; then
    printf 'NOTE: could not measure %s (exit %s); the disk figures will be omitted.\n' \
           "${work}" "${rc}" >&2
    return 0
  fi
  read -r _WORK_DIR_BYTES _WORK_DIR_INODES <<< "${out}"
  return 0
}

function pipeline_wall_seconds() {
  # Empty if `nextflow run` never started.
  local now delta
  [[ -n "${_PIPELINE_START_EPOCH:-}" ]] || return 0
  now="$(date +%s 2>/dev/null || true)"
  [[ -n "${now}" ]] || return 0
  delta=$(( now - _PIPELINE_START_EPOCH ))
  if (( delta < 0 )); then delta=0; fi    # a multi-day run can outlive an NTP step
  printf '%s' "${delta}"
}

function human_duration() {
  local total="${1:-0}"
  printf '%dh %02dm %02ds' "$(( total / 3600 ))" "$(( (total % 3600) / 60 ))" "$(( total % 60 ))"
}

function write_resource_stats() {
  # Written on success only, and before the cleanup: after the rm -rf there is nothing
  # left to measure, whereas a failed run keeps its work directory.
  #
  # Conversions are '#' comments on their own line, so `value="${line#*=}"` stays a
  # valid way to read this. An unmeasurable figure is omitted rather than written as
  # 0 - a missing key cannot be silently summed into a cost estimate.
  local gib_whole gib_frac k_whole k_frac
  [[ -n "${STATS_FILE:-}" ]] || return 0

  if ! {
    # First, so a concatenated record is attributable without its filename.
    if [[ -n "${RUN_ID:-}" ]]; then
      printf 'run_id=%s\n' "${RUN_ID}"
    fi
    if [[ -n "${REVISION:-}" ]]; then
      printf 'revision=%s\n' "${REVISION}"
    fi
    # Read by `dermatlas-http cohort analysis-workdir-stats`, which cannot infer
    # the pipeline from a relative path.
    if [[ -n "${PIPELINE_SLUG:-}" ]]; then
      printf 'pipeline_slug=%s\n' "${PIPELINE_SLUG}"
    fi
    if [[ -n "${_PIPELINE_WALL_SECONDS:-}" ]]; then
      printf '# wall_time: %s (`nextflow run` only)\n' "$(human_duration "${_PIPELINE_WALL_SECONDS}")"
      printf 'wall_time=%s\n' "${_PIPELINE_WALL_SECONDS}"
    fi
    if [[ -n "${_WORK_DIR_BYTES:-}" ]]; then
      gib_whole=$(( _WORK_DIR_BYTES / 1073741824 ))
      gib_frac=$(( (_WORK_DIR_BYTES % 1073741824) * 100 / 1073741824 ))
      printf '# disk_usage: %d.%02d GiB allocated\n' "${gib_whole}" "${gib_frac}"
      printf 'disk_usage=%s\n' "${_WORK_DIR_BYTES}"
    fi
    if [[ -n "${_WORK_DIR_INODES:-}" ]]; then
      k_whole=$(( _WORK_DIR_INODES / 1000 ))
      k_frac=$(( (_WORK_DIR_INODES % 1000) / 100 ))
      printf '# disk_inodes: %d.%d thousand\n' "${k_whole}" "${k_frac}"
      printf 'disk_inodes=%s\n' "${_WORK_DIR_INODES}"
    fi
  } > "${STATS_FILE}"
  then
    printf 'NOTE: could not write the resource stats %s.\n' "${STATS_FILE}" >&2
    return 0
  fi

  chmod 0644 "${STATS_FILE}" 2>/dev/null || true
  printf 'Wrote resource stats: %s\n' "${STATS_FILE}"
  return 0
}

function report_workdir_usage() {
  # The second half of website logging. workflow.onComplete POSTed the run
  # (the row with its status and duration); this PATCHes onto that row the
  # work-dir figures the launcher could only measure after `nextflow run`
  # returned. Fed the stats file, so the numbers on the website are the numbers
  # on disk; the server finds the row from run_id + wall_time.
  #
  # Best effort, always returns 0: the run itself succeeded, and the stats file
  # remains on disk to replay by hand with the command printed on failure.
  # dermatlas-http is not on PATH in an LSF job (Utils.groovy solves the same
  # problem for onComplete), so it is module-loaded here when needed.
  local rc=0
  [[ "${DERMATLAS_WEBSITE_LOGGING:-false}" == "true" ]] || return 0
  if [[ -z "${STATS_FILE:-}" || ! -s "${STATS_FILE}" ]]; then
    printf 'NOTE: no resource stats file; skipping the website work-dir report.\n' >&2
    return 0
  fi
  if [[ -z "${_WORK_DIR_BYTES:-}" || -z "${_WORK_DIR_INODES:-}" ]]; then
    printf 'NOTE: work directory was not measured; skipping the website work-dir report.\n' >&2
    return 0
  fi
  if ! command -v dermatlas-http >/dev/null 2>&1; then
    if ! type module >/dev/null 2>&1; then
      printf 'NOTE: dermatlas-http is not on PATH and `module` is unavailable; skipping the website work-dir report.\n' >&2
      return 0
    fi
    if ! module load "${DERMATLAS_HTTP_MODULE}"; then
      printf 'NOTE: `module load %s` failed; skipping the website work-dir report.\n' "${DERMATLAS_HTTP_MODULE}" >&2
      return 0
    fi
    if ! command -v dermatlas-http >/dev/null 2>&1; then
      printf 'NOTE: module %s loaded but dermatlas-http is still not on PATH; skipping the website work-dir report.\n' "${DERMATLAS_HTTP_MODULE}" >&2
      return 0
    fi
  fi
  printf 'Reporting work-dir usage to the Dermatlas website: %s\n' "${STATS_FILE}"
  dermatlas-http cohort analysis-workdir-stats \
    --stats-file "${STATS_FILE}" \
    --pipeline "${PIPELINE_SLUG}" \
    --api "${SELF_DESCRIBING_API}" || rc=$?
  if (( rc != 0 )); then
    printf 'NOTE: the website work-dir report failed (exit %s). Replay it with:\n' "${rc}" >&2
    printf '      dermatlas-http cohort analysis-workdir-stats --stats-file %q --pipeline %q --api %q\n' \
           "${STATS_FILE}" "${PIPELINE_SLUG}" "${SELF_DESCRIBING_API}" >&2
  fi
  return 0
}

function on_pipeline_exit() {
  # Reports nothing: by now that is workflow.onComplete's job.
  local status="${1:-0}"
  trap - EXIT INT TERM HUP         # never re-enter
  set +e                           # cleanup must not mask ${status}
  set +u                           # see on_launcher_exit
  # Before the walk below, which is not run time.
  _PIPELINE_WALL_SECONDS="$(pipeline_wall_seconds)"
  if (( status == 0 )); then
    # Before the cleanup: the last moment work/ is guaranteed to exist.
    measure_work_dir
    write_resource_stats
    report_workdir_usage
    if [[ "${DERMATLAS_CLEANUP_WORK_DIR:-true}" == "true" ]]; then
      cleanup_work_dir
    else
      _WORK_DIR_DISPOSITION="kept-opted-out"
      printf 'Keeping work directory (DERMATLAS_CLEANUP_WORK_DIR=false): %s\n' "${NXF_WORK:-unset}"
    fi
  else
    _WORK_DIR_DISPOSITION="kept-failed-run"
  fi
  # After the cleanup, never before: it must not claim "removed" for a run killed
  # part-way through its own rm -rf.
  write_completion_sentinel "${status}"
  exit "${status}"
}

############################
#### FAILURE REPORTING  ####
############################

# Covers everything up to the `trap -` before `nextflow run`. Without it a setup
# failure is silent: `nextflow run` never starts, so the pipeline's own onComplete
# reporting never runs. Two stages, because SLACK_WEBHOOK_URL is not known yet -
# until OPT-IN REPORTING sets _TRAP_CAN_SLACK=1 the trap reports to stderr only.
_LAUNCHER_LABEL="Dermatlas unmatched copy number launcher"
_TRAP_CAN_SLACK=0
_LAST_ERR_CMD=""
_LAST_ERR_LINE=""

# ERR only records context; EXIT reports, so a bare `exit 1` is reported too. The
# signal traps exit 128+n so a killed job still reaches the EXIT trap.
trap '_LAST_ERR_CMD="${BASH_COMMAND}"; _LAST_ERR_LINE="${LINENO}"' ERR
trap 'on_launcher_exit $?' EXIT
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

############################
#### RESERVED VARIABLES ####
############################
_DEFAULT_PIPELINE_SLUG="unmatched_copy_number_pipe"
_DEFAULT_SOURCE_ME="./source_me.sh"
# The environment module that puts dermatlas-http on PATH for the work-dir report.
_DEFAULT_DERMATLAS_HTTP_MODULE="dermatlas-http"
# Snapshotted before source_me.sh is read, so OPT-IN REPORTING can apply the
# submitting shell's values over it.
_ENV_WEBSITE_LOGGING="${DERMATLAS_WEBSITE_LOGGING:-}"
_ENV_SLACK_NOTIFICATIONS="${DERMATLAS_SLACK_NOTIFICATIONS:-}"
_ENV_CLEANUP_WORK_DIR="${DERMATLAS_CLEANUP_WORK_DIR:-}"
# Checked after sourcing source_me.sh; the last two only when their toggle is true.
_PIPELINE_ENV_VARS=(PROJECT_DIR COMMANDS_DIR ANALYSIS_DIR BAMS_DIR STUDY PROJECT COHORT \
                DNA_TUMOUR_LIST_ANALYSED_ALL DNA_TUMOUR_LIST_ONE_TUMOUR_PER_PATIENT_ALL \
                DNA_TUMOUR_LIST_RELATED_TUMOURS_ALL COHORT_METADATA_FILE)
_WEBSITE_ENV_VARS=(COHORT_SLUG SAMPLE_LIST_VERSION_FILE SELF_DESCRIBING_API)
_SLACK_ENV_VARS=(SLACK_WEBHOOK_URL)
PIPELINE_SLUG="${UNMATCHED_COPYNUMBER_PIPELINE_SLUG:-${_DEFAULT_PIPELINE_SLUG}}"
# Declared before any trap can fire, so `set -u` cannot turn a trap into a second,
# different failure. The empty ones are set in FILE SYSTEM SETUP. _LOCK_CONFLICT_RC
# sits outside flock's own sysexits range (64-71), so "another run holds it" cannot
# be confused with "flock itself failed".
LOCK_FILE=""
_LOCK_FD=""
_HOLDS_PIPELINE_LOCK=0
_LOCK_CONFLICT_RC=75
_WORK_DIR_DISPOSITION="kept"
_LAUNCHER_EXTRA_NOTE=""
STATS_DIR=""
STATS_FILE=""
_PIPELINE_START_EPOCH=""
_PIPELINE_WALL_SECONDS=""
_WORK_DIR_BYTES=""
_WORK_DIR_INODES=""

###########################
#### ENVIRONMENT SETUP ####
###########################

# "none" skips sourcing and relies on MANUAL ENVIRONMENT OVERRIDES below.
SOURCE_ME=${SOURCE_ME:-"${_DEFAULT_SOURCE_ME}"}
if [[ "${SOURCE_ME}" != "none" ]]; then
  check_for_source_me
  source "${SOURCE_ME}"
fi

############################
#### OPT-IN REPORTING   ####
############################

# Edit these to change every run. Opted-out runs make no network calls.
_DEFAULT_WEBSITE_LOGGING="true"        # log this run to the Dermatlas website
_DEFAULT_SLACK_NOTIFICATIONS="true"    # send a Slack message on completion
_DEFAULT_CLEANUP_WORK_DIR="true"       # delete this run's work dir after a successful run

# Precedence: submitting shell (_ENV_*) > source_me.sh > the default above, so a
# one-off run can opt out without editing this file or a shared source_me.sh:
#   export DERMATLAS_CLEANUP_WORK_DIR=false
#   bsub ... < commands/<pipeline>/run_unmatched_copynumber.sh
DERMATLAS_WEBSITE_LOGGING="${_ENV_WEBSITE_LOGGING:-${DERMATLAS_WEBSITE_LOGGING:-${_DEFAULT_WEBSITE_LOGGING}}}"
DERMATLAS_SLACK_NOTIFICATIONS="${_ENV_SLACK_NOTIFICATIONS:-${DERMATLAS_SLACK_NOTIFICATIONS:-${_DEFAULT_SLACK_NOTIFICATIONS}}}"
DERMATLAS_CLEANUP_WORK_DIR="${_ENV_CLEANUP_WORK_DIR:-${DERMATLAS_CLEANUP_WORK_DIR:-${_DEFAULT_CLEANUP_WORK_DIR}}}"

# Validated here, not where they are used: a typo must fail the launch now, not
# hours later, and onComplete must see exactly "true" or "false".
normalize_bool DERMATLAS_WEBSITE_LOGGING
normalize_bool DERMATLAS_SLACK_NOTIFICATIONS
normalize_bool DERMATLAS_CLEANUP_WORK_DIR
export DERMATLAS_WEBSITE_LOGGING DERMATLAS_SLACK_NOTIFICATIONS DERMATLAS_CLEANUP_WORK_DIR

# The module is only loaded if dermatlas-http is not already on PATH.
DERMATLAS_HTTP_MODULE="${DERMATLAS_HTTP_MODULE:-${_DEFAULT_DERMATLAS_HTTP_MODULE}}"

printf 'Dermatlas website logging:      %s\n' "${DERMATLAS_WEBSITE_LOGGING}"
printf 'Dermatlas Slack notifications:  %s\n' "${DERMATLAS_SLACK_NOTIFICATIONS}"
printf 'Work directory cleanup:         %s\n' "${DERMATLAS_CLEANUP_WORK_DIR}"

# Only now: source_me.sh has had its chance to export SLACK_WEBHOOK_URL.
_TRAP_CAN_SLACK=1

######################################
#### MANUAL ENVIRONMENT OVERRIDES ####
######################################

# The full environment contract, one commented export per variable. Uncomment
# and fill in to run without a source_me.sh (SOURCE_ME="none") or to override
# individual values - anything exported here wins over source_me.sh. To write
# your own source_me.sh instead, copy this block into a file and point
# SOURCE_ME at it.
#
# Pipeline-essential (always required):
# export PROJECT_DIR=""    # e.g. "/lustre/scratch127/.../<my_project>"; the pipeline's work/log/trace dirs are created under it
# export COMMANDS_DIR=""   # e.g. "${PROJECT_DIR}/commands"; CONFIG below defaults to living beneath it
# export ANALYSIS_DIR=""   # e.g. "${PROJECT_DIR}/analysis"; results land in ${ANALYSIS_DIR}/unmatched_copy_number (config: outdir)
# export BAMS_DIR=""       # e.g. "${PROJECT_DIR}/bams"; the config globs ${BAMS_DIR}/**bam{,.bai} (config: bam_path)
# export STUDY=""          # e.g. "6937"; part of the output file prefix (config: study_id, cohort_prefix) and the run id
# export PROJECT=""        # e.g. "3125"; part of the output file prefix and the run id
# export COHORT=""         # e.g. "PILAR_CYST"; completes the output file prefix (config: cohort_prefix)
# export DNA_TUMOUR_LIST_ANALYSED_ALL=""                # e.g. "${PROJECT_DIR}/metadata/6937_3125-analysed_all_tum.txt"; the sample universe (config: all_samples) and a subcohort
# export DNA_TUMOUR_LIST_ONE_TUMOUR_PER_PATIENT_ALL=""  # e.g. "${PROJECT_DIR}/metadata/6937_3125-one_tumour_per_patient_all_tum.txt"; a subcohort
# export DNA_TUMOUR_LIST_RELATED_TUMOURS_ALL=""         # e.g. "${PROJECT_DIR}/metadata/6937_3125-related_tumours_all_tum.txt"; a subcohort, may be empty
# export COHORT_METADATA_FILE=""   # e.g. "${PROJECT_DIR}/metadata/<cohort>_metadata.tsv"; sample id, sex and tumour/normal columns (config: metadata_manifest). Must be TSV, not .xlsx
#
# Website-essential (required only when DERMATLAS_WEBSITE_LOGGING="true"):
# export COHORT_SLUG=""              # e.g. "m10-cutaneous-mixed-tumour"; keys the analysis-log record
# export SAMPLE_LIST_VERSION_FILE="" # e.g. "${PROJECT_DIR}/metadata/VERSION"; a file holding one integer >= 1
# export SELF_DESCRIBING_API=""      # e.g. "https://team113.sanger.ac.uk/api/v1/resolve/"; dermatlas-http --api endpoint
# export DERMATLAS_HTTP_MODULE=""    # e.g. "dermatlas-http"; the module loaded to put dermatlas-http on PATH for the work-dir report
#
# Slack-essential (required only when DERMATLAS_SLACK_NOTIFICATIONS="true"):
# export SLACK_WEBHOOK_URL=""        # e.g. "https://hooks.slack.com/services/T000/B000/XXXX"
#
# The DERMATLAS_* toggles do not belong here: OPT-IN REPORTING above has already
# resolved them. Set them in your shell or source_me.sh, or edit their defaults.

################################
#### ENVIRONMENT VALIDATION ####
################################

require_env 'Regenerate source_me.sh with a dermanager version that exports these variables (or export them yourself if you git-cloned the pipeline).' \
            "${_PIPELINE_ENV_VARS[@]}"
if [[ "${DERMATLAS_WEBSITE_LOGGING}" == "true" ]]; then
  require_env 'Needed because DERMATLAS_WEBSITE_LOGGING=true. Fix source_me.sh, or set DERMATLAS_WEBSITE_LOGGING="false" in this script to opt out of website logging.' \
              "${_WEBSITE_ENV_VARS[@]}"
fi
if [[ "${DERMATLAS_SLACK_NOTIFICATIONS}" == "true" ]]; then
  require_env 'Needed because DERMATLAS_SLACK_NOTIFICATIONS=true. Fix source_me.sh, or set DERMATLAS_SLACK_NOTIFICATIONS="false" in this script to opt out of Slack notifications.' \
              "${_SLACK_ENV_VARS[@]}"
fi

###########################
#### RUN CONFIGURATION ####
###########################

# Nextflow config for this run; git-clone runs point this at their own copy.
CONFIG="${COMMANDS_DIR}/${PIPELINE_SLUG}/unmatched_copynumber.config"
# Pipeline version to run: a tag or commit hash.
REVISION="0.2.2"
# Optional. If set, RUN_ID becomes <label>_<timestamp>.
LABEL=""

###########################
#### FILE SYSTEM SETUP ####
###########################

PIPELINE_DIR="${PROJECT_DIR}/${PIPELINE_SLUG}"
mkdir -p "${PIPELINE_DIR}"

# Above the lock, so the lock file and both sentinels can always name this run.
if [[ -n "${LABEL:-}" ]]; then
  RUN_ID="$(set_run_id_from_label "${LABEL}")"
else
  RUN_ID="$(set_run_id "${STUDY:-}" "${PROJECT:-}" "${COHORT_SLUG:-}")"
fi
# Exported so the pipeline's nextflow.config and onComplete handler read it back.
export RUN_ID

# One launcher owns ${PIPELINE_DIR} at a time, taken before anything under it is
# written. Fails fast rather than queueing: two concurrent submissions of one cohort
# would share a work/ directory, and the first to finish would delete it under the
# second.
LOCK_FILE="${PIPELINE_DIR}/.lock"
acquire_pipeline_lock

export NXF_WORK="${PIPELINE_DIR}/work"
export NXF_TEMP="${PIPELINE_DIR}/tmp"
mkdir -p "${NXF_WORK}" "${NXF_TEMP}"

# One clone per revision. The shared default at ${NXF_HOME}/assets is checked out in
# place by every pull and run, so parallel runs on different revisions clobber each
# other, and a first pull of a newly published tag fails with "Cannot find revision".
export NXF_ASSETS="${PIPELINE_DIR}/clones/${REVISION}"
mkdir -p "${NXF_ASSETS}"

# Pinned because an LSF job inherits the submitter's profile, and because the
# nextflow default lives under work/, which the cleanup trap deletes. Matches
# singularity.cacheDir in the farm22 profile.
export NXF_SINGULARITY_CACHEDIR="/lustre/scratch127/casm/projects/dermatlas/singularity_images"

# Exported rather than passed as -log, so onComplete can read the path back with
# System.getenv('NXF_LOG_FILE') and name it in the Slack message.
LOG_DIR="${PIPELINE_DIR}/logs"
NXF_RUN_LOG_FILE="${LOG_DIR}/nextflow-run-${RUN_ID}.log"
NXF_PULL_LOG_FILE="${LOG_DIR}/nextflow-pull-${RUN_ID}.log"
# Consumed by the pipeline's nextflow.config.
export TRACE_DIR="${PIPELINE_DIR}/traces"
# Made here, not in the exit trap: a mkdir that fails inside a trap is a second,
# harder failure.
STATS_DIR="${PIPELINE_DIR}/stats"
STATS_FILE="${STATS_DIR}/resource-stats-${RUN_ID}.txt"
mkdir -p "${LOG_DIR}" "${TRACE_DIR}" "${STATS_DIR}"

###################################
#### EXECUTION OF THE PIPELINE ####
###################################

module load nextflow-23.10.0
module load /software/modules/ISG/singularity/3.11.4

# Nextflow conventions seems to encourage running from the pipeline directory.
cd "${PIPELINE_DIR}"

NXF_LOG_FILE="${NXF_PULL_LOG_FILE}" \
  nextflow pull "https://github.com/team113sanger/dermatlas_unmatched_copynumber_nf" -r "${REVISION}"

# From here failures are reported by workflow.onComplete (lib/Utils.groovy), so the
# launcher trap must not fire as well.
trap - EXIT ERR INT TERM HUP
trap 'on_pipeline_exit $?' EXIT
# Required: for an untrapped fatal signal bash runs the EXIT trap with $? from the last
# *completed* command, so a killed run would look successful and have its work directory
# deleted. 128+n makes it the failure it is.
trap 'exit 130' INT
trap 'exit 143' TERM
trap 'exit 129' HUP

_PIPELINE_START_EPOCH="$(date +%s)"
NXF_LOG_FILE="${NXF_RUN_LOG_FILE}" \
  nextflow run "https://github.com/team113sanger/dermatlas_unmatched_copynumber_nf" \
  -resume \
  -c "${CONFIG}" \
  -r "${REVISION}" \
  -profile farm22 \
  -work-dir "${NXF_WORK}"
