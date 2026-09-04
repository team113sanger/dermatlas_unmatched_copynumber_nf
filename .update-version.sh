#!/usr/bin/env bash
#
# .update-version.sh - set this pipeline's semantic version everywhere it is recorded.
#
# Part of the release procedure in README.md ("Cutting a release"). Sibling Dermatlas
# pipelines carry their own copy of this script; only the TARGETS table below differs.

set -euo pipefail

readonly SCRIPT_NAME="${0##*/}"
REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
readonly REPO_ROOT

# A bare X.Y.Z release version: the same shape as the git tag, and as the tag pattern
# .github/workflows/publish-assets.yml builds release bundles for.
readonly SEMVER_RE='[0-9]+\.[0-9]+\.[0-9]+'

# The files carrying the version, as "<repo-relative path>|<extended regex>".
#
# Each regex must match exactly one line and capture three groups: the text before the
# version, the version itself, and the text after. Rewriting through the capture groups
# keeps each file's own quoting and column alignment intact.
#
# '|' separates the two fields, so it cannot appear inside a pattern.
readonly -a TARGETS=(
  "assets/run_unmatched_copynumber.sh|^(REVISION=\")(${SEMVER_RE})(\")$"
  "docs/source/conf.py|^(release = ')(${SEMVER_RE})(')$"
  "nextflow.config|^([[:blank:]]*version[[:blank:]]+= ')(${SEMVER_RE})(')$"
)

# Set once main() starts rewriting; removed on any exit path.
TMP_FILE=""

cleanup() {
  [[ -n ${TMP_FILE} && -e ${TMP_FILE} ]] && rm -f "${TMP_FILE}"
  return 0
}
trap cleanup EXIT

die() {
  printf '%s: error: %s\n' "${SCRIPT_NAME}" "$*" >&2
  exit 1
}

usage() {
  cat <<EOF
Usage: ${SCRIPT_NAME} <version>
       ${SCRIPT_NAME} --help

Set this pipeline's semantic version in every file that records it.

Arguments:
  <version>     Release version as X.Y.Z, e.g. 0.4.5. No 'v' prefix and no
                pre-release or build suffix: this is the string that becomes
                the git tag.

Options:
  -h, --help    Show this message and exit.

Files updated:
$(printf '  - %s\n' "${TARGETS[@]%%|*}")

Re-running with the same version is a no-op, so the script is safe to repeat.
Nothing is written unless every file above exists and holds exactly one version
line, so a refactor that moves or renames a version cannot be silently skipped.

Run it on the release branch and commit the result; see "Cutting a release" in
README.md for the surrounding steps.
EOF
}

main() {
  local -a args=()

  while (($#)); do
    case "$1" in
      -h | --help)
        usage
        return 0
        ;;
      --)
        shift
        args+=("$@")
        break
        ;;
      -*)
        die "unknown option: $1 (try --help)"
        ;;
      *)
        args+=("$1")
        ;;
    esac
    shift
  done

  if ((${#args[@]} != 1)); then
    usage >&2
    die "expected exactly one <version> argument, got ${#args[@]}"
  fi

  local version="${args[0]}"
  if [[ ! ${version} =~ ^${SEMVER_RE}$ ]]; then
    die "not an X.Y.Z version: '${version}'"
  fi

  cd "${REPO_ROOT}"

  # Pre-flight every target before writing to any of them, so a pattern that has gone
  # stale cannot leave the tree half-updated.
  local target path pattern count
  for target in "${TARGETS[@]}"; do
    path="${target%%|*}"
    pattern="${target#*|}"
    [[ -f ${path} ]] || die "no such file: ${path}"
    count="$(grep -cE "${pattern}" "${path}" || true)"
    if ((count != 1)); then
      die "${path}: expected 1 version line, found ${count} - has the file changed shape?"
    fi
  done

  printf 'Setting version to %s\n' "${version}"

  local old
  TMP_FILE="$(mktemp)"
  for target in "${TARGETS[@]}"; do
    path="${target%%|*}"
    pattern="${target#*|}"
    old="$(sed -nE "s|${pattern}|\2|p" "${path}")"

    # Rewrite through a temp file rather than `sed -i`: portable across GNU and BSD
    # sed, and copying back preserves the target's original mode and ownership.
    sed -E "s|${pattern}|\1${version}\3|" "${path}" >"${TMP_FILE}"
    cat "${TMP_FILE}" >"${path}"

    if [[ ${old} == "${version}" ]]; then
      printf '  %-26s %s (unchanged)\n' "${path}" "${version}"
    else
      printf '  %-26s %s -> %s\n' "${path}" "${old}" "${version}"
    fi
  done
}

main "$@"
