# This script contains bash functions re-used across processes


#-------------------------------------------------------------------------------
# check_pipeline
#
# Validate the exit status of the most recently executed pipeline.
#
# This function inspects Bash's PIPESTATUS array and returns the first
# non-zero exit code encountered. It is intended for use immediately after
# a pipeline when running with `set -o pipefail`, allowing Nextflow to
# receive the correct exit status and apply retry logic for transient
# failures (e.g. OOM kills).
#
# The function also prints the full PIPESTATUS array to stderr to aid
# debugging of multi-stage pipelines.
#
# Returns:
#   0   All pipeline stages completed successfully.
#   N   Exit code of the first failed stage.
#
# Common exit codes:
#   137  Process killed with SIGKILL (often OOM)
#   143  Process terminated with SIGTERM (walltime/scheduler cancellation)
#   139  Segmentation fault
#
# Notes:
#   - Must be called immediately after the pipeline.
#   - Any intervening command will overwrite PIPESTATUS.
#
# Example:
#   set +e
#   bcftools mpileup ... \
#       | bcftools call ... \
#       | bcftools view ...
#   st=("${PIPESTATUS[@]}")
#   set -e
#   check_pipeline "${st[@]}" || exit $?
#
#-------------------------------------------------------------------------------
check_pipeline() {
    local -a st=("$@")

    echo "PIPESTATUS: ${st[*]}" >&2

    for i in "${!st[@]}"; do
        if (( st[i] != 0 )); then
            echo "Pipeline stage $((i + 1)) failed with exit code ${st[i]}" >&2
            return "${st[i]}"
        fi
    done

    return 0
}