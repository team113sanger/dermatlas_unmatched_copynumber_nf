/**
 * Small collection of static helpers used across the pipeline.
 *
 * Classes placed in lib/ are added to the Nextflow classpath automatically,
 * so these methods can be called from main.nf / nextflow.config without an
 * explicit import (e.g. `Utils.reportRun(workflow, params)`).
 */
class Utils {

    /** Minimum dermatlas-http version whose `cohort analysis-log` supports the flags we pass. */
    private static final String MIN_DERMATLAS_HTTP_VERSION = '0.6.1'

    /**
     * Parse a Nextflow CLI param into a strict boolean.
     *
     * Accepts:
     *  - Boolean: true/false
     *  - Number: 1/0
     *  - String (case-insensitive, trimmed): "true"/"false"/"1"/"0"
     *
     * If value is null (or empty string), returns defaultValue.
     * Otherwise throws IllegalArgumentException on invalid inputs.
     */
    static boolean parseCLIBool(Object givenValue, boolean defaultValue = false, String paramName = null) {
        if (givenValue == null) {
            return defaultValue
        }

        if (givenValue instanceof Boolean) {
            return (Boolean) givenValue
        }

        if (givenValue instanceof Number) {
            int n = ((Number) givenValue).intValue()
            if (n == 0) return false
            if (n == 1) return true
            throw new IllegalArgumentException(errMsg(paramName, givenValue, "expected 0 or 1"))
        }

        String s = givenValue.toString().trim()
        if (s.length() == 0) {
            throw new IllegalArgumentException(errMsg(paramName, givenValue, "empty string is not a valid, expected true/false/1/0"))
        }

        s = s.toLowerCase()
        if (s == "true" || s == "1")  return true
        if (s == "false" || s == "0") return false

        throw new IllegalArgumentException(errMsg(paramName, givenValue, "expected true/false/1/0 (case-insensitive)"))
    }

    // ---- Run-completion entry point ----

    /**
     * Single entry point for `workflow.onComplete`. Drop this into any pipeline:
     *
     *   workflow.onComplete { Utils.reportRun(workflow, params) }
     *
     * Notifications and the analysis-log call are best-effort: a failure here
     * never changes the pipeline's exit status.
     *
     * Both reporting channels are explicit opt-ins, gated on environment
     * variables exported by the pipeline's launcher (assets/run_<pipeline>.sh):
     *   DERMATLAS_SLACK_NOTIFICATIONS=true  - send a Slack message
     *   DERMATLAS_WEBSITE_LOGGING=true      - record the run in the Dermatlas
     *                                         analysis log via dermatlas-http
     * Unset means opted out, so bare `nextflow run` / test runs stay silent.
     * Stub runs (`-stub-run`, or params.is_stub) never report anywhere.
     *
     * To keep reporting decoupled from the workflow, every payload value that
     * the workflow itself does not need is read from the environment rather
     * than params: COHORT_SLUG, SAMPLE_LIST_VERSION_FILE, SELF_DESCRIBING_API,
     * SLACK_WEBHOOK_URL, RUN_ID, NXF_LOG_FILE. From `params` only
     * `is_stub`, `analysis_pipeline_slug`, `trace_file` (and an optional
     * `run_id` override) are consulted; the rest comes from `workflow.manifest`.
     */
    static void reportRun(def workflow, def params) {
        println "Pipeline completed at: ${workflow.complete}"
        if (workflow.success) {
            println "Pipeline finished successfully!"
        } else {
            System.err.println "Pipeline finished with errors!"
        }

        boolean isStub       = isStubRun(workflow, params)
        boolean slackOptIn   = envOptIn('DERMATLAS_SLACK_NOTIFICATIONS')
        boolean websiteOptIn = envOptIn('DERMATLAS_WEBSITE_LOGGING')
        String cohort   = (System.getenv('COHORT_SLUG') ?: '').trim()
        String runId    = buildRunId(params)
        // Label shown in Slack: prefer the human-readable manifest description,
        // fall back to the manifest name, then a generic default.
        String label    = (workflow.manifest?.description ?: workflow.manifest?.name ?: 'pipeline').toString()
        String version  = workflow.manifest?.version
        // Read once: both the Slack message and the analysis-log record report it.
        Integer sampleVer = readSampleListVersion(System.getenv('SAMPLE_LIST_VERSION_FILE'))

        // Slack notification (explicit opt-in; a webhook is still needed to actually send)
        String slackUrl = (System.getenv('SLACK_WEBHOOK_URL') ?: '').trim()
        if (slackOptIn && !isStub) {
            if (!slackUrl) {
                System.err.println "WARNING: DERMATLAS_SLACK_NOTIFICATIONS=true but SLACK_WEBHOOK_URL is unset; no Slack message sent"
            } else {
                boolean sent
                if (workflow.success) {
                    sent = sendSlackSuccess(slackUrl, runId, version, workflow.duration, label, cohort, sampleVer)
                } else {
                    String nxfLog = resolveNextflowLog(workflow)
                    sent = sendSlackFailure(slackUrl, runId, version, workflow.duration, label,
                                            workflow.errorReport, nxfLog, cohort, sampleVer)
                }
                if (sent) {
                    println "Slack notification sent"
                } else {
                    System.err.println "WARNING: Slack notification failed (pipeline result is unaffected)"
                }
            }
        }

        // Append a run record to the cohort analysis-log (explicit opt-in, via dermatlas-http)
        if (websiteOptIn && !isStub) {
            String apiUrl       = (System.getenv('SELF_DESCRIBING_API') ?: '').trim()
            String runStatus    = workflow.success ? 'completed' : 'failed'
            String pipelineSlug = (params.analysis_pipeline_slug ?: '').toString().trim()
            Long   durationSecs = workflow.duration?.toSeconds()
            String nxfLog       = resolveNextflowLog(workflow)
            String traceFile    = (params.trace_file ?: '').toString()
            if (reportCohortAnalysisLog(apiUrl, cohort, pipelineSlug, runStatus, sampleVer,
                                        durationSecs, version, nxfLog, traceFile)) {
                println "Reported analysis run for cohort '${cohort}' (${runStatus}) to the Dermatlas analysis log"
            } else {
                System.err.println "WARNING: Failed to report analysis run for cohort '${cohort}' to the Dermatlas analysis log (pipeline result is unaffected)"
            }
        }
    }

    // ---- Slack notification helpers ----

    /**
     * Notify Slack that the pipeline completed successfully.
     *
     * @param label human-readable pipeline name (e.g. workflow.manifest.description)
     * @param sampleListVersion sample-list version to show in the message; omitted when null
     */
    static boolean sendSlackSuccess(String webhookUrl, String runId, String version,
                                    def duration, String label, String cohortSlug = '',
                                    Integer sampleListVersion = null) {
        def ref = runRef(runId, cohortSlug)
        String msg = ":white_check_mark: *${label}* completed successfully - `${ref}`\n" +
            metaLine(version, sampleListVersion, duration)

        return postSlackMessage(webhookUrl, msg)
    }

    /**
     * Notify Slack that the pipeline failed.
     *
     * Parses workflow.errorReport to extract the failed process name and work
     * directory (if present). Includes NXF_LOG_FILE path when available.
     *
     * @param errorReport  workflow.errorReport from Nextflow (may be null)
     * @param nxfLogFile   Nextflow log path (may be null or empty)
     * @param sampleListVersion sample-list version to show in the message; omitted when null
     */
    static boolean sendSlackFailure(String webhookUrl, String runId, String version,
                                    def duration, String label, String errorReport,
                                    String nxfLogFile, String cohortSlug = '',
                                    Integer sampleListVersion = null) {
        def ref = runRef(runId, cohortSlug)
        def msg = ":octagonal_sign: *${label}* failed - `${ref}`\n" +
                  metaLine(version, sampleListVersion, duration)

        def processName = parseProcessName(errorReport)
        def workDir = parseWorkDir(errorReport)

        if (processName) {
            msg += "\nProcess: `${processName}`"
        }
        if (workDir) {
            msg += "\nWork dir: `${workDir}`"
        }
        if (nxfLogFile) {
            msg += "\nNextflow Log: `${nxfLogFile}`"
        }

        return postSlackMessage(webhookUrl, msg)
    }

    // ---- Cohort analysis-log helpers ----

    /**
     * Append one analysis-pipeline run record for a single cohort by running
     * the companion CLI:
     *
     *   dermatlas-http cohort analysis-log --api ... --cohort ... --pipeline ...
     *     --run-status ... --sample-list-version N --execution-duration SECONDS
     *     [--pipeline-version TAG] [--nextflow-log FILE] [--nextflow-performance FILE]
     *
     * The CLI owns the API schema (and its versioning drift), so this method
     * only assembles arguments. `runStatus` is one of the CLI's run-status
     * values (not_run / queued / running / completed / failed / blocked /
     * unknown); `sampleListVersion` must be an integer >= 1; `pipelineVersion`
     * is optional. The log and performance (execution trace) paths are only
     * passed when the files exist - both flags are optional to the CLI.
     *
     * Required-field problems are caught here (logged, returns false) and any
     * subprocess failure is swallowed by runAnalysisLogCli, so logging a run
     * can never change the pipeline's exit status.
     */
    static boolean reportCohortAnalysisLog(String apiUrl, String cohortSlug, String pipelineSlug,
                                           String runStatus, Integer sampleListVersion,
                                           Long executionSeconds, String pipelineVersion,
                                           String nextflowLog, String traceFile,
                                           int timeoutMs = 60_000) {
        if (!apiUrl) {
            System.err.println("WARNING: analysis-log skipped — SELF_DESCRIBING_API is not set")
            return false
        }
        if (!cohortSlug || !pipelineSlug) {
            System.err.println("WARNING: analysis-log skipped — cohort_slug and analysis_pipeline_slug are both required")
            return false
        }
        if (sampleListVersion == null || sampleListVersion < 1) {
            System.err.println("WARNING: analysis-log skipped for cohort '${cohortSlug}' — sample_list_version must be an integer >= 1 (got ${sampleListVersion})")
            return false
        }
        long durationSecs = (executionSeconds != null && executionSeconds >= 0L) ? executionSeconds : 0L
        List<String> args = [
            '--api', apiUrl,
            '--cohort', cohortSlug,
            '--pipeline', pipelineSlug,
            '--run-status', runStatus,
            '--sample-list-version', sampleListVersion.toString(),
            '--execution-duration', durationSecs.toString(),
        ]
        if (pipelineVersion) {
            args += ['--pipeline-version', pipelineVersion]
        }
        [['--nextflow-log', nextflowLog], ['--nextflow-performance', traceFile]].each { flag, path ->
            if (path) {
                def f = new File(path.toString())
                if (f.isFile()) {
                    args += [flag.toString(), f.absolutePath]
                }
            }
        }
        return runAnalysisLogCli(args, timeoutMs)
    }

    /**
     * Run `dermatlas-http cohort analysis-log <cliArgs>` in the launch environment.
     *
     * One login-shell (`bash -lc`, so `module` - a shell function set up by
     * profile scripts - is available) invocation covers everything, so the
     * PATH set by `module load` applies to the version probe and the real
     * call alike: prefer dermatlas-http already on PATH, else fall back to
     * `module load dermatlas-http`, then require the CLI to be at least
     * MIN_DERMATLAS_HTTP_VERSION before exec'ing it.
     *
     * The args ride as positional shell parameters and are never interpolated
     * into the script text, so there are no quoting/injection concerns. Exit
     * 127 is the "tool unavailable or too old" sentinel - the shell has
     * already written the specific warning to stderr in that case. Never
     * throws; returns true only when the CLI exits 0.
     */
    private static boolean runAnalysisLogCli(List<String> cliArgs, int timeoutMs) {
        // $1 = minimum version (then shifted away); $2.. = dermatlas-http args.
        String script = '''\
set -u
if ! command -v dermatlas-http >/dev/null 2>&1; then
  if command -v module >/dev/null 2>&1; then
    module load dermatlas-http >/dev/null 2>&1 || { echo "WARNING: 'module load dermatlas-http' failed; skipping analysis-log" >&2; exit 127; }
  fi
fi
command -v dermatlas-http >/dev/null 2>&1 || { echo "WARNING: dermatlas-http not on PATH and no module system found; skipping analysis-log" >&2; exit 127; }
min_version="$1"; shift
version_line="$(dermatlas-http cohort analysis-log --version 2>/dev/null || true)"
version="${version_line##* }"
case "$version" in
  [0-9]*.[0-9]*.[0-9]*) ;;
  *) echo "WARNING: could not parse dermatlas-http version from '${version_line}'; skipping analysis-log" >&2; exit 127;;
esac
printf '%s\\n%s\\n' "$min_version" "$version" | sort -V -C || { echo "WARNING: dermatlas-http ${version} is older than required ${min_version}; skipping analysis-log" >&2; exit 127; }
exec dermatlas-http cohort analysis-log "$@"
'''
        try {
            def pb = new ProcessBuilder(
                ['bash', '-lc', script, 'dermatlas-http-analysis-log', MIN_DERMATLAS_HTTP_VERSION] + cliArgs
            )
            Process p = pb.start()
            def out = new StringBuffer()
            def err = new StringBuffer()
            p.consumeProcessOutput(out, err)   // drain both streams to avoid pipe-buffer deadlock
            boolean finished = p.waitFor(timeoutMs, java.util.concurrent.TimeUnit.MILLISECONDS)
            if (!finished) {
                p.destroyForcibly()
                System.err.println("WARNING: dermatlas-http analysis-log timed out after ${timeoutMs}ms")
                return false
            }
            int code = p.exitValue()
            if (code == 0) {
                return true
            }
            if (code != 127) {
                System.err.println("WARNING: dermatlas-http analysis-log exited ${code}: ${err.toString().trim()}")
            } else if (err.length() > 0) {
                System.err.println(err.toString().trim())
            }
            return false
        } catch (Exception e) {
            System.err.println("WARNING: dermatlas-http analysis-log failed: ${e.message}")
            return false
        }
    }

    // ---- private helpers ----

    /**
     * True when this is a stub run: `-stub-run` on the CLI (workflow.stubRun,
     * Nextflow >= 22.10) OR the legacy params.is_stub flag. workflow.stubRun
     * is read defensively so older engines (the manifest allows >= 22.04.5)
     * fall back to params.is_stub alone.
     */
    private static boolean isStubRun(def workflow, def params) {
        boolean fromParams = parseCLIBool(params.is_stub, false)
        boolean fromMeta = false
        try {
            fromMeta = workflow?.stubRun as boolean
        } catch (Exception ignored) {
            // workflow.stubRun not available on this Nextflow version
        }
        return fromParams || fromMeta
    }

    /**
     * Strict opt-in read of an environment toggle. Unset means false:
     * reporting is explicit opt-in, and the run wrapper always exports these,
     * so only direct `nextflow run` invocations (tests, CI) fall through to
     * silence. Anything other than "true"/"false" (case-insensitive, trimmed)
     * is warned about and treated as false - onComplete is too late to abort,
     * and the wrapper validates the value before launch.
     */
    private static boolean envOptIn(String name) {
        String raw = System.getenv(name)
        if (raw == null) {
            return false
        }
        String s = raw.trim().toLowerCase()
        if (s == 'true')  return true
        if (s == 'false') return false
        System.err.println("WARNING: ${name}='${raw}' is not 'true' or 'false'; treating as false")
        return false
    }

    private static String errMsg(String paramName, Object givenValue, String expectation) {
        String namePart = paramName ? "--${paramName} " : ""
        return "Invalid ${namePart}${givenValue} (${expectation})"
    }

    /**
     * The id a run is referred to by in Slack messages.
     *
     * An explicit `params.run_id` wins outright. Otherwise prefer the RUN_ID
     * exported by the run wrapper - the same id that names the nextflow log,
     * execution trace and execution report, so the Slack reference matches the
     * files on disk. The last resort is assembled from the STUDY / PROJECT
     * env vars.
     */
    private static String buildRunId(def params) {
        if (params?.run_id) {
            return params.run_id.toString()
        }
        String fromEnv = (System.getenv('RUN_ID') ?: '').trim()
        if (fromEnv) {
            return fromEnv
        }
        String study    = System.getenv('STUDY') ?: ''
        String project  = System.getenv('PROJECT') ?: ''
        return "${study}_${project}"
    }

    /**
     * Read the sample-list version from a file whose contents are a single
     * integer. Returns null (with a stderr warning) if the path is unset, the
     * file is missing, or it does not contain an integer - in which case it is
     * omitted from the Slack message and the analysis-log call is skipped.
     */
    private static Integer readSampleListVersion(def sampleListVersionPath) {
        if (!sampleListVersionPath) {
            return null
        }
        def f = new File(sampleListVersionPath.toString())
        if (!f.exists()) {
            System.err.println("WARNING: sample_list_version file '${f}' not found; the version will be omitted from reporting")
            return null
        }
        try {
            return f.text.trim() as Integer
        } catch (NumberFormatException e) {
            System.err.println("WARNING: sample_list_version file '${f}' does not contain an integer; the version will be omitted from reporting")
            return null
        }
    }

    /**
     * The italic metadata line shared by the success and failure Slack messages, e.g.
     *   _Pipeline v0.4.2 | Sample list: v1 | Duration: 1h 46m 23s_
     *
     * The sample-list segment is dropped entirely when the version could not be read, so the
     * line degrades to its previous form rather than showing "vnull".
     */
    private static String metaLine(String version, Integer sampleListVersion, def duration) {
        String samplePart = sampleListVersion != null ? " | Sample list: v${sampleListVersion}" : ""
        return "_Pipeline v${version}${samplePart} | Duration: ${duration}_"
    }

    /**
     * Absolute path of this run's Nextflow log, or '' if it cannot be found.
     *
     * Our run wrappers export NXF_LOG_FILE (which also tells Nextflow where to write it).
     * Falling back to the launch directory keeps this working for anyone invoking
     * `nextflow run` directly, where that variable is unset.
     */
    private static String resolveNextflowLog(def workflow) {
        String fromEnv = System.getenv('NXF_LOG_FILE')
        if (fromEnv) {
            return fromEnv
        }
        def fallback = new File("${workflow.launchDir}/.nextflow.log")
        return fallback.isFile() ? fallback.absolutePath : ''
    }

    /**
     * Human-facing reference for a pipeline run in Slack messages.
     *
     * The cohort slug is optional (the pipeline can run without it). When it is
     * set, the run is shown as "<cohort_slug> (<run_id>)" so readers see the
     * cohort first; otherwise the bare run id is used.
     *
     * @param runId      always present (e.g. study_project_genotype)
     * @param cohortSlug optional cohort slug (null/empty when not configured)
     */
    private static String runRef(String runId, String cohortSlug) {
        return (cohortSlug != null && !cohortSlug.trim().isEmpty())
            ? "${cohortSlug.trim()} (${runId})"
            : runId
    }

    /**
     * Extract the process name from workflow.errorReport.
     * Matches: Error executing process > 'MY_PROCESS (tag)'
     * @return process name string or null
     */
    private static String parseProcessName(String errorReport) {
        if (!errorReport) return null
        def m = errorReport =~ /Error executing process > '([^']+)'/
        return m.find() ? m.group(1) : null
    }

    /**
     * Extract the work directory from workflow.errorReport.
     * Matches the line after "Work dir:" in the report.
     * @return work directory path or null
     */
    private static String parseWorkDir(String errorReport) {
        if (!errorReport) return null
        def m = errorReport =~ /(?m)Work dir:\s*\n\s*(\S+)/
        return m.find() ? m.group(1) : null
    }

    /**
     * POST a message to a Slack incoming webhook.
     *
     * Network errors and non-2xx responses are caught and logged to stderr
     * so that a notification failure never causes the pipeline itself to fail.
     *
     * @return true if the webhook accepted the message (HTTP 2xx), false otherwise
     */
    private static boolean postSlackMessage(String webhookUrl, String message, int timeoutMs = 10_000) {
        if (!webhookUrl || !message) {
            return false
        }

        try {
            def url = new java.net.URL(webhookUrl)
            def conn = (java.net.HttpURLConnection) url.openConnection()
            conn.setRequestMethod('POST')
            conn.setRequestProperty('Content-Type', 'application/json; charset=UTF-8')
            conn.setConnectTimeout(timeoutMs)
            conn.setReadTimeout(timeoutMs)
            conn.setDoOutput(true)

            // Minimal JSON-safe escaping for the message text
            def escaped = message
                .replace('\\', '\\\\')
                .replace('"', '\\"')
                .replace('\n', '\\n')

            conn.outputStream.withWriter('UTF-8') { writer ->
                writer.write("{\"text\":\"${escaped}\"}")
            }

            int code = conn.getResponseCode()
            if (code >= 200 && code < 300) {
                return true
            }
            System.err.println("WARNING: Slack webhook returned HTTP ${code}")
            return false
        } catch (SocketTimeoutException e) {
            System.err.println("WARNING: Slack notification timed out after ${timeoutMs}ms: ${e.message}")
            return false
        } catch (Exception e) {
            System.err.println("WARNING: Slack notification failed: ${e.message}")
            return false
        }
    }

}
