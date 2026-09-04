// Shared helpers for validating pipeline inputs before any work is scheduled, and for
// turning the cohort's threshold strings into the numbers the processes need.
//
// These exist so that an empty sample list or a malformed threshold string fails
// immediately with an actionable message, rather than after a day of CNVkit calling.

// Count the usable data rows in a delimited text file, ignoring blank lines and the
// first `skipLines` non-blank lines (i.e. the header).
def countDataRows(Object input, int skipLines) {
    def target = input instanceof CharSequence ? file(input) : input
    def lines = target.readLines().findAll { line -> line.trim() }
    return Math.max(lines.size() - skipLines, 0)
}

// Return the column names on the first non-blank line of a delimited file, or an empty
// list if the file has no content at all.
def readHeader(Object input, String sep) {
    def target = input instanceof CharSequence ? file(input) : input
    def first = target.readLines().find { line -> line.trim() }
    return first == null ? [] : first.split(sep, -1).collect { column -> column.trim() }
}

// Split a CNVkit `call -t` threshold string into its four log2 ratios.
//
// The string is the log2 ratio boundary for CN=0, CN=1, CN=3 and CN=4 in that order, so
// the two inner values are the loss (CN=1) and gain (CN=3) thresholds that every
// downstream step - penetrance plots and GISTIC2 alike - is keyed on.
def parseThresholds(String thresholds, String paramName) {
    def parts = thresholds.trim().split(',').collect { part -> part.trim() }
    if (parts.size() != 4) {
        error("${paramName} must be four comma-separated log2 ratios for CN=0,1,3,4 " +
              "(e.g. '-0.737,-0.322,0.263,0.485'), got '${thresholds}'.")
    }
    def numbers = parts.collect { part ->
        if (!(part ==~ /^-?\d+(\.\d+)?$/)) {
            error("${paramName} contains a value that is not a number: '${part}' (in '${thresholds}').")
        }
        return part as BigDecimal
    }
    def loss = numbers[1]
    def gain = numbers[2]
    if (loss >= 0 || gain <= 0) {
        error("${paramName} must have a negative CN=1 threshold and a positive CN=3 threshold, " +
              "got loss=${loss} gain=${gain} (in '${thresholds}').")
    }
    return [
        thresholds: parts.join(','),
        underscored: parts.join('_'),
        loss: loss,
        gain: gain,
        // GISTIC2 takes both amplitude thresholds as positive numbers.
        loss_amplitude: loss.abs()
    ]
}
