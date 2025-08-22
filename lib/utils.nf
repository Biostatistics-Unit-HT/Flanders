// Check and parse numeric intervals
def parseIntervals(intervalString) {
    if (!intervalString?.trim()) {
        error "Interval string cannot be null or empty"
    }
    
    def clean = intervalString.replaceAll(/\s+/, '')
    
    // Validate format
    if (!(clean ==~ /^[\d,-]+$/ || clean ==~ /.*[--,,].*/ || clean ==~ /^[-,].*/ || clean ==~ /.*[-,]$/ || clean ==~ /.*[-,][-,].*/ )) {
        error "Invalid interval format: '${intervalString}'"
    }
    
    def result = []
    
    try {
        clean.split(',').each { part ->
            if (part.contains('-')) {
                def (start, end) = part.split('-').collect { it.toInteger() }
                if (start > end) error "Invalid range: ${start}-${end}"
                (start..end).each { if (!result.contains(it)) result << it }
            } else {
                def num = part.toInteger()
                if (!result.contains(num)) result << num
            }
        }
        return result.sort()
    } catch (NumberFormatException e) {
        error "Invalid number in: '${intervalString}'"
    }
}