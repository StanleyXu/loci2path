setMethod(
    "show",
    signature = "eqtlSet",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat(" eQTL collected from tissue:", object@tissue, "\n")
        cat(" number of eQTLs:", length(object@snp.id), "\n")
        cat(" number of associated genes:", length(unique(object@gene)), "\n")
        invisible(NULL)
    }
)
setMethod(
    "show",
    signature = "geneSet",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat(" Number of gene sets:", length(object@gene.set), "\n")
        if (length(object@gene.set) > 0) {
            rr = range(lapply(object@gene.set, length))
            cat("   ", rr[1], "~", rr[2], " genes within sets\n")
        }
        invisible(NULL)
    }
)
