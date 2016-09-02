#!/usr/bin/env Rscript

DEPENDENCIES <- c('optparse')

#   A function to test if we have package installed
pkgTest <- function(package) {
    if(package %in% rownames(installed.packages()) == FALSE) { # check to see if a packages is available
        install.packages(package) # if not, insall it
    }
}

#   A function to install and load dependent packages
batchInstall <- function(pkgList) {
    options(repos = c(CRAN = "http://cran.rstudio.com")) # set a repo mirror, we used RStudio just because
    for(dep in pkgList) {
        pkgTest(dep) # test to see if the package is installed
    }
    lapply(X = pkgList, FUN = library, character.only = TRUE) # load the packages to be used
}

#   A function to make an argument parser
makeAruments <- function(){
    #   Create a list of options
    options <- list(
        make_option( # Option for reference fata file
            opt_str = c('-f', '--reference-fasta'),
            type = 'character',
            dest = 'reference',
            default = NULL,
            metavar = 'REFERENCE FASTA',
            help = "Reference FASTA for FOF file"
        ),
        make_option( # Option for regions file
            opt_str = c('-r', '--regions'),
            type = 'character',
            dest = 'regions',
            default = NULL,
            metavar = 'REGIONS FILE',
            help = "Regions file for FOF, must be in ANGSD format"
        ),
        make_option( # Option for output name
            opt_str = c('-o', '--output'),
            type = 'character',
            dest = 'output',
            default = paste(getwd(), 'output.fof', sep = '/'),
            metavar = 'OUTPUT',
            help = "Name of output file, defaults to '%default'"
        )
    )
    #   Create a OptionParser object using the list of arguments from above
    parser <- OptionParser(option_list = options, add_help_option = TRUE)
    return(options)
}

#   Write our output file
writeOutput <- function(data, outfile){
    write.table(x = data, file = outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste("FOF file can be found at", outfile))
}

#   Run the program here
main <- function(){
    batchInstall(DEPENDENCIES) # Install dependencies
    #   Parse and check for valid arguments
    parser <- OptionParser(option_list = makeAruments(), add_help_option = TRUE)
    args <- parse_args(object = parser)
    if (! 'reference' %in% names(x = args) | ! 'regions' %in% names(x = args)) { # Required arguments
        print_help(object = parser)
        quit(save = 'no', status = 1)
    }
    if (! file.exists(args$reference) | ! file.exists(args$regions)) { # Make sure our input files exist
        stop("Failed to find either the reference or regions file! Please check to ensure they exist...")
    }
    #   Read in data
    regions <- read.table(file = args$regions, header = FALSE, as.is = TRUE)
    regions$fasta <- replicate(n = nrow(regions), expr = args$reference)
    writeOutput(data = data.frame(regions[[1]], regions$fasta))
}

main()
