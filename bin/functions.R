#### custom R functions to be used in the pipeline

## checks existence of a variable and print error message if not defined
nf_var_check <- function(x) {
    if (!exists(x)) {
        stop(paste0("The variable'",x,"' is not defined! Make sure to check the Nextflow process inputs."))
    } else {
        print(paste0("Input variable '",x,"' = ",eval(parse(text = x))))
    }
}