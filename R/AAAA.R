## run this code at package load time


# using an environment to maintain state
classic.env <- new.env(hash = TRUE, parent = parent.frame())

#' Newhall classic (monthly inputs) constants
#' 
#' This data set contains constants used for the "classic" (monthly) Newhall model
#' @name constants.classic
#' @docType data
#' @format A list with 6 elements 
#' - `"depletion.order"` numeric length 64; depletion order
#' - `"depletion.req"` numeric length 64; depletion requirement
#' - `"temp.bins"` numeric length 24; temperature bins
#' - `"pe.bins"` numeric length 24; potential evapotransipration bins
#' - `"knorth"` matrix 31 x 12 
#' - `"ksouth"` matrix 13 x 12 
#'
#' @keywords datasets
NULL
