

#' @title Evaluate soil moisture condition
#'
#' @description Evaluate the cells of a soil moisture matrix to determine a soil moisture "condition".
#'
#' @param soil.profile numeric vector, typically derived from an 8x8 matrix
#' @param n.row integer, number of rows in the soil moisture matrix
#' @param n.col integer, number of columns in the soil moisture matrix
#'
#' @return integer representation of soil moisture condition:
#'   * dry: 1
#'   * dry/moist: 2
#'   * moist: 3
#'
#' @export
#'
moistureCondition <- function(soil.profile, n.row = 8, n.col = 8) {

  # init with empty result, to catch possible errors in logic
  res <- NA

  # sanity check: vector length
  if(length(soil.profile) != (n.row*n.col)) {
    stop('incorrect input length', call. = FALSE)
  }

  # critical positions in the soil.profile moisture matrix
  # for now, hard-coded to classic 8x8 dimensions
  #
  # [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8]
  # [1,]    1    2    3    4    5    6    7    8
  # [2,]    9   10   11   12   13   14   15   16
  # [3,]   17   18   19   20   21   22   23   24
  # [4,]   25   26   27   28   29   30   31   32
  # [5,]   33   34   35   36   37   38   39   40
  # [6,]   41   42   43   44   45   46   47   48
  # [7,]   49   50   51   52   53   54   55   56
  # [8,]   57   58   59   60   61   62   63   64

  # critical positions in the 8x8 matrix are:
  # 9, 17, 25

  # generalize to any dimension matrix
  m <- matrix(1:(n.row * n.col), nrow = n.row, byrow = TRUE)

  # column 1
  idx.1 <- m[2, 1]  # row 2
  idx.2 <- m[3, 1]  # row 3
  idx.3 <- m[4, 1]  # row 4

  # combined index
  idx <- c(idx.1, idx.2, idx.3)

  ## logic flows sequentially, and cannot result in more than one TRUE eval

  ## dry: all critical positions == 0
  if(sum(soil.profile[idx]) == 0) {
    res <- 1
  }

  ## dry/moist: some but not all critical positions > 0
  if(any(soil.profile[idx] > 0) & !all(soil.profile[idx] > 0)) {
    res <- 2
  }

  ## moist: all critical position > 0
  if(all(soil.profile[idx] > 0)) {
    res <- 3
  }

  # TODO: sanity checks
  if(is.na(res))
    stop('how did this happen?')


  # TODO: encode as factor
  return(res)

}



## reference: original implementation had an error: 2 would always replace 1

# if (sum(soil.profile[c(9,17,25)]) == 0){
#   print(1)
#   current.moisture.condition <- 1
# }
#
# if (soil.profile[9] > 0 & soil.profile[17] > 0 & soil.profile[25] > 0){
#   print(3)
#   current.moisture.condition <- 3
# }
#
# if (soil.profile[9] == 0 | soil.profile[17] == 0 | soil.profile[25] == 0 & (sum(soil.profile[c(9,17,25)]) > 0)){
#   print(2)
#   current.moisture.condition <- 2
# }



