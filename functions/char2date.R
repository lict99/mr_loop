#' @title Convert character vector to date
#'
#' @param x A character vector.
#' @return A date vector.
char2date <- function(x) {
  if (!is.character(x)) stop("Input must be a character vector.", call = FALSE)
  # remove white spaces
  x <- gsub("[ \t\r\n]", "", x)
  # spilt by "-" or "/"
  split <- strsplit(x, "-|/")
  date_char <- vapply(
    split,
    function(x) {
      # check if year is 4 digits
      if (!is.na(x[[1]]) && nchar(x[[1]]) != 4L) {
        stop("Year must be 4 digits.", x[[1]], call. = FALSE)
      }
      if (length(x) == 3L) {
        # if length is 3, then it is a full date
        sprintf(
          "%04d-%02d-%02d",
          as.integer(x[[1]]),
          as.integer(x[[2]]),
          as.integer(x[[3]])
        )
      } else if (length(x) == 2L) {
        # if length is 2, then it is a year-month date
        # set day to 15th
        sprintf(
          "%04d-%02d-15",
          as.integer(x[[1]]),
          as.integer(x[[2]])
        )
      } else {
        # if length is not 2 or 3, then it is not a valid date
        NA_character_
      }
    },
    FUN.VALUE = character(1L)
  )
  return(as.Date(date_char, format = c("%Y-%m-%d")))
}
