#### define helper method for parsing gct files ###

#' Adjust the data types for columns of a meta data frame
#' 
#' @description GCT(X) parsing initially returns data frames
#'   of row and column descriptors where all columns are of
#'   type character. This is inconvenient for analysis, so
#'   the goal of this function is to try and guess the
#'   appropriate data type for each column.
#'   
#' @param meta a data.frame
#' 
#' @details This is a low-level helper function
#'   which most users will not need to access directly
#' 
#' @return meta the same data frame with (potentially) adjusted
#'   column types.
#'   
#' @examples 
#' # meta data table with all character types
#' str(col_desc)
#' fixed <- mapGCT:::fix.datatypes(col_desc)
#' # note how some column classes have changed
#' str(fixed)
#' 
#' @keywords internal
fix.datatypes <- function(meta) {
  for (field.name in names(meta)) {
    # get the field values
    field <- meta[[field.name]]
    field.except.na <- field[!is.na(field)]
    # check if it's numeric. data may come in as a string
    # but actually contains numeric values. if so, as.numeric
    # will not result in a vector of NA values
    field.as.numeric <- suppressWarnings(as.numeric(field.except.na))
    if (!any(is.na(field.as.numeric))) {
      field <- as.numeric(field)
    }
    if (is.numeric(field)) {
      # check if it's an integer. data may be floats but
      # if we coerce to an integer and the difference from
      # original values is zero, that means data are actually
      # integers. integer conversion will return NA if there
      # are any issues.
      field.except.na <- field[!is.na(field)]
      field.as.integer <- suppressWarnings(as.integer(field.except.na))
      if (!any(is.na(field.as.integer))) {
        # integer conversion was fine, lets see if the
        # values are altered
        diffs <- field.except.na - field.as.integer
        if (all(diffs == 0)) {
          # converting to integer didn't change value,
          # set field to integer values
          field <- as.integer(field)
        }
      }
    }
    # insert back into the annotations
    meta[[field.name]] <- field
  }
  return(meta)
}


#### define helper method for parsing gct files ###

#' Check whether \code{test_names} are columns in the \code{\link{data.frame}} df
#' @param test_names a vector of column names to test
#' @param df the \code{\link{data.frame}} to test against
#' @param throw_error boolean indicating whether to throw an error if
#'   any \code{test_names} are not found in \code{df}
#' @return boolean indicating whether or not all \code{test_names} are
#'   columns of \code{df}
#' @examples 
#' check_colnames(c("pert_id", "cell_id"), col_desc)            # TRUE
#' check_colnames(c("pert_id", "foo"), col_desc, throw_error=FALSE) # FALSE, suppress error
#' 
#' @export
check_colnames <- function(test_names, df, throw_error=T) {
  # check whether test_names are valid names in df
  # throw error if specified
  diffs <- setdiff(test_names, names(df))
  if (length(diffs) > 0) {
    if (throw_error) {
      stop(paste("the following column names are not found in", deparse(substitute(df)), ":",
                 paste(diffs, collapse=" "), "\n"))
    } else {
      return(F)
    }
  } else {
    return(T)
  }
}


#' Merge two \code{\link{data.frame}}s, but where there are common fields
#' those in \code{x} are retained and those in \code{y} are dropped.
#' 
#' @param x the \code{\link{data.frame}} whose columns take precedence
#' @param y another \code{\link{data.frame}}
#' @param by a vector of column names to merge on
#' @param allow.cartesian boolean indicating whether it's ok
#'   for repeated values in either table to merge with each other
#'   over and over again.
#' @param as_data_frame boolean indicating whether to ensure
#'   the returned object is a \code{\link{data.frame}} instead of a \code{\link{data.table}}.
#'   This ensures compatibility with GCT object conventions,
#'   that is, the \code{\link{rdesc}} and \code{\link{cdesc}} slots must be strictly
#'   \code{\link{data.frame}} objects.
#'   
#' @return a \code{\link{data.frame}} or \code{\link{data.table}} object
#' 
#' @examples 
#' (x <- data.table(foo=letters[1:10], bar=1:10))
#' (y <- data.table(foo=letters[1:10], bar=11:20, baz=LETTERS[1:10]))
#' # the 'bar' column from y will be dropped on merge
#' mapGCT:::merge_with_precedence(x, y, by="foo")
#'
#' @keywords internal
#' @seealso data.table::merge
merge_with_precedence <- function(x, y, by, allow.cartesian=T,
                                  as_data_frame = T) {
  trash <- check_colnames(by, x)
  trash <- check_colnames(by, y)
  # cast as data.tables
  x <- data.table(x)
  y <- data.table(y)
  # get rid of row names
  setattr(x, "rownames", NULL)
  setattr(y, "rownames", NULL)
  common_cols <- intersect(names(x), names(y))
  y_keepcols <- unique(c(by, setdiff(names(y), common_cols)))
  y <- y[, y_keepcols, with=F]
  # if not all ids match, issue a warning
  if (!all(x[[by]] %in% y[[by]])) {
    warning("not all rows of x had a match in y. some columns may contain NA")
  }
  # merge keeping all the values in x, making sure that the
  # resulting data.table is sorted in the same order as the 
  # original object x
  merged <- merge(x, y, by=by, allow.cartesian=allow.cartesian, all.x=T)
  if (as_data_frame) {
    # cast back to a data.frame if requested
    merged <- data.frame(merged)
  }
  return(merged)
}

