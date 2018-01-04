

#' An S4 class to represent a GCT object
#' 
#' @slot mat a numeric matrix
#' @slot rid a character vector of row ids
#' @slot cid a character vector of column ids
#' @slot rdesc a \code{data.frame} of row descriptors
#' @slot cdesc a \code{data.frame} of column descriptors
#' @slot version version of the gct
#' 
#' @description The GCT_object class serves to represent annotated matrices.
#' The \code{mat} slot contains numeric matrics and the \code{rdesc} and \code{cdesc}
#' slots contain data frames with annotations about the rows and columns, respectively
setClass("GCT_object",
         representation(
           mat = "matrix",
           rid = "character",
           cid = "character",
           rdesc = "data.frame",
           cdesc = "data.frame",
           version = "character",
           src = "character"
         )
)

# define the validity for GCT_object class
setValidity("GCT_object",
            function(object) {
              # check whether dimensions of various
              # slots are in sync
              nrows <- nrow(object@mat)
              ncols <- ncol(object@mat)
              if (nrows != length(object@rid))
                "rid must be the same length as number of matrix rows"
              if (ncols != length(object@cid))
                "cid must be the same length as number of matrix columns"
              if (length(object@cid) > length(unique(object@cid)))
                "cid must be unique"
              if (length(object@rid) > length(unique(object@rid)))
                "rid must be unique"
              if (nrow(object@cdesc) != ncols)
                "cdesc must have same number of rows as matrix has columns"
              if (nrow(object@rdesc) != nrows)
                "rdesc must have same number of rows as matrix has rows"
              else
                T
            }
)

suppressMessages({
  # set method for displaying a GCT_object
  # just use the 'str' function to show its structure
  setMethod("show", signature("GCT_object"), function(object) {
    str(object)
  })
  
  # dim, nrow and ncol to display the # of rows and columns
  # for a GCT_object's matrix
  setMethod("ncol", signature("GCT_object"), function(x) {
    ncol(x@mat)
  })
  setMethod("nrow", signature("GCT_object"), function(x) {
    nrow(x@mat)
  })
  setMethod("dim", signature("GCT_object"), function(x) {
    dim(x@mat)
  })
  setMethod("range", signature("GCT_object"), function(x, na.rm=F, finite=F) {
    range(x@mat, na.rm=na.rm, finite=finite)
  })
  setMethod("max", signature("GCT_object"), function(x, na.rm=F) {
    max(x@mat, na.rm=na.rm)
  })
  setMethod("min", signature("GCT_object"), function(x, na.rm=F) {
    min(x@mat, na.rm=na.rm)
  })
  setMethod("diag", signature("GCT_object"), function(x) {
    diag(x@mat)
  })
})


# define the initialization method for the GCT_object class
setMethod("initialize",
          signature = "GCT_object",
          definition = function(.Object, mat = NULL, rdesc = NULL, cdesc = NULL,
                                version = c("#1.3", "#1.2"), src = NULL) {
            if (!is.null(mat)) {
              .Object@mat = mat
              .Object@rid = rownames(mat)
              .Object@cid = colnames(mat)
            }
            if (!is.null(rdesc)) {
              if (nrow(rdesc) != nrow(mat)) {
                stop("rdesc must have same number of rows as matrix has rows")
              }else if (all(rownames(rdesc) %in% rownames(mat))) {
                rdesc <- rdesc[rownames(mat), ]
                .Object@rdesc <- rdesc
              }else {stop("rownames of rdesc are not same as rownames of matrix")}
            }else{
              .Object@rdesc = data.frame()
            }
            
            if (!is.null(cdesc)) {
              if (nrow(cdesc) != ncol(mat)) {
                stop("cdesc must have same number of rows as matrix has columns")
              }else if (all(rownames(cdesc) %in% colnames(mat))) {
                cdesc <- cdesc[colnames(mat), ]
                .Object@cdesc <- cdesc
              }else {stop("rownames of cdesc are not same as colnames of matrix")}
            }else{
              .Object@cdesc = data.frame()
            }
            
            if (!is.null(src)) {
              .Object@src <- src
            }else {
              version <- match.arg(version)
              if (version == "#1.3") {
                if (is.null(cdesc)) {
                  stop("Column descriptors not provided, provide cdesc or set version to #1.2")
                }else{
                  .Object@version <- version
                }
              }else{
                .Object@version <- version
              }
            }
            return(.Object)
          }
)

#' Create a GCT_object for a given numeric matrix
#' 
#' @param mat a numeric matrix
#' @param rdesc a \code{data.frame} of row descriptors
#' @param cdesc a \code{data.frame} of column descriptors
#' @param version version of the gct
#' 
#' @description create a GCT_object for a given matrix and cloumn description
#' @export
to_GCT <- function(mat, cdesc=NULL, rdesc=NULL, version = NULL) {
  ds <- new("GCT_object",
            mat = mat,
            rdesc = rdesc,
            cdesc = cdesc,
            version = version)
  return(ds)
}


#' Write a GCT_object to disk in GCT format
#' 
#' @param ds the GCT_object to be written
#' @param ofile the desired output filename
#' @param precision the numeric precision at which to
#'   save the matrix. See \code{details}.
#' @param appenddim boolean indicating whether to append
#'   matrix dimensions to filename
#' @param ver the GCT version to write. See \code{details}.
#' 
#' @details Since GCT is text format, the higher \code{precision}
#'   you choose, the larger the file size.
#'   \code{ver} is assumed to be 3, aka GCT version 1.3, which supports
#'   embedded row and column metadata in the GCT file. Any other value
#'   passed to \code{ver} will result in a GCT version 1.2 file which
#'   contains only the matrix data and no annotations.
#'
#' @return NULL
#' 
#' @examples 
#' \dontrun{
#' write.gct(ds, "dataset", precision=2)
#' }
#' @family GCT parsing functions
#' @export
write_gct <- function(ds, ofile, precision=4, ver=3) {
  if (!class(ds)=="GCT_object") {
    stop("ds must be a GCT_object")
  }
  
  precision = floor(precision)
  cat(sprintf('Saving file to %s\n',ofile))
  nr <- nrow(ds@mat)
  nc <- ncol(ds@mat)
  cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
  cat(sprintf('Setting precision to %d\n',precision))
  
  # open file and write   
  if (ver==3) {
    # remove the 'id' columns
    ds@cdesc$id <- NULL
    ds@rdesc$id <- NULL
    # get the counts of meta data fields
    nrdesc = dim(ds@rdesc)[2]
    ncdesc = dim(ds@cdesc)[2]
    colkeys = colnames(ds@cdesc)
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id',colnames(ds@rdesc),ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    # line 4 + ncdesc: sample desc
    filler = 'na'
    if (ncdesc > 0) {
      for (ii in 1:ncdesc) {
        if (is.numeric(ds@cdesc[,ii])) {
          cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                      round(ds@cdesc[,ii],precision)),
                    collapse='\t'),
              file=ofile,sep='\n',append=T)  
        } else {
          cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                      ds@cdesc[,ii]),
                    collapse='\t'),
              file=ofile,sep='\n',append=T)
        }
      }
    }
    
    for (ii in 1:nr) {    
      # print rows
      cat(paste(c(ds@rid[ii],
                  ds@rdesc[ii,],
                  round(ds@mat[ii,],precision)),collapse='\t'),
          sep='\n',file=ofile,append=T)
    }
  } else {
    # assume ver 1.2 and below, ignore descriptors
    # append header
    cat(sprintf('#1.%d\n%d\t%d', ver, nr, nc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id','Description',ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    
    for (ii in 1:nr) {    
      # print rows
      cat(paste(c(ds@rid[ii],
                  ds@rdesc[ii, 2],
                  round(ds@mat[ii,],precision)),collapse='\t'),
          sep='\n',file=ofile,append=T)
    }
  }
  
  cat(sprintf('Saved.\n'))  
}


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
#' \dontrun{
#' str(cdesc_char)
#' fixed <- mapGCT:::fix.datatypes(cdesc_char)
#' str(fixed)
#' }
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
#' \dontrun{
#' check_colnames(c("column1", "column2"), df)            # TRUE
#' check_colnames(c("column1", "not_column"), df, throw_error=FALSE) # FALSE, suppress error
#' }
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
#' cmapR:::merge_with_precedence(x, y, by="foo")
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


#' Parse a GCT file into the workspace as a GCT_object
#' 
#' @param fname path to the GCT file on disk
#' @param set_annot_rownames boolean indicating whether to set the
#'   rownames on the row/column metadata data.frames. Set this to 
#'   false if the GCT file has duplicate row/column ids.
#' @param matrix_only boolean indicating whether to parse only
#'   the matrix (ignoring row and column annotations)
#'
#' @details \code{parse_gct} also supports parsing of plain text
#'   GCT files, so this function can be used as a general GCT parser.
#' 
#' @examples 
#' \dontrun{
#' ds <- parse_gct("path/to/gct/file")
#' # matrix only
#' ds <- parse_gct("path/to/gct/file", matrix_only=TRUE)
#' }
#' 
#' @family GCT parsing functions
#' @export
parse_gct <- function(fname, set_annot_rownames=F, matrix_only=F) {
  ds <- new("GCT_object", src = fname)
  
  if (!(grepl(".gct$", fname)))
    stop("A .gct file must be given")
  # get the .gct version by reading first line
  ver = scan(fname, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
  # get matrix dimensions by reading second line
  dimensions = scan(fname, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
  nrmat = dimensions[1]
  ncmat = dimensions[2]
  if (length(dimensions)==4) {
    # a #1.3 file
    message("parsing as GCT v1.3")
    nrhd <- dimensions[3]
    nchd <- dimensions[4]
  }else{
    # a #1.2 file
    message("parsing as GCT v1.2")
    nrhd <- 0
    nchd <- 0
  }
  message(paste(fname, nrmat, "rows,", ncmat, "cols,", nrhd, "row descriptors,", nchd, "col descriptors"))
  # read in header line
  header = scan(fname, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
  # construct row header and column id's from the header line
  if ( nrhd > 0 ) {
    rhd <- header[2:(nrhd+1)]
    cid <- header[-(nrhd+1):-1]
    col_offset <- 1
  }else {
    if (any(grepl("description", header, ignore.case=T))) {
      # check for presence of description column in v1.2 files
      col_offset <- 2
    } else {
      col_offset <- col_offset <- 1
    }
    rhd = NULL
    cid = header[(1+col_offset):length(header)]
  }
  # read in the next set of headers (column annotations) and shape into a matrix
  if ( nchd > 0 ) {
    header = scan(fname, what = "", nlines = nchd, skip = 3, sep = "\t", 
                  quote = NULL, quiet = TRUE)		
    header = matrix(header, nrow = nchd, 
                    ncol = ncmat + nrhd + 1, byrow = TRUE)
    # extract the column header and column descriptions
    chd = header[,1]
    cdesc = header[,-(nrhd+1):-1]
    # need to transpose in the case where there's only one column annotation
    if ( nchd == 1 )
      cdesc = t(cdesc)
  }else {
    chd = NULL
    cdesc = data.frame()
  }
  # read in the data matrix and row descriptions, shape into a matrix
  mat = scan(fname, what = "", nlines = nrmat, 
             skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
  mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + col_offset, 
               byrow = TRUE)
  # message(paste(dim(mat), collapse="\t"))
  # Extract the row id's row descriptions, and the data matrix
  rid = mat[,1]
  if ( nrhd > 0 ) {
    # need as.matrix for the case where there's only one row annotation
    rdesc = as.matrix(mat[,2:(nrhd + 1)])
    mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]), nrow = nrmat, ncol = ncmat)
  }else {
    rdesc = data.frame()
    mat = matrix(as.numeric(mat[, (1+col_offset):ncol(mat)]), nrow = nrmat, ncol = ncmat)
  }
  # assign names to the data matrix and the row and column descriptions
  # message(paste(dim(mat), collapse="\t"))
  dimnames(mat) = list(rid, cid)
  if ( nrhd > 0 ) {
    dimnames(rdesc) = list(rid,rhd)
    rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
  }
  if ( nchd > 0 ) {
    cdesc = t(cdesc)
    dimnames(cdesc) = list(cid,chd)
    cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
  }
  # assign to the GCT slots
  ds@mat = mat
  ds@rid = rownames(mat)
  ds@cid = colnames(mat)
  if (!matrix_only) {
    # return annotations as well as matrix
    ds@rdesc = fix.datatypes(rdesc)
    ds@cdesc = fix.datatypes(cdesc)
    # add id columns to rdesc and cdesc
    ds@rdesc$id <- rownames(ds@rdesc)
    ds@cdesc$id <- rownames(ds@cdesc)
  }
  return(ds)
}


#' Add annotations to a GCT_object
#' 
#' @description Given a GCT_object and either a \code{\link{data.frame}} or
#' a path to an annotation table, apply the annotations to the
#' gct using the given \code{keyfield}.
#' 
#' @param ds a GCT_object
#' @param annot a \code{\link{data.frame}} of annotations
#' @param dimension either 'row' or 'column' indicating which dimension
#'   of \code{g} to annotate
#' @param keyfield the character name of the column in \code{annot} that 
#'   matches the row or column identifiers in \code{g}
#'   
#' @return a GCT_object with annotations applied to the specified
#'   dimension
#'   
#' @examples 
#' \dontrun{
#'  ds <- parse_gct('/path/to/gct/file')
#'  ds <- annotate_gct(ds, '/path/to/annot')
#' }
#' 
#' @family GCT utilities
#' @export
annotate_gct <- function(ds, annot, dimension = "column", keyfield="id") {
  if (!class(ds)=="GCT_object") {
    stop("ds must be a GCT_object")
  }
  if (!class(annot) == "data.frame") {
    stop("annot must be a data.frame")
  }else{
    # convert to data.table
    annot <- data.table(annot)
  }
  # convert the keyfield column to id for merging
  # assumes the gct object has an id field in its existing annotations
  if (!(keyfield %in% names(annot))) {
    stop(paste("column", keyfield, "not found in annotations"))
  }
  # rename the column to id so we can do the merge
  annot$id <- annot[[keyfield]]
  if (dimension == "column") dimension <- "col"
  if (dimension == "row") {
    orig_id <- ds@rdesc$id
    merged <- merge_with_precedence(ds@rdesc, annot, by="id", allow.cartesian=T,
                                              as_data_frame=T)
    idx <- match(orig_id, merged$id)
    merged <- merged[idx, ]
    ds@rdesc <- merged
  }else if (dimension == "col") {
    orig_id <- ds@cdesc$id
    merged <- merge_with_precedence(ds@cdesc, annot, by="id", allow.cartesian=T,
                                              as_data_frame=T)
    idx <- match(orig_id, merged$id)
    merged <- merged[idx, ]
    ds@cdesc <- merged
  }else{
    stop("dimension must be either row or column")
  }
  return(ds)
}
