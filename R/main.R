

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
           version = "character"
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

# define the initialization method for the GCT_object class
setMethod("initialize",
          signature = "GCT_object",
          definition = function(.Object, mat, rdesc = NULL, cdesc = NULL,
                                version = c("#1.3", "#1.2")) {
            .Object@mat = mat
            .Object@rid = rownames(mat)
            .Object@cid = colnames(mat)
            
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
            return(.Object)
          }
)

#' Create a GCT_object for a given numeric matrix
#' 
#' @param mat a numeric matrix
#' @param rid a character vector of row ids
#' @param cid a character vector of column ids
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
#' @family GCTX parsing functions
#' @export

write_gct <- function(ds, ofile, precision=4, appenddim=F, ver=3) {
  if (!class(ds)=="GCT_object") {
    stop("ds must be a GCT_object")
  }
  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gct")
  
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


