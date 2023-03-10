# mapGCT package

**mapGCT** is a simple package to read and write gct file. It can also be used to add extra row or column descriptors to the existing GCT_object.

### Install instructions

Dependencies are listed in `DESCRIPTION`

**Installing from source**

The latest version can be installed directly from bitbucket using
`devtools::install_github("ElucidataInc/mapGCT")`.

Alternatively, R's `install.packages` function can be used at a tarball of the `mapGCT` archive. This archive can be generated by cloning this repository and doing the following:

	# make a gzip tar ball of the repo
	R CMD build mapGCT
	# makes mapGCT_x.x.x.tar.gz
	
	# check that the package is ok
	R CMD check mapGCT_x.x.x.tar.gz	

Once tarball is created, execute following command in R terminal:

	install.packages("mapGCT_x.x.x.tar.gz", type="source", repos=NULL)
	library("mapGCT")


# Code Examples

There are four main functions for parsing and writing gct files: 

  1. `to_GCT()` Converts the given matrix, column descriptions and row descriptions to a GCT_object
  2. `write_gct` write a GCT_object to disk in gct format
  3. `parse_gct` parse or load a gct file into the workspace as a GCT_object
  4. `annotate_gct` Add row or column annotations to GCT_object. Can be used to add add extra columns to gct_object

### Usages
These functions can be called as follows:
```R
# Convert to GCT_object class
ds <- to_GCT(mat = gct_mat, cdesc = col_desc, rdesc = row_desc)

# Write GCT_object to disk
write_gct(ds, "dataset.gct", precision=2)

# Read gct file in R
gct_file <- system.file("extdata", "example_n50x100.gct", package="mapGCT")
ds <- parse_gct(gct_file)

# Add column annotations
ds <- annotate_gct(ds, col_annotate) 
```