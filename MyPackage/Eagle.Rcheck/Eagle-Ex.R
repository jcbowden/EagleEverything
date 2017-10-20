pkgname <- "Eagle"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "Eagle-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('Eagle')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("AM")
### * AM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: AM
### Title: multiple-locus association mapping
### Aliases: AM

### ** Examples

  ## Not run: 
##D  
##D   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
##D   # However, the code can be run by the user. 
##D   #
##D 
##D   #-------------------------
##D   #  Example  
##D   #------------------------
##D 
##D   # read the map 
##D   #~~~~~~~~~~~~~~
##D   
##D   # File is a plain space separated text file with the first row 
##D   # the column headings
##D   complete.name <- system.file('extdata', 'map.txt', 
##D                                    package='Eagle')
##D   map_obj <- ReadMap(filename=complete.name) 
##D 
##D   # read marker data
##D   #~~~~~~~~~~~~~~~~~~~~
##D   # Reading in a PLINK ped file 
##D   # and setting the available memory on the machine for the reading of the data to 8  gigabytes
##D   complete.name <- system.file('extdata', 'geno.ped', 
##D                                      package='Eagle')
##D   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
##D  
##D   # read phenotype data
##D   #~~~~~~~~~~~~~~~~~~~~~~~
##D 
##D   # Read in a plain text file with data on a single trait and two covariates
##D   # The first row of the text file contains the column names y, cov1, and cov2. 
##D   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
##D   
##D   pheno_obj <- ReadPheno(filename=complete.name)
##D            
##D 
##D  # Performing multiple-locus genome-wide association mapping with a model 
##D  #    with no fixed effects except for an intercept. 
##D  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##D  
##D   res <- AM(trait = 'y',
##D                            fformula=c('cov1+cov2'),
##D                            map = map_obj,
##D                            pheno = pheno_obj,
##D                            geno = geno_obj, availmemGb=8)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("AM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Eagle-package")
### * Eagle-package

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Eagle-package
### Title: A short title line describing what the package does
### Aliases: Eagle-package Eagle
### Keywords: package

### ** Examples

  ## Not run: 
##D      ## Optional simple examples of the most important functions
##D      ## These can be in \dontrun{} and \donttest{} blocks.   
##D   
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Eagle-package", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("OpenGUI")
### * OpenGUI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: OpenGUI
### Title: Browser-based Graphical User Interface
### Aliases: OpenGUI

### ** Examples

## Not run: 
##D # opens a web browser 
##D OpenGUI()
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("OpenGUI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReadMap")
### * ReadMap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReadMap
### Title: Read map file
### Aliases: ReadMap

### ** Examples

# Read in  example map data from ./extdata/

# find the full location of the map data 
complete.name <- system.file('extdata', 'map.txt', package='Eagle')
  
# read in map data 
map_obj <- ReadMap(filename=complete.name) 
                               
# look at first few rows of the map file
head(map_obj)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReadMap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReadMarker")
### * ReadMarker

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReadMarker
### Title: Read marker data.
### Aliases: ReadMarker

### ** Examples

  #--------------------------------
  #  Example 1
  #-------------------------------
  #
  # Read in the genotype data contained in the text file geno.txt
  #
  # The function system.file() gives the full file name (name + full path).
  complete.name <- system.file('extdata', 'geno.txt', package='Eagle')
  # 
  # The full path and name of the file is
  print(complete.name)
  
  # Here, 0 values are being treated as genotype AA,
  # 1 values are being treated as genotype AB, 
  # and 2 values are being treated as genotype BB. 
  # 4 gigabytes of memory has been specified. 
  # The file is space separated with the rows the individuals
  # and the columns the snp loci.
  geno_obj <- ReadMarker(filename=complete.name, type='text', AA=0, AB=1, BB=2, availmemGb=4) 
   
  # view list contents of geno_obj
  print(geno_obj)

  #--------------------------------
  #  Example 2
  #-------------------------------
  #
  # Read in the allelic data contained in the PLINK ped file geno.ped
  #
  # The function system.file() gives the full file name (name + full path).
  complete.name <- system.file('extdata', 'geno.ped', package='Eagle')

  # 
  # The full path and name of the file is
  print(complete.name)
  
  # Here,  the first 6 columns are being ignored and the allelic 
  # information in columns 7 -  10002 is being converted into a reformatted file. 
  # 4 gigabytes of memory has been specified. 
  # The file is space separated with the rows the individuals
  # and the columns the snp loci.
  geno_obj <- ReadMarker(filename=complete.name, type='PLINK', availmemGb=4) 
   
  # view list contents of geno_obj
  print(geno_obj)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReadMarker", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReadPheno")
### * ReadPheno

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReadPheno
### Title: Read phenotype file
### Aliases: ReadPheno

### ** Examples

# Read in  phenotype data from ./extdata/

# find the full location of the phenotype data 
complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')

pheno_obj <- ReadPheno(filename=complete.name)
  
 ## print a couple of lines of the data file
 head(pheno_obj)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReadPheno", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SummaryAM")
### * SummaryAM

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SummaryAM
### Title: Summary of multiple locus association mapping results
### Aliases: SummaryAM

### ** Examples

 ## Not run: 
##D   # Since the following code takes longer than 5 seconds to run, it has been tagged as dontrun. 
##D   # However, the code can be run by the user. 
##D   #
##D 
##D   #---------------
##D   # read the map 
##D   #---------------
##D   #
##D   # File is a plain space separated text file with the first row 
##D   # the column headings
##D   complete.name <- system.file('extdata', 'map.txt', 
##D                                    package='Eagle')
##D   map_obj <- ReadMap(filename=complete.name) 
##D 
##D  # to look at the first few rows of the map file
##D  head(map_obj)
##D 
##D   #------------------
##D   # read marker data
##D   #------------------
##D   # Reading in a PLINK ped file 
##D   # and setting the available memory on the machine for the reading of the data to 8 gigabytes
##D   complete.name <- system.file('extdata', 'geno.ped', 
##D                                      package='Eagle')
##D   geno_obj <- ReadMarker(filename=complete.name,  type='PLINK', availmemGb=8) 
##D  
##D   #----------------------
##D   # read phenotype data
##D   #-----------------------
##D 
##D   # Read in a plain text file with data on a single trait and two fixed effects
##D   # The first row of the text file contains the column names y, cov1, and cov2. 
##D   complete.name <- system.file('extdata', 'pheno.txt', package='Eagle')
##D   
##D   pheno_obj <- ReadPheno(filename=complete.name)
##D            
##D   #-------------------------------------------------------
##D   # Perform multiple-locus genome-wide association mapping 
##D   #-------------------------------------------------------                   
##D   res <- AM(trait = 'y',
##D                            fformula=c("cov1 + cov2"),
##D                            map = map_obj,
##D                            pheno = pheno_obj,
##D                            geno = geno_obj, availmemGb=8)
##D 
##D   #-----------------------------------------
##D   # Produce additional summary information 
##D   #------------------------------------------
##D 
##D   SummaryAM(AMobj=res, pheno=pheno_obj, geno=geno_obj, map=map_obj)
##D  
## End(Not run)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SummaryAM", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
