#' Package documentation
#'
#' @name Eagle-package
#' @title Multiple-locus Genome-Wide Association Mapping
#' @docType package
#' @author Andrew W. George \email{andrew.george@@csiro.au} with lots of help from Joshua Bowden (CSIRO). 
#' @description A package for making genome-wide association mapping with multiple-locus linear mixed models routine.
#' @details  
#' @section Motivation:  Data from genome-wide association studies are analyzed, commonly, with single-locus 
#' models. That is,  analyzes are performed on a locus-by-locus basis. Multiple-locus approaches 
#' that model the association between a trait and multiple loci simultaneously are more powerful. However, 
#' these methods do not scale well with study size and many of the packages that implement these methods are not easy 
#' to use. 
#' Eagle was specifically designed to make 
#' genome-wide association mapping with multiple-locus models simple and practical.  Much effort has been 
#' devoted to making the package as easy to use as possible. As part of this effort, we 
#' developed a web-based user interface to Eagle.
#'
#' @section Assumptions:
#' \enumerate{
#' \item Individuals are diploid but they can be inbred or outbred.
#' \item The marker and phenotype data are in separate files.
#' \item The rows in the marker and phenotype file correspond to data collected on the same individuals.
#' \item  Marker loci are snps. Dominant and multi-allelic loci will need to be converted into biallelic (snp-like) loci. 
#' \item The trait is continuous and normally distributed. Eagle can handle non-normally distributed trait data but 
#'  there may a loss of power
#' to detect marker-trait associations. 
#' }
#' @section Reading in  Marker Data:
#' Eagle can handle marker data in the form of genotypes or alleles. If genotypes have been collected, then 
#' the input file must be a plain space separated text file. If allelic information have been collected on the
#' marker loci, then the input file must be in the form of a PLINK ped file. 
#' These files can be 
#' larger than the memory capacity of the machine. Other formats 
#' can also be handled via the  PLINK package and using  the \code{recode} option to 
#' convert different format marker files into ped files.  See \code{\link{ReadMarker}} for more details.
#'
#' @section Quick Start Guide:
#' At the R prompt, run \preformatted{OpenGUI()} This opens a web browser to the
#' user interface for Eagle. 
#'
#' @section Output: The aim of a genome-wide association study (GWAS) is to identify those marker loci 
#' closest to the genes that are influencing a trait. So, when the GWAS data are 
#' analysed, a set of marker loci labels are returned as the results These marker 
#'  loci are closest to the genes underlying the trait and are found while simultaneously 
#' accounting for other marker-trait associations, familial relatedness, and 
#' fixed effects such as population structure.
#' More detailed output such as the additive effect of the marker locus, its 
#' significance in the multiple-locus model ( measured by a p-value), and 
#' an estimate of the amount of variation explained by the locus can be 
#' obtained by running the summary function \code{\link{SummaryAM}}.
#'
#'@section Where to get help: 
#' A variety of different help options are available. 
#' \itemize{
#'  \item At the R prompt, type  \preformatted{library(, "Eagle")} for an overview of the package and its functions. 
#' \item For detailed help on a function called "foo" say, type  \preformatted{help("foo")} 
#' \item Visit the Eagle website at \url{http://eagle.r-forge.r-project.org/} where 
#' you can find a quick start guide, instructions on getting the most out of Eagle, 
#' tutorials, and other useful information. 
#'}
#' @keywords Association mapping, multiple-locus models, linear mixed models.
NULL




