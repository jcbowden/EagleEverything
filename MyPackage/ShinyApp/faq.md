## General 

### Can Eagle deal with inbreds or can it only handle data on outbred individuals?
Eagle can handle data recorded on inbred or outbred individuals. 

### What is the relationship between R and Eagle
Eagle is an R package, that makes heavy use of c++. 


### How should I organize my input data
Organize your input data into three files; a file with the marker data, 
a file with the phenotypic data, and a file with the map.
The marker file has rows which correspond to data on individuals and columns which correspond to data on snp. 
The phenotype file has rows which correspond to data on individuals and columns which correspond to data on traits 
and fixed effects. 
The map file has rows that correspond to snp. 
Each file is read into Eagle separately. 

### Do I need a marker map
Eagle will run without a marker map being supplied. 
Unlike classic linkage mapping, association mapping doesn't require a known marker map. 
However, for interprability of the results, it 
is best to also supply a marker map. 


### Will Eagle check for errors
Eagle does its best to capture errors, especially in the input files.
Eagle also issues error messages that we hope will be helping for solving the problem. 


## Using OpenUI

### I cannot see the file browser when I click on Choose File?
This has been a source of frustration.  The file browser has opened, it is sitting in the background. 
This behaviour of the R file browser opening behind instead of in front of the active window is a known problem which has yet to be solved. 


### When I use OpenUI, nothing happens
Make sure you are running the latest version of R along with updated versions of the packages (use `update.packages`). 


### When I use OpenUI, my web browser opens but everything is greyed out and I cannot click on any of the widgets
Go back to the window from which you issued the `OpenUI` command.  If an error has occurred, it will be printed here. 
The most likely cause is that your packages have been installed under a version of R that is different to the one being used. 
Run the following to update all your installed packages 

```r
update.packages(checkBuilt=TRUE, ask=FALSE)
```

  



## Reading in Marker Genotypes

### What types of marker data can Eagle handle?
Eagle can deal with two types of marker data; genotypic and allelic. If marker genotypes are available, then they need to be in 
a plain space separated text file where the rows are the individuals and the columns are the snps. The file should not contain row
or column headings. The genotypes can be any alphanumeric value, as long as these same alphanumeric values are used across all the loci. 
If allelic data are available, then this needs to be in PLINK ped form. 


### Am I allowed to have missing marker data
Yes. If genotype data are available, then any alphanumeric value can be used to denote missing data. If allelic data are available, 
PLINK only allows 0 or - to be missing alleles. 

### If I have missing marker genoytypes, what happens?
Eagle sets the missing marker genotypes to AB. Since Eagle assumes an additive locus model, setting a missing genotype 
to AB is equivalent to imputing a genotype that has  no effect on a trait. 
If there are a large number of missing genotypes, this will reduce the power for detecting marker-trait associations. 
A better strategy is to  impute the missing marker genotypes prior to analysis with Eagle with 
dedicated imputation software such as BEAGLE or fastPHASE. 



### Can Eagle handle huge marker data
Yes. Eagle can analyse marker data larger than the memory capacity of a computer by using out-of-memory matrix calculation. 

### I cannot get my marker data to read in
Make sure a line does not begin with a space. Also check that each line has the same number of entries. We've also encountered 
problems when transferring files from a windows system to a unix system and vise versa.  
This is because the format for a windows and a unix text file differ slightly in how they handle the end of a line. When transferring 
files between different platforms, it is good practice to use a file conversion program first, such as dos2unix and unix2dos. 


### Can Eagle handle dominant marker loci or loci with more than two alleles?
Yes. Suppose you have dominant marker data with genotype codes 0 and 1 for absence and presence, respectively. 
Then, when the marker data is being read with `Read Genotypes`, set the parameter AA to 0  and BB to 1 (or vise versa) but 
leave AB blank.   If you have a multi-allelic locus, say with 10 segregating alleles, then turn this locus into 10 
dominant loci and treat these loci as described. 


## Reading in Phenotypic Data

### Can Eagle handle missing trait data and/or missing fixed effects data?
Yes. Individuals with missing trait and/or fixed effects are removed from the analysis. 
However, only individuals whose data are being considered for analysis will be removed. Individuals with missing data not 
involved in the analysis will not be removed. 




### Can my phenotypic file contain multiple traits
Yes.

### Can my phenotypic file contain fixed effects that I may not use in an analysis
Yes

### Do I need to be concerned with the order of the rows in my phenotypic file
Yes. The rows are ordered by individual in the phenotypic file.  This same ordering of rows by individual must be followed for the 
marker file. Eagle will check to make sure the total number of rows in the genotypic and phenotypic files match. However, 
Eagle  cannot check that the ordering of the rows in both files is the same. 



## Running Eagle


### Where can I get the latest version of Eagle? 
If Eagle has not been previously installed, then start  R and at the R prompt, type

```r
install.packages("Eagle", dependencies=TRUE)
```
Eagle is dependent upon several other packages. This command will install Eagle, along with any missing packages upon which Eagle is dependent. 

If Eagle is already installed, but you you want the latest version of Eagle, then at the R prompt, type

```r
update.packages(checkBuilt=TRUE, ask=FALSE)
```

### Do I need to update Eagle if I have updated my version of R
Yes. 

We have seen some strange behaviour, especially with `OpenUI()`, when Eagle and its dependencies 
have been installed under different versions of R. When ever a newer version of R is installed, it is good practice to 
update it's packages with  

```r
update.packages(checkBuilt=TRUE, ask=FALSE)
```


### What platforms will Eagle run on? 
Eagle will run on the same platforms that R will run on which are Linux, OS X (Mac), and Windows. 




## Performing Genome-wide association mapping 

### Am I restricted to only analysing continuous traits or can I also analyse discrete/disease traits?
Yes, but there may be a loss in power for detecting marker-trait associations. 
Eagle is based on linear mixed models. Linear mixed models assume a response (or trait) is normally distributed. However, they are 
robust to violations of their assumptions. 

### Am I allowed to have interactions between my fixed effects
Yes, but using Eagle via OpenUI can only handle main effects. One solution is to add an 
extra column to the phenotypic data file that is the interaction of the two effects that are of interest. Then,  include these data 
as an extra fixed effect in the analysis. Alternately, from the R prompt, use the `AM` function and define the fixed effect part of the model 
by setting fformula. It would be best to look online for information on how to specify 
linear models in R which include interaction terms. 

### My analysis doesn't seem to be  using multiple threads/cpu
Not all parts of an Eagle analysis is parallelized. However, you should be seeing multiple threads being used at different times
throughout each iteration of the model. It is most likely that your version of R is not making use of a multi-threaded BLAS library 
such as MKL or openBLAS. Look at our installation notes on how to install R that is  multi-threaded.  




## Help

### Where can I go for help?  
* Help is available through this user interface. Just move your mouse pointer over a widget and a pop-up window will appear
  giving additional information on its use. 
* Click the Help tab at the top of a page. A sub-menu will  appear with About, FAQ, and Documents. 
* A quick start guide is available under Help/Documents. It can also be accessed by `vignette("QuickStart", "Eagle")` at the R prompt. 
* We have produced a Youtube video where we run through a complete analysis with Eagle. It can be found on Youtube by searching for 
"How to use Eagle for genome-wide association mapping".
* If all else has failed,  email us on <eaglehelp@csiro.au>




