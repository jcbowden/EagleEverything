library(shiny)
library(Eagle)
library(tcltk)

read_geno_intro <- function(){
  txt <- "
  Eagle can handle two different types of marker data; genotype data in a plain text space separated file 
  "
  return(txt)
}
read_pheno_intro <- function(){
  txt <- "adfadf"
  
  return(txt)
}  
  
    
    
shinyServer(function(input, output, session){

#  library("markdown")
#  library("knitr")

#rmdfiles <- c("faq.rmd")
#sapply(rmdfiles, knit, quiet = T)

  ##------------------------------------------
  ## Intros to pages
  ##-------------------------------------------
 
#  output$home_intro <- renderText(home_intro())
  output$read_geno_intro <- renderText(read_geno_intro())
  output$read_pheno_intro <- renderText(read_pheno_intro())
  
  ##----------------------------------------
  ##  Read marker path and file name
  ##---------------------------------------- 
  ## upload path and file name
 
  output$choose_marker_file <- renderText(NULL)
  path_to_file <- NULL
  observeEvent(input$choose_marker_file, {
   
       
       path_to_file <<- tryCatch({
            if(.Platform$OS.type=="unix"){
               path_to_file_res <- tk_choose.files()
               print(path_to_file_res)
             } else {
               path_to_file_res <- file.choose()
             }
             return (path_to_file_res)
        }, warning = function(war) {
            print(paste("Eagle::path_to_file Warning: ",war))
            path_to_file_res<-"/R/library/Eagle/shiny_app/shinydata/genoDemo.dat"
            return (path_to_file_res)
        }, error = function(err) {
            print(paste("Eagle::path_to_file Error: ",err))
            path_to_file_res<-"/R/library/Eagle/shiny_app/shinydata/genoDemo.dat"
            return (path_to_file_res)
        }, finally = {
           # path_to_file_res<-"/R/library/Eagle/shiny_app/shinydata/genoDemo.dat"
           return (path_to_file_res)
        }); # END tryCatch
 
        output$choose_marker_file <- renderText( path_to_file )
  })
  



   ## Read marker information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   geno <- NULL
   observeEvent(input$marker_go, {




     if(input$filetype == "plink"){

        withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 geno <<- ReadMarker(filename = path_to_file, type = "PLINK", availmemGb = input$memsize, quiet = FALSE)
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
        })



     }


     if(input$filetype == "text"){
        
             withCallingHandlers({
                 shinyjs::html("ReadMarker", "")
                 aa <- input$AA
                 ab <- input$AB
                 bb <- input$BB
                 missing <- input$missing
                 if(input$AA=="")
                     aa <- NULL
                 if(input$AB=="")
                     ab <- NULL
                 if(input$BB=="")
                     bb <- NULL
                 if(input$missing=="")
                     missing <- NULL
 

                 geno <<- ReadMarker(filename = path_to_file, type = "text", AA = aa, 
                            AB = ab  , BB = bb, availmemGb = input$memsize,  quiet = FALSE , missing=missing) 

              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMarker", html = m$message, add = TRUE)
             })


     }  ## end if(input$filetype == "text")

  })  ## end observeEvent




  ##-------------------------
  ## Read phenotypic data
  ##--------------------------

 ##  Read phenotypic  path and file name
  ## upload path and file name
  path_to_pheno_file <- NULL
  output$choose_pheno_file <- renderText(NULL)
  observeEvent(input$choose_pheno_file, {
    if(.Platform$OS.type=="unix"){
       path_to_pheno_file <<- tk_choose.files()
  print(path_to_pheno_file)
     } else {
       path_to_pheno_file <<- file.choose()

     }


    output$choose_pheno_file <- renderText( path_to_pheno_file )
  })




   ## Read phenotypic  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   pheno <- NULL
   observeEvent(input$pheno_go, {

   header_flag <- FALSE
   if(input$pheno_header == "yes")
      header_flag <- TRUE
   csv_flag <- FALSE
   if(input$pheno_csv == "yes")
      csv_flag <- TRUE

   pheno_missing <- input$pheno_missing
   if(input$pheno_missing=="")
      pheno_missing <- NULL




   withCallingHandlers({
                 shinyjs::html("ReadPheno", "")
                 pheno  <<- ReadPheno(filename = path_to_pheno_file, header=header_flag, csv=csv_flag, missing= pheno_missing)
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadPheno", html = m$message, add = TRUE)
       })





  })  ## end observeEvent


  ##-------------------------##
  ## Read Map                ## 
  ##-------------------------ss

 map <- NULL
 ##  Read map  path and file name
  ## upload path and file name
  path_to_map_file <- NULL
  output$choose_map_file <- renderText(NULL)
  observeEvent(input$choose_map_file, {
    if(.Platform$OS.type=="unix"){
       path_to_map_file <<- tk_choose.files()
        print(path_to_map_file)
     } else {
       path_to_map_file <<- file.choose()

     }

#    rChoiceDialogs::rchoose.files()
    output$choose_map_file <- renderText( path_to_map_file )
  })




   ## Read map  information
   ##~~~~~~~~~~~~~~~~~~~~~~~~~
   map <- NULL
   observeEvent(input$map_go, {

     csv_flag <- FALSE
     if(input$map_csv == "yes")
        csv_flag <- TRUE


   map_header_flag <- FALSE
   if(input$map_header == "yes")
      map_header_flag <- TRUE

       withCallingHandlers({
                 shinyjs::html("ReadMap", "")
                 map  <<- ReadMap(filename = path_to_map_file, csv=csv_flag, header= map_header_flag)
              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "ReadMap", html = m$message, add = TRUE)
       })





  })  ## end observeEvent


  ##-------------------
  ## Analyse Data
  ##-------------------
  traitn <- NULL
  ## gets column names of pheno file
  nms <- reactive({
     if(input$pheno_go && input$pheno_header == "yes")
        return(names(pheno))
     if(input$pheno_go && input$pheno_header == "no")
     {  ## pheno file is not named
        nms <- paste("V", 1:ncol(pheno), sep="")
        return(nms)
     } 


     })


  output$analyse_names <- renderUI({
 #     radioButtons(inputId="nmst", label=h4("Step 1: Choose trait"), 
 #              choices=nms(), inline=TRUE, selected=character(0))
  checkboxGroupInput("nmst", h4("Step 1: Choose trait"), nms(), inline=TRUE)
       
  })  ## end renderUI


  output$analyse_fnames <- renderUI({
      #nms <- names(pheno)
      checkboxGroupInput("nmsf", h4("Step 2: Choose fixed effects"), nms(), inline=TRUE)
    })  ## end renderUI

  fform <- NULL
   output$fmodel <- renderText({
        fform <<- paste(input$nmsf, collapse="+")
   })


   ## how to get traitn and effectsn from UI for later use ??????
   traitn <- reactive({input$nmst})

 res <- NULL
   observeEvent(input$analyse_go, {

       
       withCallingHandlers({
                 shinyjs::html("AM", "")
                 quietvalue <-  TRUE
                 if(input$analyse_quiet == "yes")
                    quietvalue <- FALSE
                print(quietvalue) 
                 res <<- AM(trait=input$nmst , fformula=fform , availmemGb = input$memsize , 
                            quiet = quietvalue,
                            ncpu = input$analyse_cpu, maxit = input$analyse_maxits , pheno = pheno, geno=geno, map=map) 

              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "AM", html = m$message, add = TRUE)
       })

  })  ## end observeEvent


 ##--------------------------
 ## Print findings .... 
 ##--------------------------

 ## form data frame of results 
 observeEvent(input$analyse_go, {
 dfres <- NULL
 if(length(res$Mrk)>1){
    dfres <- data.frame(snps=res$Mrk[-1], chrm=res$Chr[-1], position=res$Pos[-1])

 }

 if(is.null(dfres)){
    ## no associations found
    output$findings <- renderText("No significant associations between snp and trait were found.")
  }
  if(!is.null(dfres)){
    output$findings <- renderTable(dfres)
  }

     observeEvent(input$pvalue_go, {

      withCallingHandlers({
                 shinyjs::html("summary", "")
                  sumres <- SummaryAM(AMobj=res, pheno=pheno, geno=geno, map=map)
                  output$pvalue <- renderTable(sumres[["pvalue"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$size <- renderTable(sumres[["size"]], digits=-1, hover=TRUE, bordered=TRUE)
                  output$R <- renderTable(sumres[["R"]],  hover=TRUE, bordered=TRUE)


              },  ## end withCallingHandlers
              message = function(m) {
                 shinyjs::html(id = "summary", html = m$message, add = TRUE)
       })





  })  ## end observeEvent




}) ## end observeEvent


##----------------------
## Help - Docs
##---------------------



## FAQ html 
#   output$faq <- renderUI({
#        HTML(markdown::markdownToHTML(knit('faq.rmd', 
#              quiet = TRUE),options=c('toc'), fragment.only=TRUE))
#    })

# observeEvent(input$quickstart, {
#      RShowDoc("QuickStart", package="Eagle")
#})

  ##---------------------------------
  ## Help files   - addPopover
  ##--------------------------------

 
  addPopover(session, "dummy1", "Details", content = HTML("
Eagle can handle two types of marker genotype file; a space separated plain text file and PLINK
ped file. We assume the marker loci are snps. 
Missing marker genotypes are allowed but the 
proportion of missing genotypes is assumed to be low. 
<br><br>
The marker genotype file should not contain column names. 
We also assume that each row of data in the file corresponds to data on a different 
individual. The ordering of the rows, by individual, must be the same for the marker genotype file and phenotypic file.<br><br> 


If the file is a plain text file, then character or numeric genotypes can be used to 
denote the snp genotypes. However, Eagle needs 
to map these genotypes to its internal snp genotypes. We do this by asking the user to assign their snp genotype codes to our AA, 
AB and BB codes. If data are collected on inbred individuals, only AA and BB need be specified. 
<br> <br>

To load the marker genotype data into Eagle, follow the four steps.  Upon cliking Upload, Eagle checks the genotype file for errors, 
and recodes the genoytpes for later analysis. If the marker genotype is large (many Gbytes), this step can take several minutes. <br><br>
Output from reading in the marker genotype file will appear in the right hand-side panel. 
                                                           "), trigger = 'hover')

  
  addPopover(session, "dummy2", "Details", content = HTML("  
  Eagle assumes the phenotypic file is either a space separated or comma separated file. The rows correspond to data 
  collected on the individuals. The first row of the file can contain column headings or not.  
  The number of rows of data in the phenotypic file must equal the number of rows in the marker genotype file otherwise an error occurs.  
  Also, Eagle assumes the phenotypic data is row ordered by individual in the same way as the marker genotype data. 
 <br> <br>
Data on multiple traits and fixed effects that may or may not be used in the analysis can be included in this file.  <br> <br>
Missing values are allowed.  <br> <br>
Output from reading in the phenotypic file will appear in the right hand-side panel. 
  "), trigger = 'hover')
  
   addPopover(session, "dummy3", "Details", content = HTML("
    Eagle does not
     require a known marker map in order to analyse the data.  
     If a map file is read into Eagle, then the
     marker names are used when results are reported in  'Findings'. If a
     map file is not supplied, generic names M1, M2, ..., are
     given to the marker loci. 
      <br> <br>
     The map file can have three or four columns. If the
     map file has three columns, then it is assmed that the three
     columns are the marker locus names, the chromosome number, and the
     map position (in any units). If the map file has four columns as
     with a PLINK map file, then the columns are assumed to be the
     marker locus names, the chromosome number, the map position in
     centimorgans, and the map position in base pairs.
      <br> <br>
     Missing values are not allowed.
      <br> <br>
    The order of the marker loci in this file is assumed to be in the
     same order as the loci (or columns) in the marker data file.
  "), trigger = "hover") 



addPopover(session, "dummy4", "Details", content = HTML(paste("
    Here,  multiple_locus association mapping is performed. The analysis
     simultaneously accounts for  familial
     relatedness and nuisance fixed effects while detecting 
     multiple marker-trait associations. Unlike other association mapping methods, 
     there are no regularization parameters to be tuned, nor significance thresholds to be set. 
     <br><br>
     Output from performing the analysis is printed to the right hand panel. A table of results is printed in 'Findings'. 
     <br><br>", tags$span(style="color:red",
     "Once an analysis has been completed, a new analysis can be performed  
     by selecting a new trait in 'Step1' or different fixed effects in 'Step2' and clicking the 'Perform genome-wide analysis' button.", sep=""))


), trigger = "hover")



addPopover(session, "dummy5", "Details", content = HTML(paste("
By default, the 'best' set of snp in strongest association with the trait are reported. 
These results are given in table form. 
<br><br>
By clicking the 'Additional Summary' button on the left, two additional tables of results are shown; a table on the size and significance of the snp 
in the model and a table for the amount of phenotypic variance explained as they are added one at a time to the model.
<br><br>
There is additional computation needed to produce these extra tables. It may take a few minutes before these tables appear. 
     ", sep="")


), trigger = "hover")





 
session$onSessionEnded(stopApp)
  
})