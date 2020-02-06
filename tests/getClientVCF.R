library(resourcer)



# client: from url to GDS object

GDSFileResourceClient <- R6::R6Class(
  
  "GDSFileResourceClient",
  
  inherit = FileResourceClient,
  
  public = list(
    
    initialize = function(resource) {
      
      super$initialize(resource)
      
    },
    
    getConnection = function() {
      
      conn <- super$getConnection()
      
      if (is.null(conn)) {
        
        path <- super$downloadFile()
        
        resource <- super$getResource()
        
        resource
        
        format <- resource$format
        
        if ("VCF2GDS" == toupper(format)) {
          
          private$loadSNPRelate() # also loads gdsfmt
          
          private$.gds.file.tmp <- tempfile(fileext = ".gds")
          
          url <- super$parseURL()
          
          method <- "biallelic.only"
          
          snpfirstdim <- FALSE
          
          if (!is.null(url$query)) {
            
            if (!is.null(url$query$method)) {
              
              method <- url$query$method
              
            }
            
            if (!is.null(url$query$snpfirstdim)) {
              
              snpfirstdim <- as.logical(url$query$snpfirstdim)
              
            }
            
          }
          
          snpgdsVCF2GDS(vcf.fn = path, out.fn = private$.gds.file.tmp, method = method, snpfirstdim = snpfirstdim)
          
          path <- private$.gds.file.tmp
          
          snpgdsSummary(path)
          
        } else {
          
          loadGDSFmt()
          
        }

#     Yannick's code        
#        conn <- openfn.gds(path)
        conn <- snpgdsOpen(path)
        super$setConnection(conn)
        
      }
      
      conn
      
    },
    
    getValue = function() {
      
      self$getConnection()
      
    },
    
    close = function() {
      
      super$close()
      
      conn <- super$getConnection()
      
      if (!is.null(conn)) {
        
        closefn.gds(conn)
        
        if (!is.null(private$.gds.file.tmp) && file.exists(private$.gds.file.tmp)) {
          
          # house keeping
          
          ignore <- tryCatch(file.remove(private$.gds.file.tmp), error = function(e) {})
          
        }
        
      }
      
    }
    
  ),
  
  private = list(
    
    .gds.file.tmp = NULL,
    
    loadSNPRelate = function() {
      
      if (!require("SNPRelate")) {
        
        install.packages("SNPRelate", repos="http://R-Forge.R-project.org", dependencies = TRUE)
        
      }
      
    },
    
    loadGDSFmt = function() {
      
      if (!require("gdsfmt")) {
        
        if (!require("BiocManager")) {
          
          install.packages("BiocManager", repos = "https://cloud.r-project.org", dependencies = TRUE)
          
        }
        
        BiocManager::install("gdsfmt")
        
      }
      
    }
    
  )
  
)



# GDS resolver: from url to GDS client

GDSFileResourceResolver <- R6::R6Class(
  
  "GDSFileResourceResolver",
  
  inherit = ResourceResolver,
  
  public = list(
    
    isFor = function(x) {
      
      if (super$isFor(x)) {
        
        !is.null(findFileResourceGetter(x)) && toupper(x$format) %in% c("VCF2GDS", "GDS")
        
      } else {
        
        FALSE
        
      }
      
    },
    
    newClient = function(x) {
      
      if (self$isFor(x)) {
        
        GDSFileResourceClient$new(x)
        
      } else {
        
        NULL
        
      }
      
    }
    
  )
  
)
