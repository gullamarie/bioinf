#!/usr/bin/env Rscript

#set libPaths in R env: 
.libPaths( c( "/R/x86_64-redhat-linux-gnu-library/3.6" , .libPaths() ))

#load libraries
library(rslurm)
library(facets)

files <- as.list(list.files(path="/FACETS-INFO", pattern="*.csv", full.names=TRUE, recursive=FALSE))

facets_fun <- function(x){

        set.seed(123)
        print(x)

        rcmat <- readSnpMatrix(x)
        xx <- preProcSample(rcmat, gbuild="hg38")

        oo <- procSample(xx,cval=150)

        fit <- emcncf(oo)

        outfilename <- paste(tumor, "fit.cncf.csv", sep=".")
        print(outfilename)
        write.table(fit$cncf, outfilename, sep = ";", dec = ".", row.names = TRUE, col.names = TRUE)
}

lapply_fun <- function(files){
	lapply(files, facets_fun)
}


sjob <- slurm_call(lapply_fun, files, jobname='test_facets', submit=FALSE)
cleanup_files(sjob)


#lapply(files, function(x) {
#	sjob <- slurm_call(facets_fun, x, jobname='test_facets', submit=FALSE)
#	cleanup_files(sjob)
#})

#sjob <- slurm_call(facets_fun, files, jobname='test_facets', submit=FALSE)
#cleanup_files(sjob)


