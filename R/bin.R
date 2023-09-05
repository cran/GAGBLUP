#' @title Binning Genotypes for Dimensional Reduction
#' @description Binning the original genotypes into bins for dimensional reduction 
#' under the principle of linkage disequilibrium.
#' @param genotype a matrix for genotypes in numeric format, coded as 1, 0 and -1, with individuals in rows and markers in cols.
#' @param binvar a hyper-parameter between 0 and 1, the closer to 0, the fewer bins yields.
#' Users can choose binvar based on the required number of bins, default is 0.15.
#' @return A list with following information is returned:
#'     $bins_genotypes  binned genotypes
#'     $bins_range  start and stop of each bin 
#' @examples
#' \donttest{
#' ## load example data from GAGBLUP package
#' data(genotype)
#' ## binning genotypes
#' bins <- bin(genotype,0.2)
#'  }
#' @export

bin <- function(genotype=genotype,binvar=0.15){
  binning <- function (SNP,var_thr) {
    SNP=data.matrix(SNP)
    beg<- numeric()
    end<-numeric()
    block_mean <- numeric()
    i=1
    while(i < nrow(SNP))
    {
      var=1
      j=0
      while(var >= var_thr & i+j<nrow(SNP)){
        j=j+1
        average<- apply(SNP[i:(i+j),],2,mean)
        var <- var(average)
      }
      beg<-rbind(beg,i)
      end<-rbind(end,i+j-1)
      ifelse(j>1,average_gene<-apply(SNP[i:(i+j-1),],2,mean),average_gene<-SNP[i,])
      block_mean<- rbind(block_mean,average_gene)
      i=i+j
    } 
    range<-cbind(beg,end)
    return(list(range=range,block_mean=block_mean))
  }
  
  genotype<-scale(genotype)
  genotype <- na.omit(t(genotype))
  binvar = binvar
  
  result<-binning(SNP=genotype,var_thr=binvar)
  bins_genotypes <- t(round(result$block_mean,1))
  bins_range <- result$range
  colnames(bins_range) <- c('start','stop')
  return(list(bins_genotypes = bins_genotypes,bins_range = bins_range))
}





