#' @title Genetic algorithm assisted genomic best liner unbiased prediction for genomic selection
#' @description Performs genomic selection with genetic algorithm assisted genomic best liner unbiased prediction
#' @param genotype a matrix for genotypes in numeric format, with individuals in rows and markers in cols.
#' @param phenotype a vector of phenotype, missing (NA) values are not allowed.
#' @param fit_fun the fitness function. There are four options: fitness = "AIC"/"BIC"/"FIT"/"HAT", default is "HAT"
#' @param maxiter max number of iterations for GAGBLUP, default is 2000 
#' @param nfold the number of folds. Default is 10.
#' @param nTimes the number of independent replicates for the cross-validation. Default is 1.
#' @param seed the random number. Default is 123.
#' @param n_core the number of CPU to be used, default is 1.
#' @return A list with following information is returned:
#'         $R2  the squared pearson correlation coefficient between the true value and the predicted value,
#'         $predicted_value  the predicted value and the true value of the phenotype,
#'         $marker_selection  a vector of the selected markers, with 1 indicates selected, 0 indicates not selected.
#' @examples
#' ## load example data from GAGBLUP package
#' data(phenotype)
#' data(bins)
#' phenotype <- phenotype[1:200,3]
#' result <- gagblup(bins[1:200,],phenotype,fit_fun='HAT',maxiter=1,nfold=2,nTimes=1,seed=123,n_core=1)
#' @export

gagblup <- function(genotype,phenotype,fit_fun='HAT',maxiter=2000,nfold=10,
                    nTimes=1,seed=123,n_core=1){
  z0 <- genotype
  y0 <- phenotype
  
  set.seed(seed)
  n <- length(y0)
  foldid <- sample(rep(1:nfold,ceiling(n/nfold))[1:n])
  mixed<-function(x,y,kk,cov="qq"){
    if(cov!="qq"){
      kk<-eigen(kk,symmetric=T)
    }
    uu<-kk$vectors
    delta<-kk$values
    x<-t(uu)%*%x
    y<-t(uu)%*%y
    r<-ncol(x)
    func<-function(lambda){
      h<-1/(lambda*delta+1)
      x1<-x*sqrt(h)
      y1<-y*sqrt(h)
      xx<-crossprod(x1)
      xy<-crossprod(x1,y1)
      yy<-crossprod(y1)
      yx<-t(xy)
      dd1<-sum(log(1/h))
      yPy<-yy-yx%*%solve(xx)%*%xy
      dd2<-log(det(xx))  
      loglike<--0.5*dd1-0.5*dd2-0.5*(n-r)*log(yPy)
      return(-loglike)
    }  
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      x1<-x*sqrt(h)
      y1<-y*sqrt(h)
      xx<-crossprod(x1)
      xy<-crossprod(x1,y1)
      yy<-crossprod(y1)
      yx<-t(xy)
      cc<-solve(xx)
      yPy<-yy-yx%*%cc%*%xy
      beta<-drop(cc%*%xy)  
      s2<-drop(yPy/(n-r))
      v<-cc*s2
      stderr<-sqrt(diag(v))	
      result<-c(beta,stderr,s2)
      return(result)
    }
    par0<-1
    fit<-optim(par=par0,fn=func,method="L-BFGS-B",lower=1e-8,upper=1e8)
    lambda<-fit$par
    conv<-fit$convergence
    fn1<-fit$value
    fn0<-func(1e-8)
    lrt<-2*(fn0-fn1)
    blue<-fixed(lambda)
    beta<-blue[1:r]
    stderr<-blue[(r+1):(2*r)]
    ve<-blue[2*r+1]
    lod<-lrt/4.61
    if(lrt<1e-8){
      p<-1
    } else {
      p<-0.5*(1-pchisq(lrt,1))
    }
    va<-lambda*ve
    h2<-va/(va+ve)
    par<-data.frame(beta,stderr,va,ve,lambda,h2,conv,fn1,fn0,lrt,lod,p)
    return(par)
  }
  
  cl <- makeCluster(n_core)
  registerDoParallel(cl)
  message('Predicting by GAGBLUP...')
  rept <- NULL
  R2_all <- foreach(rept = 1:nTimes, .multicombine = TRUE, .packages=c('GA')) %dopar% {
    yy<-NULL
    for(ii in 1:nfold){
      indx1<-which(foldid==ii)
      indx2<-which(foldid!=ii)
      z1<-z0[indx2,]        
      y<-y0[indx2,drop=FALSE]
      
      n<-nrow(z1)
      m<-ncol(z1)
      
      fitness<-function(parm){
        indx<-which(parm==1)
        z<-z1[,indx]
        kk<-z%*%t(z)
        qq<-eigen(kk,symmetric=T)  
        uu<-qq$vectors
        delta<-qq$values
        x<-matrix(1,n,1)
        fit<-mixed(x,y,kk=qq,cov="qq")
        lambda<-fit$lambda
        likelihood<- -fit$fn1
        beta<-as.matrix(fit$beta)
        delta<-abs(qq$values)
        vector<-qq$vectors
        r<- y-x%*%beta
        d<-lambda*delta/(delta*lambda+1)
        mat<-t(vector)*sqrt(d)
        H<-crossprod(mat)
        eHat<-r-H%*%r
        SSE<-sum(eHat^2)
        SST<-var(r)*(n-1)
        e<-eHat/(1-diag(H))
        PRESS<-sum(e^2)
        df<-length(indx)
        
        BIC<- -(-2*likelihood+log(n)*df)
        AIC<- -(-2*likelihood+2*df)
        R2.HAT<-1-PRESS/SST
        R2.FIT<-1-SSE/SST
        
        if(fit_fun == 'AIC'){
          return(AIC)
        }else if (fit_fun == 'BIC'){
          return(BIC)
        }else if(fit_fun == 'FIT'){
          return(R2.FIT)
        }else{
          return(R2.HAT)
        }
      }
      
      
      parm <- rbinom(m,1,0.01)
      pred <- ga("binary", fitness = fitness,pmutation = 0.1,nBits = m,names = colnames(z1),
                 popSize=100, maxiter=maxiter)
      
      selection<-pred@solution[1,]
      marker<-which(selection==1)
      
      indx<-marker
      z<-z1[,indx]
      kk<-z%*%t(z)
      qq<-eigen(kk,symmetric=T)
      uu<-qq$vectors
      delta<-qq$values
      x<-matrix(1,n,1)
      fit<-mixed(x,y,kk=qq,cov="qq")
      lambda<-fit$lambda
      beta<-fit$beta
      
      z2<-z0[indx1,indx]
      zz<-z2%*%t(z)
      yhat<-beta+lambda*zz%*%solve(kk*lambda+diag(n))%*%(y-x%*%beta)
      y2<-as.matrix(y0[indx1,drop=FALSE])
      r2<-cor(y2,yhat)^2
      print(c(ii,r2))
      yy<-rbind(yy,cbind(yhat,y2))
      colnames(yy) <- c('predicted','observed')
    }
    R2<-cor(yy[,1],yy[,2])^2
    list(R2=R2,predicted_value=yy,marker_selection=selection)
  }
  stopImplicitCluster()
  stopCluster(cl)
  message('Prediction by GAGBLUP ended')
  return(R2_all)
}
