#Metric and non-metric MDPREF model
##by CW Kwan
#Mar 5, 2019
#S = ULV'
#subjects are columns, objects are rows
#


metric.mdpref<-function(x,ndim) {
  mysvd<-svd(x)
  n<-nrow(x)
  score<-mysvd$u[,1:ndim]*sqrt(n-1)
  corr<-mysvd$v[,1:ndim] %*% diag(mysvd$d[1:ndim])/sqrt(n-1)
  fitted<-score %*% t(corr)
  rownames(score)<-rownames(x)
  rownames(corr)<-colnames(x)
  colnames(score)<-colnames(corr)<-paste("dim",1:ndim,sep="")
  rownames(fitted)<-rownames(x)
  colnames(fitted)<-colnames(x)
  decomp<-list(score=score,corr=corr,fitted=fitted,d=mysvd$d)
}


mdpref<-function(x,ndim,monotone,tor,maxit) {
  p<-ncol(x)
  n<-nrow(x)
  x0<-scale(x)
  fit<-metric.mdpref(x0,ndim)
  fit0<-rbind(fit$score,fit$corr)
  if (monotone) {
    m.x<-x0
    itn<-0
    repeat {
      itn<-itn+1
      for (i in 1:p) {
        m.x[,i]<-opscale(x0[,i],fit$fitted[,i],level=2)$os
      }
      fit<-metric.mdpref(m.x,ndim)
      fit1<-rbind(fit$score,fit$corr)
      err<-sum((fit0-fit1)^2)
      if (err<tor | itn>maxit) {break}
      fit0<-fit1
    }
    fit$tpref<-m.x
    if (itn>=maxit) {print('convergence not met')}
  }
  return(fit)
}


#' biplot for mdperf
#'
#' @param x mdpref object
#' @param ... biplot() parameters
#' @export
plot.mdpref<-function(x,...) {
  ndim<-ncol(x$score)
  v<-x$d^2
  v<-round(v/sum(v)*100,2)
  for (i in 1:(ndim-1)) {
    for (j in (i+1):(ndim)) {
      xlabs<-paste("Dim",i," (",v[i],"%)")
      ylabs<-paste("Dim",j," (",v[j],"%)")
      biplot(x$score[,c(i,j)],x$corr[,c(i,j)],
               xlab=xlabs,ylab=ylabs,...)
    }
  }
}

#' summary for mdpref
#'
#' @param object mdpref object
#' @param ... summary parameters
#' @export
summary.mdpref<-function(object,...) {
  cat("Object coordinates\n")
  print(round(object$score,4))
  cat("\nSubject vectors\n")
  print(round(object$corr,4))
  ndim<-ncol(object$corr)
  cat("\nVariation explained\n")
  v<-object$d^2
  myv<-rbind(object$d,v,v/sum(v)*100,cumsum(v)/sum(v)*100)
  rownames(myv)<-c("Sing.val","Val.exp","%","Cum%")
  myv<-myv[,1:ndim]
  colnames(myv)<-colnames(object$corr)
  print(round(myv[,1:ndim],4))
}

#' MDPREF model
#'
#' Metric and non-metric MDPREF model
#'
#' @author Chi-wai Kwan
#' @param pref Preference data matrix. Rows as objects and columns as subjects. Small value for less preferable and large value for more preferable.
#' @param ndim Number of dimensions.
#' @param monotone TRUE for Kruskal monotonic transformation. FALSE for metric scale.
#' @param tor tolerance for monotonic transformation.
#' @param maxit maximum number of interactions for monotonic transformation.
#'
#' @import optiscale lattice
#'
#' @details
#' Metric scale:\cr
#' Singular decomposition is applied to data matrix S.\cr
#' S = ULV', X = U*sqrt(n-1),
#' Y = V*sqrt(L)/sqrt(n-1), n = no. of objects.\cr \cr
#' Non-metric scale:\cr
#' Kruskal monotonic transformation is applied to each subject vector by opscale(level=2,...) in optiscale package.\cr
#' Metric MDPREF is applied to the transformed data.
#'
#' @return An object of class "mdpref".
#' \item{score }{Object coordinates}
#' \item{corr }{Subject vectors}
#' \item{d }{Singular values}
#' \item{fitted }{fitted preference matrix}
#' \item{tpref }{transformed preference matrix}
#'
#' @references
#' Carroll, J. D. (1972). “Individual Differences and Multidimensional Scaling.” In Multidimensional Scaling: Theory and Applications in the Behavioral Sciences, vol. 1, edited by R. N. Shepard, A. K. Romney, and S. B. Nerlove, 105–155. New York: Seminar Press.\cr\cr
#' Kruskal, Joseph B. (1964) “Nonmetric Multidimensional Scaling: A Numerical Method.” Psychometrika 29: 115-129
#'
#' @examples
#' library(lattice)
#' library(optiscale)
#' library(cmdpref)
#' #rank order of preferences of 5 objects for 4 subjects
#' mydata
#'
#' fit<-cmdpref(pref=mydata,monotone=TRUE,maxit=100)
#' summary(fit)
#'
#' #biplot() arguments are applicable
#' plot(fit,xlim=c(-2,2),cex=1)
#' @export
cmdpref<-function(pref,ndim=2,monotone=FALSE,tor=1e-8,maxit=50) {
  value<-mdpref(x=pref,ndim=ndim,monotone=monotone,tor=tor,maxit=maxit)
  attr(value, "class") <- "mdpref"
  value
}
