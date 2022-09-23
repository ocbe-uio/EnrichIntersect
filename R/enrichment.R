#' @title Plot enrichment map
#' @description
#' Plot enrichment map through a vector (matrix) of scores and a self-defined set that summarizes a few groups of the names (rownames) of the vector (matrix)
#' 
#' @name enrichment
#' 
#' @import ggplot2
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom stats p.adjust 
#' 
#' @param x a vector (matrix) of scores to be enriched
#' @param custom.set a self-defined set that summarizes a few groups of the names (rownames) of \code{x}
#' @param alpha exponent weight of the score of ordered features. Default is \code{0} for calculating enrichment score via classic Kolmogorov-Smirnov statistic
#' @param normalize logic value to determine if normalizing enrichment scores, accounting for custom set size. Default is \code{TRUE}
#' @param permute.n number of custom-set permutations for significance testing. Default is 100
#' @param padj.method correction method, one of \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, \code{"none"}. Default is \code{"none"}
#' @param pvalue.cutoff a cutoff for both unadjusted and adjusted p-value to mark significantly enriched classes. Default is 0.05
#' @param angle angle of rotating x-axis labels. Default is 45
#' @param ... other arguments
#' 
#' @return Return a list including a matrix of (normalized) enrichment score, a matrix of corresponding p-value and ggplot object:
#' \itemize{
#' \item S - a matrix of calculated enrichment scores.
#' \item pvalue - a matrix of p-values using permuation test for the calculated enrichment scores.
#' \item g - a ggplot object for visualising the results of an enrichment analysis.
#' }
#' 
#' @references Reimand J, Isserlin R, Voisin V, et al (2019). \emph{Pathway enrichment analysis and visualization of omics data using g:profiler, gsea, cytoscape and enrichmentmap}. Nature protocols, 14:482–517.
#' 
#' @examples
#' # Data set 'cancers_drug_groups' is a list including a score dataframe with 147 drugs as rows 
#' # and 19 cancer types as columns, and a dataframe with 9 self-defined drug groups (1st column)
#' # of the 147 drugs (2nd column).
#' data(cancers_drug_groups, package = "EnrichIntersect")
#' 
#' x <- cancers_drug_groups$score
#' custom.set <- cancers_drug_groups$custom.set
#' set.seed(123)
#' enrich <- enrichment(x, custom.set, permute.n=5)
#' 
#' @export
enrichment <- function(x, custom.set, alpha=0, normalize=TRUE, permute.n=100, padj.method="none", pvalue.cutoff=0.05, angle=45, ...){
  
  if(is.matrix(x) | is.data.frame(x)){
    if(any( colSums(is.na(x))==ncol(x) ) & ncol(x)>1)
      stop("The argument 'x' matrix has some columns with all missing values!")
  }else{
    x <- as.matrix(x)
  }
  if(is.null(colnames(x)))
    colnames(x) <- "X"
  
  if(is.matrix(custom.set) | is.data.frame(custom.set)){
    if(dim(custom.set)[2] != 2)
      stop("The argument 'custom.set' has to have two columns!")
  }else{
    stop("The argument 'custom.set' has to be a matrix or dataframe!")
  }
  
  if( length(unique(custom.set[[1]])) < length(unique(custom.set[[1]])) )
    stop("The argument 'custom.set' should have more sets/groups (2nd component of 'custom.set') than unique symbols (1st component of 'custom.set')!")
    
  features <- intersect( rownames(x), custom.set[[1]] )
  groups <- unique(custom.set$group[custom.set[[1]] %in% features])
  n_groups <- length(groups)
  
  # initialize some enrichment scores and pvalues
  S <- matrix(nrow=ncol(x), ncol=n_groups)
  rownames(S) <- colnames(x); colnames(S) <- groups
  pvalue <- S_norm <- S
  
  # Initializes the progress bar
  pb <- txtProgressBar(min=0, max=ncol(x), style=3, width=50,  char="=")
  
  for(i in 1:ncol(x)){
    ## define costom sets 'myList'
    cutOff <- -Inf # this parameter remains for unordered features for the future
    myList <- x[,i]
    names(myList) <- features
    myList <- sort(myList, decreasing = TRUE)
    myList <- myList[myList>cutOff]
    n <- length(myList)
    
    for(k in 1:n_groups){ # group index
      
      idx <- rep(NA,n)
      for(j in 1:n){
        idx[j] <- names(myList[j]) %in% custom.set[custom.set$group==groups[k],1]
      }
      if(sum(idx)){
        F1 <- cumsum(abs(myList)^alpha * as.numeric(idx)) / sum(abs(myList)^alpha * as.numeric(idx))
        F2 <- cumsum(as.numeric(!idx)) / (n-sum(custom.set$group==groups[k] & features %in% names(myList)))
        
        S[i,k] <- max(F1 - F2, na.rm=T) ## NO negative, indicating no enriched drug groups for poor prediction of drug responses
      }else{
        S[i,k] <- NA
      }
      
      # permutation test
      permute_S <- rep(NA, permute.n)
      for(permute_i in 1:permute.n){
        permutationIdx <- sample(1:n,n) # n removes missing spearman's x
        myList_permute <- myList
        names(myList_permute) <- names(myList)[permutationIdx]
        for(j in 1:n){
          idx[j] <- names(myList_permute[j]) %in% custom.set[custom.set$group==groups[k],1]
        }
        if(sum(idx)){
          F1 <- cumsum(abs(myList_permute)^alpha * as.numeric(idx)) / sum(abs(myList_permute)^alpha * as.numeric(idx))
          F2 <- cumsum(as.numeric(!idx)) / (n-sum(custom.set$group==groups[k] & features %in% names(myList)))
          
          permute_S[permute_i] <- max(F1 - F2, na.rm=T)
        }
        #if(permute_i %% 100 == 0) cat("task_i=", i, "; group_k=", k, "; permute_i=",permute_i,"\n", sep="")
      }
      ## i) p-value for enrichment score
      pvalue[i,k] <- sum(permute_S >= S[i,k],na.rm=TRUE)/(permute.n-sum(is.na(permute_S)))
      
      
      ## ii) normalized enrichment score
      if( normalize ){
        S_norm <- S
        if(!is.na(S[i,k])){
          if(S[i,k]>=0){
            S_norm[i,k] <- S[i,k]/mean(permute_S[permute_S>=0], na.rm=TRUE)
          }else{
            S_norm[i,k] <- S[i,k]/mean(permute_S[permute_S<0], na.rm=TRUE)
          }
        }
        S <- S_norm
        permute_S_norm <- rep(NA, permute.n)
        for(permute_i in 1:permute.n){
          if(!is.na(permute_S[permute_i])){
            if(permute_S[permute_i]>=0){
              permute_S_norm[permute_i] <- permute_S[permute_i]/mean(permute_S[permute_S>=0], na.rm=TRUE)
            }else{
              permute_S_norm[permute_i] <- permute_S[permute_i]/mean(permute_S[permute_S<0], na.rm=TRUE)
            }
          }
        }
        pvalue[i,k] <- sum(permute_S_norm >= S_norm[i,k],na.rm=TRUE)/(permute.n-sum(is.na(permute_S_norm)))
      }
      
    }
    # Sets the progress bar to the current state
    setTxtProgressBar(pb, i)
  }
  close(pb) # Close the connection
  
  pvalue[is.na(S)] <- NA
  pvalue <- t(pvalue)
  if(! padj.method == "none"){
    pvalue <- p.adjust(pvalue, method=padj.method)
  }
  S <- t(S)
  
  # define a dataframe
  dat  <- data.frame(x = factor(rep(colnames(S), each=nrow(S))), 
                     y = rep(rownames(S), ncol(S)), 
                     ks = as.vector(S), 
                     pvalue = as.vector(pvalue))
  dat$y <- factor(dat$y, levels = levels(factor(dat$y))[c(length(unique(dat$y)):1)])
  ks.min <- min(dat$ks,na.rm=T)
  ks.max <- max(dat$ks,na.rm=T)
  dat$ks[dat$ks < 0 & !is.na(dat$ks)] <- ks.min - 0.001
  dat$ks <- dat$ks - (ks.min - 0.001)
  dat$ks[which.max(dat$ks)] <- round(sort(dat$ks,decreasing = T)[2]+0.5)
  dat$border <- rep("red", nrow(dat)); dat$border[dat$pvalue >= pvalue.cutoff] <- "gray"
  
  #globalVariables(c("y", "border", "ks"))
  y <- border <- ks <- NULL
  
  if(normalize){
    ES.name <- "Normalized\nEnrichment\nScore"
  }else{
    ES.name <- "Enrichment\nScore"
  }
  if( any(dat$pvalue<pvalue.cutoff) ){
    
    g <- ggplot(data = dat) +
      geom_point(aes(x=x, y=y, color=border, fill=pvalue, size=ks), shape=21) +
      scale_fill_gradientn(name="p-value", na.value = "black",colours=c("blue", "white"), 
                           limits=c(0,1), guide=guide_colorbar(barheight=3, barwidth=1)) +
      geom_point(aes(x=x, y=y, size=ks),color=dat$border, shape=21) +
      guides(colour = guide_legend(override.aes = list(size=5))) +
      scale_size(name = ES.name, range = c(1, 5), breaks = c(0,1,2,3), guide = guide_legend(keyheight=.8)) +
      theme(axis.text.x = element_text(size = 8, angle = angle, vjust = 1, hjust = 1)) + xlab("") + ylab("") 
    
    if(pvalue.cutoff==0.05){
      g <- g + scale_color_manual(name=NULL, values=c( "p<0.05" = "red"))
    }else{
      if(pvalue.cutoff==0.1){
        g <- g + scale_color_manual(name=NULL, values=c( "p<0.1" = "red"))
      }else{
        g <- g +  scale_color_manual(name=NULL, values=c("gray","red"), labels=paste("p",c(">=","<"),pvalue.cutoff,sep="") )
      }
    }
  }else{
    g <- ggplot(data = dat) +
      geom_point(aes(x=x, y=y, fill=pvalue, size=ks), shape=21) +
      scale_fill_gradientn(name="p-value", na.value = "black",colours=c("blue", "white"), 
                           limits=c(round(min(dat$pvalue),2),1), guide=guide_colorbar(barheight=3, barwidth=1)) +
      scale_size(name = ES.name, range = c(1, 5), breaks = c(0,1,2,3), guide = guide_legend(keyheight=.8)) +
      theme(axis.text.x = element_text(size = 8, angle = 45, vjust = 1, hjust = 1)) + xlab("") + ylab("") 
  }
  print(g)
  
  return(list(S=S, pvalue=pvalue, g=g))
}
