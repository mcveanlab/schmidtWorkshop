#' calculate BF
#'
#' @export
#'

calc.lBF <- function(
    pars = NULL,
    data.sub = NULL,
    w0 = 2,
    log10 = TRUE
) {

    if( is.null(pars) | is.null(data.sub) ) {
        stop("Missing input data.\n")
    }
    
    if (ncol(data.sub)==2) w0 <- 1;

    llk.full <- calc.llk.tree(pars, data.sub);

    q00 <- pars$p.stay + pars$p.switch * (1-pars$pi1);
    lq00 <- log(q00);
    n.trans <- nrow(data.sub)-1;
    l.p.null <- log(1-pars$pi1)+n.trans*lq00;
    l.lk.null <- sum(data.sub[,w0]);
    tmp <- c(llk.full, l.p.null+l.lk.null);
    mx <- max(tmp);
    tmp <- exp(tmp-mx);
    lBF <- mx+log(tmp[1]-tmp[2])-l.lk.null-log(1-exp(l.p.null));
    if (log10) lBF <- lBF/log(10);
    
    return(lBF);
}

#' Function to get marginal posterior on -/0/+ profile
#'
#' @param arg input 
#' @export
#' @examples
#' marginal.posterior.profile()
#'

marginal.posterior.profile <- function(
    pars = NULL,
    data.sub = NULL
) {

    ## Get forward and G matrices
    tmp <- calc.llk.tree(pars, data.sub, returnForwardMatrix=T, returnGMatrix=T);
    f <- tmp$f;
    g <- tmp$g;
    
    ## Build parents list and reverse F traversal order
    parents <- rep(0, pars$n.phenos);
    for (i in pars$t.path) parents[pars$ontology[[i]]] <- i;
    ord <- c(rev(pars$t.path), pars$terminals);

    ## Construct backward matrix
    b <- array(0, dim(f));
    b[ord[1],] <- log(pars$stat.dist);
    for (i in ord[-1]) {
        r.i <- b[parents[i],]+f[parents[i],]-g[i,];
        mx <- max(r.i);
        tmp <- mx+log(sum(exp(r.i-mx)));
        tmp <- tmp+pars$lp.switch+log(pars$stat.dist);
        tmp2 <- -pars$theta.tree+b[parents[i],]+f[parents[i],]-g[i,];
        tmp3 <- cbind(tmp2, tmp);
        mx <- apply(tmp3, 1, max);
        b[i,] <- log(rowSums(exp(tmp3-mx)))+mx;
    }
    
    ## Posteriors
    tmp <- b+f;
    mx <- apply(tmp, 1, max);
    pp <- exp(tmp-mx);
    pp <- pp/rowSums(pp);
    return(pp);
}

#' Function to calculate integrated likelihood for set of variants up tree
#'
#' @export
#'

calc.llk.tree<-function(
    pars,
    data.sub,
    returnForwardMatrix = FALSE,
    returnGMatrix = FALSE
) {

    ## Get integrated likelihood at nodes - will be overwritten for internal nodes
    mx <- apply(data.sub, 1, max);
    d <- exp(data.sub-mx);
    llk.integrated <- log(d %*% pars$stat.dist)+mx;
    if (returnGMatrix) g <- array(0, dim(d));
    
    for (i in pars$t.path) {
        emiss.node<-data.sub[i,];			# Emissions at node
        data.sub[i,]<-0;
        for (j in pars$ontology[[i]]) {
            tmp1<-cbind(data.sub[j,]-pars$theta.tree, llk.integrated[j]+pars$lp.switch);
            mx1<-apply(tmp1, 1, max);
            tmp2<-mx1+log(rowSums(exp(tmp1-mx1)));
            data.sub[i,]<-data.sub[i,]+tmp2;
            if (returnGMatrix) g[j,]<-tmp2;
        }
        data.sub[i,]<-data.sub[i,]+emiss.node;
        mx<-max(data.sub[i,]);
        llk.integrated[i]<-mx+log(sum(exp(data.sub[i,]-mx)*pars$stat.dist));
    }
    if (returnForwardMatrix) {
        if (!returnGMatrix) {
            return(data.sub);
        } else {
            return(list(f=data.sub, g=g));
        }
    } else {
        return(llk.integrated[i]);
    }
}
