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

grs.marginal.posterior <- function(
    tree,
    prior,
    p.stay,
    p.switch,
    llk.data
) {
    null.id <- which.min(abs(prior$b.grid))

    ## Get forward and G matrices
    tmp <- grs.calc.llk.tree(
        tree,prior,p.stay,p.switch,
        llk.data,
        TRUE, TRUE
    )
    
    f <- tmp$f;
    g <- tmp$g;
        
    b <- list()
    llk.tree.b <- array(0, c(nrow(tree), 2))
    colnames(llk.tree.b) <- c("LLK.0","LLK.1")
    b[[nrow(tree)]] <- list(op=prior$prior, lmx=0)

    for( i in (nrow(tree)-1):1 ) {
        
        np <- tree[i,"Par"]
        
        tmp <- log(b[[np]]$op) + log(f[[np]]$op) - log(g[[i]]$op)
        
        mx <- max(tmp)
        tmp <- exp(tmp-mx)
        val <- p.switch*sum(tmp)
        
        b.tmp <- p.stay*tmp + prior$prior*val
        mx2 <- max(b.tmp)
        
        b[[i]] <- list(
            op=b.tmp/mx2,
            lmx = mx + log(mx2) + b[[np]]$lmx + f[[np]]$lmx - g[[i]]$lmx
        )
    }

    bs.post.dec <- list()
    for (i in 1:nrow(tree)) {
        tmp <- grs.get.posterior.node_1d(
            f,
            b,
            prior, id=i, log.plot=T,plot=F, verbose=F);
        
        tmp[[4]] <- paste(tmp[[4]],collapse='-')
        tmp[[5]] <- paste(tmp[[5]],collapse='-')
        bs.post.dec[[i]] <- tmp
    
    }

    summed_llk <- c()
    for( i in 1:nrow(tree) ) {
        summed_llk[i] <- sum(f[[i]]$op * b[[i]]$op) + f[[i]]$lmx + b[[i]]$lmx   
    }

    out <- do.call(rbind,bs.post.dec)
    out$POST_ACTIVE <- as.numeric(out$POST_ACTIVE)
    out$max_b <- as.numeric(out$max_b)
    out$b_ci_lhs <- as.numeric(out$b_ci_lhs)
    out$b_ci_rhs <- as.numeric(out$b_ci_rhs)
    
    out <- cbind(tree,out)

    return(out)
}

grs.calc.lBF <- function(
    tree,
    prior,
    p.stay,
    p.switch,
    llk.data,
    eps = 1e-200
) {

    i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
    i.par <- setdiff(tree[,'ID'], i.ter);
    i.with.data <- which(tree$selectable %in% "Y")
    null.id <- which.min(abs(prior$b.grid))

    llk.full <- grs.calc.llk.tree(
        tree,prior,p.stay,p.switch,
        llk.data
    )

    ## calculate LLK under model with no active states and prior on this
    llk.full.null<-rep(0, nrow(tree));

    for (i in i.ter) {
        tmp <- llk.data[[i]]$op[null.id]
        if( tmp < eps ) tmp <- eps
        llk.full.null[i] <- log(tmp) + llk.data[[i]]$lmx
    }
    
    for (i in i.par) {
        w.d <- which(tree[,'Par']==i);
        if(! i %in% i.with.data ) {
            llk.full.null[i] <- sum(llk.full.null[w.d]) + length(w.d)*logp00;
        } else {
            j <- which(i.with.data == i)
            llk.full.null[i] <- sum(llk.full.null[w.d]) + length(w.d)*logp00 +
                log(llk.data[[j]]$op[null.id]) + llk.data[[j]]$lmx
        }
        
    }
    
    llk.full.null.mrca <- llk.full.null[nrow(tree)]+log(1-pi1);
    l.p.full.null <- log(1-pi1)+(nrow(tree)-1)*logp00;
    
    mrca <- nrow(tree);
    llk.full.mrca <- log(llk.full$llk[mrca,2])+llk.full$llk[mrca,3];
    
    ## Get BF
    tmp <- c(llk.full.mrca, llk.full.null.mrca);
    mx <- max(tmp);
    tmp2 <- mx + log(exp(tmp[1]-mx) - exp(tmp[2]-mx)) - llk.full.null.mrca;
    tmp3 <- l.p.full.null - log(1-exp(l.p.full.null));
    log10_treeBF <- (tmp2 + tmp3)/log(10)

    return(as.numeric(log10_treeBF))

}

grs.calc.llk.tree <- function(
    tree,
    prior,
    p.stat,
    p.switch,
    llk.data,
    returnF = FALSE,
    returnG = FALSE
) {

    i.ter <- tree[which(!(tree[,'ID'] %in% tree[,'Par'])),'ID'];
    i.par <- setdiff(tree[,'ID'], i.ter);
    i.with.data <- which(tree$selectable %in% "Y")
    null.id <- which.min(abs(prior$b.grid))
    if(returnG) g.surf.tree <- list();

    
    llk.tree <- array(0,c(nrow(tree),3))
    colnames(llk.tree) <- c("LLK.0","LLK.1","Max.LLK.alt")

    for ( i in i.ter ) {
        llk.tree[i,1] <- log(llk.data[[i]]$op[null.id]) + llk.data[[i]]$lmx
        llk.tree[i,2] <- sum( prior$prior * llk.data[[i]]$op )
        llk.tree[i,3] <- llk.data[[i]]$lmx
    }

    for ( i in i.par ) {
        ## find child nodes
        w.d <- which(tree[,'Par'] == i)
        tmp <- array(1, dim(llk.data[[2]]$op))
        s1 <- 0
        for( j in w.d ) {
            tmp.part <- p.stay * llk.data[[j]]$op + p.switch * llk.tree[j,2]
            tmp <- tmp * tmp.part
            if(returnG) {
                gmx <- max(tmp.part)
                g.surf.tree[[j]] <- list(
                    op = tmp.part/gmx, lmx = log(gmx) + llk.data[[j]]$lmx
                )
            }
            s1 <- s1 + llk.data[[j]]$lmx
        }
        if( i %in% i.with.data ) {
            j = which( i.with.data == i )
            tmp.part <- llk.data[[j]]$op
            tmp <- tmp * tmp.part
            s1 <- s1 + llk.data[[j]]$lmx
        }

        mx <- max(tmp)
        tmp <- tmp/mx

        llk.data[[i]] <- list(op=tmp, lmx = s1 + log(mx))
        llk.tree[i,1] <- log(llk.data[[i]]$op[null.id]) + llk.data[[i]]$lmx
        llk.tree[i,2] <- sum( prior$prior * llk.data[[i]]$op )
        llk.tree[i,3] <- llk.data[[i]]$lmx

    }
    
    if( returnF ) {
        if( ! returnG ) {
            return(list(llk=llk.tree,f=data.llk))
        } else {
            return(list(llk=llk.tree,f=llk.data, g = g.surf.tree))
        }
    } else {
        mrca <- nrow(tree)
        llk.full.mrca <- log(llk.tree[mrca,2]) + llk.tree[mrca,3]
        return(list(llk=llk.tree,llk.full.mrca=llk.full.mrca))
    }
}

grs.get.posterior.node_1d <- function(
    forward, backward, prior, id = 1, plot = FALSE, 
    return.ci = TRUE, verbose = FALSE, ci.level = 0.95, log.plot = TRUE
) {
    null.id <- which.min(abs(prior$b.grid))
        
    tmp <- forward[[id]]$op * backward[[id]]$op
    tmp <- tmp/sum(tmp)
    post.null <- tmp[null.id]
    post.active <- 1 - post.null
    tmp[null.id] <- 0
    tmp <- tmp/sum(tmp)
    mx <- arrayInd(which.max(tmp), dim(tmp))
    if (verbose) 
        cat("\nNode ", id)
    if (verbose) 
        cat("\nMax at b1 = ", prior$b.grid[mx[1]])
    if (verbose) 
        cat("\nSummed LLK = ", log(sum(forward[[id]]$op * backward[[id]]$op)) + 
            forward[[id]]$lmx + backward[[id]]$lmx)
    if (return.ci) {
        oo <- order(tmp, decreasing = T)
        cs <- cumsum(tmp[arrayInd(oo, dim(tmp))])
        w.ci <- oo[c(1, which(cs <= ci.level))]
        inds <- arrayInd(w.ci, dim(tmp))
        rg <- range(prior$b.grid[inds[, 1]])

        if (verbose) 
            cat("\nCI b1(", ci.level, ") = ", paste(rg, collapse = " - "), 
                sep = "")      
    }
    if (plot) {
        if (log.plot) 
            tmp <- log(tmp)
        plot(x = prior$b.grid, y = tmp, 
             main = paste("Node", id), xlab = "B1", ylab = "log(post)")
    }
    if (verbose) 
        cat("\n\n")

    out <- data.frame(
        max_b = prior$b.grid[mx[1]],
        summed_llk = log(sum(forward[[id]]$op * backward[[id]]$op)) + 
            forward[[id]]$lmx + backward[[id]]$lmx,
        b_ci_lhs = rg[1], 
        b_ci_rhs = rg[2],
        POST_ACTIVE = post.active)

    return(out)
}

