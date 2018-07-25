
sshhh <- function(a.package){
  suppressWarnings(suppressPackageStartupMessages(
    library(a.package, character.only=TRUE)))
}

pkgs <- c("data.table", "plotly", "Rgraphviz", "graph", "BiocGenerics", "parallel", "magrittr")
loads <- sapply(pkgs, sshhh)

###################################################################/
# Descrption: mi.univariate.subset
#
###################################################################/
mi.univariate.subset = function( )
{
  selection = c( "I251", "E780", "M206", "S9211", "M4792", "I422", "C837", "I461", "B24", "G35")
  return( mi.univariate()[ Code %in% selection ][ order( pValue )] )
}

###################################################################/
# Descrption: mi.univariate
#
###################################################################/
mi.univariate = function( originalCols = FALSE )
{
  file = 'MI_GRS_tree.rdata'
  load( file )
  data = as.data.table( MI_GRS_tree )

  if( originalCols == TRUE )
    return( data )

  setnames( data, c( "coding", "meaning", "Naffected", "BETA", "Pval" ), c( "Code", "Description", "N", "beta", "pValue") )
  return( data[, .( Code, Description, N, beta, pValue ) ][ order( pValue ) ] )
}

###################################################################/
# Descrption: ABO.univariate
#
###################################################################/
ABO.univariate = function()
{
  file = 'ABO.1df.res.rdata'
  load( file )
  data = as.data.table( data )
  return( data[ !is.na( pValue ), .( Code = coding, Description = meaning, N = COUNTS, beta = log( OR ), pValue ) ][ order( pValue ) ] )
}

###################################################################/
# Descrption: Data object with likelihood surfaces for the ABO
# SNP in the ICD-10 UK Biobank data set
###################################################################/
ABO.lk.surfs = function( ) {
    file = 'ABO.lk.surfs.rdata'
    load( file )
    return( ABO.data )
}

ABO.lk.pars = function( ) {
    file = 'ABO.lk.surfs.rdata'
    load( file )
    return( ABO.pars )
}

###################################################################/
# Descrption: Data object with likelihood surfaces for the MI GRS
# in the ICD-10 UK Biobank data set
###################################################################/
MI_GRS.lk.surfs = function( ) {
    file = 'MI_GRS_UKBB_data.rdata'
    load( file )
    return( list( llk.data = llk.data, tree = tree, prior = prior ) )
}

###################################################################
## Descrption: function to draw the likelihood surface in the GRS
## analysis
###################################################################
grs.plot.lk.surface <- function(
    code = NULL,
    tree = NULL,
    prior = NULL,
    lk.surfs = NULL
) {

    if ( is.null(code) | is.null(tree) | is.null(prior) | is.null(lk.surfs) )
        stop("missing function argument.\n")
    
    code.idx <- which(tree$coding %in% code)
    code.meaning <- tree[code.idx,'meaning']
    
    plot(
        prior$b.grid,
        lk.surfs[[code.idx]]$op,
        xlab = 'beta',
        ylab = 'likelihood (scaled)',
        pch = 19,
        col = 'black',
        bty = 'l',
        main = code.meaning
    )

}

###################################################################/
# Descrption: draw_tree_univariate
###################################################################/
draw_tree_univariate <- function(
  data,
  title            = "Univariate Analysis",
  pValueThreshold  = 1e-5,
  pValueSaturation = 1e-50
)
{
  # make sure data is correct form
  if( !is.data.table( data ) )
    data = as.data.table( data )

  if( length( setdiff( c( "Code", "pValue", "beta" ), names( data ) ) != 0 ) )
    throw( "input data requires pValue, beta and  Code columns")
  data = data[ ,.( coding = Code, Pval = pValue, BETA = beta ) ]

  # get whole tree data
  tree = mi.univariate( originalCols = TRUE )
  tree = tree[ ,.( ID, Par, coding, meaning ) ]
  tree = data[ tree, on = "coding" ]

  # convert threshold to log
  lThresh = log( pValueThreshold ) / log( 10 )
  lSat    = log( pValueSaturation) / log( 10 )
  if( lSat > lThresh )
    lWidth = 1e-5
  else
    lWidth = lThresh - lSat

  # convert data to form required by Adrian's function
  pp = data.table( tree )[ , .(  BETA = ifelse( is.na( BETA ), 0, BETA  ), logPval = ifelse( is.na( Pval), 0, ifelse( Pval < 1e-185, -185, log( Pval ) / log( 10 ) ) ) ) ]
  pp[ , effPp := ifelse( logPval > lThresh, 0, pmin( 1, ( lThresh - logPval ) / lWidth ) ) ]
  pp[ , effPp := ifelse( effPp > 0, effPp / 2 + 0.5, 0 ) ]

  pp = as.matrix( pp[ , .( ifelse( BETA < 0, effPp, 0 ), 1 - effPp, ifelse( BETA > 0, effPp, 0 ) ) ] )
  # finally drawer tree
  draw_tree( as.data.frame( tree ), pp, tree_title = title,trim_tree_pp = 0.01, measureName = "pValue", measureValueFunc = function( t, p ) return( as.data.table( t )[ , format( Pval, digits = 3 ) ] ) )$plot
}

###################################################################/
# Descrption: draw_tree
###################################################################/
draw_tree <- function(
  tree           = NULL,
  pp             = NULL,
  tree_title     = "Tree",
  only.get.stats = FALSE,
  trim_tree_pp   = NULL,
  measureName    = "PP_active",
  measureValueFunc = function( tree, pp ) return( round( 1 - pp[ ,2 ] , 2) )
)
{
  # remove tree if posterior probability is to low
  if( !is.null( trim_tree_pp ) )
  {
    tmp  <- trim_tree( tree = tree, pp = pp, pp.thr = trim_tree_pp )
    tree <- tmp$tree
    pp   <- tmp$pp
  }

  matrix <- matrix(0, ncol = nrow(tree), nrow = nrow(tree))

  for( i in 1:( nrow( tree ) - 1 ) )
  {
    p <- tree[i,"Par"]
    c <- tree[i,"ID"]
    matrix[p,c] <- 1
  }

  rownames(matrix) <- tree$ID
  colnames(matrix) <- tree$ID
  labels = tree$ID

  graph <- new("graphAM", adjMat = matrix, edgemode = 'directed')

  lGraph <- layoutGraph(graph)
  ninfo <- nodeRenderInfo(lGraph)

  node_state <- apply(pp,1,function(x) return( which.max(x) )) - 2

  node_labels <- paste(
    tree$meaning,
    "<br>",
    "State: ",node_state,
    "<br>",
    measureName, "= ", measureValueFunc( tree, pp ),
    sep = ""
  )

  nodeRI = data.frame(
    NODE    = names(ninfo$nodeX),
    PP1     = pp[,1],
    PP2     = pp[,2],
    PP3     = pp[,3],
    NODEX   = ninfo$nodeX,
    NODEY   = ninfo$nodeY,
    MEANING = tree$meaning,
    LABEL   = node_labels
  )

  col_pal_risk <- colorRampPalette( c( "white", rgb(112, 28, 28, max=255) ) )( 100 )
  col_pal_prot <- colorRampPalette( c( "white", rgb(8, 37, 103, max=255) ) )(100)

  cols1 <- map2color(pp[,3], col_pal_risk, limits = c(0,1))
  cols2 <- map2color(pp[,1], col_pal_prot, limits = c(0,1))

  bar.cols <- rep("white",nrow(tree))
  state.col <- apply(pp[,c(1,3)],1,which.max)
  bar.cols[state.col == 1] <- cols2[state.col == 1]
  bar.cols[state.col == 2] <- cols1[state.col == 2]

  cols <- rep("white",nrow(tree))
  idx <- which(node_state == 1)
  cols[idx] <- cols1[idx]
  idx <- which(node_state == -1)
  cols[idx] <- cols2[idx]

  nodeRI$COL <- bar.cols

  attrs <- list(node = list(fillcolor = 'white'), edge = list(arrowsize=0.5))

  names(cols) <- labels

  nattrs <- list(fillcolor=cols)

  nodes <- buildNodeList(graph, nodeAttrs=nattrs, defAttrs=attrs$node)
  edges <- buildEdgeList(graph)
  vv <- agopen(name="foo", nodes=nodes, edges=edges, attrs=attrs,
    edgeMode="directed")

  x <- vv
  y <- x@layoutType
  x <- graphLayout(x, y)
  ur <- upRight(boundBox(x))
  bl <- botLeft(boundBox(x))

  out <- list()
  out$nodeRI <- nodeRI

  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ## Initalize plotly
  ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p <- plotly::plot_ly()

  xlim1 <- getX(ur)*1.02
  xlim0 <- -xlim1*0.02
  xlim <- c(xlim0, xlim1)

  ## Add an axis.
  p = plotly::layout(
    p,
    title      = tree_title,
    xaxis      = list(
      title          = "",
      showgrid       = FALSE,
      showticklabels = FALSE,
      showline       = FALSE,
      zeroline       = FALSE,
      range          = xlim
    ),
    yaxis      = list(
      title          = "",
      showgrid       = FALSE,
      showticklabels = FALSE,
      showline       = FALSE,
      zeroline       = FALSE,
      range          = c(getY(bl), getY(ur))
    ),
    showlegend = FALSE
  )

  out$xlim = xlim

  ## Add the edges
  edges <- AgEdge(x)
  edges.p <- list()

  for( i in 1:length(edges) ) {

    edge <- edges[[i]]
    node.to <- edge@head
    node.from <- edge@tail

    for ( j in 1:length(splines(edge)) ) {
      z <- splines(edge)[[j]]
      points <- matrix(unlist(pointList(z)),ncol=2,byrow=TRUE)

      p <- add_trace(
        p,
        x          = points[,1],
        y          = points[,2],
        type       = "scatter",
        mode       = "lines",
        hoverinfo  = "none",
        line       = list(color = "gray"),
        showlegend = FALSE
      )
    }

    edges.p[[i]] <- points
    heads     = bezierPoints(z)
    head_from = heads[nrow(heads)-1, ]
    head_to   = heads[nrow(heads),]
}

  ## Add the nodes
  order <- order(pp[,2],decreasing=TRUE)
  p = plotly::add_trace( p,
    x          = nodeRI$NODEX[order],
    y          = nodeRI$NODEY[order],
    type       = "scatter",
    mode       = "markers",
    text       = nodeRI$LABEL[order],
    hoverinfo  = "text",
    marker     = list(
      size   = 20,
      symbol = "circle",
      color = cols[order],
      line  = list( color = "black", width = 1)
    ),
    showlegend = FALSE
  )

  out$plot <- p
  out$edges <- edges.p

  return(out)
}

###################################################################/
# Descrption: map2color
###################################################################/
map2color <- function( x, pal, limits = range(x) )
{
  return( pal[ findInterval(
    x,
    seq(limits[1],limits[2],length.out=length(pal)+1),
    all.inside=TRUE
  ) ] )
}

###################################################################/
# Descrption: trim_tree
###################################################################/
trim_tree <- function( tree = tree, pp = pp, pp.thr = 0.75 )
{
  idx <- which(pp[,2] < 1 - pp.thr)

  t2 <- tree[idx,]
  siblings <- unique(unlist(lapply(t2$ID,get_tree_siblings,tree)))
  paths_to_root <- unique(unlist(lapply(t2$ID,get_path_ids_to_root,tree)))

  nodes_to_keep <- sort(unique(c(t2$ID,siblings,paths_to_root)),decreasing=F)

  t2 <- tree[tree$ID %in% nodes_to_keep, ]
  pp2 <- pp[tree$ID %in% nodes_to_keep,]

  new_id <- 1:nrow(t2)
  new_par <- new_id[match(t2$Par,t2$ID)]

  t2$ID <- new_id
  t2$Par <- new_par

  t2[nrow(t2),'Par'] <- 0
  o <- list(tree=t2,pp=pp2)

  return(o)
}

###################################################################/
# Descrption: get_tree_siblings
###################################################################/
get_tree_siblings <- function(id,tree)
{
  par_id <- tree[ tree$ID %in% id, "Par"]
  sibling_ids <- tree[ tree$Par %in% par_id, "ID"]

  return(sibling_ids)
}

###################################################################/
# Descrption: get_path_ids_to_root
###################################################################/
get_path_ids_to_root <- function( id, tree )
{
  out     <- id
  root_id <- tree[ nrow(tree), "ID"]

  while( ! root_id %in% out )
    out <- unique( c ( out, tree[ tree$ID %in% out, "Par" ] ) )

  return(out)
}

grs.trim_tree <- function( tree = tree, pp = pp, pp.thr = 0.75 )
{
  idx <- which(pp$POST_ACTIVE >= pp.thr)

  t2 <- tree[idx,]
  siblings <- unique(unlist(lapply(t2$ID,get_tree_siblings,tree)))
  paths_to_root <- unique(unlist(lapply(t2$ID,get_path_ids_to_root,tree)))

  nodes_to_keep <- sort(unique(c(t2$ID,siblings,paths_to_root)),decreasing=F)

  t2 <- tree[tree$ID %in% nodes_to_keep, ]
  pp2 <- pp[tree$ID %in% nodes_to_keep,]

  new_id <- 1:nrow(t2)
  new_par <- new_id[match(t2$Par,t2$ID)]

  t2$ID <- new_id
  t2$Par <- new_par

  t2[nrow(t2),'Par'] <- 0
  o <- list(tree=t2,pp=pp2)

  return(o)
}

###################################################################/
# Descrption: GRS draw_tree
###################################################################/
grs.draw_tree <- function(
  tree           = NULL,
  pp             = NULL,
  tree_title     = "GRS Tree",
  only.get.stats = FALSE,
  trim_tree_pp   = NULL
) {

    if ( ! is.null( trim_tree_pp ) ) {
        tmp  <- grs.trim_tree( tree = tree, pp = pp, pp.thr = trim_tree_pp )
        tree <- tmp$tree
        pp   <- tmp$pp
    }

    matrix <- matrix(0, ncol = nrow(tree), nrow = nrow(tree))

    for( i in 1:( nrow( tree ) - 1 ) ) {
        p <- tree[i,"Par"]
        c <- tree[i,"ID"]
        matrix[p,c] <- 1
    }

    rownames(matrix) <- tree$ID
    colnames(matrix) <- tree$ID
    labels = tree$ID
    
    graph <- new("graphAM", adjMat = matrix, edgemode = 'directed')

    lGraph <- layoutGraph(graph)
    ninfo <- nodeRenderInfo(lGraph)

    node_labels <- paste(
        tree$meaning,"<br>",
        "beta: ",round(pp$max_b,2),"<br>",
        "PP: ",round(pp$POST_ACTIVE,2),"<br>",
        sep = ""
    )

    nodeRI = data.frame(
        NODE    = names(ninfo$nodeX),
        PP      = as.numeric(pp$POST_ACTIVE),
        NODEX   = ninfo$nodeX,
        NODEY   = ninfo$nodeY,
        MEANING = tree$meaning,
        LABEL   = node_labels
    )

    col_pal_risk <- colorRampPalette( c( "white", rgb(112, 28, 28, max=255) ) )( 100 )
    col_pal_prot <- colorRampPalette( c( "white", rgb(8, 37, 103, max=255) ) )(100)

    cols1 <- map2color(pp$POST_ACTIVE, col_pal_risk, limits = c(0,1))
    cols2 <- map2color(pp$POST_ACTIVE, col_pal_prot, limits = c(0,1))

    bar.cols <- rep("white",nrow(tree))
    bar.cols[pp$max_b < 0] <- cols2[pp$max_b < 0]
    bar.cols[pp$max_b > 0] <- cols1[pp$max_b > 0]
    
    nodeRI$COL <- bar.cols

    attrs <- list(node = list(fillcolor = 'white'), edge = list(arrowsize=0.5))

    cols <- nodeRI$COL
    names(cols) <- labels
    
    nattrs <- list(fillcolor=cols)
    
    nodes <- buildNodeList(graph, nodeAttrs=nattrs, defAttrs=attrs$node)
    edges <- buildEdgeList(graph)
    vv <- agopen(name="foo", nodes=nodes, edges=edges, attrs=attrs,
                 edgeMode="directed")

    x <- vv
    y <- x@layoutType
    x <- graphLayout(x, y)
    ur <- upRight(boundBox(x))
    bl <- botLeft(boundBox(x))

    out <- list()
    out$nodeRI <- nodeRI

    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ## Initalize plotly
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p <- plotly::plot_ly()
    
    xlim1 <- getX(ur)*1.02
    xlim0 <- -xlim1*0.02
    xlim <- c(xlim0, xlim1)

    ## Add an axis.
    p = plotly::layout(
        p,
        title      = tree_title,
        xaxis      = list(
            title          = "",
            showgrid       = FALSE,
            showticklabels = FALSE,
            showline       = FALSE,
            zeroline       = FALSE,
            range          = xlim
        ),
        yaxis      = list(
            title          = "",
            showgrid       = FALSE,
            showticklabels = FALSE,
            showline       = FALSE,
            zeroline       = FALSE,
            range          = c(getY(bl), getY(ur))
        ),
        showlegend = FALSE
    )

    out$xlim = xlim

    ## Add the edges
    edges <- AgEdge(x)
    edges.p <- list()
    
    for( i in 1:length(edges) ) {
        
        edge <- edges[[i]]
        node.to <- edge@head
        node.from <- edge@tail
        
        for ( j in 1:length(splines(edge)) ) {
            z <- splines(edge)[[j]]
            points <- matrix(unlist(pointList(z)),ncol=2,byrow=TRUE)
        
            p <- add_trace(
                p,
                x          = points[,1],
                y          = points[,2],
                type       = "scatter",
                mode       = "lines",
                hoverinfo  = "none",
                line       = list(color = "gray"),
                showlegend = FALSE
            )
        }
        
        edges.p[[i]] <- points
        heads     = bezierPoints(z)
        head_from = heads[nrow(heads)-1, ]
        head_to   = heads[nrow(heads),]
    }
    
    ## Add the nodes
    order <- order(nodeRI$PP,decreasing=F)

    tmp <- nodeRI[order,]
    
    p = plotly::add_trace(
        p,
        x          = nodeRI$NODEX[order],
        y          = nodeRI$NODEY[order],
        type       = "scatter",
        mode       = "markers",
        text       = nodeRI$LABEL[order],
        hoverinfo  = "text",
        marker     = list(
            size   = 20,
            symbol = "circle",
            color = nodeRI$COL[order],
            line  = list( color = "black", width = 0.5)
        ),
        showlegend = FALSE
    )

    return(p)
}
###################################################################/
# Descrption: plot.qq
#
# assume the data is normally distributed, convert to normal and
# plot against theoretical z-score
###################################################################/
plot.qq = function( data, annot.col = 'Description' )
{
    require(plotly)

    ## remove missing data
    data <- na.omit(data)

    ## convert the pValues to a normal standard error
    data[ , zData :=-qnorm( pValue )]

    ## now add theoritical z vales
    data = data[ order( zData ) ]
    data[ , zTheor := qnorm( (1:data[,.N ] - 0.5 )/data[,.N ] ) ]


    ## make plot
    plot = plot_ly( data, x = ~zTheor, y = ~zData, type = "scatter", mode = "markers", name = "data", text = data[, get( annot.col )], hoverinfo = "text" )
    plot = layout( plot,
      xaxis = list( title = "theoretical z-score", range = list( floor( data[,min( zTheor ) ] ), ceiling( data[, max( zTheor ) ] ) ) ),
      yaxis = list( title = "data z-score" ),
      title = "Q-Q Plot"
    )
    plot = add_lines( plot, x = ~zData, y = ~zData, name = "null" )

  return( plot )
}
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

    p00 <- p.stay+p.switch*(1-pi1);
    logp00 <- log(p00);
    
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

