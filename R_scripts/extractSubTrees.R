# https://stackoverflow.com/questions/33084860/sampling-subgraphs-from-different-sizes-using-igraph

library(igraph)

tmp.neigh <- V(g)[unlist(neighborhood(g,1,seed))]$name
tmp.neigh <- setdiff(tmp.neigh, seed)


library(doParallel)
cl <- makeCluster(4)
registerDoParallel(cl)


josilber <- function(size, num.rep, G) {
  score_fun <- function(vert) sum(vert$weight*vert$RWRNodeweight)/sqrt(sum(vert$RWRNodeweight^2))
  n <- length(V(G)$name)

  # Determine which nodes fall in sufficiently large connected components
  comp <- components(G)
  valid <- which(comp$csize[comp$membership] >= size)

  perm <- replicate(num.rep, {
    first.node <- sample(valid, 1)
    used <- (1:n) == first.node  # Is this node selected?
    neigh <- (1:n) %in% neighbors(G, first.node)  # Does each node neighbor our selections?
    for (iter in 2:size) {
      new.node <- sample(which(neigh & !used), 1)
      used[new.node] <- TRUE
      neigh[neighbors(G, new.node)] <- TRUE
    }
    score_fun(V(G)[used])
  })
  perm
}

josilber2 <- function(size, num.rep, G) {

  # Determine which nodes fall in sufficiently large connected components
  comp <- components(G)
  valid <- which(comp$csize[comp$membership] >= size)
  n <- length(V(G))

  # for each subgraph
  perm <- foreach (i=1:num.rep) %dopar% {
    library(igraph)

    first.node <- sample(valid, 1)
    used <- (1:n) == first.node  # Is this node selected?
    neigh <- (1:n) %in% neighbors(G, first.node)  # Does each node neighbor our selections?

    for (iter in 2:size) {
      if (sum(neigh & !used) > 0) {
        new.node <- sample(which(neigh & !used), 1)
        used[new.node] <- TRUE
        neigh[neighbors(G, new.node)] <- TRUE
      }
    }
    which(used)
  }
  if (length(valid) > 0) {  # make sure there's valid nodes to sample
    perm
  } else {
    return(NULL)
  }
}

g <- erdos.renyi.game(1000, 1/1000, directed=F)
josilber2(10, 10, g)
