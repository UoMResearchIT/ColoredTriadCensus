

#' Triad adjacency matrices, as defined in Holland & Leinhardt (1976)
triad.types <- list(
  "003"  = matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 3),
  "012"  = matrix(c(0, 0, 0, 0, 0, 0, 1, 0, 0), nrow = 3, byrow = T),
  "102"  = matrix(c(0, 0, 0, 0, 0, 1, 0, 1, 0), nrow = 3),
  "021D" = matrix(c(0, 1, 1, 0, 0, 0, 0, 0, 0), nrow = 3, byrow = T),
  "021U" = matrix(c(0, 0, 0, 1, 0, 0, 1, 0, 0), nrow = 3, byrow = T),
  "021C" = matrix(c(0, 1, 0, 0, 0, 0, 1, 0, 0), nrow = 3, byrow = T),
  "111D" = matrix(c(0, 1, 0, 0, 0, 1, 0, 1, 0), nrow = 3, byrow = T),
  "111U" = matrix(c(0, 0, 0, 1, 0, 1, 0, 1, 0), nrow = 3, byrow = T),
  "030T" = matrix(c(0, 0, 0, 1, 0, 0, 1, 1, 0), nrow = 3, byrow = T),
  "030C" = matrix(c(0, 1, 0, 0, 0, 1, 1, 0, 0), nrow = 3, byrow = T),
  "201"  = matrix(c(0, 0, 1, 0, 0, 1, 1, 1, 0), nrow = 3),
  "120D" = matrix(c(0, 1, 1, 0, 0, 1, 0, 1, 0), nrow = 3, byrow = T),
  "120U" = matrix(c(0, 0, 0, 1, 0, 1, 1, 1, 0), nrow = 3, byrow = T),
  "120C" = matrix(c(0, 1, 0, 0, 0, 1, 1, 1, 0), nrow = 3, byrow = T),
  "210"  = matrix(c(0, 1, 0, 1, 0, 1, 1, 1, 0), nrow = 3, byrow = T),
  "300"  = matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)
)

#' @importFrom Matrix diag
tr=function(mat){return(sum(diag(mat)))}

#' @importFrom gtools permutations
txtmats=function(){
  perms=permutations(3,3,1:3)
  ret=list()
  col=rep(seq(1:6),6,each=1)
  row=rep(seq(1:6),1,each=6)
  for(i in 1:36) ret[[i]]=list(perms[row[i],],perms[col[i],])
  return(ret)
}

samemat=function(mat1,mat2){
  txtmat=txtmats()
  testmats=list()
  check=F
  for(i in 1:length(txtmat)) testmats[[i]]=mat2[txtmat[[i]][[1]],txtmat[[i]][[2]]]
  for(i in 1:length(testmats)) if(all(mat1==testmats[[i]])) check=T
  return(check)
}

makeUniqueTri=function(triangle,colorCombs){
  out=1:nrow(colorCombs)
  uniqueTriangles=list()
  for(i in 1:nrow(colorCombs)){uniqueTriangles[[i]]= triangle ; diag(uniqueTriangles[[i]])=colorCombs[i,]+1}
  if(length(uniqueTriangles)>1){
    for(i in (length(uniqueTriangles)-1):1){
      for (j in (i+1):length(uniqueTriangles)){
        if(samemat(uniqueTriangles[[i]],uniqueTriangles[[j]])) out[j]=i
      }
    }
  }
  return(out)
}

#' Classification index for colored triads of a given `type`
#'
#' @param type triad type code, must be in `names(triad.types)`
#' @param ncolors max. number of colors for the graph
#'
#' @return `ncolors^3` integer vector indicating a unique identifier
#'  for each color permutation. Equivalent permutations will share the
#'  same index.
#' @importFrom gtools permutations
unique.triad.index <- function(type, ncolors, use.cache = TRUE) {

  # TODO: #4 cache e.g. using memoise, instead of shipping uniTri
  # TODO: #5 for ncolors > 3, surely it's more efficient to run makeUniqueTri
  #  for 3 colors, and then map all other permutations using that index.

  imp <- uniTri[[paste(ncolors, type, sep = " ")]]

  if (is.null(imp) | !use.cache) {
    triangle <- triad.types[[type]]
    colorCombs <- permutations(ncolors, 3, repeats.allowed = TRUE)
    imp <- makeUniqueTri(triangle = triangle, colorCombs)
  }
  return(imp)
}


remove.isomorphic <- function(triad, colorCombs, counts) {

  imp <- unique.triad.index(triad, max(colorCombs))
  stopifnot(length(counts) == length(imp))

  T003Col2 <- c()
  for (i in unique(imp))
  {
    T003Col2[length(T003Col2) + 1] <- counts[imp == i][1]
    names(T003Col2)[length(T003Col2)] <- names(counts)[imp == i][1]
  }

  for (i in unique(imp))
  {
    if (triad %in% c("003", "300")) {
      if (
        colorCombs[i, 1] == colorCombs[i, 2] |
        colorCombs[i, 1] == colorCombs[i, 3] |
        colorCombs[i, 2] == colorCombs[i, 3]
      ) {
        count <- count / 2
      }
      if (
        colorCombs[i, 1] == colorCombs[i, 2] &
        colorCombs[i, 1] == colorCombs[i, 3] &
        colorCombs[i, 2] == colorCombs[i, 3]
      ) {
        count <- count / 3
      }
    }

    if (triad == "030C")
      if (colorCombs[i, 1] == colorCombs[i, 2] &
          colorCombs[i, 1] == colorCombs[i, 3] &
          colorCombs[i, 2] == colorCombs[i, 3])
        count <- count / 3

    if (triad %in% c("102", "021D", "021U", "201", "120D", "120U"))
      if (colorCombs[i, 3] == colorCombs[i, 2])
        count <- count / 2
  }
  # print(triad)
  return(T003Col2)
}


#' @importFrom Matrix Matrix t diag
get.triad.matrices <- function(mat) {

  #Adjacency matrix
  testA=mat[,]
  #Asymmetric edges only
  testC=testA-t(testA)
  testC=matrix(pmax(0,testC),nrow=nrow(mat),byrow=F)
  #Mutual edges only
  testM=testA-testC
  #Symmetrized full matrix, and its complement
  testE=testA+t(testA)
  testE=pmin(testE,1)
  testE0=1-testE
  diag(testE0)=0

  testA=Matrix(testA)
  testC=Matrix(testC)
  testM=Matrix(testM)
  testE=Matrix(testE)
  testE0=Matrix(testE0)

  triad.matrices <- list(
    "003" = list(HT12 = testE0, HT23 = testE0, HT31 = testE0),
    "012" = list(HT12 = testE0, HT23 = testE0, HT31 = testC),
    "102" = list(HT12 = testE0, HT23 = testM, HT31 = testE0),
    "021D" = list(HT12 = testC, HT23 = testE0, HT31 = t(testC)),
    "021U" = list(HT12 = t(testC), HT23 = testE0, HT31 = testC),
    "021C" = list(HT12 = testC, HT23 = testE0, HT31 = testC),
    "111D" = list(HT12 = testC, HT23 = testM, HT31 = testE0),
    "111U" = list(HT12 = t(testC), HT23 = testM, HT31 = testE0),
    "030T" = list(HT12 = t(testC), HT23 = t(testC), HT31 = testC),
    "030C" = list(HT12 = testC, HT23 = testC, HT31 = testC),
    "201" = list(HT12 = testM, HT23 = testE0, HT31 = testM),
    "120D" = list(HT12 = testC, HT23 = testM, HT31 = t(testC)),
    "120U" = list(HT12 = t(testC), HT23 = testM, HT31 = testC),
    "120C" = list(HT12 = testC, HT23 = testM, HT31 = testC),
    "210" = list(HT12 = testM, HT23 = testM, HT31 = testC),
    "300" = list(HT12 = testM, HT23 = testM, HT31 = testM)
  )
}


#' @param mat adjacency matrix
#' @param col vector of colors for each node
#' @param color.set complete/ordered set of colors
#' @param directed whether to count directed triads
#' @return a named vector with counts of colored triads
#'
#' @importFrom Matrix Matrix t diag
#' @importFrom gtools permutations
#' @export
colored.triad.census <- function(mat, col, color.set = unique(col), directed=FALSE){

  stopifnot(all(col %in% color.set))
  stopifnot(nrow(mat) == ncol(mat))
  stopifnot(length(col) == nrow(mat))
  stopifnot(all(mat %in% c(0, 1)))

  colnum <- length(color.set)
  colorIndex <- match(col, color.set)

  # Make out-coloring matrices, evaluating the color of the nodes row-wise
  # the complete row cmat[[i]][j,] == 1 if node j has color i
  cmat <- list()
  for (i in 1:colnum) {
    cmat[[i]] <- Matrix(rep(colorIndex == i, length(col)), nrow = length(col), byrow = FALSE) * 1
    cmat[[i]][, ]
  }
  colorCombs <- permutations(colnum, 3, repeats.allowed = TRUE)

  # count triads (with redundancies)
  triad.count <- function(triad) {

    HT12 <- triad.matrices[[triad]]$HT12
    HT23 <- triad.matrices[[triad]]$HT23
    HT31 <- triad.matrices[[triad]]$HT31

    counts <- c()
    for (i in 1:nrow(colorCombs)) {
      counts[i] <- tr(
        (t(cmat[[colorCombs[i, 2]]]) * cmat[[colorCombs[i, 1]]] * HT12) %*%
        (t(cmat[[colorCombs[i, 3]]]) * cmat[[colorCombs[i, 2]]] * HT23) %*%
        (t(cmat[[colorCombs[i, 1]]]) * cmat[[colorCombs[i, 3]]] * HT31)
      )
    }
    names(counts) <- apply(colorCombs, 1, function(x) paste0(x, collapse = ""))
    names(counts) <- lapply(names(counts), function(x) paste0("T", triad, "-", x))

    remove.isomorphic(triad, colorCombs, counts)
  }

  if (directed) {
    types <- names(triad.matrices)  # use all types
  } else {
    types <- c("003", "102", "201", "300")
  }

  all_counts <- lapply(types, triad.count)
  do.call(c, all_counts)
}


