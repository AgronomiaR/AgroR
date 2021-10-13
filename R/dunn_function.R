#' Analysis: Post-hoc Dunn
#'
#' @description Perform Kruskal wallis and dunn post-hoc test
#' @param resp Vector with response
#' @param trat Numerical or complex vector with treatments
#' @param method the p-value for multiple comparisons ("none", "bonferroni", "sidak", "holm", "hs", "hochberg", "bh", "by"). The default is no adjustment for multiple comparisons
#' @param alpha Significance level of the post-hoc (\emph{default} is 0.05)
#' @param decreasing Should the order of the letters be increasing or decreasing.
#' @return Kruskal-wallis and dunn's post-hoc test returns
#' @importFrom utils capture.output
#' @importFrom dunn.test dunn.test
#' @export
#' @author Gabriel Danilo Shimizu, \email{shimizu@uel.br}
#' @examples
#' library(AgroR)
#' data(pomegranate)
#'
#' with(pomegranate, dunn(trat, WL))

dunn=function(trat,
              resp,
              method="holm",
              alpha=0.05,
              decreasing=TRUE)
  {requireNamespace("dunn.test")
  dtres <- capture.output(res <- dunn.test::dunn.test(resp, trat, method,
                                                             kw = TRUE,
                                                             altp = TRUE))
  res <- data.frame(res[-which(names(res) == "chi2")])[,c(4, 1, 2, 3)]
  names(res) <- c("Comparison", "Z", "P.unadj", "P.adj")
  vec2mat2=function (x, sep = "-"){splits <- strsplit(x, sep)
    n.spl <- sapply(splits, length)
    if (any(n.spl != 2))
      stop("Names must contain exactly one '", sep, "' each;  instead got ",
           paste(x, collapse = ", "))
    x2 <- t(as.matrix(as.data.frame(splits)))
    dimnames(x2) <- list(x, NULL)
    x2}
  multcompLetters=function (x, compare = "<",
                            threshold = alpha,
                            Letters = c(letters, LETTERS, "."),
                            reversed = decreasing)
  {x.is <- deparse(substitute(x))
    if (any(class(x) == "dist"))
      x <- as.matrix(x)
    if (!is.logical(x))
      x <- do.call(compare, list(x, threshold))
    dimx <- dim(x)
    {
      if ((length(dimx) == 2) && (dimx[1] == dimx[2])) {
        Lvls <- dimnames(x)[[1]]
        if (length(Lvls) != dimx[1])
          stop("Names requred for ", x.is)
        else {
          x2. <- t(outer(Lvls, Lvls, paste, sep = ""))
          x2.n <- outer(Lvls, Lvls, function(x1, x2) nchar(x2))
          x2.2 <- x2.[lower.tri(x2.)]
          x2.2n <- x2.n[lower.tri(x2.n)]
          x2a <- substring(x2.2, 1, x2.2n)
          x2b <- substring(x2.2, x2.2n + 1)
          x2 <- cbind(x2a, x2b)
          x <- x[lower.tri(x)]
        }
      }
      else {
        namx <- names(x)
        if (length(namx) != length(x))
          stop("Names required for ", x.is)
        x2 <- vec2mat2(namx)
        Lvls <- unique(as.vector(x2))}}
    n <- length(Lvls)
    LetMat <- array(TRUE, dim = c(n, 1), dimnames = list(Lvls, NULL))
    k2 <- sum(x)
    if (k2 == 0) {
      Ltrs <- rep(Letters[1], n)
      names(Ltrs) <- Lvls
      dimnames(LetMat)[[2]] <- Letters[1]
      return(list(Letters = Ltrs, LetterMatrix = LetMat))}
    distinct.pairs <- x2[x, , drop = FALSE]
    absorb <- function(A.) {
      k. <- dim(A.)[2]
      if (k. > 1) {
        for (i. in 1:(k. - 1)) for (j. in (i. + 1):k.) {
          if (all(A.[A.[, j.], i.])) {
            A. <- A.[, -j., drop = FALSE]
            return(absorb(A.))}
          else {
            if (all(A.[A.[, i.], j.])) {
              A. <- A.[, -i., drop = FALSE]
              return(absorb(A.))
            }
          }
        }
      }
      A.
    }
    for (i in 1:k2) {
      dpi <- distinct.pairs[i, ]
      ijCols <- (LetMat[dpi[1], ] & LetMat[dpi[2], ])
      if (any(ijCols)) {
        A1 <- LetMat[, ijCols, drop = FALSE]
        A1[dpi[1], ] <- FALSE
        LetMat[dpi[2], ijCols] <- FALSE
        LetMat <- cbind(LetMat, A1)
        LetMat <- absorb(LetMat)
      }
    }
    sortCols <- function(B) {
      firstRow <- apply(B, 2, function(x) which(x)[1])
      B <- B[, order(firstRow)]
      firstRow <- apply(B, 2, function(x) which(x)[1])
      reps <- (diff(firstRow) == 0)
      if (any(reps)) {
        nrep <- table(which(reps))
        irep <- as.numeric(names(nrep))
        k <- dim(B)[1]
        for (i in irep) {
          i. <- i:(i + nrep[as.character(i)])
          j. <- (firstRow[i] + 1):k
          B[j., i.] <- sortCols(B[j., i., drop = FALSE])
        }
      }
      B
    }
    LetMat. <- sortCols(LetMat)
    if (reversed)
      LetMat. <- LetMat.[, rev(1:ncol(LetMat.))]
    k.ltrs <- dim(LetMat.)[2]
    makeLtrs <- function(kl, ltrs = Letters) {
      kL <- length(ltrs)
      if (kl < kL)
        return(ltrs[1:kl])
      ltrecurse <- c(paste(ltrs[kL], ltrs[-kL], sep = ""),
                     ltrs[kL])
      c(ltrs[-kL], makeLtrs(kl - kL + 1, ltrecurse))
    }
    Ltrs <- makeLtrs(k.ltrs, Letters)
    dimnames(LetMat.)[[2]] <- Ltrs
    LetVec <- rep(NA, n)
    names(LetVec) <- Lvls
    for (i in 1:n) LetVec[i] <- paste(Ltrs[LetMat.[i, ]], collapse = "")
    nch.L <- nchar(Ltrs)
    blk.L <- rep(NA, k.ltrs)
    for (i in 1:k.ltrs) blk.L[i] <- paste(rep(" ", nch.L[i]),
                                          collapse = "")
    monoVec <- rep(NA, n)
    names(monoVec) <- Lvls
    for (j in 1:n) {
      ch2 <- blk.L
      if (any(LetMat.[j, ]))
        ch2[LetMat.[j, ]] <- Ltrs[LetMat.[j, ]]
      monoVec[j] <- paste(ch2, collapse = "")
    }
    InsertAbsorb <- list(Letters = LetVec, monospacedLetters = monoVec,
                         LetterMatrix = LetMat.)
    class(InsertAbsorb) <- "multcompLetters"
    InsertAbsorb}
  cldList=function (formula = NULL, data = NULL, comparison = NULL, p.value = NULL,
                    threshold = alpha, print.comp = FALSE, remove.space = TRUE,
                    remove.equal = TRUE, remove.zero = TRUE, swap.colon = TRUE,
                    swap.vs = FALSE){if (!is.null(formula)) {
      p.value = eval(parse(text = paste0("data", "$",
                                         all.vars(formula[[2]])[1])))
      comparison = eval(parse(text = paste0("data", "$",
                                            all.vars(formula[[3]])[1])))}
    Comparison = (as.numeric(p.value) <= threshold)
    if (sum(Comparison) == 0) {stop("No significant differences.", call. = FALSE)}
    if (remove.space == TRUE) {comparison = gsub(" ", "", comparison)}
    if (remove.equal == TRUE) {comparison = gsub("=", "", comparison)}
    if (remove.zero == TRUE) {comparison = gsub("0", "", comparison)}
    if (swap.colon == TRUE) {comparison = gsub(":", "-", comparison)}
    if (swap.vs == TRUE) {comparison = gsub("vs", "-", comparison)}
    names(Comparison) = comparison
    if (print.comp == TRUE) {
      Y = data.frame(Comparisons = names(Comparison), p.value = p.value,
                     Value = Comparison, Threshold = threshold)
      cat("\n", "\n")
      print(Y)
      cat("\n", "\n")}
    MCL = multcompLetters(Comparison)
    Group = names(MCL$Letters)
    Letter = as.character(MCL$Letters)
    MonoLetter = as.character(MCL$monospacedLetters)
    Z = data.frame(Group, Letter, MonoLetter)
    return(Z)}
  resp1=resp
  names(resp1)=trat
  postos=rank(resp1)
  somaposto=tapply(postos,names(postos),sum)
  N=tapply(postos,names(postos), length)
  postosmedios=somaposto/N
  media=tapply(resp1,trat,mean)
  mediana=tapply(resp1,trat,median)
  dunns=cldList(P.adj~Comparison,
          data=res)
  tabela=data.frame("group"=dunns$Group,
             "Sum Rank"=somaposto,
             "Mean Rank"=postosmedios,
             "Mean"=media,
             "Median"=mediana,
             "dunn"=dunns$Letter)
  krusk=kruskal.test(resp,trat,method=method)
  chi=krusk$statistic
  pvalor=krusk$p.value
  list("Statistic"=chi,
       "p-value"=pvalor,
       "Post-hoc"=tabela)}
