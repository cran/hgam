### Gruppe: Rest

### plot(test.m)
### summary(test.m) -> summary.hgam
### liefert zurueck Objekt der Klasse summary.hgam
### Uebergabe an print.summary.hgam
### predict(test.m, newdata = ...)
### resid(test.m)
### coef(test.m)
### fitted(test.m)
### logLik(test.m)



coef.hgam <- function(object, ...)
    return(object$coef)

print.hgam <- function(x, ...) {
    cat("\n")
    cat("High-dimensional additive model\n")
    cat("\n")
    cat("Coefficients:\n")
    print(coef(x))
    cat("\n")
    invisible(x)
}


fitted.hgam <- function(object, ...){
    predict(object, ...)
}


### link? family?
residuals.hgam <- function(object, ...) {
    return(object$y - fitted(object))
}

summary.hgam <- function(object,...){
    M <- as.matrix(residuals(object))
    SSE <- t(M) %*% M
    p <- length(object$index) - 1
    n <- length(object$y)
    MSE <- SSE/(n-(p+1))
    coef <- coef(object)
    var_beta <- diag(solve(crossprod(object$Btilde)) * as.numeric(MSE))
   	tvalue <- coef/sqrt(var_beta)


    coef.table <- cbind(coef, sqrt(var_beta), tvalue)
    dimnames(coef.table) <- list(NULL, c("Estimate", "Std.Err",
                                    "t value"))

    ret <- list(coef = coef, SSE = SSE, MSE = MSE, Residual_standard_error = sqrt(MSE))
### class
    class(ret) <- "summary.hgam"
    return(ret)

}

print.summary.hgam <- function(object, ...) {

    printCoefmat(data.frame(coef = object$coef), 
                 P.values=FALSE, has.Pvalue=FALSE)
}


### auch fuer newdata = NULL -> Lerndaten benutzen
### zulassen, dass newdata nur Teilmengen der Kovariablen
### beinhaltet
predict.hgam <- function(object, newdata = NULL, which = NULL,
                         intercept = NULL, ...)  {
    # Koeffizienten des uebergebenen hgam-Objekts herauslesen
    cf <- coef(object)

    # Falls intercept nicht vom Benutzer festgelegt wird, dann soll
    # auf gleiche Weise verfahren werden, wie bei den Lerndaten
    if(is.null(intercept) && rownames(cf)[1] != "Intercept") intercept <- FALSE
    else if(is.null(intercept) && rownames(cf)[1] == "Intercept") intercept <- TRUE

    # Falls intercept-Berucksichtigung vom Benutzer gewuenscht wird, dann muss
    # auch ein Intercept-Koeffizient berechnet worden sein, ansonsten wird
    # eine entsprechende Warnung ausgegeben und intercept auf FALSE gesetzt
    if(rownames(cf)[1] != "Intercept" && intercept){
        warning("Intercept coef missing, setting intercept to FALSE",
                call. = FALSE, immediate. = TRUE)
        intercept <- FALSE
    }

    # index-Vektor aus der Designmatrix herauslesen
    index <- attr(object$Btilde, "assign")

    # Falls which-Argument nicht angegeben wird, dann nur diese Praediktoren
    # beruecksichtigen, bei denen die Koeffizienten groesser Null sind
    index2 <- unique(index) ### unique(index[abs(cf) > 0])
    if (is.null(which)) {
        which <- index2
        which <- which[which > 0]
    } else {
        intercept <- FALSE
    }

    # Designmatrix besteht entweder aus den Lerndaten (falls keine neue Daten
    # uebergeben werden), oder aus den neuen Daten, wobei durch Benutzung der
    # "alten" Knoten die Btilde-Matrix neu berechnet wird
    if(is.null(newdata)) {
        X <- object$Btildenew(which = which, intercept = intercept)
    } else {
        X <- object$Btildenew(xnew = newdata, which = which, intercept = intercept)
    }

    # Der neue index wird nun aus dem "neuen" X herausgelesen
    indexnew <- attr(X, "assign")

    # Nur diese Koeffizienten benutzen, die zum neuen index passen
    cf <- cf[index %in% indexnew, , drop = FALSE]

    # Rueckgabe von Designmatrix * Koeffizienten
    X %*% cf
}

plot.hgam <- function(x, which = NULL, newdata = NULL, 
                      rug = TRUE, multidim = FALSE, ...) {
  
  object <- x
  
  
  if (is.null(which)) {
      which <- sort(unique(attr(object$Btilde, "assign")[abs(cf) > 0]))
      which <- which[which != 0]
  }
  
  if (multidim && length(which) == 2) {
    
  cf <- coef(object)

  if (is.null(newdata)) {
      x <- object$x
  } else {
      x <- newdata
  }
  args <- list(...)
  args$xlab <- x$name
  xlab <- args$xlab

  x1 <- x[[which[1]]]
  x2 <- x[[which[2]]]

  x1 <- seq(from = min(x1), to = max(x1), length = 100)
  x2 <- seq(from = min(x2), to = max(x2), length = 100)

  nd <- expand.grid(x1, x2)
  names(nd) <- colnames(x[, which])
  pred <- predict(object, which = which, newdata = nd)
  z <- matrix(pred, nrow = length(x1), ncol = length(x2))

#  if (is.null(args$ylim))
#      args$ylim <- range(sapply(tmp, function(x) range(x$f)))
  if (is.null(args$zlab))
      args$zlab <- "Marginal effect"
#  if (is.null(args$type))
#      args$type <- "l"

      args$theta = 30
      args$phi = 30
      args$expand = 0.5
      args$ticktype = "detailed"

      args$x <- x1
      args$y <- x2
      args$z <- z
      args$xlab <- names(x)[which[1]]
      args$ylab <- names(x)[which[2]]
      
      do.call("persp", args)
      #if (rug) rug(i$x)

  invisible(NULL)
  
  }
  else{
  cf <- coef(object)
  if (is.null(which)) {
      which <- sort(unique(attr(object$Btilde, "assign")[abs(cf) > 0]))
      which <- which[which != 0]
  }
  if (is.null(newdata)) {
      x <- object$x
  } else {
      x <- newdata
  }
  args <- list(...)
  xlab <- args$xlab
  tmp <- lapply(which, function(w) {
      xw <- x[, w, drop = TRUE]
      ox <- order(xw)
      fw <- predict(object, which = w, newdata = x)
      list(x = xw[ox], f = fw[ox], name = 
           ifelse(is.null(xlab), names(x)[w], xlab[w]), w = w)
  })

  if (is.null(args$ylim)) 
      args$ylim <- range(sapply(tmp, function(x) range(x$f)))
  if (is.null(args$ylab)) 
      args$ylab <- "Marginal effect"
  if (is.null(args$type)) 
      args$type <- "l"
  lapply(tmp, function(x) {
      args$x <- x$x
      args$y <- x$f
      args$xlab <- x$name
      do.call("plot", args)
      if (rug) rug(x$x)
  })
  invisible(NULL)
  }
}


plot.hrisk <- function(x, type = c("levelplot", "3d"), ...) {

  risk <- x$risk[[3]]
  lambda2 <- x$risk[[1]]
  lambda1 <- x$risk[[2]]

  type <- match.arg(type)
  
  if (type == "3d") {
    plot3d(lambda1, lambda2, risk, col="red") }
    else {
      levelplot(risk~lambda1*lambda2) }
}

# model.response
model.response.hgam <- function(object, ...) return(object$y)
# model.weights
model.weights.hgam <- function(object, ...) return(object$weights)

# logLik
logLik.hgam <- function(object, ...){
    y <- model.response.hgam(object)
    f <- fitted(object)
    w <- model.weights.hgam(object)
    ret <- -object$model@nloglik(y, f, w)
    class(ret) <- "logLik.hgam"
    return(ret)
}

print.logLik.hgam <- function(object){
    cat("'log Lik.'", object, "\n")
    invisible(object)
}

print.hrisk <- function(h, ...) {
    cat("\n")
    cat("    Best Lambda with method", h$type ,"\n")
    cat("\n")
    cat("    Lambda1:", h$lambda1, "\n")
    cat("    Lambda2:", h$lambda2, "\n")
    cat("\n")
    invisible(h)
}