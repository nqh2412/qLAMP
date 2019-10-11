
#pcrfit from "qpcR"
plot.pcrfit <- function(
  x, 
  type = c("all", "single", "3D", "image"),
  fitted = TRUE, 
  add = FALSE,
  col = NULL, 
  par2D = list(),
  par3D = list(),
  ...) 
{
  type <- match.arg(type)    
  object <- x
  
  print(class(x))
  
  if (class(x)[1] != "modlist") modLIST <- list(object) else modLIST <- object      
  
  ## extract cycles and fluorescence values from all curves
  allCYC <- lapply(modLIST, function(x) x$DATA[, 1])
  allFLUO <- lapply(modLIST, function(x) x$DATA[, 2])
  vecCYC <- do.call(c, allCYC)
  vecFLUO <- do.call(c, allFLUO)
  
  ## make unique cycles  
  CYC <- unique(as.numeric(vecCYC))  
  CYC <- CYC[!is.na(CYC)]    
  
  ## calculate min and max fluo values for defining ylim 
  MIN <- min(vecFLUO, na.rm = TRUE)   
  MAX <- max(vecFLUO, na.rm = TRUE)     
  
  ## length of 'modlist'
  LEN <- length(modLIST)
  ## names of 'modlist'
  NAMES <- sapply(modLIST, function(x) x$names)   
  
  ## define plotting colors
  if (is.null(col)) {
    COL <- rep(1, LEN)    
    if (class(object)[2] == "replist") COL <- rainbow(attr(object, "nlevels"))     
  } else COL <- rep(col, length.out = LEN)   
  
  ## 3D plot empty setup using par3D parameters
  if (type == "3D") {
    do.call(plot3d, modifyList(list(x = CYC, y = 1:LEN, z = MAX, type = "n", axes = FALSE, box = FALSE, xlab = "", 
                                    ylab = "", zlab = "", zlim = c(0, 1.1 * MAX)), par3D))
    do.call(axis3d, modifyList(list('x', at = pretty(CYC), cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Minute", 'x', line = 2), par3D))     
    do.call(axis3d, modifyList(list('y', at = 1:LEN, label = NAMES, cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Run", 'y', line = 2), par3D))
    do.call(axis3d, modifyList(list('z', cex = 0.5), par3D))
    do.call(mtext3d, modifyList(list("Hue", 'z', line = 2), par3D))
  }   
  
  ## standard 'all' plot empty setup
  if (type == "all" && !add) {   
    tempLIST <- modifyList(list(CYC, rep(MAX, length(CYC)), ylim = c(MIN, MAX), 
                                xlab = "Minute", ylab = "Hue", las = 1), par2D)
    tempLIST$type <- "n"
    do.call(plot, tempLIST)   
  }
  
  ## plot matrix empty setup
  if (type == "single") {
    DIM <- ceiling(sqrt(LEN))
  } 
  
  ## image plot 
  if (type == "image") {
    RUNS <- 1:length(modLIST)
    nRUNS <- length(RUNS)
    ## unique cycles
    CYCS <- unique(unlist(lapply(modLIST, function(x) x$DATA[, 1])))
    nCYCS <- length(CYCS)
    ## convert list with fluo data to matrix, fll unequal length with NA
    allLIST <- lapply(modLIST, function(x) x$DATA[, 2])
    maxCYCS <- max(sapply(allLIST, length))
    for (i in 1:length(allLIST)) allLIST[[i]] <- c(allLIST[[i]], rep(NA, maxCYCS - length(allLIST[[i]])))
    allDAT <- do.call(cbind, allLIST)
    ## image setup
    allDAT <- allDAT[, ncol(allDAT):1]
    image(allDAT, col = heat.colors(100), axes = FALSE, xlab = "Minutes", ylab = "Runs")
    axis(1, at = seq(0, 1, length.out = nCYCS), labels = CYCS)
    axis(2, at = seq(0, 1, length.out = nRUNS), labels = rev(RUNS))
  }
  
  ## iterate through all curves
  for (i in 1:LEN) {
    DATA <- modLIST[[i]]$DATA    
    DATA <- na.omit(DATA)      
    FITTED <- fitted(modLIST[[i]])       
    m <- match(CYC, DATA[, 1])
    m <- na.omit(m)
    
    ## plot 3D curves
    if (type == "3D") {
      do.call(points3d, modifyList(list(x = DATA[, 1], y = i, z = DATA[, 2], color = COL[i]), par3D))
      if (!is.null(FITTED) && fitted) do.call(lines3d, modifyList(list(x = DATA[m, 1], y = i, z = FITTED[m], color = COL[i]), par3D))      
    }
    
    ## plot 2D curves
    if (type == "all") {
      do.call(points, modifyList(list(DATA[, 1], DATA[, 2], col = COL[i]), par2D))
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[m, 1], FITTED[m], col = COL[i]), par2D)) 
    } 
    
    ## plot matrix curves
    if (type == "single") {
      NAME <- NAMES[i]
      ## color by failed fit or failed structure
      if (grepl("\\*\\*[[:alnum:]]*", NAME)) colMAIN <- "blue" 
      else if (grepl("\\*[[:alnum:]]*", NAME)) colMAIN <- "red"
      else colMAIN <- "black"
      TRY <- try(do.call(plot, modifyList(list(DATA[, 1], DATA[, 2], main = NAME, cex.main = 0.7, col.main = colMAIN, type = "p", 
                                               xlab = FALSE, ylab = FALSE, xaxt = "n", yaxt = "n", col = COL[i]), par2D)), silent = TRUE)
      if (inherits(TRY, "try-error")) next      
      if (!is.null(FITTED) && fitted) do.call(lines, modifyList(list(DATA[m, 1], FITTED[m], col = COL[i]), par2D))      
    }     
  }     
}  

#Calib from "qpcR"
calib <- function(
  refcurve, 
  predcurve = NULL, 
  thresh = "cpD2", 
  dil = NULL,
  group = NULL,
  plot = TRUE,
  conf = 0.95,
  B = 200
)
{
  if (class(refcurve)[1] != "modlist") stop("'refcurve' is not a 'modlist'!")
  if (!is.null(predcurve) & class(predcurve)[1] != "modlist") stop("'predcurve' is not a 'modlist'!")
  if (thresh != "cpD2" && !is.numeric(thresh)) stop("'thresh' must be either 'cpD2' or numeric!")
  if (is.null(dil)) stop("Please define dilutions!")
  if (!is.null(group) && (length(dil) != length(unique(group)))) stop("Supply as many dilutions as number of PCR groups in 'refcurve'!")
  
  lref <- length(refcurve)
  lpred <- length(predcurve)
  lgroup <- length(unique(group))
  dil <- log10(dil)
  COLref <- rep(rainbow(nlevels(as.factor(dil))), table(as.factor(dil)))
  COLpred <- rep(rainbow(lpred))   
  
  if(is.null(group))  {
    group <- as.factor(1:lref)
    isReps <- FALSE
  } else isReps <- TRUE
  
  LMFCT <- function(dil, ref, pred = NULL, conf) {
    linModY <- lm(ref ~ dil)
    conf.Y <- predict(linModY, interval = "confidence", level = conf)
    eff <- as.numeric(10^(-1/coef(linModY)[2]))
    FOM1 <- AIC(linModY)
    FOM2 <- AICc(linModY)
    FOM3 <- Rsq(linModY)
    FOM4 <- Rsq.ad(linModY)
    
    if (!is.null(pred)) {
      linModX <- lm(dil ~ ref)
      pred.conc <- sapply(as.numeric(pred), function(x) predict(linModX, newdata = data.frame(ref = x), interval = "confidence", level = conf))
    } else pred.conc <- NULL
    
    return(list(linModY = linModY, conf.Y = conf.Y, eff = eff, FOM1 = FOM1, FOM2 = FOM2,
                FOM3 = FOM3, FOM4 = FOM4, pred.conc = pred.conc[1, ], pred.conf = pred.conc[2:3, ]))
  }
  
  print("Calculating threshold time of reference curves...")
  flush.console()
  
  if (thresh == "cpD2") refCt <- sapply(refcurve, function(x) efficiency(x, plot = FALSE)$cpD2)
  else refCt <- as.numeric(sapply(refcurve, function(x) predict(x, newdata = data.frame(Fluo = thresh), which = "x")))   
  
  print("Calculating threshold time of prediction curves...")
  flush.console()
  
  if (!is.null(predcurve)) {
    if (thresh == "cpD2") predCt <- sapply(predcurve, function(x) efficiency(x, plot = FALSE)$cpD2)
    else predCt <- as.numeric(sapply(predcurve, function(x) predict(x, newdata = data.frame(Fluo = thresh), which = "x")))
  } else predCt <- NULL
  
  iterRef <- split(refCt, group)
  
  lmResList <- list()
  iterMat <- matrix(ncol = lgroup, nrow = B)
  
  for (i in 1:B) {
    if (isReps) selRef <- sapply(iterRef, function(x) sample(x, 1))
    else selRef <- unlist(iterRef)
    lmRes <- LMFCT(dil = dil, ref = as.numeric(selRef), pred = predCt, conf = conf)
    lmResList[[i]] <- lmRes
    iterMat[i, ] <- selRef
    
    if (plot) {
      if (i == 1) {
        plot(dil, selRef, col = c(2:7), pch = 16, cex = 1.3, xlab = "Log of DNA Copy Number", ylab = "Threshold Time", main = "Standard Calibration Curve", add = FALSE)
      } else {
        points(dil, selRef, col = c(2:7), pch = 16, cex = 1.3)
        abline(lmRes$linModY, lwd = 2)
        #Red line in the threshold graph
        lines(dil, lmRes$conf.Y[, 2], col = 2, lty = 3)
        lines(dil, lmRes$conf.Y[, 3], col = 2, lty = 3)
        
        #show the regression equation
        cf <- round(coef(lmRes$linModY), 4)
        r_sqr <- summary(lmRes$linModY)$r.squared
        thresh_data_table <- cbind.data.frame(dil, selRef)
        eq <- paste0("y = ", cf[1],
                     ifelse(sign(cf[2])==1, " + ", " - "), abs(cf[2]), "x", ". R.square = ", 
                     round(r_sqr, 4))
        
        ## printing of the equation
        mtext(eq, 3, line=0)
      }
      if (!is.null(predcurve)) {
        points(lmRes$pred.conc, predCt, pch = 15, col = COLpred, cex = 1.5)
        if (is.vector(lmRes$pred.conf)) lmRes$pred.conf <- matrix(lmRes$pred.conf, ncol = 1)
        if (!all(is.na(lmRes$pred.conc))) {
          arrows(lmRes$pred.conf[1, ], predCt, lmRes$pred.conf[2, ], predCt, code = 3, angle = 90, length = 0.1, col = "blue")
        }
      }
    }
  }
  
  summaryList <- list()
  lenRML <- 2:length(lmResList[[1]])
  
  for (i in lenRML) {
    temp <- sapply(lmResList, function(x) x[[i]])
    summaryList[[i - 1]] <- t(temp)
  }
  
  names(summaryList) <- names(lmRes[lenRML])
  
  alpha = 1 - conf
  CONFINT <- function(x, alpha = alpha) quantile(x, c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)
  
  CONF.eff <- CONFINT(summaryList$eff, alpha = alpha)
  CONF.AICc <- CONFINT(summaryList$FOM2, alpha = alpha)
  CONF.Rsq.ad <- CONFINT(summaryList$FOM4, alpha = alpha)
  
  if (!is.null(predcurve)) {
    if (nrow(summaryList$pred.conc) == 1) summaryList$pred.conc <- t(summaryList$pred.conc)
    CONF.predconc <- apply(summaryList$pred.conc, 2, function(x) CONFINT(x, alpha = alpha))
    if (!isReps) CONF.predconc <- apply(rbind(lmRes$pred.conf[1, ], lmRes$pred.conf[2, ]) , 2, function(x) CONFINT(x, alpha = alpha))
  } else {
    summaryList$pred.conc <- NULL 
    CONF.predconc <- NULL
  } 
  
  #if (plot) {
    #boxplot(as.numeric(summaryList$eff), main = "Efficiency", cex = 0.2)
    #abline(h = CONF.eff, col = 2, lwd = 2)
    #boxplot(as.numeric(summaryList$FOM2), main = "corrected AIC", cex = 0.2)
    #abline(h = CONF.AICc, col = 2, lwd = 2)
    #boxplot(as.numeric(summaryList$FOM4), main = "adjusted R-square", cex = 0.2)
    #abline(h = CONF.Rsq.ad, col = 2, lwd = 2)
    #if (!is.null(predcurve)) {
      #boxplot(summaryList$pred.conc, main = "log(conc) of predicted", cex = 0.2)
      #abline(h = CONF.predconc, col = 2, lwd = 2)
    #}
  #}
  return(list(eff = summaryList$eff, AICc = summaryList$FOM2, Rsq.ad = summaryList$FOM4, predconc = summaryList$pred.conc,
              conf.boot = list(conf.eff = CONF.eff, conf.AICc = CONF.AICc, conf.Rsq.ad = CONF.Rsq.ad, conf.predconc = CONF.predconc), thresh_data_table))
}
