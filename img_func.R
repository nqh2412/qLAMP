
#crop
img_crop <- function(img, x_1, x_2, y_1, y_2) {
  c_im1 <- imsub(img, x > width/x_1, y > height/y_1)
  c_im2 <- imsub(c_im1, x < width/x_2, y < height/y_2)
  return(c_im2)
  UseMethod('img_crop')
}
##Modifying build-in function of 'colocr'
roi_select <- function(img, threshold, shrink = 5, grow = 5, fill = 5,
                       clean = 5, tolerance = .1, n = 7) {
  UseMethod('roi_select')
}

#' @export
roi_select.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}

#' @export
roi_select.cimg <- function(img, threshold = 50, shrink = 5, grow = 5, fill = 1,
                            clean = 1, tolerance = .1, n = 7) {
  
  # check valid input
  if(!missing(threshold) & !is.numeric(threshold)) {
    stop('threshold should be a numeric >= 0 and < 100.')
  }
  if(!missing(threshold) & (threshold >= 100 | threshold < 0)) {
    stop('threshold should be a numeric >= 0 and < 100.')
  }
  
  # Ideally, I'd like to check type and value of other arguments,
  # however, I currently cannot since these arguments are optional
  #crop and rotate the image
  
  #img <- img_crop(img, x_1 = 3.5, x_2 = 1.65, y_1 = 2.3, y_2 = 5)
  img <- imrotate(img,90)
  # change image to gray scale
  img.g <- grayscale(img)
  
  # apply threshold
  img.t <- threshold(img.g, paste0(threshold, '%'))
  
  # change to pixset
  px <- as.pixset(1-img.t)
  
  # apply shrink
  px.m <- shrink(px, shrink)
  
  # apply grow
  px.m <- grow(px.m, grow)
  
  # apply fill
  px.m <- fill(px.m, fill)
  
  # apply clean
  px.m <- clean(px.m, clean)
  
  # add labels when n is provided
  labs.px <- .labels_add(px.m, tolerance = tolerance, n = n)
  attr(img, 'label') <- as.numeric(labs.px)
  
  # return object
  return(img)
}

#' @export
roi_select.list <- function(img, threshold, shrink = 5, grow = 5, fill = 1,
                            clean = 1, tolerance = .1, n = 7) {
  # get the length of the image list
  img_n <- length(img)
  
  # repeat arguments to match list length
  inputs <- list(threshold = threshold,
                 shrink = shrink,
                 grow = grow)
  
  for(i in seq_along(inputs)) {
    # lenght of argument
    input_n <- length(inputs[[i]])
    
    # use first item and return warning if not a single value or doesn't match
    # length of image list
    if(input_n != img_n & input_n != 1) {
      inputs[[i]] <- inputs[[i]][1]
      warning(paste0("Only first value in ", names(inputs)[[i]], ' will be used.'))
    }
    
    # match length of the arguments to that of the list of images
    if(input_n != img_n) {
      inputs[[i]] <- rep(inputs[[i]], img_n)
    }
  }
  
  # loop over the list of images and call roi_select
  newimgs <- list()
  for(i in 1:img_n) {
    newimgs[[i]] <- roi_select(img[[i]],
                               threshold = inputs$threshold[i],
                               shrink = inputs$shrink[i],
                               grow = inputs$grow[i],
                               fill = 1,
                               clean = 1,
                               tolerance = 0.1,
                               n = 7)
  }
  
  # return list of images
  return(newimgs)
}
#Roi_show modified to show only original and the pix set images with hightlight

roi_show <- function(img, ind = c(1,3)) {
  UseMethod('roi_show')
}


roi_show.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}


roi_show.cimg <- function(img, ind = c(1,3)) {
  
  # get labels from img
  # transform labels to cimg object
  labels <- attr(img, 'label')
  dims <- dim(grayscale(img))
  a <- array(labels, dim = dims)
  
  px <- cimg(a)
  
  # merge image
  
  plot(img,
       axes = FALSE,
       main = 'Input image')
  # pixset image
  plot(px,
       axes = FALSE,
       main = 'Read image')
  highlight(px)
  
  # return null
  invisible(NULL)
}


roi_show.list <- function(img, ind = c(1,3)) {
  
  # get the length of the image list
  img_n <- length(img)
  
  # repeat argument to match list length
  if(!is.list(ind)) {
    ind <- rep(list(ind), img_n)
  }
  
  # loop over the images of lists and call roi_show
  for(i in 1:img_n){
    roi_show(img[[i]],
             ind = ind[[i]])
  }
  
  # return null
  invisible(NULL)
}

#modified function ROI_CHECK
roi_check <- function(img, ind = c(1:3)) {
  UseMethod('roi_check')
}

roi_check.default <- function(img, ...) {
  warning(paste("img is of class",
                class(img),
                ". img should be a cimg or a list of cimg objects."))
}


roi_check.cimg <- function(img, ind = c(1:3)) {
  # get pixel intensities
  label <- attr(img, 'label')
  
  #Red channel intensities of ROIS
  r_channel <- channel(img, ind = 1)
  r_rois <- split(r_channel, label)
  str(r_rois)
  red_int <- data.frame(lapply(r_rois,mean))
  red_int <- red_int[,2:8] + (1-red_int[,1])
  #Green channel intensities of ROIS
  g_channel <- channel(img, ind = 2)
  g_rois <- split(g_channel, label)
  str(g_rois)
  gre_int <- data.frame(lapply(g_rois,mean))
  gre_int <- gre_int[,2:8] + (1-gre_int[,1])
  #Blue channel intensities of ROIS
  b_channel <- channel(img, ind = 3)
  b_rois <- split(b_channel, label)
  str(b_rois)
  blu_int <- data.frame(lapply(b_rois,mean))
  blu_int <- blu_int[,2:8] + (1-blu_int[,1])
  
  img_hsv_value <- rgb2hsv(r = as.numeric(red_int),g= as.numeric(gre_int), b= as.numeric(blu_int), maxColorValue = 1)
  img_h_value <- data.frame(img_hsv_value[1,])
  img_h_value <- data.frame(t(img_h_value))
  #img_h_value[,1] <- NULL
  return(img_h_value)
  # return null
  invisible(NULL)
}


roi_check.list <- function(img, ind = c(1,3)) {
  # get the length of the image list
  img_n <- length(img)
  a <- list()
  # repeat argument to match list length
  if(!is.list(ind)) {
    ind <- rep(list(ind), img_n)
  }
  
  # loop over the images of lists and call roi_check
  for(i in 1:img_n) {
    a[[i]] <-  roi_check(img[[i]],
                         ind = ind[[i]])
  }
  
  #name_list_table <- as.data.frame(bind_rows(name_list, .id = "order"))
  a.table <- as.data.frame(bind_rows(a, .id = "Time"))
  
  a.table$Time <- order(as.numeric(a.table$Time))
  colnames(a.table) <- c('Time', 'NC', "C1","C2","C3","C4","C5","C6")
  return(a.table)
  # return null
  invisible(NULL)
}

rlamp <- function(df,df2, standardcurve, dil2, thr_2, thr_01) {
  ml1 <- modlist(df, 1, 2:8, model = l7)
  plot1 <- plot(ml1, col = rep(1:7, each = 1))
  ml2 <- modlist(df, 1, 3:8, model = l7)
  if (standardcurve == "Standard Calibration Curve"){
    c1 <- calib(ml2, thresh = "cpD2", predcurve = NULL, dil = dil2,
                group = NULL, plot = TRUE, conf = 0.95, B = 200)
    c2 <- round(t(as.data.frame(c1[[6]])),2)
    Data <- c("DNA copy number (log)", "Threshold time")
    colnames(c2) <- c("#2","#3","#4","#5","#6","#7")
    c2<- cbind(Data, c2)
  }
  
  else{
    pred1 <- modlist(df2, 1, 3:8, model = l7)
    c1 <- calib(pred1, thresh = "cpD2", predcurve = ml2, dil = dil2,
                group = NULL, plot = TRUE, conf = 0.95, B = 200)
      c2 <- round(t(as.data.frame(c1$predconc[1,])),2)
    Data <- c("DNA copy number (log)")
    colnames(c2) <- c("#2","#3","#4","#5","#6","#7")
    c2<- cbind(Data, c2)
  }
  return(c2)
}
