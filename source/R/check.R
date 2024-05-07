check_int <- function(x, nm, min=1) {

  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be an integer"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be an integer"))
  if (!is.finite(x)) stop(paste0("ERROR: ", nm, " must be an integer"))
  x <- as.integer(x)
  if (length(min) && (x < min)) {
    stop(paste0("ERROR: ", nm, " must be >= ", min))
  }
  x
}

check_tol <- function(x, nm="tol") {

  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be a numeric value in (0, 1)"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a numeric value in (0, 1)"))
  if (!is.finite(x)) stop(paste0("ERROR: ", nm, " must be a numeric value in (0, 1)"))
  if ((x <= 0) || (x >= 1)) {
    stop(paste0("ERROR: ", nm, " must be between 0 and 1"))
  }
  NULL
}

check_sig <- function(x, nm="sig.select") {

  if (length(x) != 1) stop(paste0("ERROR: ", nm, " must be a numeric value between 0 and 1"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a numeric value between 0 and 1"))
  if (!is.finite(x)) stop(paste0("ERROR: ", nm, " must be a numeric value between 0 and 1"))
  if ((x < 0) || (x > 1)) {
    stop(paste0("ERROR: ", nm, " must be between 0 and 1"))
  }
  NULL
}

check_mat <- function(x, nm) {

  if (!is.matrix(x)) stop(paste0("ERROR: ", nm, " must be a matrix"))
  if (!is.numeric(x)) stop(paste0("ERROR: ", nm, " must be a numeric"))
  if (any(!is.finite(x))) {
    stop(paste0("ERROR: ", nm, " contains non-finite data"))
  }
  if (nrow(x) != 96) stop(paste0("ERROR: ", nm, " must have 96 rows"))
  if (ncol(x) < 1) stop(paste0("ERROR: ", nm, " has no columns"))
  NULL
}

isString <- function(x) {
  (length(x) == 1) && is.character(x)
}

check_Sig.W <- function(x, nm="Sig.W") {

  if (!length(x)) stop(paste0("ERROR: ", nm, " must be the string 'SBS' or a matrix"))
  if (isString(x)) {
    if (toupper(trimws(x)) != "SBS") {
      stop(paste0("ERROR: ", nm, " must be the string 'SBS' or a matrix"))
    }
  } else if (is.matrix(x)) {
    check_mat(x, nm)
  } else {
    stop(paste0("ERROR: ", nm, " must be the string 'SBS' or a matrix"))
  }
  NULL
}

check_BSR.object <- function(x, nm="BSR.obj") {

  if (!is.list(x)) stop(paste0("ERROR: ", nm, " must be a list returned from BayesSigRef"))
  nms   <- names(x)
  req   <- c("objects")
  tmp   <- !(req %in% nms)
  if (any(tmp)) { 
    miss <- paste0("'", req[tmp], "'")
    str  <- paste0(miss, collapse=", ")
    msg  <- paste0("ERROR: ", nm, " is missing object ", str)
    stop(msg)
  }
  NULL
}



