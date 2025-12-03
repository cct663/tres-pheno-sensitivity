# This is a very slightly modified version of the slidingwin function from climwin
# that I changed to give me more control over what windows are included.


slidingwin_custom <- function (exclude = NA, xvar, cdate, bdate, baseline, type, refday, 
          stat = "mean", func = "lin", range, cmissing = FALSE, cinterval = "day", 
          k = 0, upper = NA, lower = NA, binary = FALSE, centre = list(NULL, 
                                                                       "both"), spatial = NULL, cohort = NULL) 
{
  if (getOption("scipen") < 0) {
    current_option <- getOption("scipen")
    options(scipen = 0)
  }
  if (cmissing != FALSE && cmissing != "method1" && cmissing != 
      "method2") {
    stop("cmissing must be FALSE, 'method1' or 'method2'.")
  }
  if (type != "absolute" && type != "relative") {
    stop("type must be either absolute or relative.")
  }
  if (is.null(cohort) == TRUE) {
    cohort = lubridate::year(as.Date(bdate, format = "%d/%m/%Y"))
  }
  if (k > 0 && class(baseline)[length(class(baseline))] == 
      "coxph") {
    stop("Sorry, cross-validation is not available yet for coxph models")
  }
  if (attr(baseline, "class")[1] == "lme" && k > 0) {
    stop("Sorry, cross-validation is currently not functioning for nlme models. Consider using lme4 if possible.")
  }
  if (is.null(spatial) & length(unique(cdate)) < length(cdate)) {
    stop("Your cdate variable has repeated date measures. Do you have climate data from multiple sites? If so, you should specify the parameter `spatial`.")
  }
  if (is.null(centre[[1]]) == FALSE) {
    func = "centre"
  }
  if (is.list(xvar) == FALSE) {
    stop("xvar should be an object of type list")
  }
  if (is.null(names(xvar)) == TRUE) {
    numbers <- seq(1, length(xvar), 1)
    for (xname in 1:length(xvar)) {
      names(xvar)[xname] = paste("climate", numbers[xname])
    }
  }
  if (is.na(upper) == FALSE && is.na(lower) == FALSE) {
    combos <- expand.grid(list(upper = upper, lower = lower))
    combos <- combos[which(combos$upper >= combos$lower), 
    ]
    allcombos <- expand.grid(list(climate = names(xvar), 
                                  type = type, stat = stat, func = func, gg = c(1:nrow(combos)), 
                                  binary = binary))
    allcombos <- cbind(allcombos, combos[allcombos$gg, ], 
                       deparse.level = 2)
    binarylevel <- "two"
    allcombos$gg <- NULL
  }
  else if (is.na(upper) == FALSE && is.na(lower) == TRUE) {
    allcombos <- expand.grid(list(climate = names(xvar), 
                                  type = type, stat = stat, func = func, upper = upper, 
                                  lower = lower, binary = binary))
    binarylevel <- "upper"
  }
  else if (is.na(upper) == TRUE && is.na(lower) == FALSE) {
    allcombos <- expand.grid(list(climate = names(xvar), 
                                  type = type, stat = stat, func = func, upper = upper, 
                                  lower = lower, binary = binary))
    binarylevel <- "lower"
  }
  else if (is.na(upper) == TRUE && is.na(lower) == TRUE) {
    allcombos <- expand.grid(list(climate = names(xvar), 
                                  type = type, stat = stat, func = func))
    binarylevel <- "none"
  }
  rownames(allcombos) <- seq(1, nrow(allcombos), 1)
  combined <- list()
  for (combo in 1:nrow(allcombos)) {
    runs <- basewin_custom(exclude = exclude, xvar = xvar[[paste(allcombos[combo, 
                                                                    1])]], cdate = cdate, bdate = bdate, baseline = baseline, 
                    range = range, type = paste(allcombos[combo, 2]), 
                    refday = refday, stat = paste(allcombos[combo, 3]), 
                    func = paste(allcombos[combo, 4]), cmissing = cmissing, 
                    cinterval = cinterval, k = k, upper = ifelse(binarylevel == 
                                                                   "two" || binarylevel == "upper", allcombos$upper[combo], 
                                                                 NA), lower = ifelse(binarylevel == "two" || binarylevel == 
                                                                                       "lower", allcombos$lower[combo], NA), binary = paste(allcombos$binary[combo]), 
                    centre = centre, cohort = cohort, spatial = spatial)
    combined[[combo]] <- runs
    allcombos$DeltaAICc[combo] <- round(runs$Dataset$deltaAICc[1], 
                                        digits = 2)
    allcombos$WindowOpen[combo] <- runs$Dataset$WindowOpen[1]
    allcombos$WindowClose[combo] <- runs$Dataset$WindowClose[1]
    if (any(grepl("climate", model.frame(baseline)))) {
      if (length(which("lin" == levels(allcombos$func))) > 
          0) {
        allcombos$betaL[combo] <- round(runs$Dataset$ModelBeta[1], 
                                        digits = 2)
      }
      if (allcombos$func[1] == "centre") {
        if (centre[[2]] == "both") {
          allcombos$WithinGrpMean <- round(runs$Dataset$WithinGrpMean[1], 
                                           digits = 2)
          allcombos$WithinGrpDev <- round(runs$Dataset$WithinGrpDev[1], 
                                          digits = 2)
        }
        if (centre[[2]] == "dev") {
          allcombos$WithinGrpDev <- round(runs$Dataset$WithinGrpDev[1], 
                                          digits = 2)
        }
        if (centre[[2]] == "mean") {
          allcombos$WithinGrpMean <- round(runs$Dataset$WithinGrpMean[1], 
                                           digits = 2)
        }
      }
      if (length(which("quad" == levels(allcombos$func))) > 
          0) {
        allcombos$betaL[combo] <- round(runs$Dataset$ModelBeta[1], 
                                        digits = 2)
        allcombos$betaQ[combo] <- round(runs$Dataset$ModelBetaQ[1], 
                                        digits = 2)
      }
      if (length(which("cub" == levels(allcombos$func))) > 
          0) {
        allcombos$betaL[combo] <- round(runs$Dataset$ModelBeta[1], 
                                        digits = 2)
        allcombos$betaQ[combo] <- round(runs$Dataset$ModelBetaQ[1], 
                                        digits = 2)
        allcombos$betaC[combo] <- round(runs$Dataset$ModelBetaC[1], 
                                        digits = 2)
      }
      if (length(which("inv" == levels(allcombos$func))) > 
          0) {
        allcombos$betaInv[combo] <- round(runs$Dataset$ModelBeta[1], 
                                          digits = 2)
      }
      if (length(which("log" == levels(allcombos$func))) > 
          0) {
        allcombos$betaLog[combo] <- round(runs$Dataset$ModelBeta[1], 
                                          digits = 2)
      }
    }
  }
  allcombos <- cbind(response = colnames(model.frame(baseline))[1], 
                     allcombos)
  combined <- c(combined, combos = list(allcombos))
  if (exists("current_option")) {
    options(scipen = current_option)
  }
  return(combined)
}
