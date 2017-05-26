#-------------------------------------------------------------------------------
# Cross-correlation Actin Toca: MODULE 2: CCF computation

#-------------------------------------------------------------------------------
# Useful functions:

Count <- function(x) length(x[!is.na(x)])			 
SE <- function(x) sd(x, na.rm=TRUE)/sqrt(Count(x))	 							
# CI <- function(x) 1.96*sd(x, na.rm=TRUE)/sqrt(Count(x))     

CI <- function (x, ci = 0.95) { 							
  
  # Using T distribution; 
  # Appropriate for small samples
  # ref: Rmisc package
  # Ryan M. Hope (2013). Rmisc: Rmisc: Ryan Miscellaneous. R package version 1.5. https://CRAN.R-project.org/package=Rmisc
  # citation("Rmisc")
  
  # My adjustments from Rmisc::CI:
  # 1. output as error value (not range around mean)
  # 2. n as Count(x) instead of length(x), so that NA values do not contribute to n
  
  a <- mean(x, na.rm = TRUE)
  s <- sd(x, na.rm = TRUE)
  n <- Count(x)
  error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
  return(error)
}

DrawErrorAsPolygon <- function(x, y1, y2, tt, col = 'grey') {
  polygon(c(x[tt], rev(x[tt])), c(y1[tt], rev(y2[tt])), 
          col = col,
          border = NA)			
}

MovingAverage <- function(x, w = 5) {
  filter(x, rep(1/w, w), sides = 2)
}
#-------------------------------------------------------------------------------

n.slices 
max.lag = 10

all.lag = (-max.lag:max.lag); all.lag

ccf.toca.actin <- data.frame(matrix(NA, ncol = ncol(w.actin), nrow = 2*max.lag + 1))
ccf.toca.actin.bg <- data.frame(matrix(NA, ncol = ncol(w.actin.bg), nrow = 2*max.lag + 1))
row.names(ccf.toca.actin) <- all.lag
row.names(ccf.toca.actin.bg) <- all.lag

for(i in 1:n.spots) {
  ccf.i  <- ccf(w.toca[, i], w.actin[, i], lag.max = max.lag, na.action = na.pass, plot = FALSE) 
  ccf.i.bg <- ccf(w.toca.bg[, i], w.actin.bg[, i], lag.max = max.lag, na.action = na.pass, plot = FALSE)
  ccf.toca.actin[, i] <- ccf.i
  ccf.toca.actin.bg[, i] <- ccf.i.bg
  rm(ccf.i)
}

ccf.toca.actin
ccf.toca.actin.bg

mean.ccf.toca.actin <- apply(ccf.toca.actin, 1, mean)
mean.ccf.toca.actin.bg <- apply(ccf.toca.actin.bg, 1, mean)

ci.ccf.toca.actin  <- apply(ccf.toca.actin, 1, CI)
ci.ccf.toca.actin.bg  <- apply(ccf.toca.actin.bg, 1, CI)
#ci1 <- mean.ccf.toca.actin + ci.ccf.toca.actin
#ci2 <- mean.ccf.toca.actin - ci.ccf.toca.actin

graphics.off()

matplot(all.lag, ccf.toca.actin, type = "l", lty = 1, 
        col = "#88888850", 
        ylab = "CCF (Toca & actin)",
        ylim=c(-0.6, 1),
        xlab = "Lag [timepoints]"
        )
lines(all.lag, mean.ccf.toca.actin, add = TRUE, lwd = 5, col = "red")
# DrawErrorAsPolygon(x = all.lag, y1 = ci1, y2 = ci2, tt = all.lag, col = "#FF000030")

abline(v = 0, lty = 2)
abline(h = 0, lty  = 1)

text(-10, 0.9, pos = 4, "CCF (individual ROI)", col = "grey")
text(-10, 0.75, pos = 4, "CCF (mean)", col = "red")

text(-10, -0.55, pos = 4, "Toca -> actin")
text(10, -0.55, pos = 2, "actin -> Toca")


#  --- repeat plotting with bg-corrected data:

matplot(all.lag, ccf.toca.actin.bg, type = "l", lty = 1, 
        col = "#88888850", 
        ylab = "CCF (Toca & actin)",
        ylim=c(-0.6, 1),
        xlab = "Lag [timepoints]"
)
lines(all.lag, mean.ccf.toca.actin.bg, add = TRUE, lwd = 5, col = "red")
# DrawErrorAsPolygon(x = all.lag, y1 = ci1, y2 = ci2, tt = all.lag, col = "#FF000030")

abline(v = 0, lty = 2)
abline(h = 0, lty  = 1)

text(-10, 0.9, pos = 4, "CCF (individual ROI)", col = "grey")
text(-10, 0.75, pos = 4, "CCF (mean)", col = "red")

text(-10, -0.55, pos = 4, "Toca -> actin")
text(10, -0.55, pos = 2, "actin -> Toca")


# 
# 
# # ---------------------------------------------------------------------------
# # 6. Calculate CCFs from tip F and tip movement tables
# 
# maxlag = 20
# lag.range <- -maxlag:maxlag
# lag.in.s  <- lag.range * spt
# 
# ccf.tip.dctm <- data.frame(matrix(NA, ncol = ncol(all.move), nrow = 2*maxlag + 1))
# all.filo  <- seq_along(colnames(all.move))
# 
# 
# for (i in all.filo) {
#   ccf.i  <- ccf(tip.f[, i], all.move[, i], lag.max = 20, na.action = na.pass, plot = FALSE) 
#   ccf.tip.dctm[, i] <- ccf.i
#   rm(ccf.i, ccf.z.i)
# }
# 
# colnames(ccf.tip.dctm) <- colnames(all.move)
# row.names(ccf.tip.dctm)  <- lag.in.s