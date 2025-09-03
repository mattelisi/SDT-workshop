hablar::set_wd_to_script_path()


simulate_sdt <- function(dp, crit, nS, nN = nS, output = c("data", "table", "both")) {
  
  # input checks
  stopifnot(length(dp) == 1L, is.finite(dp),
            length(crit) == 1L, is.finite(crit),
            length(nS) == 1L, nS >= 0, is.finite(nS),
            length(nN) == 1L, nN >= 0, is.finite(nN))
  nS <- as.integer(nS); nN <- as.integer(nN)
  output <- match.arg(output)
  
  # SDT response probabilities
  pFA <- 1 - pnorm(crit)         # P(Yes | Noise)
  pH  <- 1 - pnorm(crit - dp)    # P(Yes | Signal)
  
  # simulate responses
  resp_S <- rbinom(nS, 1L, pH)   # 1="Yes"
  resp_N <- rbinom(nN, 1L, pFA)
  
  # set up simulated dataset
  df <- rbind(
    data.frame(signal_present = 1L, response = resp_S),
    data.frame(signal_present = 0L, response = resp_N)
  )
  rownames(df) <- NULL
  
  # optional 2x2 table
  H  <- sum(df$response == 1L & df$signal_present == 1L)
  M  <- sum(df$response == 0L & df$signal_present == 1L)
  FA <- sum(df$response == 1L & df$signal_present == 0L)
  CR <- sum(df$response == 0L & df$signal_present == 0L)
  
  tab <- matrix(c(H, M, FA, CR), ncol = 2,
                dimnames = list(c("Yes", "No"), c("S", "NS")))
  tab <- as.table(tab)
  
  switch(output,
         data  = df,
         table = tab,
         both  = list(data = df, table = tab,
                      params = list(dp = dp, crit = crit,
                                    pH = pH, pFA = pFA)))
}


simulate_sdt(2,1, 300, output="table")

# -----------------------------------------------------------------------------
## Multilevel example!

d <- read.csv("./data/Lisi2023_covidmisinfo.csv")
str(d)

# estimate group level d' using multilevel modelling
library(lme4)

sdt_lisi23 <- glmer(respT ~ isTrue + (isTrue|id) + (1|question), 
                    family=binomial(link="probit"),
                    data = d[d$type=="science",])

summary(sdt_lisi23)
  
  
# alternative example using Andriana's data

d <- read.csv("./data/Theodoropoulou_RR.csv")
str(d)

sdt_lisi23 <- glmer(anty ~ anti_needed + (anti_needed|ID), # + (1|scenario), 
                    family=binomial(link="probit"),
                    data = d[d$condition==1,])

summary(sdt_lisi23)



#############################################################################
# ROC curves
library(lattice)
x <- seq(0, 1, len = 500)
xy <- expand.grid(pFA = x, pH = x)
xy$dp <- with(xy, qnorm(pH) - qnorm(pFA))
contourplot(dp ~ pFA + pH, data = xy,
            subscripts = TRUE,
            panel = function(x, y, z, subscripts, ...){
              panel.levelplot(x, y, z, subscripts, at = c(0,1,2.682),
                              labels = paste("d' =", c(0,1,2.68)),
                              label.style = "mixed",
                              contour = TRUE, region = FALSE,
                              lwd=3, col="grey60")
              panel.points(0.0960, 0.9106,
                           pch = 16, col = "black", cex = 1.25)
              },
            xlab="p(FA)",
            ylab="p(H)" )


library(lattice)
# Create sequence for probabilities (excluding 0 and 1 to avoid infinite z-scores)
x <- seq(0.001, 0.999, len = 500)
xy <- expand.grid(pFA = x, pH = x)

# Calculate z-scores (inverse normal transformation)
xy$zFA <- qnorm(xy$pFA)
xy$zH <- qnorm(xy$pH)

# Calculate d' using z-scores
xy$dp <- with(xy, zH - zFA)

# Create zROC plot
contourplot(dp ~ zFA + zH, data = xy,
            subscripts = TRUE,
            panel = function(x, y, z, subscripts, ...){
              panel.levelplot(x, y, z, subscripts, at = c(0, 1, 2.682),
                              labels = paste("d' =", c(0, 1, 2.68)),
                              label.style = "mixed",
                              contour = TRUE, region = FALSE,
                              lwd = 3, col = "grey60")
              
              # Convert the original point to z-scores
              zFA_point <- qnorm(0.0960)
              zH_point <- qnorm(0.9106)
              panel.points(zFA_point, zH_point,
                           pch = 16, col = "black", cex = 1.25)
            },
            xlab = "z(FA) - z-score of False Alarm Rate",
            ylab = "z(H) - z-score of Hit Rate",
            main = "zROC Curves with d' Contours")


#############################################################################

dprime <- 2
sigma_EV <- 1
sigma_UV <- 1.6   # signal noisier -> zROC slope < 1

c_grid <- seq(-3, 6, length.out = 400)

ROC <- function(dprime, sigma) {
  FA <- 1 - pnorm(c_grid)
  H  <- 1 - pnorm((c_grid - dprime)/sigma)
  list(FA=FA, H=H,
       zFA=qnorm(pmin(pmax(FA,1e-6),1-1e-6)),
       zH =qnorm(pmin(pmax(H, 1e-6),1-1e-6)))
}

ev <- ROC(dprime, sigma_EV)
uv <- ROC(dprime, sigma_UV)

op <- par(mfrow=c(1,2), mar=c(4,4,2,1))

# ROC
plot(ev$FA, ev$H, type="l", lwd=2, col="black",
     xlab="FA", ylab="H"); lines(uv$FA, uv$H, lwd=2, col="grey40")
abline(0,1,lty=2,col="grey60")
legend("bottomright", c("EV (σ=1)","UV (σ=1.6)"), col=c("black","grey40"), lwd=2, bty="n")

# zROC
plot(ev$zFA, ev$zH, type="l", lwd=2, col="black",
     xlab="z(FA)", ylab="z(H)"); lines(uv$zFA, uv$zH, lwd=2, col="grey40")
# Fit slopes to show s
sev <- coef(lm(ev$zH ~ ev$zFA))[2]
suv <- coef(lm(uv$zH ~ uv$zFA))[2]
legend("topleft",
       legend = c(paste0("EV slope ≈ ", round(sev,3)),
                  paste0("UV slope ≈ ", round(suv,3), " (< 1)")),
       bty="n")
par(op)