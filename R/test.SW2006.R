Scafetta2006 <- function (GISS.temp = TRUE, do.MonteCarlo = TRUE, test.bp = TRUE, 
    lag = 0, stepwise = TRUE, interval = 1958:2000, same.interval = TRUE, 
    all.data = FALSE, SW06.coefs.only = FALSE, wavelets.only = FALSE, 
    bivariate = TRUE, figures = TRUE, tables = TRUE, wavelet = TRUE, 
    boundary = "reflection") 
{

  # 'lag' for lagged response in regression analysis
# same.interval <- TRUE:  1900-2000 else 1880-2000
# all.data <- FALSE: overrides same.interval

# R.E. Benestad
# R-script: to execute in R type 'source("test.SW2006.R")'
# This script will only work with an Internet connection.
# R is freely available from CRAN: http//cran.r-project.org
# (Linux, Windows, Mac)
# The R-packages (libraries) are availablefrom CRAN (contributed packages)
# 
# The script should be run with an Internet connection, and retrieves the
# necessary data over URLs. The script saves the data locally, so that it
# only needs to read from URL if not found locally.
#
# Note:
# MODWT MRA: URL http://rss.acs.unt.edu/Rdoc/library/wavelets/html/mra.html
# mra(X, filter="la8", n.levels, boundary="periodic", fast=TRUE, method="dwt")
# Newer version:
# mra(x, wf = "la8", J = 4, method = "modwt", boundary = "periodic")
#
# Scahetta & West (2005):
# D7: between 2^7 - 2^8 = 128-256 month-window -> 88-176 month-band 
#                                                (7.3-14.7 yrs).
# D8: 176-352 month frequency -> 256-512 month-window.

# Preamble -------------------------------------------------------------
    require(waveslim)
    require(lmtest)
    cmon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
        "Aug", "Sep", "Oct", "Nov", "Dec")

    # Functions used by the script: -----------------------------------------
    
    Forcings <- function() {
        #load("Debunking/data/forcings.rda")
        #data("forcings",envir=environment())
        forcing.names <- names(forcings)
        x11()
        years <- as.numeric(row.names(table(trunc(forcings$Year))))
        ny <- length(years)
        nf <- length(forcing.names)
        data <- rep(NA, ny * nf)
        dim(data) <- c(ny, nf)
        data[, 1] <- years
        forcing <- as.matrix(forcings)
        plot(range(forcing[, 1]), range(forcing[, 2:nf]), type = "n", 
            main = "GISS forcing estimates", xlab = "Time", ylab = "Forcing")
        grid()
        lines(forcings[, 1], forcing[, 2])
        lines(forcings[, 1], forcing[, 3], col = "blue")
        lines(forcings[, 1], forcing[, 4], col = "steelblue")
        lines(forcings[, 1], forcing[, 5], col = "darkblue")
        lines(forcings[, 1], forcing[, 6], col = "red")
        lines(forcings[, 1], forcing[, 7], col = "darkred")
        lines(forcings[, 1], forcing[, 8], col = "grey")
        lines(forcings[, 1], forcing[, 9], col = "pink")
        lines(forcings[, 1], forcing[, 10], col = "green")
        lines(forcings[, 1], forcing[, 11], col = "darkgreen")
        for (i in 1:ny) {
            for (ii in 2:nf) {
                iii <- is.element(trunc(forcings[, 1]), data[i, 
                  1])
                data[i, ii] <- mean(forcings[iii, ii], na.rm = TRUE)
            }
        }
        lines(data[, 1], data[, 2], lwd = 2, lty = 2)
        lines(data[, 1], data[, 3], col = "blue", lwd = 2, lty = 2)
        lines(data[, 1], data[, 4], col = "steelblue", lwd = 2, 
            lty = 2)
        lines(data[, 1], data[, 5], col = "darkblue", lwd = 2, 
            lty = 2)
        lines(data[, 1], data[, 6], col = "red", lwd = 2, lty = 2)
        lines(data[, 1], data[, 7], col = "darkred", lwd = 2, 
            lty = 2)
        lines(data[, 1], data[, 8], col = "grey", lwd = 2, lty = 2)
        lines(data[, 1], data[, 9], col = "pink", lwd = 2, lty = 2)
        lines(data[, 1], data[, 10], col = "green", lwd = 2, 
            lty = 2)
        lines(data[, 1], data[, 11], col = "darkgreen", lwd = 2, 
            lty = 2)
        colnames(data) <- forcing.names
        invisible(forcings)
    }
    test.ccf <- function() {
        t <- rnorm(100)
        S <- t
        i1 <- is.element(1:100, 1:100 + 1)
        i2 <- is.element(1:100 + 1, 1:100)
        t1 <- t[i1]
        S1 <- S[i2]
        i1 <- is.element(1:100, 1:100 - 2)
        i2 <- is.element(1:100 - 2, 1:100)
        t2 <- t[i1]
        S2 <- S[i2]
        i1 <- is.element(1:100, 1:100 - 5)
        i2 <- is.element(1:100 - 5, 1:100)
        t3 <- t[i1]
        S3 <- S[i2]
        par(mfcol = c(3, 2))
        plot(t1, type = "l")
        lines(S1, col = "red")
        plot(t2, type = "l")
        lines(S2, col = "red")
        plot(t3, type = "l")
        lines(S3, col = "red")
        ccf(t1, S1)
        grid()
        ccf(t2, S2)
        grid()
        ccf(t3, S3)
        grid()
    }
    adjust <- function(y, y.ref) {
        y <- 0.5 * sd(y.ref, na.rm = TRUE) * (y - mean(y, na.rm = TRUE))/sd(y, 
            na.rm = TRUE) + mean(y.ref, na.rm = TRUE)
        invisible(y)
    }
    lagged <- function(y, t) {
        x <- 1:length(y)
        lagged <- approx(x = x, y = y, xout = x - (t * 12))$y
        lagged
    }
    test.lagged <- function(lag = 0.5) {
        y <- cos(2 * pi * seq(120)/12)
        y.l <- lagged(y, lag)
        x11()
        plot(y, type = "l", lwd = 3, col = "grey", main = "test.lagged")
        grid()
        lines(y.l, lwd = 2, lty = 2)
    }
    T.sun <- function(S.4, D.4, D.3, Z.eq = 0.21, Z.S4 = 0.17, 
        Z.22y = 0.17, Z.11y = 0.11, t.S4 = 4.3, t.4 = 2.5, t.3 = 1.3) {
        print(paste("Z.eq=", Z.eq, "Z.S4=", Z.S4, "Z.22y=", Z.22y, 
            "Z.11y=", Z.11y, "t.S4=", t.S4, "t.4=", t.4, "t.3=", 
            t.3))
        T.sun <- Z.eq * mean(S.4, na.rm = TRUE) + Z.S4 * (lagged(S.4, 
            t.S4) - mean(S.4, na.rm = TRUE)) + Z.22y * lagged(D.4, 
            t.4) + Z.11y * lagged(D.3, t.3)
        T.sun
    }
    ma.filt <- function(x, n) {
        y <- filter(x, rep(1, n)/n)
        y
    }

# This followiung function computes the month, day, and year, given a Julian day.
# The algorithm is taken from Press et al. (1989), "Numerical Recipes 
# in Pascal", Cambridge, p. 13.
#
# This function removes the dependency to outdated packages 'chron' and
# 'date'.
#
# R.E. Benestad, met.no, Oslo, Norway 04.09.2003
# rasmus.benestad@met.no
#------------------------------------------------------------------------
    caldat <- function(julian) {
        igreg = 2299161
        julian <- trunc(julian)
        jalpha <- julian * 0
        ja <- julian * 0
        im <- (julian >= igreg)
        if (sum(im) > 0) {
            jalpha[im] <- trunc(((julian - 1867216) - 0.25)/36524.25)
            ja[im] <- julian + 1 + jalpha - trunc(0.25 * jalpha)
        }
        im <- (julian < igreg)
        if (sum(im) > 0) 
            ja[im] <- julian[im]
        jb <- ja + 1524
        jc <- trunc(6680 + ((jb - 2439870) - 122.1)/365.25)
        jd <- 365 * jc + trunc(0.25 * jc)
        je <- trunc((jb - jd)/30.6001)
        id <- jb - jd - trunc(30.6001 * je)
        mm <- je - 1
        im <- (mm > 12)
        if (sum(im) > 0) 
            mm[im] <- mm[im] - 12
        iyyy <- jc - 4715
        im <- (mm > 2)
        if (sum(im) > 0) 
            iyyy[im] <- iyyy[im] - 1
        im <- (iyyy <= 0)
        if (sum(im) > 0) 
            iyyy <- iyyy - 1
        caldat <- list(month = mm, day = id, year = iyyy)
        invisible(caldat)
    }

# This routine computes the Julian day given a month, day, and year.
# The algorithm is taken from Press et al. (1989), "Numerical Recipes 
# in Pascal", Cambridge, p. 10.
#
# This function removes the dependency to outdated packages 'chron' and
# 'date'.
#
# R.E. Benestad, met.no, Oslo, Norway 04.09.2003
# rasmus.benestad@met.no
# Bug correction. 04.02.2005: 'jy[im] <- iyyy' -> 'jy[im] <- iyyy[im]'
# 'jm[im]' <- 'mm+1' -> 'jm[im] <- mm[im]+1', 
# 'jy[im] <- iyyy-1' -> 'jy[im] <- iyyy[im]-1'
# 'jm[im] <- mm+13' -> 'jm[im] <- mm[im]+13'
# Previous warnings: 'number of items to replace is not a multiple of replacement length'
#------------------------------------------------------------------------
    julday <- function(mm, id, iyyy) {
        igreg <- 588829
        mm <- trunc(mm)
        id <- trunc(id)
        iyyy <- trunc(iyyy)
        im <- (iyyy == 0)
        if (sum(im) > 0) 
            return("There is no year zero!")
        if ((length(mm) != length(id)) | (length(mm) != length(iyyy)) | 
            (length(iyyy) != length(id))) 
            return("The vectors must have same length!")
        im <- (iyyy < 0)
        if (sum(im) > 0) 
            iyyy[im] <- iyyy[im] + 1
        jy <- mm * 0
        jm <- mm * 0
        ja <- mm * 0
        im <- (mm > 2)
        if (sum(im) > 0) {
            jy[im] <- iyyy[im]
            jm[im] <- mm[im] + 1
        }
        im <- (mm <= 2)
        if (sum(im) > 0) {
            jy[im] <- iyyy[im] - 1
            jm[im] <- mm[im] + 13
        }
        jul <- trunc(365.25 * jy) + trunc(30.6001 * jm) + id + 
            1720995
        im <- (id + 31 * (mm + 12 * iyyy) >= igreg)
        if (sum(im) > 0) {
            ja[im] <- trunc(0.01 * jy)
            jul[im] <- jul + 2 - ja[im] + trunc(0.25 * ja[im])
        }
        julday <- jul
        invisible(julday)
    }
    random.walk <- function(x) {
        s <- sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        x <- cumsum(x)
        x <- (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE) * 
            s + m
        x
    }
    monte.carlo <- function(nt = NULL, N = 300, MA.filt = FALSE, 
        tau = 22 * 12, sd.sun = 1, sd.temp = 1, m.sun = 0, m.temp = 0, 
        randomwalk = TRUE, regression = FALSE) {
        print("MONTE CARLO simulations! - please be patient...")
        ratio <- rep(NA, N)
        if (is.null(Y)) 
            nt <- 3000
        if (MA.filt) {
            for (i in 1:N) {
                sun <- rnorm(nt, sd = sd.sun, mean = m.sun)
                temp <- rnorm(nt, sd = sd.temp, mean = m.temp)
                if (randomwalk) {
                  sun <- random.walk(sun)
                  temp <- random.walk(temp)
                }
                sun.low <- ma.filt(sun, round(29.3 * 12))
                sun.high <- sun - sun.low
                sun.bp <- ma.filt(sun.high, round(14.7 * 12))
                temp.low <- ma.filt(temp, round(29.3 * 12))
                temp.high <- temp - temp.low
                temp.bp <- ma.filt(temp.high, round(14.7 * 12))
                ratio[i] <- sd(temp.bp, na.rm = TRUE)/sd(sun.bp, 
                  na.rm = TRUE)
                if (i == 1) {
                  x11()
                  plot(temp.bp, type = "l", lwd = 2, main = "Example: two stochastic series band-pass filtered")
                  lines(sun.bp, col = "grey", lwd = 2)
                }
            }
        }
        else {
            for (i in 1:N) {
                sun <- rnorm(nt, sd = sd.sun, mean = m.sun)
                temp <- rnorm(nt, sd = sd.temp, mean = m.temp)
                if (randomwalk) {
                  sun <- random.walk(sun)
                  temp <- random.walk(temp)
                }
                D.sun <- mra(sun, wf = "la8", J = floor(log(length(sun), 
                  2)), method = "modwt", boundary = boundary)
                D.temp <- mra(temp, wf = "la8", J = floor(log(length(temp), 
                  2)), method = "modwt", boundary = boundary)
                sun.bp <- D.sun$D8
                temp.bp <- D.temp$D8
                A.sun <- lm(sun.bp ~ sin(2 * pi * (1:nt)/tau) + 
                  cos(2 * pi * (1:nt)/tau))
                A.temp <- lm(temp.bp ~ sin(2 * pi * (1:nt)/tau) + 
                  cos(2 * pi * (1:nt)/tau))
                if (regression) 
                  ratio[i] <- sqrt(A.temp$coefficients[2]^2 + 
                    A.temp$coefficients[3]^2)/sqrt(A.sun$coefficients[2]^2 + 
                    A.sun$coefficients[3]^2)
                else ratio[i] <- sd(temp.bp, na.rm = TRUE)/sd(sun.bp, 
                  na.rm = TRUE)
                if (i == 1) {
                  x11()
                  plot(temp.bp, type = "n", col = "darkred", 
                    main = "Example: two stochastic series band-pass filtered")
                  points(temp, col = "red", pch = 19, cex = 0.5)
                  points(sun, col = "blue", pch = 19, cex = 0.5)
                  lines(temp.bp, lwd = 2, col = "darkred")
                  lines(sun.bp, col = "darkblue", lwd = 2)
                  lines(predict(A.sun), col = "blue", lwd = 2, 
                    lty = 2)
                  lines(predict(A.temp), col = "red", lwd = 2, 
                    lty = 2)
                }
            }
        }
        x11()
        hist(ratio, lwd = 2, breaks = seq(floor(min(ratio)), 
            ceiling(max(ratio)), length = 50), main = "A-ratios for 2 band-pass filtered stochastic series")
    }
    year2monthly <- function(year, TSI) {
        monthly <- spline(year, TSI, n = 12 * length(year))
        year2month <- list(y = monthly$y, x = monthly$x)
        invisible(year2month)
    }
    # Last N elements of the vector:
    l.N <- function(x, N) {
        n <- length(x)
        if (N >= n) {
            X <- x
        }
        else X <- x[(n - N + 1):n]
        X
    }
# Estimate the coefficients according to Scafetta & West (2005)
# assumes annual.
    est.coef <- function(TSI, TSI.year, T2m, T2m.year, tau = c(14.7, 
        29.3) * 12, plot.it = FALSE, tag1 = "", tag2 = "") {
        N <- length(T2m)
        y.TSI <- year2monthly(l.N(TSI.year, N), l.N(TSI, N))$y
        x.TSI <- year2monthly(l.N(TSI.year, N), l.N(TSI, N))$x
        D.TSI <- mra(y.TSI, wf = "la8", J = floor(log(length(y.TSI), 
            2)), method = "modwt", boundary = boundary)
        y.T2m <- year2monthly(T2m.year, T2m)$y
        x.T2m <- year2monthly(T2m.year, T2m)$x
        D.T2m <- mra(y.T2m, wf = "la8", J = floor(log(length(y.T2m), 
            2)), method = "modwt", boundary = boundary)
        if (tau[2] == 29.3 * 12) {
            TSI.bp <- D.TSI$D8
            nt.TSI <- length(TSI.bp)
            T2m.bp <- D.T2m$D8
            nt.T2m <- length(T2m.bp)
        }
        else {
            TSI.bp <- D.TSI$D7
            nt.TSI <- length(TSI.bp)
            T2m.bp <- D.T2m$D7
            nt.T2m <- length(T2m.bp)
        }
        A.TSI <- lm(TSI.bp ~ sin(2 * pi * (1:nt.TSI)/mean(tau)) + 
            cos(2 * pi * (1:nt.TSI)/mean(tau)))
        A.T2m <- lm(T2m.bp ~ sin(2 * pi * (1:nt.T2m)/mean(tau)) + 
            cos(2 * pi * (1:nt.T2m)/mean(tau)))
        # Use sin(t + P) = sin t cos P + cos t sin P
        ratio.lm <- sqrt(A.T2m$coefficients[2]^2 + A.T2m$coefficients[3]^2)/sqrt(A.TSI$coefficients[2]^2 + 
            A.TSI$coefficients[3]^2)
        ratio.sd <- sd(T2m.bp)/sd(TSI.bp)
        ratio.SW = round(2 * sqrt(2) * sd(T2m.bp), 2)/round(2 * 
            sqrt(2) * sd(TSI.bp), 2)
        if (plot.it) {
            x11()
            par(mfrow = c(3, 1), mar = c(2, 4, 3, 2) + 0.1)
            plot(x.TSI, y.TSI, pch = 19, col = "grey40", main = paste("TSI + band-passed components. Tau=", 
                tau[1]/12, "-", tau[2]/12, "years, A=", round(2 * 
                  sqrt(2) * sd(TSI.bp), 2), tag1), ylab = "TSI", 
                xlab = "Time")
            grid()
            trendfit <- predict(lm(y.TSI ~ I(x.TSI) + I(x.TSI^2) + 
                I(x.TSI^3) + I(x.TSI^4) + I(x.TSI^5)))
            polygon(c(x.TSI, reverse(x.TSI)), c(trendfit + 2 * 
                sqrt(2) * sd(TSI.bp), reverse(trendfit - 2 * 
                sqrt(2) * sd(TSI.bp))), col = "grey80", border = "grey50")
            points(x.TSI, y.TSI, pch = 19, col = "grey40")
            lines(x.TSI, TSI.bp * 1 + mean(TSI, na.rm = TRUE) + 
                0.25, col = "pink", lwd = 3)
            lines(x.TSI, predict(A.TSI) * 1 + mean(TSI, na.rm = TRUE) + 
                0.25, col = "darkred", lty = 2)
            plot(x.T2m, y.T2m, pch = 19, col = "grey40", ylab = "<T>", 
                xlab = "Time", main = paste("T2m + band-passed components. Tau=", 
                  tau[1]/12, "-", tau[2]/12, "years, A=", round(2 * 
                    sqrt(2) * sd(T2m.bp), 2), tag2))
            grid()
            trendfit <- predict(lm(y.T2m ~ I(x.T2m) + I(x.T2m^2) + 
                I(x.T2m^3) + I(x.T2m^4) + I(x.T2m^5)))
            polygon(c(x.T2m, reverse(x.T2m)), c(trendfit + 2 * 
                sqrt(2) * sd(T2m.bp), reverse(trendfit - 2 * 
                sqrt(2) * sd(T2m.bp))), col = "grey90", border = "grey70")
            if (tau[2] == 29.3 * 12) 
                polygon(c(x.T2m, reverse(x.T2m)), c(trendfit + 
                  0.06, reverse(trendfit - 0.06)), col = "grey70", 
                  border = "grey70", lwd = 3, density = 15, angle = -45)
            points(x.T2m, y.T2m, pch = 19, col = "grey40")
            lines(x.T2m, T2m.bp * 1 + mean(T2m, na.rm = TRUE) + 
                0.2, col = "steelblue", lwd = 3)
            lines(x.T2m, predict(A.T2m) * 1 + mean(T2m, na.rm = TRUE) + 
                0.2, col = "darkblue", lty = 2)
            plot(x.TSI, TSI.bp/sd(TSI.bp), type = "l", col = "pink", 
                lwd = 3, main = paste("Standardised wavelet components (", 
                  "std-ratio=", round(ratio.sd, 2), "A'(T2m)/A'(TSI)-ratio=", 
                  round(ratio.lm, 2), ")"), ylab = "dimensionless", 
                xlab = "Time")
            grid()
            text(1960, -2, paste("A'(sun)=", round(sqrt(A.TSI$coefficients[2]^2 + 
                A.TSI$coefficients[3]^2), 4), "A'(T2m)=", round(sqrt(A.T2m$coefficients[2]^2 + 
                A.T2m$coefficients[3]^2), 4), " [2*sqrt(2)*sigma=", 
                round(2 * sqrt(2) * sd(TSI.bp), 2), " (sun) &", 
                round(2 * sqrt(2) * sd(T2m.bp), 2), " (T2m)]"))
            lines(x.T2m, T2m.bp/sd(T2m.bp), col = "steelblue", 
                lwd = 3, lty = 1)
            lines(x.T2m, D.T2m$D6/sd(D.T2m$D6), lty = 3, col = "grey70")
            lines(x.T2m, D.T2m$D5/sd(D.T2m$D5), lty = 3, col = "grey70")
            lines(x.T2m, D.T2m$D7/sd(D.T2m$D7), lty = 2, col = "grey70")
        }
        ratios <- list(ratio.SW = ratio.SW, ratio.sd = ratio.sd, 
            ratio.lm = ratio.lm, D.TSI = D.TSI, D.T2m = D.T2m, 
            TSI = TSI, TSI.year = TSI.year, x.TSI = x.TSI)
        ratios
    }
    est.coefSW2006b <- function(TSI, TSI.year, T2m, T2m.year, 
        plot.it = TRUE) {
        centuries.t2m <- as.numeric(rownames(table(100 * trunc(T2m.year/100))))
        centuries.tsi <- as.numeric(rownames(table(100 * trunc(TSI.year/100))))
        NC.t2m <- length(centuries.t2m)
        NC.tsi <- length(centuries.tsi)
        coef <- matrix(rep(NA, (NC.tsi - 1) * (NC.t2m - 1)), 
            NC.tsi - 1, NC.t2m - 1)
        for (i in 1:(NC.tsi - 1)) {
            for (j in 1:(NC.t2m - 1)) {
                a <- (mean(T2m[is.element(100 * trunc(T2m.year/100), 
                  centuries.t2m[j + 1])]) - mean(T2m[is.element(100 * 
                  trunc(T2m.year/100), centuries.t2m[j])]))/(mean(TSI[is.element(100 * 
                  trunc(TSI.year/100), centuries.tsi[i + 1])]) - 
                  mean(TSI[is.element(100 * trunc(TSI.year/100), 
                    centuries.tsi[i])]))
                coef[i, j] <- round(a, 3)
            }
        }
        rownames(coef) <- paste(centuries.tsi[2:NC.tsi], "-", 
            centuries.tsi[1:(NC.tsi - 1)], sep = "")
        colnames(coef) <- paste(centuries.t2m[2:NC.t2m], "-", 
            centuries.t2m[1:(NC.t2m - 1)], sep = "")
        print(coef)
        if (plot.it) {
            x11()
            par(mfrow = c(2, 1))
            plot(TSI.year, TSI, type = "l", col = "grey70", main = "Lean (2004) TSI", 
                ylab = "W/m2", xlab = "years")
            grid()
            for (ii in 1:length(centuries.tsi)) {
                icent <- is.element(100 * trunc(TSI.year/100), 
                  centuries.tsi[ii])
                lines(c(0, 99) + centuries.tsi[ii], rep(mean(TSI[icent], 
                  na.rm = TRUE), 2), col = "pink", lwd = 3)
                text(50 + centuries.tsi[ii], mean(TSI[icent], 
                  na.rm = TRUE), round(mean(TSI[icent], na.rm = TRUE), 
                  2))
            }
            plot(T2m.year, T2m, type = "l", col = "grey70", main = "Global <T2m> GISS CTL", 
                ylab = "deg C", xlab = "model years")
            grid()
            for (ii in 1:length(centuries.t2m)) {
                icent <- is.element(100 * trunc(T2m.year/100), 
                  centuries.t2m[ii])
                lines(c(0, 99) + centuries.t2m[ii], rep(mean(T2m[icent], 
                  na.rm = TRUE), 2), col = "pink", lwd = 3)
                text(50 + centuries.t2m[ii], mean(T2m[icent], 
                  na.rm = TRUE), round(mean(T2m[icent], na.rm = TRUE), 
                  3), font = 2, cex = 0.7, srt = 60)
            }
        }
        invisible(coef)
    }
    est.coef.CTL <- function(y1, x1, y2, x2, plot.it = TRUE, 
        step = 10) {
        ny1 <- length(y1)
        ny2 <- length(y2)
        nc <- trunc(abs(ny1 - ny2)/step) - 1
        print(paste("est.coef.CTL: ny1=", ny1, " ny2=", ny2, 
            " nc=", nc))
        coef.ctl.D4 <- rep(NA, nc)
        coef.ctl.D3 <- coef.ctl.D4
        for (i in 1:nc) {
            ii <- seq(i * step, (i * step) + ny1 - 1, by = 1)
            coef.ctl.D4[i] <- est.coef(y1, x1, y2[ii], x2[ii])$ratio.SW
            coef.ctl.D3[i] <- est.coef(y1, x1, y2[ii], x2[ii], 
                tau = c(7.3, 14.7) * 12)$ratio.SW
            print(c(i, range(ii), coef.ctl.D4[i], coef.ctl.D3[i]))
        }
        coef.ctl <- rbind(coef.ctl.D4, coef.ctl.D3)
        rownames(coef.ctl) <- c("~22-year", "~11-year")
    }
    T.eq <- function(S, A = 0.3, rho = 5.67e-08) {
        T.eq <- (S * (1 - A)/(4 * rho))^0.25
        T.eq
    }

    make.table1 <- function(model, Model, text) {
        tab.results <- paste(text, round(summary(model)$coefficients[2], 
            2), " pm ", round(summary(model)$coefficients[4], 
            2), " (", round(summary(model)$coefficients[8], 2), 
            ")  & ", round(summary(Model)$coefficients[2], 2), 
            " pm ", round(summary(Model)$coefficients[5], 2), 
            " (", round(summary(Model)$coefficients[11], 2), 
            ")  & ", round(summary(Model)$coefficients[3], 2), 
            " pm ", round(summary(Model)$coefficients[6], 2), 
            " (", round(summary(Model)$coefficients[12], 2), 
            ")", sep = "")
        tab.results
    }
    make.table2 <- function(model, text) {
        TAB.results <- paste(text, round(summary(model)$coefficients[2], 
            2), " pm ", round(summary(model)$coefficients[13], 
            2), " (", round(summary(model)$coefficients[35], 
            2), ")  & ", round(summary(model)$coefficients[3], 
            2), " pm ", round(summary(model)$coefficients[14], 
            2), " (", round(summary(model)$coefficients[36], 
            2), ")  & ", round(summary(model)$coefficients[4], 
            2), " pm ", round(summary(model)$coefficients[15], 
            2), " (", round(summary(model)$coefficients[37], 
            2), ")  & ", round(summary(model)$coefficients[5], 
            2), " pm ", round(summary(model)$coefficients[16], 
            2), " (", round(summary(model)$coefficients[38], 
            2), ")  & ", round(summary(model)$coefficients[6], 
            2), " pm ", round(summary(model)$coefficients[17], 
            2), " (", round(summary(model)$coefficients[39], 
            2), ") & ", round(summary(model)$coefficients[7], 
            2), " pm ", round(summary(model)$coefficients[18], 
            2), " (", round(summary(model)$coefficients[40], 
            2), ") & ", round(summary(model)$coefficients[8], 
            2), " pm ", round(summary(model)$coefficients[19], 
            2), " (", round(summary(model)$coefficients[41], 
            2), ")", sep = "")
        TAB.results
    }

# Obtain the data ------------------------------------------------------
# Obtain the Lean TSI data... Read over the Internet

    print("START")
    print(paste("Lag=", lag))
    cmon <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", 
        "Aug", "Sep", "Oct", "Nov", "Dec")
    if (GISS.temp) {
        #data("gistemp",envir=environment())
        #load("Debunking/data/gistemp.rda")
        #attach(gistemp)
        T2m <- gistemp$T2m
        year <- gistemp$year
        level.1900 <- gistemp$level.1900
        X5.year.mean <- gistemp$X5.year.mean
        T.descr <- gistemp$T.descr
    }
    else {
        print("Get HadCRUT from URL")
        #data("crutemp",envir=environment())
        #load("Debunking/data/crutemp.rda")
        #attach(crutemp)
        gistemp <- crutemp
        T2m <- gistemp$T2m
        year <- gistemp$year
        T.descr <- gistemp$T.descr
        gistemp$level.1900 <- 0
    }
    T2m.mon <- year2monthly(year, T2m)$y
    yr.mon <- year2monthly(year, T2m)$x
    #data("AKRIM",envir=environment())
    #load("Debunking/data/AKRIM.rda")
    #data("Lean1995",envir=environment())
    #load("Debunking/data/Lean1995.rda")
    #data("Lean2004",envir=environment())
    #load("Debunking/data/Lean2004.rda")
    Lean2004.TSI <- Lean2004$X11yrCYCLE.BKGRND
    forcings <- Forcings()
    print("Combine Lean95 & AKRIM: match just one year")
    i1.match <- is.element(trunc(Lean1995$year), 1980)
    i2.match <- is.element(AKRIM$year, 1980)
    adj <- AKRIM$TSI[i2.match] - Lean1995$TSI[i1.match]
    print(paste("adj=", adj, "sum(i1.match)=", sum(i1.match), 
        " sum(i2.match)=", sum(i2.match)))
    i.less <- Lean1995$year < min(AKRIM$year)
    TSI.comb <- c(Lean1995$TSI[i.less] + adj, AKRIM$TSI)
    year.comb <- trunc(c(Lean1995$year[i.less], AKRIM$year))
    i1.match.1 <- is.element(trunc(Lean1995$year), min(AKRIM$year))
    i2.match.1 <- is.element(AKRIM$year, min(AKRIM$year))
    adj.1 <- AKRIM$TSI[i2.match.1] - Lean1995$TSI[i1.match.1]
    print(paste("adj.1=", adj.1, "sum(i1.match.1)=", sum(i1.match.1), 
        " sum(i2.match.1)=", sum(i2.match.1)))
    print(paste("year 1=", min(AKRIM$year), " No obs year 1=", 
        sum(is.element(AKRIM$akrim.dates$year, min(AKRIM$year)))))
    TSI.comb.1 <- c(Lean1995$TSI[i.less] + adj.1, AKRIM$TSI)

    # Preprocess ---------------------------------------------------------
    print("Combine Lean95 & AKRIM: match just all years")
    ii1.match <- is.element(trunc(Lean1995$year), AKRIM$year)
    ii2.match <- is.element(AKRIM$year, trunc(Lean1995$year))
    adj.all <- mean(AKRIM$TSI[ii2.match], na.rm = TRUE) - mean(Lean1995$TSI[ii1.match], 
        na.rm = TRUE)
    print(paste("adj.all=", adj.all, "sum(ii1.match)=", sum(ii1.match), 
        " sum(ii2.match)=", sum(ii2.match)))
    TSI.comb.all <- c(Lean1995$TSI[i.less] + adj.all, AKRIM$TSI)

# Decompose the solar signal using maximum overlap discrete wavelet
# transform multiresolution analysis (MODWT & MRA):

    y <- year2monthly(Lean2004$YEAR, Lean2004$X11yrCYCLE.BKGRND)$y
    x <- year2monthly(Lean2004$YEAR, Lean2004$X11yrCYCLE.BKGRND)$x
    D <- mra(y, wf = "la8", J = floor(log(length(y), 2)), method = "modwt", 
        boundary = boundary)
    y.SW <- year2monthly(year.comb, TSI.comb)$y
    x.SW <- year2monthly(year.comb, TSI.comb)$x
    D.SW <- mra(y.SW, wf = "la8", J = floor(log(length(y), 2)), 
        method = "modwt", boundary = boundary)
    y.all <- year2monthly(year.comb, TSI.comb.all)$y
    x.all <- year2monthly(year.comb, TSI.comb.all)$x
    D.all <- mra(y.all, wf = "la8", J = floor(log(length(y), 
        2)), method = "modwt", boundary = boundary)
    S.4 <- D$D9 + D$D10 + D$D11 + D$D12
    D.4 <- D$D8
    D.3 <- D$D7
    S.4.SW <- D.SW$D9 + D.SW$D10 + D.SW$D11 + D.SW$D12
    D.4.SW <- D.SW$D8
    D.3.SW <- D.SW$D7
    S.4.all <- D.all$D9 + D.all$D10 + D.all$D11 + D.all$D12
    D.4.all <- D.all$D8
    D.3.all <- D.all$D7
    #data("Mauna.Loa",envir=environment())
    #load("Debunking/data/Mauna.Loa.rda")
    SW2006 <- list(D, x, y, T.sun, Lean2004, x.SW, D.SW, y.SW, 
        x.all, D.all, y.all)
    #save(file = "Debunking/data/SW2006.rda", SW2006)

    #---------------------------------------------------------------------
    # TSI estimates:

    Y <- T.sun(S.4, D.4, D.3)
    Y.SW <- T.sun(S.4.SW, D.4.SW, D.3.SW)
    Y.all <- T.sun(S.4.all, D.4.all, D.3.all)
    Y.Zeq010 <- T.sun(S.4, D.4, D.3, Z.eq = 0.1)
    Y.Zeq000 <- T.sun(S.4, D.4, D.3, Z.eq = 0)
    Y.Z4005 <- T.sun(S.4, D.4, D.3, Z.S4 = 0.05)
    Y.Z4030 <- T.sun(S.4, D.4, D.3, Z.S4 = 0.3)
    Y.22y005 <- T.sun(S.4, D.4, D.3, Z.22y = 0.05)
    Y.22y030 <- T.sun(S.4, D.4, D.3, Z.22y = 0.3)
    Y.11y005 <- T.sun(S.4, D.4, D.3, Z.11y = 0.05)
    Y.11y030 <- T.sun(S.4, D.4, D.3, Z.11y = 0.3)
    Y.t305 <- T.sun(S.4, D.4, D.3, t.3 = 0.5)
    Y.t330 <- T.sun(S.4, D.4, D.3, t.3 = 3)
    Y.t4005 <- T.sun(S.4, D.4, D.3, t.4 = 0.5)
    Y.t4100 <- T.sun(S.4, D.4, D.3, t.4 = 10)

    # Trends:
    ii <- is.element(gistemp$year, seq(1980, 2000, by = 1))
    Z0 <- data.frame(y = gistemp$X5.year.mean[ii] - gistemp$level.1900, 
        x = gistemp$year[ii])
    trend0 <- lm(y ~ x, data = Z0)
    ii <- is.element(trunc(x), seq(1980, 2000, by = 1))
    Z.S <- data.frame(y = Y[ii], x = x[ii])
    trend1 <- lm(y ~ x, data = Z.S)
    ii <- is.element(trunc(x.SW), seq(1980, 2000, by = 1))
    Z.SW <- data.frame(y = Y.SW[ii], x = x.SW[ii])
    trend2 <- lm(y ~ x, data = Z.SW)
    Z.all <- data.frame(y = Y.all[ii], x = x.SW[ii])
    trend3 <- lm(y ~ x, data = Z.all)

    # Spectral analysis:

    T.monthly <- spline(gistemp$year, gistemp$T2m - gistemp$level.1900, 
        n = 12 * length(gistemp$year))$y

    x11()
    sp.T <- spectrum(T.monthly)
    sp.S <- spectrum(Y[100:length(Y)])
    sp.SW <- spectrum(Y.SW[100:length(Y.SW)])
    sp.D4 <- spectrum(D.4)
    sp.D3 <- spectrum(D.3)
    sp.S4 <- spectrum(S.4)
    dev.off()
    print("Use the GISS simulations to estimate tranfer fuctions")
    print("SW2005:")
    #data("GISS.GCMs",envir=environment())
    #load("Debunking/data/GISS.GCMs.rda")
    ngcms <- length(GISS.GCMs$experiments)
    coefs.D3 <- matrix(rep(NA, 6 * ngcms), 6, ngcms)
    coefs.D4 <- coefs.D3
    coefs.d4 <- coefs.D3
    coefs.d3 <- coefs.D3
    cols <- c(" ")
    if (figures) {
        x11()
        plot(c(1870, 2000), c(13, 15), type = "n", main = T.descr, 
            xlab = "year", ylab = "<T> (deg C)")
        grid()
        colours <- c("grey", "red", "blue")
        for (iv in 1:ngcms) {
            giss <- GISS.GCMs[[iv + 1]]
            head <- attr(giss, "info")
            print(c(all.data, same.interval))
            print(" ************### Estimate SW06a coefficients ###*************")
            if (!all.data) {
                if (same.interval) {
                  ix <- is.element(trunc(Lean2004$YEAR), 1900:2000)
                  IX <- is.element(giss[, 1], 1900:2000)
                  print("--- SAME INTERVAL: 1900 - 2000 ---")
                }
                else {
                  ix <- is.element(trunc(Lean2004$YEAR), 1880:2000)
                  IX <- is.element(giss[, 1], 1880:2000)
                  print("--- EXTENDED INTERVAL: 1880 - 2000 ---")
                }
            }
            else {
                ix <- is.finite(Lean2004$YEAR)
                IX <- is.finite(giss[, 1])
                print("--- ALL YEARS ---")
                print(range(Lean2004$YEAR))
                print(range(giss[, 1]))
            }
            coefs.D4[1, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 2], giss[IX, 1])$ratio.SW
            coefs.D4[2, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 3], giss[IX, 1])$ratio.SW
            coefs.D4[3, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 4], giss[IX, 1])$ratio.SW
            coefs.D4[4, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 5], giss[IX, 1])$ratio.SW
            coefs.D4[5, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 6], giss[IX, 1])$ratio.SW
            coefs.D4[6, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                rowMeans(giss[IX, 2:6]), giss[IX, 1])$ratio.SW
            coefs.D3[1, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 2], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.SW
            coefs.D3[2, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 3], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.SW
            coefs.D3[3, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 4], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.SW
            coefs.D3[4, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 5], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.SW
            coefs.D3[5, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 6], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.SW
            coefs.D3[6, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                rowMeans(giss[IX, 2:6]), giss[IX, 1], tau = c(7.3, 
                  14.7) * 12)$ratio.SW
            coefs.d4[1, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 2], giss[IX, 1])$ratio.lm
            coefs.d4[2, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 3], giss[IX, 1])$ratio.lm
            coefs.d4[3, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 4], giss[IX, 1])$ratio.lm
            coefs.d4[4, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 5], giss[IX, 1])$ratio.lm
            coefs.d4[5, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 6], giss[IX, 1])$ratio.lm
            coefs.d4[6, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                rowMeans(giss[IX, 2:6]), giss[IX, 1])$ratio.lm
            coefs.d3[1, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 2], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.lm
            coefs.d3[2, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 3], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.lm
            coefs.d3[3, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 4], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.lm
            coefs.d3[4, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 5], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.lm
            coefs.d3[5, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                giss[IX, 6], giss[IX, 1], tau = c(7.3, 14.7) * 
                  12)$ratio.lm
            coefs.d3[6, iv] <- est.coef(Lean2004.TSI[ix], Lean2004$YEAR[ix], 
                rowMeans(giss[IX, 2:6]), giss[IX, 1], tau = c(7.3, 
                  14.7) * 12)$ratio.lm
            ensemble.year <- giss[, 1]
            if (lower.case(substr(head[1], 1, 3)) == "all") {
                ensemble.mean.all <- rowMeans(giss[, 2:6])
                all.forcing <- giss[, 2:6]
            }
            if (lower.case(substr(head[1], 1, 3)) == "wel") 
                ensemble.mean.ghg <- rowMeans(giss[, 2:6])
            if (lower.case(substr(head[1], 1, 3)) == "sol") 
                ensemble.mean.tsi <- rowMeans(giss[, 2:6])
            for (ix in 2:6) lines(giss[, 1], giss[, ix], lwd = 1, 
                col = colours[iv])
            lines(giss[, 1], rowMeans(giss[, 2:6]), lwd = 3, 
                col = colours[iv])
        }
        legend(1870, 15, GISS.GCMs$experiments, col = colours, 
            bg = "grey95", lty = 1, lwd = 3)
        print("D3 (Z_11y):")
        rownames(coefs.D3) <- c(as.character(1:5), "mean")
        colnames(coefs.D3) <- cols[-1]
        print(round(coefs.D3, 2))
        print("D4 (Z_22y):")
        rownames(coefs.D4) <- c(as.character(1:5), "mean")
        colnames(coefs.D4) <- cols[-1]
        print(round(coefs.D4, 2))
        print("d3: (lm;  Z_11y)")
        rownames(coefs.d3) <- c(as.character(1:5), "mean")
        colnames(coefs.d3) <- cols[-1]
        print(round(coefs.d3, 2))
        print("d4 (lm; Z_22y):")
        rownames(coefs.d4) <- c(as.character(1:5), "mean")
        colnames(coefs.d4) <- cols[-1]
        print(round(coefs.d4, 2))
        if (SW06.coefs.only) 
            stop()
        print("SW2006b:")
        giss <- GISS.GCMs$SAT.E3oM20
        head <- attr(giss, "info")
        coefSW2006b <- est.coefSW2006b(Lean2004.TSI, Lean2004$YEAR, 
            giss[, 2], giss[, 1], plot.it = TRUE)
        x11()
        sp.ctl <- spectrum(giss[, 2])
        sp.co2 <- spectrum(Mauna.Loa$interpolated)
        dev.off()
        print("Using CTL to estimate band-pass filtered coefficients")
        coefs.CTL <- est.coef.CTL(Lean2004.TSI, Lean2004$YEAR, 
            giss[, 2], giss[, 1], plot.it = TRUE)
        print(coefs.CTL)
    }

    # -------------------------------- Wavelet analyses ----------------------------------

    if (all.data) {
        interval1 <- is.finite(ensemble.year)
        interval2 <- is.finite(forcings$Year - lag)
        interval3 <- is.finite(trunc(Lean2004$YEAR) - lag)
        interval4 <- is.finite(gistemp$year)
    }
    else {
        interval1 <- is.element(ensemble.year, 1900:2000)
        interval2 <- is.element(forcings$Year - lag, 1900:2000)
        interval3 <- is.element(trunc(Lean2004$YEAR) - lag, 1900:2000)
        interval4 <- is.element(gistemp$year, 1900:2000)
    }
    iii1 <- is.element(ensemble.year, forcings$Year - lag) & 
        is.element(ensemble.year, gistemp$year) & interval1
    iii2 <- is.element(forcings$Year - lag, ensemble.year) & 
        is.element(forcings$Year - lag, gistemp$year) & 
        interval2
    iii3 <- is.element(trunc(Lean2004$YEAR) - lag, ensemble.year) & 
        is.element(trunc(Lean2004$YEAR) - lag, gistemp$year) & 
        interval3
    iii4 <- is.element(gistemp$year, forcings$Year - lag) & 
        is.element(gistemp$year, ensemble.year) & interval4
    tag1 <- ""
    tag2 <- ""
    if (!GISS.temp) 
        tag2 <- "HadCRUT3v"
    if (figures) {
        a <- est.coef(Lean2004.TSI[iii3], Lean2004$YEAR[iii3], 
            gistemp$T2m[iii4], gistemp$year[iii4], plot.it = TRUE, 
            tag2 = tag2)
        b <- est.coef(Lean2004.TSI[iii3], Lean2004$YEAR[iii3], 
            gistemp$T2m[iii4], gistemp$year[iii4], plot.it = TRUE, 
            tau = c(7.3, 14.7) * 12, tag2 = tag2)
        jjj3 <- is.element(trunc(Lean1995$year) - lag, 1900:2000)
        jjj4 <- is.element(gistemp$year, 1900:2000)
        c <- est.coef(Lean1995$TSI[jjj3], Lean1995$year[jjj3], 
            gistemp$T2m[jjj4], gistemp$year[jjj4], plot.it = TRUE, 
            tag1 = " Lean(1995)", tag2 = tag2)
        d <- est.coef(Lean1995$TSI[jjj3], Lean1995$year[jjj3], 
            gistemp$T2m[jjj4], gistemp$year[jjj4], plot.it = TRUE, 
            tag1 = " Lean(1995)", tag2 = tag2, tau = c(7.3, 14.7) * 
                12)
        jjj3 <- is.element(trunc(Lean1995$year) - lag, 1980:2002)
        jjj4 <- is.element(gistemp$year, 1980:2002)
        e <- est.coef(AKRIM$TSI, AKRIM$year, gistemp$T2m[jjj4], 
            gistemp$year[jjj4], plot.it = TRUE, tag1 = " Lean(1995)", 
            tag2 = tag2)
        f <- est.coef(AKRIM$TSI, AKRIM$year, gistemp$T2m[jjj4], 
            gistemp$year[jjj4], plot.it = TRUE, tag1 = " Lean(1995)", 
            tag2 = tag2, tau = c(7.3, 14.7) * 12)
        par(mfrow = c(1, 1))
        plot(a$TSI.year, a$TSI, type = "b", pch = 19, main = "Lean (2004) S", 
            ylab = "TSI (W/m^2)", xlab = "year")
        lines(a$x.TSI, a$D.TSI$S10, col = "grey", lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10, col = "steelblue", 
            lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10 + a$D.TSI$D9, 
            col = "blue", lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10 + a$D.TSI$D9 + 
            a$D.TSI$D8, col = "darkblue", lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10 + a$D.TSI$D9 + 
            a$D.TSI$D8 + a$D.TSI$D7, col = "darkred", lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10 + a$D.TSI$D9 + 
            a$D.TSI$D8 + a$D.TSI$D7 + a$D.TSI$D6, col = "red", 
            lwd = 2)
        lines(a$x.TSI, a$D.TSI$S10 + a$D.TSI$D10 + a$D.TSI$D9 + 
            a$D.TSI$D8 + a$D.TSI$D7 + a$D.TSI$D6 + a$D.TSI$D5, 
            col = "darkgreen", lwd = 2)
        legend(1950, 1365.3, c("S10...", "+D10", "+D9  ", "+D8  ", 
            "+D7  ", "+D6  ", "+D5  "), lty = 1, lwd = 2, col = c("grey", 
            "steelblue", "blue", "darkblue", "darkred", "red", 
            "darkgreen"))
    }
    if (wavelet) {
        print("Standard deviations of wavelet components:")
        print(c(1, sd(c$D.TSI$D1), 2 * sqrt(2) * sd(c$D.TSI$D1)))
        wl <- c$D.TSI$D1
        print(c(2, sd(c$D.TSI$D2), 2 * sqrt(2) * sd(c$D.TSI$D2), 
            sd(wl + c$D.TSI$D2), 2 * sqrt(2) * sd(wl + c$D.TSI$D2)))
        wl <- wl + c$D.TSI$D2
        print(c(3, sd(c$D.TSI$D3), 2 * sqrt(2) * sd(c$D.TSI$D3), 
            sd(wl + c$D.TSI$D3), 2 * sqrt(2) * sd(wl + c$D.TSI$D3)))
        wl <- wl + c$D.TSI$D3
        print(c(4, sd(c$D.TSI$D4), 2 * sqrt(2) * sd(c$D.TSI$D4), 
            sd(wl + c$D.TSI$D4), 2 * sqrt(2) * sd(wl + c$D.TSI$D4)))
        wl <- wl + c$D.TSI$D4
        print(c(5, sd(c$D.TSI$D5), 2 * sqrt(2) * sd(c$D.TSI$D5), 
            sd(wl + c$D.TSI$D5), 2 * sqrt(2) * sd(wl + c$D.TSI$D5)))
        wl <- wl + c$D.TSI$D5
        print(c(6, sd(c$D.TSI$D6), 2 * sqrt(2) * sd(c$D.TSI$D6), 
            sd(wl + c$D.TSI$D6), 2 * sqrt(2) * sd(wl + c$D.TSI$D6)))
        wl <- wl + c$D.TSI$D6
        print(c(7, sd(c$D.TSI$D7), 2 * sqrt(2) * sd(c$D.TSI$D7), 
            sd(wl + c$D.TSI$D7), 2 * sqrt(2) * sd(wl + c$D.TSI$D7)))
        wl <- wl + c$D.TSI$D7
        print(c(8, sd(c$D.TSI$D8), 2 * sqrt(2) * sd(c$D.TSI$D8), 
            sd(wl + c$D.TSI$D8), 2 * sqrt(2) * sd(wl + c$D.TSI$D8)))
        wl <- wl + c$D.TSI$D8
        print(c(9, sd(c$D.TSI$D9), 2 * sqrt(2) * sd(c$D.TSI$D9), 
            sd(wl + c$D.TSI$D9), 2 * sqrt(2) * sd(wl + c$D.TSI$D9)))
        wl <- wl + c$D.TSI$D9
        print(c(7, sd(c$D.TSI$D10), 2 * sqrt(2) * sd(c$D.TSI$D10), 
            sd(wl + c$D.TSI$D10), 2 * sqrt(2) * sd(wl + c$D.TSI$D10)))
    }
    if (wavelets.only) 
        stop()

    #Regression analysis:---------------------------------------------------------------

    yr.co2 <- as.numeric(rownames(table(Mauna.Loa$year)))
    n.co2 <- length(yr.co2)
    co2 <- rep(NA,n.co2)
    for (i in 1:n.co2) {
      ii <- is.element(Mauna.Loa$year,yr.co2[i])
      co2[i] <- mean(Mauna.Loa$interpolated[ii])
    }
    
    print("Regression:")
    i1 <- is.element(ensemble.year, trunc(Lean2004$YEAR))
    i2 <- is.element(trunc(Lean2004$YEAR), ensemble.year)
    i5 <- is.element(year.comb, ensemble.year[i1])
    ii1 <- is.element(ensemble.year, trunc(Lean2004$YEAR) - lag) & 
        is.element(ensemble.year, yr.co2 - lag) & is.element(ensemble.year, 
        gistemp$year)
    ii2 <- is.element(trunc(Lean2004$YEAR) - lag, ensemble.year) & 
        is.element(trunc(Lean2004$YEAR) - lag, yr.co2 - 
            lag) & is.element(trunc(Lean2004$YEAR) - lag, gistemp$year)
    ii3 <- is.element(yr.co2 - lag, trunc(Lean2004$YEAR) - 
        lag) & is.element(yr.co2 - lag, ensemble.year) & 
        is.element(yr.co2 - lag, gistemp$year)
    ii4 <- is.element(gistemp$year, trunc(Lean2004$YEAR) - lag) & 
        is.element(gistemp$year, yr.co2 - lag) & is.element(gistemp$year, 
        ensemble.year)
    ii5 <- is.element(year.comb, ensemble.year[ii1])

    # Trend analysis: 1900-2000
    j1 <- is.element(gistemp$year, 1900:2000)
    j2 <- is.element(gistemp$year, 1900:2000)
    lin.trend.fit.0 <- lm(gistemp$T2m[j1] - mean(gistemp$T2m[j1]) ~ 
        gistemp$year[j1])
    if (figures) {
        plot(gistemp$year[j1], gistemp$T2m[j1], ylim = c(-0.5, 
            1), type = "b", pch = 19)
        abline(lin.trend.fit.0, lty = 2)
        grid()
        lin.trend.fit.1 <- lm(ensemble.mean.all[j2] - mean(ensemble.mean.all[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.all - mean(ensemble.mean.all), 
            type = "b", col = "red")
        abline(lin.trend.fit.1, lty = 2, col = "red")
        lin.trend.fit.2 <- lm(ensemble.mean.ghg[j2] - mean(ensemble.mean.ghg[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.ghg - mean(ensemble.mean.ghg), 
            type = "b", col = "darkgreen")
        abline(lin.trend.fit.2, lty = 2, col = "darkgreen")
        lin.trend.fit.3 <- lm(ensemble.mean.tsi[j2] - mean(ensemble.mean.tsi[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.tsi - mean(ensemble.mean.tsi), 
            type = "b", col = "blue")
        abline(lin.trend.fit.3, lty = 2, col = "blue")
        legend(1880, 1, c(round(lin.trend.fit.0$coefficients[2] * 
            100, 2), round(lin.trend.fit.1$coefficients[2] * 
            100, 2), round(lin.trend.fit.2$coefficients[2] * 
            100, 2), round(lin.trend.fit.3$coefficients[2] * 
            100, 2)), col = c("black", "red", "darkgreen", "blue"), 
            lty = 1)
    }
    print(paste("Obs. warming over 100 years:", round(lin.trend.fit.0$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.0)$coefficients[4] * 
        100, 2), "K"))
    print(paste("Warming in 'all' over 100 years:", round(lin.trend.fit.1$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.1)$coefficients[4] * 
        100, 2), "K"))
    print(paste("Warming in 'ghg' over 100 years:", round(lin.trend.fit.2$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.2)$coefficients[4] * 
        100, 2), "K"))
    print(paste("Warming in 'sol' over 100 years ('all'):", round(lin.trend.fit.3$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        100, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.1$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.1)$coefficients[4]/summary(lin.trend.fit.1)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))
    print(paste("Warming in 'sol' over 100 years ('obs'):", round(lin.trend.fit.3$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        100, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.0$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.0)$coefficients[4]/summary(lin.trend.fit.0)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))
    print(paste("Warming in 'sol' over 100 years ('ghg'):", round(lin.trend.fit.3$coefficients[2] * 
        100, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        100, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.2$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.2)$coefficients[4]/summary(lin.trend.fit.2)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))

    # Trend analysis: 1980-2002

    j1 <- is.element(gistemp$year, 1980:2003)
    j2 <- is.element(gistemp$year, 1980:2003)
    lin.trend.fit.0 <- lm(gistemp$T2m[j1] - mean(gistemp$T2m[j1]) ~ 
        gistemp$year[j1])
    if (figures) {
        plot(gistemp$year[j1], gistemp$T2m[j1], ylim = c(-0.5, 
            1), type = "b", pch = 19)
        abline(lin.trend.fit.0, lty = 2)
        grid()
        lin.trend.fit.1 <- lm(ensemble.mean.all[j2] - mean(ensemble.mean.all[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.all - mean(ensemble.mean.all), 
            type = "b", col = "red")
        abline(lin.trend.fit.1, lty = 2, col = "red")
        lin.trend.fit.2 <- lm(ensemble.mean.ghg[j2] - mean(ensemble.mean.ghg[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.ghg - mean(ensemble.mean.ghg), 
            type = "b", col = "darkgreen")
        abline(lin.trend.fit.2, lty = 2, col = "darkgreen")
        lin.trend.fit.3 <- lm(ensemble.mean.tsi[j2] - mean(ensemble.mean.tsi[j2]) ~ 
            ensemble.year[j2])
        lines(ensemble.year, ensemble.mean.tsi - mean(ensemble.mean.tsi), 
            type = "b", col = "blue")
        abline(lin.trend.fit.3, lty = 2, col = "blue")
        legend(1980, 1, c(round(lin.trend.fit.0$coefficients[2] * 
            24, 2), round(lin.trend.fit.1$coefficients[2] * 24, 
            2), round(lin.trend.fit.2$coefficients[2] * 24, 2), 
            round(lin.trend.fit.3$coefficients[2] * 24, 2)), 
            col = c("black", "red", "darkgreen", "blue"), lty = 1)
    }
    print(paste("Obs. warming over 1980-2003:", round(lin.trend.fit.0$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.0)$coefficients[4] * 
        24, 2), "K"))
    print(paste("Warming in 'all' over 1980-2003:", round(lin.trend.fit.1$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.1)$coefficients[4] * 
        24, 2), "K"))
    print(paste("Warming in 'ghg' over 1980-2003:", round(lin.trend.fit.2$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.2)$coefficients[4] * 
        24, 2), "K"))
    print(paste("Warming in 'sol' over 100 years ('all'):", round(lin.trend.fit.3$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        24, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.1$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.1)$coefficients[4]/summary(lin.trend.fit.1)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))
    print(paste("Warming in 'sol' over 100 years ('obs'):", round(lin.trend.fit.3$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        24, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.0$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.0)$coefficients[4]/summary(lin.trend.fit.0)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))
    print(paste("Warming in 'sol' over 100 years ('ghg'):", round(lin.trend.fit.3$coefficients[2] * 
        24, 2), "+-", round(summary(lin.trend.fit.3)$coefficients[4] * 
        24, 2), "K ->", round(100 * lin.trend.fit.3$coefficients[2]/lin.trend.fit.2$coefficients[2], 
        2), "+-", round(100 * sqrt((summary(lin.trend.fit.2)$coefficients[4]/summary(lin.trend.fit.2)$coefficients[2])^2 + 
        (summary(lin.trend.fit.3)$coefficients[4]/summary(lin.trend.fit.3)$coefficients[2])^2), 
        2), "%"))
    print(paste("Estimated Z_22y=", round(a$ratio.sd, 2)))
    iv1 <- is.element(ensemble.year, forcings$Year - lag) & 
        is.element(ensemble.year, gistemp$year) & is.element(ensemble.year, 
        interval)
    iv2 <- is.element(forcings$Year - lag, ensemble.year) & 
        is.element(forcings$Year - lag, gistemp$year) & 
        is.element(forcings$Year - lag, interval)
    iv4 <- is.element(gistemp$year, forcings$Year - lag) & 
        is.element(gistemp$year, ensemble.year) & is.element(gistemp$year, 
        interval)
    print(c(sum(i1), sum(i2), sum(i5)))
    print(c(sum(ii1), sum(ii2), sum(ii3), sum(ii4)))
    print(c(sum(iii1), sum(iii2), sum(iii3), sum(iii4)))
    print(c(sum(iv1), sum(iv2), sum(iv4)))
    gis.res <- ensemble.mean.all[i1] - ensemble.mean.ghg[i1]
    gis.res <- gis.res - mean(gis.res)
    tsi <- Lean2004$X11yrCYCLE.BKGRND[i2]
    calibrate.all <- data.frame(y = ensemble.mean.all[i1] - mean(ensemble.mean.all[i1], 
        na.rm = TRUE), x = tsi)
    calibrate.ghg <- data.frame(y = ensemble.mean.ghg[i1] - mean(ensemble.mean.ghg[i1], 
        na.rm = TRUE), x = tsi)
    calibrate.res <- data.frame(y = gis.res, x = tsi)
    calibrate.tsi <- data.frame(y = ensemble.mean.tsi[i1] - mean(ensemble.mean.tsi[i1], 
        na.rm = TRUE), x = tsi)
    calibrate.SW <- data.frame(y = ensemble.mean.tsi[i1] - mean(ensemble.mean.tsi[i1], 
        na.rm = TRUE), x = TSI.comb[i5])
    Gis.res <- ensemble.mean.all[ii1] - ensemble.mean.ghg[ii1]
    Gis.res <- Gis.res - mean(Gis.res)
    tsi <- Lean2004$X11yrCYCLE.BKGRND[ii2]
    cal.interval <- trunc(Lean2004$YEAR[ii2])
    Calibrate.obs <- data.frame(y = gistemp$T2m[ii4] - mean(gistemp$T2m[ii4], 
        na.rm = TRUE), x1 = 0.7/4 * tsi, x2 = 5.35 * log(co2[ii3]))
    Calibrate.all <- data.frame(y = ensemble.mean.all[ii1] - 
        mean(ensemble.mean.all[ii1], na.rm = TRUE), x1 = 0.7/4 * 
        tsi, x2 = 5.35 * log(co2[ii3]))
    Calibrate.ghg <- data.frame(y = ensemble.mean.ghg[ii1] - 
        mean(ensemble.mean.ghg[ii1], na.rm = TRUE), x1 = 0.7/4 * 
        tsi, x2 = 5.35 * log(co2[ii3]))
    Calibrate.res <- data.frame(y = Gis.res, x1 = 0.7/4 * tsi, 
        x2 = 5.35 * log(co2[ii3]))
    Calibrate.tsi <- data.frame(y = ensemble.mean.tsi[ii1] - 
        mean(ensemble.mean.tsi[ii1], na.rm = TRUE), x1 = 0.7/4 * 
        tsi, x2 = 5.35 * log(co2[ii3]))
    Calibrate.SW <- data.frame(y = gistemp$T2m[ii4] - mean(gistemp$T2m[ii4], 
        na.rm = TRUE), x1 = 0.7/4 * TSI.comb[ii5], x2 = 5.35 * 
        log(co2[ii3]))
    CALIBRATE.obs <- data.frame(y = gistemp$T2m[iii4] - mean(gistemp$T2m[iii4], 
        na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.all <- data.frame(y = ensemble.mean.all[iii1] - 
        mean(ensemble.mean.all[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2], x3 = forcings$O3[iii2], 
        x4 = forcings$StratH2O[iii2], x5 = forcings$LandUse[iii2], 
        x6 = forcings$SnowAlb[iii2], x7 = forcings$StratAer[iii2], 
        x8 = forcings$BC[iii2], x9 = forcings$ReflAer[iii2], 
        x10 = forcings$AIE[iii2])
    CALIBRATE.all1 <- data.frame(y = all.forcing[iii1, 1] - mean(all.forcing[iii1, 
        1], na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.all2 <- data.frame(y = all.forcing[iii1, 2] - mean(all.forcing[iii1, 
        2], na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.all3 <- data.frame(y = all.forcing[iii1, 3] - mean(all.forcing[iii1, 
        3], na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.all4 <- data.frame(y = all.forcing[iii1, 4] - mean(all.forcing[iii1, 
        4], na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.all5 <- data.frame(y = all.forcing[iii1, 5] - mean(all.forcing[iii1, 
        5], na.rm = TRUE), x1 = forcings$Solar[iii2], x2 = forcings$W.M_GHGs[iii2], 
        x3 = forcings$O3[iii2], x4 = forcings$StratH2O[iii2], 
        x5 = forcings$LandUse[iii2], x6 = forcings$SnowAlb[iii2], 
        x7 = forcings$StratAer[iii2], x8 = forcings$BC[iii2], 
        x9 = forcings$ReflAer[iii2], x10 = forcings$AIE[iii2])
    CALIBRATE.ghg <- data.frame(y = ensemble.mean.ghg[iii1] - 
        mean(ensemble.mean.ghg[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2], x3 = forcings$O3[iii2], 
        x4 = forcings$StratH2O[iii2], x5 = forcings$LandUse[iii2], 
        x6 = forcings$SnowAlb[iii2], x7 = forcings$StratAer[iii2], 
        x8 = forcings$BC[iii2], x9 = forcings$ReflAer[iii2], 
        x10 = forcings$AIE[iii2])
    CALIBRATE.tsi <- data.frame(y = ensemble.mean.tsi[iii1] - 
        mean(ensemble.mean.tsi[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2], x3 = forcings$O3[iii2], 
        x4 = forcings$StratH2O[iii2], x5 = forcings$LandUse[iii2], 
        x6 = forcings$SnowAlb[iii2], x7 = forcings$StratAer[iii2], 
        x8 = forcings$BC[iii2], x9 = forcings$ReflAer[iii2], 
        x10 = forcings$AIE[iii2])
    CALIBRATE.obs.s.co2 <- data.frame(y = gistemp$T2m[iii4] - 
        mean(gistemp$T2m[iii4], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2])
    CALIBRATE.all.s.co2 <- data.frame(y = ensemble.mean.all[iii1] - 
        mean(ensemble.mean.all[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2])
    CALIBRATE.ghg.s.co2 <- data.frame(y = ensemble.mean.ghg[iii1] - 
        mean(ensemble.mean.all[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2])
    CALIBRATE.tsi.s.co2 <- data.frame(y = ensemble.mean.tsi[iii1] - 
        mean(ensemble.mean.all[iii1], na.rm = TRUE), x1 = forcings$Solar[iii2], 
        x2 = forcings$W.M_GHGs[iii2])
    CALIBRATE.obs.s.co2.1958 <- data.frame(y = gistemp$T2m[iv4] - 
        mean(gistemp$T2m[iv4], na.rm = TRUE), x1 = forcings$Solar[iv2], 
        x2 = forcings$W.M_GHGs[iv2])
    CALIBRATE.all.s.co2.1958 <- data.frame(y = ensemble.mean.all[iv1] - 
        mean(ensemble.mean.all[iv1], na.rm = TRUE), x1 = forcings$Solar[iv2], 
        x2 = forcings$W.M_GHGs[iv2])
    CALIBRATE.ghg.s.co2.1958 <- data.frame(y = ensemble.mean.ghg[iv1] - 
        mean(ensemble.mean.all[iv1], na.rm = TRUE), x1 = forcings$Solar[iv2], 
        x2 = forcings$W.M_GHGs[iv2])
    CALIBRATE.tsi.s.co2.1958 <- data.frame(y = ensemble.mean.tsi[iv1] - 
        mean(ensemble.mean.all[iv1], na.rm = TRUE), x1 = forcings$Solar[iv2], 
        x2 = forcings$W.M_GHGs[iv2])
    print("=======================================================")
    print("Multiple-regressnio with only S and CO2:")
    print("1880--2002")
    print(summary(lm(CALIBRATE.obs.s.co2)))
    print(summary(lm(CALIBRATE.all.s.co2)))
    print(summary(lm(CALIBRATE.ghg.s.co2)))
    print(summary(lm(CALIBRATE.tsi.s.co2)))
    print("1950--2000")
    print(summary(lm(CALIBRATE.obs.s.co2.1958)))
    print(summary(lm(CALIBRATE.all.s.co2.1958)))
    print(summary(lm(CALIBRATE.ghg.s.co2.1958)))
    print(summary(lm(CALIBRATE.tsi.s.co2.1958)))
    print("=======================================================")
    attr(calibrate.all, "Description") <- c("GISS GCM all forcings", 
        "Regression: TSI")
    attr(calibrate.ghg, "Description") <- c("GISS GCM GHG", "Regression:  TSI")
    attr(calibrate.res, "Description") <- c("GISS GCM residual (GHG-all)", 
        "Regression:  TSI")
    attr(calibrate.tsi, "Description") <- c("GISS GCM solar", 
        "Regression: TSI")
    attr(calibrate.res, "Description") <- c("Observations", "Regression: TSI from SW2006b")
    attr(Calibrate.obs, "Description") <- c("Observations", "Multiple regression: C02 & TSI")
    attr(Calibrate.all, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression: C02 & TSI")
    attr(Calibrate.ghg, "Description") <- c("GISS GCM GHG", "Multiple regression: C02 & TSI")
    attr(Calibrate.res, "Description") <- c("GISS GCM residual (GHG-all)", 
        "Multiple regression: C02 & TSI")
    attr(Calibrate.tsi, "Description") <- c("GISS GCM solar", 
        "Multiple regression: C02 & TSI")
    attr(Calibrate.res, "Description") <- c("Observations", "Multiple regression: C02 & TSI from SW2006b")
    attr(CALIBRATE.all, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.all1, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.all2, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.all3, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.all4, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.all5, "Description") <- c("GISS GCM all forcings", 
        "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.obs, "Description") <- c("Observations", "Multiple regression  with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.ghg, "Description") <- c("GISS GCM GHG", "Multiple regression with GISS forcings: C02, TSI, aerosols, O3, ++")
    attr(CALIBRATE.tsi, "Description") <- c("GISS GCM GHG", "Multiple regression  with GISS forcings: C02, TSI, aerosols, O3, ++")

    # Prepare forcing data for input to prediction experiments:
    Pred.obs.tsi <- Calibrate.obs
    Pred.obs.tsi$x1[] <- mean(Calibrate.obs$x1, na.rm = TRUE)
    Pred.obs.co2 <- Calibrate.obs
    Pred.obs.co2$x2[] <- mean(Calibrate.obs$x2, na.rm = TRUE)
    Pred.all.tsi <- Calibrate.all
    Pred.all.tsi$x1[] <- mean(Calibrate.all$x1, na.rm = TRUE)
    Pred.all.co2 <- Calibrate.all
    Pred.all.co2$x2[] <- mean(Calibrate.all$x2, na.rm = TRUE)
    Pred.res.tsi <- Calibrate.res
    Pred.res.tsi$x1[] <- mean(Calibrate.res$x1, na.rm = TRUE)
    Pred.res.co2 <- Calibrate.res
    Pred.res.co2$x2[] <- mean(Calibrate.res$x2, na.rm = TRUE)
    Pred.tsi.tsi <- Calibrate.tsi
    Pred.tsi.tsi$x1[] <- mean(Calibrate.tsi$x1, na.rm = TRUE)
    Pred.tsi.co2 <- Calibrate.tsi
    Pred.tsi.co2$x2[] <- mean(Calibrate.tsi$x2, na.rm = TRUE)
    Pred.ghg.tsi <- Calibrate.ghg
    Pred.ghg.tsi$x1[] <- mean(Calibrate.ghg$x1, na.rm = TRUE)
    Pred.ghg.co2 <- Calibrate.ghg
    Pred.ghg.co2$x2[] <- mean(Calibrate.ghg$x2, na.rm = TRUE)
    Pred.SW.tsi <- Calibrate.SW
    Pred.SW.tsi$x1[] <- mean(Calibrate.SW$x1, na.rm = TRUE)
    Pred.SW.co2 <- Calibrate.SW
    Pred.SW.co2$x2[] <- mean(Calibrate.SW$x2, na.rm = TRUE)
    PRED.OBS.tsi <- CALIBRATE.obs
    PRED.OBS.tsi$x1[] <- mean(CALIBRATE.obs$x1, na.rm = TRUE)
    PRED.OBS.co2 <- CALIBRATE.obs
    PRED.OBS.co2$x2[] <- mean(CALIBRATE.obs$x2, na.rm = TRUE)
    PRED.ALL.tsi <- CALIBRATE.all
    PRED.ALL.tsi$x1[] <- mean(CALIBRATE.all$x1, na.rm = TRUE)
    PRED.ALL.co2 <- CALIBRATE.all
    PRED.ALL.co2$x2[] <- mean(CALIBRATE.all$x2, na.rm = TRUE)
    PRED.GHG.tsi <- CALIBRATE.ghg
    PRED.GHG.tsi$x1[] <- mean(CALIBRATE.ghg$x1, na.rm = TRUE)
    PRED.GHG.co2 <- CALIBRATE.ghg
    PRED.GHG.co2$x2[] <- mean(CALIBRATE.ghg$x2, na.rm = TRUE)
    PRED.TSI.tsi <- CALIBRATE.tsi
    PRED.TSI.tsi$x1[] <- mean(CALIBRATE.tsi$x1, na.rm = TRUE)
    PRED.TSI.co2 <- CALIBRATE.tsi
    PRED.TSI.co2$x2[] <- mean(CALIBRATE.tsi$x2, na.rm = TRUE)
    PRED.OBS.sol.only <- CALIBRATE.obs
    PRED.OBS.sol.only$x2[] <- mean(CALIBRATE.obs$x2, na.rm = TRUE)
    PRED.OBS.sol.only$x3[] <- mean(CALIBRATE.obs$x3, na.rm = TRUE)
    PRED.OBS.sol.only$x4[] <- mean(CALIBRATE.obs$x4, na.rm = TRUE)
    PRED.OBS.sol.only$x5[] <- mean(CALIBRATE.obs$x5, na.rm = TRUE)
    PRED.OBS.sol.only$x6[] <- mean(CALIBRATE.obs$x6, na.rm = TRUE)
    PRED.OBS.sol.only$x7[] <- mean(CALIBRATE.obs$x7, na.rm = TRUE)
    PRED.OBS.sol.only$x8[] <- mean(CALIBRATE.obs$x8, na.rm = TRUE)
    PRED.OBS.sol.only$x9[] <- mean(CALIBRATE.obs$x9, na.rm = TRUE)
    PRED.OBS.sol.only$x10[] <- mean(CALIBRATE.obs$x10, na.rm = TRUE)
    PRED.OBS.ghg.only <- PRED.OBS.sol.only
    PRED.OBS.ghg.only$x1[] <- mean(CALIBRATE.obs$x1, na.rm = TRUE)
    PRED.OBS.ghg.only$x2 <- CALIBRATE.obs$x2
    PRED.ALL.sol.only <- CALIBRATE.all
    PRED.ALL.sol.only$x2[] <- mean(CALIBRATE.all$x2, na.rm = TRUE)
    PRED.ALL.sol.only$x3[] <- mean(CALIBRATE.all$x3, na.rm = TRUE)
    PRED.ALL.sol.only$x4[] <- mean(CALIBRATE.all$x4, na.rm = TRUE)
    PRED.ALL.sol.only$x5[] <- mean(CALIBRATE.all$x5, na.rm = TRUE)
    PRED.ALL.sol.only$x6[] <- mean(CALIBRATE.all$x6, na.rm = TRUE)
    PRED.ALL.sol.only$x7[] <- mean(CALIBRATE.all$x7, na.rm = TRUE)
    PRED.ALL.sol.only$x8[] <- mean(CALIBRATE.all$x8, na.rm = TRUE)
    PRED.ALL.sol.only$x9[] <- mean(CALIBRATE.all$x9, na.rm = TRUE)
    PRED.ALL.sol.only$x10[] <- mean(CALIBRATE.all$x10, na.rm = TRUE)
    PRED.ALL.ghg.only <- PRED.ALL.sol.only
    PRED.ALL.ghg.only$x1[] <- mean(CALIBRATE.all$x1, na.rm = TRUE)
    PRED.ALL.ghg.only$x2 <- CALIBRATE.all$x2
    NT <- length(Lean2004$X11yrCYCLE.BKGRND)
    PRED.obs.tsi <- data.frame(x1 = Lean2004$X11yrCYCLE.BKGRND, 
        x2 = rep(mean(Calibrate.obs$x2, na.rm = TRUE), NT))
    PRED.SW.tsi <- data.frame(x1 = Lean2004$X11yrCYCLE.BKGRND, 
        x2 = rep(mean(Calibrate.obs$x2, na.rm = TRUE), NT))
    PRED.all.tsi <- data.frame(x1 = Lean2004$X11yrCYCLE.BKGRND, 
        x2 = rep(mean(Calibrate.all$x2, na.rm = TRUE), NT))
    Obs.tsi <- lm(y ~ x1 + x2, data = Calibrate.obs)
    obs.tsi <- lm(y ~ x1, data = Calibrate.obs)
    Obs.SW <- lm(y ~ x1 + x2, data = Calibrate.SW)
    obs.SW <- lm(y ~ x, data = calibrate.SW)
    all.tsi <- lm(y ~ x, data = calibrate.all)
    All.tsi <- lm(y ~ x1 + x2, data = Calibrate.all)
    ghg.tsi <- lm(y ~ x, data = calibrate.ghg)
    Ghg.tsi <- lm(y ~ x1 + x2, data = Calibrate.ghg)
    res.tsi <- lm(y ~ x, data = calibrate.res)
    Res.tsi <- lm(y ~ x1 + x2, data = Calibrate.res)
    sol.tsi <- lm(y ~ x, data = calibrate.tsi)
    Sol.tsi <- lm(y ~ x1 + x2, data = Calibrate.tsi)
    OBS.forcings <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.obs)
    if (stepwise) 
        OBS.forcings <- step(OBS.forcings, trace = 0)
    DW.obs <- dwtest(OBS.forcings)
    ALL.forcings <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all)
    if (stepwise) 
        ALL.forcings <- step(ALL.forcings, trace = 0)
    ALL.forcings1 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all1)
    if (stepwise) 
        ALL.forcings1 <- step(ALL.forcings1, trace = 0)
    ALL.forcings2 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all2)
    if (stepwise) 
        ALL.forcings2 <- step(ALL.forcings2, trace = 0)
    ALL.forcings3 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all3)
    if (stepwise) 
        ALL.forcings3 <- step(ALL.forcings3, trace = 0)
    ALL.forcings4 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all4)
    if (stepwise) 
        ALL.forcings4 <- step(ALL.forcings4, trace = 0)
    ALL.forcings5 <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.all5)
    if (stepwise) 
        ALL.forcings5 <- step(ALL.forcings5, trace = 0)
    GHG.forcings <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.ghg)
    if (stepwise) 
        GHG.forcings <- step(GHG.forcings, trace = 0)
    TSI.forcings <- lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + 
        x8 + x9 + x10, data = CALIBRATE.tsi)
    if (stepwise) 
        TSI.forcings <- step(TSI.forcings, trace = 0)
    pred.obs.tsi <- predict(Obs.tsi, newdata = Pred.obs.tsi)
    pred.obs.co2 <- predict(Obs.tsi, newdata = Pred.obs.co2)
    pred.SW.tsi <- predict(Obs.SW, newdata = Pred.SW.tsi)
    pred.SW.co2 <- predict(Obs.SW, newdata = Pred.SW.co2)
    pred.all.tsi <- predict(All.tsi, newdata = Pred.all.tsi)
    pred.all.co2 <- predict(All.tsi, newdata = Pred.all.co2)
    pred.res.tsi <- predict(Res.tsi, newdata = Pred.res.tsi)
    pred.res.co2 <- predict(Res.tsi, newdata = Pred.res.co2)
    pred.tsi.tsi <- predict(Sol.tsi, newdata = Pred.tsi.tsi)
    pred.tsi.co2 <- predict(Sol.tsi, newdata = Pred.tsi.co2)
    pred.ghg.tsi <- predict(Ghg.tsi, newdata = Pred.ghg.tsi)
    pred.ghg.co2 <- predict(Ghg.tsi, newdata = Pred.ghg.co2)
    pred.obs.TSI <- predict(Obs.tsi, newdata = PRED.obs.tsi)
    pred.all.TSI <- predict(All.tsi, newdata = PRED.all.tsi)
    pred.SW.TSI <- predict(All.tsi, newdata = PRED.SW.tsi)
    pre.OBS.TSI <- predict(OBS.forcings, newdata = PRED.OBS.tsi)
    pre.OBS.CO2 <- predict(OBS.forcings, newdata = PRED.OBS.co2)
    pre.ALL.TSI <- predict(ALL.forcings, newdata = PRED.ALL.tsi)
    pre.ALL.CO2 <- predict(ALL.forcings, newdata = PRED.ALL.co2)
    pre.GHG.TSI <- predict(GHG.forcings, newdata = PRED.GHG.tsi)
    pre.GHG.CO2 <- predict(GHG.forcings, newdata = PRED.GHG.co2)
    pre.TSI.TSI <- predict(TSI.forcings, newdata = PRED.TSI.tsi)
    pre.TSI.CO2 <- predict(TSI.forcings, newdata = PRED.TSI.co2)
    pre.OBS.sol.only <- predict(OBS.forcings, newdata = PRED.OBS.sol.only)
    pre.ALL.sol.only <- predict(ALL.forcings, newdata = PRED.ALL.sol.only)
    pre.OBS.ghg.only <- predict(OBS.forcings, newdata = PRED.OBS.ghg.only)
    pre.ALL.ghg.only <- predict(ALL.forcings, newdata = PRED.ALL.ghg.only)
    T.e <- T.eq(Lean2004$X11yrCYCLE.BKGRND)
    pred.obs.TSI <- pred.obs.TSI - mean(pred.obs.TSI[is.element(trunc(Lean2004$YEAR), 
        1885:1914)], na.rm = TRUE)
    pred.SW.TSI <- pred.SW.TSI - mean(pred.SW.TSI[is.element(trunc(Lean2004$YEAR), 
        1885:1914)], na.rm = TRUE)
    pred.all.TSI <- pred.all.TSI - mean(pred.all.TSI[is.element(trunc(Lean2004$YEAR), 
        1885:1914)], na.rm = TRUE)
    T.e <- T.e - mean(T.e[is.element(trunc(Lean2004$YEAR), 1885:1914)], 
        na.rm = TRUE)
    iv <- is.element(gistemp$year, 1901:2000)
    iiv <- is.element(ensemble.year, 1901:2000)
    linear.obs <- data.frame(y = gistemp$T2m[iv] - mean(gistemp$T2m[iii4], 
        na.rm = TRUE), x = 1901:2000)
    linear.all <- data.frame(y = ensemble.mean.all[iiv] - mean(ensemble.mean.all[iii1], 
        na.rm = TRUE), x = 1901:2000)
    linear.tsi <- data.frame(y = ensemble.mean.tsi[iiv] - mean(ensemble.mean.tsi[iii1], 
        na.rm = TRUE), x = 1901:2000)
    iv <- is.element(gistemp$year, 1980:2000)
    iiv <- is.element(ensemble.year, 1980:2000)
    linear.all2 <- data.frame(y = ensemble.mean.all[iiv] - mean(ensemble.mean.all[iii1], 
        na.rm = TRUE), x = 1980:2000)
    linear.tsi2 <- data.frame(y = ensemble.mean.tsi[iiv] - mean(ensemble.mean.tsi[iii1], 
        na.rm = TRUE), x = 1980:2000)
    TREND.OBS <- predict(lm(y ~ x, data = linear.obs))
    TREND.ALL <- predict(lm(y ~ x, data = linear.all))
    TREND.TSI <- predict(lm(y ~ x, data = linear.tsi))
    TREND.TSI2 <- predict(lm(y ~ x, data = linear.tsi2))
    twentiethC <- is.element(trunc(Lean2004$YEAR), 1900:1999)
    lin.trend.obs <- data.frame(y = pred.obs.TSI[twentiethC], 
        x = 1900:1999)
    lin.trend.SW <- data.frame(y = pred.SW.TSI[twentiethC], x = 1900:1999)
    lin.trend.all <- data.frame(y = pred.all.TSI[twentiethC], 
        x = 1900:1999)
    lin.trend.Te <- data.frame(y = T.e[twentiethC], x = 1900:1999)
    trend.obs <- predict(lm(y ~ x, data = lin.trend.obs))
    trend.SW <- predict(lm(y ~ x, data = lin.trend.SW))
    trend.all <- predict(lm(y ~ x, data = lin.trend.all))
    trend.Te <- predict(lm(y ~ x, data = lin.trend.Te))
    print("Details about the regression models:")
    print(anova(OBS.forcings))
    print(anova(ALL.forcings))

    # Same exercise with the regression analysis for the GISS forcing data:

    print("1901--2000 trend estimates for multiple regression models:")
    TwentiethC <- is.element(trunc(forcings$Year[iii1]), 
        1900:1999)
    LIN.trend.obs <- data.frame(y = pre.OBS.TSI[TwentiethC], 
        x = 1900:1999)
    LIN.trend.all <- data.frame(y = pre.ALL.TSI[TwentiethC], 
        x = 1900:1999)
    LIN.trend.all <- data.frame(y = pre.ALL.TSI[TwentiethC], 
        x = 1900:1999)

    # pred.obs.tsi
    # Bivariate regression model:
    if (bivariate) {
        LIN.trend.exp.obs <- data.frame(y = pred.obs.tsi, x = ensemble.year[ii1])
        LIN.trend.exp.all <- data.frame(y = pred.all.tsi, x = ensemble.year[ii1])
        LIN.trend.ghg.obs <- data.frame(y = pred.obs.co2, x = ensemble.year[ii1])
        LIN.trend.ghg.all <- data.frame(y = pred.all.co2, x = ensemble.year[ii1])
    }
    else {
      # Multiple regression model:
        LIN.trend.exp.obs <- data.frame(y = pre.OBS.sol.only[TwentiethC], 
            x = 1900:1999)
        LIN.trend.exp.all <- data.frame(y = pre.ALL.sol.only[TwentiethC], 
            x = 1900:1999)
        LIN.trend.ghg.obs <- data.frame(y = pre.OBS.ghg.only[TwentiethC], 
            x = 1900:1999)
        LIN.trend.ghg.all <- data.frame(y = pre.ALL.ghg.only[TwentiethC], 
            x = 1900:1999)
    }
    print("Regression:")
    c.obs <- summary(lm(y ~ x, data = linear.obs))$coefficients
    c.all <- summary(lm(y ~ x, data = linear.all))$coefficients
    c.tsi <- summary(lm(y ~ x, data = linear.tsi))$coefficients
    c.all2 <- summary(lm(y ~ x, data = linear.all2))$coefficients
    c.tsi2 <- summary(lm(y ~ x, data = linear.tsi2))$coefficients
    
    plot(linear.obs$x, linear.obs$y, type = "b", pch = 19)
    grid()
    abline(lm(y ~ x, data = linear.obs), lty = 2)
    lines(linear.all$x, linear.all$y, col = "red")
    lines(linear.tsi$x, linear.tsi$y, col = "blue")
    abline(lm(y ~ x, data = linear.all), lty = 2, col = "red")
    abline(lm(y ~ x, data = linear.tsi), lty = 2, col = "blue")
    
    plot(LIN.trend.obs$x, LIN.trend.obs$y, type = "b", pch = 19, 
        ylim = c(-0.5, 2))
    grid()
    points(LIN.trend.all$x, LIN.trend.all$y, col = "red")
    lines(LIN.trend.exp.obs$x, LIN.trend.exp.obs$y, col = "blue")
    lines(LIN.trend.exp.all$x, LIN.trend.exp.all$y, col = "darkgreen")
    lines(LIN.trend.ghg.obs$x, LIN.trend.ghg.obs$y, col = "darkblue")
    lines(LIN.trend.ghg.all$x, LIN.trend.ghg.all$y, col = "darkred", 
        lty = 2)
    
    print("Masking solar forcing")
    TREND.year <- 1901:2000
    TREND.obs <- predict(lm(y ~ x, data = LIN.trend.obs))
    TREND.all <- predict(lm(y ~ x, data = LIN.trend.all))
    c.obs.co2.only <- summary(lm(y ~ x, data = LIN.trend.ghg.obs))$coefficients
    c.all.co2.only <- summary(lm(y ~ x, data = LIN.trend.ghg.all))$coefficients
    c.obs.sol.only <- summary(lm(y ~ x, data = LIN.trend.exp.obs))$coefficients
    c.all.sol.only <- summary(lm(y ~ x, data = LIN.trend.exp.all))$coefficients
    dT.0S.obs <- round(100 * c.obs.co2.only[2]/c.obs[2], 4)
    dT.0S.all <- round(100 * c.all.co2.only[2]/c.all[2], 4)
    dT.tsi <- round(100 * c.tsi[2]/c.all[2], 4)
    dT.tsi2 <- round(100 * c.tsi2[2]/c.all2[2], 4)
    dT.dS.obs <- round(100 * c.obs.sol.only[2]/c.obs[2], 2)
    dT.dS.all <- round(100 * c.all.sol.only[2]/c.all[2], 2)
    dT.0S.obs.err <- 100 * sqrt((c.obs.sol.only[4]/c.obs.sol.only[2])^2 + 
        (c.obs[4]/c.obs[2])^2)
    dT.0S.all.err <- 100 * sqrt((c.all.sol.only[4]/c.all.sol.only[2])^2 + 
        (c.all[4]/c.all[2])^2)
    dT.tsi.err <- 100 * sqrt((c.tsi[4]/c.tsi[2])^2 + (c.all[4]/c.all[2])^2)
    dT.tsi2.err <- 100 * sqrt((c.tsi2[4]/c.tsi2[2])^2 + (c.all2[4]/c.all2[2])^2)
    print(paste("Obs. ", round(min(TREND.obs, na.rm = TRUE), 
        4), "-", round(max(TREND.obs, na.rm = TRUE), 4), "which is", 
        dT.0S.obs, " pm ", dT.0S.obs.err, "%"))
    print(paste("All. ", round(min(TREND.all, na.rm = TRUE), 
        4), "-", round(max(TREND.all, na.rm = TRUE), 4), "which is", 
        dT.0S.all, " pm ", dT.0S.all.err, "%"))
    print(paste("1901-2000 change in observations: ", paste(diff(range(TREND.obs))[1])))
    print(paste("1901-2000 change in observations: ", paste(diff(range(TREND.all))[1])))
    print("Solar forcing only")
    TREND.obs.sol.only <- predict(lm(y ~ x, data = LIN.trend.exp.obs), 
        newdata = LIN.trend.exp.obs)
    TREND.all.sol.only <- predict(lm(y ~ x, data = LIN.trend.exp.all), 
        newdata = LIN.trend.exp.all)
    print(paste("Obs. ", round(min(TREND.obs.sol.only, na.rm = TRUE), 
        4), "-", round(max(TREND.obs.sol.only, na.rm = TRUE), 
        4), "which is", dT.dS.obs, " pm ", dT.0S.obs.err, "%"))
    print(paste("All. ", round(min(TREND.all.sol.only, na.rm = TRUE), 
        4), "-", round(max(TREND.all.sol.only, na.rm = TRUE), 
        4), "which is", dT.dS.all, " pm ", dT.0S.all.err, "%"))
    print("Trend in 'sol'")
    print(paste("'sol' 1901-2000: ", round(min(TREND.TSI, na.rm = TRUE), 
        4), "-", round(max(TREND.TSI, na.rm = TRUE), 4), "which is", 
        dT.tsi, " pm ", dT.tsi.err, "% compared to 'all'"))
    print(paste("'sol' 1980-2000: ", round(min(TREND.TSI2, na.rm = TRUE), 
        4), "-", round(max(TREND.TSI2, na.rm = TRUE), 4), "which is", 
        dT.tsi2, " pm ", dT.tsi2.err, "% compared to 'all'"))

    #Numbers, tables, etc.:---------------------------------------------------------------
    if (tables) {
        print("###################### Results ##############################")
        print(c(round(100 * trend1$coefficients[2]/trend0$coefficients[2]), 
            round(100 * trend2$coefficients[2]/trend0$coefficients[2]), 
            round(100 * trend3$coefficients[2]/trend0$coefficients[2])))
        print(c(round(100 * summary(trend1)$coefficients[4]/trend0$coefficients[2], 
            1), round(100 * summary(trend2)$coefficients[4]/trend0$coefficients[2], 
            1), round(100 * summary(trend3)$coefficients[4]/trend0$coefficients[2], 
            1)))
        print(summary(res.tsi))
        print(summary(sol.tsi))
        print(summary(all.tsi))
        print(summary(Res.tsi))
        print(summary(Sol.tsi))
        print(summary(All.tsi))
        tab.results <- rep("  ", 6)
        tab.results[1] <- make.table1(obs.tsi, Obs.tsi, "obs. & ")
        tab.results[2] <- make.table1(all.tsi, All.tsi, "All. & ")
        tab.results[3] <- make.table1(sol.tsi, Sol.tsi, "Sol. & ")
        tab.results[4] <- make.table1(res.tsi, Res.tsi, "Res. & ")
        tab.results[5] <- make.table1(ghg.tsi, Ghg.tsi, "GHG. & ")
        tab.results[6] <- make.table1(obs.SW, Obs.SW, "SW   & ")
        print(paste("lag=", lag))
        print(tab.results, quote = FALSE)

        #Table for the regression analysis based on GISS forcings:
        TAB.results <- rep("  ", 9)
        TAB.results[1] <- make.table2(OBS.forcings, "obs. & ")
        TAB.results[2] <- make.table2(ALL.forcings, "all. & ")
        TAB.results[3] <- make.table2(ALL.forcings1, "all1 & ")
        TAB.results[4] <- make.table2(ALL.forcings2, "all2 & ")
        TAB.results[5] <- make.table2(ALL.forcings3, "all3 & ")
        TAB.results[6] <- make.table2(ALL.forcings4, "all4 & ")
        TAB.results[7] <- make.table2(ALL.forcings5, "all5 & ")
        TAB.results[8] <- make.table2(GHG.forcings, "ghg.& ")
        TAB.results[9] <- make.table2(TSI.forcings, "sol.& ")
        if (!stepwise) 
            print(TAB.results, quote = FALSE)
        else {
            print("Obs:")
            print(summary(OBS.forcings))
            print("All:")
            print(summary(ALL.forcings))
            print("All1:")
            print(summary(ALL.forcings1))
            print("All2:")
            print(summary(ALL.forcings2))
            print("All3:")
            print(summary(ALL.forcings3))
            print("All4:")
            print(summary(ALL.forcings4))
            print("All5:")
            print(summary(ALL.forcings5))
            print("GHG:")
            print(summary(GHG.forcings))
            print("TSI:")
            print(summary(TSI.forcings))
        }
        print(DW.obs)
    }

    #Graphics:---------------------------------------------------------------
    t.rng <- range(c(gis.res, ensemble.mean.all[i1] - mean(ensemble.mean.all[i1], 
        na.rm = TRUE), ensemble.mean.tsi[i1] - mean(ensemble.mean.tsi[i1], 
        na.rm = TRUE), gistemp$T2m[ii4] - mean(gistemp$T2m[ii4], 
        na.rm = TRUE)), na.rm = TRUE)
    if (figures) {
        x11()
        plot(gistemp$year[iii4], gistemp$T2m[iii4] - mean(gistemp$T2m[iii4], 
            na.rm = TRUE), main = "Reconstruction of <T> by linear models", 
            pch = 19, col = "red", xlab = "year", ylab = "Temperature anomaly (K)", 
            ylim = t.rng, type = "l", lwd = 5,
             sub = "Reconstruction by multiple regression & GISS forcings")
        lines(ensemble.year, ensemble.mean.all - mean(ensemble.mean.all), 
            col = "blue", lwd = 5)
        grid()
        lines(gistemp$year[iii4], predict(OBS.forcings), lwd = 1, 
            lty = 1, col = "darkred", type = "b", pch = 19)
        lines(ensemble.year[iii1], predict(ALL.forcings), lwd = 1, 
            lty = 1, col = "darkblue", type = "b", pch = 19)
        lines(gistemp$year[ii4], predict(Obs.tsi) + mean(gistemp$T2m[ii4], 
            na.rm = TRUE), lwd = 1, lty = 1, col = "pink", type = "b")
        lines(ensemble.year[ii4], predict(All.tsi) + mean(ensemble.mean.all[ii1], 
            na.rm = TRUE) - mean(ensemble.mean.all), lwd = 1, 
            lty = 1, col = "steelblue", type = "b")
        points(gistemp$year[ii4], predict(Obs.tsi) + mean(gistemp$T2m[ii4], 
            na.rm = TRUE), col = "pink", pch = 21)
        points(ensemble.year[ii4], predict(All.tsi) + mean(ensemble.mean.all[ii1], 
            na.rm = TRUE) - mean(ensemble.mean.all), col = "steelblue", 
            pch = 21)
        lines(ensemble.year, ensemble.mean.tsi - mean(ensemble.mean.tsi), 
            col = "grey", lty = 2)
        grid()
        legend(min(gistemp$year[iii4], na.rm = TRUE), max(gistemp$T2m[iii4] - 
            mean(gistemp$T2m[iii4], na.rm = TRUE), na.rm = TRUE) * 
            0.85, c("Obs  ", "all  ", "Eq1 - obs ", "Eq1 - all  ", 
            "Eq2 - obs   ", "Eq2 - all   ", "solar   "), col = c("red", 
            "blue", "pink", "steelblue", "darkred", "darkblue", 
            "grey"), lwd = c(5, 5, 1, 1, 1, 1, 1), lty = c(1, 
            1, 0, 0, 0, 0, 2), pch = c(26, 26, 21, 21, 20, 20, 
            26), bg = "grey95", cex = 0.8)
        print("HERE2")

        x11()
        plot(gistemp$year, gistemp$T2m - gistemp$level.1900, 
            pch = 17, col = "grey", xlab = "Year",
            ylab = "global mean temperature anomaly (Celsius)", 
            main = "Measured and calculated global mean temperature", 
            xlim = c(1900, 2000))
        polygon(c(1900, 2000, 2000, 1900, 1900), c(0, 0, 0.6, 
            0.6, 0), col = "grey92", border = "grey80")
        points(gistemp$year, gistemp$T2m - gistemp$level.1900, 
            pch = 17, col = "grey")
        grid()
        lines(Z0$x, predict(trend0), lty = 2, col = "grey")
        lines(Z.S$x, predict(trend1), lty = 2)
        lines(Z.SW$x, predict(trend2), lty = 2, col = "grey40")
        lines(Z.all$x, predict(trend3), lty = 2, col = "grey20")
        lines(gistemp$year, gistemp$T2m - gistemp$level.1900, 
            lty = 3, col = "grey")
        lines(gistemp$year, gistemp$X5.year.mean - gistemp$level.1900, 
            lwd = 5, col = "grey")
        lines(x, T.sun(S.4, D.4, D.3), lwd = 2)
        lines(x.SW, T.sun(S.4.SW, D.4.SW, D.3.SW), lty = 2, col = "grey40")
        lines(x.all, T.sun(S.4.all, D.4.all, D.3.all), lty = 5, col="grey20")
        legend(1970, 0, c("SW06a             ", "Lean2004            ", 
            "adjust all years", T.descr), col = c("grey40", "black", 
            "grey20", "grey"), lwd = c(1, 2, 1, 5),
               lty = c(2, 1, 5, 1), bg = "grey95", cex = 0.8)
    }

    #----- Auxiliary figures:--------------------------------------------------------

    print("Auxiliary figures:")
    x11()
    par(col.axis = "white", col.lab = "white")
    plot(c(0, 1), c(0, 1), type = "n")
    text(0.5, 0.5, "Auxiliary Figures", cex = 1.7)
    x11()
    ar1 <- ccf(gistemp$T2m[ii4], Lean2004$X11yrCYCLE.BKGRND[ii2])
    plot(ar1, sub = paste("Obs <T>:", lag, "-year delay"), lwd = 5, 
        main = "Observed <T> & Lean (2004) S", ylab = "Cross-correlation", 
        xlab = "Lag (years)")
    grid()
    x11()
    ar2 <- ccf(ensemble.mean.all[ii1], Lean2004$X11yrCYCLE.BKGRND[ii2])
    plot(ar2, sub = paste("'all' <T>:", lag, "-year delay"), 
        lwd = 5, main = "'all' <T> & Lean (2004) S", ylab = "Cross-correlation", 
        xlab = "Lag (years)")
    grid()
    x11()
    ar3 <- ccf(ensemble.mean.ghg[ii1], Lean2004$X11yrCYCLE.BKGRND[ii2])
    plot(ar3, sub = paste("'GCM' <T>:", lag, "-year delay"), 
        lwd = 5, main = "'GHG' <T> & Lean (2004) S", ylab = "Cross-correlation", 
        xlab = "Lag (years)")
    grid()
    x11()
    ar4 <- ccf(ensemble.mean.tsi[ii1], Lean2004$X11yrCYCLE.BKGRND[ii2])
    plot(ar4, sub = paste("GCM <T>:", lag, "-year delay"), lwd = 5, 
        main = "'sol' <T> & Lean (2004) S", ylab = "Cross-correlation", 
        xlab = "Lag (years)")
    grid()
    x11()
    ar4 <- ccf(y.SW, y)
    plot(ar4, sub = paste("S&W <T>:", lag, "-year delay"), lwd = 5, 
        main = "'sol' <T> & Lean (2004) S", ylab = "Cross-correlation", 
        xlab = "Lag (years)")
    grid()
    x11()
    plot(ensemble.year[i1], ensemble.mean.tsi[i1] - mean(ensemble.mean.tsi[i1], 
        na.rm = TRUE), type = "l", ylim = t.rng, main = "<T> response to TSI", 
        lwd = 2, ylab = "<T>", xlab = "year")
    grid()
    lines(ensemble.year[i1], ensemble.mean.all[i1] - mean(ensemble.mean.all[i1], 
        na.rm = TRUE), col = "blue", lwd = 2)
    lines(ensemble.year[i1], ensemble.mean.all[i1] - ensemble.mean.ghg[i1], 
        col = "red", lwd = 2)
    legend(1880, 0.45, c("Solar only", "All - GHG", "All"), col = c("black", 
        "red", "blue"), lwd = 2, bg = "grey95")
    x11()
    plot(gistemp$year[ii4], gistemp$T2m[ii4] - mean(gistemp$T2m[ii4], 
        na.rm = TRUE), main = "Decomposition of <T>", pch = 19, 
        col = "grey", xlab = "year", ylab = "Temperature anomaly (K)", 
        ylim = t.rng, type = "l", lwd = 3)
    grid()
    lines(gistemp$year[ii4], pred.obs.co2, col = "grey70", lty = 2)
    lines(gistemp$year[ii4], pred.all.co2, col = "lightblue", 
        lty = 2)
    lines(gistemp$year[ii4], pred.res.co2, col = "pink", lty = 2)
    lines(gistemp$year[ii4], pred.tsi.co2, col = "grey30", lty = 2)
    lines(gistemp$year[ii4], pred.obs.tsi, col = "grey70", lty = 1)
    lines(gistemp$year[ii4], pred.all.tsi, col = "lightblue", 
        lty = 1)
    lines(gistemp$year[ii4], pred.res.tsi, col = "pink", lty = 1)
    lines(gistemp$year[ii4], pred.tsi.tsi, col = "grey30", lty = 1)
    lines(Lean2004$YEAR[ii2], T.eq(Lean2004$X11yrCYCLE.BKGRND[ii2]) - 
        mean(T.eq(Lean2004$X11yrCYCLE.BKGRND[ii2]), na.rm = TRUE), 
        lty = 3, col = "red")
    legend(1955, 0.45, c("Obs", "all", "solar", "GHG-solar", 
        "T_e"), col = c("grey70", "lightblue", "pink", "grey20", 
        "red"), lwd = c(3, rep(1, 4)), lty = c(rep(1, 4), 3), 
        bg = "grey95")
    legend(1955, -0.37, c("CO2", "TSI"), lty = c(1, 2), lwd = c(1, 
        1), cex = 0.8, bg = "grey95")
    x11()
    plot(range(Lean2004$YEAR, na.rm = TRUE), range(c(pred.obs.TSI, 
        pred.all.TSI, T.e), na.rm = TRUE), main = "TSI reconstructions", 
        lwd = 3, col = "green", type = "n", ylab = "<T> anomaly", 
        xlab = "Time", sub = "multiple regression")
    grid()
    lines(Lean2004$YEAR, pred.obs.TSI, col = "blue", lwd = 2)
    lines(Lean2004$YEAR, pred.all.TSI, col = "red", lwd = 2)
    lines(Lean2004$YEAR, pred.SW.TSI, col = "steelblue", lwd = 1)
    lines(Lean2004$YEAR, T.e, col = "black", lty = 2)
    lines(Lean2004$YEAR[twentiethC], trend.obs, col = "blue", 
        lty = 3)
    lines(Lean2004$YEAR[twentiethC], trend.SW, col = "lightblue", 
        lty = 3)
    lines(Lean2004$YEAR[twentiethC], trend.all, col = "red", 
        lty = 3)
    lines(Lean2004$YEAR[twentiethC], trend.Te, col = "black", 
        lty = 3)
    dT.obs <- round(max(trend.obs) - min(trend.obs), 2)
    dT.SW <- round(max(trend.SW) - min(trend.SW), 2)
    dT.all <- round(max(trend.all) - min(trend.all), 2)
    dT.Te <- round(max(trend.Te) - min(trend.Te), 2)
    grid()
    legend(min(Lean2004$YEAR, na.rm = TRUE), max(c(pred.obs.TSI, 
        pred.all.TSI, T.e), na.rm = TRUE), c(paste("Obs.", dT.obs), 
        paste("S(SW06a)", dT.SW), paste("All", dT.all), paste("T_e", 
            dT.Te)), col = c("blue", "steelblue", "red", "black"), 
        lwd = c(2, 1, 2, 1), lty = c(1, 1, 2, 2), bg = "grey95")
    x11()
    plot(gistemp$year[iii4], gistemp$T2m[iii4] - mean(gistemp$T2m[iii4], 
        na.rm = TRUE), main = "Decomposition of <T> by ignoring CO2", 
        pch = 19, col = "red", xlab = "year", ylab = "Temperature anomaly (K)", 
        ylim = t.rng, type = "l", lwd = 1, sub = "Reconstruction by multiple regression & GISS forcings")
    lines(ensemble.year, ensemble.mean.all - mean(ensemble.mean.all), 
        col = "blue", lwd = 1)
    grid()
    lines(gistemp$year[iii4], predict(OBS.forcings), lty = 2, 
        lwd = 1, col = "darkred", type = "b", pch = 19)
    lines(ensemble.year[iii1], predict(ALL.forcings), lty = 2, 
        lwd = 1, col = "darkblue", type = "b", pch = 19)
    lines(ensemble.year[iii1], pre.OBS.CO2, col = "black", lwd = 2)
    lines(ensemble.year[iii1], pre.ALL.CO2, col = "grey", lwd = 2)
    lines(ensemble.year[iii1], pre.GHG.CO2, col = "steelblue", 
        lwd = 2)
    lines(ensemble.year[iii1], pre.TSI.CO2, col = "darkgreen", 
        lwd = 2)
    grid()
    legend(min(gistemp$year[iii4], na.rm = TRUE), max(gistemp$T2m[iii4] - 
        mean(gistemp$T2m[iii4], na.rm = TRUE), na.rm = TRUE) * 
        0.85, c("observations    ", "all             ", "Obs. eq.2         ", 
        "All eq.2          ", "Obs. mean(GHG)  ", "All mean(GHG)   ", 
        "GHG mean(GHG)     ", "solar mean(GHG)   "), col = c("red", 
        "blue", "darkred", "darkblue", "black", "grey", "steelblue", 
        "darkgreen"), lwd = c(1, 1, 1, 1, 2, 2, 2, 2), lty = c(1, 
        1, 2, 2, 1, 1, 1, 1), , pch = c(26, 26, 19, 19, 26, 26, 
        26, 26), bg = "grey95", cex = 0.6)
    x11()
    plot(gistemp$year[iii4], gistemp$T2m[iii4] - mean(gistemp$T2m[iii4], 
        na.rm = TRUE), main = "Decomposition of <T> by ignoring TSI", 
        pch = 19, col = "red", xlab = "year", ylab = "Temperature anomaly (K)", 
        ylim = t.rng, type = "l", lwd = 1, sub = "Reconstruction by multiple regression & GISS forcings")
    lines(ensemble.year, ensemble.mean.all - mean(ensemble.mean.all), 
        col = "blue", lwd = 1)
    grid()
    lines(gistemp$year[iii4], predict(OBS.forcings), lty = 2, 
        lwd = 1, col = "darkred", type = "b", pch = 19)
    lines(ensemble.year[iii1], predict(ALL.forcings), lty = 2, 
        lwd = 1, col = "darkblue", type = "b", pch = 19)
    lines(ensemble.year[iii1], pre.OBS.TSI, col = "black", lwd = 2)
    lines(ensemble.year[iii1], pre.ALL.TSI, col = "grey", lwd = 2)
    lines(ensemble.year[iii1], pre.GHG.TSI, col = "steelblue", 
        lwd = 2)
    lines(ensemble.year[iii1], pre.TSI.TSI, col = "darkgreen", 
        lwd = 2)
    grid()
    legend(min(gistemp$year[iii4], na.rm = TRUE), max(gistemp$T2m[iii4] - 
        mean(gistemp$T2m[iii4], na.rm = TRUE), na.rm = TRUE) * 
        0.85, c("observations    ", "all             ", "Obs. eq.2       ", 
        "All eq.2          ", "Obs. mean(F_s)  ", "All mean(F_s)   ", 
        "GHG mean(F_s)   ", "solar mean(F_s)   "), col = c("red", 
        "blue", "darkred", "darkblue", "black", "grey", "steelblue", 
        "darkgreen"), lwd = c(1, 1, 1, 1, 2, 2, 2, 2), lty = c(1, 
        1, 2, 2, 1, 1, 1, 1), , pch = c(26, 26, 19, 19, 26, 26, 
        26, 26), bg = "grey95", cex = 0.6)
    x11()
    plot(range(Lean2004$X11yrCYCLE.BKGRND[ii2]), t.rng, type = "n", 
        col = "grey", xlab = "TSI", ylab = "<T>", main = paste(T.descr, 
            "response to delta TSI"), sub = paste(min(trunc(Lean2004$YEAR[ii2])), 
            "-", max(trunc(Lean2004$YEAR[ii2]))))
    grid()
    points(Lean2004$X11yrCYCLE.BKGRND[ii2], gistemp$T2m[ii4] - 
        mean(gistemp$T2m[ii4], na.rm = TRUE), pch = 19, col = "grey", 
        cex = 1.2)
    points(Lean2004$X11yrCYCLE.BKGRND[i2], ensemble.mean.all[i1] - 
        mean(ensemble.mean.all[i1], na.rm = TRUE), pch = 19, 
        col = "black")
    points(Lean2004$X11yrCYCLE.BKGRND[i2], gis.res, pch = 2, 
        cex = 0.6, col = "grey70")
    points(Lean2004$X11yrCYCLE.BKGRND[i2], ensemble.mean.tsi[i1] - 
        mean(ensemble.mean.tsi[i1], na.rm = TRUE), col = "grey10")
    abline(obs.tsi, col = "grey", lty = 1)
    abline(all.tsi, col = "black", lty = 1)
    abline(res.tsi, col = "grey70", lty = 3)
    abline(sol.tsi, col = "grey10", lty = 2)
    legend(1366.5, -0.35, c("Solar ", "All-GHG  ", "All  ", "obs  "), 
        col = c("grey10", "grey70", "black", "grey"), pch = c(21, 
            2, 19, 19), lty = c(2, 3, 1, 1), cex = 0.75, bg = "grey95")
    x11()
    plot(x, Y, type = "l", lty = 1, lwd = 3, xlim = c(1900, 2000), 
        ylim = c(-0.05, 0.45), main = "Sensitivity to parameters", 
        xlab = "Time", ylab = "T")
    grid()
    legend(1900, 0.45, c("Lean (2004)       ", "SW06a        ", 
        "all         ", "Z_S4        "), lwd = c(3, 3, 3, 1), 
        lty = c(1, 2, 1, 2), cex = 0.8, col = c("black", "grey30", 
            "grey60", "black"), bg = "grey95")
    lines(x.all, Y.all, col = "grey60", lwd = 3)
    lines(x.SW, Y.SW, lty = 2, col = "grey30", lwd = 3)
    lines(x, Y.Zeq010, lty = 2, col = "grey40")
    lines(x, Y.Zeq000, lty = 2, col = "grey40")
    lines(x, Y.Z4005, lty = 2)
    lines(x, Y.Z4030, lty = 2)
    lines(x, Y.22y005, lty = 3, col = "grey40")
    lines(x, Y.22y030, lty = 3, col = "grey40")
    lines(x, Y.11y005, lty = 3, col = "grey40")
    lines(x, Y.11y030, lty = 3, col = "grey40")
    lines(x, Y.t305, lty = 3, col = "grey40")
    lines(x, Y.t330, lty = 3, col = "grey40")
    lines(x, Y.t4005, lty = 3, col = "grey40")
    lines(x, Y.t4100, lty = 3, col = "grey40")
    x11()
    plot(x, y, type = "l", lwd = 3, ylim = c(1363.5, 1367.5), 
        xlim = c(1900, 2000), xlab = "Year", ylab = "TSI reconstruction", 
        main = "Total Solar Irradiance")
    grid()
    legend(1900, 1367.5, c("Lean (2004)       ", "SW06a        ", 
        "All years        ", "1st year    "), lwd = c(3, 2, 2, 
        1), lty = c(1, 1, 2, 1), cex = 0.8, col = c("black", 
        "grey30", "grey70", "grey50"), bg = "grey95")
    lines(year.comb, TSI.comb, pch = 19, col = "grey30", lwd = 2)
    n.akrim <- length(AKRIM$year)
    for (i in 1:n.akrim) {
        lines(rep(AKRIM$year[i], 2), AKRIM$TSI[i] + c(-1, 1) * 
            AKRIM$TSI.sd[i], col = "grey80", lty = 1)
        lines(rep(AKRIM$year[i], 2) + c(-0.3, 0.3), rep(AKRIM$TSI[i] + 
            AKRIM$TSI.sd[i], 2), col = "grey80", lty = 1)
        lines(rep(AKRIM$year[i], 2) + c(-0.3, 0.3), rep(AKRIM$TSI[i] - 
            AKRIM$TSI.sd[i], 2), col = "grey80", lty = 1)
        points(AKRIM$year[i], AKRIM$TSI[i], pch = 19, col = "grey30")
    }
    lines(year.comb, TSI.comb.all, pch = 19, col = "grey70", 
        lty = 2, lwd = 2)
    lines(year.comb, TSI.comb.1, pch = 21, col = "grey30", lwd = 1)
    x11()
    plot(1/(sp.T$freq * 12), sp.T$spec, log = "y", xlim = c(0, 
        60), ylim = c(1e-15, 100), type = "l", main = "Spectral analysis", 
        lwd = 2, xlab = "Time scale", ylab = "Power")
    grid()
    lines(1/(sp.S$freq * 12), sp.S$spec, col = "grey30", lty = 2)
    lines(1/(sp.SW$freq * 12), sp.SW$spec, col = "grey60", lty = 2)
    lines(1/(sp.ctl$freq), sp.ctl$spec, col = "grey80", lty = 2)
    lines(1/(sp.co2$freq), sp.co2$spec, col = "black", lty = 3)
    legend(30, 1e-11, c("<T>         ", "Lean (2004)        ", 
        "SW06a        ", "CTL      ", "CO2       "), col = c("black", 
        "grey30", "grey60", "grey80", "black"), lty = c(1, 2, 
        2, 2, 3), lwd = c(2, 1, 1, 1, 1), cex = 0.8, bg = "grey95")
    x11()
    plot(1/(sp.D4$freq * 12), sp.D4$spec, log = "y", xlim = c(0, 
        60), type = "l", main = "Spectral analysis: test band-passed components", 
        lwd = 2, xlab = "Time scale", ylab = "Power TSI")
    grid()
    lines(1/(sp.D3$freq * 12), sp.D3$spec, col = "grey30", lwd = 2)
    lines(1/(sp.S4$freq * 12), sp.S4$spec, col = "grey60", lty = 2)
    legend(30, 1e-11, c("D4 (~22yr) ", "D3 (~11yr)", "S4 (secular)"), 
        col = c("black", "grey30", "grey60"), lty = c(1, 1, 2), 
        lwd = c(2, 2, 1), bg = "grey95")
    x11()
    hist(c(coefSW2006b), main = "Z-coefficients from GISS CTL after SW06b", 
        lwd = 3)
    grid()
    if (do.MonteCarlo) {
        monte.carlo(nt = length(T2m.mon), sd.sun = sd(TSI.comb, 
            na.rm = TRUE), sd.temp = sd(T2m.mon, na.rm = TRUE), 
            m.sun = mean(TSI.comb, na.rm = TRUE), m.temp = mean(T2m.mon, 
                na.rm = TRUE))
        grid()
        title(sub = "Monte-Carlo simulations (random walk)")
        monte.carlo(nt = length(T2m.mon), sd.sun = sd(TSI.comb, 
            na.rm = TRUE), sd.temp = sd(T2m.mon, na.rm = TRUE), 
            m.sun = mean(TSI.comb, na.rm = TRUE), m.temp = mean(T2m.mon, 
                na.rm = TRUE), randomwalk = FALSE)
        grid()
        title(sub = "Monte-Carlo simulations (white noise)")
        monte.carlo(nt = length(T2m.mon), sd.sun = sd(TSI.comb, 
            na.rm = TRUE), sd.temp = sd(T2m.mon, na.rm = TRUE), 
            m.sun = mean(TSI.comb, na.rm = TRUE), m.temp = mean(T2m.mon, 
                na.rm = TRUE), MA.filt = TRUE, randomwalk = FALSE)
        grid()
        title(sub = "Monte-Carlo simulations (MA-filter)")
    }
    x11()
    plot(yr.co2, co2, type = "l", lwd = 3, 
        main = "Mauna Loa annual mean CO2", ylab = "ppmv", xlab = "time", 
        sub = "http://cdiac.ornl.gov/ftp/trends/co2/maunaloa.co2")
    grid()
    x11()
    plot(Obs.tsi$residual, main = "residuals", lwd = 3, col = "green", 
        type = "l")
    grid()
    lines(All.tsi$residual, col = "blue")
    lines(Res.tsi$residual, col = "red")
    lines(Sol.tsi$residual, col = "black")
    grid()
    matrix <- cbind(forcings$Solar, forcings$W.M_GHGs, 
        forcings$O3, forcings$StratH2O, forcings$LandUse, 
        forcings$SnowAlb, forcings$StratAer, forcings$BC, 
        forcings$ReflAer, forcings$AIE)
    udv <- svd(matrix)
    print("SVD:")
    print(round(100 * udv$d^2/sum(udv$d^2)))
    print(sqrt((0.01/0.06)^2 + (0.1/0.35)^2))
}
