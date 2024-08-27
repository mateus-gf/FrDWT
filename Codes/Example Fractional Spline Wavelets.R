#===========================================#
#  Example Fractional splines wavelets in R #
#===========================================#

rm(list = ls())
cat("\014")

#Sources

SDir <- "D:\\FractionalWavelets\\"
#Functions 

source(paste0(SDir, "functions.R"))

alpha <- 0
J     <- 10
N     <- 2^J 
type  <- "+ortho"
x     <- simulate.doppler(n=2^J)
x     <- simulate.chirps(n = 2^J)

plot.ts(x, main="Example Chirps Signal")
grid()
dev.off()
#Get filters 

#filt <- get_FFT_filters(N = N, alpha = alpha, type = type)

w    <- Frac.dwt(signal = x, J = J,type = type, alpha = alpha, out.mat=TRUE)

plot.dwt(w,type = type,alpha=alpha)

rec  <- Frac.idwt(w = w$w, J = J, type = type, alpha = alpha)

plot.ts(x, main = "Original (B) x Reconstructed (R)")
lines(rec, col='red', lty='dashed')
grid()

#Sum os squared errors

sum((x-rec)^2)

#Relative error as calculated in Bruce and Gao
relative.error(signal = x, approximation = rec)
SMAPE(actual = x, predicted = rec)
#Check in decibells

dev.off()

#=====================================#
#  Wavelet Scalogram and Periodogram  #
#=====================================#

w  <- Frac.dwt(signal = x, J = J,type = type, alpha = alpha, out.mat=TRUE)

plot.dwt(w.mat = w,type = type, alpha = alpha)

#Periodograma de wavelets 
Id <- frac.dwt.periodogram(w = w,J=J)
S  <- frac.dwt.scalogram(w = w, I=Id, J = J)

plot.ts(Id$I.w)

#============================#
#  Multiresolution analysis  #
#============================#

mra.list <- frac.MRA(w = w, J = J,type = type,alpha = alpha)

for(i in 1:(J+1)){
  plot.ts(mra.list[[i]], ylab='Reconstruction', main = "MRA of signal")
  Sys.sleep(1)
}

lines(x, col='red', lty='dashed')

mrd.list <- frac.MRD(w = w, J = J,type = type,alpha = alpha)

for(i in 1:(J+1)){
  plot.ts(mrd.list[[i]], ylab='Reconstruction', main = "MRD of signal")
  Sys.sleep(1)
}

dev.off()

#==============================#
# Example of Wavelet Shrinkage #
#==============================#

set.seed(123)

type='*ortho'
alpha=1.5

cont <- x+rnorm(2^J)/20

plot.ts(cont, main='Contaminated series', ylab='', col='darkblue')
grid()

wd   <- Frac.dwt(signal = cont, J = J,type = type, alpha = alpha, out.mat=TRUE)
recs <- wavelet.shrinkage(w = wd$w,J = J, type = type, alpha = alpha, reconstruct = T)

plot.ts(recs$shrinked)

plot.ts(x, ylim=c(-1.5,1.5), main = 'Wavelet Shrinkage procedure')
lines(recs$shrinked, col='red',lty='dashed')
grid()
dev.off()

#Assess the quality of wavelet shirnkage
Metrics::rmse(actual = x, predicted = recs$shrinked)
relative.error(signal = x, approximation = recs$shrinked)

#===============================#
# Bootstrap confidence interval #
#===============================#

#Example of bootstrapping with FrWT
#This is only ilustrative. Better performance can be achieved using a Fortran or C code to perform bootstrap

bootstrap <- bootstrap.ci(signal = cont, N.rep = 100, type = type, J = J,alpha=alpha)

plot.ts((recs$shrinked), main="Bootstrap CI using FrWT Shrinkage", col='blue', lwd=1, ylab='Estimated signal',
        ylim = c(-2,2))
xd <- c(0:(2^J-1), rev(0:(2^J-1)))
yd <- c(bootstrap$lower.ci, rev(bootstrap$upper.ci))
polygon(xd,yd, col='red', lty='dashed', lwd=0.5, border='red')
legend(0, -1.1, legend=c('Estimated', 'Bootstrap CI'), col=c('blue', 'red'), border = c("blue", 'red'), fill =c("blue", 'red'))
grid()

dev.off()

#=====================================#
# MATLAB example used by Blu & Unser  #
#=====================================#

alpha <- 0.25
type  <- '*ortho'
#Generate a 1/f-like noise

noise <- simulate.fractal.noise(J = J, alpha = (alpha+1)/2, type = type,seed=123,beta = 1+alpha)

par(mfrow=c(1,2))
plot.ts(noise, main='Synthetized 1/f-noise', col='blue')
grid()

M  <- 2^J
y  <- 1/M * abs(fft(x))^2
y  <- c(y[(1+M/2):M], y[1:(M/2)])

nu <- abs(seq(from = -M/2, to = M/2-1, by=1))/M

#Plot in Decibells
plot(x=nu, 10*log(y, base = 10), 
     xlab = 'Normalized Frequency',
     ylab = 'Fourier amplitude (dB)', 
     main = 'DSP of the synthetized noise', type='l',
     log  = 'x') 
grid()

dev.off()

#GPH estimate for d

fracdiff::fdGPH((noise))$d-1

#noise     <- simulate.chirps(n = 2^J)

w    <- Frac.dwt(signal = noise, J = J,type = type, alpha = alpha, out.mat=TRUE)

plot.ts(noise)

plot.dwt(w,type = type,alpha=alpha)


mra.list <- frac.MRA(w = w, J = J,type = type,alpha = alpha)

for(i in 1:(J+1)){
  plot.ts(mra.list[[i]], ylab='Reconstruction', main = "MRA of signal")
  Sys.sleep(1)
}

lines(noise, col='red', lty='dashed')

dev.off()
#==========================#
#  Analysis of pink noise  #
#==========================#

clean <- white.signal(signal = noise, J=J, type = type, alpha = alpha)

#Comparison of the series 

par(mfrow=c(1,2))
plot.ts(noise, col='magenta', main='Simulated 1/f noise', ylab='')
grid()
plot.ts(clean, col='blue', main='Whitened 1/f noise with FrWT', ylab='')
grid()
dev.off()

#Spectral density

M  <- 2^J
y  <- 1/M * abs(fft(noise))^2
y  <- c(y[(1+M/2):M], y[1:(M/2)])

yw <- 1/M * abs(fft(clean))^2
yw <- c(yw[(1+M/2):M], yw[1:(M/2)])
nu <- abs(seq(from = -M/2, to = M/2-1, by=1))/M

par(mfrow=c(1,2))

#Plot in Decibells
plot(x=nu, 10*log(y, base = 10), 
     xlab = 'Normalized Frequency',
     ylab = 'Fourier amplitude (dB)', 
     main = 'DSP of the 1/f synthetized noise', type='l',
     log  = 'x', col='magenta') 
grid()

plot(x=nu, 10*log(yw, base = 10), 
     xlab = 'Normalized Frequency',
     ylab = 'Fourier amplitude (dB)', 
     main = 'DSP of the whitened series', type='l',
     log  = 'x', col='blue', ylim=c(-33, 60)) 
grid()
#Autocorrelations 

par(mfrow=c(1,2))
acf(clean,60, main='Whitened ACF')
grid()
pacf(clean,60, main='Whitened PACF')
grid()

dev.off()

#Generating a white noise

test <- sim.pink.noise(J = J, alpha = alpha, type = type, seed = 123)

plot.ts(test)

acf(test)
pacf(test)

#Jarque-Bera test for normality 

tseries::jarque.bera.test(x = test)

#=====================================#
#    Generate a FWN with the FrWT     #
#=====================================#

fwn.sim <- sim.FWN(J = J,alpha = alpha,type = type, seed = 123)

par(mfrow=c(1,2))
plot.ts(fwn.sim, main='Synthetized FWN', col='blue')
grid()

M  <- 2^J
y  <- 1/M * abs(fft(fwn.sim))^2
y  <- c(y[(1+M/2):M], y[1:(M/2)])

#Plot in Decibells
plot(x=nu, 10*log(y, base = 10), 
     xlab = 'Normalized Frequency',
     ylab = 'Fourier amplitude (dB)', 
     main = 'DSP of the 1/f synthetized noise', type='l',
     log  = 'x', col='magenta') 
grid()

#Autocorrelation function

par(mfrow=c(1,2))
acf(fwn.sim,60, main='FWN ACF')
grid()
pacf(fwn.sim,60, main='FWN PACF')
grid()

fracdiff::fracdiff(fwn.sim)$d
alpha

