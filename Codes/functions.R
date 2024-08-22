library(moments)

#==========================#
#  Wavelet Energy Entropy  #
#==========================#

wavelet.energy.entropy <- function(w, J){
  #Computes wavelet energy entropy
  # Calculates the energy per level as a proportion of total energy.
  # This is p[j] for each level [j]
  # Then WEE is - sum(p[j]*log(p[j], base=2))
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  if((round(J)==J)==F){stop("J must be integer!")}
  if(is.list(w)){stop("Input must be a vector!")}
  M = length(w)
  if(is.pow.two(M)==FALSE){stop("Input must be a power of two.")}
  
  E <- sum(w^2)
  
  tab <- generate.dyadic(J = J)
  
  p <- rep(x = NA, nrow(tab)-1)
  
  for(i in 1:(nrow(tab)-1)){
    p[i] <- sum(w[(tab[i,"Start"]:tab[i,"End"])]^2)/E
  }
  
  
  WEE <- -sum(p*log(p, base = 2)) 
  return(WEE)
}

#==========================#
#  Relative error funciton #
#==========================#

relative.error <- function(signal, approximation){
  #Reference: 
  #Bruce, Gao 'Applied Wavelet Analysis with S-PLUS', 1996
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  return(sqrt(sum((signal-approximation)^2))/sqrt(sum(signal^2)))
}

SMAPE <- function(actual, predicted){
  #Function to calculate SMAPE
  #Following Martinez et al (2017) but without multiplying by '100'
  #We do so to match other metrics used
  #Author: Mateus Gonzalez de Freitas Pinto, 2022
  if(length(actual)!=length(predicted)){stop("Length of predcited and actual series must be the same")}
  N <- length(actual) 
  return((1/N)*sum(abs(actual-predicted)/((abs(actual)+abs(predicted))/2)))
}


#============================#
# Calculate LP norm function #
#============================#

lp_norm <- function(x,p){
  stopifnot(is.numeric(x) || is.complex(x), is.numeric(p))
  if(p<=0){stop("p must be positive!")}
  if(typeof(x)!='double' & typeof(x) != 'vector'){
    stop("Input must be a vector or double")
  }
  
  if(p==1){
    norm <- sum(abs(x))
  }else if(p==2){
    norm <- sqrt(sum(x^2))
  }else if(p==Inf){
    norm <- max(abs(x))
  }else{
    norm <- sum(abs(x)^p)^(1/p)
  }
  
  return(norm)
}

#=====================#
# Simulated Signals   #
#=====================#

simulate.doppler <- function(n=1024){
  #Computes wavelet energy entropy
  # Simulate doppler signals following S+Wavelets
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
    t = seq(0,1,length.out = n)
    sqrt(t * (1 - t)) * sin((2 * pi * 1.05)/(t + 0.05))
}

simulate.chirps <- function(n=1024){
  #Computes wavelet energy entropy
  # Simulate chirps signal, following S+Wavelets
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  x <- 1e-05 + seq(from = -1, to = 1, length = n + 1)[1:n]
  y <- sin(pi/x)
  return(y)
}

#====================================#
# Check if signal is a power of two  #
#====================================#

is.pow.two <- function(x){
  # Check if it is power of two
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  if(x%%2==0){return(TRUE)}else{return(FALSE)}
}

#==========================#
#  Give names to wavelets  #
#==========================#

.give.name <- function(type){
  # Internal use for beautiful plots
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  nm <- gsub("\\-", "Anticausal ", gsub("\\+", "Causal ",gsub("\\*", "Symmetric ",type)))
  nm <- gsub("ortho", "Orthonormal Frac. Spline Wavelets",nm)
  nm <- gsub("dual", "Frac. Dual Spline Wavelets",nm)
  nm <- gsub("B-spline", " Frac. B-spline",nm)
  return(nm)
}

#=========================#
# MATLAB-compatible IFFT  #
#=========================#

ifft <- function(x){
  #A MATLAB-compatible Inverse FFT
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  return(fft(x, inverse = TRUE)/length(x))
}

#================#
#  SINC function #
#================#

sinc <- function(x){
  # Sinc(x)= sin(pi*x)/(pi*x)
  # Author: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  if(is.numeric(x)==F){stop("Numeric-only inputs allowed")}
  if(x!=0){
    y=sin(pi*x)/(pi*x)
  }else if(x==0){
    y=1
  }
  
  return(y)
}


#====================================#
# Fractional Splines autocorrelation # 
#====================================#

fractsplineautocorr <- function(alpha, nu, M=100){
  #Frequency domain computation of fractional spline autocorrelation
  #n=100 is default number of terms of the summation for computing the autocorrelation. Can be changed.
  if(alpha<=-0.5){stop("The autocorrelation of the fractional splines exists only for degrees strictly larger than -0.5!")}
  
  S    <- rep(0, length(nu))
  err  <- c()
  err0 <- c()
  
  for(m in -M:M){
    S = S+abs(sapply(X = as.vector(nu+m), sinc))^(2*alpha+2) 
  }
  
  U <- 2/(2*alpha+1)/M^(2*alpha+1);
  U <- U-1/M^(2*alpha+2);
  U <- U+(alpha+1)*(1/3+2*nu*nu)/M^(2*alpha+3);
  U <- U-(alpha+1)*(2*alpha+3)*nu*nu/M^(2*alpha+4);
  U <- U*abs(sin(pi*nu)/pi)^(2*alpha+2);
  A <- S + U 
  
}


#=====================#
# Get FFT FWT filters #
#=====================#

get_FFT_filters <- function(N, alpha, type){
  #+ causal
  #- anticausal
  #* symmetric
  
  if(alpha<=-0.5){stop("The autocorrelation of the fractional splines exists only for degrees strictly larger than -0.5!")}
  if(is.pow.two(N)==FALSE){stop("Signal must be of length 2^J, for a natural J!")}
  
  nu <- seq(0, 1-1/N, by=1/N)
  
  A  <-fractsplineautocorr(alpha=alpha,nu=nu)
  A2 <- c(A,A)
  A2 <- A2[seq(1,length(A2), by=2)]
  
  if(length(grep("ortho", type)>0)){
    if(substr(type,1,1)=="+"){
      lowa=sqrt(2)*((1+exp(-2*1i*pi*nu))/2)^(alpha+1)*sqrt(A/A2)
    }else if(substr(type,1,1)=="-"){
      lowa=sqrt(2)*((1+exp(2*1i*pi*nu))/2)^(alpha+1)*sqrt(A/A2)
    }else if(substr(type,1,1)=="*"){
      lowa=sqrt(2)*abs((1+exp(-2*as.complex(1i)*pi*nu))/2)^(alpha+1)*sqrt(A/A2)
    }else{
      stop("Unknown filter prefix!")
      }
    
    higha=exp(-2*1i*pi*nu)*lowa
    higha=Conj(higha[c((N/2 + seq(1, N/2)), seq(1, N/2))])
    
    lows=lowa
    highs=higha
    filters <- list(FFTanalysisfilters  = list(lowa=lowa, higha=higha),
                    FFTsynthesisfilters = list(lows=lows, highs=highs))
  }else{
    #semi-orthonormal and spline filters
    if(substr(type,1,1)=='*'){
      lowa=sqrt(2)*abs((1+exp(-2*1i*pi*nu))/2)^(alpha+1)
    }else if(substr(type,1,1)=='+'){
      lowa=sqrt(2)*((1+exp(-2*1i*pi*nu))/2)^(alpha+1)
    }else if(substr(type,1,1)=='-'){
      lowa=sqrt(2)*((1+exp(2*1i*pi*nu))/2)^(alpha+1)
    }else{
      stop("Unknown filter prefix!")
    }
    
    higha=exp(-2*1i*pi*nu)*lowa*A
    higha=Conj(higha[c((N/2 + seq(1, N/2)), seq(1, N/2))])
    lows=lowa*A/A2
    highs=higha/(A2*A[c((N/2 + seq(1, N/2)), seq(1, N/2))])
    
    if(length(grep("dual",type))>0){#Dual spline wavelets
      filters <- list(FFTanalysisfilters  = list(lowa=lowa, higha=higha),
                      FFTsynthesisfilters = list(lows=lows, highs=highs))
    }else if(length(grep("B-spline", type))>0){#B-spline wavelets
      filters <- list(FFTsynthesisfilters  = list(lows=lowa, highs=higha),
                      FFTanalysisfilters = list(lowa=lows, higha=highs))
    }else{
      stop("Unknown wavelet")
    }
  }

  return(filters)
}


#=========================================#
#  Discrete Fractional Wavelet Transform  # 
#        and also its inverse DWT         #
#=========================================#

#In BLU & UNSER this is referred to as wavelet anaylsis and synthesis

Frac.dwt <- function(signal, J=log(length(signal), 2), type, alpha, out.mat=TRUE){
  #Fast Fourier based implementation of wavelet transform, based on Blu and Unser's Matlab code
  #Computes the wavelet transform of a signal x using the Fourier method. 
  #It uses periodic boudnary condition 
  #The legnth of the signal is 2^N
  #Frequency response filters are results of the get_FFT_filters
  
  #signal = input signal, size 2^N = length of signal 
  #J: depth of decomposition: J wavelet bands (wav1 to wavJ)+lowpass lowJ
  #type = families of wavelet used. In this case, must be same names used in get_FFT_filters
  
  #Output: 
  # w = wav1, wav2, ... wavJ, lowJ
  # If list = TRUE, then atribute coefficients to a list. Else, to a vector
  # The latter is better for the algorithm of IDFWT
  
  #Based on original MATLAB code of Thierry Blu, October 1999
  # 	Biomedical Imaging Group, EPFL, Lausanne, Switzerland.
  # 	This software is downloadable at http://bigwww.epfl.ch/
  #References:
  # 	[1] M. Unser and T. Blu, "Fractional splines and wavelets," 
  # 	SIAM Review, Vol. 42, No. 1, pp. 43--67, January 2000.
  # 	[2] M. Unser and T. Blu, "Construction of fractional spline wavelet bases," 
  # 	Proc. SPIE, Wavelet Applications in Signal and Image Processing VII,
  #     Denver, CO, USA, 19-23 July, 1999, vol. 3813, pp. 422-431. 
  # 	[3] T. Blu and M. Unser, "The fractional spline wavelet transform: definition and 
  #	implementation," Proc. IEEE International Conference on Acoustics, Speech, and 
  #	Signal Processing (ICASSP'2000), Istanbul, Turkey, 5-9 June 2000, vol. I, pp. 512-515 .
  
  #Author of R version: Mateus Gonzalez de Freitas Pinto
  # Used for the following articles: 
  # [1] PINTO, M. G. DE F. and CHIANN, C. "Long-memory parameter estimation based on fractional spline wavelets"
  # Digital Signal processing Vol. 133, March 2023, pp. 103836-1-12
  # [2] PINTO, M. G. DE F. and CHIANN, C. "A maximum-likelihood-based approach to estimate the long memory parameter using fractional spline wavelets"
  # Signal Processing Vol. 222, 2024, pp. 109518
  #
  #Sao Paulo, Brazil, 2022
  
  if(is.data.frame(signal)){stop("Signal must be either a vector or numeric format!")}
  if(is.pow.two(length(signal))==FALSE){stop("Length of signal must be a power of two!")}
  if(any(is.character(signal))){stop('Signal must numeric')}
  
  signal <- as.numeric(as.character(signal))
  
  if(any(is.na(signal))){stop('Missing values in input signal')}
  
  
  M = length(signal)
  
  filt <- get_FFT_filters(N = M, alpha = alpha, type = type)
  
  if((length(filt$FFTanalysisfilters$lowa)==M & length(filt$FFTanalysisfilters$higha)==M &
      length(filt$FFTsynthesisfilters$lows==M & length(filt$FFTsynthesisfilters$highs)==M))==FALSE){
    stop("Length of filters and length of signal must match!")
  }
  
  X <- fft(z = signal)
  
  G <- filt$FFTanalysisfilters$lowa
  H <- filt$FFTanalysisfilters$higha
  w <- c()

  w.mat <- matrix(nrow = J+1, ncol = M)
  
  for(j in 1:J){
     Y <- G*X
     Z <- H*X
     Y <- 0.5*(Y[seq(1, M/2,by=1)]+Y[M/2+(seq(1, M/2,by=1))])
     Z <- 0.5*(Z[seq(1, M/2,by=1)]+Z[M/2+(seq(1, M/2,by=1))])
     z <- ifft(x = Z)
     w <- c(w,z)
     w.mat[j,1:length(z)] <- z
     M <- M/2
     X <- Y
     G <- G[seq(1, length(G), by=2)]
     H <- H[seq(1, length(H), by=2)]
     
   }
    
  w <- c(w,ifft(X))
  w.mat[J+1,1] <- ifft(X)
  w <- Re(w)
  w.mat <- Re(w.mat)
  
  out.list <- list(w=w, w.mat=w.mat)
  
  if(out.mat==TRUE){
    return(out.list)
  }else{
    return(w)
  }
}



Frac.idwt <- function(w, J=log(length(w), 2), type, alpha){
  if(is.list(w)){stop("Input must be a vector!")}
  
  M = length(w)
  if(is.pow.two(M)==FALSE){stop("Input must be a power of two.")}
  
  filt <- get_FFT_filters(N = M, alpha=alpha, type=type)
  
  if((length(filt$FFTanalysisfilters$lowa)==M & length(filt$FFTanalysisfilters$higha)==M &
      length(filt$FFTsynthesisfilters$lows==M & length(filt$FFTsynthesisfilters$highs)==M))==FALSE){
    stop("Length of filters and length of signal must match!")
  }
  
  G <- Conj(filt$FFTsynthesisfilters$lows)
  H <- Conj(filt$FFTsynthesisfilters$highs)
  
  M <- M/2^J

  y <- w[length(w)+(-M+1):0]
  w <- w[1:(length(w)-M)]
  
  Y <- fft(y)
  
  for(j in seq(J, 1, by=-1)){
    z <-  w[length(w)+((-M+1):0)]
    w <-  w[1:(length(w)-M)]
    Z <- fft(z)
    M <- 2*M
    
    H1 <- H[seq(1,length(H), by=(2^(j-1)))]
    G1 <- G[seq(1,length(G), by=(2^(j-1)))]
    
    Y0 <- G1[seq(1,M/2)]*Y + H1[seq(1,M/2)]*Z
    Y1 <- G1[M/2+seq(1, M/2)]*Y+H1[M/2+seq(1, M/2)]*Z
    
    Y  <- c(Y0,Y1)
    
  }
  
  x <- Re(ifft(Y))
  return(x)
}

#===========================#
#   Plot wavelet crystals   #
#===========================#

plot.dwt <- function(w.mat, main="Fractional Splines DWT Coefficients", ylab="Resolution level", xlab="Translate",first.level=0,type, alpha){
  #Plot the coefficients of the discrete fractional wavelet transform, given the w
  #w.mat is either the matrix form of the fractional DWT or the list output of function Frac.dwt
  #type is the type of wavelet one uses to perform the FDWT
  #alpha is the order of the FDWT
  #We plot following Daubechies's notation in "Ten lectures on wavelets":
  # the smoother coefficients are in the highest level, whereas the details are in the lower levels
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  
  if(is.list(w.mat)){
    if(is.null(w.mat$w.mat)){
      stop("Format of the input is wrong")
    }else{
      w.mat <- w.mat$w.mat
    }
  }else if(is.vector(w.mat)){stop("Input must be a matrix")}
  
  
  if(is.pow.two(ncol(w.mat))==F){stop("Input must be a power of two!")}
  
  #NotPlotVal  <- 0.005
  NotPlotVal  <- -Inf
  J           <- log(ncol(w.mat),2)
  levels      <- J 
  nlevels     <- levels - first.level
  n           <- 2^(levels - 1)
  
  nm  <- .give.name(type)
  sub <- paste0(nm, " w/ ", expression(alpha),"=", alpha)
  
  plot(c(0, 0, n, n), c(0, nlevels + 1, nlevels + 1, 0), type = "n", 
       xlab = xlab, ylab = ylab, main = main, yaxt = "n", 
       xaxt = "n", sub = sub)

  yll <- (levels - 1):first.level+1
  axis(2, at = 1:(nlevels), labels = rev(yll))
  
  axx <- c(0, 2^(levels - 3), 2^(levels - 2), 2^(levels - 2) + 2^(levels - 3), 2^(levels - 1))
  
  myxx   <- 1:n
  height <- 1
  axr    <- NULL
  
  my <- 0
  for (i in (((levels - 1):first.level)+1)) {
    y  <- w.mat[i,]
    y  <- y[which(!(is.na(y)))]
    my <- max(c(my, abs(y)))
  }
  
  shift <- 1
for (i in (((levels - 1):first.level)+1) ) {
    y  <- w.mat[nrow(w.mat)-i,]
    y  <- y[which(!(is.na(y)))]
    n  <- 2^i
    xplot <- myxx
    ly <- length(y)
    if (my == 0) {
      y <- rep(0, length(y))
    }else{y <- (0.5 * y)/my}
    
    axr <- c(axr, my)
    if (max(abs(y)) > NotPlotVal){ 
      segments(xplot, height, xplot, height + y)}
    if (i != first.level) {
      x1 <- myxx[seq(1, n - 1, 2)]
      x2 <- myxx[seq(2, n, 2)]
      myxx <- (x1 + x2)/2
      
      height <- height + 1
    }
  }
  
}


#==============================#
#   Multiresolution Analysis   #
#==============================#

generate.dyadic <- function(J){
  #Auxiliary function: generate start and end values to separate in levels the DWT
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  if((round(J)==J)==F){stop("J must be integer!")}
  powers <- rev(sapply(seq(0,J-1,by=1), function(x) 2^x))

  coef <- c(paste0('d',1:(J)),paste0('c',J))
  posi <- 2^J-powers
  posi <- c(posi, 2^J)
  tab <- data.frame(cbind(coef, posi), stringsAsFactors = F)
  tab$posi <- as.numeric(as.character(tab$posi))
  colnames(tab) <- c("Crystal", "End")
  tab$Start     <- 1
  for(i in 2:(nrow(tab))){
    tab[i,"Start"] <- tab[i-1, "End"]+1
  }
  
  
  
  tab <- tab[,c("Crystal", "Start", "End")]
  return(tab)
}

frac.MRA <- function(w, J, type, alpha){
  #Multiresulution analysis
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  if(is.pow.two(length(w$w))==F){stop("Input must be a power of two!")}
  if(missing(type) | missing(alpha)){stop("Alpha and type must be specified")}
  tab <- generate.dyadic(J)
  
  w.ywd <- w$w
  
  MRA <- list()
  
  w0 <- w.ywd
  w0[1:(2^J-1)] <- 0
  
  S0 <- Frac.idwt(w = w0, J = J, type = type, alpha = alpha)
  
  MRA[[1]] <- S0
  for(j in rev(1:J)){
    w.aux <- w.ywd
    w.aux[1:tab[j,]$Start-1] <- 0
    assign(x = paste0("D",j-1),Frac.idwt(w = w.aux,J = J, type = type, alpha = alpha))
    MRA[[J-j+2]] <- get(paste0("D",j-1))
  }
  

  names(MRA) <- c("S0",paste0("D", (1:J)-1))
  return(MRA)
}


#===================================#
#   Multiresolution Decomposition   #
#===================================#

frac.MRD <- function(w, J, type, alpha){
  #Returns multiresolution decomposition, with each of the components of the signal
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  
  if(is.list(w)){
    w <- w$w
  }else if(is.matrix(w)){stop("w must be a vector")}
  
  if(is.pow.two(length(w))==F){stop("Input must be a power of two!")}
  if(missing(type) | missing(alpha)){stop("Alpha and type must be specified")}
  tab <- generate.dyadic(J)
  
  w.ywd <- w
  
  MRD <- list()
  
  w0 <- w.ywd
  w0[1:(2^J-1)] <- 0
  
  S0 <- Frac.idwt(w = w0, J = J, type = type, alpha = alpha)
  
  MRD[[1]] <- S0
  for(j in rev(1:J)){
    w.aux <- w.ywd
    seq   <- 1:length(w.ywd)
    seq   <- seq[which(!(seq %in% (tab[j,'Start']:tab[j,'End']) ))]
    w.aux[seq] <- 0
    assign(x = paste0("D",j-1),Frac.idwt(w = w.aux,J = J, type = type, alpha = alpha))
    MRD[[J-j+2]] <- get(paste0("D",j-1))
  }
  
  
  names(MRD) <- c("S0",paste0("D", (1:J)-1))
  return(MRD)
}

#==========================#
#   Wavelet thresholding   # 
#==========================#

wavelet.shrinkage <- function(w, J, type, alpha, policy='hard', thresh='MAD', reconstruct=TRUE, shrink.rule='universal',
                              scale.rules='all',p, constant=1){
  if(scale.rules=='each'){
    #For rules scales each, call function as each
    
    obj <- each.threshold(w = w,J = J, type = type, alpha = alpha, thresh=thresh,policy = policy, reconstruct = reconstruct, constant = constant)
    
    return(obj)
    
  }else if(scale.rules=='all'){
    if(is.list(w)){
      if(is.null(w$w)){stop("Input must be a vector!")}
      else{
        w=w$w
      }
    }
    
    if(thresh=='LP'){if(missing(p)){stop("LP norm chosen with unspecified 'p' ! ")}}
    
    if(is.pow.two(length(w))==F){stop("Input must be a power of two!")}
    if(missing(type) | missing(alpha)){
      if(reconstruct == TRUE){
        stop("Alpha and type must be specified")
      }
    }
    
    if(shrink.rule=='universal'){
      if(thresh=='MAD' | thresh == 'mad'){
        lambda <- constant*stats::mad(w, center = 0)*sqrt(2*log(length(w))) 
      }else if(thresh=='SD' | thresh == 'sd'){
        lambda <- constant*stats::sd(w)*sqrt(2*log(length(w))) 
      }else if(thresh=='LP' | thresh=='lp'){
        lambda <- constant*lp_norm(w,p=p)*sqrt(2*log(length(w))) 
      }else if(thresh=='MEDIAN' | thresh=='median'){
        lambda <- constant*median(w)*sqrt(2*log(length(w)))
      }else if(thresh=='UNITY' | thresh=='unity'){
        lambda <- constant*sqrt(2*log(length(w)))
      }else{
        stop('Thresh options are: MAD, SD, LP, UNITY, MEDIAN')
      }
      
    }
    
    if(policy=='hard'){
      w[which(abs(w)<=lambda)] <- 0
    }else if(policy=='soft'){
      w[which(abs(w)<=lambda)] <- 0
      w[which(abs(w)>lambda)] <- sign(w[which(abs(w)>lambda)])*(abs(w[which(abs(w)>lambda)]) - lambda)
    }else{
      stop("Policy must be either 'hard' or 'soft'!")
    }
    
    if(reconstruct==FALSE){
      return(w)
    }else if(reconstruct==TRUE){
      shrinked <- Frac.idwt(w = w, J = J,type = type, alpha = alpha)
      ans      <- list(w=w, shrinked=shrinked)
      return(ans)
    }
  }else{
    stop("Scales rules must be either each or all")
  }
}

each.threshold <- function(w, J, type, alpha, policy='hard', thresh='MAD', reconstruct=TRUE, shrink.rule='universal',
                           scale.rules='all',p, constant=1){
  tab <- generate.dyadic(J)
  if(is.list(w)){
    if(is.null(w$w)){stop("Input must be a vector!")}
    else{
      w=w$w
    }
  }
  
  if(thresh=='LP'){if(missing(p)){stop("LP norm chosen with unspecified 'p' ! ")}}
  
  if(is.pow.two(length(w))==F){stop("Input must be a power of two!")}
  if(missing(type) | missing(alpha)){
    if(reconstruct == TRUE){
      stop("Alpha and type must be specified")
    }
  }
  
  wd <- rep(NA, length(w))
  
  for(i in 1:nrow(tab)){
    aux <- w[tab[i,'Start']:tab[i,'End']]
    
    if(thresh=='MAD' | thresh == 'mad'){
      lambda <- constant*stats::mad(aux, center = 0)*sqrt(2*log(length(aux))) 
    }else if(thresh=='SD' | thresh == 'sd'){
      lambda <- constant*stats::sd(aux)*sqrt(2*log(length(aux))) 
    }else if(thresh=='LP' | thresh=='lp'){
      lambda <- constant*lp_norm(aux,p=p)*sqrt(2*log(length(aux))) 
    }else if(thresh=='MEDIAN' | thresh=='median'){
      lambda <- constant*median(aux)*sqrt(2*log(length(aux)))
    }else if(thresh=='UNITY' | thresh=='unity'){
      lambda <- constant*sqrt(2*log(length(aux)))
    }else{
      stop('Thresh options are: MAD, SD, LP')
    }
    
   
    if(policy=='hard'){
      aux[which(abs(aux)<=lambda)] <- 0
    }else if(policy=='soft'){
      aux[which(abs(aux)<=lambda)] <- 0
      aux[which(abs(aux)>lambda)] <- sign(aux[which(abs(aux)>lambda)])*(abs(aux[which(abs(aux)>lambda)]) - lambda)
    }else{
      stop("Policy must be either 'hard' or 'soft'!")
    }
    
    wd[tab[i,'Start']:tab[i,'End']] <- aux
    
  }
  
  if(reconstruct==FALSE){
    return(wd)
  }else if(reconstruct==TRUE){
    shrinked <- Frac.idwt(w = wd, J = J,type = type, alpha = alpha)
    ans      <- list(w=wd, shrinked=shrinked)
    return(ans)
  }else{
    stop("reconstruct must be either TRUE or FALSE")
  }
}




wavelet.shrinkage.reverse <- function(w, J, type, alpha, levels, n.levels=J, policy='hard', thresh='MAD', reconstruct=TRUE, shrink.rule='universal',
                                      scale.rules='all',p,constant=1){
  if(is.list(w)){
    if(is.null(w$w)){stop("Input must be a vector!")}
    else{
      w=w$w
    }
  }
  
  if(thresh=='LP'){if(missing(p)){stop("LP norm chosen with unspecified 'p' ! ")}}
  
  if(is.pow.two(length(w))==F){stop("Input must be a power of two!")}
  if(missing(type) | missing(alpha)){
    if(reconstruct == TRUE){
      stop("Alpha and type must be specified")
    }
  }
  
  if(shrink.rule=='universal'){
    if(thresh=='MAD' | thresh == 'mad'){
      lambda <- constant*stats::mad(w, center = 0)*sqrt(2*log(length(w))) 
    }else if(thresh=='SD' | thresh == 'sd'){
      lambda <- constant*stats::sd(w)*sqrt(2*log(length(w))) 
    }else if(thresh=='LP' | thresh=='lp'){
      lambda <- constant*lp_norm(w,p=p)*sqrt(2*log(length(w))) 
    }else{
      stop('Thresh options are: MAD, SD, LP')
    }
    
  }
  
  if(policy=='hard'){
    w[which(abs(w)>lambda)] <- 0
  }else if(policy=='soft'){
    w[which(abs(w)>lambda)] <- 0
    w[which(abs(w)<=lambda)] <- sign(w[which(abs(w)>lambda)])*(abs(w[which(abs(w)>lambda)]) - lambda)
  }else{
    stop("Policy must be either 'hard' or 'soft'!")
  }
  
  if(reconstruct==FALSE){
    return(w)
  }else if(reconstruct==TRUE){
    shrinked <- Frac.idwt(w = w, J = J,type = type, alpha = alpha)
    ans <- list(w=w, shrinked=shrinked)
    return(ans)
  }
}

#=====================#
#  Plot MATLAB demo   #
#=====================#

simulate.fractal.noise <- function(J,alpha, beta, type,seed){
  # Generation of 1/f-like noise using a fractional
  #	spline wavelet transform.
  # Adapted from T. Blu and M. Unser, 1999
  #
  #	Example of parameters : J=12, alpha=0.33, beta=1.33
  # Author of R version: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  
  w <- c()
  M <- 2^J 
  s <- M/2
  
  if(missing(beta)){
    beta=alpha+1
  }
  if(missing(seed)==F){
    set.seed(seed)
  }
  for(j in 1:J){
    w <- c(w, sqrt(12)*((runif(n = s, min = 0, max = 1))-0.5)*2^(j*beta))
    s <- s/2
  }
  
  w <- c(w, sqrt(12)*((runif(n = 2*s, min = 0, max = 1))-0.5)*2^(J*beta))
  
  x <- Frac.idwt(w = w, J = J, type = type, alpha = alpha)
  return(x)
}

#=========================#
#  Simulate a test-noise  #
#=========================#


sim.pink.noise <- function(J,alpha, beta, type,seed){
  # Generation of 1/f-like noise using a fractional
  #	spline wavelet transform.
  # Adapted from T. Blu and M. Unser, 1999
  #
  #	Example of parameters : J=12, alpha=0.33, beta=1.33
  # Author of R version: Mateus Gonzalez de Freitas Pinto
  # Sao Paulo, Brazil, 2022
  w <- c()
  M <- 2^J 
  s <- M/2
  
  if(missing(beta)){
    beta=alpha+1
  }
  if(missing(seed)==F){
    set.seed(seed)
  }
  for(j in 1:J){
    w <- c(w, sqrt(12)*((runif(n = s, min = 0, max = 1))-0.5)*2^(j*beta))
    s <- s/2
  }
  
  w <- c(w, sqrt(12)*((runif(n = 2*s, min = 0, max = 1))-0.5)*2^(J*beta))
  
  x <- Frac.idwt(w = w, J = J, type = type, alpha = alpha)
  return(x)
}

#========================#
#   Simulate a Frac WN   #
#========================#

sim.FWN <- function(J,alpha, type,seed){
  # Generation of fracional white noise using a fractional
  #	spline wavelet transform.
  # Adapted from T. Blu and M. Unser, 1999
  #
  #	Example of parameters : J=12, alpha=0.33, beta=1.33
  
  w <- c()
  M <- 2^J 
  s <- M/2
  beta=alpha
  
  if(missing(seed)==F){
    set.seed(seed)
  }
  for(j in 1:J){
    w <- c(w, sqrt(12)*((runif(n = s, min = 0, max = 1))-0.5)*2^(j*beta))
    s <- s/2
  }
  
  w <- c(w, sqrt(12)*((runif(n = 2*s, min = 0, max = 1))-0.5)*2^(J*beta))
  
  x <- Frac.idwt(w = w, J = J, type = type, alpha = alpha)
  return(x)
}

#============================#
# Energy wavelet calculation #
#============================#

compute.energy <- function(w, J){
  #As in Bruce, Gao 'Applied Wavelet Analysis with S-PLUS', 1996
  #A little bit compatible to S+wavelets functions
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  if((round(J)==J)==F){stop("J must be integer!")}
  if(is.list(w)){
    if(is.null(w$w)){
      stop("Format of w is wrong")
    }else{
      w <- w$w
    }
  }else if(is.matrix(w)){
    stop("Input must be a vector or the list with vector and matrix")
  }
  
  if(is.pow.two(length(w))==F){stop("w must be a power of two!")}
  tab <- generate.dyadic(J)
  
  energy <- rep(NA, nrow(tab))
  for(i in 1:nrow(tab)){
    energy[i] <- sum(w[tab[i,'Start']:tab[i,'End']]^2)
  }
  
  relative.energy  <- energy/sum(energy)
  vec              <- c(paste0('d',1:J), paste0('c',J)) 
  table            <- cbind(vec,energy, relative.energy)
  table            <- as.data.frame(table)
  colnames(table)  <- c("Level", "Energy", "Rel.Energy")
  table$Energy     <- as.numeric(as.character(table$Energy))
  table$Rel.Energy <- as.numeric(as.character(table$Rel.Energy))
  return(table)
}

energy.plot <- function(energy.table, type, alpha, from0.to1=F){
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  if(missing(alpha) | missing(type)){
    stop("Please supply type and alpha!")
  }
  nm  <- .give.name(type)
  sub <- paste0(nm, " w/ ", expression(alpha),"=", alpha)
  
  if(from0.to1==T){
    barplot(energy.table$Rel.Energy,names.arg = energy.table$Level, cex.names =0.75, cex.axis = 0.75,
            main = 'Energy by decomposition level', horiz = TRUE, sub = sub, xlab='Energy (100%)', ylab='Crystal', xlim=c(0,1))
  }else{
    barplot(energy.table$Rel.Energy,names.arg = energy.table$Level, cex.names =0.75, cex.axis = 0.75,
            main = 'Energy by decomposition level', horiz = TRUE, sub = sub, xlab='Energy (100%)', ylab='Crystal')
  }
}

#======================================#
#  Wavelet Scalogram and Periodogram   #
#======================================#

frac.dwt.periodogram <- function(w,J){
  if(is.list(w)){
    w <- w$w
    }else{
      stop("w must be a list")
    }
  
  tab <- generate.dyadic(J = J)
  
  I.mat <- matrix(data = NA, nrow=length(w), ncol = nrow(tab)-1)
  I.w   <- NULL
  
  for(i in 1:ncol(I.mat)){
    I.mat[,i] <- c((w[(tab[i,"Start"]:tab[i,"End"])])^2, rep(NA, 2^J - length(w[(tab[i,"Start"]:tab[i,"End"])])   ))
    I.w       <- c(I.w, (w[(tab[i,"Start"]:tab[i,"End"])])^2)
    
  }
  I.mat <- as.matrix(I.mat)
  
  I.list <- list(I.mat = I.mat, I.w=I.w)
  
  return(I.list)
  
}

# Wavelet Scalogram

frac.dwt.scalogram <- function(w, I, J, plot = TRUE, perc=T){
  if(missing(w) & missing(I)){stop("You need to input either w or I")}
  if(missing(J)){stop("You need to input J!")}
  if(is.logical(plot)==FALSE){stop("plot is a logical variable")}
  
  if(missing(I)){
    I <- frac.dwt.periodogram(w = w, J=J)
  }else{
    I <- I$I.mat
  }
  
  
  S <- rep(NA, J)
  
  if(ncol(I)!=J){stop("There's a problem in the periodogram matrix.")}
  
  for(i in 1:ncol(I)){
    S[i] <- sum(I[which(!(is.na(I[,i]))),i])
  }
 
  if(perc == T){
    S <- S/sum(S)
  }
  
  if(plot==T){
    plot(x = 1:J, y=S, main = 'FrDWT Scalogram', 
         pch = 16, col = 'blue', type = 'b', cex = 1.25, 
         xlab= "Decomposition levels", ylab='S(j)')
  }
  
  S <- cbind(1:J, S)
  
  S <- as.data.frame(S)
  colnames(S) <- c("j", "Energy")
  return(S) 
}

#=========================================#
#  Wavelet bootstrap confidence interval  #
#=========================================#

resample.dwt <- function(w, J, alpha, type){
  #Function to resample the DWT coefficients for bootstrap
  #Internal use only
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  tab <- generate.dyadic(J)
  wn  <- rep(NA, length(w))
  for(i in 1:nrow(tab)){
    wd <- w[tab[i,'Start']:tab[i,'End']]
    wn[tab[i,'Start']:tab[i,'End']] <- sample(wd, replace = T)
  }
  
  new.signal <- Frac.idwt(w = w, J = J, type = type, alpha = alpha)
  return(new.signal)
}

bootstrap.ci <- function(signal, N.rep = 1000, type, alpha, J, probs=c(0.05, 0.95)){
  #Based on the proposed procedure of Bruce and Gao 1996
  #Compatible to the S-PLUS wavelet proposal of bootstrap CI 
  #Author: Mateus Gonzalez de Freitas Pinto
  #Sao Paulo, Brazil, 2022
  stats.fun <- function(x) quantile(x, probs = probs)
  w         <- Frac.dwt(signal = signal, J = J, type = type, alpha = alpha, out.mat = F)
  estimate  <- wavelet.shrinkage(w = w, J = J,type = type,alpha = alpha, reconstruct = T)$shrinked
  
  resid.dwt <- Frac.dwt(signal = signal-estimate,J = J, type = type, alpha = alpha, out.mat = F)
  
  samples <- matrix(0, length(signal), N.rep)
  
  for(i in 1:N.rep){
    new.resid     <- resample.dwt(w = resid.dwt, J = J,alpha = alpha, type = type)
    new.sign      <- new.resid+estimate
    wd            <- Frac.dwt(signal = new.sign, J = J, type = type, alpha = alpha, out.mat = F)
    new.estimate  <- wavelet.shrinkage(w = wd, J = J,type = type,alpha = alpha, reconstruct = T)$shrinked
    samples[,i]   <- new.estimate
  }
  
  mat <- apply(samples, 1, stats.fun)
  lower.ci <- mat[1,]
  upper.ci <- mat[2,]
  
  list.boot <- list(signal=signal, lower.ci=lower.ci, upper.ci=upper.ci, bootstrap.samples=samples)
  return(list.boot)
}

