#Arguments
# nsims: number of simulations to run (default is 1000)
# boots: number of bootstrap samples (default is 0, which doesn't run bootstrap analysis)
# MCsamples: number of monte carlo samples (default is 0, which does run Monte Carlo analysis)
# alpha: level of test, acceptable type I error rate (default is 0.05)
# N: sample size to test
# corm: correlation among repeated measurements of mediator
# cory: correlation among repeated measurements of outcome
# aest: a-path estimates standardized by SD(M)
# best: b-path estimate standardized by SD(M) and SD(Y)
# cest: c-path estimate standardized by SD(Y)
# dest: d-path estimate standardized by SD(M) and SD(Y)

WSmedpower <- function(nsims=1000, boots=0, MCsamples=0, alpha=0.05, N, corm, cory, aest, best, cest, dest){
  sigm1y1 <- best + dest
  sigm2y2 <- best - dest
  sigm1y2 <- corm*(best - dest)
  sigm2y1 <- corm*(best + dest)
  
  cormat <- matrix(c(1,       corm,    sigm1y1, sigm1y2,
                     corm,    1,       sigm2y1, sigm2y2,
                     sigm1y1, sigm2y1, 1,       cory, 
                     sigm1y2, sigm2y2, cory,    1),
                   nrow = 4, ncol = 4)
  dets <- rep(NA, times = 3)
  for (i in 2:4){
    dets[i-1] <- det(cormat[1:i, 1:i]) 
  }
  
  if(all(dets > 0)){
    Results <- data.frame(matrix(data = NA, nrow = nsims, ncol = 3))
    names(Results) <- c("JointSig", "Bootstrap", "MonteCarlo")
    for(i in 1:nsims){
      M1M2Y1Y2 <- data.frame(matrix(rnorm(n = N*4), ncol = 4)%*%chol(cormat))
      names(M1M2Y1Y2) <- c("M1", "M2", "Y1", "Y2")
      M1M2Y1Y2$M2 <- M1M2Y1Y2$M2-aest
      M1M2Y1Y2$Y2 <- M1M2Y1Y2$Y2-cest 
      contrast <- matrix(c( 1,  1,  0, 
                            -1,  1,  0,
                            0,  0,  1,
                            0,  0, -1), nrow = 4, byrow = TRUE) 
      MDiffMSumYDiff <- data.frame(as.matrix(M1M2Y1Y2)%*%contrast)
      names(MDiffMSumYDiff) <- c("Mdiff", "Msum", "Ydiff")
      MSumUncentered <- MDiffMSumYDiff$Msum 
      MDiffMSumYDiff$Msum <- MDiffMSumYDiff$Msum -mean(MDiffMSumYDiff$Msum)
      
      apath <- mean(MDiffMSumYDiff[,1])
      seapath <- sd(MDiffMSumYDiff[,1])/sqrt(N)
      apathsig <- 2*pt(abs(apath/seapath), df = N-1, lower.tail = FALSE)
      bmodel <- lm(Ydiff~Mdiff+Msum, data = MDiffMSumYDiff)
      betas <- coef(bmodel)
      bpath <- betas[2]
      sebpath <- summary(bmodel)$coef[2,2]
      bpathsig <- summary(bmodel)$coef[2,4]

      
      #Joint Significance
      JointSigReject <- ((apathsig < alpha)&(bpathsig < alpha))
      Results[i,1] <- JointSigReject
      
      #Bootstrap CI
      if(boots > 0){
             confidence <- 1- alpha
      bootResults <- vector(length = boots)
      for (j in 1:boots){
        rows <- sample(1:N, size = N, replace = TRUE)
        bootdat <- MDiffMSumYDiff[rows, ]
        bootapath <- mean(bootdat$Mdiff)
        bootmodel <- lm(Ydiff~Mdiff+Msum, data = bootdat)
        bootbpath <- coef(bootmodel)[2]
        bootResults[j] <- bootapath*bootbpath 
      }
      bootResults <- sort(bootResults)
      LowerCI = (1-confidence)/2
      UpperCI = 1 - LowerCI
      BootLCI = bootResults[ceiling(boots*LowerCI)] #something weird going on with ceiling function 25 --> 26
      BootUCI = bootResults[ceiling(boots*UpperCI)]
      BootReject = ((BootLCI > 0)|(BootUCI < 0))
      Results[i,2] <- BootReject 
      }

      
      #Monte Carlo CI
      if(MCsamples > 0){
      asamp <- rnorm(n = MCsamples, mean = apath, sd = seapath)
      bsamp <- rnorm(n = MCsamples, mean = bpath, sd = sebpath)
      absamp <- asamp*bsamp
      absamp <- sort(absamp)
      MCLCI <- absamp[ceiling(MCsamples*LowerCI)]
      MCUCI <- absamp[ceiling(MCsamples*UpperCI)]
      MCReject <- ((MCLCI >0)|(MCUCI <0))
      Results[i,3] <- MCReject
      }

    }
    
    testpower <- colMeans(Results)
    return(testpower)
    
  }
  if(all(dets <= 0)){
    print("Impossible Parameter Combinations")
  }
  
}


WSmedpower(nsims = 5000, boots = 0, MCsamples = 0, alpha = 0.05, N = 100, corm = .22, cory = .32, aest = -.44, best = -.34, cest = .39, dest = 0)



