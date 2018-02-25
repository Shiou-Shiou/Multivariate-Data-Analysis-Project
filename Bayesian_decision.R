# Final Project : Shiou-Shiou Deng
# the true probabilities pij were used as the prior means 
priorfunction <- function(scenario){
	if (scenario == "A") {
		#the prior probabilities of the scenario A 
		prior <- matrix(c(0.04, 0.1, 0.16, 0.22, 
							0.08, 0.14, 0.2, 0.26,
							0.12, 0.18, 0.24, 0.3, 
							0.16, 0.22, 0.28, 0.34), nrow = 4, byrow=TRUE)	
	}
	if (scenario == "B"){
		# the prior probabilities of the scenario B
		prior <- matrix(c(0.02, 0.05, 0.08, 0.11, 
				0.04, 0.07, 0.1, 0.13,
				0.06, 0.09, 0.12, 0.15, 
				0.08, 0.11, 0.14, 0.17), nrow = 4, byrow=TRUE)
	}
	if (scenario == "C"){
		# the prior probabilities of the scenario C 
		prior <- matrix(c(0.10, 0.25, 0.40, 0.55, 
              0.2, 0.35, 0.50, 0.65,
              0.3, 0.45, 0.60, 0.75, 
              0.4, 0.55, 0.70, 0.85), nrow = 4, byrow=TRUE)
	}
	if (scenario == "D"){
		# the prior probabilities of the scenario D  
		prior <- matrix(c(0.44, 0.50, 0.56, 0.62, 
               0.48, 0.54, 0.60, 0.66,
               0.52, 0.58, 0.64, 0.70, 
               0.56, 0.62, 0.68, 0.74), nrow = 4, byrow=TRUE)
	}
	if (scenario == "E"){
		# the prior probabilities of the scenario E  
		prior <- matrix(c(0.08, 0.09, 0.10, 0.11, 
              0.18, 0.19, 0.20, 0.21,
              0.28, 0.29, 0.30, 0.31, 
              0.29, 0.30, 0.31, 0.41), nrow = 4, byrow=TRUE)
	}
	if (scenario == "F"){
		# the prior probabilities of the scenario F  
		prior <- matrix(c(0.12, 0.16, 0.44, 0.50, 
              0.13, 0.18, 0.45, 0.52,
              0.14, 0.20, 0.46, 0.54, 
              0.15, 0.22, 0.47, 0.55), nrow = 4, byrow=TRUE)
	}
	if (scenario == "G"){
		# the prior probabilities of the scenario G  
		prior <- matrix(c(0.01, 0.04, 0.6, 0.10, 
              0.02, 0.1, 0.15, 0.30,
              0.03, 0.15, 0.30, 0.50, 
              0.04, 0.20, 0.45, 0.80), nrow = 4, byrow=TRUE)
	}
	return(prior)
}
# sensitivity analysis
priorchoose <- function(adjust, prior){
	if (adjust == "none") prior <- prior
	if (adjust == "subtract"){
		# substract 1-5 % points from the true DLT probabilities
		prior <- prior - matrix(runif(16, 0.01, 0.05), nrow = 4, ncol = 4)
		prior <- ifelse(prior < 0.01, 0.01, prior)
	}
	if (adjust == "add"){
		# add 1- 5 % points from the true DLT probabilities
		prior <- prior + matrix(runif(16, 0.01, 0.05), nrow = 4, ncol = 4)	
	}
	if (adjust == "both"){
	# add or substract 1-5 % points from the true DLT probabilities
		prior <- prior + matrix(runif(16, c(-0.05,-0.01), c(0.01, 0.05)), nrow = 4, ncol = 4)
		prior <- ifelse(prior < 0.01, 0.01, prior)
	}
	return(prior)
}
## Part one: escalate the dose of one agent until the first toxicity or the highest dose combination
firstToxicity <- function(patients, doses, prior){
    toxicity <- 0
    i <- 1
    j <- 1
    while (toxicity == 0 & i <= 4 & j <= 4){
        # toxicity is bernoulli distribution
        p <- rbinom(n=1, size = 1, prob = prior[i,j])
        # if toxicity is observed, then start to Bayescian method.
        if (p == 1) {
            toxicity <- 1
            patients[i,j] <- patients[i,j] + 1
            doses[i,j] <- 1 
        } else {
        # if toxicity is not observed, then try the higher dose unless the higest combination is reached
            # if agentA and agentB are smaller than 4, randomly choose one to increase
            patients[i,j] <- patients[i,j] + 1 
            if(i < 4 & j < 4) ifelse(runif(1)>=0.5, i <- i + 1, j <- j + 1) 
            # if one is already 4, then increase the other one
            if(i == 4 & j < 4) j <- j + 1 
            if(i < 4 & j == 4) i <- i + 1
            # reach the highest combination
            if(i == 4 & j == 4) patients[i,j] <- 1
        }
    }
    # give the n patients to find the first toxicity
    return(cbind(patients,doses,c(i,j,0,0)))
} 
updateDistriburion <- function(tij, n = 1, tp = p, a0, b0){
    # tij is the patient treated as dose i and dose j and p is toxicity or not.
    # n patients were treated, and tp of them experienced DLTS
    if ( n %% 1 != 0 | tp %% 1 != 0 | n <= 0 | tp < 0) stop("N and tp should be positive integers.")
    # no toxicity, update bij
    for (i in 1:tij[1]){
        for (j in 1:tij[2]){
            b0[i,j] <- b0[i,j] + (n - tp)          
        }
    }
    # toxicity, update aij
    for (i in tij[1]:4){
        for (j in tij[2]:4){
            a0[i,j] <- a0[i,j] + tp           
        }
    }
    return(cbind(a0, b0))
} 
# decide use skip rule or not 
skipchoose <- function(skip, a0, b0, alpha, theta, eta, tij){
	if (skip== TRUE){
		ratio <- a0/(a0+b0)
		expectation <- (-(alpha + eta))*(theta * pbeta(theta, a0, b0)- ratio * pbeta(theta, a0+1, b0))- eta*(ratio-theta)
		# Dose skipping is allowed
		tij <- which(expectation == max(expectation), arr.ind = TRUE)
	} else {
		ratio <- a0/(a0+b0)
		expectation <- (-(alpha + eta))*(theta * pbeta(theta, a0, b0)- ratio * pbeta(theta, a0+1, b0))- eta*(ratio-theta)
		# Dose skipping is not allowed
		doseless <- as.vector(expectation[1:tij[1], 1:tij[2]])
		#plus 1 to agentA or agentB
		doseplusA <- ifelse(tij[1] <= 3, expectation[tij[1]+1, tij[2]], -Inf)
		doseplusB <- ifelse(tij[2] <= 3, expectation[tij[1], tij[2]+1], -Inf)
		# find the remain dose has the max value
		tij <- which(expectation == max(c(doseless, doseplusA, doseplusB)), arr.ind = TRUE)  
	}
	return(tij)
}
MTD <- function(trails, scenario = "A", adjust = "none", skip = TRUE, nmin = 10, nmax = 50, theta = 0.2, delta = 0.05, r1 = 0.5, r2 = 0.95, eta = 1, alpha = 1.2, C = 4){
	if ( nmin %% 1 != 0 | nmax %% 1 != 0 | trails %% 1 != 0) stop("Nmin, nmax, and trails should be positive integers.")
	if ( theta < 0 | delta < 0  | r1 < 0 | r2 < 0 | theta > 1 | delta > 1  | r1 > 1 |  r2 > 1) stop("Theta, delta, r1, and r2 should be between 0 and 1.")
	if ( is.logical(skip) == FALSE) stop("Skip should be logical value.")
	if ( eta < 0 | alpha < 0) stop("Eta and alpha should be positive.")
	if ( is.numeric(C) == FALSE) stop("C should be number.")
    set.seed(267)
    noMTD <- numeric(1)
    MTD <- matrix(0, nrow=4, ncol=4)    
    samplesize <- numeric(1)
	prior <- priorfunction(scenario) 
	prior <- priorchoose(adjust, prior = prior)
    for (i in 1:trails){
		# (aij, bij) pairs such that aij = c*uij and bij*(1-uij) where c = 4 
        a0 <- C * prior 
        b0 <- C * (1 - prior)
        patients <- matrix(0, nrow=4, ncol=4)
        doses <- matrix(0, nrow=4, ncol=4)
        resultT <- firstToxicity(patients, doses, prior= prior)
        # the first toxicity
        patients <- resultT[,1:4]
        doses <- resultT[,5:8]
        tij <- resultT[1:2,9]
        # update the distribution
        resultP <- updateDistriburion(tij=tij, tp=1, a0=a0, b0=b0) 
        a0 <- resultP[,1:4]
        b0 <- resultP[,5:8] 
        # start to apply the bayesian discion theorey
        repeat{
            ntp <- sum(apply(patients, 2, sum))
            # check the S3 rules
            if ((1 - pbeta(theta + delta, a0[1,1], b0[1,1])) > r1 & ntp >= nmin){
                noMTD <- noMTD + 1
                break 
            } 
            # check the S4 rules
            higherdose <- as.vector(1 - pbeta(theta + delta, a0[tij[1]:4, tij[2]:4], b0[tij[1]:4,tij[2]:4]))
            if (!any(is.na(min(higherdose[2:length(higherdose)]))) & min(higherdose[2:length(higherdose)]) > r2 & ntp >= nmin) {
                MTD[tij[1], tij[2]] <- MTD[tij[1], tij[2]] + 1
                break
            }
            # check S1 and S2 rules            
            if (ntp == nmax & ntp >= nmin ) {
                MTD[tij[1], tij[2]] <- MTD[tij[1], tij[2]] + 1
                break
            }
            # Use the OSLA strategy
            # check Dose skipping is allowed or not and return the next dose
            tij <- skipchoose(skip=skip, a0=a0, b0=b0, alpha=alpha, eta = eta, theta = theta, tij=tij)
			# add 1 patient
			patients[tij[1],tij[2]] <- patients[tij[1],tij[2]] + 1
			# check if the patient have toxicity or not
            p <- rbinom(n=1, size=1, prob = prior[tij[1],tij[2]])
            doses[tij[1],tij[2]] <- doses[tij[1],tij[2]] + p 
			# update aij and bij
            resultP <- updateDistriburion(tij = tij, tp = p, a0=a0, b0=b0) 
            a0 <- resultP[,1:4]
            b0 <- resultP[,5:8]   
        }  
        samplesize <- samplesize + ntp
    }    
    averagesize <- samplesize/trails
	MTD <- MTD / trails
	noMTD <- noMTD / trails
    return(cbind(MTD, c(averagesize, noMTD,0,0)))
}
