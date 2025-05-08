
install.packages("bayesplot")
library(bayesplot)
library(spAbundance)
library(coda)
library(stars)
library(ggplot2)
library(readxl)
library(sf)

dataset <- read.csv("data.csv")
dataset

#creating abundance and detection covariates formulas and standarize them 

y <- dataset[11:13]


primunite <- as.numeric(factor(dataset$ï..unite.primaire))
view(primunite)
emveg <- scale(dataset$emveg)
emveg
wt <- dataset$Wt
wt
surfarea <- scale(dataset$Area..meters.square.)
surfarea
year <- as.numeric(factor(dataset$Year))
year
abund.cov <- data.frame(primunite = primunite, emveg = emveg, wt = wt, surfarea = surfarea, year = year)





effort <- scale(dataset[29:31])
effort
julianday <- dataset[32:34]
julianday
heure <- dataset[26:28]
heure

det.cov <- list(effort = effort, julianday = julianday, heure = heure)



datalist <- list(y=y, det.covs=det.cov, abund.covs=abund.cov)


# Format with explicit specification of inits for alpha and beta
# with three detection parameters and two abundance parameters
# (including the intercept).
inits <- list(alpha = 0,
              beta = 0,
              kappa = 0.5,
              sigma.sq.mu = 0.5,
              N = apply(y, 1, max, na.rm = TRUE))#to set the initial values to the largest observed count at each site


priors <- list(alpha.normal = list(mean = 0, var = 2.72),
               beta.normal = list(mean = 0, var = 100),
               kappa.unif = c(0, 100),
               sigma.sq.mu.ig = list(0.1, 0.1))#The prior you are using for the random effects variance (σ 2) is an 
#Inverse Gamma distribution with parameters α=0.1 (shape) and β=0.1 (scale). This is a very vague prior, allowing for 
#a wide range of values and potentially leading to the large uncertainty in your random effect variance.

priors.1 <- list(alpha.normal = list(mean = 0, var = 2.72),
                 beta.normal = list(mean = 0, var = 100),
                 kappa.unif = c(0, 100),
                 sigma.sq.mu.ig = list(2, 2))# The prior for the random effects variance (σ²) is an Inverse Gamma distribution 
# with shape α = 2 and scale β = 2. This is a mildly informative prior with mean = 2. 
# It reflects moderate prior belief about the scale of the random effect variance, 
# allowing for some shrinkage without being overly restrictive.

batch.length <- 50
n.batch <- 3100
# Total number of MCMC samples per chain
batch.length * n.batch


tuning <- list(beta = 0.5, alpha = 0.5, kappa = 0.5, beta.star = 0.5)
# accept.rate = 0.43 by default, so we do not specify it.


n.burn <- 10000
n.thin <- 10
n.chains <- 3


#------------------------------------------------------------------newdataset with no replicates of sites-----------------------

#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------- Abundances Estimation part without Habitat covariates--------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------




#---------------------N-Mixture Models; Poisson and NB distributions with informative prior in sigma -----------------------------

out.i1 <- NMix(abund.formula = ~(1 | primunite) + year, 
               det.formula = ~ julianday, data = datalist, inits = inits, priors = priors.1, 
               tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
               verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(out.i1)

out.i11 <- NMix(abund.formula = ~(1 | primunite) + year, 
                det.formula = ~ effort, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i11)


out.i12 <- NMix(abund.formula = ~(1 | primunite) + year, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i12)

out.i13 <- NMix(abund.formula = ~(1 | primunite), 
                det.formula = ~ effort, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i13)


out.i14 <- NMix(abund.formula = ~(1 | primunite), 
                det.formula = ~ heure + effort, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i14)

out.i15 <- NMix(abund.formula = ~(1 | primunite) + year, 
                det.formula = ~ heure + poly(heure,1), data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i15)

out.i16 <- NMix(abund.formula = ~(1 | primunite) + year, 
                det.formula = ~ effort + heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i16)

out.i17 <- NMix(abund.formula = ~(1 | primunite) + year, 
                det.formula = ~ effort + heure + poly(heure,1), data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i17)

out.i11nb <- NMix(abund.formula = ~(1 | primunite) + year, 
                  det.formula = ~ effort, data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i11nb)


out.i12nb <- NMix(abund.formula = ~(1 | primunite) + year, 
                  det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i12nb)

out.i13nb <- NMix(abund.formula = ~(1 | primunite), 
                  det.formula = ~ effort, data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i13nb)


out.i14nb <- NMix(abund.formula = ~(1 | primunite), 
                  det.formula = ~ heure + effort, data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i14nb)

out.i15nb <- NMix(abund.formula = ~(1 | primunite) + year, 
                  det.formula = ~ heure + poly(heure,1), data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i15nb)

out.i16nb <- NMix(abund.formula = ~(1 | primunite) + year, 
                  det.formula = ~ effort + heure, data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i16nb)

out.i17nb <- NMix(abund.formula = ~(1 | primunite) + year, 
                  det.formula = ~ effort + heure + poly(heure,1), data = datalist, inits = inits, priors = priors.1, 
                  tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                  verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary(out.i17nb)


#---------------------------Posterior N-samples for abundance estimation and 95% CI--------------------------------------------------------

N_samples.i1 <- out.i1$N.samples
summary(N_samples.i1)
N_means.i1 <- colMeans(as.matrix(N_samples.i1))
total_N_abundance_mean.i1 <- sum(N_means.i1)
total_N_abundance_mean.i1
# Convert coda object to matrix
N_matrix.i1 <- as.matrix(N_samples.i1)
# Sum across sites for each posterior sample
total_N_abundance_samples.i1 <- rowSums(N_matrix.i1)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i1 <- mean(total_N_abundance_samples.i1)
N.samples_credible_interval.i1 <- quantile(total_N_abundance_samples.i1, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i1, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i1, "\n")


N_samples.i11 <- out.i11$N.samples
summary(N_samples.i11)
N_means.i11 <- colMeans(as.matrix(N_samples.i11))
total_N_abundance_mean.i11 <- sum(N_means.i11)
total_N_abundance_mean.i11
# Convert coda object to matrix
N_matrix.i11 <- as.matrix(N_samples.i11)
# Sum across sites for each posterior sample
total_N_abundance_samples.i11 <- rowSums(N_matrix.i11)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i11 <- mean(total_N_abundance_samples.i11)
N.samples_credible_interval.i11 <- quantile(total_N_abundance_samples.i11, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i11, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i11, "\n")

N_samples.i12 <- out.i12$N.samples
summary(N_samples.i12)
N_means.i12 <- colMeans(as.matrix(N_samples.i12))
total_N_abundance_mean.i12 <- sum(N_means.i12)
total_N_abundance_mean.i12
# Convert coda object to matrix
N_matrix.i12 <- as.matrix(N_samples.i12)
# Sum across sites for each posterior sample
total_N_abundance_samples.i12 <- rowSums(N_matrix.i12)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i12 <- mean(total_N_abundance_samples.i12)
N.samples_credible_interval.i12 <- quantile(total_N_abundance_samples.i12, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i12, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i12, "\n")

N_samples.i13 <- out.i13$N.samples
summary(N_samples.i13)
N_means.i13 <- colMeans(as.matrix(N_samples.i13))
total_N_abundance_mean.i13 <- sum(N_means.i13)
total_N_abundance_mean.i13
# Convert coda object to matrix
N_matrix.i13 <- as.matrix(N_samples.i13)
# Sum across sites for each posterior sample
total_N_abundance_samples.i13 <- rowSums(N_matrix.i13)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i13 <- mean(total_N_abundance_samples.i13)
N.samples_credible_interval.i13 <- quantile(total_N_abundance_samples.i13, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i13, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i13, "\n")

N_samples.i14 <- out.i14$N.samples
summary(N_samples.i14)
N_means.i14 <- colMeans(as.matrix(N_samples.i14))
total_N_abundance_mean.i14 <- sum(N_means.i14)
total_N_abundance_mean.i14
# Convert coda object to matrix
N_matrix.i14 <- as.matrix(N_samples.i14)
# Sum across sites for each posterior sample
total_N_abundance_samples.i14 <- rowSums(N_matrix.i14)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i14 <- mean(total_N_abundance_samples.i14)
N.samples_credible_interval.i14 <- quantile(total_N_abundance_samples.i14, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i14, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i14, "\n")

N_samples.i15 <- out.i15$N.samples
summary(N_samples.i15)
N_means.i15 <- colMeans(as.matrix(N_samples.i15))
total_N_abundance_mean.i15 <- sum(N_means.i15)
total_N_abundance_mean.i15
# Convert coda object to matrix
N_matrix.i15 <- as.matrix(N_samples.i15)
# Sum across sites for each posterior sample
total_N_abundance_samples.i15 <- rowSums(N_matrix.i15)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i15 <- mean(total_N_abundance_samples.i15)
N.samples_credible_interval.i15 <- quantile(total_N_abundance_samples.i15, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i15, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i15, "\n")

N_samples.i16 <- out.i16$N.samples
summary(N_samples.i16)
N_means.i16 <- colMeans(as.matrix(N_samples.i16))
total_N_abundance_mean.i16 <- sum(N_means.i16)
total_N_abundance_mean.i16
# Convert coda object to matrix
N_matrix.i16 <- as.matrix(N_samples.i16)
# Sum across sites for each posterior sample
total_N_abundance_samples.i16 <- rowSums(N_matrix.i16)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i16 <- mean(total_N_abundance_samples.i16)
N.samples_credible_interval.i16 <- quantile(total_N_abundance_samples.i16, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i16, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i16, "\n")

N_samples.i17 <- out.i17$N.samples
summary(N_samples.i17)
N_means.i17 <- colMeans(as.matrix(N_samples.i17))
total_N_abundance_mean.i17 <- sum(N_means.i17)
total_N_abundance_mean.i17
# Convert coda object to matrix
N_matrix.i17 <- as.matrix(N_samples.i17)
# Sum across sites for each posterior sample
total_N_abundance_samples.i17 <- rowSums(N_matrix.i17)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i17 <- mean(total_N_abundance_samples.i17)
N.samples_credible_interval.i17 <- quantile(total_N_abundance_samples.i17, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i17, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i17, "\n")

N_samples.i11nb <- out.i11nb$N.samples
summary(N_samples.i11nb)
N_means.i11nb <- colMeans(as.matrix(N_samples.i11nb))
total_N_abundance_mean.i11nb <- sum(N_means.i11nb)
total_N_abundance_mean.i11nb
# Convert coda object to matrix
N_matrix.i11nb <- as.matrix(N_samples.i11nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i11nb <- rowSums(N_matrix.i11nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i11nb <- mean(total_N_abundance_samples.i11nb)
N.samples_credible_interval.i11nb <- quantile(total_N_abundance_samples.i11nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i11nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i11nb, "\n")

N_samples.i12nb <- out.i12nb$N.samples
summary(N_samples.i12nb)
N_means.i12nb <- colMeans(as.matrix(N_samples.i12nb))
total_N_abundance_mean.i12nb <- sum(N_means.i12nb)
total_N_abundance_mean.i12nb
# Convert coda object to matrix
N_matrix.i12nb <- as.matrix(N_samples.i12nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i12nb <- rowSums(N_matrix.i12nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i12nb <- mean(total_N_abundance_samples.i12nb)
N.samples_credible_interval.i12nb <- quantile(total_N_abundance_samples.i12nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i12nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i12nb, "\n")

N_samples.i13nb <- out.i13nb$N.samples
summary(N_samples.i13nb)
N_means.i13nb <- colMeans(as.matrix(N_samples.i13nb))
total_N_abundance_mean.i13nb <- sum(N_means.i13nb)
total_N_abundance_mean.i13nb
# Convert coda object to matrix
N_matrix.i13nb <- as.matrix(N_samples.i13nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i13nb <- rowSums(N_matrix.i13nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i13nb <- mean(total_N_abundance_samples.i13nb)
N.samples_credible_interval.i13nb <- quantile(total_N_abundance_samples.i13nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i13nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i13nb, "\n")

N_samples.i14nb <- out.i14nb$N.samples
summary(N_samples.i14nb)
N_means.i14nb <- colMeans(as.matrix(N_samples.i14nb))
total_N_abundance_mean.i14nb <- sum(N_means.i14nb)
total_N_abundance_mean.i14nb
# Convert coda object to matrix
N_matrix.i14nb <- as.matrix(N_samples.i14nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i14nb <- rowSums(N_matrix.i14nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i14nb <- mean(total_N_abundance_samples.i14nb)
N.samples_credible_interval.i14nb <- quantile(total_N_abundance_samples.i14nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i14nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i14nb, "\n")

N_samples.i15nb <- out.i15nb$N.samples
summary(N_samples.i15nb)
N_means.i15nb <- colMeans(as.matrix(N_samples.i15nb))
total_N_abundance_mean.i15nb <- sum(N_means.i15nb)
total_N_abundance_mean.i15nb
# Convert coda object to matrix
N_matrix.i15nb <- as.matrix(N_samples.i15nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i15nb <- rowSums(N_matrix.i15nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i15nb <- mean(total_N_abundance_samples.i15nb)
N.samples_credible_interval.i15nb <- quantile(total_N_abundance_samples.i15nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i15nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i15nb, "\n")

N_samples.i16nb <- out.i16nb$N.samples
summary(N_samples.i16nb)
N_means.i16nb <- colMeans(as.matrix(N_samples.i16nb))
total_N_abundance_mean.i16nb <- sum(N_means.i16nb)
total_N_abundance_mean.i16nb
# Convert coda object to matrix
N_matrix.i16nb <- as.matrix(N_samples.i16nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i16nb <- rowSums(N_matrix.i16nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i16nb <- mean(total_N_abundance_samples.i16nb)
N.samples_credible_interval.i16nb <- quantile(total_N_abundance_samples.i16nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i16nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i16nb, "\n")

N_samples.i17nb <- out.i17nb$N.samples
summary(N_samples.i17nb)
N_means.i17nb <- colMeans(as.matrix(N_samples.i17nb))
total_N_abundance_mean.i17nb <- sum(N_means.i17nb)
total_N_abundance_mean.i17nb
# Convert coda object to matrix
N_matrix.i17nb <- as.matrix(N_samples.i17nb)
# Sum across sites for each posterior sample
total_N_abundance_samples.i17nb <- rowSums(N_matrix.i17nb)
# Calculate mean and credible interval
total_N.samples_abundance_mean.i17nb <- mean(total_N_abundance_samples.i17nb)
N.samples_credible_interval.i17nb <- quantile(total_N_abundance_samples.i17nb, c(0.025, 0.975))
cat("Total Abundance Mean:", total_N.samples_abundance_mean.i17nb, "\n")
cat("95% Credible Interval:", N.samples_credible_interval.i17nb, "\n")



#---------------detection probability from the best model-----------------------------

# Extract posterior samples for detection intercept
post_samples <- out.i12$alpha.samples

# Extract posterior samples for detection intercept
logit_p <- out.i12$alpha.samples[, 1]  # assuming the intercept is the first column

# Define inverse logit function
inv_logit <- function(x) exp(x) / (1 + exp(x))

# Transform to probability scale
p_samples <- inv_logit(logit_p)

# Calculate mean and 95% credible interval
mean_p <- mean(p_samples)
ci_p <- quantile(p_samples, probs = c(0.025, 0.975))

# Print results
cat("Mean detection probability:", round(mean_p, 2), "\n")
cat("95% Credible Interval:", round(ci_p[1], 2), "-", round(ci_p[2], 2), "\n")

#------------------------------------------------------#




#--------------------------------------- PPC and AIC of new models without julianday different priors and without year as well ------


ppc.outi11 <- ppcAbund(out.i11, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi11)
ppc.outi12 <- ppcAbund(out.i12, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi12)
ppc.outi13 <- ppcAbund(out.i13, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi13)
ppc.outi14 <- ppcAbund(out.i14, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi14)
ppc.outi15 <- ppcAbund(out.i15, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi15)
ppc.outi16 <- ppcAbund(out.i16, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi16)
ppc.outi17 <- ppcAbund(out.i17, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi17)

aicouti11 <- waicAbund(out.i11)
aicouti11
aicouti12 <- waicAbund(out.i12)
aicouti12
aicouti13 <- waicAbund(out.i13)
aicouti13
aicouti14 <- waicAbund(out.i14)
aicouti14
aicouti15 <- waicAbund(out.i15)
aicouti15
aicouti16 <- waicAbund(out.i16)
aicouti16
aicouti17 <- waicAbund(out.i17)
aicouti17

ppc.outi11nb <- ppcAbund(out.i11nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi11nb)
ppc.outi12nb <- ppcAbund(out.i12nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi12nb)
ppc.outi13nb <- ppcAbund(out.i13nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi13nb)
ppc.outi14nb <- ppcAbund(out.i14nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi14nb)
ppc.outi15nb <- ppcAbund(out.i15nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi15nb)
ppc.outi16nb <- ppcAbund(out.i16nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi16nb)
ppc.outi17nb <- ppcAbund(out.i17nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.outi17nb)

aicouti11nb <- waicAbund(out.i11nb)
aicouti11nb
aicouti12nb <- waicAbund(out.i12nb)
aicouti12nb
aicouti13nb <- waicAbund(out.i13nb)
aicouti13nb
aicouti14nb <- waicAbund(out.i14nb)
aicouti14nb
aicouti15nb <- waicAbund(out.i15nb)
aicouti15nb
aicouti16nb <- waicAbund(out.i16nb)
aicouti16nb
aicouti17nb <- waicAbund(out.i17nb)
aicouti17nb



#-----------------------------------------relationship covariates with detection----------------------------------------------




datalist$det.covs$effort[is.na(datalist$det.covs$effort)] <- mean(datalist$det.covs$effort, na.rm = TRUE)
if (ncol(detection_probs_effort) != length(effort_range)) {
  detection_probs_effort <- t(detection_probs_effort)  # Transpose if needed
}

# Extract posterior samples
effort_samples <- out.i13$alpha.samples  # Detection coefficients
# Generate effort values
effort_range <- seq(min(datalist$det.covs$effort), max(datalist$det.covs$effort), length.out = 100)
# Compute detection probability for each effort value
detection_probs_effort <- apply(effort_samples, 1, function(a) plogis(a[1] + a[2] * effort_range))
# Plot
plot(effort_range, rowMeans(detection_probs_effort), type = "l", col = "blue", lwd = 2,
     xlab = "Survey Effort", ylab = "Detection Probability",
     main = "Effect of Effort on Detection")

polygon(c(effort_range, rev(effort_range)), 
        c(apply(detection_probs_effort, 2, quantile, probs = 0.975), 
          rev(apply(detection_probs_effort, 2, quantile, probs = 0.025))),
        col = rgb(0, 0, 1, 0.2), border = NA)


# Generate a sequence of heure values
heure_range <- seq(min(datalist$det.covs$heure, na.rm = TRUE), 
                   max(datalist$det.covs$heure, na.rm = TRUE), 
                   length.out = 100)

# Compute detection probability for each heure value
heure_detection_samples <- sapply(heure_range, function(h) {
  plogis(out.i12$alpha.samples[, "(Intercept)"] + out.i12$alpha.samples[, "heure"] * h)
})

heure_detection_samples <- t(heure_detection_samples)  # Fix the dimension issue
print(dim(heure_detection_samples))  # Should be (100, 43500)
mean_detection_heure <- apply(heure_detection_samples, 1, mean)  # Mean for each time step
lower_CI_heure <- apply(heure_detection_samples, 1, quantile, probs = 0.025)
upper_CI_heure <- apply(heure_detection_samples, 1, quantile, probs = 0.975)


print(length(mean_detection_heure))  # Should be 100
print(length(lower_CI_heure))  # Should be 100
print(length(upper_CI_heure))  # Should be 100



plot(heure_range, mean_detection_heure, type = "l", col = "blue", lwd = 2,
     xlab = "Time of Day (Heure)", ylab = "Detection Probability",
     ylim = c(min(lower_CI_heure), max(upper_CI_heure)),
     main = "Effect of Heure on Detection Probability")

polygon(c(heure_range, rev(heure_range)), 
        c(upper_CI_heure, rev(lower_CI_heure)), 
        col = rgb(0, 0, 1, 0.2), border = NA)







# Extract posterior samples for heure from the model
heure_samples <- out.i13$alpha.samples  # Ensure the model contains heure

# Generate heure values across the observed range
heure_range <- seq(min(datalist$det.covs$heure), max(datalist$det.covs$heure), length.out = 100)

# Compute detection probability for each heure value across posterior samples
detection_probs_heure <- apply(heure_samples, 1, function(a) plogis(a[1] + a[2] * heure_range))

# Compute 95% credible intervals
lower_CI_heure <- apply(detection_probs_heure, 1, quantile, probs = 0.025)
upper_CI_heure <- apply(detection_probs_heure, 1, quantile, probs = 0.975)

# Check that lengths match
print(length(heure_range))  # Should be 100
print(length(lower_CI_heure))  # Should be 100
print(length(upper_CI_heure))  # Should be 100

stopifnot(length(heure_range) == length(lower_CI_heure), 
          length(heure_range) == length(upper_CI_heure))

# Plot mean detection probability
plot(heure_range, rowMeans(detection_probs_heure), type = "l", col = "red", lwd = 2,
     xlab = "Hour of Day", ylab = "Detection Probability",
     main = "Effect of Hour on Detection")

# Add 95% CI shading
polygon(c(heure_range, rev(heure_range)), 
        c(upper_CI_heure, rev(lower_CI_heure)),
        col = rgb(1, 0, 0, 0.2), border = NA)

#-------------------------------------EXTRAPOLATION----------------------------------------------------------------------------


# Load necessary package
library(MASS)  # For uncertainty propagation (mvrnorm)

# Define known values
total_water_surface_morocco <- 730   # Replace with actual value
sum_water_surface_study_sites <- 46.75   # Replace with actual value

# Compute extrapolation metric
extrapolation_metric <- total_water_surface_morocco / sum_water_surface_study_sites

# Brood abundance estimate from your study
mean_brood_abundance <- 24.5  # Replace with your estimated mean
lower_CI <- 21                # Replace with lower 95% CI
upper_CI <- 34                # Replace with upper 95% CI

# Approximate standard deviation from CI
sd_brood_abundance <- (upper_CI - lower_CI) / (2 * 1.96)  

# Simulate national estimates using uncertainty propagation
set.seed(123)  # For reproducibility
n_simulations <- 10000  # Number of simulations

# Simulate abundance based on estimated uncertainty
simulated_brood_abundance <- rnorm(n_simulations, mean_brood_abundance, sd_brood_abundance)

# Apply extrapolation metric
simulated_national_population <- simulated_brood_abundance * extrapolation_metric

# Compute summary statistics
national_mean <- mean(simulated_national_population)
national_CI <- quantile(simulated_national_population, c(0.025, 0.975))  # 95% CI

# Display results
cat("Estimated breeding population at national level:", round(national_mean), "\n")
cat("95% CI:", round(national_CI[1]), "-", round(national_CI[2]), "\n")



#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------- effect of habitat covariates on the abundance  --------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------------------


mod.1 <- NMix(abund.formula = ~(1 | primunite) + year + emveg, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "Poisson", n.omp.threads = 1, 
              verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary (mod.1)

mod.2 <- NMix(abund.formula = ~(1 | primunite) + year + wt, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.2)

mod.3 <- NMix(abund.formula = ~(1 | primunite) + year + surfarea, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.3)

mod.4 <- NMix(abund.formula = ~(1 | primunite) + year + emveg + wt, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.4)
mod.5 <- NMix(abund.formula = ~(1 | primunite) + year + emveg + surfarea, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.5)
mod.6 <- NMix(abund.formula = ~(1 | primunite) + year + wt + surfarea, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.6)

mod.7 <- NMix(abund.formula = ~(1 | primunite) + year + wt + surfarea + emveg, 
              det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
              tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
              family = "Poisson", n.omp.threads = 1, verbose = TRUE, 
              n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.7)

library(dplyr)
library(ggplot2)

# Step 1: Extract beta samples from the model
beta_samples <- mod.7$beta.samples

# Step 2: Compute mean and 95% CI for each covariate
summary_stats <- apply(beta_samples, 2, function(x) {
  c(mean = mean(x), lower = quantile(x, 0.025), upper = quantile(x, 0.975))
})

# Step 3: Convert to data frame
effects_df <- as.data.frame(t(summary_stats))
effects_df$term <- rownames(effects_df)
rownames(effects_df) <- NULL

# Step 4: Drop 'year' if needed
effects_df <- effects_df %>% filter(term != "year")

# Step 5: Add readable labels and significance info
effects_df <- effects_df %>%
  mutate(
    Significant = ifelse(`lower.2.5%` > 0 | `upper.97.5%` < 0, "Yes", "No"),
    term_label = case_when(
      term == "(Intercept)" ~ "Dam",
      term == "wtLagoon" ~ "Lagoon",
      term == "wtlake" ~ "Lake",
      term == "wtmarsh" ~ "Marsh",
      term == "wtPond" ~ "Pond",
      term == "wtriver" ~ "River",
      term == "wtSaltpan" ~ "Saltpan",
      term == "emveg" ~ "Emergent vegetation",
      term == "surfarea" ~ "Surface area",
      TRUE ~ term
    )
  )
# Desired order from bottom to top
desired_order <- c(
  "Emergent vegetation",
  "Surface area",
  "Lagoon",
  "Lake",
  "Marsh",
  "Pond",
  "River",
  "Saltpan",
  "Dam"
)

# Convert to factor with specified order
effects_df$term_label <- factor(effects_df$term_label, levels = desired_order)

# Step 6: Create the plot
ggplot(effects_df, aes(y = term_label, mean, x = mean, color = Significant)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = `lower.2.5%` , xmax = `upper.97.5%`), height = 0.0) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Yes" = "tomato1", "No" = "gray60")) +
  coord_cartesian(xlim = c(-25, 13)) +
  labs(
    title = "Posterior estimates of abundance covariates",
    x = "Effect on log-abundance",
    y = NULL
  ) +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.2)
  )
#----------------------------------------------------------------------------------------------------

ppc.mod.1 <- ppcAbund(mod.1, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.1)
ppc.mod.2 <- ppcAbund(mod.2, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.2)
ppc.mod.3 <- ppcAbund(mod.3, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.3)
ppc.mod.4 <- ppcAbund(mod.4, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.4)
ppc.mod.5 <- ppcAbund(mod.5, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.5)
ppc.mod.6 <- ppcAbund(mod.6, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.6)
ppc.mod.7 <- ppcAbund(mod.7, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.7)

aicmod.1 <- waicAbund(mod.1)
aicmod.1
aicmod.2 <- waicAbund(mod.2)
aicmod.2
aicmod.3 <- waicAbund(mod.3)
aicmod.3
aicmod.4 <- waicAbund(mod.4)
aicmod.4
aicmod.5 <- waicAbund(mod.5)
aicmod.5
aicmod.6 <- waicAbund(mod.6)
aicmod.6
aicmod.7 <- waicAbund(mod.7)
aicmod.7


mod.1nb <- NMix(abund.formula = ~(1 | primunite) + year + emveg, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, family = "NB", n.omp.threads = 1, 
                verbose = TRUE, n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

summary (mod.1nb)

mod.2nb <- NMix(abund.formula = ~(1 | primunite) + year + wt, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.2nb)

mod.3nb <- NMix(abund.formula = ~(1 | primunite) + year + surfarea, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.3nb)

mod.4nb <- NMix(abund.formula = ~(1 | primunite) + year + emveg + wt, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.4nb)
mod.5nb <- NMix(abund.formula = ~(1 | primunite) + year + emveg + surfarea, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.5nb)
mod.6nb <- NMix(abund.formula = ~(1 | primunite) + year + wt + surfarea, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.6nb)

mod.7nb <- NMix(abund.formula = ~(1 | primunite) + year + wt + surfarea + emveg, 
                det.formula = ~ heure, data = datalist, inits = inits, priors = priors.1, 
                tuning = tuning, n.batch = n.batch, batch.length = batch.length, 
                family = "NB", n.omp.threads = 1, verbose = TRUE, 
                n.report = 400, n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)
summary(mod.7nb)

ppc.mod.1nb <- ppcAbund(mod.1nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.1nb)
ppc.mod.2nb <- ppcAbund(mod.2nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.2nb)
ppc.mod.3nb <- ppcAbund(mod.3nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.3nb)
ppc.mod.4nb <- ppcAbund(mod.4nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.4nb)
ppc.mod.5nb <- ppcAbund(mod.5nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.5nb)
ppc.mod.6nb <- ppcAbund(mod.6nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.6nb)
ppc.mod.7nb <- ppcAbund(mod.7nb, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.mod.7nb)

aicmod.1nb <- waicAbund(mod.1nb)
aicmod.1nb
aicmod.2nb <- waicAbund(mod.2nb)
aicmod.2nb
aicmod.3nb <- waicAbund(mod.3nb)
aicmod.3nb
aicmod.4nb <- waicAbund(mod.4nb)
aicmod.4nb
aicmod.5nb <- waicAbund(mod.5nb)
aicmod.5nb
aicmod.6nb <- waicAbund(mod.6nb)
aicmod.6nb
aicmod.7nb <- waicAbund(mod.7nb)
aicmod.7nb

#---------------------------------BEST MODEL IS MOD.7----------------------------------------------------------------------


#-------------------------------------------------explore the effect of emergent vegetation on abundance--------------------


# Extract posterior samples
intercept_samples <- mod.7$beta.samples[, "(Intercept)"]
emveg_samples <- mod.7$beta.samples[, "emveg"]

# Standardized range used in model
emveg_range <- seq(
  min(mod.7$X[, "emveg"], na.rm = TRUE),
  max(mod.7$X[, "emveg"], na.rm = TRUE),
  length.out = 100
)

# Compute abundance samples (lambda) for each emveg value
lambda_samples <- sapply(emveg_range, function(em) {
  exp(intercept_samples + emveg_samples * em)
})

# Compute posterior summaries
mean_abundance <- apply(lambda_samples, 2, mean)
lower_CI <- apply(lambda_samples, 2, quantile, probs = 0.025)
upper_CI <- apply(lambda_samples, 2, quantile, probs = 0.975)

# Unstandardize emveg for x-axis labeling (if needed)
emveg_mean <- mean(datalist$abund.covs$emveg, na.rm = TRUE)
emveg_sd <- sd(datalist$abund.covs$emveg, na.rm = TRUE)
emveg_unscaled <- emveg_range * emveg_sd + emveg_mean


# Create the plotting dataframe
plot_df <- data.frame(
  emveg = emveg_range,
  mean = mean_abundance,
  lower = lower_CI,
  upper = upper_CI
)

# Generate the plot
ggplot(plot_df, aes(x = emveg, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lavender", alpha = 4) +
  geom_line(color = "blue", size = 0.8) +
  labs(
    x = "Emergent vegetation (standardized)",
    y = "Expected abundance",
    title = "Effect of Emergent Vegetation on Abundance"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey90")
  )

#----------------------------wt----------------------------------------------------------


wt_types <- c("Dam", "Lagoon", "Lake", "Marsh", "Pond", "River", "Saltpan")

# Mapping wetland type to beta column names
wt_effect_names <- c(
  "Dam" = NA,  # baseline level
  "Lagoon" = "wtLagoon",
  "Lake" = "wtlake",
  "Marsh" = "wtmarsh",
  "Pond" = "wtPond",
  "River" = "wtriver",
  "Saltpan" = "wtSaltpan"
)


# Set values for other covariates
emveg_val <- 0      # standardized
surfarea_val <- 0   # standardized
year_val <- 0       # standardized

lambda_matrix <- matrix(NA, nrow = nrow(beta_samples), ncol = length(wt_types))
colnames(lambda_matrix) <- wt_types

for (i in seq_along(wt_types)) {
  wt <- wt_types[i]
  wt_col <- wt_effect_names[wt]
  
  # Apply wetland type effect
  wt_term <- if (is.na(wt_col)) {
    0
  } else {
    beta_samples[, wt_col]
  }
  
  log_lambda <- beta_samples[, "(Intercept)"] +
    year_val * beta_samples[, "year"] +
    emveg_val * beta_samples[, "emveg"] +
    surfarea_val * beta_samples[, "surfarea"] +
    wt_term
  
  lambda_matrix[, i] <- exp(log_lambda)
}


abundance_summary <- data.frame(
  Wetland = wt_types,
  Mean = apply(lambda_matrix, 2, mean),
  Lower = apply(lambda_matrix, 2, quantile, probs = 0.025),
  Upper = apply(lambda_matrix, 2, quantile, probs = 0.975)
)

# Mark significant wetland types manually
abundance_summary$Significant <- ifelse(abundance_summary$Wetland %in% c("Marsh", "Pond", "Saltpan"), "Yes", "No")


ggplot(abundance_summary, aes(x = Wetland, y = Mean, fill = Significant)) +
  geom_bar(stat = "identity", color = "black", width = 0.6) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 0.8) +
  scale_fill_manual(values = c("Yes" = "#F8766D", "No" = "grey70")) +
  scale_y_log10()+
  labs(
    x = "Wetland Type",
    y = "Expected Abundance (λ)",
    fill = "Significant",
    title = "Effect of Wetland Type on Expected Abundance"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )

# Make sure the original wetland type column is present
lambda_df$wt <- datalist$abund.covs$wt

ggplot(lambda_df, aes(x = wt, y = lambda, fill = wt)) +
  geom_boxplot(
    alpha = 0.7,
    outlier.shape = 16,
    outlier.size = 1,
    outlier.alpha = 0.4
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  scale_y_log10() +
  labs(
    x = "Wetland Type",
    y = "Posterior Mean of Expected Abundance (λ)",
    title = "Expected Abundance by Original Wetland Types (Log Scale)"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none")




#---------------------------------------------------------------------------------------------------------------------------

# Create an empty matrix to store abundance estimates
lambda_samples <- matrix(NA, nrow = nrow(mod.7$beta.samples), ncol = length(emveg_range))

# Compute abundance estimates iteratively
for (i in 1:length(emveg_range)) {
  lambda_samples[, i] <- exp(intercept_samples + emveg_samples * emveg_range[i])
}

# Compute mean and credible intervals
mean_abundance <- apply(lambda_samples, 2, mean)
lower_CI <- apply(lambda_samples, 2, quantile, probs = 0.025)
upper_CI <- apply(lambda_samples, 2, quantile, probs = 0.975)

plot(emveg_range, mean_abundance, type = "l", col = "blue", lwd = 2,
     xlab = "Emergent Vegetation (Standardized)", ylab = "Estimated Abundance",
     main = "Effect of Emergent Vegetation on Abundance",
     ylim = c(min(lower_CI), max(upper_CI)))

# Add credible interval shading
polygon(c(emveg_range, rev(emveg_range)), 
        c(lower_CI, rev(upper_CI)),  
        col = rgb(0, 0, 1, 0.2), border = NA)

# Add mean line
lines(emveg_range, mean_abundance, col = "blue", lwd = 2)

print(range(mean_abundance))
print(range(lower_CI))
print(range(upper_CI))

hist(lambda_samples, breaks = 50, main = "Distribution of Predicted Abundance",
     xlab = "Predicted Abundance", col = "lightblue", border = "black")


plot(emveg_range, mean_abundance, type = "l", col = "blue", lwd = 2,
     xlab = "Emergent Vegetation (Standardized)", ylab = "Estimated Abundance",
     main = "Effect of Emergent Vegetation on Abundance",
     ylim = c(0.01, max(upper_CI)), log = "y")  # Log-scale Y-axis

# Add credible interval shading
polygon(c(emveg_range, rev(emveg_range)), 
        c(lower_CI, rev(upper_CI)),  
        col = rgb(0, 0, 1, 0.2), border = NA)

# Add mean line
lines(emveg_range, mean_abundance, col = "blue", lwd = 2)

# Scatter plot of raw data points
plot(datalist$abund.covs$emveg, apply(mod.1$N.samples, 2, mean), 
     xlab = "Emergent Vegetation (%)", ylab = "Observed Abundance",
     main = "Observed Abundance vs. Emergent Vegetation", 
     col = "black", pch = 16)




# Extract posterior samples
beta_samples_wt <- mod.2$beta.samples  # Extract abundance coefficients
intercept <- beta_samples_wt[, "(Intercept)"]  # Intercept (reference wetland type)

# Identify wetland type columns in the model output
wt_columns <- grep("^wt", colnames(beta_samples_wt), value = TRUE)  # Select wetland type variables
wt_effects <- beta_samples_wt[, wt_columns, drop = FALSE]  # Get coefficients for each wetland type

# Compute abundance for each wetland type
lambda_samples_wt <- exp(sweep(wt_effects, 1, intercept, "+"))

# Compute mean and credible intervals
mean_abundance_wt <- apply(lambda_samples_wt, 2, mean)
lower_CI_wt <- apply(lambda_samples_wt, 2, quantile, probs = 0.025)
upper_CI_wt <- apply(lambda_samples_wt, 2, quantile, probs = 0.975)

# Set up wetland type names
wt_names <- sub("wt", "", colnames(lambda_samples_wt))  # Remove "wt" prefix

# Define reasonable upper limit (e.g., 1.5x the 97.5% percentile of non-outliers)
y_limit <- min(max(upper_CI_wt), max(mean_abundance_wt) * 1.5)


# Define colors for each wetland type
wt_colors <- c("red", "blue", "green", "purple", "orange", "brown", "pink")  # Adjust as needed

# Create bar plot and capture bar midpoints
bar_midpoints <- barplot(mean_abundance_wt, names.arg = wt_names, col = wt_colors,
                         ylim = c(0.01, max(upper_CI_wt) * 1.5), log = "y", 
                         xlab = "Wetland Type", ylab = "Estimated Abundance (log scale)", 
                         main = "Effect of Wetland Type on Abundance", border = "black")

# Compute correct CI positioning: Ensure midpoint aligns with the bar top
ci_offset <- (upper_CI_wt - lower_CI_wt) / 2  # Half the CI range to center it

# Prevent disappearing CIs:
lower_CI_fixed <- pmax(mean_abundance_wt - ci_offset, 0.01)  # Ensure lower bound is within log-scale range
upper_CI_fixed <- pmin(mean_abundance_wt + ci_offset, max(upper_CI_wt) * 1.5)  # Keep upper bound within plot

# Adjust arrows to center the CI around the bar top
arrows(x0 = bar_midpoints, y0 = lower_CI_fixed, 
       x1 = bar_midpoints, y1 = upper_CI_fixed,
       angle = 90, code = 3, length = 0.05, col = "black", lwd = 2)

