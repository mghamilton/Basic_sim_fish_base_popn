#Matthew Hamilton
#26/04/19

#Basic simulation of base population from founders (one sampled founder population, no fixed effects,
#no half sibs, no inbreeding, no inbreeding depression, no missing data, 
#same number of indivduals per family etc)

###############################################################################################
#Inputs
###############################################################################################

n_fams          <- 20 #set number of full-sib families in base population
n_indiv_per_fam <- 50 #number of individuals per family

sim_reps <- 100 #number of times simulation is repeated

#true genetic parameters 
true_additive_var_1    <- 0.3
true_additive_var_2    <- 0.5
true_additive_cor_1_2  <- 0.8

#Note that dominance is not fitted in the model for asreml analyses below as the additive and dominance effects cannot be partitioned in this case (i.e. no half sibs)
#True parameter values can be entered to examine their impact on estimates of additive effects
#Enter a very small value in the place of zero (e.g. 0.0000001)
true_dominance_var_1   <- 0.2 
true_dominance_var_2   <- 0.0000001
true_dominance_cor_1_2 <- 0.0000001

true_residual_var_1    <- 0.5
true_residual_var_2    <- 0.3
true_residual_cor_1_2  <- 0.8

###############################################################################################
#load packages
###############################################################################################
library(asreml)
library(nadiv)

###############################################################################################
#Generate pedigree file
###############################################################################################

#founders 
ped <- data.frame(INDIV = 1:(2*n_fams),
                                    SIRE = NA,
                                    DAM = NA)

#base population (G0)
tmp <- data.frame(INDIV =  (max(ped$INDIV) + 1):(max(ped$INDIV) + n_fams * n_indiv_per_fam),
                  SIRE  =  rep(1:n_fams,              each = n_indiv_per_fam),
                  DAM   =  rep((n_fams+1):(2*n_fams), each = n_indiv_per_fam))

ped <- rbind(ped, tmp)
rm(tmp)

ped$INDIV <- as.factor(ped$INDIV) 
ped$SIRE  <- as.factor(ped$SIRE)
ped$DAM   <- as.factor(ped$DAM)

#covariances from correlations
true_additive_covar_1_2  <- true_additive_cor_1_2  * sqrt(true_additive_var_1  * true_additive_var_2)
true_dominance_covar_1_2 <- true_dominance_cor_1_2 * sqrt(true_dominance_var_1 * true_dominance_var_2)
true_residual_covar_1_2  <- true_residual_cor_1_2  * sqrt(true_residual_var_1  * true_residual_var_2)

# Define covariance matrices for random effects
Ga <- matrix(c(true_additive_var_1,  true_additive_covar_1_2,  true_additive_covar_1_2,  true_additive_var_2),  2, 2)
Gd <- matrix(c(true_dominance_var_1, true_dominance_covar_1_2, true_dominance_covar_1_2, true_dominance_var_2), 2, 2)
R <- matrix(c(true_residual_var_1,  true_residual_covar_1_2,  true_residual_covar_1_2,  true_residual_var_2),  2, 2)

#generate additive and dominance relationship matrices
#http://finzi.psych.upenn.edu/library/nadiv/html/warcolak.html
A_matrix <- makeA(ped)
A_matrix_inv <- asreml.Ainverse(dat[,1:3])$ginv  

D_matrix     <- makeD(ped)
D_matrix_inv <- D_matrix$listDinv
D_matrix     <- D_matrix$D

###############################################################################################
#Begin simulations
###############################################################################################
dat <- NULL
all_out <- NULL

for(rep in 1:sim_reps) {
  
###############################################################################################
#Generate phenotypes for base population (two traits from random from a multivariate normal distribution)
###############################################################################################

rm(dat)  
  
#True additive values for individuals
samp_a <- grfx(nrow(ped), G = Ga, incidence = A_matrix)
colnames(samp_a) <- c("TRAIT_1_TRUE_ADDITIVE",  "TRAIT_2_TRUE_ADDITIVE")
dat <- cbind(ped,samp_a)

#True dominance values for individuals
samp_d <- grfx(nrow(ped), G = Gd, incidence = D_matrix)
colnames(samp_d) <- c("TRAIT_1_TRUE_DOMINANCE", "TRAIT_2_TRUE_DOMINANCE")
dat <- cbind(dat,samp_d)

#True residual values for individuals
samp_r <- grfx(nrow(ped), G = R, incidence = Diagonal(nrow(ped))) # incidence = identity matrix
colnames(samp_r) <- c("TRAIT_1_TRUE_RESIDUAL", "TRAIT_2_TRUE_RESIDUAL")
dat <- cbind(dat, samp_r)

#Phenotypes for individuals
dat$TRAIT_1_PHENOTYPE <- dat$TRAIT_1_TRUE_ADDITIVE + dat$TRAIT_1_TRUE_DOMINANCE + dat$TRAIT_1_TRUE_RESIDUAL
dat[is.na(dat[,"SIRE"]),"TRAIT_1_PHENOTYPE"] <- NA #remove phenotype data for founders
dat$TRAIT_2_PHENOTYPE <- dat$TRAIT_2_TRUE_ADDITIVE + dat$TRAIT_2_TRUE_DOMINANCE + dat$TRAIT_2_TRUE_RESIDUAL
dat[is.na(dat[,"SIRE"]),"TRAIT_2_PHENOTYPE"] <- NA #remove phenotype data for founders

###############################################################################################
#Run ASReml and extract parameter estimates
###############################################################################################

add_start_values=c(true_additive_cor_1_2,true_additive_var_1,true_additive_var_2)
names(add_start_values)=c(true_additive_cor_1_2,true_additive_var_1,true_additive_var_2)

dom_start_values=c(true_dominance_cor_1_2,true_dominance_var_1,true_dominance_var_2)
names(dom_start_values)=c(true_dominance_cor_1_2,true_dominance_var_1,true_dominance_var_2)  

res_start_values <- c(true_residual_var_1,true_residual_covar_1_2,true_residual_var_2)
names(res_start_values) <- c(true_residual_var_1,true_residual_covar_1_2,true_residual_var_2)

model <- "asreml (cbind(TRAIT_1_PHENOTYPE, TRAIT_2_PHENOTYPE) ~ trait,
                  random = ~ ped(INDIV):  corgh(trait, init = add_start_values),
                  ginverse = list(INDIV=A_matrix_inv),
                  rcov   = ~ units:us(trait,init=res_start_values),
                  data=dat,
                  workspace=64e06,
                  maxiter=100,"

#get inital values in first rep of simulation
  if(rep == 1) {
    gammas <- eval(parse(text=paste(model,"start.values = TRUE)")))$gammas.table
    gammas[gammas[,"Constraint"] == "U","Constraint"] <- "P" #overwrite constraints 
  } 
dat.asr <- eval(parse(text=paste(model,"R.param = gammas)")))

#Undertake further iterations using parameter estimates from the intial model as starting values (this is to ensure that parameter estimates have stabilised/converged)
for (i in 1:10) {
  loglik <- summary(dat.asr)$loglik #get log likelihood value from last repetition
  dat.asr <- update(dat.asr) 
  if(signif(loglik,8) == signif(summary(dat.asr)$loglik,8)) {break} #break when likelihood value equal to previous repetition
  if(i==10) {stop} #if after 10 attempts it has not convergered then stop executing script
}

#Estimate heritability and SE
h2_est_1 <- nadiv:::pin(dat.asr, h2 ~ V2 / (V2 + V5))
h2_est_2 <- nadiv:::pin(dat.asr, h2 ~ V3 / (V3 + V7))

#Output variance components
var_comp <- summary(dat.asr)$varcomp
var_comp <- var_comp[-4,] # remove R!variance row
rownames(var_comp) <- c("add_cor_est_1_2","add_var_est_1","add_var_est_2","res_var_est_1","res_cov_est_1_2","res_var_est_2")

############################################################################################
#Generate files to output
############################################################################################

nrow <- nrow(var_comp) #Count number of parameters estimated

#Create summary table of variance components
out <- data.frame(replicate=rep)

#Rearrange variance component data so that it is on one line
for(i in 1:nrow){
  out2 <- data.frame(replicate=rep,Gamma=var_comp[i,1],Comp=var_comp[i,2], SE=var_comp[i,3], Z=var_comp[i,4],Con=var_comp[i,5])
  colnames(out2) <-  paste(c("",rep(rownames(var_comp)[i],5)),c("",rep("_",5)),colnames(out2), sep="")   #Add model term suffix (first 5 characters only) to each colname
  colnames(out2) <-  strsplit(colnames(out2), split="_Comp")  #remove "_Comp" suffix
  
  out <- merge(out, out2, by=c("replicate"), sort = TRUE)      
} 
rm(out2)

#get residual correlations, SEs and Z
  out$res_cor_est_1_2     <- out$res_cov_est_1_2 / (sqrt(out$res_var_est_1) * sqrt(out$res_var_est_2))
  if(out$res_cov_est_1_2_Con == "Fixed") {
    out$res_cor_est_1_2_SE <- 0
    out$res_cor_est_1_2_Z  <- 0
  } else {
    out$res_cor_est_1_2_SE  <- out$res_cor_est_1_2 / out$res_cov_est_1_2_Z
    out$res_cor_est_1_2_Z   <- out$res_cov_est_1_2_Z
  }  
  
  out$h2_est_1 <- h2_est_1$Estimate
  out$h2_est_2 <- h2_est_2$Estimate
  out$h2_est_1_se <- h2_est_1$SE
  out$h2_est_2_se <- h2_est_2$SE
  
  all_out <- rbind(all_out, out)
  rm(out)
} #end simulations

########################################################################################## 
#Basic graphical outputs
########################################################################################## 

hist(all_out$add_var_est_1)
abline(v=true_additive_var_1,col="red")

hist(all_out$add_var_est_2)
abline(v=true_additive_var_2,col="red")

hist(all_out$add_cor_est_1_2)
abline(v=true_additive_cor_1_2,col="red")

########################################################################################## 
#Output data to file
########################################################################################## 

write.table(all_out, file = "all_out.csv",append = FALSE, quote = FALSE, sep = ",",
            na = "NA", dec = ".", row.names = FALSE, col.names = TRUE) 
getwd()