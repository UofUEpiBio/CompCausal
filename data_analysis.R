###################
## load packages ##
###################
library(tidyverse)
library(splines)
library(gridExtra)
library(rms)

source(paste0(getwd(), "/HelperFunction.R"))
source(paste0(getwd(), "/singleindexmodelfunctions.R"))
source(paste0(getwd(), "/est_s_t_y.R"))
source(paste0(getwd(), "/est_exchange.R"))
source(paste0(getwd(), "/SIDR_Ravinew.R"))
source(paste0(getwd(), "/SIDRnew.R"))
source(paste0(getwd(), "/SensIAT_sim_outcome_modeler_mave.R"))

###############
## load data ##
###############

df <- readRDS(paste0(getwd(), "/TOIB3/TOIB_aggregated.rds"))
R <- df$isRCT
t <- df$isTropical
Y <- df$womac_12m
M <- ifelse(is.na(df$womac_12m), 0, 1)
n= length(t)

## X
X <- df %>% dplyr::select(c(age, womac_bq, expectationb, ChronicPainb))
## complete cases
complete_indx <- complete.cases(X)
X <- X[complete_indx, ]
R <- R[complete_indx]
t <- t[complete_indx]
Y <- Y[complete_indx]
M <- M[complete_indx]
n <- length(t)

data <- data.frame(cbind(Y, M, R, t))
colnames(data) <- c("Y", "M", "R", "t")

## relevel covariates
X$expectationb <- ifelse(X$expectationb %in% c(0, 1), "Much/A little worse",
                         ifelse(X$expectationb==2, "About the same", "A little/much better"))
X$expectationb <- relevel(factor(X$expectationb), ref="Much/A little worse")
X$ChronicPainb <- ifelse(X$ChronicPainb %in% c(0, 1), "1-2", "3-4")
X$ChronicPainb <- relevel(factor(X$ChronicPainb), ref="1-2")

################
## estimation ##
################

gamma1 <- seq(-2, 2, length=101)
gamma0 <- seq(-2, 2, length=101)

seed <- 1000096

topical_WOMAC_12m <- fit_one_analysis(trt_val=1, fold=fold, seed=rand_seed, gamma=seq(-2, 2, by=0.5), IF_output=TRUE,
                                      single_index_method="norm1coef", method="optim", kernel="dnorm")
oral_WOMAC_12m    <- fit_one_analysis(trt_val=0, fold=fold, seed=rand_seed, gamma=seq(-2, 2, by=0.5), IF_output=TRUE,
                                      single_index_method="norm1coef", method="optim", kernel="dnorm")

r_topical <- data.frame(est=topical_WOMAC_12m$est_trunc, var=topical_WOMAC_12m$var_trunc, 
                        est_R0=topical_WOMAC_12m$est_trunc_R0, var_R0=topical_WOMAC_12m$var_trunc_R0, 
                        est_R1=topical_WOMAC_12m$est_trunc_R1, var_R1=topical_WOMAC_12m$var_trunc_R1, 
                        t=t.topical, t_R0=t_R0.topical, t_R1=t_R1.topical)

r_oral <- data.frame(est=oral_WOMAC_12m$est_trunc, var=oral_WOMAC_12m$var_trunc, 
                     est_R0=oral_WOMAC_12m$est_trunc_R0, var_R0=oral_WOMAC_12m$var_trunc_R0, 
                     est_R1=oral_WOMAC_12m$est_trunc_R1, var_R1=oral_WOMAC_12m$var_trunc_R1, 
                     t=t.oral, t_R0=t_R0.oral, t_R1=t_R1.oral)

est_plot <- data.frame(est=c(r_topical$est, r_oral$est, 
                             r_topical$est_R0, r_oral$est_R0, 
                             r_topical$est_R1, r_oral$est_R1), 
                       treatment=c(rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0)), 
                                   rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0)), 
                                   rep("Topical", each=length(gamma1)), rep("Oral", each=length(gamma0))), 
                       lowerCI=c(r_topical$est-r_topical$t*sqrt(r_topical$var), 
                                 r_oral$est-r_oral$t*sqrt(r_oral$var), 
                                 r_topical$est_R0-r_topical$t_R0*sqrt(r_topical$var_R0), 
                                 r_oral$est_R0-r_oral$t_R0*sqrt(r_oral$var_R0), 
                                 r_topical$est_R1-r_topical$t_R1*sqrt(r_topical$var_R1), 
                                 r_oral$est_R1-r_oral$t_R1*sqrt(r_oral$var_R1)), 
                       upperCI=c(r_topical$est+r_topical$t*sqrt(r_topical$var), 
                                 r_oral$est+r_oral$t*sqrt(r_oral$var), 
                                 r_topical$est_R0+r_topical$t_R0*sqrt(r_topical$var_R0), 
                                 r_oral$est_R0+r_oral$t_R0*sqrt(r_oral$var_R0), 
                                 r_topical$est_R1+r_topical$t_R1*sqrt(r_topical$var_R1), 
                                 r_oral$est_R1+r_oral$t_R1*sqrt(r_oral$var_R1)), 
                       gamma=c(gamma1, gamma0, gamma1, gamma0, gamma1, gamma0), 
                       group=c(rep("CC", each=length(gamma1)+length(gamma0)), 
                               rep("PPS", each=length(gamma1)+length(gamma0)), 
                               rep("RCT", each=length(gamma1)+length(gamma0))))
est_plot$treatment <- factor(est_plot$treatment, levels=c("Topical", "Oral"))
est_plot$group <- factor(est_plot$group, levels=c("CC", "PPS", "RCT"))

## plot for estimates
plot1 <- ggplot(est_plot[which(est_plot$treatment=="Topical"), ], aes(x=gamma, y=est, color=group))+
  geom_line(linewidth = 1.2)+
  geom_line(aes(x=gamma, y=lowerCI), alpha=0.3, linewidth = 1.2, linetype = "dashed")+
  geom_line(aes(x=gamma, y=upperCI), alpha=0.3, linewidth = 1.2, linetype = "dashed")+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_color_manual(values=c("black", "blue", "green4"))+
  scale_x_continuous(breaks=seq(-2, 2, 1), limits = c(-2.05, 2.3)) +
  scale_y_continuous(breaks=seq(30, 55, 5), limits = c(30, 55)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))
plot1_f <-   plot1+
  annotate("text", x = 2, y = 40, label = as.character(expression(E*"["*Y(1)-Y[0]*"]")), parse = TRUE, size = 4)+
  annotate("text", x = 2, y = 43, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~R==0*"]")), parse = TRUE, size = 4, color="blue")+
  annotate("text", x = 2, y = 36.4, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~R==1*"]")), parse = TRUE, size = 4, color="green4")

plot2 <- plot1 + est_plot[which(est_plot$treatment=="Oral"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("text", x = 2, y = 41, label = as.character(expression(E*"["*Y(0)-Y[0]*"]")), parse = TRUE, size = 4)+
  annotate("text", x = 2, y = 44.8, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~R==0*"]")), parse = TRUE, size = 4, color="blue")+
  annotate("text", x = 2, y = 37, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~R==1*"]")), parse = TRUE, size = 4, color="green4")


png("est_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()

####################
## Counterfactual ##
####################

Counter_topical <- counterfactual(Y=Y, Y0=X$womac_bq, M=M, R=R, X=X, t=t, trt=1, 
                                  gamma=gamma1, est=topical_WOMAC_12m$est_trunc, 
                                  est_R0=topical_WOMAC_12m$est_trunc_R0, 
                                  est_R1=topical_WOMAC_12m$est_trunc_R1)

Counter_oral <- counterfactual(Y=Y, Y0=X$womac_bq, M=M, R=R, X=X, t=t, trt=0, 
                               gamma=gamma0, est=oral_WOMAC_12m$est_trunc, 
                               est_R0=oral_WOMAC_12m$est_trunc_R0, 
                               est_R1=oral_WOMAC_12m$est_trunc_R1)

Counter.plot.dat <- tibble(Data=c(Counter_topical$EY_t0, Counter_oral$EY_t0, Counter_topical$EY_t, Counter_oral$EY_t, 
                                  Counter_topical$EY_R0_t0, Counter_oral$EY_R0_t0, Counter_topical$EY_R0_t, Counter_oral$EY_R0_t, 
                                  Counter_topical$EY_R1_t0, Counter_oral$EY_R1_t0, Counter_topical$EY_R1_t, Counter_oral$EY_R1_t), 
                           Gamma=c(gamma1, gamma0, gamma1, gamma0, gamma1, gamma0, gamma1, gamma0, 
                                   gamma1, gamma0, gamma1, gamma0), 
                           Treatment=factor(c(rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0)), 
                                              rep("Topical", length(gamma1)), rep("Oral", length(gamma0))), 
                                            levels=c("Topical", "Oral")), 
                           Estimate=factor(c(rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0)), 
                                             rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0)), 
                                             rep("Induced", length(gamma1)), rep("Induced", length(gamma0)), 
                                             rep("Observed", length(gamma1)), rep("Observed", length(gamma0))), 
                                           levels=c("Observed", "Induced")), 
                           group=factor(c(rep("CC", each=2*length(gamma1)+2*length(gamma0)), 
                                          rep("PPS", each=2*length(gamma1)+2*length(gamma0)), 
                                          rep("RCT", each=2*length(gamma1)+2*length(gamma0))), 
                                        levels=c("CC", "PPS", "RCT")))

Counter.plot.data_split <- split(Counter.plot.dat, f=Counter.plot.dat$Treatment)

Counter.plot1 <- ggplot(data=Counter.plot.data_split$Topical, 
                        aes(x=Gamma, y=Data, color=group, linetype = Estimate)) +
  geom_line() +
  scale_color_manual(values=c("black", "blue", "green4"))+
  scale_x_continuous(breaks=seq(-2, 2, 1), limits = c(-2.1, 2.1)) +
  scale_y_continuous(breaks=seq(25, 55, 5), limits = c(25, 55)) +
  facet_wrap(~Treatment, scales = "fixed")+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months") +
  theme_bw()+
  scale_linetype_manual(values= c("dashed", "solid"), guide = "none")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 16, face="bold"),
        axis.title.x = element_text(size = 18, face="bold"),
        axis.text.y = element_text(size = 16, face="bold"),
        axis.title.y = element_text(size = 18, face="bold"),
        strip.text.x = element_text(size = 16, face="bold"))

Counter.plot1_f <- Counter.plot1+
  annotate("segment", x = -0.5, xend = 0, y = 34,  yend = 34, linetype = "dashed", color="black")+
  annotate("text",  x = 0.2, y = 34, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==1*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = -0.5, xend = 0, y = 32.5,  yend = 32.5, linetype = "solid", color="black")+
  annotate("text",  x = 0.2, y = 32.5, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==0*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = -0.5, xend = 0, y = 31,  yend = 31, linetype = "dashed", color="blue")+
  annotate("text",  x = 0.2, y = 31, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==1* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = -0.5, xend = 0, y = 29.5,  yend = 29.5, linetype = "solid", color="blue")+
  annotate("text",  x = 0.2, y = 29.5, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==0* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = -0.5, xend = 0, y = 28,  yend = 28, linetype = "dashed", color="green4")+
  annotate("text",  x = 0.2, y = 28, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==1* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")+
  annotate("segment", x = -0.5, xend = 0, y = 26.5,  yend = 26.5, linetype = "solid", color="green4")+
  annotate("text",  x = 0.2, y = 26.5, label = as.character(expression(E*"["*Y(1)-Y[0]~'|'~T==0* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")

Counter.plot2 <- Counter.plot1 + Counter.plot.data_split$Oral+
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x = -0.5, xend = 0, y = 34,  yend = 34, linetype = "dashed", color="black")+
  annotate("text",  x = 0.2, y = 34, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==0*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = -0.5, xend = 0, y = 32.5,  yend = 32.5, linetype = "solid", color="black")+
  annotate("text",  x = 0.2, y = 32.5, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==1*"]")), parse = TRUE, size = 3, hjust = 0)+
  annotate("segment", x = -0.5, xend = 0, y = 31,  yend = 31, linetype = "dashed", color="blue")+
  annotate("text",  x = 0.2, y = 31, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==0* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = -0.5, xend = 0, y = 29.5,  yend = 29.5, linetype = "solid", color="blue")+
  annotate("text",  x = 0.2, y = 29.5, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==1* "," ~R==0*"]")), parse = TRUE, size = 3, hjust = 0, color="blue")+
  annotate("segment", x = -0.5, xend = 0, y = 28,  yend = 28, linetype = "dashed", color="green4")+
  annotate("text",  x = 0.2, y = 28, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==0* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")+
  annotate("segment", x = -0.5, xend = 0, y = 26.5,  yend = 26.5, linetype = "solid", color="green4")+
  annotate("text",  x = 0.2, y = 26.5, label = as.character(expression(E*"["*Y(0)-Y[0]~'|'~T==1* "," ~R==1*"]")), parse = TRUE, size = 3, hjust = 0, color="green4")

png("Counter_Plot.png", width = 10, height = 7, units = 'in',res = 600)
grid.arrange(Counter.plot1_f,Counter.plot2, ncol=2)
dev.off()

#######################
## Treatment effects ##
#######################

diff <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_R0 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
diff_R1 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var_R0 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))
var_R1 <- matrix(0, nrow=length(gamma1), ncol=length(gamma0))

fold_idx <- split(seq_along(topical_WOMAC_12m$fold_index_pain),
                  topical_WOMAC_12m$fold_index_pain)
topical_IF_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc))
topical_IF_R0_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc_R0))
topical_IF_R1_mat <- t(do.call(rbind, topical_WOMAC_12m$IF_trunc_R1))
scale_n <- 1/n

for(g_0 in seq_along(gamma0)){
  
  diff[, g_0] <- topical_WOMAC_12m$est_trunc-oral_WOMAC_12m$est_trunc[g_0]
  diff_R0[, g_0] <- topical_WOMAC_12m$est_trunc_R0-oral_WOMAC_12m$est_trunc_R0[g_0]
  diff_R1[, g_0] <- topical_WOMAC_12m$est_trunc_R1-oral_WOMAC_12m$est_trunc_R1[g_0]
  
  IF_diff    <- topical_IF_mat - oral_WOMAC_12m$IF_trunc[[g_0]]
  IF_R0_diff <- topical_IF_R0_mat - oral_WOMAC_12m$IF_trunc_R0[[g_0]]
  IF_R1_diff <- topical_IF_R1_mat - oral_WOMAC_12m$IF_trunc_R1[[g_0]]
  
  var_temp    <- vapply(fold_idx, function(id) colVars(IF_diff[id, ]),    numeric(101))
  var_R0_temp <- vapply(fold_idx, function(id) colVars(IF_R0_diff[id, ]), numeric(101))
  var_R1_temp <- vapply(fold_idx, function(id) colVars(IF_R1_diff[id, ]), numeric(101))
  
  var[, g_0]    <- rowMeans(var_temp)    * scale_n
  var_R0[, g_0] <- rowMeans(var_R0_temp) * scale_n
  var_R1[, g_0] <- rowMeans(var_R1_temp) * scale_n
  
}

lowerCI <- diff-diff_t*sqrt(var)
lowerCI_R0 <- diff_R0-diff_t_R0*sqrt(var_R0)
lowerCI_R1 <- diff_R1-diff_t_R1*sqrt(var_R1)

upperCI <- diff+diff_t*sqrt(var)
upperCI_R0 <- diff_R0+diff_t_R0*sqrt(var_R0)
upperCI_R1 <- diff_R1+diff_t_R1*sqrt(var_R1)

## contour plot
png("contour_est_s_t_y.png", width = 10, height = 10, units = 'in',res = 300)
par(cex.lab = 2, cex.axis = 1.5, mgp = c(2.5, 1, 0))
contour(gamma1, gamma0, diff, nlevels = 15 , xlab = substitute(paste(bold(gamma[1]))), ylab = substitute(paste(bold(gamma[0]))), 
        lty="solid", lwd=1, col="black", labcex = 1.5, axes = F)
axis(1, at = seq(-5, 5, 1), labels = round(seq(-5, 5, 1), 1))
axis(2, at = seq(-5, 5, 1), labels = round(seq(-5, 5, 1), 1))
contour(gamma1, gamma0, diff, levels = 0, lwd = 3, add = T, col = "red", drawlabels = FALSE)
contour(gamma1, gamma0, lowerCI, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
contour(gamma1, gamma0, upperCI, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
#points(0, 0, pch = 3, col = "black", cex = 2, lwd = 2)
dev.off()

png("contour_est_R0_s_t_y.png", width = 10, height = 10, units = 'in',res = 300)
par(cex.lab = 2, cex.axis = 1.5, mgp = c(2.5, 1, 0))
contour(gamma1, gamma0, diff_R0, nlevels = 15 , xlab = substitute(paste(bold(gamma[1]))), ylab = substitute(paste(bold(gamma[0]))), 
        lty="solid", lwd=1, col="black", labcex = 1.5, axes = F)
axis(1, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
axis(2, at = seq(-2, 2, 1), labels = round(seq(-2, 2, 1), 1))
contour(gamma1, gamma0, diff_R0, levels = 0, lwd = 3, add = T, col = "red", drawlabels = FALSE)
contour(gamma1, gamma0, lowerCI_R0, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
contour(gamma1, gamma0, upperCI_R0, lwd=3,levels = 0, add = T, col = "blue", drawlabels = FALSE)
#points(0, 0, pch = 3, col = "black", cex = 2, lwd = 2)
dev.off()


####################
## Exchangability ##
####################

gamma1_extend <- seq(-2.5, 2.5, length=101)
gamma0_extend <- seq(-2.5, 2.5, length=101)

fit_one_analysis <- function(trt_val, fold, seed, gamma) {
  # try mave first
  out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                     t=t, trt=trt_val, gamma=gamma, fold=fold,
                     IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim",
                     single_index_method="norm1coef", use_mave=TRUE, seed=seed),
             silent = TRUE)
  
  if (inherits(out, "try-error")) {
    # fallback cumSIR
    out <- try(est_psi(Y=Y, M=M, R=R, X=X,
                       t=t, trt=trt_val, gamma=gamma, fold=fold,
                       IF_output=TRUE, simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim",
                       single_index_method="norm1coef", use_mave=FALSE, seed=seed),
               silent = TRUE)
  }
  
  if (inherits(out, "try-error") || is.null(out)) {
    return(NULL)
  }else{
    return(out)
  }
}

topical_WOMAC_12m_extend <- fit_one_analysis(trt_val=1, fold=5, seed=2000056, gamma=gamma1_extend)
oral_WOMAC_12m_extend <- fit_one_analysis(trt_val=0, fold=5, seed=2000056, gamma=gamma0_extend)


topical_WOMAC_12m_exchange <- est_psi_exchange(Y=Y, Y0, M=M, Y0=X$womac_bq, R=R, X=X, t=t, trt=1, 
                                               gamma=gamma1_extend, fold=5, seed=2000056, IF_output=TRUE, 
                                               simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                                               single_index_method="norm1coef", use_mave=TRUE)
oral_WOMAC_12m_exchange <- est_psi_exchange(Y=Y, M=M, Y0=X$womac_bq, R=R, X=X, t=t, trt=0, 
                                            gamma=gamma0_extend, fold=5, seed=2000056, IF_output=TRUE, 
                                            simple_trunc=FALSE, quant=NULL, kernel="dnorm", method="optim", 
                                            single_index_method="norm1coef", use_mave=TRUE)

est_plot_combine <- data.frame(est=c(topical_WOMAC_12m_extend$est_trunc, oral_WOMAC_12m_extend$est_trunc, 
                                     topical_WOMAC_12m_extend$est_trunc_R0, oral_WOMAC_12m_extend$est_trunc_R0), 
                               treatment=c(rep("Topical", each=length(gamma1_extend)), rep("Oral", each=length(gamma0_extend)), 
                                           rep("Topical", each=length(gamma1_extend)), rep("Oral", each=length(gamma0_extend))), 
                               gamma=c(gamma1_extend, gamma0_extend, gamma1_extend, gamma0_extend), 
                               group=c(rep("CC", each=length(gamma1_extend)+length(gamma0_extend)), rep("PPS", each=length(gamma1_extend)+length(gamma0_extend))))
est_plot_combine$treatment <- factor(est_plot_combine$treatment, levels=c("Topical", "Oral"))
est_plot_combine$group <- factor(est_plot_combine$group, levels=c("CC", "PPS"))

plot1 <- ggplot(est_plot_combine[which(est_plot_combine$treatment=="Topical"&est_plot_combine$group=="CC"), ], aes(x=gamma, y=est))+
  geom_line(linewidth = 1.2)+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_x_continuous(breaks=seq(-2.5, 2.5, 1), limits = c(-2.5, 2.5)) +
  scale_y_continuous(breaks=seq(33, 48, 5), limits = c(33, 48)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))

plot1_f <- plot1+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which(gamma1_extend==0)], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which(gamma1_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-0.5))], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-1))], 
           yend = topical_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma1_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[1]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[1]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[1]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

plot2 <- plot1 + est_plot_combine[which(est_plot_combine$treatment=="Oral"&est_plot_combine$group=="CC"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which(gamma0_extend==0)], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which(gamma0_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-0.5))], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-1))], 
           yend = oral_WOMAC_12m_exchange$est_trunc[which.min(abs(gamma0_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[0]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[0]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[0]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

png("combine_exchange_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()

plot1 <- ggplot(est_plot_combine[which(est_plot_combine$treatment=="Topical"&est_plot_combine$group=="PPS"), ], aes(x=gamma, y=est))+
  geom_line(linewidth = 1.2)+
  facet_grid(.~`treatment`, scales = "fixed")+
  scale_x_continuous(breaks=seq(-2.5, 2.5, 1), limits = c(-2.5, 2.5)) +
  scale_y_continuous(breaks=seq(33, 48, 5), limits = c(33, 48)) +
  theme_bw()+
  labs(x=expression(gamma[1]), y = "WOMAC pain scores at 12 months")+
  theme(legend.position = "none",
        axis.text.x = element_text(size = 18, face="bold"),
        axis.title.x = element_text(size = 20, face="bold"),
        axis.text.y = element_text(size = 18, face="bold"),
        axis.title.y = element_text(size = 20, face="bold"),
        strip.text.x = element_text(size = 18, face="bold"))

plot1_f <- plot1+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which(gamma1_extend==0)], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which(gamma1_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-0.5))], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-1))], 
           yend = topical_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma1_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[1]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[1]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[1]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

plot2 <- plot1 + est_plot_combine[which(est_plot_combine$treatment=="Oral"&est_plot_combine$group=="PPS"), ]+ 
  labs(x=expression(gamma[0]), y=NULL)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which(gamma0_extend==0)], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which(gamma0_extend==0)], linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-0.5))], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-0.5))], linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("segment", x=-2.5, xend=2.5,  
           y = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-1))], 
           yend = oral_WOMAC_12m_exchange$est_trunc_R0[which.min(abs(gamma0_extend-1))], linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("segment", x = 0.5, xend = 1, y = 36,  yend = 36, linetype = "dashed", color="cyan", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 36, label = as.character(expression(gamma[0]*"'"==0)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 35,  yend = 35, linetype = "dashed", color="coral", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 35, label = as.character(expression(gamma[0]*"'"==0.5)), parse = TRUE, size = 6, hjust = 0)+
  annotate("segment", x = 0.5, xend = 1, y = 34,  yend = 34, linetype = "dashed", color="blueviolet", linewidth = 1.1)+
  annotate("text",  x = 1.2, y = 34, label = as.character(expression(gamma[0]*"'"==1)), parse = TRUE, size = 6, hjust = 0)

png("combine_exchange_R0_s_t_y.png", width = 15, height = 10, units = 'in',res = 300)
grid.arrange(plot1_f,plot2,ncol=2)
dev.off()





