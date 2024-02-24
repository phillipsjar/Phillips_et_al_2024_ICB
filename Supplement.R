# supplemental analyses


#test an ordinal model of responsiveness
rm(list=ls())

load(file = "data/full_data.RData")

#first do our best to fit expected model type: ordinal data
library(ordinal)
mod_clmm <- clmm(response_ordered ~ habitat_factor*lung_factor*treatment*
                   temp_difference + (1|Run) + (1|Code) + (1|genus) + 
                   (1|population) , data = full_data)
summary(mod_clmm)

mod_clmm2 <- update(mod_clmm, ~.-habitat_factor:lung_factor:treatment:temp_difference)
summary(mod_clmm2)

mod_clmm3 <- update(mod_clmm2, ~.-habitat_factor:lung_factor:temp_difference)
mod_clmm4 <- update(mod_clmm2, ~.-habitat_factor:lung_factor:treatment)
mod_clmm5 <- update(mod_clmm3, ~.-habitat_factor:lung_factor:treatment)

anova(mod_clmm, mod_clmm2, mod_clmm3, mod_clmm4, mod_clmm5)
summary(mod_clmm5)
summary(phylo_mod5)
a <- coef(mod_clmm5)[1:4]

# very evenly spaced, justifying a non-ordinal approach

rm(list=ls())
load(file = "data/lung_matrix.RData")

lung_matrix2 <- lung_matrix[which(lung_matrix$genus != "Phrynomantis"),]
lung_matrix2 <- lung_matrix2[which(lung_matrix2$genus != "Strongylopus"),]

library(phyr)
library(ape)
tree <- read.tree("data/vid_tree.tre")

phylo_mod <- phyr::communityPGLMM(sums ~ habitat*treatment + offset(log(num_periods)) + (1|run) +
                                     (1|genus__) + (1 | genus__@population), cov_ranef = list(genus = tree), 
                                   data = lung_matrix2, family = "poisson", add.obs.re = FALSE)
glmer_mod <- glmer(sums ~ habitat*treatment + offset(log(num_periods)) + (1|run) +
                     (1|genus) + (1|population)
                   , data = lung_matrix, family = poisson(link = "log"), 
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))


summary(phylo_mod)
summary(glmer_mod)


phylo_mod_genus <- phyr::communityPGLMM(sums ~ genus*treatment + offset(log(num_periods)) + (1|run) +
                                    (1|genus__) + (1 | genus__@population), cov_ranef = list(genus = tree), 
                                  data = lung_matrix2, family = "poisson", add.obs.re = FALSE)

glmer_mod_genus <- glmer(sums ~ genus*treatment + offset(log(num_periods)) + (1|run) +
                           (1|population)
                         , data = lung_matrix2, family = poisson(link = "log"), 
                         glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))

summary(phylo_mod_genus)
summary(glmer_mod_genus)

Anova(glmer_mod_genus, type = "3")

#phylogenetic models very similar

















