#### trying some modeling approaches
setwd("/Users/jack/Desktop/Research/ICB_Phys/physiology/physiology_space")

rm(list=ls())

load(file = "data/full_data.RData")

# primary responsiveness model

library(car)
library(lme4)
library(lmerTest)

library(evolvability)
library(ape)
library(car)

tree <- read.tree("data/resp_tree.tre") #Portik et al tree previously trimmed and names replaced with genera
plot(tree)                              #check to make sure tree makes sense
A <- Matrix::Matrix(ape::vcv(tree), sparse = TRUE) # create matrix of relatedness
colnames(A) <- rownames(A) <- tree$tip.label       # fix row and column names

phylo_mod <- evolvability::Almer(response_numeric ~ habitat_factor*lung_factor*treatment*temp_difference +
                     (1|Run) + (1|Code) + (1|genus) + (1|population)
                   , data = full_data, A = list(genus = A))
summary(phylo_mod)
car::Anova(phylo_mod, type = "3") #same coefficients significant

phylo_mod2 <- update(phylo_mod, ~.-habitat_factor:lung_factor:treatment:temp_difference)
car::Anova(phylo_mod2, type = "3")
#test if four-way interaction should be included or not
anova(phylo_mod, phylo_mod2)

#check next group of models (without each non-significant 3-way interaction and without both)
phylo_mod3 <- update(phylo_mod2, ~.-habitat_factor:lung_factor:temp_difference)
phylo_mod4 <- update(phylo_mod2, ~.-habitat_factor:lung_factor:treatment)
phylo_mod5 <- update(phylo_mod3, ~.-habitat_factor:lung_factor:treatment)

#test if adding factors supported
anova(phylo_mod, phylo_mod2, phylo_mod3, phylo_mod4, phylo_mod5)

#most parsimonious wins, all remaining factors are significant or involved in higher order interactions
summary(phylo_mod5) 
car::Anova(phylo_mod5, type = "3")

rm(list = c("phylo_mod", "phylo_mod2", "phylo_mod3", "phylo_mod4"))

#post hoc tests
library(emmeans)

# estimated marginal means of responsiveness at temp_diff = 0 for treatment:lung groups

emm_lung <- emmeans(phylo_mod5, specs = list( ~ treatment:lung_factor:temp_difference),
                    at = list(temp_difference = 0))
#summary for custom contrasts:
lung_em <- summary(contrast(emm_lung, "consec", simple = "each", adjust = "bonferroni"))
lung_em$`simple contrasts for treatment`   #test treatment effects within lunged and lungless
lung_em$`simple contrasts for lung_factor` #test lung effects within high and low treatments

# custom contrast to be compared 
lung_treatment <- contrast(emm_lung, list("lunged_low_high" = c(-1,1,0,0), "lungless_low_high" = c(0,0,-1,1)), adjust = "bonferroni")
pairs(lung_treatment) #compare differences at high and low for lunged and lungless

#same process for habitat comparisons:
emm_habitat <- emmeans(phylo_mod5, specs = list( ~ treatment:habitat_factor:temp_difference),
                       at = list(temp_difference = 0))

habitat_em <- summary(contrast(emm_habitat, "consec", simple = "each", adjust = "bonferroni"))
habitat_em$`simple contrasts for treatment`   #test treatment effects within streams and ponds
habitat_em$`simple contrasts for habitat_factor` #test habitat effects within high and low treatments

# custom contrast to be compared 
habitat_treatment <- contrast(emm_habitat, list("habitated_low_high" = c(-1,1,0,0), "habitatless_low_high" = c(0,0,-1,1)), adjust = "bonferroni")
pairs(habitat_treatment) #compare differences at high and low for habitated and habitatless

# compare trends (slopes) as well, not just point comparisons at a given temp_diff
#lungs:treatment
emt_lung<- emtrends(phylo_mod5,  ~ treatment:lung_factor, var = "temp_difference")
lung_emt <- summary(contrast(emt_lung, "consec", simple = "each", adjust = "bonferroni"))
lung_emt$`simple contrasts for treatment`
lung_emt$`simple contrasts for lung_factor`
lung_treat_slope <- contrast(emt_lung, list("lunged_low_high" = c(-1,1,0,0), "lungless_low_high" = c(0,0,-1,1)), adjust = "bonferroni")
pairs(lung_treat_slope)


emt_hab <- emtrends(phylo_mod5,  ~ treatment:habitat_factor, var = "temp_difference")
hab_emt <- summary(contrast(emt_hab, "consec", simple = "each", adjust = "bonferroni"))
hab_emt$`simple contrasts for treatment`
hab_emt$`simple contrasts for habitat_factor`

hab_treat_slope <- contrast(emt_hab, list("pond_low_high" = c(-1,1,0,0), "stream_low_high" = c(0,0,-1,1)), adjust = "bonferroni")
pairs(hab_treat_slope)


######################################################################################################

rm(list=ls())

load(file = "data/lung_matrix.RData")

#remove two genera with very few observations/individuals for genus-level analyses
lung_matrix2 <- lung_matrix[which(lung_matrix$genus != "Phrynomantis"),]
lung_matrix2 <- lung_matrix2[which(lung_matrix2$genus != "Strongylopus"),]


glmer_mod <- glmer(sums ~ habitat*treatment + offset(log(num_periods)) + (1|run) +
                      (1|genus) + (1|population)
                    , data = lung_matrix, family = poisson(link = "log"), 
                    glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(glmer_mod)
car::Anova(glmer_mod, type = "III")

em <- emmeans(glmer_mod, ~ treatment:habitat)
contrast(em, "consec", simple = "each", adjust = "bonferroni")


glmer_mod_genus <- glmer(sums ~ genus*treatment + offset(log(num_periods)) + (1|run) +
                     (1|population)
                   , data = lung_matrix2, family = poisson(link = "log"), 
                   glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(glmer_mod_genus)
car::Anova(glmer_mod_genus, type = "III")

em2 <- emmeans(glmer_mod_genus, ~ genus*treatment)

a <- summary(contrast(em2, "consec", simple = "each", adjust = "bonferroni"))
b <- a$`simple contrasts for treatment`; rm(a)

save(b, file = "data/emm_genus.RData") #save emms for Fig. 3


###############################################################################
#other models

#testing whether climbing behaviors differ across treatments and taxa

rm(list=ls())

load(file = "data/full_data.RData")

#only keep relevant rows with climbing taxa:
heleodata <- full_data[which(full_data$Genus == "Heleophryne" | full_data$Genus == "Hadromophryne"),]
heleodata <- heleodata[which(heleodata$response != "dead"),] #remove non-living observations
rm(full_data) #cleanup

heleodata$OOW <- as.factor(heleodata$Out_of_water)

#full model with all interactions
heleo_mod_full <- lme4::glmer(Out_of_water ~ treatment*lung_factor +
                          (1|population) + (1|Code),
                        data = heleodata,
                        family = binomial())                                                                          
summary(heleo_mod_full)
car::Anova(heleo_mod_full)

#simplified model with no interaction
heleo_mod_simp <- glmer(Out_of_water ~ treatment+lung_factor +
                          (1|population) + (1|Code),
                        data = heleodata,
                        family = binomial())      

summary(heleo_mod_simp)
car::Anova(heleo_mod_simp)

#model comparison
anova(heleo_mod_full, heleo_mod_simp)

##############################################

#test whether surfacing behaviors differ between treatments in Schismaderma

rm(list=ls())

schis_data <- read.csv("data/schismaderma_vids.csv") #raw data 

#reformat data for modeling
schis_data_stack <- as.data.frame(matrix(ncol = 10, 
                                  nrow = nrow(schis_data)*20))

colnames(schis_data_stack) <- c(colnames(schis_data)[1:7], "vid_number", "time", "surface")

schis_data_stack$vid_number <- rep(c(1,2,3,4), each = nrow(schis_data)*5)
schis_data_stack$time <- rep(c(1,2,3,4,5), times = nrow(schis_data)*4)

for(i in 1:7){
  schis_data_stack[,i] <- rep(rep(schis_data[,i], each = 5), times = 4)
}

for(i in 1:dim(schis_data_stack)[1]){
  schis_data_stack$surface[i] <- schis_data[,match(paste("X",schis_data_stack$vid_number[i],"_",
                                            schis_data_stack$time[i], sep = ""), 
                                            colnames(schis_data))][which(schis_data$Code== schis_data_stack$Code[i])]}


mod_null <-  glmer(as.factor(surface) ~ 1 + (1|Code),
                   data = schis_data_stack,
                   family = binomial()) 
summary(mod_null)

mod_full <-  glmer(as.factor(surface) ~ 1 + treatment + (1|Code),
                   data = schis_data_stack,
                   family = binomial()) 
summary(mod_full)

#model comparison
anova(mod_full, mod_null)









