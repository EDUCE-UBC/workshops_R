##############################################################
# load packages
##############################################################
# Suite of packages for data manipulation and visualization
library(tidyverse)
# puts output of statistical models in a nice data frame
library(broom) 
# Split, process, and recombine data
library(plyr)
# Fit linear and generalized linear mixed-effects models
library(lme4)
# various functions, including Anova
library(car) 
# least-square means
library(lsmeans)
## generalized linear model fitter
## Also has quine data  set
library(MASS)

# Data set libraries
## Fruitfly longevity, size, and sexual activity
library(faraway)
## Life expectancy, GDP and population by country
library(gapminder)  
## A variety of data sets; we will use the plasma data
library(HSAUR3) 

##############################################################
# Experimental design
##############################################################
#### Exercise: Experimental design

#1. Discuss with following in pairs
    #- What are some advantages and disadvantages  of using a balanced experimental design? 
    #- Give an example of when a balanced design might not be possible.
    #-  There are 3 undergraduates assisting you with your experiment that assess the addiction potential of Saturday morning cartoons in rats. You need to run the experiments every Saturday, but one of your undergraduate assistants can only help out 2 Saturdays a month, while the other two undergraduate assistants can be there every Saturday. Rat behaviour is sensitive to handler. What should you do?
  
#2. *True or false.* A completely randomized design offers no control for lurking variables (a variable that is not included as an explanatory or response variable in the analysis).

##############################################################
# 1-way ANOVA with 2 groups
##############################################################
###############################
# explore fruit fly data
###############################
# Read in data from faraway package
data("fruitfly")

# Create categorical groups and subset data
fruitfly_2groups <- fruitfly %>%
  # Convert factors to character variables
  mutate_if(is.factor, as.character) %>% 
  # Create 2-level activity variable
  mutate(activity = ifelse(activity %in% c("isolated", "one", "many"), "no", "yes")) %>%
  # Create 2-level size variable (for later ANOVA)
  mutate(thorax = ifelse(thorax <= 0.8, "short", "long")) %>% 
  # Subset to equal group sizes for activity and size
  group_by(activity, thorax) %>% 
  sample_n(20)

# get number of rows and columns
dim(fruitfly_2groups)
# view first 6 records of data
head(fruitfly_2groups)
# see summary of data
summary(fruitfly_2groups)

# plot raw data points for each group as a transparent grey/black point
# overlay mean as a red diamond
ggplot(fruitfly_2groups, aes(x = activity, y = longevity)) +
  geom_jitter(position = position_jitter(0.15),
              alpha = 0.7) +
  stat_summary(fun.y = mean,
               geom = "point",
               shape = 18,
               size = 4,
               color="red") +
  xlab("Sexually Active") +
  ylab("Longevity (days)")

###############################
# ANOVA in R: 2 groups
###############################
# create an ANOVA "model" object
fruitfly_2groups_model <- aov(longevity ~ activity,
                              data = fruitfly_2groups)

# view output of aov() as a nice dataframe using tidy() from the broom package
tidy(fruitfly_2groups_model)

### EXERCISE: 1-way ANOVA
# 1. Using ANOVA, test if fruit fly longevity is effected by size (as measured by thorax length). What are your null and alternate hypotheses? What can you conclude from these results?

##############################################################
# 1-way ANOVA with 2 groups
##############################################################
###############################
# Re-explore the data
###############################
# Create categorical groups and subset data
fruitfly_3groups <- fruitfly %>%
  # Convert factors to character variables
  mutate_if(is.factor, as.character) %>% 
  # Create 3-level activity variable
  mutate(activity = ifelse(activity %in% c("isolated", "one", "many"), "none", activity)) %>%
  # Subset to equal group sizes for activity
  group_by(activity) %>% 
  sample_n(25)

# get number of rows and columns
dim(fruitfly_3groups)
# view first 6 records of data
head(fruitfly_3groups)
# see summary of data
summary(fruitfly_3groups)

# re-order factors to make them show up how we would like them on the plot
# instead of alphabetically (default R behaviour)
fruitfly_3groups$activity <- factor(fruitfly_3groups$activity,
                                    levels = c("none","low","high"))

# plot raw data points for each group as a transparent grey/black point
# overlay mean as a red diamond
ggplot(fruitfly_3groups, aes(x = activity, y = longevity)) +
  geom_jitter(position = position_jitter(0.15),
              alpha = 0.7) +
  stat_summary(fun.y = mean,
               geom = "point",
               shape = 18,
               size = 4,
               color = "red") +
  xlab("Sexual Activity") +
  ylab("Longevity (days)")

###############################
# ANOVA in R: > 2 groups
###############################
# create an ANOVA "model" object
fruitfly_3groups_model <- aov(longevity ~ activity,
                              data = fruitfly_3groups)

# view output of aov() as a nice dataframe using tidy() from the broom package
tidy(fruitfly_3groups_model)

###############################
# Assess which groups differ
###############################
# pairwise t-tests to observe group differences
tidy(pairwise.t.test(fruitfly_3groups$longevity,
                     fruitfly_3groups$activity,
                     p.adjust.method = "bonferroni",
                     pool.sd = TRUE,
                     paired = FALSE))

# Tukey's HSD test to observe group differences
tidy(TukeyHSD(fruitfly_3groups_model, "activity"))

##############################################################
# 2-way ANOVA with 2 groups
##############################################################
###############################
# Re-explore the data
###############################
# get number of rows and columns
dim(fruitfly_2groups)
# view first 6 records of data
head(fruitfly_2groups)
# see summary of data
summary(fruitfly_2groups)

# re-order factors to make them show up how we would like them on the plot
# instead of alphabetically (default R behaviour)
fruitfly_2groups$thorax <- factor(fruitfly_2groups$thorax,
                                  levels = c("short","long"))
# plot strip charts of longevity, grouped by sexual activity
# and colored by thorax length
ggplot(fruitfly_2groups,
       aes(x = activity, y = longevity, color = thorax)) +
  stat_summary(fun.y = mean,
               geom = "point",
               shape = 5,
               size = 4,
               position = position_dodge(0.5)) +
  geom_jitter(position = position_dodge(0.5), alpha = 0.3) +
  scale_color_manual(values=c("black", "dodgerblue3")) +
  xlab("Sexual Activity") +
  ylab("Longevity (days)")

###############################
# ANOVA in R: 2 variables with 2 groups
###############################
# create an ANOVA "model" object
fruitfly_2var_model <- aov(longevity ~ activity + thorax,
                           data = fruitfly_2groups)

# view output of aov() as a nice dataframe using tidy() from the broom package
tidy(fruitfly_2var_model)

###############################
# 2-way ANOVA with 2 groups including an interaction term
###############################
# plot to investigate possible interaction effect of sexual
# activity and thorax length on longevity
ggplot(fruitfly_2groups,
       aes(x = activity, y = longevity, color = thorax)) +
  stat_summary(fun.y = mean,
               geom = "point",
               shape = 18,
               size = 3) +
  stat_summary(fun.y = mean,
               geom = "line",
               aes(group = thorax)) +
  scale_color_manual(values=c("black", "dodgerblue3")) +
  xlab("Sexual Activity") +
  ylab("Longevity (days)")

# create an ANOVA "model" object
fruitfly_2var_model2 <- aov(longevity ~ activity * thorax,
                            data = fruitfly_2groups)

# view output of aov() as a nice dataframe using tidy() from the broom package
tidy(fruitfly_2var_model2)

#### Exercise: ANOVA
#Determine whether the following statements are *true or false*?
  
#1. ANOVA tests the null hypothesis that the sample means are	all	equal?
  
#2. We use ANOVA to compare the variances of the population?
  
#3. A one-way ANOVA is equivalent to a *t*-test when there are 2 groups to be compared.	

#4. In rejecting the null hypothesis, one can conclude that all the population means are different from one another?
  
###*Questions courtesy of Dr. Gabriela Cohen Freue's DSCI 562 course (UBC)*

##############################################################
# Linear regression
##############################################################
###############################
# Explore the gapminder data
###############################
#View the data frame
gapminder

#Plot life expectancy by year
gapminder %>%
  ggplot(aes(x = year, y = lifeExp)) +
  geom_point() +
  labs(y="Life expectancy (yrs)")

#Plot life expectancy by per-capita GCP
gapminder %>%
  ggplot(aes(x = gdpPercap, y = lifeExp)) +
  geom_point() +
  labs(y="Life expectancy (yrs)", x="Per-capita GDP")

###############################
# Linear models
###############################
gapminder %>%
  ggplot(aes(x = year, y = lifeExp)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 0.033, colour = "green") +
  geom_abline(intercept = -575, slope = 0.32, colour = "purple") +
  geom_hline(aes(yintercept = mean(lifeExp)), colour = "red") +
  geom_vline(aes(xintercept = mean(year)), colour = "blue") +
  labs(y="Life expectancy (yrs)")

### Exercise: Best fit lines
### At your table, discuss the following questions.

#1. Which line best describes the data?

#2. The red one is a horizontal line at the overall mean life expectancy. It seems a reasonable model, but what is missing?
  
###############################
# Simple linear regression
###############################
# create a lm "model" object
lifeExp_model1 <- lm(lifeExp ~ year, 
                     data = gapminder)

# view output of lm() as using summary()
summary(lifeExp_model1)

# Plot this fit
gapminder %>%
  ggplot(aes(x = year, y = lifeExp)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, colour = "green") +
  labs(y="Life expectancy (yrs)")

# Assess the model: residuals
# Set plot frame to 2 by 2
par(mfrow=c(2,2))
# Create diagnostic plots
plot(lifeExp_model1)

### Exercise: Linear models
### At your table, discuss the following questions.

#1. Looking at the summary plots above, do you feel that our model can be extrapolated to a much wider year range? Why or why not?
  
#2. Fit a linear model of life expectancy as a function of per-capita gdp. Using the summary table and diagnostic plots, discuss whether or not you think this is a good fit for these data.

###############################
# CAUTION!
###############################
#Anscombe data sets 
absc <- with(datasets::anscombe,
             tibble(X=c(x1,x2,x3,x4),
                    Y=c(y1,y2,y3,y4),
                    anscombe_quartet=gl(4,nrow(datasets::anscombe))
             )
)

#the same LM
absc %>%
  ggplot(aes(x=X,y=Y,group=anscombe_quartet)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_continuous(limits = c(0,20)) +
  facet_wrap(~anscombe_quartet)

###############################
# Multiple linear regression
###############################
# create a lm "model" object
lifeExp_model2 <- lm(lifeExp ~ year + gdpPercap, 
                     data = gapminder)

# view output of lm() as using summary()
summary(lifeExp_model2)

#  Assess model
par(mfrow=c(2,2))
plot(lifeExp_model2)

###############################
# Transforming predictors
###############################
# Sine example
gapminder %>%
  ggplot(aes(x = sin(gdpPercap), y = lifeExp)) +
  geom_point()

### Exercise: Tranforming predictors

#1. Find a function that makes the plot more “linear” and fit a model of life expectancy as a function of the transformed per-capita gdp. Is it a better model?
###  Go back to your original gdpPercap vs. lifeExp plot and think about what function creates a similar trend.

###############################
# log transformation
###############################
lifeExp_model3a <- lm(lifeExp ~ year + log(gdpPercap),
                      data = gapminder)

# equivalent to
gapminder <- gapminder %>%
  mutate(log_gdp = log(gdpPercap))

lifeExp_model3b <- lm(lifeExp ~ year + log_gdp,
                      data = gapminder)
# view results
summary(lifeExp_model3a)
summary(lifeExp_model3b)

###############################
# polynomial transformation
###############################
lifeExp_model4 <- lm(lifeExp ~ year + poly(gdpPercap),
                     data = gapminder)
summary(lifeExp_model4)

#or
lifeExp_model5 <- lm(lifeExp ~ year + gdpPercap + I(gdpPercap^2),
                     data = gapminder)
summary(lifeExp_model5)

###Exercise: Multiple linear regression

#1. So far, we have worked with lifeExp as our independent variable. Now, in small groups, try to produce a model of population (pop) using one or more of the variables available in gapminder.

##############################################################
# interactions and ANCOVA
##############################################################
#Add continent to the gapminder plot
gapminder %>%
  ggplot(aes(x = year, y = lifeExp, colour = continent)) +
  geom_point() +
  geom_smooth(method = "lm", se= FALSE) +
  labs(y="Life expectancy (yrs)")

#add an interaction term
lifeExp_model6 <- lm(lifeExp ~ year*continent,
                     data = gapminder)
summary(aov(lifeExp_model6))

# subet the data and fit again
gapminder %>%
  filter(continent %in% c("Oceania","Europe")) %>%
  lm(lifeExp ~ year*continent, data = .) %>% 
  aov() %>% 
  summary()

##############################################################
# LME
##############################################################
###############################
# explore esoph data
###############################
p <- ggplot(esoph, aes(ncontrols, ncases, group=agegp, colour=agegp)) +
  geom_jitter(height=0.25) +
  scale_colour_discrete("Age Group") +
  ylab("Number of Cases") + xlab("Number of Controls")
p

# Many LMs
p + geom_smooth(method="lm", se=FALSE, size=0.5)

# but small sample size
esoph %>%
  group_by(agegp) %>%
  dplyr::summarise(n=length(ncases)) %>%
  as.data.frame

###############################
# fitting LME
###############################
esoph_model <- lmer(ncases ~ ncontrols + (ncontrols | agegp),
                    data=esoph)
summary(esoph_model)

#extract line slopes and intercepts
tidy(esoph_model, "ran_modes")
#or
coef(esoph_model)

## Plot
ggplot(esoph, aes(ncontrols, ncases, group=agegp, colour=agegp)) +
  geom_jitter(height=0.25) +
  geom_abline(aes(intercept=intercept, slope=slope, colour=agegp)) +
  scale_colour_discrete("Age Group") +
  ylab("Number of Cases") + xlab("Number of Controls")

#### Exercise: LME

#1. Using the `sleepstudy` dataset, fit an LME on Reaction against Days, grouped by Subject.

#2. What is the intercept and slope of subject #310 in the model from question 1?

#3. CHALLENGE. Using the Teams dataset from the Lahman package, fit a model on runs (`R`) from the variables 'walks' (`BB`) and 'Hits' (`H`), grouped by team (`teamID`).
    #- *Hint*: wrap the scale function around each predictor variable.

##############################################################
# GLM
##############################################################
###############################
# logistic (binomial)
###############################
# Load and format data
ucb <- as.data.frame(UCBAdmissions) %>% 
  dplyr::rename(sex=Gender)

# Fit GLM binomial
ucb_model <- glm(Admit ~ sex * Dept, 
                 data = ucb, 
                 family = binomial, 
                 weights = Freq)

summary(ucb_model)

#check significance of covariates
Anova(ucb_model)

#can you drop a covariate?
drop1(ucb_model, test = "Chisq")

#fit success probabilities
ucb_model_sum <- lsmeans(ucb_model, ~ sex + Dept, type = "response")
ucb_model_sum
    # grouped by  department
      summary(ucb_model_sum, by = "Dept")
    # odds ratio
      contrast(ucb_model_sum, "pairwise", by = "Dept")

###Exercise: Logistic GLM

#1.In the plasma data (from the HSAUR3 package), use logistic regression to estimate the probabilities of ESR > 20, given the level of fibrinogen in the blood.
      
#2. Using the womensrole data set from the HSAUR3 package, try to fit a logistic regression to the agreement with the statement, given the years of education and the respondent’s sex (also attributed as gender in these data).
      
###############################
# count data (poisson)
###############################
#assuming variance equal
polyps_model1 <- glm(number ~ treat + age,
                           data = polyps,
                           family = poisson)
summary(polyps_model1)

#not assuming variance equal
polyps_model2 <- glm(number ~ treat + age, 
                     data = polyps, 
                     family = quasipoisson)
summary(polyps_model2)

# investigate differnce between 2 levels
polyps_model2_sum <- lsmeans(polyps_model2, ~ treat, 
                             type = "response")

contrast(polyps_model2_sum, "pairwise")

###Exercise: Poisson GLM

#1. Check which covariates have a significant effect on the response in the model fitted with the Poisson family and with the quasi-Poisson family and compare the results. What do you observe?

###############################
# Negative binomial
###############################
quine_model <- glm.nb(Days ~ Sex * (Age + Eth * Lrn), 
                      data = quine)
## equivalent to
## quine_model <- glm.nb(Days ~ Sex * Age + Sex * Eth * Lrn, data = quine)
summary(quine_model)

#check significance
Anova(quine_model)

#can we drop a covariate?
drop1(quine_model, test = "Chisq")

#difficult to interpret relationships or compare groups directly
##using 3 terms...
quine_model_sum1 <- lsmeans(quine_model, ~ Sex + Eth + Lrn, type = "response")

summary(quine_model_sum1, by = c("Sex", "Eth"))
summary(quine_model_sum1, by = c("Eth", "Lrn"))

## or even using 2 terms
quine_model_sum2 <- lsmeans(quine_model, ~ Sex + Age, type = "response")

summary(quine_model_sum2, by = "Sex")
summary(quine_model_sum2, by = "Age")

###Exercise: Quasi-poisson vs. negative binomial GLM

#1. Fit the above model with a Quasi-Poisson family and check the over-dispersion in that fit. Is there a difference in the significance of any terms compared to the NB model? Would a Poisson model be appropriate as well?