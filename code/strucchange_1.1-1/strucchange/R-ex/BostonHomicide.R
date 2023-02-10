### Name: BostonHomicide
### Title: Youth Homicides in Boston
### Aliases: BostonHomicide
### Keywords: datasets

### ** Examples

data(BostonHomicide)

fm <- glm(homicides ~ population + season, data = BostonHomicide,
          family = poisson)
anova(fm, test = "F")



