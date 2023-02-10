setwd("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")

load("../output/simOP110.01.01.Rdata")

View(rs[[1]])
a <- rs[[1]][rs[[1]]>300]-300
b <- rs[[1]][rs[[1]]<301]
mean(a)
sd(a)
length(b)
