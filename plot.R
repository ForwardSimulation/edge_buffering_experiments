library(ggplot2)
library(readr)

args = commandArgs(trailing=TRUE)

data = read_delim(args[1], delim=" ")

p = ggplot(aes(x=N,y=time),data=data) + 
    geom_point(aes(col=method,shape=as.factor(tsimplify))) +
    xlab("Diploid pop size") +
    ylab("Run time (seconds)")
ggsave("timings.pdf")
                                                   
