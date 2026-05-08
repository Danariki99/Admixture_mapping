args = commandArgs(trailingOnly=TRUE)
library("data.table")

o<-fread(args[1])
i<-fread(args[2])

q<-merge(o, i, by="ID")

q$P_ratio<-q$P.y/q$P.x

q$P_delta<-log10(q$P.y)-log10(q$P.x)

write.table(q, args[3], row.names=FALSE, quote=FALSE)
