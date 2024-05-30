
args = commandArgs(trailingOnly=TRUE)
#print(args[1])
#print(args[2])
tt <- readRDS(args[1])
write.table(tt,file=args[2],sep=",", row.names=FALSE, quote=FALSE)
#print(tt)