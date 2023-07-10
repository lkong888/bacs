a=read.delim(file("stdin"),skip=1)
a[,7] =a[,7]/a$Length *10**3/sum(a[,7]/a$Length *10**3) *10**6
write.table(a,stdout(),sep='\t',quote=F,row.names=F)

#rpk=a[,7]/a$Length *10**3
#TPM=rpk/sum(rpk) *10**6