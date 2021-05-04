first2<-read.table('/projectnb/bf528/users/swiss_cheese2/project_4/data/freq4.txt',row.names=NULL)
second2<-read.table('/projectnb/bf528/users/swiss_cheese2/project_4/data/freq5.txt',row.names=NULL)
third2<-read.table('/projectnb/bf528/users/swiss_cheese2/project_4/data/freq6.txt',row.names=NULL)

plot(ecdf(first2$row.names))

plot(ecdf(second2$row.names))


plot(ecdf(third2$row.names))

plot(ecdf(first2$row.names),
     col="green",xlab="Barcode Count",
ylab="Cumulative Proportion",
main="Distribution of Barcode Counts")
lines(ecdf(second2$row.names),
      col="blue")
lines(ecdf(third2$row.names),
      col="red")


plot(ecdf(first2$row.names),
     col="green",xlab="Barcode Count",
     ylab="Cumulative Proportion",
     main="Distribution of Barcode Counts-SRR3879604_1_bc.fastq")


plot(ecdf(second2$row.names),
     col="blue",xlab="Barcode Count",
     ylab="Cumulative Proportion",
     main="Distribution of Barcode Counts-SRR3879605_1_bc.fastq")

plot(ecdf(third2$row.names),
     col="red",xlab="Barcode Count",
     ylab="Cumulative Proportion",
     main="Distribution of Barcode Counts-SRR3879606_1_bc.fastq")

