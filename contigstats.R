lengths = read.table('~/length4a.txt') # include appropriate path to this file
# gives a quick view of the distribution of contig lengths
quantile(lengths[,1])
mean(lengths[,1])
# calculates how many total base pairs were assembled
sum(lengths[,1])
# calculates how many base pairs were assembled into contigs longer than 500bp for Community 4a
sum(lengths[1:474,])
