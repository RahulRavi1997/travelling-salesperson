wd = "/home/local/ZOHOCORP/rahul-8808/Downloads/tsp"
numberOfNodes = 5000
testcase = 10


nodes = matrix(c(runif(numberOfNodes*2)*100), ,ncol = 2)
inputfile = paste(wd, "/input/testcase_" ,testcase, sep="")
write.table(nodes ,file = inputfile, sep = " ", row.names = FALSE, col.names = FALSE)

mydata <- read.table(inputfile)
print(mydata)
