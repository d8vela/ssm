library(gplots)
#cat("heatmap.r","\n")

cmd_args = commandArgs();
args = commandArgs(TRUE)
data <- read.table(paste(args[1],args[2],sep = ''),header=T,na.strings = c("na","NONE","NOTE"),check.names=FALSE)
row.names(data) <- data$position
data <- data[,2:length(colnames(data))]
data = t(data)

pairs.breaks <- seq(-5, 5, by=0.125)

data_matrix <- data.matrix(log2(data))

data2 <- read.table(paste(args[1],args[2],sep = ''),header=T,check.names=FALSE)
data2 <- data2[,2:length(colnames(data2))]
data2 = t(data2)

labels <- matrix(data=' ',nrow=dim(data_matrix)[1],ncol=dim(data_matrix)[2])
for (i in 1:dim(labels)[1]) { 
  for (j in 1:dim(labels)[2]) { 
	if (data2[i,j] == 'NOTE')
    		{labels[i,j] <- '\''}
	else if (data2[i,j] == 'NONE')
    		{labels[i,j] <- '.'}
	else
		{labels[i,j] <- NA }
 } 
} 

for (x in 1:nchar(args[3])) { 
	f <- substring(args[3], x:x, x:x)
		if(f == 'A')
			{labels[1,x] <- 'A'}
		if(f == 'V')
			{labels[2,x] <- 'V'}
		if(f == 'I')
			{labels[3,x] <- 'I'}
		if(f == 'L')
			{labels[4,x] <- 'L'}
		if(f == 'M')
			{labels[5,x] <- 'M'}
		if(f == 'F')
			{labels[6,x] <- 'F'}
		if(f == 'Y')
			{labels[7,x] <- 'Y'}
		if(f == 'W')
			{labels[8,x] <- 'W'}
		if(f == 'S')
			{labels[9,x] <- 'S'}
		if(f == 'T')
			{labels[10,x] <- 'T'}
		if(f == 'N')
			{labels[11,x] <- 'N'}
		if(f == 'Q')
			{labels[12,x] <- 'Q'}
		if(f == 'R')
			{labels[13,x] <- 'R'}
		if(f == 'H')
			{labels[14,x] <- 'H'}
		if(f == 'K')
			{labels[15,x] <- 'K'}
		if(f == 'D')
			{labels[16,x] <- 'D'}
		if(f == 'E')
			{labels[17,x] <- 'E'}
		if(f == 'C')
			{labels[18,x] <- 'C'}
		if(f == 'G')
			{labels[19,x] <- 'G'}
		if(f == 'P')
			{labels[20,x] <- 'P'}
	}


data_matrix[is.na(data_matrix)] <- 0

postscript(paste(args[1],'heatmap_',args[2], '.eps', sep = ''),horizontal = TRUE,onefile = TRUE, paper = "special", height = 10, width = 10)
#pdf(paste(args[1],'heatmap_',args[2], '.pdf', sep = ''), width=12, height=7)
#data_heatmap <- heatmap.2(data_matrix, Rowv=TRUE, Colv=TRUE, cellnote=labels,na.color=par("bg"),notecol="white",trace=c("none"), dendrogram = c("both"), 
data_heatmap <- heatmap.2(data_matrix, Rowv=FALSE, Colv=FALSE, cellnote=labels,na.color=par("bg"),notecol="white",trace=c("none"), dendrogram = c("none"), 
key=FALSE, col = redgreen(79), breaks=pairs.breaks[1:80],symbreaks=20,
sepwidth=c(0.001, 0.001),
sepcolor="white",
colsep=1:ncol(data_matrix),
rowsep=1:nrow(data_matrix),
scale=c("none"), 
margins=c(3,3),
keysize = 2.0,
xlab='',
ylab='',
font.main=15,
lmat=rbind( c(0, 3), c(2,1), c(0,4) ), 
lhei=c(0.5, 3, 1 ),
lwid=c(0.5, 6 ),
tracecol="cyan",
na.rm=TRUE,
notecex=1.5,
cexRow=1.4,
cexCol=1.0
)

title(main = list(args[2], cex=2, font=1.5))
dev.off()
