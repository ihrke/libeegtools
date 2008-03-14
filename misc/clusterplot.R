fname = "clustertest2.txt";
dm = read.table(fname, header=F);
d = matrix(as.numeric(as.matrix(dm)),nrow=length(dm[,1]), ncol=length(dm[1,]));
#d = as.matrix(dm);

postscript("test.ps", width=6, height=6, horizontal=F, paper="special");
fit <- hclust(as.dist(d), method="ward")
plot(fit) # display dendogram
groups <- cutree(fit, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
clusters = rect.hclust(fit, k=2, border="red")
dev.off();
