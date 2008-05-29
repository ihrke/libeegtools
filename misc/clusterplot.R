#fname = "clustdat_test_cleaned_ttjj.raw_tw.tsv";
#fname = "clustdat_test_cleaned_ttjj.raw_euclid.tsv";

#fname1 = "clustdat_artclusters2.raw_euclid.tsv"
#fname2 = "clustdat_artclusters2.raw_tw.tsv"

fname1 = "clustdat_tt_filtered.raw_euclid.tsv"
fname2 = "clustdat_tt_filtered.raw_tw.tsv"


dm = read.table(fname1, header=F);
d = matrix(as.numeric(as.matrix(dm)),nrow=length(dm[,1]), ncol=length(dm[1,]));

dm2 = read.table(fname2, header=F);
d2 = matrix(as.numeric(as.matrix(dm2)),nrow=length(dm2[,1]), ncol=length(dm2[1,]));


#---------------------------------------------------------------------------------------
postscript(sprintf("cplot_%s.eps",fname1), width=6, height=6, horizontal=F, paper="special");
fit1 <- hclust(as.dist(d), method="ward")
plot(fit1) # display dendogram
groups1 <- cutree(fit1, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
clusters = rect.hclust(fit1, k=2, border="red")
dev.off();


postscript(sprintf("heatmapplot_%s.eps",fname1), width=6, height=6, horizontal=F, paper="special");
heatmap(d);
dev.off();

#---------------------------------------------------------------------------------------

postscript(sprintf("cplot_%s.eps",fname2), width=6, height=6, horizontal=F, paper="special");
fit2 <- hclust(as.dist(d2), method="ward")
plot(fit2) # display dendogram
groups2 <- cutree(fit2, k=2) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
clusters = rect.hclust(fit2, k=2, border="red")
dev.off();


postscript(sprintf("heatmapplot_%s.eps",fname2), width=6, height=6, horizontal=F, paper="special");
heatmap(d2);
dev.off();
