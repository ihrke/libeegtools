# This is an evaluation script to accompany progs/t_dtw.
#  It clusters electrodes according to the difference of the
#  warppaths.


source("functions.R");
source("channels.R");

cluster_channels_with_pathdist <- function(fname){
  
  dm = read.table(fname, header=F);
  d = matrix(as.numeric(as.matrix(dm)),nrow=length(dm[,1]), ncol=length(dm[1,]));
  
  printf("writing file '%s'\n", sprintf("channelclustmap_%s.eps",fname));
  postscript(sprintf("channelclustmap_%s.eps",fname), width=12, height=12, horizontal=F, paper="special");

  fit1 <- hclust(as.dist(d), method="ward")
  plot(fit1, labels=channels_64eog.name) # display dendogram
  #groups1 <- cutree(fit1, k=2) # cut tree into 5 clusters
                                        # draw dendogram with red borders around the 5 clusters
  clusters = rect.hclust(fit1, k=2, border="red")

  # ---------------
  
  heatmap(d, labRow=channels_64eog.name,labCol=channels_64eog.name);

  # -------------
  
  cols = c("red", "green", "blue", "magenta", "black", "yellow");
  # loop over num_clust
  for( this.k in 1:6 ){
    clusters <- pam(as.dist(d), this.k); #cutree(fit1, k=this.k)
    plot.new();
    plot.window(c(-100, 100), c(-100, 100));
    title(main=sprintf("K = %i", this.k));
    printf("tk = %i\n", this.k);
    for( k2 in 1:this.k ){
#      printf("  k2=%i\n", k2);
 #     print(clusters);
      c = which(clusters$clustering==k2);
      text( channels_64eog.y[c], channels_64eog.x[c], labels=channels_64eog.name[c], col=cols[k2] );
    }
  }
  
  dev.off();
}
