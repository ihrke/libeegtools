errbarplot <- function (x, y = NULL, uiw, liw = uiw, add=F, ..., sfrac = 0.01)  {
  if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (is.null(y)) {
    if (is.null(x)) 
      stop("both x and y NULL")
    y <- as.numeric(x)
    x <- seq(along = x)
  }
  ui <- y + uiw;
  li <- y - liw;
  if(add)
    points(x, y, ylim = range(c(y, ui, li)), ...)
  else
    plot(x, y,...)# ylim = range(c(y, ui, li)), ...)

  smidge <- diff(par("usr")[1:2]) * sfrac
  segments(x, li, x, ui,...)
  x2 <- c(x, x)
  ul <- c(li, ui)
  segments(x2 - smidge, ul, x2 + smidge, ul, ...)
  invisible(list(x = x, y = y))
}

sem <- function(y){
  return(sd(y)/(sqrt(length(y))));
}


#d <- read.table("results/old/res_artdat_wgn.txt", header=T);
#d <- read.table("results/res_paramscan_denoisestat_beta0.0_100_samples2499.txt", header=T);
d <- read.table("data/res_paramscan_denoisestat_beta0.5_100_samples2499_snr-10.txt", header=T);
#"res_paramscan_denoisestat_beta2.0_100_sample2499_snr-10.txt"
# HIER noch unbedingt eine Condition mit 2^J langen Signalen nehmen!
#postscript(file="results/sigext.eps", horizontal=T);
postscript("results/sigext.eps", width = 13.0, height = 7.0,
                horizontal = FALSE, onefile = FALSE, paper = "special");
par(mfrow=c(1,2), bty="n",mai=c(1.0,1,0.2,0.2));
boxplot(avg.SNR. ~ e, data=d,
        main = "",
        xlab = "Extension Scheme",
        ylab = "SNR (db)",
	cex.axis=1.5,
	lwd=2.0,
	cex.lab=1.5);
boxplot(avg.RMSE. ~ e, data=d,
        xlab = "Extension Scheme",
        ylab = "RMSE",
	cex.axis=1.5,
	lwd=2.0,
	cex.lab=1.5
);
dev.off();

#postscript(file="results/eta.eps", horizontal=T);
postscript("results/eta.eps", width = 13.0, height = 7.0,
                horizontal = FALSE, onefile = FALSE, paper = "special");
par(mfrow=c(1,2), bty="n", mai=c(1.0,1,0.2,0.2));
boxplot(avg.SNR. ~ t, data=d,
        main = "",
        xlab = "Threshold",
        ylab = "SNR (db)",
	cex.axis=1.5,
	lwd=2.0,
	cex.lab=1.5,
        names=c("hard", "soft"));
boxplot(avg.RMSE. ~ t, data=d,
        xlab = "Threshold",
        ylab = "RMSE",	cex.axis=1.5,
	lwd=2.0,
	cex.lab=1.5,
        names=c("hard", "soft"));
dev.off();

#postscript(file="results/L.eps", horizontal=T);
postscript("results/L.eps", width = 13.0, height = 13.0,
                horizontal = FALSE, onefile = FALSE, paper = "special");
par(mfrow=c(2,2),bty="n", mai=c(0.7,1,0.2,0.2));

#postscript("results/L.eps", width = 16.0, height = 8.0,
#                horizontal = FALSE, onefile = FALSE, paper = "special");
#par(mfrow=c(1,2), bty="n", mai=c(0.8,1,0.2,0.2));

meanify_L<- function(d, func){
  m = c(); v = c(); m2=c(); v2=c();
  for(l in levels(as.factor(d$L))){
    l = as.numeric(l);
    m = c(m, mean(d$avg.SNR.[d$L==l & d$f==func]));
    v = c(v, sem(d$avg.SNR.[d$L==l & d$f==func]));
    m2 = c(m2, mean(d$avg.RMSE.[d$L==l & d$f==func]));
    v2 = c(v2, sem(d$avg.RMSE.[d$L==l & d$f==func]));
  }
  return(data.frame(m, v, m2, v2));
}

plot_L<- function(d, thisylim=c(-10,5), leg=F, rmse=F){
  m1 = meanify_L(d, "conventional");
  m2 = meanify_L(d, "ti");
  m3 = meanify_L(d, "sureshrink");
  m4 = meanify_L(d, "heursure");
  
  colors = c("black", "red", "green", "brown");
  pchs = c(15:19);
  lwds = 2.0;
  ltys = c(1,1,1,1);
  mycex = 2.0;
  if(rmse){
	m1$m = m1$m2;
	m1$v = m1$v2;
	m2$m = m2$m2;
	m2$v = m2$v2;
	m3$m = m3$m2;
	m3$v = m3$v2;
	m4$m = m4$m2;
	m4$v = m4$v2;
  }

  errbarplot(m1$m, NULL, m1$v, col=colors[1],
             xlab = "L", type="b",
             ylab = "SNR (db)", 
	     ylim=thisylim,
	     cex.axis=1.8,
	     cex.lab=2.0,
	     pch=pchs[1], 
	     lwd=lwds,
             lty=ltys[1],
	     cex=mycex
	     );
  errbarplot(m2$m, NULL, m2$v, add=T,  type='b',col=colors[2], pch=pchs[2], lwd=lwds,lty=ltys[2],cex=mycex);
  errbarplot(m3$m, NULL, m3$v, add=T, type='b', col=colors[3], pch=pchs[3], lwd=lwds,lty=ltys[3],cex=mycex);
  errbarplot(m4$m, NULL, m4$v, add=T,  type='b',col=colors[4], pch=pchs[4], lwd=lwds,lty=ltys[4],cex=mycex);

  if(leg){
    legend("topright", c("conventional", "ti", "sureshrink", "heursure"), lty=ltys,col=colors,pch=pchs, cex=mycex, bty="n");
  }
}

plot_L(d, thisylim=c(-10,5), T);
#plot_L(d, thisylim=c(3,17), F, T);
d2 <- read.table("data/res_paramscan_denoisestat_beta1.0_100_samples2499_snr-10.txt", header=T);
d3 <- read.table("data/res_paramscan_denoisestat_beta1.5_100_samples2499_snr-10.txt", header=T);
d4 <- read.table("data/res_paramscan_denoisestat_beta2.0_100_samples2499_snr-10.txt", header=T);
plot_L(d2, thisylim=c(-10,-2));
plot_L(d3, thisylim=c(-15,-4), T);
plot_L(d4, thisylim=c(-10,-3));

dev.off();




postscript("results/theta.eps", width = 13.0, height = 7.0,
                horizontal = FALSE, onefile = FALSE, paper = "special");
par(mfrow=c(1,2),bty="n");
#th <- read.table("results/res_theta.txt", header=T);
th <- read.table("data/res_timewarp_beta1.0_100_samples2499_snr-10.txt", header=T);

z1 = matrix(th$avg.SNR., nrow=11,ncol=11);
z2 = matrix(th$avg.RMSE., nrow=11,ncol=11);


persp(seq(0,1,0.1), seq(0,1,0.1), z1,
      xlab = '', ylab = '', zlab = '', main = NULL, sub = NULL,
      theta = -15, phi = 30, r = 50, d = 1, scale = TRUE, expand = 1,
      col = "green", border = NULL, ltheta = -135, lphi = 0, shade = 0.75,
      box = TRUE, axes = TRUE, nticks = 5, ticktype = "detailed");

persp(seq(0,1,0.1), seq(0,1,0.1), z2,
      xlab = '', ylab = '', zlab = '', main = NULL, sub = NULL,
      theta = 105, phi = 20, r = 50, d = 1, scale = TRUE, expand = 1,
      col = "green", border = NULL, ltheta = -135, lphi = 0, shade = 0.75,
      box = TRUE, axes = TRUE, nticks = 5, ticktype = "detailed");


dev.off();




#errbarplot(m, NULL, v, labels=l);

#mRMSE = matrix(NA, nrow=15, ncol=3);
#mSNR = matrix(NA, nrow=15, ncol=3);

#map = c('conventional', 'ti', 'sureshrink');
#for(i in 1:3){
#  mRMSE[,i]=c(
#        mean(d$avg.RMSE.[d$t=='h' & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$t=='s' & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$e==0 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$e==1 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==0 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==1 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==2 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==3 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==4 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==5 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==6 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==7 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==8 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==9 & d$f==map[i]]),
#        mean(d$avg.RMSE.[d$L==10 & d$f==map[i]])
#  );
#  mSNR[,i]=c(
#        mean(d$avg.SNR.[d$t=='h' & d$f==map[i]]),
#        mean(d$avg.SNR.[d$t=='s' & d$f==map[i]]),
#        mean(d$avg.SNR.[d$e==0 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$e==1 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==0 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==1 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==2 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==3 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==4 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==5 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==6 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==7 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==8 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==9 & d$f==map[i]]),
#        mean(d$avg.SNR.[d$L==10 & d$f==map[i]])
 # );
#  
#}##

#colors = c(rep('red', 2), #t
#  rep('blue', 2), #e
#  rep('green', 11)); #L
#h = barplot(mRMSE, names.arg=map,beside=T, col=colors);
#title("RMSE, parameter t, e, L");
#x11()
##par(ask=T);
#h = barplot(mSNR, names.arg=map,beside=T, col=colors);
#title("SNR, parameter t, e, L");#


