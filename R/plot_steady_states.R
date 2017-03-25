plot_steady_states <- function(ss, png=NULL) {
  
  # Load libraries
  require(reshape2)
  require(dplyr)
  require(coRal)
  
  # For runs that did NOT reach a steady state of positive growth, set all results to zero
  ss[which(ss$gr <= 0), c("gr", "sh")] <- NA#, "hl", "sl", "ee", "pl")] <- NA
  # For runs where steady state value of S:H > 0.5, set S:H to 0.5
  ss[which(ss$sh > 0.5), "sh"] <- 0.5
  
  # Reshape steady state results into matrices for plotting cont
  sh <- acast(ss, L~N, value.var="sh")
  gr <- acast(ss, L~N, value.var="gr")

  # Build plotting function
  imagef <- function(r, bin, main, col) {
    par(tcl=-0.2, cex.main=1, cex.axis=0.75, mar=c(2,2,3,3), mgp=c(0,0.1,0))
    if (any(r<0, na.rm=T) & any(r>0, na.rm=T)) {
      nb <- rev(seq(0, -1*min(r[which(r!=0)], na.rm=T)+bin, bin))*-1  # negative portion of breaks
      pb <- seq(0, max(r[which(r!=0)], na.rm=T)+bin, bin)
      breaks <- c(nb, pb[-1])
    } else {
      breaks <- seq(min(r[which(r!=0)], na.rm=T), max(r, na.rm=T)+bin, bin)
    }
    ncolors <<- length(breaks) - 1
    if (col[[1]]=="grayscale") {
      colors <- gray.colors(n=ncolors, start=0.9, end=0.1)
    } else if (length(table(breaks>0))>1) {
      dir <- which.max(table(breaks>0))
      nn <- max(table(breaks>0))*2
      rc1 <- colorRampPalette(colors=c(col[[1]], "gray90"), space="Lab")(nn/2)
      rc1 <- rc1[-(length(rc1))]
      rc2 <- colorRampPalette(colors=c("gray90", col[[2]]), space="Lab")((nn/2)+1)[-1]
      colors <- c(rc1, rc2)
      if (dir==1) {
        colors <- colors[1:ncolors]
      } else {
        colors <- colors[(length(colors)-ncolors+1):length(colors)]
      }
    } else if (breaks[1] < 0) {
      colors <- colorRampPalette(colors=c(col[[1]], "gray90"), space="Lab")(ncolors)
    } else {
      colors <- colorRampPalette(colors=c("gray90", col[[2]]), space="Lab")(ncolors)
    }
    image(t(r), breaks=breaks, col=colors, xaxt="n", yaxt="n")
    #title(main=main, line=0.5, adj=-0.5, outer=F, xpd=T)
    conv <- diff(grconvertX(0:1, 'inches', 'user')) # length of one inch in user coordinates
    text(x=par("usr")[1]-par("mai")[2]*conv/1.1, y=par("usr")[4]+par("mai")[3]/2.5*conv, adj=c(0,0), 
         labels=main, xpd=NA, font=2, cex=1)
    mtext(side=1, text=expression("[DIN]"~(mol~L^{-1})), line=1, cex=0.75)
    axis(side=1, at=seq(0,1,0.125), 
         labels=quantile(seq(as.numeric(first(colnames(r))),
                             as.numeric(last(colnames(r))),
                             length.out=100),
                         probs=seq(0,1,0.125)))
    mtext(side=2, text=expression("Light" ~ (mol~photons~m^{2}~d^{-1})), line=1, cex=0.75)
    axis(side=2, at=seq(0,1,0.1), 
         labels=quantile(seq(as.numeric(first(rownames(r))),
                             as.numeric(last(rownames(r))),
                             length.out=100),
                         probs=seq(0,1,0.1)))
    par(new=T, mar=c(3,2,1,0))
    
    # Plot color key
    w <<- 0.05
    xl <<- par("usr")[2] * (1 + w)
    yb <<- par("usr")[3]
    h <<- (par("usr")[4] - par("usr")[3])/ncolors
    for(i in 1:ncolors) {
      rect(xl, yb+h*(i-1), xl+w, yb+h*i, col=colors[i], xpd=NA)
    }
    axis(side=4, at=seq(par("usr")[3], par("usr")[4], h), labels=round(as.numeric(as.character(breaks)), 2), 
         pos=xl+w, padj=0.5, hadj=0, las=1, tck=0)
  }
  
  # Plot steady states (save to .png if filename given to png argument)
  if (!is.null(png)) png(png, width=7, height=3.5, units="in", res=300)
  par(mfrow=c(1,2), mar=c(3,3,3,3))
  imagef(gr, bin=0.01, main="A. Specific growth", col="grayscale")
  text(x=xl+w, y=(yb+h*ncolors)*1.05, labels=expression(d^{-1}), xpd=T)
  print(c(xl+w, yb+h*ncolors))
  #imagef(hl, bin=0.25, main="B. Host growth limitation", col=list("#67001F","#053061"))
  imagef(sh, bin=0.05, main="B. Symbiont:host biomass", col="grayscale")
  text(x=xl+w, y=(yb+h*ncolors)*1.05, labels="S:H", xpd=T)
  #imagef(pl, bin=0.5, main="D. Photosynthesis limitation", col=list("#67001F", "gray30"))
  if (!is.null(png)) dev.off()
}
