##plotting sequence. Can be turned into a function
##Plot the potential circles
#Circle by circle
#Version 3 with plotting the gene names on the blocks.

plot_circles12 <-function(circles,def1,def2, outputname,org1,org2,size_cutoff) { #for circles with multiple blocks, circles_12 is a list with dataframes containing block info
  require(RColorBrewer)

  #png(filename = "1-2-_potential_circles_Synchro_m%02d.png",
  #    width = 1000, height = 1000, units = "px", pointsize = 12,
  #    bg = "white")
  #Test:
  #i=2
  #circle<-circles1[[1]][[1]][["1-2-"]][[i]]
  #def1<-annotations[[1]]
  #def2<-annotations[[2]]
  #org1<-organisms[1]
  #org2<-organisms[2]

  
  pdf(file = outputname,
      pointsize = 10,paper="a4", bg = "white")
  
  pal <- brewer.pal(9, "Set1")
  
    for (i in 1:length(circles)) {
      #skip for test:
      circle<-circles[[i]]
      
      no<-1:dim(circle)[1]
      organism2<-circle[order(circle$start2),] %>% select(contains("2")) #for some plots, the order has to be sorted
      left2<-vector(mode = "numeric",length = length(no))
      right2<-vector(mode = "numeric",length = length(no))
      left1<-vector(mode = "numeric",length = length(no))
      right1<-vector(mode = "numeric",length = length(no))
      length2 <- max(circle$stop2)-min(circle$start2)
      length1 <- max(circle$stop1)-min(circle$start1)
      largestLength = max (length2, length1)
    
    if(length2 > length1) {
      for (n in no) {
        left2[n]  <- 1 - (max(circle[, "stop2"])-circle[n, "start2"])/largestLength
        right2[n] <- 1 - (max(circle[, "stop2"])-circle[n, "stop2"])/largestLength
        left1[n]  <- (1-((length2-length1)/2)/largestLength) - (max(circle[, "stop1"])-circle[n, "start1"])/largestLength
        right1[n] <- (1-((length2-length1)/2)/largestLength) - (max(circle[, "stop1"])-circle[n, "stop1"])/largestLength
      }
      
    } else if(length2 < length1) {
      for (n in no){
      left1[n]  <- 1 - (max(circle[, "stop1"])-circle[n, "start1"])/largestLength
      right1[n] <- 1 - (max(circle[, "stop1"])-circle[n, "stop1"])/largestLength
      left2[n]  <- (1-((length1-length2)/2)/largestLength) - (max(circle[, "stop2"])-circle[n, "start2"])/largestLength
      right2[n] <- (1-((length1-length2)/2)/largestLength) - (max(circle[, "stop2"])-circle[n, "stop2"])/largestLength
      }
    }
    
    plot.new()
    height1 = 0.7
    height2 = 0.3
    
    #Plot the blocks
    for (n in no) {
      segments(left1[n], height1, right1[n], height1, lwd=3) #plot the organism1 blocks
      segments(left2[n], height2, right2[n], height2, lwd=3) #plot the organism2 blocks
    }
    for (n in no) {
      polygon(c(left1[n], right2[n], left2[n], right1[n]), c(height1-0.01, height2+0.01, height2+0.01, height1-0.01), col = "lightgray", border = pal[n],fillOddEven=F)
      
    }#For 1-2-, both are inversed

    left2_sort<-left2 %>%sort
    right2_sort<-right2 %>%sort
    
    #Add the breakpoint intervals
    for (n in no[-length(no)]) {
      nobph <-(circle[n+1, "start1"]-circle[n, "stop1"])<=0
      if(!nobph) {
      segments(right1[n], height1+0.05, left1[n+1], height1+0.05, lwd=3, col="darkgray") #circle breakpoint
      segments(right1[n], height1, right1[n], height1+0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
      segments(left1[n+1], height1, left1[n+1], height1+0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
      #label the breakpoints with size
      text((right1[n]+left1[n+1])/2, 0.9, paste(format(circle[n+1, "start1"]-circle[n, "stop1"], nsmall=1, big.mark=","),"bps"), col="gray50",
           cex = 1, srt=0, font=2)
      }
      nobpm <-(circle[n+1, "start2"]-circle[n, "stop2"])<=0
      if(!nobpm) {
      segments(right2_sort[n], height2-0.05, left2_sort[n+1], height2-0.05, lwd=3, col="darkgray") #open circle breakpoint
      segments(right2_sort[n], height2, right2_sort[n], height2-0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
      segments(left2_sort[n+1], height2, left2_sort[n+1], height2-0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
      #label the breakpoints with size
      text((right2_sort[n]+left2[n+1])/2, 0.1, paste(format(organism2[n+1, "start2"]-organism2[n, "stop2"], nsmall=1, big.mark=","),"bps"),col="gray50",
           cex = 1, srt=0, font=2)
      }
    }
    
    #Get all gene info in blocks
    genes <- list()
    genes[["blocks1"]] <-sapply(no, function(n) {
      def1[def1$chr==circle[1,"chr1"] & #genes in: same chromosome
                        def1$start>=circle[n,"start1"] & def1$end<=circle[n,"stop1"] & #the correct interval
                        def1$size>=size_cutoff,] #pick only genes about the cutoff
    }, simplify = F) %>% reduce(bind_rows) #dplyr bind_rows appears to be the only row binder that works
    
    genes[["breakpoints1"]] <-sapply(no[-length(no)], function(n) {
      def1[def1$chr==circle[1,"chr1"] & #genes in: same chromosome
              def1$start>circle[n,"stop1"] & def1$end<circle[n+1,"start1"] & #the correct interval
              def1$size>=size_cutoff,]
    }, simplify = F) %>% reduce(bind_rows)
    
    genes[["blocks2"]] <-sapply(no, function(n) {
      def2[def2$chr==organism2[1,"chr2"] & #genes in: same chromosome
              def2$start>= organism2[n,"start2"] & def2$end<=organism2[n,"stop2"] & #the correct interval
              def2$size>=size_cutoff,]
    }, simplify = F) %>% reduce(bind_rows) #dplyr bind_rows appears to be the only row binder that works
    
    genes[["breakpoints2"]] <-sapply(no[-length(no)], function(n) {
      def2[def2$chr==organism2[1,"chr2"] & #genes in: same chromosome
              def2$start>organism2[n,"stop2"] & def2$end<organism2[n+1,"start2"] & #the correct interval
              def2$size>=size_cutoff,]
    }, simplify = F) %>% reduce(bind_rows)
    
    
              ##Add organism1 genes to blocks
    if (dim(genes[["blocks1"]])[1]>0) {
      genenames<-genes[["blocks1"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["blocks1"]])[1]) {
        gene_l<-(genes[["blocks1"]][r,"start"]-min(circle[, "start1"]))/largestLength
        gene_r<-(genes[["blocks1"]][r,"end"]-min(circle[, "start1"]))/largestLength
        polygon(c(gene_l+left1[1], gene_l+left1[1], gene_r+left1[1],gene_r+left1[1]),c(height1-0.01,height1+0.01,height1+0.01,height1-0.01), height1,col="white",border = "black")
        positions[r]<-((gene_l+gene_r)/2)+min(left1)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height1, length(positions)),genenames, col.label="darkgreen", cex.label = 0.6)
      }#text(((gene_l+gene_r)/2)+min(left1), height1+0.07,genes[["blocks1"]][r,"symbol"], offset = 0.05, col="black", cex = .6, srt=90)

      #add organism1 genes in breakpoints
    if(dim(genes[["breakpoints1"]])[1]>0) {
      genenames<-genes[["breakpoints1"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["breakpoints1"]])[1]) {
      gene_l<-(genes[["breakpoints1"]][r,"start"]-min(circle[, "start1"]))/largestLength
      gene_r<-(genes[["breakpoints1"]][r,"end"]-min(circle[, "start1"]))/largestLength
      polygon(c(gene_l+left1[1], gene_l+left1[1], gene_r+left1[1],gene_r+left1[1]),c(height1+0.04,height1+0.06,height1+0.06,height1+0.04), height1,col="white",border = "darkgray")
      positions[r]<-((gene_l+gene_r)/2)+min(left1)
          }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height1+0.05, length(positions)),genenames, col.label="black", cex.label = 0.6)
      }#text(((gene_l+gene_r)/2)+min(left1), height1+0.08,genes[["breakpoints1"]][r,"symbol"], offset = 0.05, col="darkgray", cex = .6, srt=90)

    
    ##Add organism2 genes to blocks
    if(dim(genes[["blocks2"]])[1]>0) {
      genenames<-genes[["blocks2"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["blocks2"]])[1]) {
      gene_l<-(genes[["blocks2"]][r,"start"]-min(circle[, "start2"]))/largestLength
      gene_r<-(genes[["blocks2"]][r,"end"]-min(circle[, "start2"]))/largestLength
      polygon(c(gene_l+min(left2), gene_l+min(left2), gene_r+min(left2),gene_r+min(left2)),c(height2-0.01,height2+0.01,height2+0.01,height2-0.01), height1,col="white",border = "black")
      positions[r]<-((gene_l+gene_r)/2)+min(left2)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height2, length(positions)),genenames, col.label="purple", cex.label = 0.6)
      }#text(((gene_l+gene_r)/2)+min(left2), height2-0.07,genes[["blocks2"]][r,"symbol"], offset = 0.05, col="black", cex = .6, srt=90)
    
    #add organism2 genes in breakpoints
    if(dim(genes[["breakpoints2"]])[1]>0) {
      genenames<-genes[["breakpoints2"]][["symbol"]]
      positions<-numeric(length = length(genenames))      
    for (r in 1:dim(genes[["breakpoints2"]])[1]) {
      gene_l<-(genes[["breakpoints2"]][r,"start"]-min(circle[, "start2"]))/largestLength
      gene_r<-(genes[["breakpoints2"]][r,"end"]-min(circle[, "start2"]))/largestLength
      polygon(c(gene_l+min(left2), gene_l+min(left2), gene_r+min(left2),gene_r+min(left2)),c(height2-0.04,height2-0.06,height2-0.06,height2-0.04), height1,col="white",border = "darkgray")
      positions[r]<-((gene_l+gene_r)/2)+min(left2)
          }#text(((gene_l+gene_r)/2)+min(left2), height2,genes[["breakpoints2"]][r,"symbol"], offset = 0.05, col="darkgray", cex = .6, srt=90)
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height2-0.05, length(positions)),genenames, col.label="black", cex.label = 0.6)
      
      }
    
       #label species and chromosome number
    text(0.5, 1, paste(org1,"\n Chr",circle[1,"chr1"],": ",format(min(circle$start1), nsmall=1, big.mark=","),"-",format(max(circle$stop1), nsmall=1, big.mark=",")),col="black",cex = .8, font=3)
    text(0.5, 0, paste(org2,"\n Chr",circle[1,"chr2"],": ",format(min(circle$start2), nsmall=1, big.mark=","),"-",format(max(circle$stop2), nsmall=1, big.mark=",")),col="black",cex = .8, font=3)
    }#loop complete, next circle plot
  dev.off()
} #function complete

plot_circles21 <-function(circles,def1,def2,outputname,org1,org2,size_cutoff) { #for circles with multiple blocks, circles_12 is a list with dataframes containing block info
  require(RColorBrewer)
#  png(filename = "2+1+_potential_circles_Synchro_m%02d.png",
#      width = 1000, height = 1000, units = "px", pointsize = 12,
#      bg = "white")
  #Test:
  #i=2
  #circle<-circles1[[1]][[1]][["2+1+"]][[i]]
  #def1<-annotations[[1]]
  #def2<-annotations[[2]]
  #org1<-organisms[1]
  #org2<-organisms[2]
  
  pdf(file = outputname,
      pointsize = 10,paper="a4", bg = "white")
    pal <- brewer.pal(9, "Set1")
  
  for (i in 1:length(circles)) {
    #don't execute next line for test
    circle<-circles[[i]]
    no<-1:dim(circle)[1]
    organism2<-circle[order(circle$start2),] %>% select(contains("2")) #for some plots, the order has to be sorted
    left2<-vector(mode = "numeric",length = length(no))
    right2<-vector(mode = "numeric",length = length(no))
    left1<-vector(mode = "numeric",length = length(no))
    right1<-vector(mode = "numeric",length = length(no))
    length2 <- max(circle$stop2)-min(circle$start2)
    length1 <- max(circle$stop1)-min(circle$start1)
    largestLength = max (length2, length1)
    
    if(length2 > length1) {
      for (n in no) {
        left2[n]  <- 1 - (max(circle[, "stop2"])-circle[n, "start2"])/largestLength
        right2[n] <- 1 - (max(circle[, "stop2"])-circle[n, "stop2"])/largestLength
        left1[n]  <- (1-((length2-length1)/2)/largestLength) - (max(circle[, "stop1"])-circle[n, "start1"])/largestLength
        right1[n] <- (1-((length2-length1)/2)/largestLength) - (max(circle[, "stop1"])-circle[n, "stop1"])/largestLength
      }
      
    } else if(length2 < length1) {
      for (n in no){
        left1[n]  <- 1 - (max(circle[, "stop1"])-circle[n, "start1"])/largestLength
        right1[n] <- 1 - (max(circle[, "stop1"])-circle[n, "stop1"])/largestLength
        left2[n]  <- (1-((length1-length2)/2)/largestLength) - (max(circle[, "stop2"])-circle[n, "start2"])/largestLength
        right2[n] <- (1-((length1-length2)/2)/largestLength) - (max(circle[, "stop2"])-circle[n, "stop2"])/largestLength
      }
    }
    
    plot.new()
    height1 = 0.7
    height2 = 0.3

    ##Difference until here.
    
    #Plot the blocks
    for (n in no) {
      segments(left1[n], height1, right1[n], height1, lwd=3) #plot the organism1 blocks
      segments(left2[n], height2, right2[n], height2, lwd=3) #plot the organism2 blocks
    }

    for (n in no) {
      polygon(c(left1[n], left2[n], right2[n], right1[n]), c(height1-0.01, height2+0.01, height2+0.01, height1-0.01), col = "lightgray", border = pal[n],fillOddEven=F)
    }
    
    
    left2_sort<-left2 %>%sort
    right2_sort<-right2 %>%sort
    
    #Add the breakpoint intervals
    
    #Add the breakpoint intervals
    for (n in no[-length(no)]) {
      nobph <-(circle[n+1, "start1"]-circle[n, "stop1"])<=0
      if(!nobph) {
        segments(right1[n], height1+0.05, left1[n+1], height1+0.05, lwd=3, col="darkgray") #circle breakpoint
        segments(right1[n], height1, right1[n], height1+0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
        segments(left1[n+1], height1, left1[n+1], height1+0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
        #label the breakpoints with size
        text((right1[n]+left1[n+1])/2, 0.9, paste(format(circle[n+1, "start1"]-circle[n, "stop1"], nsmall=1, big.mark=","),"bps"), col="gray50",
             cex = 1, srt=0, font=2)
      }
      nobpm <-(organism2[n+1, "start2"]-organism2[n, "stop2"])<=0
      if(!nobpm) {
        segments(right2_sort[n], height2-0.05, left2_sort[n+1], height2-0.05, lwd=3, col="darkgray") #open circle breakpoint
        segments(right2_sort[n], height2, right2_sort[n], height2-0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
        segments(left2_sort[n+1], height2, left2_sort[n+1], height2-0.05, lwd=2, lty = "dotted",col="darkgray") #circle breakpoint
        
        #label the breakpoints with size
        text((right2_sort[n]+left2_sort[n+1])/2, 0.1, paste(format(organism2[n+1, "start2"]-organism2[n, "stop2"], nsmall=1, big.mark=","),"bps"),col="gray50",
             cex = 1, srt=0, font=2)
      }
    }
    
    
    #Get all gene info in blocks
    genes <- list()
    genes[["blocks1"]] <-sapply(no, function(n) {
      def1[def1$chr==circle[1,"chr1"] & #genes in: same chromosome
              def1$start>=circle[n,"start1"] & def1$end<=circle[n,"stop1"] & #the correct interval
              def1$size>=size_cutoff,] #for now only pick the genes
    }, simplify = F) %>% reduce(bind_rows) #dplyr bind_rows appears to be the only row binder that works
    
    genes[["breakpoints1"]] <-sapply(no[-length(no)], function(n) {
      def1[def1$chr==circle[1,"chr1"] & #genes in: same chromosome
              def1$start>circle[n,"stop1"] & def1$end<circle[n+1,"start1"] & #the correct interval
              def1$size>=size_cutoff,] #for now only pick the genes
    }, simplify = F) %>% reduce(bind_rows)
    
    genes[["blocks2"]] <-sapply(no, function(n) {
      def2[def2$chr==organism2[1,"chr2"] & #genes in: same chromosome
              def2$start>=organism2[n,"start2"] & def2$end<=organism2[n,"stop2"] & #the correct interval
              def2$size>=size_cutoff,] #for now only pick the genes
    }, simplify = F) %>% reduce(bind_rows) #dplyr bind_rows appears to be the only row binder that works
    
    genes[["breakpoints2"]] <-sapply(no[-length(no)], function(n) {
      def2[def2$chr==organism2[1,"chr2"] & #genes in: same chromosome
              def2$start>organism2[n,"stop2"] & def2$end<organism2[n+1,"start2"] & #the correct interval
              def2$size>=size_cutoff,] #for now only pick the genes
    }, simplify = F) %>% reduce(bind_rows)
    
    
    ##Add organism1 genes to blocks
    
    if (dim(genes[["blocks1"]])[1]>0) {
      genenames<-genes[["blocks1"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["blocks1"]])[1]) {
      gene_l<-(genes[["blocks1"]][r,"start"]-min(circle[, "start1"]))/largestLength
      gene_r<-(genes[["blocks1"]][r,"end"]-min(circle[, "start1"]))/largestLength
      polygon(c(gene_l+left1[1], gene_l+left1[1], gene_r+left1[1],gene_r+left1[1]),c(height1-0.01,height1+0.01,height1+0.01,height1-0.01), height1,col="white",border = "black")
      positions[r]<-((gene_l+gene_r)/2)+min(left1)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height1, length(positions)),genenames, col.label="darkgreen", cex.label = 0.6)
      }
    
    #add organism1 genes in breakpoints
    if(dim(genes[["breakpoints1"]])[1]>0) {
      genenames<-genes[["breakpoints1"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["breakpoints1"]])[1]) {
      gene_l<-(genes[["breakpoints1"]][r,"start"]-min(circle[, "start1"]))/largestLength
      gene_r<-(genes[["breakpoints1"]][r,"end"]-min(circle[, "start1"]))/largestLength
      polygon(c(gene_l+left1[1], gene_l+left1[1], gene_r+left1[1],gene_r+left1[1]),c(height1+0.04,height1+0.06,height1+0.06,height1+0.04), height1,col="white",border = "darkgray")
      positions[r]<-((gene_l+gene_r)/2)+min(left1)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height1+0.05, length(positions)),genenames, col.label="black", cex.label = 0.6) 
      }
    
    ##Add organism2 genes to blocks
    if (dim(genes[["blocks2"]])[1]>0) {  
      genenames<-genes[["blocks2"]][["symbol"]]
      positions<-numeric(length = length(genenames))
    for (r in 1:dim(genes[["blocks2"]])[1]) {
      gene_l<-(genes[["blocks2"]][r,"start"]-min(circle[, "start2"]))/largestLength
      gene_r<-(genes[["blocks2"]][r,"end"]-min(circle[, "start2"]))/largestLength
      polygon(c(gene_l+min(left2), gene_l+min(left2), gene_r+min(left2),gene_r+min(left2)),c(height2-0.01,height2+0.01,height2+0.01,height2-0.01), height1,col="white",border = "black")
      positions[r]<-((gene_l+gene_r)/2)+min(left2)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height2, length(positions)),genenames, col.label="purple", cex.label = 0.6)
      }
    
    #add organism2 genes in breakpoints
    if(dim(genes[["breakpoints2"]])[1]>0) {
      genenames<-genes[["breakpoints2"]][["symbol"]]
      positions<-numeric(length = length(genenames))         
    for (r in 1:dim(genes[["breakpoints2"]])[1]) {
      gene_l<-(genes[["breakpoints2"]][r,"start"]-min(circle[, "start2"]))/largestLength
      gene_r<-(genes[["breakpoints2"]][r,"end"]-min(circle[, "start2"]))/largestLength
      polygon(c(gene_l+min(left2), gene_l+min(left2), gene_r+min(left2),gene_r+min(left2)),c(height2-0.04,height2-0.06,height2-0.06,height2-0.04), height1,col="white",border = "darkgray")
      positions[r]<-((gene_l+gene_r)/2)+min(left2)
    }
      positions<-positions[genenames!=""]
      genenames<-genenames[genenames!=""]
      addTextLabels(positions, rep(height2-0.05, length(positions)),genenames, col.label="black", cex.label = 0.6)
      }

    
    #label species and chromosome number
    text(0.5, 1, paste(org1,"\n Chr",circle[1,"chr1"],": ",format(min(circle$start1), nsmall=1, big.mark=","),"-",format(max(circle$stop1), nsmall=1, big.mark=",")),col="black",cex = .8, font=3)
    text(0.5, 0, paste(org2,"\n Chr",circle[1,"chr2"],": ",format(min(circle$start2), nsmall=1, big.mark=","),"-",format(max(circle$stop2), nsmall=1, big.mark=",")),col="black",cex = .8, font=3)

  }#loop complete, next circle plot
  dev.off()
} #function complete
