as.dataframe.synbuilder<-function(df=Hg38_mm10, gen1=hg38,gen2=mm10){
  df<-data.frame("ID"=df[grepl(">",df)],
                 "ID1"=df[5:length(df)][grepl(gen1,df[5:length(df)])],
                 "ID2"=df[5:length(df)][grepl(gen2,df[5:length(df)])])
  df<-data.frame("ID"=gsub(">","",df$ID),
                 "ID1"=gsub(" .*","",df$ID1),
                 "ID2"=gsub(" .*","",df$ID2),
                 "chr1"=str_extract(df$ID1,"(?<=\\.chr)[:alnum:]+(?=\\:)"),
                 "start1"=str_extract(df$ID1,"(?<=\\:)[:alnum:]+(?=-)") %>% as.numeric,
                 "stop1"=str_extract(df$ID1,"(?<=\\-)[:alnum:]+(?= )") %>% as.numeric,
                 "direction1"=str_extract(df$ID1,"(?<= )[:graph:]$"),
                 "chr2"=str_extract(df$ID2,"(?<=\\.chr)[:alnum:]+(?=\\:)"),
                 "start2"=str_extract(df$ID2,"(?<=\\:)[:alnum:]+(?=-)") %>% as.numeric,
                 "stop2"=str_extract(df$ID2,"(?<=\\-)[:alnum:]+(?= )") %>% as.numeric,
                 "direction2"=str_extract(df$ID2,"(?<= )[:graph:]$"))
  return(df)
}

find_circles<-function(df=MAMU) {
  #Test:
  #df<-random_blocks_HOMU[[365]]
  #options(warn = 0)
  synteny_dir<-sapply(2:nrow(df), function(i) {
     link<- if (df[i,"chr2"]==df[i-1,"chr2"] & df[i,"chr1"]==df[i-1,"chr1"]) { #both blocks are on the same chromosome
      if    (df[i-1,"direction2"]=="+" & df[i,"direction2"]=="+" & df$start2[i-1]<df$start2[i]& #it's not possible to have any additional blocks for 12 synteny
             !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i-1] & df[,"start2"] < df$start2[i])) {"1+2+"}  
      else if (df[i-1,"direction2"]=="-" & df[i,"direction2"]=="-" & df$start2[i-1]<df$start2[i]&
               !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i-1] & df[,"start2"] < df$start2[i])) {"1-2-"}      
      else if (df[i-1,"direction2"]=="+" & df[i,"direction2"]=="-" & df$start2[i-1]<df$start2[i]&
               !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i-1] & df[,"start2"] < df$start2[i])) {"1+2-"}
      else if (df[i-1,"direction2"]=="-" & df[i,"direction2"]=="+" & df$start2[i-1]<df$start2[i]&
               !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i-1] & df[,"start2"] < df$start2[i])) {"1-2+"}
#for 21 combinations, there are more options that need to be tested
      #2+1+
      else if (df[i-1,"direction2"]=="+" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
              !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i] & df[,"start2"] < df$start2[i-1])) {"2+1+"} #nothing between the blocks -->OK, or 

      else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks in the breakpoint in org2 that pairs with upstream of org1 block1 coordinates
              !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i])& #and no blocks downstream of org1 block2
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & 
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & #stuffs in the breakpoint have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "+")) { #and correct direction
              ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] ==  #the coordinates upstream of block1 have to be sorted
              sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) &             
                all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                      df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]),"2+1+m2","rj_2+1+_2")}

                else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
                !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are no blocks upstream of org1 block1 coordinates
                any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & 
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & #stuffs have to be on the same chromosome
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "+")) #and correct direction
                {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] ==    #the coordinates downstream of block2 have to be sorted
                       sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"])) &
                   all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                         df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]), "2+1+m1","rj_2+1+_2")} 
        
      else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates, same chromosome in org2
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & 
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & #stuffs have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "+")) #and correct direction   
              {ifelse(!any(max(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]) > #there can not be any org2 blocks starting from  block2
              min(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) &   #that are bigger than blocks extending from block1
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] == #coordinates upstream of block1 have to be sorted
              sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) & 
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] == #and the coordinates upstream of block2 have to be sorted
              sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"])) & 
                all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                      df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]) &
                all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                      df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]), "2+1+m3","rj_2+1+_3")}
      
    #2-1+
    else if (df[i-1,"direction2"]=="-" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
             !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i] & df[,"start2"] < df$start2[i-1])) {"2-1+"} #nothing between the blocks -->OK, or 
        else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
                any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
                !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & 
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & #stuffs have to be on the same org1 chromosome
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "+"))        
         {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] ==  #the coordinates upstream of block1 have to be sorted
                       sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"]))&
                   all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                         df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]),"2-1+m1","rj_2-1+_1")}
        
            else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
            !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are no blocks upstream of org1 block1 coordinates
             any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) & 
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & 
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & #stuffs have to be on the same chromosome
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "-"))  #and correct direction 
            {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] ==    #the coordinates downstream of block2 have to be sorted
                   sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]))&
               all(df[(i+1):(i+length(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                     df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]), "2-1+m2","rj_2-1+_2")}

    else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="+" & df$start2[i-1]>df$start2[i] &
            any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
            any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i])   & #and blocks downstream of org1 block2 coordinates
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"])   & 
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   & #stuffs have to be on the same org1 chromosome
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"direction2"]=="+") &   #  blocks starting from  block1 are + in org2   
            all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"direction2"]=="-")) #  blocks starting from  block2 are - in org2
          {ifelse(!max(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]) > #there can not be any org2 blocks starting from  block2
                    min(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"]) &   #that are bigger than blocks extending from block1
               all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] == #and the coordinates upstream of block1 have to be sorted
               sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) & 
               all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] == #and the coordinates upstream of block2 have to be sorted
               sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"])) & 
               all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                     df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]) &
               all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                     df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]), "2-1+m3","rj_2-1+_3")}
    
      #2-1-
      else if (df[i-1,"direction2"]=="-" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
               !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i] & df[,"start2"] < df$start2[i-1])) {"2-1-"} #nothing between the blocks -->OK, or 
      else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are no blocks upstream of org1 block1 coordinates
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"])   &#stuffs have to be on the same chromosome 
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "-")) #and correct direction 
            {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] ==    #the coordinates downstream of block2 have to be sorted
                     sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]))&
                 all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]), "2-1-m1","rj_2-1-_1")}
       
      else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
              !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i])  &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"])   &#stuffs have to be on the same chromosome 
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "-")) #and correct direction 
              {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] ==  #the coordinates upstream of block1 have to be sorted
                     sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) &
                 all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]),"2-1-m2","rj_2-1-_2")}

      else if(df[i-1,"direction2"]=="-" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & #stuffs have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]=="-"))   #  blocks in -
              {ifelse(!max(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]) > #there can not be any org2 blocks starting from  block2
                      min(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"]) &   #that are bigger than blocks extending from block1
                 all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] == #and the coordinates upstream of block1 have to be sorted
                 sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) & 
                 all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] == #and the coordinates upstream of block2 have to be sorted
                 sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"])) & 
                 all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]) &
                 all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]), "2+1+m3","rj_2-1-_3")}
      
      #2+1-
      else if (df[i-1,"direction2"]=="+" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
               !any(df[,"chr2"] == df$chr2[i] & df[,"start2"] > df$start2[i] & df[,"start2"] < df$start2[i-1])) {"2+1-"} #nothing between the blocks -->OK, or 
      else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are no blocks upstream of org1 block1 coordinates
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & #stuffs have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "+")) #downstream they should be +
               {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] ==    #the coordinates downstream of block2 have to be sorted
                     sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]))&
                     all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]), "2+1-","rj_2+1-_1")}
       
      else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
              !any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & #stuffs have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"])   &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"direction2"]== "-")) #upstream they should be -              
              {ifelse(all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] ==  #the coordinates upstream of block1 have to be sorted
                     sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"]))&
                 all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]),"2+1-m2","rj_2+1-_2")}
       
      else if(df[i-1,"direction2"]=="+" & df[i,"direction2"]=="-" & df$start2[i-1]>df$start2[i] &
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1]) & #if there are blocks upstream of org1 block1 coordinates
              any(df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i]) &
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr1"]== df[i,"chr1"]) & #stuffs have to be on the same chromosome
              all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1],"chr2"]== df[i,"chr2"]) & 
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"direction2"]=="-") &   #  blocks starting from  block1 are + in org2   
                all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"direction2"]=="+")) #  blocks starting from  block2 are - in org2              
               {ifelse(!max(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"]) > #there can not be any org2 blocks starting from  block2
                 min(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"]) &   #that are bigger than blocks extending from block1
                 all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"] == #and the coordinates upstream of block1 have to be sorted
                 sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start2"])) & 
                 all(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"] == #and the coordinates upstream of block2 have to be sorted
                 sort(df[df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start2"])) & 
                 all(df[(i+1):(i+length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"])),"start1"] == #and the org1 coordinates have to be directly downstream of block 2
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1> df$start1[i],"start1"]) &
                 all(df[(i-length(df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"])):(i-1),"start1"] == #and the org1 coordinates have to be directly upstream of block 1
                       df[df$chr1==df$chr1[i] & df$chr2==df$chr2[i] & df$start2> df$start2[i] & df$start2< df$start2[i-1] & df$start1< df$start1[i-1],"start1"]), "2+1-m3","rj_2+1-_3")}
    } #end if
    #print(i)
    #print(link)
    link
  })
  
  synteny_dir[sapply(synteny_dir,is.null)]<-NA
  synteny_dir<-c(NA,unlist(synteny_dir))
  synteny_dir %>% table
  
  circles <-list() #list of circles returned, add the raw calls for synteny
  df<-cbind(df,synteny_dir)
  circles_12_block1<-df[(which(grepl("^1-2-",df$synteny_dir))-1),]  # ignore the last column
  circles_12_block2<-df[which(grepl("^1-2-",df$synteny_dir)),]
  
  if ( dim(circles_12_block1)[1]>0 & dim(circles_12_block2)[1]>0 ) { #if there are no pairs found at all, skip sequence
    circles_12_block1[,"Number"]<-1:dim(circles_12_block1)[1]
    circles_12_block2[,"Number"]<-1:dim(circles_12_block2)[1]
    circles_12<-full_join(select(circles_12_block1, -"synteny_dir"),select(circles_12_block2, -starts_with(c("chr","direction"))),by="Number")
    colnames(circles_12)<-gsub("\\.x","_b1",colnames(circles_12))
    colnames(circles_12)<-gsub("\\.y","_b2",colnames(circles_12))
    #In this, there are no blocks used in two circle calls
    #intersect(circles_12$ID_b1,circles_12$ID_b2)
    
    duplicates<-(circles_12$ID_b1 %in% intersect(circles_12$ID_b1,circles_12$ID_b2)|circles_12$ID_b2 %in% intersect(circles_12$ID_b1,circles_12$ID_b2))
    duplicated_block<-rep(intersect(circles_12$ID_b1,circles_12$ID_b2),each=2)
    
    n=0
    circles_12[,"Duplicated_blocks"]<-sapply(seq_along(duplicates), function(i) {
      if (duplicates[i]) {
        n <<- n+1 #you can update a variable outside of the apply environment by using <<
        paste(duplicated_block[n])}
      else {"Unique"}
    })
    
    ##
    circles_12<-transform(circles_12, 
                          bp_size=circles_12$start2_b2-circles_12$stop2_b1,
                          total_size2=circles_12$stop2_b2-circles_12$start2_b1,
                          total_size1=circles_12$stop1_b2-circles_12$start1_b1)
    circles_12<-relocate(circles_12,"Number")
    

    #Identify if there are blocks flanking with the correct synteny, that could be part of the circle
    circles_12_OK<-lapply(1:dim(circles_12)[1], function(i) {
      #Are there any correct blocks upstream
        n<-which(df$start1==circles_12[i,"start1_b1"] &
                   df$chr1 == circles_12[i,"chr1"] &
                   df$chr2 == circles_12[i,"chr2"])#the block2 position in df
        if(!n == 1 & df[n,"start2"] != min(df[df$chr2==df[n,"chr2"],"start2"])) {       
        link<-if(df[n-1,"start2"]==sort(df[df$chr2==circles_12[i,"chr2"],"start2"])[which(sort(df[df$chr2==circles_12[i,"chr2"],"start2"])==df[n,"start2"])-1]) df[n-1,] else df[FALSE,]
        link<-link[link[,"direction2"]==circles_12[i,"direction2"],] 
        } else {link<-df[FALSE,]}
        link_up<-link 
        
        n<-which(df$start1==circles_12[i,"start1_b2"] &
                   df$chr1 == circles_12[i,"chr1"] &
                   df$chr2 == circles_12[i,"chr2"])#the block2 position in df 
        if(!n == dim(df)[1] & df[n,"start2"] != max(df[df$chr2==df[n,"chr2"],"start2"])) {
        link<-if(df[n+1,"start2"]==sort(df[df$chr2==circles_12[i,"chr2"],"start2"])[which(sort(df[df$chr2==circles_12[i,"chr2"],"start2"])==df[n,"start2"])+1]) df[n+1,] else df[FALSE,]
        link<-link[link[,"direction2"]==circles_12[i,"direction2"],] 
        } else {link<-df[FALSE,]}
        
        df[df$ID %in% unlist(c(link_up$ID,circles_12[i,c("ID_b1","ID_b2")],link$ID)),]
      })

    
    circles[["1-2-"]] <- circles_12_OK
    circles[["1-2-"]][sapply(circles[["1-2-"]], function(x) !dim(x)[1]>0)]<- NULL
    circles[["1-2-"]]<-if(!length(circles[["1-2-"]])>0) "No circles found" else circles[["1-2-"]]
    
  } else {circles[["1-2-"]]<-"No circles found"}#end 1-2- circle sequence
  
  #begin 2+1+ circle sequence
  circles_21_block1<-df[(which(grepl("^2\\+1\\+",df$synteny_dir))-1),]  # ignore the last column
  circles_21_block2<-df[which(grepl("^2\\+1\\+",df$synteny_dir)),]
  
  if ( dim(circles_21_block1)[1]>0 & dim(circles_21_block2)[1]>0 ) { #if there are no pairs found at all, skip sequence
    
    circles_21_block1[,"Number"]<-1:dim(circles_21_block1)[1]
    circles_21_block2[,"Number"]<-1:dim(circles_21_block2)[1]
    circles_21<-full_join(select(circles_21_block1, -"synteny_dir"),select(circles_21_block2, -starts_with(c("chr","direction"))),by="Number")
    colnames(circles_21)<-gsub("\\.x","_b1",colnames(circles_21))
    colnames(circles_21)<-gsub("\\.y","_b2",colnames(circles_21))
    duplicates<-(circles_21$ID_b1 %in% intersect(circles_21$ID_b1,circles_21$ID_b2)|circles_21$ID_b2 %in% intersect(circles_21$ID_b1,circles_21$ID_b2))
    duplicated_block<-rep(intersect(circles_21$ID_b1,circles_21$ID_b2),each=2)
    
    
    
    n=0
    circles_21[,"Duplicated_blocks"]<-sapply(seq_along(duplicates), function(i) {
      if (duplicates[i]) {
        n <<- n+1 #you can upate a variable outside of the apply environment by using <<
        paste(duplicated_block[n])}
      else {"Unique"}
    })
    
    circles[["2+1+"]] <- lapply(1:dim(circles_21)[1], function(i) { 
      df[df$chr2==circles_21[i,"chr2"] & df$start2 >= circles_21[i,"start2_b2"] & df$start2 <= circles_21[i,"start2_b1"],]})
    
    circles[["2+1+"]][sapply(circles[["2+1+"]], function(x) !dim(x)[1]>0)]<- NULL
    circles[["2+1+"]]<-if(!length(circles[["2+1+"]])>0) "No circles found" else circles[["2+1+"]]
    
  } else {circles[["2+1+"]]<-"No circles found"} #end 2+1+ sequence
  
  return(list(circles,table(synteny_dir)))
}

##Script for extraction of synteny blocks from Synteny Portal data, randomization and plotting

###################FUNCTIONS END - SCRIPT STARTS

setwd("L:/synteny-circle project/R")
library(tidyverse)
library(ggpubr)
library(data.table)
library(seqinr)
library(plyr)
library(ggpubr)
library(rlist)
options(stringsAsFactors=FALSE)
select<-dplyr::select
reduce<-purrr::reduce
#When you'd like to save the R environment:
#save.image(file = "210625_R_workspace_analysis_of_circules.RData")
#when you'd like to load it:
#load(file = ".Rdata")

HOMU<-readLines("SynBuilder/hg38_mm10.txt")
HOMA<-readLines("SynBuilder/hg38_rheMac3.txt")
MAMU<-readLines("SynBuilder/rheMac3_mm10.txt")
#added 29/7
HOBO<-readLines("SynBuilder/hg38_bosTau8.txt") #Bos taurus
HOCF<-readLines("SynBuilder/hg38_canFam3.txt") #Canis lupus familiaris
HODR<-readLines("SynBuilder/hg38_danRer10.txt") #Danio rerio
HORN<-readLines("SynBuilder/hg38_m5.txt") #Rattus norvegicus
HOMD<-readLines("SynBuilder/hg38_monDom5.txt") #Monodelphis domestica
HOOA<-readLines("SynBuilder/hg38_oviAri3.txt") #Ovis aries 
HOPT<-readLines("SynBuilder/hg38_panTro6.txt") #Pan troglodytes
#Following organisms are not assembled to chromomsome, so they were not analysed:
#HOSB<-readLines("SynBuilder/hg38_saiBol1.txt") #Saimiri boliviensis boliviensis
#HOTS<-readLines("SynBuilder/hg38_tarSyr2.txt") #tarsius syrichta
#HOXT<-readLines("SynBuilder/hg38_xenTro7.txt") #Xenopus tropicalis 

organisms<-c(reference="Homo sapiens",HOMU="Mus musculus",HOMA="Macaca mulatta",HOBO="Bos taurus",HOCF="Canis lupus familiaris",
             HODR="Danio rerio",HORN="Rattus norvegicus",HOMD="Monodelphis domestica",HOOA="Ovis aries",
             HOPT="Pan troglodytes") #HOSB="Saimiri boliviensis boliviensis",HOTS="Tarsius syrichta",HOXT="Xenopus tropicalis"
ensemble_name<-c("hsapiens","mmusculus","mmulatta","btaurus","clfamiliaris", "drerio","rnorvegicus","mdomestica","oaries","ptroglodytes")

assemly_to_org<-c(Hg38="Homo sapiens",mm10="Mus Musculus",rheMac3="Macaca mulatta",bosTau8="Bos taurus",canFam3="Canis lupus familiaris",
                  danRer10="Danio rerio",rn5="Rattus norvegicus",monDom5="Monodelphis domestica",oviAri3="Ovis aries",
                  panTro6="Pan troglodytes")

dfs<-list(
HOMU=as.dataframe.synbuilder(df=HOMU,gen1="hg38",gen2="mm10"),
HOMA=as.dataframe.synbuilder(df=HOMA,gen1="hg38",gen2="rheMac3"),
HOBO=as.dataframe.synbuilder(df=HOBO,gen1="hg38",gen2="bosTau8"),
HOCF=as.dataframe.synbuilder(df=HOCF,gen1="hg38",gen2="canFam3"),
HODR=as.dataframe.synbuilder(df=HODR,gen1="hg38",gen2="danRer10"),
HORN=as.dataframe.synbuilder(df=HORN,gen1="hg38",gen2="rn5"),
HOMD=as.dataframe.synbuilder(df=HOMD,gen1="hg38",gen2="monDom5"),
HOOA=as.dataframe.synbuilder(df=HOOA,gen1="hg38",gen2="oviAri3"),
HOPT=as.dataframe.synbuilder(df=HOPT,gen1="hg38",gen2="panTro6")
#HOSB=as.dataframe.synbuilder(df=HOSB,gen1="hg38",gen2="saiBol1"),
#HOTS=as.dataframe.synbuilder(df=HOTS,gen1="hg38",gen2="tarSyr2"),
#HOXT=as.dataframe.synbuilder(df=HOXT,gen1="hg38",gen2="xenTro7")
)

#export as csv to check if the find_circles function is removing wrong stuff
#for (x in names(dfs)) {
#  write.table(dfs[[x]],file=paste0(x,"_allblocks_Synbuilder.txt"),sep = "\t",quote=FALSE,col.names=T,row.names=F)
#}
rm(list=names(organisms[-1]))

noblocks<-sapply(dfs, dim)[1,] # number of synteny blocks
names(noblocks)<-organisms[-1]
#run seperately:
MAMU=as.dataframe.synbuilder(df=MAMU,gen1="rheMac3",gen2="mm10")

##Finding the circles
circles<-lapply(dfs, find_circles)
circles1<-circles

##gather all circle syntenies in one data frame and export
df_12<-lapply(seq_along(circles), 
       function(i) circles[[i]][[1]]$'1-2-' %>% reduce(rbind)) %>% reduce(rbind)
df_12$ID<-rep(1:(dim(df_12)[1]/2),each=2)
df_12<-df_12[order(df_12$chr1,df_12$ID, df_12$start1),]
df_12["org2"]<-gsub("\\..*","",df_12$ID2)

df_21<-lapply(seq_along(circles), 
              function(i) circles[[i]][[1]]$'2+1+' %>% reduce(rbind)) %>% reduce(rbind)
df_21$ID<-rep(1:(dim(df_21)[1]/2),each=2)
df_21<-df_21[order(df_21$chr1,df_21$ID, df_21$start1),]
df_21["org2"]<-gsub("\\..*","",df_21$ID2)

write.table(df_12,file="Circles12.txt",sep = "\t",quote=FALSE,col.names=T,row.names=F)
write.table(df_21,file="Circles21.txt",sep = "\t",quote=FALSE,col.names=T,row.names=F)

#BiocManager::install("chromPlot")
circles_coordinates_12<- data.frame(Chrom=aggregate(chr1~ID, unique, data=df_12)[,2],
           Start=aggregate(start1~ID, min, data=df_12)[,2],
           End=aggregate(stop1~ID, max, data=df_12)[,2],
           Chrom_org2=aggregate(chr2~ID, unique, data=df_12)[,2],
           Start_org2=aggregate(start2~ID, min, data=df_12)[,2],
           End_org2=aggregate(stop2~ID, max, data=df_12)[,2],
           Direction=aggregate(direction2~ID, unique, data=df_12)[,2],
           Size_Human_mb=round((aggregate(stop1~ID, max, data=df_12)[,2]-aggregate(start1~ID, min, data=df_12)[,2])/1000000,digits=2),
           Size_org2_mb=round((aggregate(stop2~ID, max, data=df_12)[,2]-aggregate(start2~ID, min, data=df_12)[,2])/1000000,digits=2),
           Assembly=aggregate(org2~ID, unique, data=df_12)[,2],
           Breakpoint_start_human=aggregate(stop1~ID, min, data=df_12)[,2],
           Breakpoint_end_human=aggregate(start1~ID, max, data=df_12)[,2],
           Breakpoint_size_human_kb=round((aggregate(start1~ID, max, data=df_12)[,2]-aggregate(stop1~ID, min, data=df_12)[,2])/1000,digits=2),
           Breakpoint_start_org2=aggregate(stop2~ID, min, data=df_12)[,2],
           Breakpoint_end_org2=aggregate(start2~ID, max, data=df_12)[,2],
           Breakpoint_size_org2_kb=round((aggregate(start2~ID, max, data=df_12)[,2]-aggregate(stop2~ID, min, data=df_12)[,2])/1000,digits=2))

#If blocks are overlapping, correct with 0
circles_coordinates_12[circles_coordinates_12$Breakpoint_size_org2_kb<0,"Breakpoint_size_org2_kb"]<-0

circles_coordinates_21<- data.frame(Chrom=aggregate(chr1~ID, unique, data=df_21)[,2],
                                    Start=aggregate(start1~ID, min, data=df_21)[,2],
                                    End=aggregate(stop1~ID, max, data=df_21)[,2],
                                    Chrom_org2=aggregate(chr2~ID, unique, data=df_21)[,2],
                                    Start_org2=aggregate(start2~ID, min, data=df_21)[,2],
                                    End_org2=aggregate(stop2~ID, max, data=df_21)[,2],
                                    Direction=aggregate(direction2~ID, unique, data=df_21)[,2],
                                    Size_Human_mb=round((aggregate(stop1~ID, max, data=df_21)[,2]-aggregate(start1~ID, min, data=df_21)[,2])/1000000,digits=2),
                                    Size_org2_mb=round((aggregate(stop2~ID, max, data=df_21)[,2]-aggregate(start2~ID, min, data=df_21)[,2])/1000000,digits=2),
                                    Assembly=aggregate(org2~ID, unique, data=df_21)[,2],
                                    Breakpoint_start_human=aggregate(stop1~ID, min, data=df_21)[,2],
                                    Breakpoint_end_human=aggregate(start1~ID, max, data=df_21)[,2],
                                    Breakpoint_size_human_kb=round((aggregate(start1~ID, max, data=df_21)[,2]-aggregate(stop1~ID, min, data=df_21)[,2])/1000,digits=2),
                                    Breakpoint_start_org2=aggregate(stop2~ID, min, data=df_21)[,2],
                                    Breakpoint_end_org2=aggregate(start2~ID, max, data=df_21)[,2],
                                    Breakpoint_size_org2_kb=round((aggregate(start2~ID, max, data=df_21)[,2]-aggregate(stop2~ID, min, data=df_21)[,2])/1000,digits=2))

circles_coordinates_21[circles_coordinates_21$Breakpoint_size_org2_kb<0,"Breakpoint_size_org2_kb"]<-0

circles_coordinates<-rbind(circles_coordinates_12,circles_coordinates_21)
circles_coordinates[,"Group"]<-assemly_to_org[match(circles_coordinates$Assembly,names(assemly_to_org))]
circles_coordinates<-circles_coordinates[order(circles_coordinates$Chrom,circles_coordinates$Start,circles_coordinates$Group),]

write.table(circles_coordinates,file="circles_coordinates_nosizeselection.txt",sep = "\t",quote=FALSE,col.names=T,row.names=F)
write.table(circles_coordinates,file="circles_coordinates_wodanio.txt",sep = "\t",quote=FALSE,col.names=T,row.names=F)


###PLOTTING THE ENTIRE CHROMOSOME

library(viridisLite)
library(chromPlot)
data(hg_gap) #load the gaps, including centromere, telomere...
#chromPlot(gaps=hg_gap[!hg_gap=="Y",], segment=circles_coordinates_12,noHist=TRUE,chr=c(1:22,"X"))
#chromPlot(gaps=hg_gap[!hg_gap=="Y",], segment=circles_coordinates_21,noHist=TRUE,chr=c(1:22,"X"))
Order_speciesname<-c(a="Pan_t",b="Macca_m",c="Mus_m",d="Rattus_n",e="Ovis_a",f="Bos_t",g="Canis_l",h="Monodelphis_d",i="Danio_r")
pal=c(a="black",b="#757d79",c="#0072b2",d="#56b4e9",e="#e69d00",f="#f0e442",g="#18e442",h="#cc79a7",i="#d55e00")


pdf(file = "export/211026_circlesynteny_plot_all_chr_all.pdf",
    7, 6) # 7 and 6 inches 
chromPlot(gaps=hg_gap[!hg_gap=="Y",], segment=circles_coordinates,noHist=TRUE,chr=c(1:22,"X"), 
          colSegments=pal[names(sort(Order_speciesname))])
dev.off()

?chromPlot()

pdf(file = "export/211026_circlesynteny_plot_all_chr_all_below150kb.pdf",
    7, 6) # 7 and 6 inches 
chromPlot(gaps=hg_gap[!hg_gap=="Y",], segment=circles_coordinates[circles_coordinates$Breakpoint_size_human_kb<=150,],noHist=TRUE,chr=c(1:22,"X"), 
          colSegments=pal[names(sort(Order_speciesname))])
dev.off()

#without zebrafish
Order_speciesname2<-Order_speciesname[!grepl("i",names(Order_speciesname))]
c(a="Pan_t",b="Macca_m",c="Mus_m",d="Rattus_n",e="Ovis_a",f="Bos_t",g="Canis_l",h="Monodelphis_d")
pal2=c(a="#0072b2",b="#56b4e9",c="#229439",d="#18e442",e="#e69d00",f="#f0e442",g="#d55e00",h="black")


pdf(file = "export/211026_circlesynteny_plot_all_chr_wodanio_rerio.pdf",
    4.5, 10) # 7 and 6 inches 
chromPlot(gaps=hg_gap[!hg_gap=="Y",], segment=circles_coordinates[circles_coordinates$Group!="Danio rerio",],noHist=TRUE,chr=c(1:22,"X"), 
          colSegments=pal2[names(sort(Order_speciesname2))])
dev.off()



###BEST RECIPROCAL HIT TO TEST CIRCLES AND PLOTTING ALL THE INDIVIDUAL CIRCLES
##Preparing for plotting

annotations<-sapply(ensemble_name, function(x) {
  link<-read.table(file=paste0("biomart_maps/", x, ".tsv"),sep = "\t", header = T)
  link<-link[order(link$chr,link$start),]
}, simplify = F, USE.NAMES = T)


##plot the insertion patterns
source("circle_plot_functions.R")

for (i in 1:length(circles)) {
  plot_circles12(circles=circles[[i]][[1]][["1-2-"]], def1=annotations[[1]],def2=annotations[[i+1]],
                 outputname=paste0("export/211123_10kbcutoff/",names(organisms[i+1]), "_1-2-potential_circles_SynBuilder.pdf"), 
                 org1=organisms[1],
                 org2=organisms[i+1],size_cutoff=10) #50 kb cutoff for plotting
}

for (i in 1:length(circles)) {
  plot_circles21(circles=circles[[i]][[1]][["2+1+"]], def1=annotations[[1]],def2=annotations[[i+1]],
                 outputname=paste0("export/211123_10kbcutoff/",names(organisms[i+1]), "_2+1+potential_circles_SynBuilder.pdf"), 
                 org1=organisms[1],
                 org2=organisms[i+1],size_cutoff=10) #50 kb cutoff for plotting
}


###TO CONFIRM THE SYNTENY, HOMOLOGY OF ALL GENES WERE INVESTIGATED WITH MMSEQ2 (command line)
###FIRST STEP. LOAD THE BEST RECIRPOCAL HITS GENERATED WITH MMSEQ2

##ADD ONLY GENES WITH BEST RECIPROCAL BLAST HITS IN THE AREA OF THE
##CIRCLE PATTERN

#fastafile<-gzfile("../mmseq/proteomes/Homo_sapiens.GRCh38.pep.all.fa.gz")
#header<-grep(">",fastafile, value=T)
#close(fastafile)

dfs$HOOA[dfs$HOOA$chr1==15 & dfs$HOOA$start1>=(82385406-1000) & dfs$HOOA$stop1<=(85170013+1000),]

#What genes are inside the circles?
#The file has been manually annotated with "id" and is loaded again
insertion_patterns<-read.table(file="Circle_insertion_patterns_backtoR.csv", sep=";",header = T,dec = ",")

#get the genes into a list for each circle pattern
genes_list<-list()
for (n in seq_along(ensemble_name[-1])) {
link0<-annotations[[ensemble_name[-1][n]]]
genes_list[[ensemble_name[-1][n]]]<-
  lapply(1:dim(insertion_patterns)[1], function(i) {
    print(n)
  if(insertion_patterns[i,"Assembly"]==names(assemly_to_org[-1][n])) {
  link<-link0[link0$chr==insertion_patterns[i,"Chrom_org2"] &
              ((link0$start>=insertion_patterns[i,"Start_org2"] &
              link0$end<=insertion_patterns[i,"Breakpoint_start_org2"]) |
              (link0$start>=insertion_patterns[i,"Breakpoint_end_org2"] &
              link0$end<=insertion_patterns[i,"End_org2"])
                  )  ,]
  if(dim(link)[1]>0) link[order(link$start),] else NULL
  }
})
names(genes_list[[ensemble_name[-1][n]]])<-insertion_patterns$Pattern_id
}
#remove NULLs
genes_list<-lapply(genes_list, function(x) x[!sapply(x, is.null)])

#what is in insertion pattern 22?
insertion_patterns[insertion_patterns$Pattern_id=="22",]
lapply(genes_list, function(x) x["22"])
cow22<-genes_list$btaurus[["22"]]
cow22[cow22$size_kb>=100,]

#Insertion pattern 66 is the most difficult. Check what genes are inside
#note: some genes span blocks and gaps. Fx Bos taurus has ELF1 in block and gap and was plotted anyways to show the homology to human
sheep66<-genes_list$mdomestica[["66"]]
sheep66[sheep66$size_kb>=100,]
#Check if there are genes in block and gap in rat pattern 66:
rat66<-insertion_patterns[insertion_patterns$Pattern_id=="66" & insertion_patterns$Organism2 == "Rattus norvegicus",]
annotations$rnorvegicus[annotations$rnorvegicus$start>=rat66$Start_org2 & annotations$rnorvegicus$end<=rat66$End_org2 & annotations$rnorvegicus$chr==rat66$Chrom_org2,] %>% .[,"symbol"]

#Check pattern 74
cow74<-genes_list$btaurus[["74"]]

##Plot missing circle insertion pattern 66 in BT
plot_missing <-list()
#bos taurus
plot_missing[[1]]<-dfs$HOBO[dfs$HOBO$chr1==15 & dfs$HOBO$start1>=82536244 & dfs$HOBO$stop1<=85169524,] 
plot_missing[[2]]<-dfs$HOBO[dfs$HOBO$chr2==21 & dfs$HOBO$start2>=(22698129) & dfs$HOBO$stop2<=25519818,]

insertion_patterns[insertion_patterns$Pattern_id=="66",]

blocks_cow66<-dfs$HOBO[dfs$HOBO$chr1==15 & dfs$HOBO$start1>=82536244 & dfs$HOBO$stop1<=85169524,8:10]
blocks_human66<-dfs$HOBO[dfs$HOBO$chr1==15 & dfs$HOBO$start1>=82536244 & dfs$HOBO$stop1<=85169524,4:6]

lapply(1:dim(blocks_cow)[1], function(i) { 
  annotations$btaurus[annotations$btaurus$chr==blocks_cow[i,1] & annotations$btaurus$start>=blocks_cow[i,2] & annotations$btaurus$start<=blocks_cow[i,3],]
})
lapply(1:dim(blocks_human)[1], function(i) { 
  annotations$hsapiens[annotations$hsapiens$chr==blocks_human[i,1] & annotations$hsapiens$start>=blocks_human[i,2] & annotations$hsapiens$start<=blocks_human[i,3],]
})

##Where is human block 4 in cow?
#look for genes
where_is_sec11<-annotations$btaurus[annotations$btaurus$wikigene_name=="SEC11A",]
dfs$HOBO[which(dfs$HOBO$chr2==where_is_sec11$chr & dfs$HOBO$start2<where_is_sec11$start &dfs$HOBO$stop2>where_is_sec11$start),]

##Block 74 dog, cow, sheep
pattern74<-insertion_patterns[insertion_patterns$Pattern_id=="74",]
human_genes_pattern74<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern74$Chrom) & 
                       ((annotations$hsapiens$start>=min(pattern74$Start) & 
                       annotations$hsapiens$end<=max(pattern74$Breakpoint_start)) | 
                       (annotations$hsapiens$start>=min(pattern74$Breakpoint_end) & 
                          annotations$hsapiens$end<=max(pattern74$End))),]
human_genes_pattern74_block1<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern74$Chrom) & 
                                                     annotations$hsapiens$start>=min(pattern74$Start) & 
                                                         annotations$hsapiens$end<=max(pattern74$Breakpoint_start),]

#look up pubmed hits per intersecting gene
require(rentrez)

#in block1
intersecting_genes_pattern74_block1<-Reduce(intersect, list(genes_list$clfamiliaris$`74`$symbol,genes_list$btaurus$`74`$symbol,genes_list$oaries$`74`$symbol,human_genes_pattern74_block1$symbol))
#there are no intersecting genes with all three organisms and human in block1:
#check individually:
for (y in c("clfamiliaris","btaurus","oaries")) {
print(y)
  link<-intersect(genes_list[[y]]$`74`$symbol,human_genes_pattern74_block1$symbol)
  if(length(link>0)) {
    pubmed_hits<-lapply(link, function(x) entrez_search(db="pubmed",term=paste0(x,"[All]")))
  pubmed_counts<-sapply(pubmed_hits, function(x) x$count)
  names(pubmed_counts)<-link
  print(pubmed_counts %>% sort(., decreasing = T))
  }else print("no intersect")
  
}
#->NCOR1,TRPV2,UBB,TTC19,PIGL

#check the gene order for sheep
annotations$oaries[annotations$oaries$symbol %in%
                     genes_list$oaries$`74`$symbol[!genes_list$oaries$`74`$symbol==""],] %>% .[order(.["start"]),]

#intersecting genes with all three organisms and human
intersecting_genes_pattern74<-Reduce(intersect, list(genes_list$clfamiliaris$`74`$symbol,genes_list$btaurus$`74`$symbol,genes_list$oaries$`74`$symbol,human_genes_pattern74$symbol)) %>% .[!.==""]
pubmed_hits<-lapply(intersecting_genes_pattern74, function(x) entrez_search(db="pubmed",term=paste0(x,"[All]")))
pubmed_counts<-sapply(pubmed_hits, function(x) x$count)
names(pubmed_counts)<-intersecting_genes_pattern74
pubmed_counts %>% sort(., decreasing = T)
#->TOP3A  SREBF1    FLCN  ALKBH5    FLII ALDH3A1   MAPK7    RAI1   SHMT1    ULK2 ALDH3A2

annotations$hsapiens[annotations$hsapiens$symbol %in% intersecting_genes_pattern74[!intersecting_genes_pattern74==""],] %>% .[order(.["start"]),]

#Check insertion pattern 77:
pattern77<-insertion_patterns[insertion_patterns$Pattern_id=="77",]
human_genes_pattern77<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern77$Chrom) & 
                                              ((annotations$hsapiens$start>=min(pattern77$Start) & 
                                                  annotations$hsapiens$end<=max(pattern77$Breakpoint_start)) | 
                                                 (annotations$hsapiens$start>=min(pattern77$Breakpoint_end) & 
                                                    annotations$hsapiens$end<=max(pattern77$End))),]
human_genes_pattern77_block1<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern77$Chrom) & 
                                                     annotations$hsapiens$start>=min(pattern77$Start) & 
                                                     annotations$hsapiens$end<=max(pattern77$Breakpoint_start),]

#all intersecting genes
intersecting_genes_pattern77<-Reduce(intersect, list(genes_list$btaurus$`77`$symbol,genes_list$oaries$`77`$symbol,human_genes_pattern77$symbol)) %>% .[!.==""]
#intersecting genes in block1
intersecting_genes_pattern77_block1<-Reduce(intersect, list(genes_list$btaurus$`77`$symbol,genes_list$oaries$`77`$symbol,human_genes_pattern77_block1$symbol)) %>% .[!.==""]
#in block1, sheep
annotations$oaries[(annotations$oaries$symbol %in% human_genes_pattern77_block1$symbol),]
#sort by pubmed hits for all genes, smaller than 10kb
human_genes_pattern77<-annotations$hsapiens[(annotations$hsapiens$symbol %in% intersecting_genes_pattern77) & annotations$hsapiens$size_kb>=10,]$symbol
pubmed_hits<-lapply(intersecting_genes_pattern77, function(x) entrez_search(db="pubmed",term=paste0(x,"[All]")))
pubmed_counts<-sapply(pubmed_hits, function(x) x$count)
names(pubmed_counts)<-intersecting_genes_pattern77
pubmed_counts %>% sort(., decreasing = T)
pubmed_counts<-data.frame(symbol=names(pubmed_counts),count=pubmed_counts)
pubmed_counts<-merge(annotations$hsapiens[annotations$hsapiens$symbol %in% pubmed_counts$symbol,],pubmed_counts,by="symbol")
pubmed_counts[order(pubmed_counts$count, decreasing = T),] %>% .[.["size_kb"]>=10,c("symbol","start","size_kb")]

##check insertion pattern 64
#all intersect
pattern64<-insertion_patterns[insertion_patterns$Pattern_id=="64",]
human_genes_pattern64<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern64$Chrom) & 
                                              ((annotations$hsapiens$start>=min(pattern64$Start) & 
                                                  annotations$hsapiens$end<=max(pattern64$Breakpoint_start)) | 
                                                 (annotations$hsapiens$start>=min(pattern64$Breakpoint_end) & 
                                                    annotations$hsapiens$end<=max(pattern64$End))),]

intersecting_genes_pattern64<-Reduce(intersect, list(#genes_list$btaurus$`64`$symbol,
                                                     genes_list$oaries$`64`$symbol,
                                                     genes_list$mdomestica$`64`$symbol,
                                                     toupper(genes_list$mmusculus$`64`$symbol),
                                                     toupper(genes_list$rnorvegicus$`64`$symbol),
                                                     human_genes_pattern64$symbol)) %>% .[!.==""]
sapply(genes_list, function(x) x$`64`)

#only three genes are found in all->"FBXO22"  "TMEM266" "ISL2"
#cow and sheep:
pattern64_cow_and_sheep<-intersect(genes_list$btaurus$`64`$symbol,
          genes_list$oaries$`64`$symbol) %>% .[!.==""]
#mouse and rat:
pattern64_mouse_and_rat<-intersect(toupper(genes_list$mmusculus$`64`$symbol),
                                   toupper(genes_list$rnorvegicus$`64`$symbol)) %>% .[!.==""]

pubmed_hits<-lapply(pattern64_mouse_and_rat, function(x) entrez_search(db="pubmed",term=paste0(x,"[All]")))
pubmed_counts<-sapply(pubmed_hits, function(x) x$count)
names(pubmed_counts)<-pattern64_mouse_and_rat
pubmed_counts %>% sort(., decreasing = T)
pubmed_counts<-data.frame(symbol=names(pubmed_counts),count=pubmed_counts)
pubmed_counts<-merge(annotations$hsapiens[annotations$hsapiens$symbol %in% pubmed_counts$symbol,],pubmed_counts,by="symbol")
pubmed_counts[order(pubmed_counts$count, decreasing = T),] %>% .[,c("symbol","start","size_kb")]



#what happend with dog?
dog64<-dfs$HOCF[dfs$HOCF$chr1==unique(pattern64$Chrom) & 
           dfs$HOCF$start1>=(min(pattern64$Start)-5000000) &
           dfs$HOCF$stop1<=(max(pattern64$End)),]

#->blocks containing genes from dog are now in different chromosomes


#where are the genes in block2 in the cow genome?
where_is_DNAJA4<-annotations$btaurus[annotations$btaurus$wikigene_name=="DNAJA4",]
where_is_DNAJA4<-dfs$HOBO[which(dfs$HOBO$chr2==where_is_DNAJA4$chr & dfs$HOBO$start2<where_is_DNAJA4$start &dfs$HOBO$stop2>where_is_DNAJA4$start),]
cow64<-insertion_patterns[insertion_patterns$Pattern_id=="64" & insertion_patterns$Organism2=="Bos taurus",]
cow64_extended_upstream<-dfs$HOBO[dfs$HOBO$chr2==cow64$Chrom_org2 &
                                    dfs$HOBO$start2>=where_is_DNAJA4$start2 &
                                    dfs$HOBO$stop2<=cow64$End_org2,]
plot_missing[[3]]<-cow64_extended_upstream

lapply(1:dim(cow64_extended_upstream)[1], function(i) annotations$btaurus[annotations$btaurus$chr==cow64_extended_upstream[i,"chr2"] &
                                                                          annotations$btaurus$start>=cow64_extended_upstream[i,"start2"] &
                                                                            annotations$btaurus$end<=cow64_extended_upstream[i,"stop2"], "symbol"] )
#some of the genes are in the block upstream
#human block 2
human64_block2<-annotations$hsapiens[annotations$hsapiens$chr>=cow64_extended_upstream[3,"chr1"] &
  annotations$hsapiens$start>=cow64_extended_upstream[3,"start1"] &
    annotations$hsapiens$end<=cow64_extended_upstream[3,"stop1"],"symbol"] 

cow64_blockup<-annotations$btaurus[annotations$btaurus$chr>=cow64_extended_upstream[1,"chr2"] &
                                       annotations$btaurus$start>=cow64_extended_upstream[1,"start2"] &
                                       annotations$btaurus$end<=cow64_extended_upstream[1,"stop2"],"symbol"] %>% .[!.==""]
intersect(human64_block2,cow64_blockup)

##check insertion pattern 75
#all intersect
pattern75<-insertion_patterns[insertion_patterns$Pattern_id=="75",]
human_genes_pattern75<-annotations$hsapiens[annotations$hsapiens$chr==unique(pattern75$Chrom) & 
                                              ((annotations$hsapiens$start>=min(pattern75$Start) & 
                                                  annotations$hsapiens$end<=max(pattern75$Breakpoint_start)) | 
                                                 (annotations$hsapiens$start>=min(pattern75$Breakpoint_end) & 
                                                    annotations$hsapiens$end<=max(pattern75$End))),]
plot_missing[[4]]<-dfs$HORN[dfs$HORN$chr1==unique(pattern75$Chrom) &
        dfs$HORN$start1>=min(pattern75$Start) &
        dfs$HORN$stop1<=max(pattern75$End),]

#where is NF1 in rat?
annotations$rnorvegicus[grepl("Nf1",annotations$rnorvegicus$symbol),]
#in the linked area, linked to another human block
#take the human coordinates and find all rat blocks for plotting
plot_missing[[5]]<-dfs$HORN[dfs$HORN$chr2==unique(plot_missing[[4]]$chr2) &
                              dfs$HORN$start2>=min(plot_missing[[4]]$start2) &
                              dfs$HORN$stop2<=max(plot_missing[[4]]$stop2),]
plot_missing[[6]]<-plot_missing[[5]][-dim(plot_missing[[5]])[1],]

plot_missing[[7]]<-dfs$HORN[dfs$HORN$chr1==unique(plot_missing[[5]]$chr1) &
                              dfs$HORN$start1>=min(plot_missing[[5]]$start1) &
                              dfs$HORN$stop1<=max(plot_missing[[5]]$stop1),]


#plot missing rn circles
plot_circles21(circles=plot_missing[4:7],
               def1=annotations$hsapiens,def2=annotations$rnorvegicus,
               outputname=paste0("export/","Rattus_norvigicus_degenerated_insertion_pattern75", "_2+1+potential_circles_SynBuilder.pdf"), 
               org1=organisms["reference"],
               org2=organisms["HORN"],size_cutoff=0) #50 kb cutoff for plotting



#plot missing bt circles
plot_circles21(circles=plot_missing[3],
               def1=annotations$hsapiens,def2=annotations$btaurus,
               outputname=paste0("export/","Bos_taurus_degenerated_insertion_pattern64", "_2+1+potential_circles_SynBuilder.pdf"), 
               org1=organisms["reference"],
               org2=organisms["HOBO"],size_cutoff=0) #50 kb cutoff for plotting
plot_circles12(circles=plot_missing[3],
               def1=annotations$hsapiens,def2=annotations$btaurus,
               outputname=paste0("export","Bos_taurus_degenerated_insertion_pattern64", "_1-2-potential_circles_SynBuilder.pdf"), 
               org1=organisms["reference"],
               org2=organisms["HOBO"],size_cutoff=0) #50 kb cutoff for plotting


#Check the probability of getting circles randomly, compared to 
##Number of circles found in the native syntenies
##Finding the circles
no_circles<-lapply(circles,function(x) lapply(x[[1]], length))

#write.table(circles_21,file="circles_21_Synchro.txt",sep = "\t",quote=FALSE,col.names=T,row.names=F)

#how big a portion of the ungulate genomes could have been derived from circles
insertion_patterns$Assembly %>% unique %>% sort


genomesize<-c("bosTau8"=2670028162,"oviAri3"=2619037772,"canFam3"=2410976875,
              "mm10"=2730855475, "monDom5"=3598443077,"panTro6"=3050398082,"rheMac3"=2969988180,"rn5"=	2909698938)

sapply(names(genomesize), function(x) {
  link<-insertion_patterns[insertion_patterns$Assembly==x,"Size_org2_mb"] %>% sum
  paste0(format(round(100*link/(genomesize[x]/(10^6)),digits = 1),nsmall=1),"%")
}, USE.NAMES = T)



###RANDOMIZED SYNTENY BLOCKS
######Completely randomized org2 genome

random_blocks<-lapply(1:length(dfs), function(i) {
  link<-dfs[[i]]
  lapply(1:10000, function(dummy) {
    #randomize everything with a "2" in the column name
    cbind(link %>% select(-contains("2")),
          link %>% select(contains("2")) %>% .[sample(1:dim(link)[1], dim(link)[1], replace=F),])
  })
})

##Call the circles
random_test<-lapply(1:length(dfs), function(i) {
  print(names(dfs[i]))
  lapply(1:10000, function(n) {
    find_circles(random_blocks[[i]][[n]])
  })
})
names(random_test)<-names(dfs)

#for single organism
#random_test_MAMU<-lapply(1:10000, function(i) {
#  print(i)
#  find_circles(random_blocks_MAMU[[i]])
#})

#list the number of circle synteny pairs found for 1-2- and 2+1+
number_of_circles <- list()
number_of_circles[["1-2-"]]<-sapply(names(random_test),function(x) {
  sapply(seq_along(random_test[[x]]), function(i) random_test[[x]][[i]][[1]][["1-2-"]][!(random_test[[x]][[i]][[1]][["1-2-"]]=="No circles found")] %>% length) }, USE.NAMES = T, simplify = F)
number_of_circles[["2+1+"]]<-sapply(names(random_test),function(x) {
  sapply(seq_along(random_test[[x]]), function(i) random_test[[x]][[i]][[1]][["2+1+"]][!(random_test[[x]][[i]][[1]][["2+1+"]]=="No circles found")] %>% length) }, USE.NAMES = T, simplify = F)

stats<-list()
stats[["1-2-"]]<-sapply(names(number_of_circles[["1-2-"]]), function(x) {
  link<-table(Number=number_of_circles[["1-2-"]][[x]]) %>% data.frame %>% transform(., Number=Number %>% as.character %>% as.numeric)
  link<-link %>% mutate( ToHighlight = ifelse(link$Number >=  no_circles[[x]][["1-2-"]], "yes", "no" ))
  link}, simplify = F)
stats[["2+1+"]]<-sapply(names(number_of_circles[["2+1+"]]), function(x) {
  link<-table(Number=number_of_circles[["2+1+"]][[x]]) %>% data.frame %>% transform(., Number=Number %>% as.character %>% as.numeric)
  link<-link %>% mutate( ToHighlight = ifelse(link$Number >= no_circles[[x]][["2+1+"]] , "yes", "no" ))
  link}, simplify = F)

maxno<-sapply(names(stats), function(x) sapply(names(stats[[x]]), function(y) stats[[x]][[y]][["Freq"]] %>% max )) %>% max
maxsyn<-sapply( names(stats), function(x) sapply(names(stats[[x]]), function(y) stats[[x]][[y]][["Number"]] %>% max))

plot_hist<-sapply(names(stats), function(n1) sapply(names(stats[[n1]]), function(n2) {
  ggplot(stats[[n1]][[n2]])+aes(x=Number, y=Freq, fill=ToHighlight) + geom_bar(stat = "identity") +
    scale_x_continuous(breaks = 0:(maxsyn[n2,]%>% max),limits = c(-0.5,(maxsyn[n2,]%>% max)+0.5))  +xlab("") + ylab ("") + ylim(c(0,maxno))+ 
    scale_fill_manual(values = c( "yes"="black", "no"="white" ), guide = FALSE ) +
    theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))
}, simplify = F, USE.NAMES = T), USE.NAMES = T, simplify = F)

ggarrange(plotlist=c(plot_hist[["1-2-"]],plot_hist[["2+1+"]])[c(rbind(1:length(plot_hist[[1]]),(length(plot_hist[[1]])+1):(length(plot_hist[[1]])*2)))],
          labels=paste(names(plot_hist), rep(organisms[-1],each=2),", unscrambled =", no_circles %>% unlist),
          font.label = list(size = 8, color = "red"),vjust=-0.2, hjust=-0.1,
          ncol = length(plot_hist),nrow = length(plot_hist[[1]]))
# Annotate the figure by adding a common labels
annotate_figure(random_hist_c,left = text_grob("Occurence of circle synteny", color = "black", rot = 90)) +
  theme(plot.margin = margin(1,0.5,1,0.5, "cm")) #,
#top = text_grob("Potential circles found in 10,000 randomized synteny maps", color = "black", face = "bold", size = 14))



