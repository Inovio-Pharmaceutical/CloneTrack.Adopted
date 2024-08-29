# ID clones with 2 fold increase per CloneTrack paper

##Read in periphery data for CloneTrack Analysis
## Expected column headers: 
##    SID - Subject Identifier
##    VISIT - Sample Visit Name
##    Clones - Clone Count
##    CDR3.nt - Nucleic acid sequence of CDR3 region
##    CDR3.aa - Amino Acid Sequence of CDR3 region

library(tidyverse)

pbmc.source.data.df <- 
  bind_rows(readRDS(
    "CloneTrack/RRP001.TcrBeta.PBMC.Master.21Feb2024.RDS"), .id = "SID_VISIT") %>% 
  mutate(
    SID_VISIT = str_remove(SID_VISIT, "-"), 
    SID = str_remove(SID_VISIT, "_.*"), 
    VISIT = str_remove(SID_VISIT, ".*_")
  ) %>% 
  relocate(SID, VISIT) %>% 
  filter(V.MEAN >= 70, J.MEAN >= 70, # Vendor suggested QC criteria
         !str_detect(CDR3.aa, "_|\\*") # Select only productive TCR sequences
  ) %>%
  select(SID, VISIT, Clones, CDR3.nt, CDR3.aa) 

##Calculate total read count for each sample
##    lib.sizes - total read count for each sample
lib.sizes <- 
  pbmc.source.data.df %>% 
  group_by(SID, VISIT) %>% 
  summarise(lib.size = sum(Clones))

#Prep dataframe for analysis

fold.master.df <- 
  
  ## Start with baseline visits(Day0 or Screening)
  left_join(by = c("SID", "CDR3.aa", "CDR3.nt"),
            left_join(
              by = c("SID", "VISIT"),
              pbmc.source.data.df %>% 
                filter(str_detect(VISIT, "Day0|Screening")) %>% 
                rename(Base.Clones = Clones), 
              lib.sizes
            ) %>% 
              rename(Base.Lib = lib.size) %>% 
              select(SID, Base.Clones, Base.Lib, CDR3.aa, CDR3.nt),
            ## Join non-baseline visit data to dataframe
            left_join(
              pbmc.source.data.df %>% 
                filter(!str_detect(VISIT, "Day0|Screening")) %>% 
                rename(VISIT.Clones = Clones), 
              lib.sizes
            ) %>% 
              rename(VISIT.Lib = lib.size) %>% 
              select(SID, VISIT, VISIT.Clones, VISIT.Lib, CDR3.aa, CDR3.nt)
            
  ) %>% 
  mutate(
    # (T/F) is visit ratio 2x more than Baseline read
    First.Cut.2Fold = (2*(Base.Clones/Base.Lib)) < (VISIT.Clones/VISIT.Lib),
    # Actual visit/Base fold change
    First.Cut.value = (VISIT.Clones/VISIT.Lib)/(Base.Clones/Base.Lib)) %>% 
  filter(First.Cut.2Fold == "TRUE") %>% 
  mutate(
    Base.Lib.Fish = (Base.Lib/2) - Base.Clones, #Baseline ratio for fisher test - adopted from orig. CloneTrack 
    Visit.Lib.Fish = (VISIT.Lib - VISIT.Clones),# Visit ration for fisher test 
    Target = paste(SID, VISIT, CDR3.nt, sep = "_") #used to align SID, VISIT, CDR3.nt for testing
  )

#Custom function to calculate fisher exact pvalue - adopted from CloneTrack

my.fish.fun <- 
  function(temp){
    rbind(
      c(temp["Base.Clones"] %>% unlist,temp["VISIT.Clones"]%>% unlist), 
      c(temp["Base.Lib.Fish"]%>% unlist,temp["Visit.Lib.Fish"]%>% unlist)
    ) %>% unname %>% 
      stats::fisher.test() %>% 
      .$p.value
  }

#Prep for function loop application
targets <- fold.master.df$Target %>% unique


#empty df for fisher results
output <- 
  data.frame(Target = character, 
             P.Val = numeric)

#calc fisher exact pvalues - adopted from CloneTrack
for(t in 1: length(targets)){
  output <- 
    rbind(output, 
          cbind(
            Target = targets[t],
            P.Val = fold.master.df %>% 
              filter(Target == targets[t]) %>% 
              my.fish.fun()
          ))
}

#Merge pvalues with original data
fold.output <- 
  full_join(by = "Target",fold.master.df,output) %>% 
  select(-Target)


# Calculate adjusted pvalues - adopted from CloneTrack
## Padj = "bonferroni" reported.pvalue * number of test for Subject-visit combination

test.count <- 
  left_join(by = c("SID", "CDR3.aa", "CDR3.nt"),
            pbmc.master %>% 
              filter(str_detect(VISIT, "Day0|Screening")) %>% 
              select(-VISIT),
            
            pbmc.master %>% 
              filter(!str_detect(VISIT, "Day0|Screening"))
  ) %>%  
  na.omit %>% 
  group_by(SID, VISIT) %>% 
  tally() %>% 
  rename(TCR.Overlap = n)

fold.output.adj <- 
  left_join(by = c("SID", "VISIT"), fold.output, test.count) %>% 
  mutate(P.Val = as.numeric(P.Val),
         P.ADJ = P.Val * TCR.Overlap) %>% 
  filter(P.ADJ <= 0.05) # value adopted from CloneTrack

# filter CloneTrack Significant TCR sequences seen in EOS tissue samples
##  EOS - Week47 or Week52 tissue samples
##  Read in tissue TCR data
##     Expected column headers: 
##        SID - Subject Identifier
##        VISIT - Sample Visit Name
##        Clones - Clone Count
##        CDR3.nt - Nucleic acid sequence of CDR3 region
##        CDR3.aa - Amino Acid Sequence of CDR3 region

tissue.master <- 
  bind_rows(readRDS("RRP001.TcrBeta.Master.Jan2024.RDS"), .id = "SID_VISIT") %>% 
  mutate(
    SID_VISIT = str_remove(SID_VISIT, "-"), 
    SID = str_remove(SID_VISIT, "_.*"), 
    VISIT = str_remove(SID_VISIT, ".*_")
  ) %>% 
  relocate(SID, VISIT) %>% 
  filter(
    !str_detect(CDR3.aa, "_|\\*"),  #Select only productive TCR sequences
    V.MEAN >= 70, J.MEAN >= 70 # Vendor suggested QC criteria
  ) %>% 
  filter(str_detect(VISIT, "Week52|Week47")) %>% #Select EOS visits only
  select(SID, Clones, CDR3.nt, CDR3.aa) %>% 
  rename(Tissue.EOS.Clones = Clones)


eos.tissue.CloneTrack.pbmc.df<- 
  left_join(by = c("SID", "CDR3.nt", "CDR3.aa"),
            fold.output.adj, 
            tissue.master) %>% 
  na.omit #remove records that do not have a match



