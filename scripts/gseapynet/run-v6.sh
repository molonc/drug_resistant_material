

#### This was for the January version of Farhia's thesis, with 4 (or 5) columns
#COMPS=("pathway9" "ressens8" "pathway4Tvs4U_2" "ressens18" "fitness56")
#TITLES=("SA609 Rx X7:cloneA vs. UnRx X7:cloneH" "SA1035 Rx X8cloneH vs. UnRx X8:cloneE" "SA535:cisplatin Rx X9:cloneST vs. UnRx X9:cloneJ" "SA535:CX5461 Rx X8:cloneU vs. UnRx X8:cloneJ" "SA535:CX5461:Rx X6:cloneR vs. SA535:UnRx X6:cloneJ")
#DATASETS=("SA609" "SA1035" "SA535_cisplatin" "SA535_CX5461" "SA535_CX5461")


# this is for pathway evolution
## COMPS=("comp10" "comp11" "comp12" "comp13" "comp14" "comp15" "comp16" "comp17")

# this is for the other pathway
#COMPS=("pathway1" "pathway5" "pathway9" "pathway2" "pathway4Tvs4U_1" "pathway6" "pathway12" "pathway3" "pathway4Tvs4U_2" "pathway4" "pathway7" "pathway8")
#TITLES=("SA609 UTTTT vs. UT all clones" "SA609 UTTTT R vs. UT R" "SA609 UTTTT R vs. UUUUU H" "SA1035 UTTTT vs. UT all clones" "SA1035 UTTTT G,H vs. UUUUU E" "SA1035 UTTTT G,H vs. UT G,H" "SA535_CX5461 UXXXX U vs. UUUUU Q" "SA535_cisplatin UUTTTT vs. UUT all clones" "SA535_cisplatin UUTTTT S,T vs. UUUUUU J,Q" "SA535_CX5461 UXXXX vs. UX all clones" "SA535_cisplatin UUTTTT S,T vs. UUT S,T" "SA535_CX5461 UXXXX U vs. UX R")

#COMPS=("pathway12" "pathway4Tvs4U_2")
#TITLES=("SA535_CX5461 UXXXX U vs. UUUUU Q" "SA535_cisplatin UUTTTT S,T vs. UUUUUU J,Q")

#COMPS=("pathway9" "ressens7" "ressens8" "pathway4Tvs4U_2" "ressens18")
#TITLES=("SA609 UTTTT R vs. UUUUU H" "SA1035 UTTTT G vs. UUUUU E" "SA1035 UTTTT H vs. UUUUU E" "SA535 cisplatin UUTTTT S_T vs. UUUUUU J" "SA535 CX5461 UXXXX U vs. UUUUU J")
#DATASETS=("SA609" "SA1035" "SA535_cisplatin" "SA535_CX5461")



### Started 13 Feb 2021, de pathways at different time points, combines clones
## SA609
COMPS=("ressenstime1" "ressenstime2" "ressenstime3" "ressenstime4")
TITLES=("SA609 Rx X4:cloneAB vs. UnRx X4:cloneCH" "SA609 Rx X5:cloneAB vs. UnRx X5:cloneCH" "SA609 Rx X6:cloneAB vs. UnRx X6:cloneCH" "SA609 Rx X7:cloneAB vs. UnRx X7:cloneCH")
DATASETS=("SA609" "SA609" "SA609" "SA609")

# SA1035
COMPS=("ressenstime6" "ressenstime7" "ressenstime8")
TITLES=("SA1035 Rx X6:cloneGH vs. UnRx X6:cloneDE" "SA1035 Rx X7:cloneGH vs. UnRx X7:cloneDE" "SA1035 Rx X8:cloneGH vs. UnRx X8:cloneDE" )
DATASETS=("SA1035" "SA1035" "SA1035")

# SA535_cisplatin
COMPS=("ressenstime9" "ressenstime11" "ressenstime12")
TITLES=("SA535:cisplatin Rx X6:cloneRST vs. UnRx X6:cloneJQ" "SA535:cisplatin Rx X8:cloneRST vs. UnRx X8:cloneJQ" "SA535:cisplatin Rx X9:cloneRST vs. UnRx X9:cloneJQ" )
DATASETS=("SA535_cisplatin" "SA535_cisplatin" "SA535_cisplatin")

# SA535_CX5461
COMPS=("ressenstime13" "ressenstime14" "ressenstime15")
TITLES=("SA535:CX5461 Rx X5:cloneRU vs. UnRx X5:cloneJQ" "SA535:CX5461 Rx X6:cloneRU vs. UnRx X6:cloneJQ" "SA535:CX5461 Rx X8:cloneRU vs. UnRx X8:cloneJQ" )
DATASETS=("SA535_CX5461" "SA535_CX5461" "SA535_CX5461")


####################
####################
### Started 19 Feb 2021, de pathways at different time points, single clones
## SA609
COMPS=("ressenstime21" "ressenstime22" "ressenstime23" "ressenstime24")
TITLES=("SA609 Rx X4:cloneA vs. UnRx X4:cloneH" "SA609 Rx X5:cloneA vs. UnRx X5:cloneH" "SA609 Rx X6:cloneA vs. UnRx X6:cloneH" "SA609 Res(X7:A) vs. Sen(X7:H)")
DATASETS=("SA609" "SA609" "SA609" "SA609")

# SA1035
COMPS=("ressenstime27" "ressenstime28")
TITLES=("SA1035 Rx X7:cloneH vs. UnRx X7:cloneE" "SA1035 Res(X8:H) vs. Sen(X8:E)" )
DATASETS=("SA1035" "SA1035")

# SA535_cisplatin
COMPS=("ressenstime29" "ressenstime31" "ressenstime32" "ressenstime33")
TITLES=("SA535:Cisplatin Res(X6:R) vs. Sen(X6:J)" "SA535:Cisplatin Rx X8:cloneST vs. UnRx X8:cloneJ" "SA535:Cisplatin Res(X9:ST) vs. Sen(X9:J)" "SA535:Cisplatin Res(X10:ST) vs. Sen(X9:J)")
DATASETS=("SA535_cisplatin" "SA535_cisplatin" "SA535_cisplatin" "SA535_cisplatin")

# SA535_CX5461
COMPS=("ressenstime34" "ressenstime35" "ressenstime36" "ressenstime37" "ressenstime38" "ressenstime39" "ressenstime40")
TITLES=("SA535:CX5461 Rx X5:cloneR vs. UnRx X5:cloneQ" "SA535:CX5461 Res(X6:R) vs. Sen(X6:J)" "SA535:CX5461 Rx X6:cloneR vs. UnRx X6:cloneQ" "SA535:CX5461 Res(X8:R) vs. Sen(X8:J)" "SA535:CX5461 Rx X8:cloneR vs. UnRx X8:cloneQ" "SA535:CX5461 Res(X8:U) vs. Sen(X8:J)" "SA535:CX5461 Rx X8:cloneU vs. UnRx X8:cloneQ")
DATASETS=("SA535_CX5461" "SA535_CX5461" "SA535_CX5461" "SA535_CX5461" "SA535_CX5461" "SA535_CX5461" "SA535_CX5461")


##########################
### The networks for the figures in thesis
COMPS=("ressenstime24" "ressenstime28" "ressenstime29" "ressenstime32" "ressenstime35" "ressenstime39" "ressenstime37")
TITLES=("SA609 Res(X7:A) vs. Sen(X7:H)" "SA1035 Res(X8:H) vs. Sen(X8:E)" "SA535:Cisplatin Res(X6:R) vs. Sen(X6:J)" "SA535:Cisplatin Res(X9:ST) vs. Sen(X9:J)" "SA535:CX5461 Res(X6:R) vs. Sen(X6:J)" "SA535:CX5461 Res(X8:U) vs. Sen(X8:J)" "SA535:CX5461 Res(X8:R) vs. Sen(X8:J)")
DATASETS=("SA609" "SA1035" "SA535_cisplatin" "SA535_cisplatin" "SA535_CX5461" "SA535_CX5461" "SA535_CX5461")
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "4.5")
FIGH=("3" "3" "3" "3" "3" "3" "4")


##########################
### 18 May 2021: networks for resistant vs sensitive, all samples
COMPS=("allsamples1" "allsamples2" "allsamples3" "allsamples4" "allsamples5" "allsamples6" "allsamples7")
TITLES=("SA609 Res(all:A) vs. Sen(all:CDGH)" "SA609 Res(all:B) vs. Sen(all:CDGH)" "SA609 Res(all:A) vs. Sen(all:B)" "SA1035 Res(all:H) vs. Sen(all:E)" "SA535:Cisplatin Res(all:ST) vs. Sen(all:JQ)" "SA535:CX5461 Res(all:U) vs. Sen(all:JQ)" "SA535 Res(all:R) vs. Sen(all:JQ)")
DATASETS=("SA609" "SA609" "SA609" "SA1035" "SA535_cisplatin" "SA535_CX5461" "SA535_CX5461")
FIGW=("3.25" "3.25" "3.25" "3.25"  "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3")


##########################
### 7 Jun 2021: resistant vs sensitive for all samples combined
## Now including ver into the data sets
COMPS=("ihc1" "ihc3" "ihc5")
TITLES=("SA609 Res(all:A) vs. Sen(all:H)" "SA1035 Res(all:H) vs. Sen(all:E)" "SA535 Res(all:A) vs. Sen(all:G)")
DATASETS=("SA609-v6" "SA1035-v6" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")

##########################
### 9 Jun 2021: resistant vs sensitive by time 
## Now including ver into the data sets
COMPS=("rstmm1" "rstmm2" "rstmm3" "rstmm4" "rstmm5" "rstmm6" "rstmm7" "rstmm8" "rstmm9" "rstmm10" "rstmm11")
TITLES=("SA609 Res(X4:A) vs. Sen(X4:H)" "SA609 Res(X5:A) vs. Sen(X5:H)" "SA609 Res(X6:A) vs. Sen(X6:H)" "SA609 Res(X7:A) vs. Sen(X7:H)" "SA1035 Res(X7:H) vs. Sen(X7:E)" "SA1035 Res(X8:H) vs. Sen(X8:E)" "SA535 Res(X6:AG) vs. Sen(X6:G)" "SA535 Res(X7:AEJ) vs. Sen(X7:G)" "SA535 Res(X8:A) vs. Sen(X8:G)" "SA535 Res(X9:A) vs. Sen(X9:G)" "SA535 Res(X10:A) vs. Sen(X9:G)")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA1035-v6" "SA1035-v6" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3")


## clone unaware SA609
COMPS=("tutmm01" "tutmm02" "tutmm03" "tutmm04" "uutmm01" "uutmm02" "uutmm03" "uutmm04" "hutmm01" "hutmm02" "hutmm03")
TITLES=("SA609 Rx X4 vs. Un X3" "SA609 Rx X5 vs. Un X3" "SA609 Rx X6 vs. Un X3" "SA609 Rx X7 vs. Un X3" "SA609 Un X4 vs. Un X3" "SA609 Un X5 vs. Un X3" "SA609 Un X6 vs. Un X3" "SA609 Un X7 vs. Un X3" "SA609 Ho X5 vs. Un X3" "SA609 Ho X6 vs. Un X3" "SA609 Ho X7 vs. Un X3")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" )
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3")


## clone unaware SA1035
COMPS=("tutmm10" "tutmm11" "tutmm12" "tutmm13" "uutmm10" "uutmm11" "uutmm12" "uutmm13")
TITLES=("SA1035 Rx X5 vs. Un X4" "SA1035 Rx X6 vs. Un X4" "SA1035 Rx X7 vs. Un X4" "SA1035 Rx X8 vs. Un X4" "SA1035 Un X5 vs. Un X4" "SA1035 Un X6 vs. Un X4" "SA1035 Un X7 vs. Un X4" "SA1035 Un X8 vs. Un X4")
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6" )
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3")

## clone unaware SA535
COMPS=("tutmm21" "tutmm22" "tutmm23" "tutmm24" "tutmm25" "uutmm21" "uutmm22" "uutmm23" "uutmm24" "hutmm21" "hutmm22" "hutmm23" "hutmm24")
TITLES=("SA535 Rx X6 vs. Un X5" "SA535 Rx X7 vs. Un X5" "SA535 Rx X8 vs. Un X5" "SA535 Rx X9 vs. Un X5" "SA535 Rx X10 vs. Un X5" "SA535 Un X6 vs. Un X5" "SA535 Un X7 vs. Un X5"  "SA535 Un X8 vs. Un X5"  "SA535 Un X9 vs. Un X5" "SA535 Ho X7 vs. Un X5" "SA535 Ho X8 vs. Un X5"  "SA535 Ho X9 vs. Un X5"  "SA535 Ho X10 vs. Un X5")
DATASETS=("SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" )
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3")

### The new DE with scran normalization before edgeR
## clone unaware SA535
COMPS=("scrande_SA609_1" "scrande_SA609_2" "scrande_SA609_3" "scrande_SA609_4" "scrande_SA1035_1" "scrande_SA1035_2" "scrande_SA535_1" "scrande_SA535_2" "scrande_SA535_3")
TITLES=("SA609 Res(X4:A) vs. Sen(X4:H)" "SA609 Res(X5:A) vs. Sen(X5:H)" "SA609 Res(X6:A) vs. Sen(X6:H)" "SA609 Res(X7:A) vs. Sen(X7:H)" "SA1035 Res(X7:H) vs. Sen(X7:E)" "SA1035 Res(X8:H) vs. Sen(X8:E)" "SA535 Res(X8:A) vs. Sen(X8:G)" "SA535 Res(X9:A) vs. Sen(X9:G)" "SA535 Res(X10:A) vs. Sen(X9:G)")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA1035-v6" "SA1035-v6" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3")

### No significant pathways for scrande_SA609_20


### New scran edgeR with  
COMPS=("scrande_SA609_20")  # "scrande_SA609_22" "scrande_SA535_4" "scrande_SA535_5")
TITLES=("SA609 Res(X5:A) vs. Hol(X5:B)")  # "SA609 Hol(X5:B) vs. Sen(X5:H)" "SA535 Res(X6:AG) vs. Sen(X6:G)" "SA535 Res(X7:AEJ) vs. Sen(X7:G)")
DATASETS=("SA609-v6")  # "SA609-v6" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25")  # "3.25" "3.25" "3.25")
FIGH=("3")  # "3" "3" "3")

### New scran edgeR with  
COMPS=("scrande_SA609_23" "scrande_SA609_24")  # "scrande_SA609_22" "scrande_SA535_4" "scrande_SA535_5")
TITLES=("SA609 Hol(X6:R) vs. Hol(X6:H)" "SA609 Hol(X7:R) vs. Hol(X7:H)")  # "SA609 Hol(X5:B) vs. Sen(X5:H)" "SA535 Res(X6:AG) vs. Sen(X6:G)" "SA535 Res(X7:AEJ) vs. Sen(X7:G)")
DATASETS=("SA609-v6" "SA609-v6")  # "SA609-v6" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25")  # "3.25" "3.25" "3.25")
FIGH=("3" "3")  # "3" "3" "3")


COMPS=("scran-nocl-SA609-1" "scran-nocl-SA609-2" "scran-nocl-SA609-3" "scran-nocl-SA609-4"  "scran-nocl-SA609-12" "scran-nocl-SA609-13" "scran-nocl-SA609-14" "scran-nocl-SA609-22" "scran-nocl-SA609-23" "scran-nocl-SA609-24")
TITLES=("SA609 Rx X4 vs. UnRx X4" "SA609 Rx X5 vs. UnRx X5" "SA609 Rx X6 vs. UnRx X6" "SA609 Rx X7 vs. UnRx X7" "SA609 Rx X5 vs. RxH X5" "SA609 Rx X6 vs. RxH X6" "SA609 Rx X7 vs. RxH X7" "SA609 RxH X5 vs. UnRx X5" "SA609 RxH X6 vs. UnRx X6" "SA609 RxH X7 vs. UnRx X7")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6" "SA609-v6")
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3" "3")


COMPS=("scran-nocl-SA535-1" "scran-nocl-SA535-2" "scran-nocl-SA535-3" "scran-nocl-SA535-4" "scran-nocl-SA535-5" "scran-nocl-SA535-12" "scran-nocl-SA535-13" "scran-nocl-SA535-14" "scran-nocl-SA535-15" "scran-nocl-SA535-22" "scran-nocl-SA535-23" "scran-nocl-SA535-24" "scran-nocl-SA535-25")
TITLES=("SA535 Rx X6 vs. UnRx X6" "SA535 Rx X7 vs. UnRx X7" "SA535 Rx X8 vs. UnRx X8" "SA535 Rx X9 vs. UnRx X9" "SA535 Rx X10 vs. UnRx X9" "SA535 Rx X7 vs. RxH X7" "SA535 Rx X8 vs. RxH X8" "SA535 Rx X9 vs. RxH X9" "SA535 Rx X10 vs. RxH X10" "SA535 RxH X7 vs. UnRx X7" "SA535 RxH X8 vs. UnRx X8" "SA535 RxH X9 vs. UnRx X9" "SA535 RxH X10 vs. UnRx X9")
DATASETS=("SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7" )
FIGW=("3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3" "3")


COMPS=("scran-nocl-SA609-22" "scran-nocl-SA609-23" "scran-nocl-SA609-24")
TITLES=("SA609 RxH X5 vs. UnRx X5" "SA609 RxH X6 vs. UnRx X6" "SA609 RxH X7 vs. UnRx X7")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")
## make this one a bit larger
#COMPS=()
#TITLES=("SA535:CX5461 Res(X8:R) vs. Sen(X8:J)")
#DATASETS=("SA535_CX5461")


COMPS=("scran-nocl-SA1035-1" "scran-nocl-SA1035-2" "scran-nocl-SA1035-3" "scran-nocl-SA1035-4")
TITLES=("SA1035 Rx X5 vs. UnRx X5" "SA1035 Rx X6 vs. UnRx X6" "SA1035 Rx X7 vs. UnRx X7" "SA1035 Rx X8 vs. UnRx X8" )
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25" "3.25" "3.25")
FIGH=("3" "3" "3" "3")


COMPS=("scran-nocl-SA1035-22" "scran-nocl-SA1035-23" "scran-nocl-SA1035-24")
TITLES=("SA1035 RxH X6 vs. UnRx X6" "SA1035 RxH X7 vs. UnRx X7" "SA1035 RxH X8 vs. UnRx X8" )
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


COMPS=("scran-nocl-SA1035-12" "scran-nocl-SA1035-13" "scran-nocl-SA1035-14")
TITLES=("SA1035 Rx X6 vs. RxH X6" "SA1035 Rx X7 vs. RxH X7" "SA1035 Rx X8 vs. RxH X8" )
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


COMPS=("scrande_SA535_52" "scrande_SA535_16" "scrande_SA535_36")
TITLES=("Pt5 Hol(X7:AEJ) vs. Sen(X7:G)" "Pt5 Hol(X8:D) vs. Sen(X8:G)" "Pt5 Res(X10:E) vs. Sen(X9:G)" )
DATASETS=("SA535_cisplatin-v7" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


COMPS=("scrande_SA609_20" "scrande_SA609_25" "scrande_SA609_26")
TITLES=("SA609 X5 Rx:A vs. RxH:B" "SA609 X6 Rx:A vs. RxH:A" "SA609 X7 Rx:A vs. RxH:A")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")
## make this one a bit large


COMPS=("scrande_SA609_31" "scrande_SA535_76" "scrande_SA535_78")
TITLES=("Pt4 Hol(X5:A) vs. Sen(X5:H)" "Pt5 Hol(X8:A) vs. Sen(X8:G)" "Pt5 Hol(X10:A) vs. Sen(X10:G)")
DATASETS=("SA609-v6" "SA535_cisplatin-v7" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


COMPS=("scrande_SA609_11")
TITLES=("Pt4 Res(X4:B) vs. Sen(X4:H)")
DATASETS=("SA609-v6")
FIGW=("3.25")
FIGH=("3")


COMPS=("scrande_SA609_22" "scrande_SA609_23" "scrande_SA609_24")
TITLES=("Pt4 Res(X5:B) vs. Sen(X5:H)" "Pt4 Res(X6:A) vs. Sen(X6:H)" "Pt4 Res(X7:A) vs. Sen(X7:H)")
DATASETS=("SA609-v6" "SA609-v6" "SA609-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")

COMPS=("scrande_SA609_31" "scrande_SA535_78")
TITLES=("Pt4 Hol(X5:A) vs. Sen(X5:H)" "Pt5 Hol(X10:A) vs. Sen(X10:G)")
DATASETS=("SA609-v6" "SA535_cisplatin-v7")
FIGW=("3.25" "3.25")
FIGH=("3" "3")

COMPS=("scrande_SA535_16")
TITLES=("Pt5 Hol(X8:D) vs. Sen(X8:G)")
DATASETS=("SA535_cisplatin-v7")
FIGW=("3.25")
FIGH=("3")

COMPS=("scrande_SA535_78")
TITLES=("Pt5 Hol(X10:A) vs. Sen(X10:G)")
DATASETS=("SA535_cisplatin-v7")
FIGW=("3.25")
FIGH=("3")

COMPS=("scrande_SA535_36")
TITLES=("Pt5 Res(X10:E) vs. Sen(X9:G)" )
DATASETS=("SA535_cisplatin-v7")
FIGW=("3.25")
FIGH=("3")


COMPS=("scrande_SA1035_41" "scrande_SA1035_51" "scrande_SA1035_52")
TITLES=("Pt6 Res(X5:D) vs. Sen(X5:A)" "Pt6 Res(X7:G) vs. Sen(X7:E)" "Pt6 Res(X7:G) vs. Sen(X7:E)")
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")

COMPS=("scrande_SA1035_36" "scrande_SA1035_34")
TITLES=("Pt6 Hol(X7:H) vs. Sen(X7:E)" "Pt6 Hol(X8:H) vs. Sen(X8:E)")
DATASETS=("SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25")
FIGH=("3" "3")


COMPS=("scrande_SA1035_3")
TITLES=("Pt6 Res(X6:G) vs. Sen(X6:G)")
DATASETS=("SA1035-v6")
FIGW=("3.25")
FIGH=("3")


COMPS=("scrande_SA1035_32" "scrande_SA1035_12" "scrande_SA1035_22")
TITLES=("Pt6 Hol(X6:B) vs. Sen(X6:G)" "Pt6 Hol(X7:G) vs. Sen(X7:E)" "Pt6 Hol(X8:G) vs. Sen(X8:E)" )
DATASETS=("SA1035-v6" "SA1035-v6" "SA1035-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


COMPS=("scrande_SA1035_41" "scrande_SA609_22" "scrande_SA609_31")
TITLES=("Pt6 Res(X5:B) vs. Sen(X5:A)" "Pt4 Res(X5:B) vs. Sen(X5:H)" "Pt4 Hol(X5:A) vs. Sen(X5:H)")
DATASETS=("SA1035-v6" "SA609-v6" "SA609-v6")
FIGW=("3.25" "3.25" "3.25")
FIGH=("3" "3" "3")


SET="HALLMARK"
#SET="KEGG"
VER="v6"

#STHR=0.0001
STHR=0.01
#STHR=0.01
#STHR=0.05
#SCOL="PValue"
SCOL="FDR"

## Setting on MINFC 0.25 and PFDR 0.05
#MINFC="0.0"
MINFC="0.25"
#MINFC="0.5"

#PFDR="0.5"
PFDR="0.05"
#PFDR="0.01"

## NOTE Feb 12 2021: it was run with PValue as sig-col, not FDR


if [ "$SET" == "HALLMARK" ]; then
    GMTFILE="h.all.v7.0.symbols.gmt"
else     
    GMTFILE="KEGG_2019_Human"
fi

i=0
for COMP in ${COMPS[@]}; do
  echo $COMP
  echo $i
  echo ${TITLES[$i]}
  for DATASET in ${DATASETS[@]}; do
    # echo "ls $DATASET-$VER/${COMP}_*_logfc_results.csv"
    echo "ls $DATASET/${COMP}_*_logfc_results.csv"
    ls $DATASET/${COMP}_*_logfc_results.csv
    if ls $DATASET/${COMP}_*_logfc_results.csv 1> /dev/null 2>&1; then
    	echo $DATASET
    	python3 gseapynet_Viki.py \
  	  --csv $DATASET/${COMP}_*_logfc_results.csv \
  	  --sig $STHR \
  	  --sig-col $SCOL \
  	  --minfc $MINFC \
  	  --minfc-col logFC \
  	  --pathway-fdr $PFDR \
  	  --gene-col gene_symbol \
  	  --png $DATASET/${COMP}-${DATASET}-${SET}-${SCOL}-${STHR}-MINFC-${MINFC}-PFDR-${PFDR}-pathway.pdf \
  	  --gmt $GMTFILE \
  	  --title "${TITLES[$i]}" \
  	  --figw "${FIGW[$i]}" \
  	  --figh "${FIGH[$i]}"
  	else 
  	  echo "$DATASET does not exist"
  	fi  
  done
  i=$((i+1))
done

##exit

