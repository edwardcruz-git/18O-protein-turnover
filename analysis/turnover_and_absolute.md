General relationship between turnover, dynamicity, and abundance
================
Edward Cruz
2026-05-13

``` r
rm(list=ls(all=T))
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(fgsea)
library(stringr)
library(ggplot2)
library(Biostrings)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: generics

    ## 
    ## Attaching package: 'generics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     explain

    ## The following objects are masked from 'package:base':
    ## 
    ##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
    ##     setequal, union

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     combine

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
    ##     unsplit, which.max, which.min

    ## Loading required package: S4Vectors

    ## Loading required package: stats4

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     first, rename

    ## The following object is masked from 'package:utils':
    ## 
    ##     findMatches

    ## The following objects are masked from 'package:base':
    ## 
    ##     expand.grid, I, unname

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     collapse, desc, slice

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: XVector

    ## Loading required package: Seqinfo

    ## 
    ## Attaching package: 'Biostrings'

    ## The following object is masked from 'package:base':
    ## 
    ##     strsplit

``` r
library(Peptides)
library(lsa)
```

    ## Loading required package: SnowballC

``` r
library(patchwork)

knitr::opts_chunk$set(fig.path = "figures/turnover_and_absolute/")
```

``` r
#FASTA file used to search data
frog_fasta <- readAAStringSet("Files/FASTA/XENLA_longest-transcript_10p1_GFY-Fix.fasta")

#Converting to a dataframe for merging
frog_sizes <- data.frame(
  Protein_ID = names(frog_fasta),
  Sequence = as.character(frog_fasta),
  stringsAsFactors = FALSE )
row.names(frog_sizes) <- NULL
frog_sizes["Protein_ID"] <- sapply(strsplit(frog_sizes$`Protein_ID`, split = " "),
                                   "[", 1)

#FASTA file used to search data
fly_fasta <- readAAStringSet("Files/Reference/Dmel_Uniprot_Proteome_092722.fasta")

#Converting to a dataframe for merging
fly_sizes <- data.frame(
  Protein_ID = names(fly_fasta),
  Sequence = as.character(fly_fasta),
  stringsAsFactors = FALSE )
row.names(fly_sizes) <- NULL
fly_sizes["Protein_ID"] <- sapply(strsplit(fly_sizes$`Protein_ID`, split = " "),
                                  "[", 1)
fly_sizes["Protein_ID"] <- sapply(strsplit(fly_sizes$Protein_ID, split = "\\|"),
                                  "[", 2)
#---- Calculating MW ----------

frog_sizes["Sequence"] <- str_remove_all(frog_sizes$Sequence,
                                         "(?i)[^ARNDCEQGHILKMFPSTWYV]")
fly_sizes["Sequence"] <- str_remove_all(fly_sizes$Sequence,
                                        "(?i)[^ARNDCEQGHILKMFPSTWYV]")

frog_sizes["kDa"] <- sapply(frog_sizes$Sequence, function(seq){
    return(mw(seq)/1000) })
fly_sizes["kDa"] <- sapply(fly_sizes$Sequence, function(seq){
    return(mw(seq)/1000) })

frog_sizes <- frog_sizes %>% select(Protein_ID, kDa)
fly_sizes <- fly_sizes %>% select(Protein_ID, kDa)
```

``` r
XLA_T0_Abs <- read.csv("Files/AbsC/XLA_RT0_Abs-Concentrations.csv")
XLA_T0_Abs <- XLA_T0_Abs[-c(3:4)]

Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405) - 161
Frog_T5_light.min <- Frog_T5_all.min[c(1:2, 4:5, 7:8, 10:11)]

Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394) - 150
Frog_T6_light.min <- Frog_T6_all.min[c(1:2, 4:5, 7:8, 10:11)]

avg_light_time <- (Frog_T5_light.min + Frog_T6_light.min)/2
frog_hpf_light <- avg_light_time + ((150+161)/2)

avg_heavy_time <- (Frog_T5_all.min + Frog_T6_all.min)/2
frog_hpf_heavy <- avg_heavy_time + ((150+161)/2)

#------------------------------------------

XLA_fits_df <- read.csv("Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv")
head(XLA_fits_df)
```

    ##                   Protein_ID   XLA_Gene Human_Gene         kD    HL_Hrs
    ## 1  XBmRNA10569|XBXL10_1g5638    ccnb1.S      CCNB1 0.05817981 0.1985646
    ## 2 XBmRNA35507|XBXL10_1g19214  ccnb1.2.L      CCNB1 0.05236107 0.2206306
    ## 3 XBmRNA34126|XBXL10_1g18541  LOC398156      PTTG1 0.05055962 0.2284917
    ## 4 XBmRNA40199|XBXL10_1g21682  ccnb1.2.S      CCNB1 0.04860139 0.2376980
    ## 5   XBmRNA1720|XBXL10_1g1100    ccna2.L      CCNA2 0.03471313 0.3327977
    ## 6 XBmRNA20588|XBXL10_1g10943    ccna1.S      CCNA1 0.03404127 0.3393661
    ##   Estimate      T5_kD      T6_kD deltaBIC      T0_L        N1        N3
    ## 1 interval 0.05196905 0.06439057 59.78006 1.7147274 1.0893752 1.1090663
    ## 2    point         NA 0.05236107 19.02400 0.5715232 0.3167322 0.4122321
    ## 3    point         NA 0.05055962 53.52402 0.7593631 0.2741591 0.3558655
    ## 4    point 0.04860139         NA 16.69609 0.8122257 0.5816428 0.5083152
    ## 5    point         NA 0.03471313 50.29648 0.3719585 0.2490067 0.3498724
    ## 6 interval 0.02838204 0.03970050 51.62111 2.1104195 1.5523841 1.7135697
    ##          N4        N7         N9        N11          N12     T0_H        O1
    ## 1 1.2411665 1.8280703 0.92944429 0.06468336 2.346668e-02 9.479349 0.6087783
    ## 2 0.4270137 1.8856945 1.67675998 1.13389467 1.576150e+00 6.415712 0.3157785
    ## 3 0.4436930 1.9537817 1.95445907 0.69835989 1.560319e+00 9.303156 0.5306161
    ## 4 0.3789010 1.0024334 1.50784145 1.53953923 1.669101e+00 6.061506 0.5321921
    ## 5 0.3315833 1.4536714 1.20658237 1.33815842 2.699167e+00 8.538092 1.2239487
    ## 6 1.7205897 0.8069156 0.07777375 0.01829269 5.503508e-05 8.399209 1.3292753
    ##           O2        O3        O4         O5        O7         O9       O10
    ## 1 0.00000000 0.2046160 0.1690556 0.15700874 0.1380140 0.06531466 0.1457002
    ## 2 0.05760827 0.4069183 0.6042256 0.44871175 0.5156761 0.49982006 0.2410755
    ## 3 0.00000000 0.1234161 0.2365997 0.18184335 0.1569108 0.00000000 0.0000000
    ## 4 0.20740221 0.4029011 0.3367262 0.26181593 0.7716691 0.50081417 0.4388957
    ## 5 0.00000000 0.2604006 0.2245385 0.15593400 0.2258017 0.17048824 0.1304815
    ## 6 0.24473731 0.3145507 0.2879517 0.08938614 0.1435368 0.00000000 0.0000000
    ##          O11       O12 mClass Number_Pep_Fits
    ## 1 0.03216316 0.0000000    Deg               2
    ## 2 0.74512284 0.7493515    Deg               1
    ## 3 0.18012809 0.2873295    Deg               1
    ## 4 0.48312743 1.0029505    Deg               1
    ## 5 0.07031518 0.0000000    Deg               1
    ## 6 0.11168100 0.0796721    Deg               3

``` r
frog_dyn_df <- XLA_fits_df[c("Protein_ID", "Human_Gene", "kD",
                             "HL_Hrs", "mClass")]

frog_dyn_df["Dyn"] <- apply(XLA_fits_df[10:17], 1, function(x){

  flat_line <- rep(1, length(x))
  cosine_distance <- 1-cosine(flat_line, x)[1,]

  return(cosine_distance) })

frog_dyn_df <- merge(frog_dyn_df, XLA_T0_Abs, by.x="Protein_ID", by.y="Ind_Protein_ID")
frog_dyn_df <- merge(frog_dyn_df, frog_sizes, by="Protein_ID")
```

``` r
frog_dyn_df["Plot_iBAQ"] <- log10(frog_dyn_df$iBAQ_Conc*1000)
frog_dyn_df["Plot_Dyn"] <- log10(frog_dyn_df$Dyn)
frog_dyn_df["Plot_HL"] <- log10(frog_dyn_df$HL_Hrs)
frog_dyn_df["Plot_kD"] <- log10(abs(frog_dyn_df$kD))
frog_dyn_df <- frog_dyn_df %>% arrange(desc(abs(kD)))

head(frog_dyn_df)
```

    ##                   Protein_ID Human_Gene         kD    HL_Hrs mClass        Dyn
    ## 1  XBmRNA10569|XBXL10_1g5638      CCNB1 0.05817981 0.1985646    Deg 0.15135278
    ## 2 XBmRNA35507|XBXL10_1g19214      CCNB1 0.05236107 0.2206306    Deg 0.14394443
    ## 3 XBmRNA34126|XBXL10_1g18541      PTTG1 0.05055962 0.2284917    Deg 0.16720048
    ## 4 XBmRNA40199|XBXL10_1g21682      CCNB1 0.04860139 0.2376980    Deg 0.09807862
    ## 5   XBmRNA1720|XBXL10_1g1100      CCNA2 0.03471313 0.3327977    Deg 0.21787494
    ## 6 XBmRNA20588|XBXL10_1g10943      CCNA1 0.03404127 0.3393661    Deg 0.22785139
    ##      iBAQ_Conc      kDa Plot_iBAQ   Plot_Dyn    Plot_HL   Plot_kD
    ## 1 4.388700e-07 43.87020 -3.357664 -0.8200096 -0.7020981 -1.235228
    ## 2 1.925095e-06 44.62239 -2.715548 -0.8418051 -0.6563343 -1.280991
    ## 3 6.530381e-06 20.65477 -2.185061 -0.7767625 -0.6411296 -1.296196
    ## 4 2.428176e-06 44.67334 -2.614720 -1.0084256 -0.6239745 -1.313351
    ## 5 2.069420e-07 46.90004 -3.684151 -0.6617927 -0.4778196 -1.459506
    ## 6 5.024590e-06 46.81479 -2.298899 -0.6423483 -0.4693315 -1.467994

``` r
Dmel_T0_Abs <- read.csv("Files/AbsC/Dmel_T0_Abs-Concentrations.csv")
Dmel_T0_Abs <- Dmel_T0_Abs[-c(3:4)]

#------------------------------------------

#Ultimately - collected these minutes apart so used one vector
Fly_T6.time.min <- c(257, 287, 320, 351, 386, 441, 500,
                     565, 682, 868, 995, 1218, 1514) - 257
light_time <- c(0, Fly_T6.time.min[c(4, 8, 11, 13)])

Fly_T6_hpf_time <- ((Fly_T6.time.min + 257) + (Fly_T6.time.min + 117))/2
Fly_T6_hpf_time <- Fly_T6_hpf_time / 60
Fly_T6_hpf_light_time <- c(Fly_T6_hpf_time[1], Fly_T6_hpf_time[c(4, 8, 11, 13)])

#------------------------------------------

fly_blast_fits <-
  read.csv("Files/Fits/Dyn-Model/Fly_O18_T6-7_FinalFits_ALL.csv")
head(fly_blast_fits)
```

    ##   Protein_ID Gene_Symbol         kD    HL_Hrs Estimate      T6_kD      T7_kD
    ## 1     Q9VFD5      CG6966 0.01539471 0.7504170    point 0.01539471         NA
    ## 2     Q9VMA3         cup 0.01486216 0.7773063 interval 0.01538270 0.01434163
    ## 3     Q8IPM1         srl 0.01408818 0.8200105    point         NA 0.01408818
    ## 4     Q9VED4        bard 0.01264941 0.9132802    point         NA 0.01264941
    ## 5     P09085         cad 0.01217044 0.9492226    point         NA 0.01217044
    ## 6     Q9VNG1      CG2182 0.01191980 0.9691817 interval 0.01182441 0.01201520
    ##    deltaBIC     T0_L       N13       N14       N15       N16     T0_H       O1
    ## 1  7.832369 1.430443 0.3653358 0.6719770 1.3041889 1.2280553 3.365396 2.688617
    ## 2 10.631890 3.149509 0.6416273 0.4164807 0.3952220 0.3971613 3.778549 2.663961
    ## 3 11.704248 1.706303 0.5332616 1.4637543 0.9769972 0.3196836 3.141349 2.952435
    ## 4 23.447389 2.822755 1.7116034 0.2114468 0.1230994 0.1310951 4.380005 2.508002
    ## 5 19.866379 2.215896 1.0709636 0.5444190 0.5127122 0.6560089 3.627582 2.951559
    ## 6 18.057844 1.141686 1.0725223 1.0495665 0.9212109 0.8150147 3.075938 2.791074
    ##         O2        O3        O4        O5        O6        O7        O8
    ## 1 1.144322 0.8776167 0.7173222 0.4833720 0.6023776 0.4890450 0.5963080
    ## 2 1.149432 0.9644864 0.6122611 0.5266244 0.4758346 0.5489199 0.4854475
    ## 3 1.405408 1.1108608 0.6885891 0.7144997 0.3689344 0.5798764 0.5538892
    ## 4 1.952674 1.7729541 0.7960269 0.4390828 0.2705379 0.2474848 0.1318265
    ## 5 2.063225 1.0804270 0.5648673 0.5653617 0.3289622 0.3461073 0.3201948
    ## 6 1.871093 1.2094890 0.7646547 0.6387289 0.4852980 0.4144466 0.4046215
    ##          O9       O10       O11       O12 mClass
    ## 1 0.5719438 0.5362716 0.3889168 0.5384916    Deg
    ## 2 0.4306688 0.4159514 0.4550532 0.4928102    Deg
    ## 3 0.3652615 0.4694481 0.2848861 0.3645625    Deg
    ## 4 0.1570708 0.1219197 0.1053457 0.1170699    Deg
    ## 5 0.2565727 0.2858838 0.3429187 0.2663379    Deg
    ## 6 0.3697731 0.3911873 0.2936905 0.2900057    Deg

``` r
fly_dyn_df <- fly_blast_fits[c("Protein_ID", "Gene_Symbol", "kD",
                               "HL_Hrs", "mClass")]

fly_dyn_df["Dyn"] <- apply(fly_blast_fits[9:13], 1, function(x){

  flat_line <- rep(1, length(x))
  cosine_distance <- 1-cosine(flat_line, x)[1,]

  return(cosine_distance) })

fly_dyn_df <- merge(fly_dyn_df, Dmel_T0_Abs, by.x="Protein_ID", by.y="Ind_Protein_ID")
fly_dyn_df <- merge(fly_dyn_df, fly_sizes, by="Protein_ID")
```

``` r
fly_dyn_df["Plot_iBAQ"] <- log10(fly_dyn_df$iBAQ_Conc*1000)
fly_dyn_df["Plot_Dyn"] <- log10(fly_dyn_df$Dyn)
fly_dyn_df["Plot_HL"] <- log10(fly_dyn_df$HL_Hrs)
fly_dyn_df["Plot_kD"] <- log10(abs(fly_dyn_df$kD))
fly_dyn_df <- fly_dyn_df %>% arrange(desc(abs(kD)))

head(fly_dyn_df)
```

    ##   Protein_ID Gene_Symbol         kD    HL_Hrs mClass         Dyn    iBAQ_Conc
    ## 1     Q9VFD5      CG6966 0.01539471 0.7504170    Deg 0.074725848 1.406290e-06
    ## 2     Q9VMA3         cup 0.01486216 0.7773063    Deg 0.320167710 2.895060e-04
    ## 3     Q8IPM1         srl 0.01408818 0.8200105    Deg 0.115780501 7.180895e-07
    ## 4     Q9VED4        bard 0.01264941 0.9132802    Deg 0.325019851 9.597112e-06
    ## 5     P09085         cad 0.01217044 0.9492226    Deg 0.157655254 2.395412e-05
    ## 6     Q9VNG1      CG2182 0.01191980 0.9691817    Deg 0.006752819 3.286632e-06
    ##         kDa  Plot_iBAQ   Plot_Dyn     Plot_HL   Plot_kD
    ## 1  73.78584 -2.8519250 -1.1265292 -0.12469732 -1.812628
    ## 2 125.67200 -0.5383425 -0.4946225 -0.10940781 -1.827918
    ## 3 118.78416 -3.1438214 -0.9363646 -0.08618057 -1.851145
    ## 4  52.86154 -2.0178594 -0.4880901 -0.03939597 -1.897930
    ## 5  45.72605 -1.6206197 -0.8022916 -0.02263195 -1.914694
    ## 6  77.06133 -2.4832490 -2.1705149 -0.01359479 -1.923731

``` r
dyn_ibaq_plots <- function(df, gs_col) {
  
  ribo_df <- df[(grepl("^rpl\\d", df[[gs_col]],
                              ignore.case = TRUE)) &
                         (!grepl("like", df[[gs_col]],
                                 ignore.case = TRUE)),]
  int_df <- df[df$Protein_ID == "XBmRNA21384|XBXL10_1g11382",]
  
  # print(ribo_df)
  # print(median(ribo_df$HL_Hrs))
  
  p1 <- ggplot() +
    geom_point(data=df,
               aes(y=Plot_iBAQ, x=Plot_Dyn),
               shape=1, size=3, alpha=0.2, color="black") +
    
    geom_point(data=ribo_df,
               aes(y=Plot_iBAQ, x=Plot_Dyn),
               size=3, color="#E69F00") +
    geom_point(data=ribo_df,
               aes(y=Plot_iBAQ, x=Plot_Dyn),
               shape=1, size=3, color="black") +

    stat_density_2d(data = df, aes(y=Plot_iBAQ, x=Plot_Dyn),
                    color = "#CC79A7", linewidth = 2, bins = 3) +

    geom_point(data=int_df,
               aes(y=Plot_iBAQ, x=Plot_Dyn),
               size=3, color="#56B4E9") +
    geom_point(data=int_df,
               aes(y=Plot_iBAQ, x=Plot_Dyn),
               shape=1, size=3, color="black") +

    labs(x=expression("log"[10]*"( Dynamicity )"),
         y=expression("log"[10]*"( [μM] )")) +
    theme_bw() +
    theme(axis.text=element_text(size=22,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          axis.title = element_text(size=20),
          legend.position = "none",
          panel.border = element_rect(linewidth=1.5),
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          aspect.ratio = 1,
          plot.margin = margin(t = 10, r = 0, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(-6,2), xlim = c(-5,0))

  #-----------------------------
  
  deg_df <- df[df$mClass == "Deg",] %>% arrange(desc(kD))
  deg_df["Deg_Rank"] <- seq(1, nrow(deg_df))
  
  non_deg_ribo <- ribo_df[!ribo_df$Protein_ID %in% deg_df$Protein_ID,]
  # non_deg_ribo["Deg_Rank"] <- seq(nrow(deg_df) + 240, nrow(deg_df) + nrow(non_deg_ribo)+239)
  non_deg_ribo["Deg_Rank"] <- seq(nrow(deg_df) + 100, nrow(deg_df) + nrow(non_deg_ribo)+99)
    
  med_ibaq_line <- deg_df %>% arrange(Deg_Rank) %>%
    # Create a bin for every 10 ranks (1-10, 11-20, etc.)
    # Subtracting 1 ensures rank 10 stays in the first group
    group_by(Rank_Group = (Deg_Rank - 1) %/% 200) %>%
    summarize(
      Avg_Plot_iBAQ = median(Plot_iBAQ, na.rm = TRUE),
      Median_Rank = median(Deg_Rank, na.rm = TRUE)
    ) %>%
    ungroup() %>%
    select(Median_Rank, Avg_Plot_iBAQ)


  p2 <- ggplot() +
    annotate("rect", xmin = nrow(deg_df), xmax = 100000,
             ymin = -1000, ymax = 1000,
             fill = "grey80", alpha = 0.8) +
    
    geom_point(data=deg_df,
               aes(y=Plot_iBAQ, x=Deg_Rank),
               shape=1, size=3, alpha=0.3) +
    geom_line(data=med_ibaq_line,
              aes(x=Median_Rank, y=Avg_Plot_iBAQ, group=1),
              linewidth=2, color="red3") +
    
    geom_vline(xintercept=nrow(deg_df), linetype="dashed",
               linewidth=1) +

    geom_point(data=deg_df[deg_df$Protein_ID %in% ribo_df$Protein_ID,],
               aes(y=Plot_iBAQ, x=Deg_Rank),
               size=3, color="#E69F00") +
    geom_point(data=deg_df[deg_df$Protein_ID %in% ribo_df$Protein_ID,],
               aes(y=Plot_iBAQ, x=Deg_Rank),
               shape=1, size=3, color="black") +
    geom_point(data=non_deg_ribo,
               aes(y=Plot_iBAQ, x=Deg_Rank),
               size=3, color="#E69F00") +
    geom_point(data=non_deg_ribo,
               aes(y=Plot_iBAQ, x=Deg_Rank),
               shape=1, size=3, color="black") +

    geom_point(data=deg_df[deg_df$Protein_ID == "XBmRNA21384|XBXL10_1g11382",],
               aes(y=Plot_iBAQ, x=Deg_Rank),
               size=3, color="#56B4E9") +
    geom_point(data=deg_df[deg_df$Protein_ID == "XBmRNA21384|XBXL10_1g11382",],
               aes(y=Plot_iBAQ, x=Deg_Rank),
               shape=1, size=3, color="black") +    
    
    scale_color_manual(values=c("#5E4FA2", "grey80")) +    
        
    labs(x="Deg. speed (Rank 1 = fastest)",
         y=expression("log"[10]*"( [μM] )")) +
    theme_bw() +
    theme(axis.text=element_text(size=22,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          axis.title = element_text(size=20),
          panel.border = element_rect(linewidth=1.25),
          legend.position = "none",
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
          aspect.ratio = 1,
          plot.margin = margin(t = 10, r = 25, b = 10,
                               l = 10, unit = "pt")) +
    coord_cartesian(ylim=c(-6,2),
                    # xlim=c(0, nrow(deg_df) + nrow(non_deg_ribo)+240)) +
                    xlim=c(0, nrow(deg_df) + nrow(non_deg_ribo)+100)) +    
    # scale_x_continuous(breaks = c(1, 1000, 2000, 4000, 6000))
    scale_x_continuous(breaks = c(1, 1000, 3000, 5000))

  return(p1 | p2)
  
  }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=8, height=5, res=300)

dyn_ibaq_plots(frog_dyn_df, "Human_Gene")
```

![](figures/turnover_and_absolute/Relationship%20between%20turnover%20and%20turnover-1.png)<!-- -->

``` r
dyn_ibaq_plots(fly_dyn_df, "Gene_Symbol")
```

![](figures/turnover_and_absolute/Relationship%20between%20turnover%20and%20turnover-2.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
sessionInfo()
```

    ## R version 4.5.3 (2026-03-11 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26200)
    ## 
    ## Matrix products: default
    ##   LAPACK version 3.12.1
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United States.utf8 
    ## [2] LC_CTYPE=English_United States.utf8   
    ## [3] LC_MONETARY=English_United States.utf8
    ## [4] LC_NUMERIC=C                          
    ## [5] LC_TIME=English_United States.utf8    
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] patchwork_1.3.2     lsa_0.73.4          SnowballC_0.7.1    
    ##  [4] Peptides_2.4.6      Biostrings_2.78.0   Seqinfo_1.0.0      
    ##  [7] XVector_0.50.0      IRanges_2.44.0      S4Vectors_0.48.1   
    ## [10] BiocGenerics_0.56.0 generics_0.1.4      ggplot2_4.0.3      
    ## [13] stringr_1.6.0       fgsea_1.36.2        dplyr_1.2.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] stringi_1.8.7       lattice_0.22-9      digest_0.6.39      
    ##  [4] magrittr_2.0.4      evaluate_1.0.5      grid_4.5.3         
    ##  [7] RColorBrewer_1.1-3  fastmap_1.2.0       Matrix_1.7-5       
    ## [10] scales_1.4.0        isoband_0.3.0       codetools_0.2-20   
    ## [13] cli_3.6.5           rlang_1.1.7         crayon_1.5.3       
    ## [16] cowplot_1.2.0       withr_3.0.3         yaml_2.3.12        
    ## [19] otel_0.2.0          tools_4.5.3         parallel_4.5.3     
    ## [22] BiocParallel_1.44.0 fastmatch_1.1-8     vctrs_0.7.1        
    ## [25] R6_2.6.1            lifecycle_1.0.5     MASS_7.3-65        
    ## [28] pkgconfig_2.0.3     pillar_1.11.1       gtable_0.3.6       
    ## [31] data.table_1.18.4   glue_1.8.0          Rcpp_1.1.1-1.1     
    ## [34] xfun_0.58           tibble_3.3.1        tidyselect_1.2.1   
    ## [37] rstudioapi_0.19.0   knitr_1.51          farver_2.1.2       
    ## [40] htmltools_0.5.9     labeling_0.4.3      rmarkdown_2.31     
    ## [43] compiler_4.5.3      S7_0.2.1
