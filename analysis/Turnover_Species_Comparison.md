Comparing Fly and Frog Protein Turnover
================
Edward Cruz
2025-05-08

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
library(tidyr)
library(stringr)
library(ggplot2)
library(MASS)
```

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
knitr::opts_chunk$set(fig.path = "figures/Turnover_Species_Comparison/")
```

``` r
frog_min_kd <- log(2)/(4244*3)
frog_max_hl <- log(2)/frog_min_kd/60
```

``` r
fly_fits <- read.csv("Files/Fits/Dyn-Model/Fly_O18_T6-7_FinalFits_ALL.csv")
frog_fits <-read.csv("Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv")

head(fly_fits)
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
head(frog_fits)
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
frog_fits_only <- frog_fits[frog_fits$Estimate=="interval",]
frog_fits_only <- frog_fits_only[c("Protein_ID", "XLA_Gene", "kD",
                                   "HL_Hrs", "mClass")]
colnames(frog_fits_only)[5] <- "Frog_Class"

fly_fits_only <- fly_fits[fly_fits$Estimate=="interval",]
fly_fits_only <- fly_fits_only[c("Protein_ID", "Gene_Symbol", "kD",
                                 "HL_Hrs")]
colnames(fly_fits_only)[2] <- "Fly_Gene"

#----------------------------------

sc_df <- read.csv("Data/XLA-Dmel_orthologue_pairs.csv")
colnames(sc_df) <- c("Fly_ID", "Frog_ID")
sc_df["Fly_ID"] <- sapply(strsplit(sc_df$Fly_ID, split = "\\|"), "[", 2)

sc_df <- merge(sc_df, fly_fits_only,
               by.x = "Fly_ID", by.y="Protein_ID")
sc_df <- merge(sc_df, frog_fits_only,
               by.x = "Frog_ID", by.y="Protein_ID",
               suffixes = c("_Fly", "_Frog"))

sc_df <- sc_df %>%
  dplyr::select(Frog_ID, Fly_ID, XLA_Gene, Fly_Gene, everything()) %>%
  arrange(desc(kD_Frog))

head(sc_df) 
```

    ##                      Frog_ID Fly_ID      XLA_Gene Fly_Gene      kD_Fly
    ## 1 XBmRNA36839|XBXL10_1g19916 P29993  LOC108714625     Itpr 0.001027758
    ## 2 XBmRNA69786|XBXL10_1g37002 Q9VQ91       tdrkh.L     papi 0.001306100
    ## 3 XBmRNA56526|XBXL10_1g30025 Q9VND2     dclre1a.L     Snm1 0.001604518
    ## 4  XBmRNA11752|XBXL10_1g6259 Q9VNA8      donson.L       hd 0.002088953
    ## 5   XBmRNA3736|XBXL10_1g2050 P56175        rcl1.L     Rtc1 0.001000166
    ## 6 XBmRNA45716|XBXL10_1g24428 Q9VRJ0  LOC108717192      Gen 0.001413296
    ##   HL_Hrs_Fly      kD_Frog HL_Hrs_Frog Frog_Class
    ## 1  11.240446 0.0034152334    3.382625        Deg
    ## 2   8.845000 0.0027468404    4.205724        Deg
    ## 3   7.199952 0.0011060341   10.444934        Deg
    ## 4   5.530260 0.0008580523   13.463577        Deg
    ## 5  11.550534 0.0007876533   14.666927        Deg
    ## 6   8.174121 0.0006915136   16.706039        Deg

``` r
paired_kd_comp <- data.frame(kD=c(sc_df$kD_Fly,
                                  sc_df$kD_Frog),
                             Group=c(rep("Fly", nrow(sc_df)),
                                     rep("Frog", nrow(sc_df))))
# paired_kd_comp <- data.frame(kD=c(fly_fits$kD,
#                                   frog_fits$kD),
#                              Group=c(rep("Fly", nrow(fly_fits)),
#                                      rep("Frog", nrow(frog_fits))))
paired_kd_comp["Plot_kD"] <- paired_kd_comp$kD
paired_kd_comp[paired_kd_comp$kD > 0.001, "Plot_kD"] <- 0.001

#Define exact bin edges to guarantee the bars and the outline align perfectly
bw <- 0.000025
brks <- seq(-0.001, 0.002 + bw, by = bw)

#Calculate half-life x-intercepts (kd = ln(2) / minutes)
kd_1day <- log(2) / (24 * 60)
kd_12hr <- log(2) / (12 * 60)
kd_6hr  <- log(2) / (6 * 60)

comp_kd_pdf <- ggplot(paired_kd_comp,
                      aes(x = Plot_kD, fill = Group,
                          color = Group)) +
  # after_stat(density) normalizes the area to 1
  # position = "identity" prevents stacking so they overlap naturally
  geom_histogram(aes(y = after_stat(count / sum(count))), 
                 breaks = brks, position = "identity", 
                 alpha=0.7, color=NA) +
  
  # ONLY draws an outline for "Deg"
  geom_step(aes(y = after_stat(count / sum(count))),
            stat = "bin", breaks = brks, direction = "mid",
            linewidth = 1, show.legend = FALSE) +
  
  geom_vline(xintercept = c(kd_1day, kd_12hr, kd_6hr),
             linetype = "dashed", color = "black", 
             linewidth = 1) +
  
  scale_fill_manual(values = c("Fly" = "#CC79A7", "Frog" = "#E69F00")) +
  scale_color_manual(values = c("Fly" = "#CC79A7", "Frog" = "#E69F00")) +

  # LABELS & THEME: Matching your reference image
  labs(x = expression("k"["d"]~"(min"^{-1}*")"), y = "Density") +
  theme_bw() +
  theme(axis.text.y=element_text(size=20,colour="black"),
        axis.text.x=element_text(size=18,colour="black", angle=45,
                                 hjust=1),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=18),
        panel.border = element_rect(linewidth=1.5),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25),
        plot.margin = margin(t = 10, r = 25, b = 10,
                             l = 10, unit = "pt")) +
  coord_cartesian(xlim=c(0, 0.001)) +
  scale_x_continuous(breaks=c(0, 0.00025, 0.0005, 0.00075, 0.001),
                     labels=c("0", "0.00025", "0.0005",
                              "0.00075", "0.001"))

median(paired_kd_comp[paired_kd_comp$Group=="Frog",]$kD)
```

    ## [1] 1.742766e-05

``` r
median(paired_kd_comp[paired_kd_comp$Group=="Fly",]$kD)
```

    ## [1] 0.0004989113

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=6, height=3.5, res=300)

comp_kd_pdf
```

![](figures/Turnover_Species_Comparison/Histograms%20comparing%20all%20frog%20and%20fly%20single-copy%20half-lives-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# Fit a robust linear model forcing a zero intercept (y = mx)
# In this formula, '0 +' tells R to remove the intercept term.
robust_ortho <- rlm(HL_Hrs_Frog ~ 0 + HL_Hrs_Fly, data = sc_df)

# Extract just the slope coefficient
robust_ortho <- coef(robust_ortho)[1]
print(paste("The robust slope (m) is:", round(robust_ortho, 1)))
```

    ## [1] "The robust slope (m) is: 7.4"

``` r
p1 <- ggplot() +
  geom_point(data=sc_df, aes(x=HL_Hrs_Fly, y=HL_Hrs_Frog),
             size=3, shape=1, stroke=1) +
  geom_abline(slope=1, intercept=0, color="black",
              linetype="dashed", linewidth=1) +
  geom_abline(slope=robust_ortho, intercept=0, color="red2",
              linetype="dashed", linewidth=1) +

  labs(x="Fly half-life (Hrs)", y="Frog half-life (Hrs)") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=20),
        panel.border = element_rect(linewidth=1.25),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none",
        plot.margin = margin(t = 5, r = 10, b = 5,
                             l = 5, unit = "pt")) +
  coord_cartesian(xlim=c(0,frog_max_hl), ylim=c(0,frog_max_hl))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4, height=4, res=300)

p1
```

![](figures/Turnover_Species_Comparison/Biplot%20of%20frog%20and%20fly%20half-lives-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
hl_ratio <- sc_df
hl_ratio["Ratio_HL"] <- hl_ratio$HL_Hrs_Frog / hl_ratio$HL_Hrs_Fly
hl_ratio[hl_ratio$Ratio_HL > 20, "Ratio_HL"] <- 20

comp_violin <- ggplot(hl_ratio,
                      aes(x = "", y = Ratio_HL)) +
  geom_violin(fill = "grey70") +
  labs(x=paste("Median:",
               round(median(hl_ratio$Ratio_HL), 1)),
       y="Frog / fly half-life FC") +

  stat_summary(fun = median, geom = "point", shape = 23,
               size = 3, fill = "black") +  # Add median
  stat_summary(fun.data = function(y) { data.frame(y = median(y),
                                                   ymin = quantile(y, 0.25),
                                                   ymax = quantile(y, 0.75))},
               geom = "errorbar", width = 0.1) +  # Add IQR
  theme_bw() +
  theme(#aspect.ratio = 1,
        plot.title = element_text(size=20,colour="black"),
        axis.text=element_text(size=20,colour="black"),
        axis.title.x=element_text(size=20,colour="black"),
        axis.title.y=element_text(size=20,colour="black"),
        panel.border = element_rect(linewidth=1.25))
  # coord_cartesian(ylim=c(0,100))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=2.55, height=4, res=300)

comp_violin
```

![](figures/Turnover_Species_Comparison/Violin%20plot%20of%20single-copy%20pair%20fold-change-1.png)<!-- -->

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
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] MASS_7.3-65   ggplot2_4.0.3 stringr_1.6.0 tidyr_1.3.2   dplyr_1.2.0  
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] vctrs_0.7.1        cli_3.6.5          knitr_1.51         rlang_1.1.7       
    ##  [5] xfun_0.58          stringi_1.8.7      otel_0.2.0         purrr_1.2.2       
    ##  [9] generics_0.1.4     S7_0.2.1           labeling_0.4.3     glue_1.8.0        
    ## [13] htmltools_0.5.9    scales_1.4.0       rmarkdown_2.31     grid_4.5.3        
    ## [17] evaluate_1.0.5     tibble_3.3.1       fastmap_1.2.0      yaml_2.3.12       
    ## [21] lifecycle_1.0.5    compiler_4.5.3     RColorBrewer_1.1-3 pkgconfig_2.0.3   
    ## [25] rstudioapi_0.19.0  farver_2.1.2       digest_0.6.39      R6_2.6.1          
    ## [29] tidyselect_1.2.1   pillar_1.11.1      magrittr_2.0.4     withr_3.0.3       
    ## [33] gtable_0.3.6       tools_4.5.3
