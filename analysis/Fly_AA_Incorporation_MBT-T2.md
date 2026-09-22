Fly AA Incorporation Gast T2
================
Edward Cruz
2026-04-24

``` r
rm(list=ls(all=T))
library(tidyr)
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
library(splines)
library(stringr)
library(ggplot2)
library(patchwork)
library(data.table)
```

    ## 
    ## Attaching package: 'data.table'

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

``` r
library(minpack.lm)
library(scam)
```

    ## This is scam 1.2-22.

``` r
library(parallel)

knitr::opts_chunk$set(fig.path = "figures/Fly_AA_Incorporation_MBT-T2/")
num_cores <- detectCores() - 1  # Save one core for system stability
```

``` r
aa.profiles <- read.csv("Files/Reference/AA_Profiles.csv")
ivanov_T2_AA <- read.csv("Data/Dmel_AAA/2026_03_16_Fly_AA_2.csv")

head(aa.profiles)
```

    ##   AA       Name   Formula     MW Nitrogen Carbon Oxygen Hydrogen Sulfur
    ## 1  A    alanine   C3H7NO2  89.10        1      3      2        7      0
    ## 2  R   arginine C6H14N4O2 174.20        4      6      2       14      0
    ## 3  N asparagine  C4H8N2O3 132.12        2      4      3        8      0
    ## 4  D  aspartate   C4H7NO4 133.11        1      4      4        7      0
    ## 5  C   cysteine  C3H7NO2S 121.16        1      3      2        7      1
    ## 6  E  glutamate   C5H9NO4 147.13        1      5      4        9      0

``` r
head(ivanov_T2_AA)
```

    ##         AA     Isotope T0A_1 T0A_2 T0B_1 T0B_2 T1_1 T1_2 T2_1 T2_2 T3_1
    ## 1  alanine  C12 PARENT     1     1     1     1    1    1    1    1    1
    ## 2  alanine O18-label-1     0     0     0     0    0    0    0    0    0
    ## 3  alanine O18-label-2     0     0     0     0    0    0    0    0    0
    ## 4 arginine  C12 PARENT     1     1     1     1    1    1    1    1    1
    ## 5 arginine O18-label-1     0     0     0     0    0    0    0    0    0
    ## 6 arginine O18-label-2     0     0     0     0    0    0    0    0    0
    ##        T3_2      T4_1 T4_2    T5_1      T5_2      T6_1     T6_2     T7_1
    ## 1 0.9261990 0.9400660    1 0.80414 0.9177050 0.8440580 0.828058 0.416063
    ## 2 0.0738005 0.0599343    0 0.19586 0.0822954 0.1259860 0.171942 0.442817
    ## 3 0.0000000 0.0000000    0 0.00000 0.0000000 0.0299557 0.000000 0.141120
    ## 4 1.0000000 1.0000000    1 1.00000 1.0000000 1.0000000 1.000000 0.509250
    ## 5 0.0000000 0.0000000    0 0.00000 0.0000000 0.0000000 0.000000 0.490750
    ## 6 0.0000000 0.0000000    0 0.00000 0.0000000 0.0000000 0.000000 0.000000
    ##       T7_2     T8_1     T8_2    T9_1     T9_2     T10_1    T10_2    T11_1
    ## 1 0.390624 0.527875 0.584823 0.67442 0.528575 0.3940670 0.420092 0.363563
    ## 2 0.413118 0.472125 0.290907 0.32558 0.288633 0.4052390 0.365469 0.387962
    ## 3 0.196258 0.000000 0.124270 0.00000 0.182792 0.2006940 0.214440 0.248474
    ## 4 0.543816 0.598866 0.655618 1.00000 0.536893 0.3481840 0.421541 0.392224
    ## 5 0.456184 0.401134 0.344382 0.00000 0.463107 0.5703130 0.578459 0.607776
    ## 6 0.000000 0.000000 0.000000 0.00000 0.000000 0.0815035 0.000000 0.000000
    ##      T11_2    T12_1    T12_2
    ## 1 0.273674 0.315786 0.313694
    ## 2 0.443658 0.549001 0.520217
    ## 3 0.282667 0.135213 0.166089
    ## 4 0.359279 0.377675 0.369665
    ## 5 0.640721 0.498746 0.529556
    ## 6 0.000000 0.123579 0.100779

``` r
MBT_AA_matrix <- lapply(unique(ivanov_T2_AA$AA),
                       function(AA){
  one_letter_code <- aa.profiles[aa.profiles$Name==AA, "AA"]
  sub_AA <- ivanov_T2_AA[ivanov_T2_AA$AA == AA,]

  tp_rows <- lapply(paste0("T", seq(0,12)), function(tp){

    out_df <- data.frame(matrix(NA, nrow = 1, ncol = 5 ))
    colnames(out_df) <- c("C12 PARENT", "O18-label-1", "O18-label-2",
                          "O18-label-3", "O18-label-4")

    col_regex <- paste0(tp, "[_|AB]")
    col_bind <- sub_AA[grepl(col_regex, colnames(sub_AA) )]

    col_meds <- apply(col_bind, 1, median)
    col_meds <- col_meds / sum(col_meds)

    out_df[1:length(col_meds)] <- col_meds
    out_df <- cbind(data.frame(Name=AA, AA=one_letter_code,
                               Total=length(col_meds)-1,
                               Sample=tp),
                    out_df)

    return(out_df) }) %>% bind_rows
  
  return(tp_rows) }) %>% bind_rows

MBT_AA_matrix <- MBT_AA_matrix[!MBT_AA_matrix$Sample %in% c("T4", "T10"),]
MBT_AA_matrix <- MBT_AA_matrix[!MBT_AA_matrix$Name == "histidine",]
```

``` r
measured_AAs <- unique(MBT_AA_matrix$AA)

no_Oxy_AA_all <- c("A", "R", "C", "G", "V",
                   "H", "I", "L", "K",
                   "M", "F", "P", "W",
                   "T", "Y") 
no_Oxy_AA <- c(no_Oxy_AA_all[no_Oxy_AA_all %in% measured_AAs], "X")

one_Oxy_AA_all <- c("Q", "S", "N") 
one_Oxy_AA <- c(one_Oxy_AA_all[one_Oxy_AA_all %in% measured_AAs], "U")

two_Oxy_AA_all <- c("D", "E")
two_Oxy_AA <- c(two_Oxy_AA_all[two_Oxy_AA_all %in% measured_AAs], "Z")
```

``` r
generate_unknown_meds <- function(AA_group, no_labels, AA_name) {

  sub_matrix <- MBT_AA_matrix[MBT_AA_matrix$AA %in% AA_group,]

  tp_meds <- lapply(unique(sub_matrix$Sample), function(tp){
    
    tp_sub <- sub_matrix[sub_matrix$Sample == tp,]
    tp_med <- apply(tp_sub[5:9], 2, median)

    return(data.frame(t(tp_med))) }) %>%
    bind_rows
  colnames(tp_meds) <- c("C12 PARENT", "O18-label-1", "O18-label-2",
                         "O18-label-3", "O18-label-4")
    
  tp_meds <- cbind(data.frame(Name=rep(AA_name, nrow(tp_meds)),
                              AA=rep(AA_name, nrow(tp_meds)),
                              Total=rep(no_labels, nrow(tp_meds)),
                              Sample=unique(sub_matrix$Sample)),
                   tp_meds)
  
  return(tp_meds) }  

matrix_med_X <- generate_unknown_meds(no_Oxy_AA_all, 2, "X")
matrix_med_U <- generate_unknown_meds(one_Oxy_AA_all, 3, "U")
matrix_med_Z <- generate_unknown_meds(two_Oxy_AA_all, 4, "Z")

MBT_AA_matrix <- rbind(MBT_AA_matrix, matrix_med_X,
                       matrix_med_U, matrix_med_Z)
MBT_AA_matrix$Name <- sub("^(\\w)", "\\U\\1", MBT_AA_matrix$Name, perl = TRUE)
head(MBT_AA_matrix)
```

    ##      Name AA Total Sample C12 PARENT O18-label-1 O18-label-2 O18-label-3
    ## 1 Alanine  A     2     T0  1.0000000  0.00000000  0.00000000          NA
    ## 2 Alanine  A     2     T1  1.0000000  0.00000000  0.00000000          NA
    ## 3 Alanine  A     2     T2  1.0000000  0.00000000  0.00000000          NA
    ## 4 Alanine  A     2     T3  0.9630997  0.03690026  0.00000000          NA
    ## 6 Alanine  A     2     T5  0.8609223  0.13907767  0.00000000          NA
    ## 7 Alanine  A     2     T6  0.8360581  0.14896402  0.01497785          NA
    ##   O18-label-4
    ## 1          NA
    ## 2          NA
    ## 3          NA
    ## 4          NA
    ## 6          NA
    ## 7          NA

``` r
#Time in minutes of timepoints collected, spline time is coarser
AA_T2.time <- c(0, 20, 40, 57, 95, 118, 180, 297, 421, 717, 1204)
```

``` r
build_C12_curve <- function(AA, aa_matrix, times,
                            interp_times = seq(0, 1204, by = 1),
                            onset_tol = 0.999) {
  sub.data <- aa_matrix[aa_matrix$AA == AA, ]
  sub.data["time_num"] <- as.numeric(gsub("T", "", sub.data$Sample))
  sub.data <- sub.data %>% arrange(time_num)

  tk <- times; yk <- sub.data$`C12 PARENT`

  # --- Isotonic fit to smooth curve ---
  iso_fit <- isoreg(tk, -yk)
  y_iso <- -iso_fit$yf
  
  #-------------------------------------------------

  n_pts <- length(unique(tk))
  k_use <- min(n_pts - 1, 8)
  
  fit <- scam(y ~ s(t, bs = "mpd", k = k_use),
              data = data.frame(t = tk, y = y_iso))
  pred <- pmin(pmax(predict(fit,
                            newdata = data.frame(t = interp_times)), 0),
               1)

  # plot(tk, yk, ylim=c(0,1), xlim=c(0,1210))
  # lines(tk, y_iso, col="blue")
  # lines(interp_times, pred, col="red")

  c_df <- data.frame(Type=rep("C12", length(pred)),
                     Time=interp_times, Value=pred)
  
  return(c_df) }


build_labeled_curve <- function(tk, y_label, name,
                                interp_times = seq(0, 1204, by = 1),
                                k_max = 6) {
  
  n_pts <- length(unique(tk))
  k_use <- min(n_pts - 1, k_max)

  w <- rep(1, length(tk))
  w[tk == 0] <- 1000    
  
  # Unconstrained smooth: allows rise AND fall (no monotone constraint)
  fit  <- scam(y ~ s(t, k = k_use),
               data = data.frame(t = tk, y = y_label),
               weights = w)
  pred <- as.numeric(predict(fit, newdata = data.frame(t = interp_times)))

  pred <- pmax(pred, 0) # floor at 0, no upper/monotone constraint

  c_df <- data.frame(Type=rep(name, length(pred)),
                     Time = interp_times, Value = pred)

  # plot(tk, y_label, ylim=c(0,1), xlim=c(0,1210))
  # lines(interp_times, pred, col="green")
  # print(y_label)
    
  return(c_df) }

#------------------------------------------

all_splines.merge <- apply(unique(MBT_AA_matrix[1:3]), 1,
                           function(row) {

  AA <- row["AA"]
  no_labels <- as.numeric(row["Total"])

  #Fit C12 data as isotonic regression then fit a spline through points
  C12_data <- build_C12_curve(AA, MBT_AA_matrix, AA_T2.time)

  #Fit spline bounded above 0 for labels
  label_data <- lapply(seq(no_labels), function(set){
    label_name <- paste0("O18-label-", set)
    label_data <- MBT_AA_matrix[MBT_AA_matrix$AA==AA, label_name]

    ind_label_data <- build_labeled_curve(AA_T2.time,
                                          label_data,
                                          label_name)

    return(ind_label_data) }) %>% bind_rows

  #Merge data together
  all_labels <- rbind(C12_data, label_data)

  #Set total sum at each t to 1
  cum_df <- lapply(seq(0, 1204, by = 1), function(t) {
    sub_t <- all_labels[all_labels$Time==t,]
    
    c12_row <- sub_t[sub_t$Type=="C12",]
    sub_t <- sub_t[!sub_t$Type=="C12",]
    
    leftover <- 1 - c12_row$Value
    sub_t["Value"] <- (sub_t$Value / sum(sub_t$Value)) * leftover
    sub_t <- rbind(c12_row, sub_t)
    sub_t[is.na(sub_t$Value), "Value"] <- 0
    
    return(sub_t) })
  cum_df <- do.call(rbind, cum_df)

  # sub_data <- MBT_AA_matrix[MBT_AA_matrix$AA==AA,2:(5+no_labels)]
  # sub_data["Total"] <- NULL
  # sub_data["Sample"] <- AA_T2.time
  # 
  # p1 <- ggplot() +
  #   geom_line(data=cum_df, aes(x=Time, y=Value,
  #                          color=Type),
  #             size=2) +
  #   geom_point(data=sub_data, aes(x=Sample, y=`C12 PARENT`),
  #              color="black", size=3) +
  #   geom_point(data=sub_data, aes(x=Sample, y=`O18-label-1`),
  #              color="green", size=3) +
  #   geom_point(data=sub_data, aes(x=Sample, y=`O18-label-2`),
  #              color="blue", size=3) +
  #   theme_bw() +
  #   scale_color_manual(values=c("red", "green", "blue", "purple", "pink"))
  # 
  # print(p1)

  cum_df["AA"] <- rep(AA, nrow(cum_df))
  cum_df["Labels"] <- rep(no_labels, nrow(cum_df))
  cum_df <- cum_df[c("AA", "Labels", "Time", "Type", "Value")]  
  
  return(cum_df) }) %>% bind_rows()

head(all_splines.merge)
```

    ##           AA Labels Time        Type Value
    ## 1...1      A      2    0         C12     1
    ## 11000...2  A      2    0 O18-label-1     0
    ## 1206...3   A      2    0 O18-label-2     0
    ## 2...4      A      2    1         C12     1
    ## 2411...5   A      2    1 O18-label-1     0
    ## 1207...6   A      2    1 O18-label-2     0

``` r
all_spline.p <- lapply(no_Oxy_AA, function(ind_AA){

  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA, ]
  sub.data["time_num"] <- as.numeric(gsub("T", "", sub.data$Sample))
  sub.data <- sub.data %>% arrange(time_num)
  sub.data["Time"] <- AA_T2.time/60 + (190/60)

  spline.data <- all_splines.merge[all_splines.merge$AA==ind_AA,]
  spline.data["Time"] <- spline.data$Time/60 + (190/60)

  AA_spline.p <- ggplot() +
    geom_line(data=spline.data[spline.data$Type=="O18-label-1",],
              aes(x=Time, y=Value),
              color="#009E73", size=1) +
    geom_line(data=spline.data[spline.data$Type=="O18-label-2",],
              aes(x=Time, y=Value),
              color="#B45100", size=1) +    
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=Value),
              color="black", size=3, alpha=0.4) +
    
    
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(190/60,(1204+190)/60), ylim=c(0,1))
  
  return(AA_spline.p) })
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once per session.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
all_spline.p <- c(all_spline.p,
                  lapply(one_Oxy_AA, function(ind_AA){

  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA, ]
  sub.data["time_num"] <- as.numeric(gsub("T", "", sub.data$Sample))
  sub.data <- sub.data %>% arrange(time_num)
  sub.data["Time"] <- AA_T2.time/60 + (190/60)

  spline.data <- all_splines.merge[all_splines.merge$AA==ind_AA,]
  spline.data["Time"] <- spline.data$Time/60 + (190/60)

  AA_spline.p <- ggplot() +
    geom_line(data=spline.data[spline.data$Type=="O18-label-1",],
              aes(x=Time, y=Value),
              color="#009E73", size=1) +
    geom_line(data=spline.data[spline.data$Type=="O18-label-2",],
              aes(x=Time, y=Value),
              color="#B45100", size=1) +    
    geom_line(data=spline.data[spline.data$Type=="O18-label-3",],
              aes(x=Time, y=`Value`),
              color="#0072B2", size=1) +
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=Value),
              color="black", size=3, alpha=0.4) +
    
    
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(190/60,(1204+190)/60), ylim=c(0,1))
  
  return(AA_spline.p) }) ) 
  
all_spline.p <- c(all_spline.p,
                  lapply(two_Oxy_AA, function(ind_AA){

  sub.data <- MBT_AA_matrix[MBT_AA_matrix$AA == ind_AA, ]
  sub.data["time_num"] <- as.numeric(gsub("T", "", sub.data$Sample))
  sub.data <- sub.data %>% arrange(time_num)
  sub.data["Time"] <- AA_T2.time/60 + (190/60)

  spline.data <- all_splines.merge[all_splines.merge$AA==ind_AA,]
  spline.data["Time"] <- spline.data$Time/60 + (190/60)

  AA_spline.p <- ggplot() +
    geom_line(data=spline.data[spline.data$Type=="O18-label-1",],
              aes(x=Time, y=Value),
              color="#009E73", size=1) +
    geom_line(data=spline.data[spline.data$Type=="O18-label-2",],
              aes(x=Time, y=Value),
              color="#B45100", size=1) +    
    geom_line(data=spline.data[spline.data$Type=="O18-label-3",],
              aes(x=Time, y=`Value`),
              color="#0072B2", size=1) +
    geom_line(data=spline.data[spline.data$Type=="O18-label-4",],
              aes(x=Time, y=`Value`),
              color="#CC79A7", size=1) +
    geom_line(data=spline.data[spline.data$Type=="C12",],
              aes(x=Time, y=Value),
              color="black", size=3, alpha=0.4) +
    
    
    geom_point(data=sub.data, aes(x=Time, y=`C12 PARENT`),
               color="black", size=3) +
    labs(x="Hours after Labeling", y="Relative Abundance",
         title=unique(sub.data$Name)) +
    theme_bw() +
    theme(axis.text=element_text(size=21,colour="black"),
          plot.title = element_text(size=24, hjust=0.5),
          # axis.title = element_text(size=22),
          axis.title = element_blank(),                    
          panel.border = element_rect(linewidth=2),
          panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
          aspect.ratio = 1,
          legend.position = "none") +
    coord_cartesian(xlim=c(190/60,(1204+190)/60), ylim=c(0,1))
  
  return(AA_spline.p) }) )          

names(all_spline.p) <- c(no_Oxy_AA, one_Oxy_AA, two_Oxy_AA)
```

``` r
# #Uncomment for new image!
# tiff("gRaPh_Me.tiff", units="in",
#      width=13, height=12, res=300)

(all_spline.p[["A"]] | all_spline.p[["R"]] | all_spline.p[["G"]] | all_spline.p[["V"]]) /
  (all_spline.p[["I"]] | all_spline.p[["L"]] | all_spline.p[["K"]] | all_spline.p[["M"]]) /
  (all_spline.p[["F"]] | all_spline.p[["P"]] | all_spline.p[["W"]] | all_spline.p[["T"]]) /
  (all_spline.p[["Y"]] | all_spline.p[["Y"]] | all_spline.p[["Y"]] | all_spline.p[["Y"]])
```

![](figures/Fly_AA_Incorporation_MBT-T2/Supp%20Figure%20of%20all%20amino%20acids%20-%20Part%201-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=15, height=3, res=300)

all_spline.p[["Q"]] | all_spline.p[["S"]] | all_spline.p[["N"]] | all_spline.p[["D"]] | all_spline.p[["E"]] 
```

![](figures/Fly_AA_Incorporation_MBT-T2/Supp%20Figure%20of%20all%20amino%20acids%20-%20Part%202-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=13/4*3, height=3, res=300)

all_spline.p[["X"]] | all_spline.p[["U"]] | all_spline.p[["Z"]]
```

![](figures/Fly_AA_Incorporation_MBT-T2/Median%20plots%20of%20AAs-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
AA_matrix_fit <- data.frame(Time=seq(0, 1204, by = 1))

AA_matrix_fit <- cbind(AA_matrix_fit,
                       do.call(cbind, lapply(unique(all_splines.merge$AA),
                                             function(AA){
                                               
  sub_data <- all_splines.merge[all_splines.merge$AA==AA&
                                  all_splines.merge$Type=="C12",]
  AA_labels <- unique(sub_data$Labels)

  if (AA %in% no_Oxy_AA) { p_ao <- sub_data$Value^(1/2)  }
  else if (AA %in% one_Oxy_AA) { p_ao <- sub_data$Value^(2/3)  }
  else if (AA %in% two_Oxy_AA) { p_ao <- sub_data$Value^(3/4)  }  
  
  out_col <- data.frame(AA_Name=p_ao)
  colnames(out_col) <- AA
  
  return(out_col) })))

matrix_header <- colnames(AA_matrix_fit)[2:length(colnames(AA_matrix_fit))]
AA_matrix_fit <- AA_matrix_fit %>% as.matrix

head(AA_matrix_fit)
```

    ##      Time A R N D E L I K F P S V Y Q W T M G X U Z
    ## [1,]    0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [2,]    1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [3,]    2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [4,]    3 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [5,]    4 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
    ## [6,]    5 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

``` r
filter_raw <- function(csv.file) {

  sn.cols <- c("126 Sn", "127n Sn", "127c Sn", "128n Sn", "128c Sn", "129n Sn",
               "129c Sn", "130n Sn", "130c Sn", "131n Sn", "131c Sn", "132n Sn",
               "132c Sn", "133n Sn", "133c Sn", "134n Sn", "134c Sn", "135n Sn")

  #Selecting columns needed for filtering and data analysis
  raw.columns <- c('Protein ID', 'Parsimony', 'Trimmed Peptide', sn.cols,
                   'Theo m/z', 'z', 'Isolation m/z')

  
  cat("Reading",  tail(unlist(strsplit(csv.file, "/")), n = 1), "...\n\n")  
    
  data.df <- fread(input = csv.file, select = raw.columns) %>%
    as.data.frame
  
  cat("Initial rows:\t\t\t\t\t", nrow(data.df), "\n")
  
  #'*Parsimony, REV sequences, contaminants, and oxM*
  data.df <- data.df[data.df$Parsimony %in% c("R", "U"),] #Drop NA Parsimony
  data.df <- data.df[!grepl("#", data.df$`Protein ID`),] #Remove reverse sequences
  data.df <- data.df[!grepl('contaminant', data.df$`Protein ID`),] #Remove contaminants
  data.df <- data.df[!grepl("\\*", data.df$`Trimmed Peptide`),] #Remove oxM
  
  cat("After Initial Parsimony/REV/Contam/oxM Filter:\t", nrow(data.df), "\n")

  data.df <- data.df[data.df$z==2 | data.df$z==3, ]

  cat("After Charge State Filter:\t\t\t", nrow(data.df), "\n")  
      
  #'*M0 Filter*  
  #Calculate theoretical/observed masses + "error"
  # ... Then select rows within error range
  data.df <- data.df[((data.df['z']*data.df['Isolation m/z']) >
                        ((data.df['z']*data.df['Theo m/z']) - 0.1)) &
                       ((data.df['z']*data.df['Isolation m/z']) <
                          ((data.df['z']*data.df['Theo m/z']) + 0.1)), ]
  
  cat("After M0 Filter:\t\t\t\t", nrow(data.df), "\n")  

  #'*Missed Cleavages*
  #Remove peptides with missed cleavages
  # ... Contains KR after removing last AA
  data.df <- data.df[!(grepl("K",
                        str_sub(data.df$`Trimmed Peptide`, end=-2))),]
  data.df <- data.df[!(grepl("R",
                        str_sub(data.df$`Trimmed Peptide`, end=-2))),]
  cat("After removing Missed Cleavages:\t\t", nrow(data.df), "\n")
  
  #'*TMTpro Signal*
  #Signal per channel
  data.df["Sn_per_channel"] <- apply(data.df, 1, function(x) {
    return(all(as.numeric(x[sn.cols]) >= 4)) })
  data.df <- data.df[data.df$Sn_per_channel==TRUE,]

  #Total signal
  data.df["sum_sn"] <- rowSums(data.df[sn.cols])
  data.df <- data.df[(data.df$z==2 & data.df$sum_sn>257)|
                       (data.df$z==3 & data.df$sum_sn>529),]  
    
  cat("After TMTpro Signal Filter:\t\t\t", nrow(data.df), "\n")  
  
  #Removing unnecessary columns and renaming data
  data.df[c("Parsimony", "Sn_per_channel", "z",
            "Theo m/z", "Isolation m/z")] <- NULL

  #Summing signal across peptides
  data.df <- data.df %>%
    group_by(`Protein ID`, `Trimmed Peptide`) %>%
    dplyr::summarise(across(where(is.numeric), \(x) sum(x, na.rm = TRUE)),
                     .groups = "drop") %>%
    as.data.frame()  

  #Max signal within peptide group
  # data.df <- data.df %>% group_by(`Trimmed Peptide`) %>%
  #   filter(sum_sn == max(sum_sn)) %>% # ... Select row with max signal
  #   ungroup

  cat("Final Unique Peptides:\t\t\t\t", nrow(data.df), "\n\n\n")

  return(data.df) }
  
# filter_raw("Data/Dmel_O18/Dmel-O18-T6_MS3-RTS_18plex_Fractionated_Merge.csv")
```

``` r
#'*MBT Timeseries Replicate 1*
Dmel_T6.df <- filter_raw("Data/Dmel_O18/Dmel-O18-T6_MS3-RTS_18plex_Fractionated_Merge.csv")
```

    ## Reading Dmel-O18-T6_MS3-RTS_18plex_Fractionated_Merge.csv ...
    ## 
    ## Initial rows:                     246006 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    191795 
    ## After Charge State Filter:            182536 
    ## After M0 Filter:              152241 
    ## After removing Missed Cleavages:      132949 
    ## After TMTpro Signal Filter:           96566 
    ## Final Unique Peptides:                52176

``` r
colnames(Dmel_T6.df) <- c("Protein_ID", "Peptide",
                          "T0A", "T0B", "O1", "O2", "O3", "O4", "O5",
                          "O6", "O7", "O8", "O9", "O10", "O11", "O12",
                          "N13", "N14", "N15", "N16", "sum_sn")
Dmel_T6.df["Protein_ID"] <- sapply(strsplit(Dmel_T6.df$Protein_ID, split = "\\|"),
                                   "[", 2)

#'*MBT Timeseries Replicate 2*
Dmel_T7.df <- filter_raw("Data/Dmel_O18/Dmel-O18-T7_MS3-RTS_18plex_Fractionated_Merge.csv")
```

    ## Reading Dmel-O18-T7_MS3-RTS_18plex_Fractionated_Merge.csv ...
    ## 
    ## Initial rows:                     237602 
    ## After Initial Parsimony/REV/Contam/oxM Filter:    177149 
    ## After Charge State Filter:            169742 
    ## After M0 Filter:              142722 
    ## After removing Missed Cleavages:      126437 
    ## After TMTpro Signal Filter:           92947 
    ## Final Unique Peptides:                50908

``` r
colnames(Dmel_T7.df) <- c("Protein_ID", "Peptide",
                          "T0A", "T0B", "O1", "O2", "O3", "O4", "O5",
                          "O6", "O7", "O8", "O9", "O10", "O11", "O12",
                          "N13", "N14", "N15", "N16", "sum_sn")
Dmel_T7.df["Protein_ID"] <- sapply(strsplit(Dmel_T7.df$Protein_ID, split = "\\|"),
                                   "[", 2)
```

``` r
#Exponential decay function with two parameters
decay_fcn <- function(t, k1, k2) { k1* exp(-k2*t) }

#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CH]", "X", peptide)

  pep_table <- table(strsplit(mod_pep,"")) #Splitting to structure matrix

  #Creating empty matrix of the entire amino acid table
  peptide_matrix <- matrix(0, nrow = nrow(AA_matrix_fit),
                           ncol = ncol(AA_matrix_fit) - 1)  
  colnames(peptide_matrix) <- colnames(AA_matrix_fit)[-1]

  #Replacing AA column with number of amino acids
  for (i in names(pep_table)) { peptide_matrix[,i] <- as.numeric(pep_table[i]) }

  calc_per_AA <- AA_matrix_fit[,-1] ^ peptide_matrix
  pep_deg.tot <- apply(calc_per_AA, 1, prod)
  pep_deg.tot <- pep_deg.tot / pep_deg.tot[1]  

  # Fitting using nlsLM which is less limited for starting values
  #   decay_fcn is an exponential decay function with two parameters defined above
  # Note: You need to constrain k1 or else it will find solutions above 1, which are meaningless.
  #       I've set the max k1 to be the unlabeled fraction at T0 which is reasonable.
  #       Can't use 1 -> too strong of a local minimum -> fixes k1 = 1
  decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 1204, by = 1), 1, k2),
                     start = list(k2 = 0.04),
                     lower = c(k2 = 0),
                     upper = c(k2 = Inf),
                     control = list(maxiter = 1000))
  fit.values <- as.numeric(coef(decay.fit)) #Extracting fit  

  # print(peptide)
  # plot(seq(0, 1204, by = 1)/60, pep_deg.tot,
  #      xlim=c(0,20), ylim=c(0,1), pch=16,
  #      ylab="Probability of Light Synthesis",
  #      xlab="Time (min)")
  # lines(seq(0,1204,1)/60, decay_fcn(seq(0,1204,1), 1,
  #                                   fit.values[1]),
  #       col="red")
  # print(decay.fit)
  # invokeRestart("abort")
  
  return(data.frame(t(c(1, fit.values)))) }

#All unique identified peptides
mbt_fly_peps <- unique(c(Dmel_T6.df$Peptide, Dmel_T7.df$Peptide))

#Applying function to fit each peptide
mbt_fly_pep_fits <- data.frame(Peptide = mbt_fly_peps)

#------------------------------------------------------------

#Setup cluster
cl <- makeCluster(num_cores)

#Export necessary variables/functions to the cluster nodes
clusterExport(cl,
              varlist=c("mbt_fly_peps", "peptide_AA_fit", "decay_fcn",
                        "AA_matrix_fit"))
invisible(clusterEvalQ(cl, library(minpack.lm)))

mbt_fly_pep_fits <- cbind(mbt_fly_pep_fits,
                          do.call(rbind, parLapply(cl, mbt_fly_peps, peptide_AA_fit)))
colnames(mbt_fly_pep_fits)[2:3] <- c("Pep_k1", "Pep_k2")

stopCluster(cl) #Stop Cluster

#------------------------------------------------------------

head(mbt_fly_pep_fits)
```

    ##       Peptide Pep_k1      Pep_k2
    ## 1 ANALEMVNEIR      1 0.015914141
    ## 2     ESYIAYR      1 0.008395398
    ## 3     LAAYYAK      1 0.005535717
    ## 4     LDIDPDK      1 0.009292738
    ## 5    LGGWPLIK      1 0.007442714
    ## 6     LPNLLYK      1 0.008213647

``` r
filtered_fly_T6.df <- merge(Dmel_T6.df,
                            mbt_fly_pep_fits, by="Peptide") %>%
  arrange(Protein_ID) %>% dplyr::select(Protein_ID, everything())

filtered_fly_T7.df <- merge(Dmel_T7.df,
                            mbt_fly_pep_fits, by="Peptide") %>%
  arrange(Protein_ID) %>% dplyr::select(Protein_ID, everything())

write.csv(mbt_fly_pep_fits,
          "Files/Fits/Fly_O18-MBT-TS_Theo-Peptide-Decay.csv",
          row.names = FALSE)
write.csv(filtered_fly_T6.df,
          "Data/Dmel_O18/Fly_O18-T6_Filtered_Data.csv",
          row.names = FALSE)
write.csv(filtered_fly_T7.df,
          "Data/Dmel_O18/Fly_O18-T7_Filtered_Data.csv",
          row.names = FALSE)
```

``` r
min_pep_k2 <- mbt_fly_pep_fits %>%
  filter(min(Pep_k2)==Pep_k2)

max_pep_k2 <- mbt_fly_pep_fits %>%
  filter(max(Pep_k2)==Pep_k2)

med_pep_k2 <- mbt_fly_pep_fits %>%
  mutate(diff = abs(Pep_k2 - median(Pep_k2))) %>%
  filter(diff == min(diff)) %>%
  dplyr::select(-diff)
med_pep_k2 <- med_pep_k2[1,]

example_fits <- rbind(min_pep_k2, med_pep_k2, max_pep_k2)
example_fits["Type"] <- c("Min", "Median", "Max")

ex.fit_data <- lapply(example_fits$Peptide, function(peptide) {
  
  ex.fit <- example_fits[example_fits$Peptide==peptide,]

  curve_values <- decay_fcn(seq(0,1204,1), ex.fit$Pep_k1, ex.fit$Pep_k2)  

  example.df <- data.frame(Time=seq(0,1204,1)/60,
                           Fit_values=curve_values,
                           Type=rep(ex.fit$Type, length(curve_values)))
  
  return(example.df) })

ex.fit_data <- do.call(rbind, ex.fit_data) 

ex.p <- ggplot() +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Min",], 
             aes(x=Time, y=Fit_values), color="#CC79A7", size=1.5) +
  geom_line(data=ex.fit_data[ex.fit_data$Type=="Median",], 
             aes(x=Time, y=Fit_values), color="#009E73", size=1.5) +
  geom_hline(yintercept=1, linetype="dashed", linewidth=1) +
  theme_bw() +
  labs(x="Hours after Labeling", y="Probability of light peptide") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,20))
```

``` r
round(log(2)/example_fits$Pep_k2,0)
```

    ## [1] 278  55  21

``` r
# #Uncomment for new image!
# tiff("Graphs/Dmel_AA_Incorporation/Fly_T2_Example_Fits.tiff", units="in",
#      width=4.5, height=4.5, res=300)

ex.p
```

![](figures/Fly_AA_Incorporation_MBT-T2/Example%20fits%20plots-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
example_fits
```

    ##                                   Peptide Pep_k1      Pep_k2   Type
    ## 1                               YGGGGGGGR      1 0.002493333    Min
    ## 2                              LTSSILDDCR      1 0.012705695 Median
    ## 3 LGLGIDEDEPMTTDDAQSAGDAPSLVEDTEDASHMEEVD      1 0.033231847    Max

``` r
#Function to fit all peptide decay
peptide_AA_fit <- function(peptide) {
  # Z not used
  # CHLA -> X
  # YN -> U
  mod_pep <- gsub("[CH]", "X", peptide)

  pep_table <- table(strsplit(mod_pep,"")) #Splitting to structure matrix

  #Creating empty matrix of the entire amino acid table
  peptide_matrix <- matrix(0, nrow = nrow(AA_matrix_fit),
                           ncol = ncol(AA_matrix_fit) - 1)  
  colnames(peptide_matrix) <- colnames(AA_matrix_fit)[-1]

  #Replacing AA column with number of amino acids
  for (i in names(pep_table)) { peptide_matrix[,i] <- as.numeric(pep_table[i]) }

  calc_per_AA <- AA_matrix_fit[,-1] ^ peptide_matrix
  pep_deg.tot <- apply(calc_per_AA, 1, prod)
  pep_deg.tot <- pep_deg.tot / pep_deg.tot[1]  

  # # Fitting using nlsLM which is less limited for starting values
  # #   decay_fcn is an exponential decay function with two parameters defined above
  # # Note: You need to constrain k1 or else it will find solutions above 1, which are meaningless.
  # #       I've set the max k1 to be the unlabeled fraction at T0 which is reasonable.
  # #       Can't use 1 -> too strong of a local minimum -> fixes k1 = 1
  # decay.fit <- nlsLM(pep_deg.tot ~ decay_fcn(seq(0, 1204, by = 1), 1, k2),
  #                    start = list(k2 = 0.04),
  #                    lower = c(k2 = 0),
  #                    upper = c(k2 = Inf),
  #                    control = list(maxiter = 1000))
  # fit.values <- as.numeric(coef(decay.fit)) #Extracting fit  

  # print(peptide)
  # plot(seq(0, 1204, by = 1)/60, pep_deg.tot,
  #      xlim=c(0,20), ylim=c(0,1), pch=16,
  #      ylab="Probability of Light Synthesis",
  #      xlab="Time (min)")
  # lines(seq(0,1204,1)/60, decay_fcn(seq(0,1204,1), 1,
  #                                   fit.values[1]),
  #       col="red")
  # print(decay.fit)
  # invokeRestart("abort")
  
  return(pep_deg.tot) }


median_pep_residues <- data.frame(Time = seq(0, 1204, by = 1)/60,
                                  Fit_values = peptide_AA_fit("LTSSILDDCR"))
longest_pep_residues <- data.frame(Time = seq(0, 1204, by = 1)/60,
                                   Fit_values = peptide_AA_fit("YGGGGGGGR"))

res.p <- ggplot() +
  geom_point(data=longest_pep_residues, 
             aes(x=Time, y=Fit_values), color="#CC79A7", size=2,
             shape=1) +
  geom_point(data=median_pep_residues, 
             aes(x=Time, y=Fit_values), color="#009E73", size=2,
             shape=1) +

  theme_bw() +
  labs(x="Hours after Labeling", y="Probability of light peptide") +
  theme_bw() +
  theme(axis.text=element_text(size=21,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=22),
        panel.border = element_rect(linewidth=2),
        panel.grid.major = element_line(color = "grey70", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey70", linewidth = 0.25),
        aspect.ratio = 1,
        legend.position = "none") +
  coord_cartesian(ylim=c(0,1), xlim=c(0,20))

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=4, height=4, res=300)

res.p
```

![](figures/Fly_AA_Incorporation_MBT-T2/Plots%20for%20supplemental%20figure-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
median_pep_residues[median_pep_residues$Fit_values < 0.5,][1,]
```

    ##        Time Fit_values
    ## 63 1.033333  0.4980982

``` r
longest_pep_residues[longest_pep_residues$Fit_values < 0.5,][1,]
```

    ##         Time Fit_values
    ## 224 3.716667  0.4989351

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
    ## [1] parallel  splines   stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ## [1] scam_1.2-22       minpack.lm_1.2-4  data.table_1.18.4 patchwork_1.3.2  
    ## [5] ggplot2_4.0.3     stringr_1.6.0     dplyr_1.2.0       tidyr_1.3.2      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Matrix_1.7-5       gtable_0.3.6       compiler_4.5.3     tidyselect_1.2.1  
    ##  [5] scales_1.4.0       yaml_2.3.12        fastmap_1.2.0      lattice_0.22-9    
    ##  [9] R6_2.6.1           labeling_0.4.3     generics_0.1.4     knitr_1.51        
    ## [13] tibble_3.3.1       pillar_1.11.1      RColorBrewer_1.1-3 rlang_1.1.7       
    ## [17] stringi_1.8.7      xfun_0.58          S7_0.2.1           otel_0.2.0        
    ## [21] cli_3.6.5          mgcv_1.9-4         withr_3.0.3        magrittr_2.0.4    
    ## [25] digest_0.6.39      grid_4.5.3         rstudioapi_0.19.0  nlme_3.1-169      
    ## [29] lifecycle_1.0.5    vctrs_0.7.1        evaluate_1.0.5     glue_1.8.0        
    ## [33] farver_2.1.2       rmarkdown_2.31     purrr_1.2.2        tools_4.5.3       
    ## [37] pkgconfig_2.0.3    htmltools_0.5.9
