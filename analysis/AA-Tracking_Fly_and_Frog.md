Amino acid tracking
================
Edward Cruz
2026-04-23

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

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     expand

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
library(stringr)
library(ggbreak)
```

    ## ggbreak v0.1.7 Learn more at https://yulab-smu.top/

    ## If you use ggbreak in published research, please cite the following
    ## paper:
    ## 
    ## S Xu, M Chen, T Feng, L Zhan, L Zhou, G Yu. Use ggbreak to effectively
    ## utilize plotting space to deal with large datasets and outliers.
    ## Frontiers in Genetics. 2021, 12:774846. doi: 10.3389/fgene.2021.774846

``` r
library(patchwork)

knitr::opts_chunk$set(fig.path = "figures/AA-Tracking_Fly_and_Frog/")
```

``` r
Fly_heavy_time <- c(257, 287, 320, 351, 386, 441, 500,
                     565, 682, 868, 995, 1218, 1514) - 257
Fly_light_time <- c(0, Fly_heavy_time[c(4, 8, 11, 13)])

Fly_heavy_hpf_time <- ((Fly_heavy_time + 257) + (Fly_heavy_time + 117))/2
Fly_heavy_hpf_time <- Fly_heavy_hpf_time / 60

Fly_light_hpf_time <- c(Fly_heavy_hpf_time[1], Fly_heavy_hpf_time[c(4, 8, 11, 13)])
```

``` r
Frog_T5_all.min <- c(161, 221, 281, 401, 521, 641,
                     881, 1605, 1961, 2918, 4405) - 161
Frog_T5_light.min <- Frog_T5_all.min[c(1:2, 4:5, 7:8, 10:11)]
Frog_T6_all.min <- c(150, 210, 270, 390, 510, 635,
                     930, 1591, 1950, 2910, 4394) - 150
Frog_T6_light.min <- Frog_T6_all.min[c(1:2, 4:5, 7:8, 10:11)]

Frog_heavy_time <- (Frog_T5_all.min + Frog_T6_all.min)/2
Frog_light_time <- (Frog_T5_light.min + Frog_T6_light.min)/2

Frog_heavy_hpf_time <- (Frog_heavy_time + ((161+150)/2))/60
Frog_light_hpf_time <- (Frog_light_time + ((161+150)/2))/60
```

``` r
XLA_2C_fits <- read.csv("Files/Fits/Dyn-Model/XLA_O18_T5-6_FinalFits_ALL.csv")
Fly_gast_fits <- read.csv("Files/Fits/Dyn-Model/Fly_O18_T6-7_FinalFits_ALL.csv")
```

``` r
fly_yolk_set <- c("P02843", "P02844", "P06607")
names(fly_yolk_set) <- c("Yp1", "Yp2", "Yp3")

frog_yolk_set <- c("XBmRNA35014|XBXL10_1g18990", "XBmRNA39801|XBXL10_1g21498",
                   "XBmRNA35013|XBXL10_1g18989", "XBmRNA35012|XBXL10_1g18988",
                   "XBmRNA67639|XBXL10_1g35803")
names(frog_yolk_set) <- c("vtgb1.L", "vtgb1.S", "vtga1.L", "vtga2.L", "serpina")
```

``` r
aa_vector <- c("A", "C", "D", "E", "F", "G", "H", "I", 
               "K", "L", "M", "N", "P", "Q", "R", "S", 
               "T", "V", "W", "Y")
emb_essential_AAs <- c("F", "W", "I", "L", "V", "K",
                       "H", "T", "M", "C", "Y")
```

``` r
XLA_RT0_Abs <- read.csv("Files/AbsC/XLA_RT0_Abs-Concentrations.csv")
XLA_N12_Abs <- read.csv("Files/AbsC/XLA_N12_Abs-Concentrations.csv")
Dmel_T0_Abs <- read.csv("Files/AbsC/Dmel_T0_Abs-Concentrations.csv")
Dmel_N16_Abs <- read.csv("Files/AbsC/Dmel_N16_Abs-Concentrations.csv")

XLA_abs_conc <- merge(XLA_RT0_Abs[1:2], XLA_N12_Abs[c(1:2,4)], by="Ind_Protein_ID",
      suffixes=c("_Start", "_Late"))
Dmel_abs_conc <- merge(Dmel_T0_Abs[1:2], Dmel_N16_Abs[c(1:2,4)], by="Ind_Protein_ID",
      suffixes=c("_Start", "_Late"))

XLA_abs_conc <- XLA_abs_conc[XLA_abs_conc$Ind_Protein_ID %in% 
                               c(XLA_2C_fits$Protein_ID, frog_yolk_set),]
Dmel_abs_conc <- Dmel_abs_conc[Dmel_abs_conc$Ind_Protein_ID %in%
                                 c(Fly_gast_fits$Protein_ID, fly_yolk_set),]
```

``` r
GB_NYS_df <- read.csv("Data/XLA_Norm/XLA-O18_GB-FinalDataset.csv")
XLA_yolk_trend <- GB_NYS_df[GB_NYS_df$Protein_ID %in% frog_yolk_set,]

XLA_yolk_FC <- XLA_yolk_trend$T8 / XLA_yolk_trend$T0
names(XLA_yolk_FC) <- XLA_yolk_trend$Protein_ID

modify_frog_abs <- function(df){

  nonyolk_df <- df[!df$Ind_Protein_ID %in% frog_yolk_set,]

  nonyolk_fits <- XLA_2C_fits[XLA_2C_fits$Protein_ID %in% nonyolk_df$Ind_Protein_ID,]
  
  nonyolk_fits["lightFC"] <- nonyolk_fits$N12 / nonyolk_fits$T0_L
  nonyolk_fits["heavyFC"] <- nonyolk_fits$O12 / nonyolk_fits$T0_H
  nonyolk_fits[nonyolk_fits$heavyFC > 1, "heavyFC"] <- 1

  nonyolk_fits["AbsCategory"] <- rep("CHANGE", nrow(nonyolk_fits))
  nonyolk_fits[nonyolk_fits$lightFC > 1, "AbsCategory"] <- "Inc"
  nonyolk_fits[nonyolk_fits$lightFC < 1, "AbsCategory"] <- "Dec"

  nonyolk_final <- lapply(nonyolk_df$Ind_Protein_ID, function(x){
    fit_row <- nonyolk_fits[nonyolk_fits$Protein_ID == x,]
    abs_row <- nonyolk_df[nonyolk_df$Ind_Protein_ID == x,]
    start_conc <- abs_row$iBAQ_Conc_Start
    
    T_end <- tail(Frog_light_time, 1)

    if (fit_row$mClass == "Deg") { surv <- exp(-fit_row$kD * T_end) }
    else { surv <- 1 }
    
    conc_deg    <- start_conc * (1 - surv)
    conc_needed <- max(start_conc * fit_row$lightFC - start_conc * surv, 0)

    outline <- data.frame(Protein_ID = x, StartConc = start_conc,
                          Conc_Deg = conc_deg,
                          Conc_Synth = conc_needed,
                          Sequence = abs_row$Sequence)

    return(outline) }) %>% bind_rows

  #------------------------------------------
  
  yolk_df <- df[df$Ind_Protein_ID %in% frog_yolk_set,]

  yolk_final <- lapply(yolk_df$Ind_Protein_ID, function(yolk_id){
    row <- yolk_df[yolk_df$Ind_Protein_ID == yolk_id,]
    
    start_conc <- row$iBAQ_Conc_Start
    nospin_FC <- XLA_yolk_FC[yolk_id]
    
    conc_deg <- start_conc-(start_conc * nospin_FC)
    conc_needed <- 0
    
    outline <- data.frame(Protein_ID = yolk_id, StartConc = start_conc,
                          Conc_Deg = conc_deg,
                          Conc_Synth = conc_needed,
                          Sequence = row$Sequence)
    
    return(outline) }) %>% bind_rows

  out_df <- rbind(yolk_final, nonyolk_final)
  row.names(out_df) <- NULL
  
  return(out_df) }

XLA_prot_concs <- modify_frog_abs(XLA_abs_conc)
head(XLA_prot_concs)
```

    ##                   Protein_ID   StartConc    Conc_Deg  Conc_Synth
    ## 1 XBmRNA35012|XBXL10_1g18988 0.235544077 0.020676666 0.000000000
    ## 2 XBmRNA35013|XBXL10_1g18989 0.332823344 0.031244357 0.000000000
    ## 3 XBmRNA35014|XBXL10_1g18990 0.471484027 0.025625924 0.000000000
    ## 4 XBmRNA39801|XBXL10_1g21498 0.074227212 0.007298034 0.000000000
    ## 5 XBmRNA67639|XBXL10_1g35803 1.100507096 0.195450749 0.000000000
    ## 6                 XBgroup104 0.003225629 0.000000000 0.001879298
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Sequence
    ## 1           MKGIVLALLLALAGSERTHIEPVFSESKISVYNYEAVILNGFPESGLSRAGIKINCKVEISAYAQRSYFLKIQSPEIKEYNGVWPKDPFTRSSKLTQALAEQLTKPARFEYSNGRVGDIFVADDVSDTVANIYRGILNLLQVTIKKSQDVYDLQESSVGGICHTRYVIQEDKRGDQIRIIKSTDFNNCQDKVSKTIGLELAEFCHSCKQLNRVIQGAATYTYKLKGRDQGTVIMEVTARQVLQVTPFAERHGAATMESRQVLAWVGSKSGQLTPPQIQLKNRGNLHYQFASELHQMPIHLMKTKSPEAQAVEVLQHLVQDTQQHIREDAPAKFLQLVQLLRASNFENLQALWKQFAQRTQYRRCLLDALPMAGTVDCLKFIKQLIHNEELTTQEAAVLITFAMRSARPGQRNFQISADLVQDSKVQKYSTVHKAAILAYGTMVRRYCDQLSSCPEHALEPLHELAAEAANKGHYEDIALALKALGNAGQPESIKRIQKFLPGFSSSADQLPVRIQTDAVMALRNIAKEDPRKVQEILLQIFMDRDVRTEVRMMACLALFETRPGLATVTAIANVAARESKTNLQLASFTFSQMKALSKSSVPHLEPLAAACSVALKILNPSLDNLGYRYSKVMRVDTFKYNLMAGAAAKVFIMNSANTMFPVFILAKFREYTSLVENDDIEIGIRGEGIEEFLRKQNIQFANFPMRKKISQIVKSLLGFKGLPSQVPLISGYIKLFGQEIAFTELNKEVIQNTIQALNQPAERHTMIRNVLNKLLNGVVGQYARRWMTWEYRHIIPTTVGLPAELSLYQSAIVHAAVNSDVKVKPTPSGDFSAAQLLESQIQLNGEVKPSVLVHTVATMGINSPLFQAGIEFHGKVHAHLPAKFTAFLDMKDRNFKIETPPFQQENHLVEIRAQTFAFTRNIADLDSARKTLVVPRNNEQNILKKHFETTGRTSAEGASMMEDSSEMGPKKYSAEPGHHQYAPNINSYDACTKFSKAGVHLCIQCKTHNAASRRNTIFYQAVGEHDFKLTMKPAHTEGAIEKLQLEITAGPKAASKIMGLVEVEGTEGEPMDETAVTKRLKMILGIDESRKDTNETALYRSKQKKKNKIHNRRLDAEVVEARKQQSSLSSSSSSSSSSSSSSSSSSSSSSSSSPSSSSSSSYSKRSKRREHNPHHQRESSSSSSQEQNKKRNLQENRKHGQKGMSSSSSSSSSSSSSSSSSSSSSSSSSSSSEENRPHKNRQHDNKQAKMQSNQHQQKKNKFSESSSSSSSSSSSEMWNKKKHHRNFYDLNFRRTARTKGTEHRGSRLSSSSESSSSSSESAYRHKAKFLGDKEPPVLVVTFKAVRNDNTKQGYQMVVYQEYHSSKQQIQAYVMDISKTRWAACFDAVVVNPHEAQASLKWGQNCQDYKINMKAETGNFGNQPALRVTANWPKIPSKWKSTGKVVGEYVPGAMYMMGFQGEYKRNSQRQVKLVFALSSPRTCDVVIRIPRLTVYYRALRLPVPIPVGHHAKENVLQTPTWNIFAEAPKLIMDSIQGECKVAQDQITTFNGVDLASALPENCYNVLAQDCSPEMKFMVLMRNSKESPNHKDINVKLGEYDIDMYYSADAFKMKINNLEVSEEHLPYKSFNYPTVEIKKKGNGVSLSASEYGIDSLDYDGLTFKFRPTIWMKGKTCGICGHNDDESEKELQMPDGSVAKDQMRFIHSWILPAESCSEGCNLKHTLVKLEKAIATDGAKAKCYSVQPVLRCAKGCSPVKTVEVSTGFHCLPSDVSLDLPEGQIRLEKSEDFSEKVEAHTACSCETSPCAA
    ## 2  MRGIILALLLAIAGSERTQIEPVFSESKTSVYNYEAVILNGFPESGLSRAGIKINCKVEISAYAQRSYFLKIQSPELKEYNGVWPKDPFTRSSKLTQALAEQLTKPARFQYSNGRVGDIFVSDDVSDTVLNIYRGILNLLQVTIKKSQNVYDLQEPSIGGICHTRYVLQEDSRGDRISIIKATDFNNCHEKVFKGIGFELTESCESCKQFHRNIRGTATYTYKLKGRDQGSVIMEVTARQVLQFTPFAERNGAATMESRQILVWTGSKSGPIHPPQIQLKNRGNLHYQFASELHQMPFHLMKTKSPEAQAVEVLQHLMKDTQQQIREDTPVKFLQLIQLLRSSDFEDLQALWKQFAQRTQYRRWLLDAIPMAGTVDCLKFVKQLIHNEELNPYEAAVTITLALRSARPSQRAFQISTDFVQDSKVQKYSTVHKAAILAYGTMVKKYCDQLSSCPEQALEPLHDLAAEAADKGHAEIIALALKALGNAGQRESIKRIQKFLPGFSSSAYQLPVRIQTDAVMALRNIAKQDPQKVQEILLQIFMDRDVRTEVRMMACLALFETKPHLATVTTIANVAARESKTNLQLASFTYSQLKALAKSSVPHLEPLAAACSVALKILNPSLDNLGYRYSKVMRVDTFKYNLMAGGAAKVYVMNSANTMFPVFILAKFREYLAGVESDILEIGIRGEGIEEILRKQNIQFADFPMRKKISQILKTLLGFKGLPSQMPLISGYVKVLGQEIAFTEVNKEVIQNIIRALNEPAERHTMIRNILNKFLGGVAGHYSQALMTGEYRYLIPTTVGLPAEFSLYHSAIVNAAVNSDIKVKPSPSKDFSVAQLMESQIQLNADVNPSFMFYKVGTMGINSRLIQAGFEFHGKVHARLPAKFTAYLDMKDRNFKIETPPCQQENHLVELRAQTFAVTRNIEDLDAARKNLVLPRNSEQNILKKQFQSTGRSSAEGASMMEDSSEMGPRKYSTEPGHPQFAPNINSYDACTKLSKAGVHFCIQCKTHNAASHRNTFLYHIIGEHEFKLTMKPAHTEGAIEKLQLEIAAGPKAASKLMGLIEVEETEGEKMGDSAVQKRLKLILGIDDSVMNTNETAVLRSKKSKKERKEHKRPHDAEVVEAKRQQSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSSSSSSYRDRPNRRRQKSGQQGESSSSSSQKQDKHKPQENRKHGQKAQSSSSSSSSSSSSSSSSSSSSSSRSSSSSSSESSPSRSHNKQHENEQRKKQPKQNQQNQQNKNKYGSSSSSSGSSSSSSSSEMWTKKKHHRHFYDLHFRRAARTKDKKQRGSQWSSSSESSSSSSESANTNRRKTNFLGDKESPVFVATFRAVLNDNTKQGYQMVVYQDQHSSKQQIQAYVMDISKSRWATCLNAVVNNQYEAQASLKWGQNCQDYKINFKAQTGNFGNQPALKVTANWPKIPSKLKSTGKYVAEYVPGAMYMMGFQGEYKRNSQRQVKLVFALSSPRTCDVVIQTPRVTVYYRALRLPVPIPVGHHTKENVLQTPTWNIFAEAPQIIMDSLQGECKVAQDQITTFNGVDLVSALPENCYHVLAQDCSPEMKFMVLMRNSKESPNHKAINVKLGTYDIDMHYSSEMLKMKINGIELSEERLPYKSFEDPTVELKKKGNGVSLSAPEYGINSLDYDGLTFKLKAAIWMKGKTCGICGHNDDENEKELQMPDGSVAKGPMSFIHSWILPAESCSEGCNLKRSLVQLEKEIDGAMAKCYSVQPVLRCAKGCSPLKTEKVSTGFHCLPTDASLDLPTSQTRLEKSEDFSESVDAHTACSCETSQCAA
    ## 3 MRGIILALLLALAGSEKSQYEPFFSESKTYVYNYEGIILNGIPENGLARSGIKLNCKAEISGYAQRSYMLKIRDPEIKEFNGLWPKDPFTRATKLSQILAEHLTRPVKFEYRNGQVGNIFAAEDVSETALNIQRGILNLLQLTIKKSQNVYDLQEDSISGVCHTKYTIQEDKKAERLTVMKSRDFDRCESRVEKKIGVAYLETCTQCKKKVQHLRGSATFTYKLQNKDQGALIVEASAKQVYRYAPYYELYGAAVMEARQILTLAETKSNRVSESQVQLTKRGGLNYHFWQEEHQMPFHLLKTKNAESQVVETLQHLVQNNQELVHGESTSKFFELVQLLRTTNQENIQAIWKQFSDRPQHRHWLVNAIPMAGSVDALKFIKQMIHNQDMTHQEAATAILMAMHYTRASQRAVEIAADLVTDSRVKNAPSLYKSALLAYGSLVYKYCVDSQSCPEAALQPLHEFAAEAANNAHQEEITLALKAIGNAGQPASIKRIQKFLPGFSSGASRLPVRLQTDAVMALRNIAIKDPRKVQEILLQLFMDRDVHPEVRMITCVALFETNPGIAIVTAVANAVLRESKDNMQLASFTYSHIKSLSKSMMPEFYNLASACRVALNIMNPVFDQLSFRYSKAFHVDTFNYPMMGGLALNVFFINSPNTAFPVFALIKVRECSFGASTDVFEVALRAEGLQEVIRKQNIQFSEFSMHKKIVKILKALMAYKEAPANVPLVSASIKLLGQEVNFNQATRDTIQNAMKLLNEPAERHTMVKKILNKLLNGVTGQWSHPMLVAEARHIVPTCVGLPLEIATHSAAVANAALNVDMKVSPSPSNDFSLSQLLESNIQLHTDMSPSVLIHLVNTMGVNTPLIQSGLEFHGKVRIQTPAKFTAKLDMKDKNFKIETPPCQREVELVTVRHQMFAVTRNVGQQDSEKRMLVMPKQHGPSTMEQHFEKSGRTSAESASMMEVSSEMGIRNKNAQHQHILNNVESYRYCFRLSKLNINACLQSKLNTVQGDRRSFLYENIGEHEIKFLLKPTQNEPAIEKLQVEITAGPRASSKIIALAEVEEKEGENLDESAVQKRLKIILGINEEEYSARNSSLVRKASKKRKVKKTKDAGETDVLEAERKWQQSSSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSKSSSSSSSSSPQRNSRNKGIKVDVLSRVSKRQQEKNKDTHGQTGRNKHRASSSESSSSESSSSESSRSSESSSSSSSESSSSSESSRRHGKSGKQQQGDREQQREKQTHQNKQKQQKHRSRHPQDPSQMPPQLYKLKFQRPYKQQYKGMSGQSESSESSSSSSSSSSSSESARSPHHKKFLGDKKSPTVVITAQAVRNDNVKQGYQVILYTDDTTRPKMQAYLVDITGKSRWRACAEAFVRNPHEAQAELKWGQNCKQYKMEFIAQTGIFGNQPALRLRVNWPRIPSKLKSVGEEVGEFLPGAAYMMGLNTKFQRNPSKQAKFIFSLSSPRTMNAVIQVPWATCYYQAINLPFAVTAAYHLHSSDIQAPTWNIFANSYSSMVDRLSGECTMSSDQIKTFNEVKFNYTLPGSCFHVLAQDCSTELKFLVLAKKSPESSELRDINIKLADYDIDMYTSSEEIRLKVNGLEVPKEKLPYIASTNPIVEIKMEEGALVLTAPVFGIEQLYFNGKLTKLRTAVWMRGGTCGICGRHDDETEREYQMPDGYAARDENSFTHSWILPEEPCSGTCKLQHTSVKLEKTIHIDDQKSKCQSVRPVLRCVKGCSPTSVTPVTVGFHCLPFDASTDMLESQNIMGKTEDLVETVDAHTACSCEALKCAA
    ## 4                  MRGIILALLLALAGCEKSQYEPFFSESKTYVYNYEGIILNGIPENGLARSGIKLNCKVELSGYAQRSYMLKIRNPEIKEFNGLWPKDPFTRATKLSQILAEQLTRPVKFEYRNGQVGNIFAPEDVSETVLNIQRGILNMLQVTIKKSQNVYDLQEDSISGVCHTKYTVQEDKKAERLIVIKSRDFENCESRVEKKIGVAYMETCTQCKKRAKHLQGSATFTYKLQNKDQGALIVEASAKQVYRYTPYYELYGAAVMEARQILTLAETKSSRNSESQVQLTKRGGLNYHFWQEVHQMPFHLLKTKNAESQVVETLQHLVQNNQELVHGESTSKFFQLVQLLRATDQENIQAIWKQFSDRPQHRHWLVNAITMAGSVDALKFIKQMIHNQDMTHQEAATAILMAMHYTRASQRVVEIAADLLTDSKVQNSPTLYKSAMLAYGSLVYKYCVDSQSCPEAALQPLHEFAAEAANKAHQEDITLALKAIGNAGQPASIKRIQKFLPGFSSGASRFPVRVQTDAVMALRNIAKKDPRKVQEVLLQLFMDRDVHPEVRMITCVALFETNPGIAIVTAVANAVLRESKDNLQLASFTYSHIKSLTKSTMPEFHNLASACRVALNILNPVFDQLSFRYSKAFHVDTFSYPMMGGLALNVFFINSPRTAFPVFALLKVRESFAGASTDVLEVALRAEGLQEVIRKQNIQFSEYSMRKKIVKILKALLAYKEAPANVPLVSASIKMLGQEVNFNQVTRDTIQNAMKLLNEPAERHTVVKRVLNKLLSGATGQWSHPMLLTEARHIIPTCVGLPLEIATHSAAVATTALNVDMKVTPSPSNDFSLSQLLESNIQLHTDMSPSVLIHVVETMGVNTPLIQSGLEFHGKLRIQTPVKFTAKLDMKDKNFKMEFPPCQREVEILTARHQMFAVTRNVGQLDSEKRMLVVPKQHGPSTLEQHFEKSGRTSAESASMMEVSSEMGIRKNAYSGGGAQHQRIPNHVQSYHFCLRLSKLNINACLQSKLNTAHGDMQSFLYEQIGEHETKVILKPVQNDAAIEKLQLEITTGPRASSKIIALTEVEVKEGEQLDESAVQKRLKIILGINEEEYRARNRSLVGKESKKAKIKKNKGAGETDVLEAGRNWKQSSSSSSSSPSSSSSSSSSSSSSSSSSSSSSSSSESSSSSSQRNRRNHGRKVDILSRVSQQAMNKETHGQTGRQKQRDSSSSSSSSQSSNKESSSSSSSSSSESSSSSESSRGKHQQNSKQQTRQNKQKQQKHHGSHPQDPSQMPPQMYKLRFQRPYKQQYKGMSSQSSESSSSSSSSSSSESARSQHHNKFLGDKKPPTVVITAQAVRNDNVKQGYQVILYTDDTTRPKMQAYVVDITGKSRWRACAEAYIRNPHEAQAELKWGQNCKQYKMEFIAQTGIFGNQPALRFKMNWPRIPSKLKSAGEEVAEFLPGAAYMMGFNTKNQRNPSKQAKIILSLSSPRTMNAVIQVPWATCYYQAVNLPFAVTSAYHLHSSDIQAPTWNIFADTYSSLVDRFSGECTMSRDQIKTFNEVKFNYTLPGSCFHVLAQDCSSELKFLVMAKKSPESSELRDINIKLADYDIDMYTSSEEIRLKVNGLEVPKEKLPYIASADPAVEIKMEEKAVILSAPDFGIEQLYYNGKLTKVRTAVWMRGSTCGICGQHDDETEREYQMPGGHAARDENRFTHSWILPEEPCSGTCKLQHTSVKLERTVHIDDQESKCYSVRPVLRCVKGCSPTSVTPVTVGFHCLPADASTDLVESQNIMGKTEDLVDTVDAHTACSCEALKCAA
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      MHLLVYLSLFFALALASVTEISLDNKHRHRHEQQGHHDSAKHGHQKDKQQQEQIKNDEGKLTKEEKILSEENSDFSVNLFNQLSTESKRSPRKNIFFSPISISAAFYMLALGAKSETHQQILKGLSFNKKKLSESQVHEAFKRLIEDSNNPMKAHQFTIGNALFVEQTVNILKGFEENVKHYYQAGVFPMNFKDPDNAKKQLNNYVKDKTHGVIQEMIRELDSNTEMVLVNYVLYKGEWANNFNPTLTQKSLFSVDKNTNVTVQMMNRLGLYRTYQDDDCKIIELPYKNDTAMLLVVPQLGKIQELVLTSKLINHWYESLATSIVDLYMPTFSISGKVVLKDTLRKMGISDIFTDKADLTGISEQIKLKVSMASHNAVLNVNEFGTEAVGATSAQASPTKLFPPFLIDSPFLVMIYSRTLGSQLFMGKVMDPTNAQ
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       MFSTSAKIVKPNGEKPDEFESGISQALLELEMNSDLKAQLRELNITAAKEIEVGAGRKAIIIFVPVPQLKSFQKIQVRLVRELEKKFSGKHVVFIAQRRILPKPTRKSRTKNKQKRPRSRTLTAVHDAILEDLVYPSEIVGRRIRVKLDGSRLIKVHLDKAQQNNVEHKVVETFSGVYKKLTGKDVVFEFPEFQL

``` r
modify_fly_abs <- function(df){

  fits_df <- Fly_gast_fits[Fly_gast_fits$Protein_ID %in% df$Ind_Protein_ID,]

  fits_df["lightFC"] <- fits_df$N16 / fits_df$T0_L
  fits_df["heavyFC"] <- fits_df$O12 / fits_df$T0_H
  fits_df[fits_df$heavyFC > 1, "heavyFC"] <- 1
    
  fits_df["AbsCategory"] <- rep("CHANGE", nrow(fits_df))
  fits_df[fits_df$lightFC > 1, "AbsCategory"] <- "Inc"
  fits_df[fits_df$lightFC < 1, "AbsCategory"] <- "Dec"

  final_df <- lapply(df$Ind_Protein_ID, function(x){

    fit_row <- fits_df[fits_df$Protein_ID == x,]
    abs_row <- df[df$Ind_Protein_ID == x,]
    start_conc <- abs_row$iBAQ_Conc_Start

    T_end <- tail(Fly_light_time, 1)

    if (fit_row$mClass == "Deg") { surv <- exp(-fit_row$kD * T_end) }
    else { surv <- 1 }
    
    conc_deg    <- start_conc * (1 - surv)
    conc_needed <- max(start_conc * fit_row$lightFC - start_conc * surv, 0)
    
    outline <- data.frame(Protein_ID = x, StartConc = start_conc,
                          Conc_Deg = conc_deg,
                          Conc_Synth = conc_needed,
                          Sequence = abs_row$Sequence)

    return(outline) }) %>% bind_rows

  final_df <- final_df %>% arrange(desc(StartConc))
  
  return(final_df) }

Dmel_prot_concs <- modify_fly_abs(Dmel_abs_conc)
head(Dmel_prot_concs)
```

    ##   Protein_ID  StartConc    Conc_Deg   Conc_Synth
    ## 1     P02843 0.09930386 0.017996817 0.0000000000
    ## 2     P06607 0.05240221 0.010267455 0.0008521061
    ## 3     P02844 0.04581583 0.009591635 0.0009946195
    ## 4     P08736 0.03983611 0.014140260 0.0123430050
    ## 5     Q9W334 0.02782544 0.006610322 0.0035318179
    ## 6     Q24560 0.02518482 0.009813133 0.0069138194
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Sequence
    ## 1                         MNPMRVLSLLACLAVAALAKPNGRMDNSVNQALKPSQWLSGSQLEAIPALDDFTIERLENMNLERGAELLQQVYHLSQIHHNVEPNYVPSGIQVYVPKPNGDKTVAPLNEMIQRLKQKQNFGEDEVTIIVTGLPQTSETVKKATRKLVQAYMQRYNLQQQRQHGKNGNQDYQDQSNEQRKNQRTSSEEDYSEEVKNAKTQSGDIIVIDLGSKLNTYERYAMLDIEKTGAKIGKWIVQMVNELDMPFDTIHLIGQNVGAHVAGAAAQEFTRLTGHKLRRVTGLDPSKIVAKSKNTLTGLARGDAEFVDAIHTSVYGMGTPIRSGDVDFYPNGPAAGVPGASNVVEAAMRATRYFAESVRPGNERSFPAVPANSLQQYKQNDGFGKRAYMGIDTAHDLEGDYILQVNPKSPFGRNAPAQKQSSYHGVHQAWNTNQDSKDYQ
    ## 2                                            MMSLRICLLATCLLVAAHASKDASNDRLKPTKWLTATELENVPSLNDITWERLENQPLEQGAKVIEKIYHVGQIKHDLTPSFVPSPSNVPVWIIKSNGQKVECKLNNYVETAKAQPGFGEDEVTIVLTGLPKTSPAQQKAMRRLIQAYVQKYNLQQLQKNAQEQQQQLKSSDYDYTSSEEAADQWKSAKAASGDLIIIDLGSTLTNFKRYAMLDVLNTGAMIGQTLIDLTNKGVPQEIIHLIGQGISAHVAGAAGNKYTAQTGHKLRRITGLDPAKVLSKRPQILGGLSRGDADFVDAIHTSTFAMGTPIRCGDVDFYPNGPSTGVPGSENVIEAVARATRYFAESVRPGSERNFPAVPANSLKQYKEQDGFGKRAYMGLQIDYDLRGDYILEVNAKSPFGQRSPAHKQAAYHGMHHAQN
    ## 3                      MNPLRTLCVMACLLAVAMGNPQSGNRSGRRSNSLDNVEQPSNWVNPREVEELPNLKEVTLKKLQEMSLEEGATLLDKLYHLSQFNHVFKPDYTPEPSQIRGYIVGERGQKIEFNLNTLVEKVKRQQKFGDDEVTIFIQGLPETNTQVQKATRKLVQAYQQRYNLQPYETTDYSNEEQSQRSSSEEQQTQRRKQNGEQDDTKTGDLIVIQLGNAIEDFEQYATLNIERLGEIIGNRLVELTNTVNVPQEIIHLIGSGPAAHVAGVAGRQFTRQTGHKLRRITALDPTKIYGKPEERLTGLARGDADFVDAIHTSAYGMGTSQRLANVDFFPNGPSTGVPGADNVVEATMRATRYFAESVRPGNERNFPSVAASSYQEYKQNKGYGKRGYMGIATDFDLQGDYILQVNSKSPFGRSTPAQKQTGYHQVHQPWRQSSSNQGSRRQ
    ## 4 MGKEKIHINIVVIGHVDSGKSTTTGHLIYKCGGIDKRTIEKFEKEAQEMGKGSFKYAWVLDKLKAERERGITIDIALWKFETAKYYVTIIDAPGHRDFIKNMITGTSQADCAVLIVAAGTGEFEAGISKNGQTREHALLAFTLGVKQLIVGVNKMDSSEPPYSEARYEEIKKEVSSYIKKIGYNPAAVAFVPISGWHGDNMLEPSTNMPWFKGWKVERKEGNADGKTLIDALDAILPPARPTDKALRLPLQDVYKIGGIGTVPVGRVETGVLKPGTVVVFAPANITTEVKSVEMHHEALQEAVPGDNVGFNVKNVSVKELRRGYVAGDSKANPPKGAADFTAQVIVLNHPGQIANGYTPVLDCHTAHIACKFAEIKEKVDRRSGKTTEENPKFIKSGDAAIVNLVPSKPLCVEAFQEFPPLGRFAVRDMRQTVAVGVIKAVNFKDASGGKVTKAAEKATKGKK
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                               MDKPVVWARVMKVLGRTGSQGQCTQVKVEFLGEQNRQIIRNVKGPVREGDILTLLESEREARRLR
    ## 6                 MREIVHIQAGQCGNQIGAKFWEIISDEHGIDATGAYHGDSDLQLERINVYYNEASGGKYVPRAVLVDLEPGTMDSVRSGPFGQIFRPDNFVFGQSGAGNNWAKGHYTEGAELVDSVLDVVRKEAESCDCLQGFQLTHSLGGGTGSGMGTLLISKIREEYPDRIMNTYSVVPSPKVSDTVVEPYNATLSVHQLVENTDETYCIDNEALYDICFRTLKLTTPTYGDLNHLVSLTMSGVTTCLRFPGQLNADLRKLAVNMVPFPRLHFFMPGFAPLTSRGSQQYRALTVPELTQQMFDAKNMMAACDPRHGRYLTVAAIFRGRMSMKEVDEQMLNIQNKNSSYFVEWIPNNVKTAVCDIPPRGLKMSATFIGNSTAIQELFKRISEQFTAMFRRKAFLHWYTGEGMDEMEFTEAESNMNDLVSEYQQYQEATADEDAEFEEEQEAEVDEN

``` r
calculate_AAs <- function(row) {

  #Splits sequence into individual characters
  chars <- strsplit(row["Sequence"], split = "")[[1]]
  chars <- chars[chars %in% aa_vector]

  #Setting order 
  counts <- table(factor(chars, levels = aa_vector))

  AA_needed <- counts * as.numeric(row["Conc_Synth"])
  names(AA_needed) <- paste(names(AA_needed), "Need", sep="_")
  
  AA_gen <- counts * as.numeric(row["Conc_Deg"])
  names(AA_gen) <- paste(names(AA_gen), "Gen", sep="_")
    
  outline <- data.frame(t(c(AA_needed, AA_gen)))
  outline <- cbind(data.frame(Protein_ID=as.character(row["Protein_ID"])),
                   outline)
    
  return(outline) }

XLA_AA_tracking <- apply(XLA_prot_concs, 1, calculate_AAs) %>% bind_rows
Dmel_AA_tracking <- apply(Dmel_prot_concs, 1, calculate_AAs) %>% bind_rows

head(XLA_AA_tracking)
```

    ##                   Protein_ID     A_Need C_Need     D_Need     E_Need     F_Need
    ## 1 XBmRNA35012|XBXL10_1g18988 0.00000000      0 0.00000000 0.00000000 0.00000000
    ## 2 XBmRNA35013|XBXL10_1g18989 0.00000000      0 0.00000000 0.00000000 0.00000000
    ## 3 XBmRNA35014|XBXL10_1g18990 0.00000000      0 0.00000000 0.00000000 0.00000000
    ## 4 XBmRNA39801|XBXL10_1g21498 0.00000000      0 0.00000000 0.00000000 0.00000000
    ## 5 XBmRNA67639|XBXL10_1g35803 0.00000000      0 0.00000000 0.00000000 0.00000000
    ## 6                 XBgroup104 0.02067228      0 0.01315509 0.03006877 0.01879298
    ##       G_Need      H_Need     I_Need     K_Need     L_Need      M_Need
    ## 1 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.000000000
    ## 2 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.000000000
    ## 3 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.000000000
    ## 4 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.000000000
    ## 5 0.00000000 0.000000000 0.00000000 0.00000000 0.00000000 0.000000000
    ## 6 0.01691368 0.007517192 0.02631017 0.04322385 0.03382736 0.003758596
    ##       N_Need     P_Need     Q_Need     R_Need     S_Need     T_Need     V_Need
    ## 1 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    ## 2 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    ## 3 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    ## 4 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    ## 5 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
    ## 6 0.01127579 0.01691368 0.01879298 0.02818947 0.02255158 0.01503438 0.03570666
    ##   W_Need      Y_Need    A_Gen     C_Gen     D_Gen    E_Gen    F_Gen     G_Gen
    ## 1      0 0.000000000 3.018794 0.6409768 1.4266902 2.398494 1.343984 1.8195470
    ## 2      0 0.000000000 4.311722 0.9685752 2.1558608 3.499368 2.030883 2.8744811
    ## 3      0 0.000000000 3.587629 0.7944035 1.5631811 3.254492 1.614433 2.0500736
    ## 4      0 0.000000000 1.021725 0.2262391 0.4597761 0.861168 0.430584 0.5984388
    ## 5      0 0.000000000 4.690817 0.1954507 4.2999154 5.081718 4.495366 3.7135633
    ## 6      0 0.003758596 0.000000 0.0000000 0.0000000 0.000000 0.000000 0.0000000
    ##       H_Gen     I_Gen     K_Gen    L_Gen     M_Gen     N_Gen     P_Gen    Q_Gen
    ## 1 1.0958635 1.9436070 2.8740571 2.874057 0.9718035 1.7781936 1.5507503 2.171050
    ## 2 1.4684849 2.9682142 4.4054548 4.467943 1.4997293 2.5932819 2.4058157 3.624346
    ## 3 1.2812960 2.3575846 3.4594992 3.767010 1.2812960 2.2550810 1.9219440 2.844477
    ## 4 0.3794978 0.6568231 0.9487444 1.050917 0.3649017 0.5984388 0.5400545 0.853870
    ## 5 2.9317605 4.6908168 7.6225773 8.990732 2.9317605 5.8635210 2.7363098 4.886267
    ## 6 0.0000000 0.0000000 0.0000000 0.000000 0.0000000 0.0000000 0.0000000 0.000000
    ##      R_Gen    S_Gen     T_Gen     V_Gen      W_Gen     Y_Gen
    ## 1 1.860900 4.135334 1.8195470 2.2951104 0.26879671 1.0751868
    ## 2 2.749504 6.811270 2.7182593 3.0931916 0.37493232 1.7184398
    ## 3 2.331959 5.381443 2.3575846 2.8444771 0.33313696 1.2812960
    ## 4 0.642227 1.474203 0.7079093 0.8465719 0.09487444 0.3794978
    ## 5 1.954507 6.840775 4.8862675 5.4726196 0.39090140 2.5408591
    ## 6 0.000000 0.000000 0.0000000 0.0000000 0.00000000 0.0000000

``` r
#Sum of all amino acids needed for new protein synthesis + generated from deg
frog_AA_sums <- apply(XLA_AA_tracking[2:41], 2, sum)

#Amino acids needed to synthesize new protein
frog_AA_needed <- frog_AA_sums[1:20]
names(frog_AA_needed) <- aa_vector

#Amino acids provided from degradation
frog_AA_deg <- frog_AA_sums[21:40]
names(frog_AA_deg) <- aa_vector

#Amino acids provided from yolk
frog_yolk_contribution <-
  apply(XLA_AA_tracking[XLA_AA_tracking$Protein_ID %in%
                           frog_yolk_set, 22:41], 2, sum)
names(frog_yolk_contribution) <- aa_vector

#Changing vector to only be amount provided by degradation w/out yolk
frog_AA_deg <- frog_AA_deg - frog_yolk_contribution
```

``` r
#Sum of all amino acids needed for new protein synthesis + generated from deg
fly_AA_sums <- apply(Dmel_AA_tracking[2:41], 2, sum)

#Amino acids needed to synthesize new protein
fly_AA_needed <- fly_AA_sums[1:20]
names(fly_AA_needed) <- aa_vector

#Amino acids provided from degradation
fly_AA_deg <- fly_AA_sums[21:40]
names(fly_AA_deg) <- aa_vector

#Amino acids provided from yolk
fly_yolk_contribution <-
  apply(Dmel_AA_tracking[Dmel_AA_tracking$Protein_ID %in%
                           fly_yolk_set, 22:41], 2, sum)
names(fly_yolk_contribution) <- aa_vector

#Changing vector to only be amount provided by degradation w/out yolk
fly_AA_deg <- fly_AA_deg - fly_yolk_contribution
```

Citations for frog nucleotides:

<https://www.pnas.org/doi/10.1073/pnas.51.1.139> Gurdon and WIckens
table

Justification:

To estimate the total metabolic cost of nucleic acids during our
timecourse, we derived the net biosynthesis requirements from
established literature. The Xenopus laevis embryo begins with roughly 12
pg of DNA, which expands to nearly 1 µg by late embryonic stages (Gurdon
and Wickens, 1989; Gurdon and Brown, 1963). Additionally, studies
utilizing anucleolate mutants that lack the ability to synthesize new
rRNA demonstrate that the wild-type embryo synthesizes approximately 6
µg of net rRNA by the early swimming stages (Brown & Gurdon, 1964).
Thus, the embryo requires the net synthesis of roughly 7 µg of the major
nucleic acids by NF35. We treat the net mass accumulation of minor RNA
species (mRNA, tRNA, etc) as negligible, as their combined net expansion
(\< 0.2 µg) is functionally subsumed by our upper-bound estimates for
DNA and rRNA, which were derived from NF41 embryos. Furthermore,
mitochondrial RNA remains static within our analysis window (Chase and
David 1972). Given an average elemental nitrogen fraction of 15.5% for
mixed biological nucleotides, the synthesis of 7 µg of net nucleic acid
requires almost 1.1 µg of atomic nitrogen.

Likewise, we repeated the same approximation in fly from established
literature for the net biosynthesis requirements. We begin our
experimental window at roughly 3 hours post-fertilization, a stage
comprising approximately 6,000 diploid nuclei (Zalokar and Erk, 1976).
Derived from a genome size of 180 megabases (Adams, et al. 2000), the
diploid genome mass contains 0.36 pg DNA which corresponds to 2.2 ng of
starting DNA. At the point of hatching, the embryo contains roughly
50,000 cells (Lehner et al., 2001). However, the fly embryo is known to
enter the endocycle, where specific tissues (salivary gland, midgut,
hindgut, Malpighian tubules, and dorsal cells) replicate their genomes
an average of two additional times without cleavage (Smith and
Orr-Weaver, 1991). Based on anatomical proportions, we estimate these
tissues comprise roughly 10% of the embryo. Thus, by hatching, the
embryo contains an approximate total of 23.9 ng of DNA, resulting in a
net synthesis of 21.7 ng during our experimental window. We treat the
net mass accumulation of bulk RNA species (rRNA, mRNA, tRNA, etc.) as
negligible, as they are known not to appreciably increase during this
developmental window (Anderson and Lengyel, 1979; Becker et al., 2018).
Utilizing the 16.8% mass fraction of elemental nitrogen for the
Drosophila genome, the synthesis of 21.8 ng of net DNA requires
approximately 3.6 ng of atomic nitrogen.

To allocate the atomic nitrogen demand, we modeled the free amino acid
pool using a mass balance approach. First, an intermediate pool was
established by summing yolk-derived and proteome-degraded amino acids,
then subtracting the demands of de novo protein synthesis. Assuming
rapid transaminase-driven equilibration, each amino acid’s contribution
was weighted by its intermediate molarity and atomic nitrogen count.
These proportional nitrogen demands were then converted back into amino
acid molarities, yielding the final physiological amino acid demands for
nucleic acid synthesis.

Math for fly nucleotides:

- Genome size is 180 megabases (Adams 2000) \> 978 megabases per
  picogram

- 6,000 2C nuclei ~ 3.1 hpf (Zalokar and Erk 1976) \>
  ((180x2/978)\*6000)/1000 \[1\] 2.208589

> 2.2 ng DNA at 3.1 hpf

- 50,000 cells at hatching (Lehner, et al. 2001)

- 5,000 endocycle based on anatomical proportions -\> 8C

- 45,000 normal -\> 2C

- ((180*2/978)*45000)/1000 -\> 16.6 ng from 2C cells

- ((180*8/978)*5000)/1000 -\> 7.4 ng from endocycle cells

> 23.9 ng DNA at hatching

- Thus, 21.7 ng DNA made during experimental window.

Atomic nitrogen: \>
((((180*2/978)*45000)/1000+((180*8/978)*5000)/1000)-(((180*2/978)*6000)/1000))\*.168
\[1\] 3.648589

``` r
aa_profiles <- read.csv("Files/Reference/AA_Profiles.csv")

aa_nitrogens <- aa_profiles$Nitrogen
names(aa_nitrogens) <- aa_profiles$AA
aa_nitrogens <- aa_nitrogens[aa_vector]

head(aa_nitrogens)
```

    ## A C D E F G 
    ## 1 1 1 1 1 1

``` r
#Embryo needs 1 ug of nitrogen
#Atomic nitrogen ~14.007 ug/umol
#Cell volume is 1 uL
frog_nuc_acid_nitrogen <- (((7*0.155)/14.007)/1)*1000

#Amount leftover after protein synthesis is taken into account
frog_AA_remaining <- (frog_yolk_contribution + frog_AA_deg) - frog_AA_needed

#Mass action -> Each AA donates based on no. of nitrogens
frog_N_pool <- frog_AA_remaining * aa_nitrogens
frog_N_pool <- frog_N_pool / sum(frog_N_pool)

#AAs used for nucleic acid synthesis
frog_AA_nuc_synthesis <- frog_N_pool * frog_nuc_acid_nitrogen
frog_AA_nuc_synthesis <- frog_AA_nuc_synthesis / aa_nitrogens

final_frog_AA_leftover <- frog_AA_remaining - frog_AA_nuc_synthesis

sum(final_frog_AA_leftover * aa_nitrogens)
```

    ## [1] 53.21301

``` r
fly_nuc_acid_nitrogen <- ((((180*2/978)*45000)/1000+((180*8/978)*5000)/1000)-
                            (((180*2/978)*6000)/1000))*.168
#Convert to grams and divide by g/mol
fly_nuc_acid_nitrogen <- (fly_nuc_acid_nitrogen*10e-10)/(14.01)
#Divide by 10nL (in L) to get Molar
fly_nuc_acid_nitrogen <- fly_nuc_acid_nitrogen / (1e-8)
fly_nuc_acid_nitrogen <- fly_nuc_acid_nitrogen * 1000

#Amount leftover after protein synthesis is taken into account
fly_AA_remaining <- (fly_yolk_contribution + fly_AA_deg) - fly_AA_needed

#Mass action -> Each AA donates based on no. of nitrogens
fly_N_pool <- fly_AA_remaining * aa_nitrogens
fly_N_pool <- fly_N_pool / sum(fly_N_pool)

#AAs used for nucleic acid synthesis
fly_AA_nuc_synthesis <- fly_N_pool * fly_nuc_acid_nitrogen
fly_AA_nuc_synthesis <- fly_AA_nuc_synthesis / aa_nitrogens

final_fly_AA_leftover <- fly_AA_remaining - fly_AA_nuc_synthesis

sum(final_fly_AA_leftover * aa_nitrogens)
```

    ## [1] 49.03092

``` r
# your vectors: frog_AA_needed, frog_AA_deg, frog_AA_yolk
# your essential list:
# emb_essential_AAs <- c("F","W","I","L","V","K","H","T","M","C","Y")

# ---- order: essentials first, then the rest (keep remaining in alphabetical order) ----
aa_all <- sort(unique(names(frog_AA_needed)))               # assumes all three share names
aa_levels <- c(emb_essential_AAs, setdiff(aa_all, emb_essential_AAs))

# data for the stacked supply bar (Deg + Yolk)
df_supply <- tibble(
  AA = factor(aa_levels, levels = aa_levels),
  Deg  = as.numeric(frog_AA_deg[aa_levels]),
  Yolk = as.numeric(frog_yolk_contribution[aa_levels])) |>
  pivot_longer(c(Deg, Yolk),
               names_to = "component",
               values_to = "value")

# data for the single needed bar, plus essential flag
df_needed <- tibble(
  AA = factor(aa_levels, levels = aa_levels),
  Synth = as.numeric(frog_AA_needed[aa_levels]),
  Nuc = as.numeric(frog_AA_nuc_synthesis[aa_levels]),
  essential = if_else(AA %in% emb_essential_AAs, "Essential", "Other")) |>
  pivot_longer(c(Synth, Nuc),
               names_to = "component",
               values_to = "value")

# numeric x positions to separate the two bars per AA
df_supply$x <- as.numeric(df_supply$AA) + 0.22
df_needed$x <- as.numeric(df_needed$AA) - 0.22

p1 <- ggplot() +
  # ----- stacked Deg+Yolk bar -----
  geom_col(data = df_supply,
           aes(x = x, y = value, fill = component),
           width = 0.44, color="black") +
  scale_fill_manual(values = c(Deg = "#0072B2", Yolk = "#F0E442"), name = NULL) +

  # start a fresh fill scale so the next layer's colors don't clash
  ggnewscale::new_scale_fill() +

  # ----- single Needed bar with essential-specific colors -----
  geom_col(data = df_needed,
           aes(x = x, y = value, fill = component),
           width = 0.44, color="black") +
  scale_fill_manual(values = c(Synth = "grey90", Nuc = "grey60"),
                    name = NULL) +

  # ----- axes & theme -----
  scale_x_continuous(breaks = seq_along(aa_levels), labels = aa_levels) +
  labs(x = "Amino acid (AA)", y = "Concentration (mM)") +
  theme_bw() +  #minimal(base_size = 12)
  theme(axis.text=element_text(size=24,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=24),
        panel.border = element_rect(linewidth=2),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)) 
  # coord_cartesian(ylim=c(0,25))
  # scale_y_continuous(position = "right")

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=7.5, height=4.5, res=300)

p1
```

![](figures/AA-Tracking_Fly_and_Frog/Frog%20supply%20and%20demand%20barplot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
# your vectors: frog_AA_needed, frog_AA_deg, frog_AA_yolk
# your essential list:
# emb_essential_AAs <- c("F","W","I","L","V","K","H","T","M","C","Y")

# ---- order: essentials first, then the rest (keep remaining in alphabetical order) ----
aa_all <- sort(unique(names(fly_AA_needed)))               # assumes all three share names
aa_levels <- c(emb_essential_AAs, setdiff(aa_all, emb_essential_AAs))

# data for the stacked supply bar (Deg + Yolk)
df_supply <- tibble(
  AA = factor(aa_levels, levels = aa_levels),
  Deg  = as.numeric(fly_AA_deg[aa_levels]),
  Yolk = as.numeric(fly_yolk_contribution[aa_levels])) |>
  pivot_longer(c(Deg, Yolk),
               names_to = "component",
               values_to = "value")

# data for the single needed bar, plus essential flag
df_needed <- tibble(
  AA = factor(aa_levels, levels = aa_levels),
  Synth = as.numeric(fly_AA_needed[aa_levels]),
  Nuc = as.numeric(fly_AA_nuc_synthesis[aa_levels]),
  essential = if_else(AA %in% emb_essential_AAs, "Essential", "Other")) |>
  pivot_longer(c(Synth, Nuc),
               names_to = "component",
               values_to = "value")

# numeric x positions to separate the two bars per AA
df_supply$x <- as.numeric(df_supply$AA) + 0.22
df_needed$x <- as.numeric(df_needed$AA) - 0.22

p1 <- ggplot() +
  # ----- stacked Deg+Yolk bar -----
  geom_col(data = df_supply,
           aes(x = x, y = value, fill = component),
           width = 0.44, color="black") +
  scale_fill_manual(values = c(Deg = "#0072B2", Yolk = "#F0E442"), name = NULL) +

  # start a fresh fill scale so the next layer's colors don't clash
  ggnewscale::new_scale_fill() +

  # ----- single Needed bar with essential-specific colors -----
  geom_col(data = df_needed,
           aes(x = x, y = value, fill = component),
           width = 0.44, color="black") +
  scale_fill_manual(values = c(Synth = "grey90", Nuc = "grey60"),
                    name = NULL) +

  # ----- axes & theme -----
  scale_x_continuous(breaks = seq_along(aa_levels), labels = aa_levels) +
  labs(x = "Amino acid (AA)", y = "Concentration (mM)") +
  theme_bw() +  #minimal(base_size = 12)
  theme(axis.text=element_text(size=24,colour="black"),
        plot.title = element_text(size=24, hjust=0.5),
        axis.title = element_text(size=24),
        panel.border = element_rect(linewidth=2),
        legend.position = "none",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
        panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)) +
  coord_cartesian(ylim=c(0,25)) +
  scale_y_continuous(position = "right")


# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=7.5, height=4.5, res=300)

p1
```

![](figures/AA-Tracking_Fly_and_Frog/Fly%20supply%20and%20demand%20barplot-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

Mesh plot as concentration

``` r
mesh_calculator <- function(row, light_cols, heavy_cols, l_time, h_time) {
  
  start_conc <- as.numeric(row["StartConc"])
  conc_deg <- as.numeric(row["Conc_Deg"])
  conc_synth <- as.numeric(row["Conc_Synth"])  

  light_trend <- as.numeric(row[light_cols])
  light_FC <- light_trend / light_trend[1]
  kd <- as.numeric(row["kD"])
  start_conc <- as.numeric(row["StartConc"])

  #-------------------------------------------------
  
  if (row["mClass"] == "Deg") { 
    surv_L  <- exp(-kd * l_time)    
    surv_H  <- exp(-kd * h_time) }
  else {
    surv_L <- rep(1, length(l_time))    
    surv_H <- rep(1, length(h_time)) }
  
  protein_degraded <- start_conc * surv_H  

  #-------------------------------------------------

  protein_made <- pmax(start_conc * light_FC - start_conc * surv_L, 0)

  outline <- data.frame(t(c(protein_made, protein_degraded)))
  colnames(outline) <- c(names(row[light_cols]), names(row[heavy_cols]))
  outline <- cbind(data.frame(Protein_ID = as.character(row["Protein_ID"])),
                   outline)
  
  return(outline) }


XLA_mesh_conc <- merge(XLA_prot_concs[1:4],
                       XLA_2C_fits[!XLA_2C_fits$Protein_ID %in% frog_yolk_set,],
                       by="Protein_ID")

XLA_mesh_conc <- apply(XLA_mesh_conc, 1, function(row){
  mesh_calculator(row, 13:20, 21:31, Frog_light_time, Frog_heavy_time) }) %>% bind_rows

Fly_mesh_conc <- merge(Dmel_prot_concs[1:4], Fly_gast_fits, by="Protein_ID")


Fly_mesh_conc <- apply(Fly_mesh_conc, 1, function(row){
  mesh_calculator(row, 12:16, 17:29, Fly_light_time, Fly_heavy_time)}) %>% bind_rows

head(XLA_mesh_conc)
```

    ##   Protein_ID T0_L           N1           N3           N4           N7
    ## 1 XBgroup104    0 0.000000e+00 7.707123e-06 0.000000e+00 0.000000e+00
    ## 2 XBgroup106    0 1.121476e-06 2.681220e-06 4.721364e-06 6.755482e-06
    ## 3 XBgroup109    0 3.945277e-06 9.251651e-06 3.659456e-05 3.175371e-05
    ## 4 XBgroup111    0 1.619357e-07 1.495262e-06 1.849185e-06 3.688551e-06
    ## 5 XBgroup112    0 1.337436e-07 2.687515e-08 7.083485e-09 1.805828e-07
    ## 6 XBgroup113    0 0.000000e+00 5.502687e-05 2.231594e-05 7.430808e-06
    ##             N9          N11          N12         T0_H           O1           O2
    ## 1 9.271818e-05 0.000000e+00 1.879298e-03 3.225629e-03 3.225629e-03 3.225629e-03
    ## 2 1.066707e-05 3.413793e-05 5.202872e-05 1.984376e-05 1.984376e-05 1.984376e-05
    ## 3 4.784126e-05 1.158917e-04 3.811711e-04 3.185393e-04 3.176232e-04 3.167096e-04
    ## 4 6.023563e-06 8.998699e-06 1.450575e-05 1.807641e-06 1.795039e-06 1.782525e-06
    ## 5 1.611136e-07 6.232077e-07 2.172957e-06 2.804193e-06 2.789708e-06 2.775297e-06
    ## 6 6.864126e-05 0.000000e+00 0.000000e+00 2.848291e-04 2.848291e-04 2.848291e-04
    ##             O3           O4           O5           O7           O9          O10
    ## 1 3.225629e-03 3.225629e-03 3.225629e-03 3.225629e-03 3.225629e-03 3.225629e-03
    ## 2 1.984376e-05 1.984376e-05 1.984376e-05 1.984376e-05 1.984376e-05 1.984376e-05
    ## 3 3.148905e-04 3.130818e-04 3.112461e-04 3.072750e-04 2.972282e-04 2.921709e-04
    ## 4 1.757759e-06 1.733336e-06 1.708755e-06 1.656282e-06 1.527805e-06 1.465430e-06
    ## 5 2.746700e-06 2.718396e-06 2.689804e-06 2.628410e-06 2.475903e-06 2.400669e-06
    ## 6 2.848291e-04 2.848291e-04 2.848291e-04 2.848291e-04 2.848291e-04 2.848291e-04
    ##            O11          O12
    ## 1 3.225629e-03 3.225629e-03
    ## 2 1.984376e-05 1.984376e-05
    ## 3 2.790322e-04 2.598274e-04
    ## 4 1.310477e-06 1.102066e-06
    ## 5 2.210047e-06 1.944086e-06
    ## 6 2.848291e-04 2.848291e-04

``` r
head(Fly_mesh_conc)
```

    ##   Protein_ID T0_L          N13          N14          N15          N16
    ## 1 A0A021WW64    0 0.000000e+00 4.786152e-07 1.157323e-06 9.184661e-07
    ## 2 A0A023GRW3    0 9.128538e-09 3.528000e-07 5.897489e-07 8.258937e-07
    ## 3 A0A0B4JCU3    0 3.868542e-07 3.234142e-06 2.903742e-06 5.165446e-06
    ## 4 A0A0B4JCZ0    0 0.000000e+00 2.048577e-06 9.683103e-06 1.738359e-05
    ## 5 A0A0B4JCZ1    0 3.248238e-07 4.635734e-07 4.353160e-07 3.928566e-07
    ## 6 A0A0B4JCZ3    0 4.930416e-05 9.271944e-05 8.902317e-05 7.317788e-05
    ##           T0_H           O1           O2           O3           O4           O5
    ## 1 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06
    ## 2 1.413181e-06 1.367754e-06 1.319470e-06 1.275665e-06 1.227954e-06 1.156559e-06
    ## 3 1.510468e-05 1.462445e-05 1.411382e-05 1.365038e-05 1.314541e-05 1.238937e-05
    ## 4 8.306108e-06 8.214621e-06 8.115148e-06 8.022802e-06 7.919802e-06 7.760610e-06
    ## 5 4.846549e-07 4.744919e-07 4.635585e-07 4.535173e-07 4.424417e-07 4.255811e-07
    ## 6 8.431519e-05 8.240434e-05 8.035238e-05 7.847135e-05 7.640047e-05 7.325611e-05
    ##             O6           O7           O8           O9          O10          O11
    ## 1 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06 4.533839e-06
    ## 2 1.084579e-06 1.010455e-06 8.895636e-07 7.264418e-07 6.326028e-07 4.961972e-07
    ## 3 1.162661e-05 1.084053e-05 9.557096e-06 7.822184e-06 6.822227e-06 5.365648e-06
    ## 4 7.593397e-06 7.413347e-06 7.099947e-06 6.628767e-06 6.325140e-06 5.825261e-06
    ## 5 4.082079e-07 3.898879e-07 3.589589e-07 3.147610e-07 2.877517e-07 2.458117e-07
    ## 6 7.002678e-05 6.663361e-05 6.093487e-05 5.286159e-05 4.797270e-05 4.045663e-05
    ##            O12
    ## 1 4.533839e-06
    ## 2 3.594589e-07
    ## 3 3.900982e-06
    ## 4 5.222232e-06
    ## 5 1.994304e-07
    ## 6 3.226707e-05

``` r
GB_NYS_time <- c(176, 286, 420, 807, 1488, 1920, 2930, 3641, 4406)

XLA_yolk_conc <- apply(merge(XLA_prot_concs[1:4], XLA_yolk_trend,
                              by="Protein_ID"), 1,
  function(row) {
    start_conc <- as.numeric(row["StartConc"])
    conc_deg <- as.numeric(row["Conc_Deg"])
    heavy_average <- as.numeric(row[5:13])
    heavy_average <- heavy_average/mean(heavy_average)

    # Maximum relative drop in raw vector + scaling factor
    max_rel_drop <- heavy_average[1] - tail(heavy_average, 1)
    H_scaling_factor <- conc_deg / max_rel_drop
    
    # For every point, calculate how far it has dropped from the start,
    # multiply by the scaling factor, and subtract it from your starting concentration.
    absolute_heavy_pool <- start_conc - ((heavy_average[1] - 
                                            heavy_average) * H_scaling_factor)
    names(absolute_heavy_pool) <- names(row[5:13])
    absolute_heavy_pool <- cummin(absolute_heavy_pool) #fixes to a strict decline

    outline <- data.frame(t(absolute_heavy_pool))
    outline <- cbind(data.frame(Protein_ID = as.character(row["Protein_ID"])),
                     outline)

    return(outline) }) %>% bind_rows

head(XLA_yolk_conc)
```

    ##                   Protein_ID         T0         T1         T2         T3
    ## 1 XBmRNA35012|XBXL10_1g18988 0.23554408 0.23367572 0.23367572 0.23350936
    ## 2 XBmRNA35013|XBXL10_1g18989 0.33282334 0.32708783 0.32674908 0.32607923
    ## 3 XBmRNA35014|XBXL10_1g18990 0.47148403 0.47148403 0.47148403 0.47148403
    ## 4 XBmRNA39801|XBXL10_1g21498 0.07422721 0.07379854 0.07144474 0.07144474
    ## 5 XBmRNA67639|XBXL10_1g35803 1.10050710 1.09595049 1.05362335 1.03133796
    ##           T4         T5         T6         T7         T8
    ## 1 0.23017110 0.22736725 0.22321686 0.22067047 0.21486741
    ## 2 0.31933049 0.31877662 0.31246044 0.30590618 0.30157898
    ## 3 0.45937848 0.45077272 0.44638691 0.44168769 0.44168769
    ## 4 0.07138531 0.06999798 0.06999798 0.06724045 0.06692918
    ## 5 1.03133796 1.03133796 0.94997866 0.87469532 0.87469532

``` r
XLA_mesh_conc_avg <- apply(XLA_mesh_conc[!XLA_mesh_conc$Protein_ID %in%
                                    frog_yolk_set, 2:20], 2, sum) %>%
  as.data.frame()
colnames(XLA_mesh_conc_avg) <- "Trend"
XLA_mesh_conc_avg["Time"] <- c(Frog_light_hpf_time, Frog_heavy_hpf_time)
XLA_mesh_conc_avg["Type"] <- c(rep("New", length(Frog_light_hpf_time)),
                               rep("Start", length(Frog_heavy_hpf_time)))

XLA_yolk_conc_avg <- apply(XLA_yolk_conc[2:10], 2, sum) %>%
  as.data.frame

colnames(XLA_yolk_conc_avg) <- "Trend"
XLA_yolk_conc_avg["Time"] <- GB_NYS_time/60
XLA_yolk_conc_avg["Type"] <- "Yolk"

XLA_mesh_conc_avg <- rbind(XLA_mesh_conc_avg, XLA_yolk_conc_avg)

head(XLA_mesh_conc_avg)
```

    ##           Trend      Time Type
    ## T0_L 0.00000000  2.591667  New
    ## N1   0.01870066  3.591667  New
    ## N3   0.01753008  6.591667  New
    ## N4   0.02313820  8.591667  New
    ## N7   0.02733216 15.091667  New
    ## N9   0.03581946 26.633333  New

``` r
Fly_yolk_conc_avg<- apply(Fly_mesh_conc[Fly_mesh_conc$Protein_ID %in%
                                    fly_yolk_set,7:19], 2, sum) %>%
  as.data.frame()
colnames(Fly_yolk_conc_avg) <- "Trend"
Fly_yolk_conc_avg["Time"] <- Fly_heavy_hpf_time
Fly_yolk_conc_avg["Type"] <- "Yolk"

Fly_mesh_conc_avg <- apply(Fly_mesh_conc[!Fly_mesh_conc$Protein_ID %in%
                                    fly_yolk_set,2:19], 2, sum) %>%
  as.data.frame()
colnames(Fly_mesh_conc_avg) <- "Trend"
Fly_mesh_conc_avg["Time"] <- c(Fly_light_hpf_time, Fly_heavy_hpf_time)
Fly_mesh_conc_avg["Type"] <- c(rep("New", length(Fly_light_hpf_time)),
                               rep("Start", length(Fly_heavy_hpf_time)))
Fly_mesh_conc_avg <- rbind(Fly_mesh_conc_avg, Fly_yolk_conc_avg)
```

``` r
mesh_plot <- function(df, xlab="Hours Post-Fertilization", ylab="Protein",
                      fun_breaks, approx_x_min, approx_x_max) {
  
  time_grid <- seq(approx_x_min, approx_x_max, length.out = 100)
  
  new_interp <- approx(df[df$Type == "New",]$Time,
                       df[df$Type == "New",]$Trend,
                       xout = time_grid)$y

  start_interp <- approx(df[df$Type == "Start",]$Time,
                         df[df$Type == "Start",]$Trend,
                         xout = time_grid)$y
  
  yolk_interp <- approx(df[df$Type == "Yolk",]$Time,
                        df[df$Type == "Yolk",]$Trend,
                        xout = time_grid)$y

  # Stacked total must not exceed the initial (t0) total. Start and Yolk are
  # fitted/monotone curves; New is the observed, noise-prone term, so any
  # overshoot is attributed to (and trimmed from) New.
  total0     <- new_interp[1] + start_interp[1] + yolk_interp[1]
  headroom   <- pmax(total0 - start_interp - yolk_interp, 0)
  new_interp <- pmin(new_interp, headroom)  
  
  interp_df <- data.frame(Time=rep(time_grid, 3),
                          Trend=c(new_interp, start_interp, yolk_interp),
                          Type=c(rep("New", length(time_grid)),
                                 rep("Start", length(time_grid)),
                                 rep("Yolk", length(time_grid))))
  # interp_df["Type"] <- factor(interp_df$Type, levels=c("New", "Start", "Yolk"))
  interp_df["Type"] <- factor(interp_df$Type, levels=c("Yolk", "New", "Start"))
  
  # interp_df <- data.frame(Time=rep(time_grid, 2),
  #                         Trend=c(new_interp, start_interp),
  #                         Type=c(rep("New", length(time_grid)),
  #                                rep("Start", length(time_grid))))
  # interp_df["Type"] <- factor(interp_df$Type, levels=c("New", "Start"))
  
  p1 <- ggplot(interp_df, aes(x = Time, y = Trend, fill = Type)) +
    geom_area(position = "stack", color="black", linewidth=0.5) +
    # scale_fill_manual(values = c("New" = "#E69F00", "Start" = "#999999")) +
    scale_fill_manual(values = c("New" = "#E69F00", "Start" = "#999999",
                                 "Yolk"="#56B4E9")) +
    labs(x = xlab, y = ylab) +
    theme_bw() +
    theme(axis.text= element_text(size=21,colour="black"),
          # axis.text.y = element_blank(),
          # axis.ticks = element_blank(),
          plot.title = element_text(size=24, hjust=0.5),
          axis.title = element_text(size=22),
          panel.border = element_rect(linewidth=2),
          legend.position = "none",
          panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
          panel.grid.minor = element_line(color = "grey90", linewidth = 0.25)) +
    scale_x_continuous(breaks=fun_breaks)

  return(p1) }

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=5, height=4.5, res=300)

mesh_plot(XLA_mesh_conc_avg, xlab="Hours Post-Fertilization (16°C)",
          ylab="Protein Concentration", fun_breaks=c(3, 24, 48, 72),
          approx_x_min=3, approx_x_max=72)
```

![](figures/AA-Tracking_Fly_and_Frog/Mesh%20plots%20for%20frog%20and%20fly%20as%20concentration%20over%20time-1.png)<!-- -->

``` r
mesh_plot(Fly_mesh_conc_avg, xlab="Hours Post-Fertilization (22°C)",
          ylab="Protein Concentration", fun_breaks=c(3, 10, 24),
          approx_x_min=3.2, approx_x_max=24)
```

![](figures/AA-Tracking_Fly_and_Frog/Mesh%20plots%20for%20frog%20and%20fly%20as%20concentration%20over%20time-2.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

# Mesh plot as protein amount

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

frog_sizes <- frog_sizes %>% dplyr::select(Protein_ID, kDa)
fly_sizes <- fly_sizes %>% dplyr::select(Protein_ID, kDa)
```

``` r
mesh_ug_conv <- function(row, data_cols, volume) {
  
  size <- as.numeric(row["kDa"])
  conc_vec <- row[data_cols]

  ug_vec <- as.numeric(conc_vec) * volume * size
  names(ug_vec) <- names(conc_vec)
  
  outline <- data.frame(t(ug_vec))
  outline <- cbind(data.frame(Protein_ID = as.character(row["Protein_ID"])),
                   outline)

  return(outline) }

frog_mesh_ug <- apply(merge(frog_sizes, XLA_mesh_conc[!XLA_mesh_conc$Protein_ID %in% frog_yolk_set,],
                            by="Protein_ID"),
                      1, function(x) { mesh_ug_conv(x, 3:21, 1) }) %>% bind_rows

fly_mesh_ug <- apply(merge(fly_sizes, Fly_mesh_conc, by="Protein_ID"),
      1, function(x) { mesh_ug_conv(x, 3:20, 0.01) }) %>% bind_rows

frog_mesh_yolk_ug <- apply(merge(frog_sizes, XLA_yolk_conc, by="Protein_ID"), 1,
                           function(x) { mesh_ug_conv(x, 3:11, 1) }) %>% bind_rows
```

``` r
XLA_mesh_ug_avg <- apply(frog_mesh_ug[!frog_mesh_ug$Protein_ID %in%
                                        frog_yolk_set, 2:20], 2, sum) %>%
  as.data.frame()
colnames(XLA_mesh_ug_avg) <- "Trend"
XLA_mesh_ug_avg["Time"] <- c(Frog_light_hpf_time, Frog_heavy_hpf_time)
XLA_mesh_ug_avg["Type"] <- c(rep("New", length(Frog_light_hpf_time)),
                             rep("Start", length(Frog_heavy_hpf_time)))

XLA_yolk_ug_avg <- apply(frog_mesh_yolk_ug[2:10], 2, sum) %>%
  as.data.frame

colnames(XLA_yolk_ug_avg) <- "Trend"
XLA_yolk_ug_avg["Time"] <- GB_NYS_time/60
XLA_yolk_ug_avg["Type"] <- "Yolk"

XLA_mesh_ug_avg <- rbind(XLA_mesh_ug_avg, XLA_yolk_ug_avg)
```

``` r
Fly_yolk_ug_avg<- apply(fly_mesh_ug[fly_mesh_ug$Protein_ID %in%
                                    fly_yolk_set,7:19], 2, sum) %>%
  as.data.frame()
colnames(Fly_yolk_ug_avg) <- "Trend"
Fly_yolk_ug_avg["Time"] <- Fly_heavy_hpf_time
Fly_yolk_ug_avg["Type"] <- "Yolk"

Fly_mesh_ug_avg <- apply(fly_mesh_ug[!fly_mesh_ug$Protein_ID %in%
                                    fly_yolk_set,2:19], 2, sum) %>%
  as.data.frame()
colnames(Fly_mesh_ug_avg) <- "Trend"
Fly_mesh_ug_avg["Time"] <- c(Fly_light_hpf_time, Fly_heavy_hpf_time)
Fly_mesh_ug_avg["Type"] <- c(rep("New", length(Fly_light_hpf_time)),
                             rep("Start", length(Fly_heavy_hpf_time)))
Fly_mesh_ug_avg <- rbind(Fly_mesh_ug_avg, Fly_yolk_ug_avg)
```

``` r
# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=10, height=4.5, res=300)

mesh_plot(XLA_mesh_ug_avg, xlab="Hours post-fertilization (16°C)",
          ylab="Protein amount (μg)", fun_breaks=c(3, 24, 48, 72),
          approx_x_min=3, approx_x_max=72) |
mesh_plot(Fly_mesh_ug_avg, xlab="Hours post-fertilization (22°C)",
          ylab="Protein amount (μg)", fun_breaks=c(3, 10, 24),
          approx_x_min=3.2, approx_x_max=24) +
  theme(axis.title.y = element_blank())
```

![](figures/AA-Tracking_Fly_and_Frog/Mesh%20plots%20for%20frog%20and%20fly%20as%20amount%20over%20time-1.png)<!-- -->

``` r
# dev.off() #Uncomment for new image!
```

``` r
frog_base_mesh <- mesh_plot(XLA_mesh_ug_avg, xlab="Hours post-fertilization (16°C)",
                            ylab="Protein amount (μg)", fun_breaks=c(3, 24, 48, 72),
                            approx_x_min=3, approx_x_max=72)

frog_top_mesh <- frog_base_mesh +
  coord_cartesian(ylim = c(290, 325)) +
  scale_y_continuous(breaks=c(240, 280, 320)) +  
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        plot.margin = margin(t = 10, r = 5, b = 0,
                             l = 10, unit = "pt"))

frog_bottom_mesh <- frog_base_mesh +
  coord_cartesian(ylim = c(-2.45, 60)) +
  scale_y_continuous(breaks=c(0, 20, 40, 60)) +
  theme(plot.margin = margin(t = 0, r = 5, b = 10,
                             l = 10, unit = "pt"))

frog_split_compiled <- (frog_top_mesh / frog_bottom_mesh) + 
  plot_layout(heights = c(20, 80), axis_titles = "collect")


fly_base_mesh <- mesh_plot(Fly_mesh_ug_avg, xlab="Hours post-fertilization (22°C)",
          ylab="Protein amount (μg)", fun_breaks=c(3, 10, 24),
          approx_x_min=3.2, approx_x_max=24) +
  theme(plot.margin = margin(t = 10, r = 10, b = 10,
                             l = 5, unit = "pt"),
        axis.title.y = element_blank())

# #Uncomment for new image!
# tiff("graph.tiff", units="in",
#      width=9.5, height=4.5, res=300)

frog_split_compiled | fly_base_mesh
```

![](figures/AA-Tracking_Fly_and_Frog/Splitting%20frog%20data%20to%20show%20more%20nonyolk%20part-1.png)<!-- -->

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
    ##  [1] patchwork_1.3.2     ggbreak_0.1.7       stringr_1.6.0      
    ##  [4] Peptides_2.4.6      Biostrings_2.78.0   Seqinfo_1.0.0      
    ##  [7] XVector_0.50.0      IRanges_2.44.0      S4Vectors_0.48.1   
    ## [10] BiocGenerics_0.56.0 generics_0.1.4      ggplot2_4.0.3      
    ## [13] tidyr_1.3.2         dplyr_1.2.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] yulab.utils_0.2.4  rappdirs_0.3.4     ggplotify_0.1.3    stringi_1.8.7     
    ##  [5] digest_0.6.39      magrittr_2.0.4     evaluate_1.0.5     grid_4.5.3        
    ##  [9] RColorBrewer_1.1-3 fastmap_1.2.0      ggnewscale_0.5.2   purrr_1.2.2       
    ## [13] aplot_0.2.9        scales_1.4.0       cli_3.6.5          rlang_1.1.7       
    ## [17] crayon_1.5.3       withr_3.0.3        yaml_2.3.12        otel_0.2.0        
    ## [21] tools_4.5.3        vctrs_0.7.1        R6_2.6.1           gridGraphics_0.5-1
    ## [25] lifecycle_1.0.5    fs_2.1.0           ggfun_0.2.0        pkgconfig_2.0.3   
    ## [29] pillar_1.11.1      gtable_0.3.6       glue_1.8.0         Rcpp_1.1.1-1.1    
    ## [33] xfun_0.58          tibble_3.3.1       tidyselect_1.2.1   rstudioapi_0.19.0 
    ## [37] knitr_1.51         farver_2.1.2       htmltools_0.5.9    labeling_0.4.3    
    ## [41] rmarkdown_2.31     compiler_4.5.3     S7_0.2.1
