DiffBind_analysis_onlyPeakCount
================
Guandong Shang
2022-05-29

有些时候，你想要其他的软件，比如说DESeq2或者edgeR等去分析ChIP-seq和ATAC-seq的数据。他们会要求一个表达量矩阵，这里我演示如何得到feature
count矩阵，尽管这些内容在前面的两个文档中都已经提过了。

这里用的是CIM0d-SIM8d的ATAC-seq数据

``` bash
tree rawdata/{bam,peak}
## rawdata/bam
## ├── C01.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C01.rm_organelle.bam
## ├── C01.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C01.rm_organelle.bam.bai
## ├── C02.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C02.rm_organelle.bam
## ├── C02.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C02.rm_organelle.bam.bai
## ├── C31.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C31.rm_organelle.bam
## ├── C31.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C31.rm_organelle.bam.bai
## ├── C32.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C32.rm_organelle.bam
## ├── C32.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C32.rm_organelle.bam.bai
## ├── C51.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C51.rm_organelle.bam
## ├── C51.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C51.rm_organelle.bam.bai
## ├── C52.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C52.rm_organelle.bam
## ├── C52.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/C52.rm_organelle.bam.bai
## ├── s01.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s01.rm_organelle.bam
## ├── s01.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s01.rm_organelle.bam.bai
## ├── s02.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s02.rm_organelle.bam
## ├── s02.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s02.rm_organelle.bam.bai
## ├── s31.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s31.rm_organelle.bam
## ├── s31.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s31.rm_organelle.bam.bai
## ├── s32.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s32.rm_organelle.bam
## ├── s32.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s32.rm_organelle.bam.bai
## ├── s61.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s61.rm_organelle.bam
## ├── s61.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s61.rm_organelle.bam.bai
## ├── s62.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s62.rm_organelle.bam
## ├── s62.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s62.rm_organelle.bam.bai
## ├── s81.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s81.rm_organelle.bam
## ├── s81.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s81.rm_organelle.bam.bai
## ├── s82.rm_organelle.bam -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s82.rm_organelle.bam
## └── s82.rm_organelle.bam.bai -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/bam/s82.rm_organelle.bam.bai
## rawdata/peak
## ├── C01_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C01_peaks.narrowPeak
## ├── C02_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C02_peaks.narrowPeak
## ├── C31_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C31_peaks.narrowPeak
## ├── C32_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C32_peaks.narrowPeak
## ├── C51_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C51_peaks.narrowPeak
## ├── C52_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/C52_peaks.narrowPeak
## ├── s01_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s01_peaks.narrowPeak
## ├── s02_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s02_peaks.narrowPeak
## ├── s31_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s31_peaks.narrowPeak
## ├── s32_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s32_peaks.narrowPeak
## ├── s61_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s61_peaks.narrowPeak
## ├── s62_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s62_peaks.narrowPeak
## ├── s81_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s81_peaks.narrowPeak
## └── s82_peaks.narrowPeak -> /home/sgd/project/202005/WLY_Total_Fig_202005/ATAC_CIM_SIM/rawdata/peak/s82_peaks.narrowPeak
## 
## 0 directories, 42 files
```

# 前期准备

``` r
library(DiffBind)
## Loading required package: GenomicRanges
## Loading required package: stats4
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
##     union, unique, unsplit, which.max, which.min
## Loading required package: S4Vectors
## 
## Attaching package: 'S4Vectors'
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## Loading required package: IRanges
## Loading required package: GenomeInfoDb
## Loading required package: SummarizedExperiment
## Loading required package: MatrixGenerics
## Loading required package: matrixStats
## 
## Attaching package: 'MatrixGenerics'
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## Attaching package: 'Biobase'
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
##  >>> DiffBind 3.6.1

library(BiocParallel)

library(magrittr)
## 
## Attaching package: 'magrittr'
## The following object is masked from 'package:GenomicRanges':
## 
##     subtract
```

``` r
bam_files <- list.files("rawdata/bam",pattern = "bam$", full.names = T)
peak_files <- list.files("rawdata/peak", full.names = T)

bam_files
##  [1] "rawdata/bam/C01.rm_organelle.bam" "rawdata/bam/C02.rm_organelle.bam"
##  [3] "rawdata/bam/C31.rm_organelle.bam" "rawdata/bam/C32.rm_organelle.bam"
##  [5] "rawdata/bam/C51.rm_organelle.bam" "rawdata/bam/C52.rm_organelle.bam"
##  [7] "rawdata/bam/s01.rm_organelle.bam" "rawdata/bam/s02.rm_organelle.bam"
##  [9] "rawdata/bam/s31.rm_organelle.bam" "rawdata/bam/s32.rm_organelle.bam"
## [11] "rawdata/bam/s61.rm_organelle.bam" "rawdata/bam/s62.rm_organelle.bam"
## [13] "rawdata/bam/s81.rm_organelle.bam" "rawdata/bam/s82.rm_organelle.bam"

peak_files
##  [1] "rawdata/peak/C01_peaks.narrowPeak" "rawdata/peak/C02_peaks.narrowPeak"
##  [3] "rawdata/peak/C31_peaks.narrowPeak" "rawdata/peak/C32_peaks.narrowPeak"
##  [5] "rawdata/peak/C51_peaks.narrowPeak" "rawdata/peak/C52_peaks.narrowPeak"
##  [7] "rawdata/peak/s01_peaks.narrowPeak" "rawdata/peak/s02_peaks.narrowPeak"
##  [9] "rawdata/peak/s31_peaks.narrowPeak" "rawdata/peak/s32_peaks.narrowPeak"
## [11] "rawdata/peak/s61_peaks.narrowPeak" "rawdata/peak/s62_peaks.narrowPeak"
## [13] "rawdata/peak/s81_peaks.narrowPeak" "rawdata/peak/s82_peaks.narrowPeak"
```

``` r
tissue <- rep(c("CIM0d", "CIM3d", "CIM5d", "CIM7d", "SIM3d", "SIM6d", "SIM8d"), each = 2)

# 注意这里ATAC-seq没有control相关的几列
sample_info <- data.frame(
  SampleID = paste(tissue, c("R1", "R2"), sep = "_") ,
  
  Tissue = tissue,
  
  Replicate = 1:2,
  
  bamReads = bam_files, 
  
  Peaks = peak_files, 
  
  PeakCaller = "narrow"
  
) 

sample_info
##    SampleID Tissue Replicate                         bamReads
## 1  CIM0d_R1  CIM0d         1 rawdata/bam/C01.rm_organelle.bam
## 2  CIM0d_R2  CIM0d         2 rawdata/bam/C02.rm_organelle.bam
## 3  CIM3d_R1  CIM3d         1 rawdata/bam/C31.rm_organelle.bam
## 4  CIM3d_R2  CIM3d         2 rawdata/bam/C32.rm_organelle.bam
## 5  CIM5d_R1  CIM5d         1 rawdata/bam/C51.rm_organelle.bam
## 6  CIM5d_R2  CIM5d         2 rawdata/bam/C52.rm_organelle.bam
## 7  CIM7d_R1  CIM7d         1 rawdata/bam/s01.rm_organelle.bam
## 8  CIM7d_R2  CIM7d         2 rawdata/bam/s02.rm_organelle.bam
## 9  SIM3d_R1  SIM3d         1 rawdata/bam/s31.rm_organelle.bam
## 10 SIM3d_R2  SIM3d         2 rawdata/bam/s32.rm_organelle.bam
## 11 SIM6d_R1  SIM6d         1 rawdata/bam/s61.rm_organelle.bam
## 12 SIM6d_R2  SIM6d         2 rawdata/bam/s62.rm_organelle.bam
## 13 SIM8d_R1  SIM8d         1 rawdata/bam/s81.rm_organelle.bam
## 14 SIM8d_R2  SIM8d         2 rawdata/bam/s82.rm_organelle.bam
##                                Peaks PeakCaller
## 1  rawdata/peak/C01_peaks.narrowPeak     narrow
## 2  rawdata/peak/C02_peaks.narrowPeak     narrow
## 3  rawdata/peak/C31_peaks.narrowPeak     narrow
## 4  rawdata/peak/C32_peaks.narrowPeak     narrow
## 5  rawdata/peak/C51_peaks.narrowPeak     narrow
## 6  rawdata/peak/C52_peaks.narrowPeak     narrow
## 7  rawdata/peak/s01_peaks.narrowPeak     narrow
## 8  rawdata/peak/s02_peaks.narrowPeak     narrow
## 9  rawdata/peak/s31_peaks.narrowPeak     narrow
## 10 rawdata/peak/s32_peaks.narrowPeak     narrow
## 11 rawdata/peak/s61_peaks.narrowPeak     narrow
## 12 rawdata/peak/s62_peaks.narrowPeak     narrow
## 13 rawdata/peak/s81_peaks.narrowPeak     narrow
## 14 rawdata/peak/s82_peaks.narrowPeak     narrow
```

``` r
# Reading in peaksets
dba <- dba(sampleSheet = sample_info)
## CIM0d_R1 CIM0d    1 narrow
## CIM0d_R2 CIM0d    2 narrow
## CIM3d_R1 CIM3d    1 narrow
## CIM3d_R2 CIM3d    2 narrow
## CIM5d_R1 CIM5d    1 narrow
## CIM5d_R2 CIM5d    2 narrow
## CIM7d_R1 CIM7d    1 narrow
## CIM7d_R2 CIM7d    2 narrow
## SIM3d_R1 SIM3d    1 narrow
## SIM3d_R2 SIM3d    2 narrow
## SIM6d_R1 SIM6d    1 narrow
## SIM6d_R2 SIM6d    2 narrow
## SIM8d_R1 SIM8d    1 narrow
## SIM8d_R2 SIM8d    2 narrow
dba
## 14 Samples, 27602 sites in matrix (29941 total):
##          ID Tissue Replicate Intervals
## 1  CIM0d_R1  CIM0d         1     20971
## 2  CIM0d_R2  CIM0d         2     20350
## 3  CIM3d_R1  CIM3d         1     21916
## 4  CIM3d_R2  CIM3d         2     18774
## 5  CIM5d_R1  CIM5d         1     24260
## 6  CIM5d_R2  CIM5d         2     22883
## 7  CIM7d_R1  CIM7d         1     22482
## 8  CIM7d_R2  CIM7d         2     24312
## 9  SIM3d_R1  SIM3d         1     21663
## 10 SIM3d_R2  SIM3d         2     22669
## 11 SIM6d_R1  SIM6d         1     21591
## 12 SIM6d_R2  SIM6d         2     21227
## 13 SIM8d_R1  SIM8d         1     26310
## 14 SIM8d_R2  SIM8d         2     25667
```

``` r
# Counting reads
# ATAC-seq没有input，所以不需要做graylist
# 由于ATAC-seq相对来说也是narrowPeak，所以我这里将summit关闭

# 同时注意，因为我们这里只想要raw count做后面的分析
# 而不是矫正过后的count
# 所以这里的score要设置为DBA_SCORE_READS
dba <- dba.count(dba, summits = FALSE, score = DBA_SCORE_READS)

dba
## 14 Samples, 27602 sites in matrix:
##          ID Tissue Replicate    Reads FRiP
## 1  CIM0d_R1  CIM0d         1  6970957 0.75
## 2  CIM0d_R2  CIM0d         2  5142481 0.74
## 3  CIM3d_R1  CIM3d         1  4952533 0.75
## 4  CIM3d_R2  CIM3d         2  3730200 0.64
## 5  CIM5d_R1  CIM5d         1  8334734 0.73
## 6  CIM5d_R2  CIM5d         2  9016566 0.70
## 7  CIM7d_R1  CIM7d         1 28397538 0.75
## 8  CIM7d_R2  CIM7d         2 25294498 0.75
## 9  SIM3d_R1  SIM3d         1 22335798 0.68
## 10 SIM3d_R2  SIM3d         2 26180816 0.60
## 11 SIM6d_R1  SIM6d         1  7758298 0.73
## 12 SIM6d_R2  SIM6d         2  6721092 0.74
## 13 SIM8d_R1  SIM8d         1 21225290 0.65
## 14 SIM8d_R2  SIM8d         2 15985456 0.65
```

# 提取和输出结果

``` r
peakInfo <- dba.peakset(dba, bRetrieve = TRUE)
names(peakInfo) <- paste("ATAC_", names(peakInfo))

# 注意到这里的值都是整数，而不像我们之前的文档里面有小数点
# 因为我们提取的是rawcount
peakInfo
## GRanges object with 27602 ranges and 14 metadata columns:
##               seqnames            ranges strand |  CIM0d_R1  CIM0d_R2  CIM3d_R1
##                  <Rle>         <IRanges>  <Rle> | <numeric> <numeric> <numeric>
##       ATAC_ 1     Chr1         1422-3756      * |       212       153       266
##       ATAC_ 2     Chr1         6494-6860      * |        31         7        44
##       ATAC_ 3     Chr1         8527-8978      * |        72        53        61
##       ATAC_ 4     Chr1         9381-9877      * |        44        28        70
##       ATAC_ 5     Chr1       13375-15027      * |        17        14        22
##           ...      ...               ...    ... .       ...       ...       ...
##   ATAC_ 27598     Chr5 26966551-26967607      * |       198       114       158
##   ATAC_ 27599     Chr5 26969183-26969925      * |        76        48        77
##   ATAC_ 27600     Chr5 26970635-26971311      * |       110        60        76
##   ATAC_ 27601     Chr5 26972192-26972813      * |       134        73        61
##   ATAC_ 27602     Chr5 26975252-26975432      * |         9        10        15
##                CIM3d_R2  CIM5d_R1  CIM5d_R2  CIM7d_R1  CIM7d_R2  SIM3d_R1
##               <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
##       ATAC_ 1       153       318       364      1470      1284      1047
##       ATAC_ 2        31        55        81       167       134       125
##       ATAC_ 3        51       101       121       300       336       215
##       ATAC_ 4        45       102        86       353       303       289
##       ATAC_ 5        30        68        85       421       407       463
##           ...       ...       ...       ...       ...       ...       ...
##   ATAC_ 27598       105       179       211       639       625       456
##   ATAC_ 27599        45       101       102       355       312       247
##   ATAC_ 27600        53       112       129       440       381       220
##   ATAC_ 27601        50       118       106       357       369       249
##   ATAC_ 27602        11        36        44       105       108        78
##                SIM3d_R2  SIM6d_R1  SIM6d_R2  SIM8d_R1  SIM8d_R2
##               <numeric> <numeric> <numeric> <numeric> <numeric>
##       ATAC_ 1       958       324       257       862       657
##       ATAC_ 2       130        56        56       147        94
##       ATAC_ 3       208       108        93       290       186
##       ATAC_ 4       223        73        55       154       114
##       ATAC_ 5       521       268       247       454       323
##           ...       ...       ...       ...       ...       ...
##   ATAC_ 27598       481       215       195       549       395
##   ATAC_ 27599       233       132       109       287       207
##   ATAC_ 27600       260       122       105       259       183
##   ATAC_ 27601       304       143       102       312       260
##   ATAC_ 27602       100        34        21        55        40
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths

# 这里的矩阵就可以作为一些软件的输入了
peakCount_mt <- as.matrix(mcols(peakInfo))
head(peakCount_mt)
##         CIM0d_R1 CIM0d_R2 CIM3d_R1 CIM3d_R2 CIM5d_R1 CIM5d_R2 CIM7d_R1 CIM7d_R2
## ATAC_ 1      212      153      266      153      318      364     1470     1284
## ATAC_ 2       31        7       44       31       55       81      167      134
## ATAC_ 3       72       53       61       51      101      121      300      336
## ATAC_ 4       44       28       70       45      102       86      353      303
## ATAC_ 5       17       14       22       30       68       85      421      407
## ATAC_ 6        3        3        9        8       23       24      145      129
##         SIM3d_R1 SIM3d_R2 SIM6d_R1 SIM6d_R2 SIM8d_R1 SIM8d_R2
## ATAC_ 1     1047      958      324      257      862      657
## ATAC_ 2      125      130       56       56      147       94
## ATAC_ 3      215      208      108       93      290      186
## ATAC_ 4      289      223       73       55      154      114
## ATAC_ 5      463      521      268      247      454      323
## ATAC_ 6      106       76       23       14       63       46
```

``` r
mcols(peakInfo) <- NULL
peakInfo$feature_id <- names(peakInfo)

# 这里的peakInfo可以用作注释软件的输入
peakInfo
## GRanges object with 27602 ranges and 1 metadata column:
##               seqnames            ranges strand |  feature_id
##                  <Rle>         <IRanges>  <Rle> | <character>
##       ATAC_ 1     Chr1         1422-3756      * |     ATAC_ 1
##       ATAC_ 2     Chr1         6494-6860      * |     ATAC_ 2
##       ATAC_ 3     Chr1         8527-8978      * |     ATAC_ 3
##       ATAC_ 4     Chr1         9381-9877      * |     ATAC_ 4
##       ATAC_ 5     Chr1       13375-15027      * |     ATAC_ 5
##           ...      ...               ...    ... .         ...
##   ATAC_ 27598     Chr5 26966551-26967607      * | ATAC_ 27598
##   ATAC_ 27599     Chr5 26969183-26969925      * | ATAC_ 27599
##   ATAC_ 27600     Chr5 26970635-26971311      * | ATAC_ 27600
##   ATAC_ 27601     Chr5 26972192-26972813      * | ATAC_ 27601
##   ATAC_ 27602     Chr5 26975252-26975432      * | ATAC_ 27602
##   -------
##   seqinfo: 5 sequences from an unspecified genome; no seqlengths
```
