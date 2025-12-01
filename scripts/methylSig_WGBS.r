library(methylSig)
argv<-commandArgs(T)
sampleInfoFilePath<-argv[1]
sampleSheet<-read.table(sampleInfoFilePath, header = TRUE, colClasses = c("character", "character"))
# name   condition
# SRR13171471    control
# SRR13171472    control
# SRR13171473    test
# SRR13171474    test
files<-paste0("mapped/",sampleSheet$name,".WGBS_filter.bitmapperbs_CpG.bedGraph")
bsseq_stranded = bsseq::read.bismark(
    files = files,
    colData = data.frame(row.names = sampleSheet$name),
    rmZeroCov = FALSE,
    strandCollapse = FALSE
)
# bsseq_stranded$Type<-c("case","case","control","control")
bsseq_stranded$Type<-sampleSheet$condition
bs = filter_loci_by_coverage(bsseq_stranded,min_count = 5,max_count = 500)
windowed_bs = tile_by_windows(bs = bs,win_size = 10000)
diff_fit_simple = diff_dss_fit(
    bs = bs,
    design = bsseq::pData(bs),
    formula = as.formula('~ Type'))
simple_contrast = matrix(c(0,1),ncol = 1)
diff_simple_gr = diff_dss_test(
    bs = bs,
    diff_fit = diff_fit_simple,
    contrast = simple_contrast,
    methylation_group_column = 'Type',
    methylation_groups = c('case' = 'case', 'control' = 'control'))
write.table(diff_simple_gr, file = paste0("text/diff_simple_gr.txt"),
        row.names = FALSE, quote = FALSE, sep = "\t")