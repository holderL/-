附页R代码
# 载入需要的R包
library(ShortRead)
library(GenomicRanges)
library(Rsamtools)
library(DESeq2)

# 读取测序数据
fastqFiles <- c("sample1.fastq.gz", "sample2.fastq.gz")
reads <- readFastq(fastqFiles)

# 数据预处理：质量控制、去除接头序列
filteredReads <- qualityFilter(reads, cutoff=20, maxN=0, maxConsecN=0, maxLen=100)
filteredReads <- trimLRPatterns(filteredReads, c("TTAGGG"))

# 参考基因组比对
indexFile <- "hg19.fa"
bamFile <- "aligned.bam"
align(index=indexFile, readfile=filteredReads, output=bamFile)

# 标记PCR重复序列，去除PCR重复
sortedFile <- "sorted.bam"
sortBam(bamFile, sortedFile)
dedupFile <- "deduplicated.bam"
deduplicateBam(sortedFile, dedupFile, tag="RG", remove=TRUE)

# 差异表达分析
counts <- featureCounts(dedupFile, annot.ext="hg19.gtf", isPairedEnd=FALSE)
colData <- data.frame(condition=c("control", "treatment"))
dds <- DESeqDataSetFromMatrix(countData=counts, colData=colData, design=~condition)
dds <- DESeq(dds)
results <- results(dds, contrast=c("condition", "treatment", "control"))
# 安装和加载所需的包
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("methylKit")

library(methylKit)

# 读取测序数据
setwd("/path/to/your/data")
files <- list.files(pattern=".txt.gz$")
sampleNames <- sapply(strsplit(files, "[._]"), `[`, 1)

# 创建methylRawList对象
myobj <- read.metharray.exp(files, sampleNames = sampleNames)

# 过滤掉控制位、质控位和低甲基化位点
myobj.filtered <- filterData(myobj, type = "Illumina")

# 将甲基化水平转换为百分比
myobj.norm <- preprocessRaw(myobj.filtered)
myobj.final <- normalizeData(myobj.norm)

# 将样本按组分配
group <- c(rep("A", 3), rep("B", 3))
myobj.final$Sample_Group <- factor(group)

# 比较组之间的甲基化水平差异
meth.diff <- calculateDiffMeth(myobj.final, group1 = "A", group2 = "B", difference = 25, adjust.method = "fdr")

# 提取差异甲基化位点
dmrs <- getMethylDiff(meth.diff, difference = 25, qvalue = 0.01)

# 导出结果
write.table(dmrs, file = "dmrs.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
# 加载所需的包
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(GenomicAlignments)

# 获取参考基因组序列
hg19 <- BSgenome.Hsapiens.UCSC.hg19

# 建立比对索引
indexFa(hg19)

# 读取比对索引
hg19_idx <- readIndex(hg19)

# 将序列比对到基因组上
alignments <- readGAlignments("aligned_reads.bam", use.names = TRUE)
aligned <- suppressWarnings(coverage(hg19_idx, alignments))

# 获取甲基化水平信息
meth <- methylRaw("aligned_reads.bam", genome = hg19, sample.id = c("M0_1", "M0_2", "M0_3", "M1_1", "M1_2", "M1_3", "M2_1", "M2_2", "M2_3"))
meth_bismark <- bismark(bam = meth, genome = hg19)

# 对甲基化水平信息进行过滤和归一化
meth_filt <- filterByCoverage(meth_bismark)
meth_norm <- normalizeCoverage(meth_filt)

# 比较甲基化水平差异
dmr <- DMRfinder(meth_norm, mc.p = 0.05, mc.min = 5, min.length = 100)
# 基于Bismark进行甲基化比对和差异分析
library(bismark)

# 设定基因组参考序列路径
ref_genome <- "/path/to/genome.fa"

# 设定Bismark的安装路径
bismark_dir <- "/path/to/bismark/"

# 设定Bowtie2的安装路径
bowtie2_dir <- "/path/to/bowtie2/"

# 设定样本数据的路径和前缀
m0_mc_prefix <- "/path/to/M0_MC_"
m0_unt_prefix <- "/path/to/M0_UNT_"
m1_mc_prefix <- "/path/to/M1_MC_"
m1_unt_prefix <- "/path/to/M1_UNT_"
m2_mc_prefix <- "/path/to/M2_MC_"
m2_unt_prefix <- "/path/to/M2_UNT_"

# 设定输出文件夹路径
out_dir <- "/path/to/output/"

# 设定Bismark参数
bismark_params <- "--bowtie2 --multicore 8"

# 对M0样本进行甲基化比对
bismark_genome <- bismark_genome_folder(ref_genome, bismark_dir = bismark_dir, bowtie2_dir = bowtie2_dir)
m0_mc_bam <- bismark(m0_mc_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)
m0_unt_bam <- bismark(m0_unt_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)

# 对M1样本进行甲基化比对
m1_mc_bam <- bismark(m1_mc_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)
m1_unt_bam <- bismark(m1_unt_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)

# 对M2样本进行甲基化比对
m2_mc_bam <- bismark(m2_mc_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)
m2_unt_bam <- bismark(m2_unt_prefix, bismark_genome, output_dir = out_dir, params = bismark_params)

# 合并重复样本
m0_mc_bam_merged <- merge(bam_file = c(m0_mc_bam, c("_1_bismark_bt2.bam", "_2_bismark_bt2.bam", "_3_bismark_bt2.bam")), output_dir = out_dir)
m0_unt_bam_merged <- merge(bam_file = c(m0_unt_bam, c("_1_bismark_bt2.bam", "_2_bismark_bt2.bam", "_3_bismark_bt2.bam")), output_dir = out_dir)
m1_mc_bam_merged <- merge(bam_file = c(m1_mc_bam, c("_1_bismark_bt2.bam", "_2_bismark_bt2.bam", "_3_bismark_bt2.bam")), output_dir = out_dir)
m1_unt_bam_merged <- merge(bam_file = c(m1_unt_bam, c("_1_bismark_bt2.bam", "_2_bismark_bt2
# 导入肿瘤共培养组甲基化芯片数据
shared_culture <- read.csv("shared_culture.csv")

# 导入单独培养组甲基化芯片数据
single_culture <- read.csv("single_culture.csv")
# 加载limma包
library(limma)

# 创建设计矩阵
design <- model.matrix(~0 + c(rep("Shared culture", 30), rep("Single culture", 30)))

# 合并数据
all_data <- cbind(shared_culture, single_culture)

# 执行差异分析
fit <- lmFit(all_data, design)
fit <- eBayes(fit)

# 提取差异甲基化位点
diff_genes <- topTable(fit, coef=1, adjust.method="BH", sort.by="B", number=Inf)
# 加载随机森林包
library(randomForest)

# 创建响应变量
response <- c(rep("M0", 30), rep("M1", 30), rep("M2", 30))

# 合并数据
all_data <- rbind(shared_culture, single_culture)

# 执行随机森林建模
rf_model <- randomForest(all_data, response, ntree=500, importance=TRUE)
