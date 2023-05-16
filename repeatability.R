library(pheatmap)
# 
### 检验重复样本数据的相似性
### https://www.cnblogs.com/zhanmaomao/p/12266551.html
save_pheatmap_pdf <- function(x, filename, width = 7, height = 7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width = width, height = height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
rc <- read.table("../hisat2/total_rc.txt", header = FALSE, sep = " ")
rc_ercc <- read.table("../hisat2_ERCC/total_rc.txt", header = FALSE, sep = " ")

total_rc <- merge(rc, rc_ercc, by = "V1")

total_rc$ratio <- total_rc$V2.x / total_rc$V2.y

total_rc$sf <- total_rc$ratio / total_rc$ratio[3]


## reture gene length and count after featurecounts
fc_out <- function(df) {
    df <- df[, c(1, 7:ncol(df))]
    df <- df[which(rowSums(df[, 2:ncol(df)]) > 0), ]
    return(df)
}

df <- read.table("fc_gene.txt", header = TRUE, sep = "\t") ## reads count

df <- fc_out(df)
names(df)[2:ncol(df)] <-
   c("rpm_asy_1", "rpm_asy_2", "rpm_0_1", "rpm_0_2", "rpm_160_2",
   "rpm_300_1", "rpm_300_2", "rpm_40_1", "rpm_40_2", "rpm_80_1", "rpm_80_2")

df_2 <- df[apply(df[, 2:ncol(df)], 1, function(x) sum(x >= 1) > 10), ]

df_3 <- scale(cor(log2(df_2[, 2:ncol(df)] + 1)))

ph <- pheatmap(df_3)

save_pheatmap_pdf(ph, "cor.pdf")

## calculate the RPM ----
rpm <- function(df, total_rc){
    # for (i in total_rc$V1) {
    #     df[, i] <- df[, i] / total_rc[total_rc$V1 == i, 2] * 1000000
    # }
    for (i in 2:ncol(df)) {
        df[, i] <- df[, i] / total_rc[i - 1, 2] * 1000000 * total_rc[i - 1, 5]
    }
    return(df)
}

# total_rc <- read.table("../total_rc.txt", header = FALSE, sep = " ")


df_rpm <- rpm(df, total_rc)

df_rpm_2 <- df_rpm[
    apply(df_rpm[, 2:ncol(df_rpm)], 1, function(x) sum(x >= 1) > 10), ]

df_rpm_3 <- scale(cor(log2(df_rpm_2[, 2:ncol(df_rpm)] + 1)))

ph <- pheatmap(df_rpm_3)

save_pheatmap_pdf(ph, "cor_rpm.pdf")