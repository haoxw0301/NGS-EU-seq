library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
# pwd: /NFSdata01/YBT/RNA-seq/L1_subfamily_diff/Cellcycle_EU/map/fc_L1
## filter out the featureCounts results
fc_out <- function(df) {
    df <- df[, c(1, 7:ncol(df))]
    df <- df[which(rowSums(df[, 2:ncol(df)]) > 0), ]

    return(df)
}

# load data ----
rc <- read.table("../hisat2/total_rc.txt", header = FALSE, sep = " ")
rc_ercc <- read.table("../hisat2_ERCC/total_rc.txt", header = FALSE, sep = " ")

total_rc <- merge(rc, rc_ercc, by = "V1")

total_rc$ratio <- total_rc$V2.x / total_rc$V2.y
total_rc$sf <- total_rc$ratio / total_rc$ratio[3]

df <- read.table("fc_L1.txt", header = TRUE, sep = "\t")

df <- fc_out(df)

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

df_rpm <- rpm(df, total_rc)

## 根据样本修改列名
names(df_rpm)[2:ncol(df_rpm)] <-
   c("rpm_asy_1", "rpm_asy_2", "rpm_0_1", "rpm_0_2", "rpm_160_2",
   "rpm_300_1", "rpm_300_2", "rpm_40_1", "rpm_40_2", "rpm_80_1", "rpm_80_2")

## 根据样本修改列名
## clean data ----
df2 <- df_rpm %>%
    transmute(L1 = Geneid,
        rpm_0 = (rpm_0_1 + rpm_0_2) / 2,
        rpm_40 = (rpm_40_1 + rpm_40_2) / 2,
        rpm_80 = (rpm_80_1 + rpm_80_2) / 2,
        rpm_160 = rpm_160_2,
        rpm_300 = (rpm_300_1 + rpm_300_2) / 2,
        rpm_asy = (rpm_asy_1 + rpm_asy_2) / 2
    )
## annotation files ----
young_L1 <- read.table("~/haoxw/reference/mm10/L1xn/young_L1.txt")
A_L1 <- read.table("~/haoxw/reference/mm10/L1xn/A_L1_FC.txt")
genic_L1 <- read.table("~/haoxw/reference/mm10/L1xn/genic_L1_FC.txt")
proximal_L1 <- read.table("~/haoxw/reference/mm10/L1xn/proximal_itg_L1_FC.txt")

df2$age <- "old"
df2$age[df2$L1 %in% young_L1$V1] <- "young"

df2$compartment <- "B"
df2$compartment[df2$L1 %in% A_L1$V1] <- "A"

df2$region <- "distal_intergenic"
df2$region[df2$L1 %in% genic_L1$V1] <- "genic"
df2$region[df2$L1 %in% proximal_L1$V1] <- "proximal_intergenic"

df2$region <- factor(df2$region,
    levels = c("genic", "proximal_intergenic", "distal_intergenic"),
    labels = c("genic", "proximal intergenic", "distal intergenic"),
    ordered = TRUE)

df3 <- melt(df2, id.vars = c("L1", "age", "compartment", "region"),
    variable.name = "stage", value.name = "rpm")
    
## 根据样本修改列名
df3$stage <- factor(df3$stage,
    levels = c("rpm_0", "rpm_40", "rpm_80", "rpm_160", "rpm_300", "rpm_asy"),
    labels = c("0", "40", "80", "160", "300", "asynaptic"),
    ordered = TRUE)


## plot ----
### single expression level ----
pdf("L1_rpm_violin.pdf", width = 20, height = 6)
ggplot(df3, aes(x = stage, y = rpm, fill = age)) +
    geom_violin(scale = "width") +
    geom_boxplot(width = 0.1, color = "black", alpha = 0.5,
        position = position_dodge(0.9),
        lwd = 0.5, outlier.colour = NA) +
    facet_grid(compartment ~ region) +
    scale_fill_manual(values = c("orange", "#9933CC")) +
    scale_y_continuous(trans = "log2") +
    # ylim(0, 2) +
    labs(x = "Stage", y = "RPM", fill = "") +
    theme_bw(base_rect_size = 1.3) +
    theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 20)
    )
dev.off()

### total expression level ----
df3_total <- df3 %>%
    group_by(stage, age, compartment, region) %>%
    summarise(rpm = sum(rpm))
pdf("L1_rpm_col.pdf", width = 20, height = 5)
ggplot(df3_total, aes(x = stage, y = rpm, fill = age)) +
    geom_col(width = 0.5, position = position_dodge(), color = "grey") +
    facet_grid(compartment ~ region) +
    scale_fill_manual(values = c("orange", "#9933CC")) +
    theme_bw(base_rect_size = 1.3) +
    labs(x = "Stage", y = "RPM", fill = "") +
    theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 20)
    )
dev.off()

pdf("L1_rpm_total.pdf", width = 4, height = 4)
ggplot(df3_total, aes(x = stage, y = rpm, fill = age)) +
    geom_col(width = 0.6, position = position_dodge(), color = NA) +
    # facet_grid(compartment ~ region) +
    scale_fill_manual(values = c("orange", "#9933CC")) +
    theme_classic(base_rect_size = 1.3) +
    labs(x = "", y = "RPM", fill = "") +
    theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1,
            size = 18, color = "black"))
dev.off()
### number of transcribed L1 ----
df3_number <- df3[which(df3$rpm > 0), ] %>%
    group_by(stage, age, compartment, region) %>%
    summarise(n = length(rpm))

pdf("L1_number_col.pdf", width = 20, height = 5)
ggplot(df3_number, aes(x = stage, y = n, fill = age)) +
    geom_col(width = 0.5, position = position_dodge(), color = "grey") +
    facet_grid(compartment ~ region) +
    scale_fill_manual(values = c("orange", "#9933CC")) +
    theme_bw(base_rect_size = 1.3) +
    labs(x = "Stage", y = "# of L1", fill = "") +
    theme(
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 16, color = "black"),
        legend.text = element_text(size = 18),
        strip.text = element_text(size = 20)
    )
dev.off()