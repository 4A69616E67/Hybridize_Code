library(circlize)
chr_size_file <- "../../../Xenopus/XTrXLa.chrom.sizes"
TeLs_coverage_file <- "../WGS/mapping/GTeLs_7.1M.RPKM.bedgraph"
P53Ls_coverage_file <- "../WGS/mapping/GP53eLsM_12.1M.RPKM.bedgraph"
# chromosome szie file must cotain two columns
chr_size <- read.table(chr_size_file, colClasses = c("character", "numeric"))
chr_size <- data.frame("chr" = chr_size$V1, "start" = 1, "end" = chr_size$V2)
chr_size$chr <- factor(levels = chr_size$chr, chr_size$chr)

## read genomeic coverage file, 1M resolution
TeLs_wgs <- read.table(TeLs_coverage_file, colClasses = c("character", "numeric", "numeric", "numeric"))
P53Ls_wgs <- read.table(P53Ls_coverage_file, colClasses = c("character", "numeric", "numeric", "numeric"))
## modify the column name
names(TeLs_wgs) <- c("chr", "start", "end", "score")
names(P53Ls_wgs) <- c("chr", "start", "end", "score")
## different colours used in different type of chromosome
TeLs_wgs$color <- "#8FCF69"
TeLs_wgs$color[TeLs_wgs$chr %in% paste("Chr", c(1:8, "9_10"), "S", sep = "")] <- "#EFC165"
TeLs_wgs$color[TeLs_wgs$chr %in% paste("Chr", c(1:8, "9_10"), "L", sep = "")] <- "#75B5F0"

P53Ls_wgs$color <- "#8FCF69"
P53Ls_wgs$color[P53Ls_wgs$chr %in% paste("Chr", c(1:8, "9_10"), "S", sep = "")] <- "#EFC165"
P53Ls_wgs$color[P53Ls_wgs$chr %in% paste("Chr", c(1:8, "9_10"), "L", sep = "")] <- "#75B5F0"

###### generate data which will be plot (plot all chromosome in this case)
plot_chr <- c(paste("Chr", 1:10, sep = ""), paste("Chr", c(1:8, "9_10"), "S", sep = ""), paste("Chr", c(1:8, "9_10"), "L", sep = ""))
plot_chr_szie <- chr_size[chr_size$chr %in% plot_chr, ]
plot_chr_szie$chr <- factor(levels = plot_chr, plot_chr_szie$chr)

plot_TeLs_wgs <- TeLs_wgs[TeLs_wgs$chr %in% plot_chr, ]
plot_P53Ls_wgs <- P53Ls_wgs[P53Ls_wgs$chr %in% plot_chr, ]

circos.clear() # clear the circos settings
## circos setting
circos.par("track.height" = 0.1, "start.degree" = 90, "cell.padding" = c(0, 0, 0, 0), "points.overflow.warning" = FALSE)
circos.genomicInitialize(plot_chr_szie, sector.names = gsub("Chr", "", plot_chr)) # draw chromosome without "Chr" prefix

# draw tracks one
circos.track(plot_chr_szie$chr, ylim = c(-0.4, 0), bg.border = NA, bg.col = "white")
circos.trackLines(
    sectors = plot_TeLs_wgs$chr, x = (plot_TeLs_wgs$start + plot_TeLs_wgs$end) / 2,
    y = -plot_TeLs_wgs$score, type = "s", area = T, col = plot_TeLs_wgs$color, lwd = 0.5, baseline = "top"
)
# draw tracks two
circos.track(plot_chr_szie$chr, ylim = c(-0.4, 0), bg.border = NA, bg.col = "white")
circos.trackLines(
    sectors = plot_P53Ls_wgs$chr, x = (plot_P53Ls_wgs$start + plot_P53Ls_wgs$end) / 2,
    y = -plot_P53Ls_wgs$score, type = "s", area = T, col = plot_P53Ls_wgs$color, lwd = 0.5, baseline = "top"
)