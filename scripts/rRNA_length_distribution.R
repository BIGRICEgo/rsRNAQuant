suppressPackageStartupMessages(library(ggforce))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(patchwork))
library(jsonlite)


option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "GSE id"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for plots and data"),
  make_option(c("-s", "--samples"), type = "character", help = "Comma-separated list of sample names"),
  make_option(c("-g", "--groups"), type = "character", help = "String representation of the groups dictionary")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## Get input parameters
GSE_ID = opt$input
BASIC_OUTPUT_DIR <- opt$output
SAMPLES <- strsplit(opt$samples, ",")[[1]]
GROUPS <- fromJSON(opt$groups)
OUTPUT_DIR <- file.path(BASIC_OUTPUT_DIR, "plots")

sample_names <- SAMPLES
case_samples <- GROUPS$case
control_samples <- GROUPS$control

## Load Files
matrix_file_path <- file.path(BASIC_OUTPUT_DIR, "stats",paste0(GSE_ID, "_rRNA_expression_matrix.txt"))
expr <- read.delim(matrix_file_path, header = TRUE, sep = "\t")

##### Length Distribution 
out <- expr %>%
  group_by(Length) %>%
  summarise(across(all_of(sample_names), sum, na.rm=TRUE)) %>%
  ungroup()

out <- as.data.frame(out)

out$Class <- "rRNA"
out$sum_case <- rowSums(out[, case_samples, drop=FALSE])
out$sum_control <- rowSums(out[, control_samples, drop=FALSE])
out <- out[, c("Class", "Length", sample_names, "sum_case", "sum_control")]
out <- out[order(out$Length), ]
write.table(out, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_rRNA_Length_distribution_RPM.txt"), sep="\t", quote=FALSE, row.names=FALSE)

plot_data <- out[ ,c("Class","Length","sum_case","sum_control")]
data_plot <- pivot_longer(
  plot_data,
  cols = c("sum_case", "sum_control"),
  names_to = "Type",
  values_to = "RPM"
) %>%
  mutate(Type = recode(Type, sum_case = "case", sum_control = "control"))

data_plot <- as.data.frame(data_plot)

## Plot
p1 <- ggplot(data_plot, aes(x = Length, y = RPM, fill = Type)) + 
  geom_bar(stat = "identity",position = "dodge") + 
  #facet_grid(. ~ Type) +
  labs(x = "length", y = "RPM", title = "") + 
  theme(axis.line = element_line()
        , panel.grid.minor = element_blank() 
        , panel.grid.major = element_blank()
        , panel.background = element_rect(fill = "transparent", colour = NA)
        , plot.background = element_rect(fill = "transparent", colour = NA)
        , panel.border = element_rect(fill = "transparent")
        , text = element_text(size = 10, color = "black")
        , axis.text.x = element_text(size = 7, color = "black")
        , axis.text.y = element_text(size = 10, color = "black")					
  ) +
  scale_fill_brewer(palette = "Paired") +
  theme(legend.position = "right") +
  guides(fill = guide_legend(nrow = 2, bycol = TRUE)) +
  labs(fill = "") 

ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_rRNA_Length_Distribution.pdf",sep=""), plot = p1,width = 10, height = 3, dpi = 300)
ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_rRNA_Length_Distribution.png",sep=""), plot = p1,width = 10, height = 3, dpi = 300)

message(paste0("Saved rRNA Length Distribution to ", OUTPUT_DIR))






##### HollowPiechart
out <- expr %>%
  group_by(Subtype) %>%
  summarise(across(all_of(sample_names), sum, na.rm=TRUE)) %>%
  ungroup()

out <- as.data.frame(out)


out$sum_case <- rowSums(out[, case_samples, drop=FALSE])
out$sum_control <- rowSums(out[, control_samples, drop=FALSE])
out <- out[, c("Subtype", sample_names, "sum_case", "sum_control")]

## Order
subtype_order <- c("5S","5.8S","18S","28S","45S","12S","16S")
out <- out[order(subtype_order), ]
out$Subtype <- factor(out$Subtype,levels =subtype_order )

write.table(out, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_rRNA_subtype_RPM.txt"), sep="\t", quote=FALSE, row.names=FALSE)


## Plot
data_plot_case <- out[ ,c("Subtype","sum_case")]
data_plot_case$percentage <- round(data_plot_case$sum_case/sum(data_plot_case$sum_case)* 100, 2)
data_plot_case$label <- paste0(data_plot_case$Subtype, " ","(",data_plot_case$percentage, "%",")")
data_plot_case$label <- factor(data_plot_case$label, levels = data_plot_case$label[order(match(data_plot_case$Subtype, subtype_order))])

# Colors
my_colors <- c("#a6cee3","#2078b4","#b2df8a","#32a02e",
               "#fc9b9a","#e3191a","#febe6f","#ff8000","#cab2d5")
names(my_colors) <- data_plot_case$label


p2 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+labs(title="case")+
  geom_arc_bar(data=data_plot_case,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=percentage,fill=label), color = "white")+
  scale_fill_manual(values = my_colors)


data_plot_control <- out[ ,c("Subtype","sum_control")]
data_plot_control$percentage <- round(data_plot_control$sum_control/sum(data_plot_control$sum_control)* 100, 2)
data_plot_control$label <- paste0(data_plot_control$Subtype, " ","(",data_plot_control$percentage, "%",")")
data_plot_control$label <- factor(data_plot_control$label, levels = data_plot_control$label[order(match(data_plot_control$Subtype, subtype_order))])

# Colors
my_colors <- c("#a6cee3","#2078b4","#b2df8a","#32a02e",
               "#fc9b9a","#e3191a","#febe6f","#ff8000","#cab2d5")
names(my_colors) <- data_plot_control$label

p3 <- ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+
  xlab("")+ylab('')+labs(title="control")+
  geom_arc_bar(data=data_plot_control,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=percentage,fill=label), color = "white")+
  scale_fill_manual(values = my_colors)


combined_plot <- p2 + p3 + plot_layout(ncol = 2)
ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_rRNA_subtype_HollowPieChart.pdf",sep=""), plot = combined_plot, width = 12, height = 5, dpi = 300)
ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_rRNA_subtype_HollowPieChart.png",sep=""), plot = combined_plot, width = 12, height = 5, dpi = 300)
message(paste0("Saved rRNA subtype HollowPie Chart to ", OUTPUT_DIR))





##### End Proportion
out <- expr %>%
  group_by(Fragment) %>%
  summarise(across(all_of(sample_names), sum, na.rm=TRUE)) %>%
  ungroup()

out <- as.data.frame(out)

out$sum_case <- rowSums(out[, case_samples, drop=FALSE])
out$sum_control <- rowSums(out[, control_samples, drop=FALSE])
out <- out[, c("Fragment", sample_names, "sum_case", "sum_control")]

# Order
End_order <- c("5_end", "inner", "3_end")
out <- out[order(End_order), ]
out$Fragment <- factor(out$Fragment,levels = End_order)
write.table(out, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_rRNA_End_RPM.txt"), sep="\t", quote=FALSE, row.names=FALSE)


plot_data <- out[ ,c("Fragment","sum_case","sum_control")]
data_plot <- pivot_longer(
  plot_data,
  cols = c("sum_case", "sum_control"),
  names_to = "Type",
  values_to = "RPM"
) %>%
  mutate(Type = recode(Type, sum_case = "case", sum_control = "control"))

data_plot <- as.data.frame(data_plot)


## Plot
p4 <- ggplot(data = data_plot, aes(x="", y=RPM, fill=Fragment))+ 
  geom_bar(stat = "identity",
           position = "fill")+ 
  scale_y_continuous(expand = expansion(mult=c(0.01,0.1)),
                     labels = scales::percent_format())+ 
  theme(panel.background = element_blank(),
        axis.line = element_line(),
        legend.position = "right",
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank())+ 
  labs(x=NULL,y="End (%)",title="rRNA Fragement End Proportion") + 
  scale_fill_manual(values = c("#493f8a", "#b0b0ca", "#599f89" ) )+ 
  facet_grid(. ~ Type) 

ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID,  "_rRNA_End_StackBarPlot.pdf",sep=""), plot = p4 ,width = 6, height = 5, dpi = 300)
ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID,  "_rRNA_End_StackBarPlot.png",sep=""), plot = p4 ,width = 6, height = 5, dpi = 300)

message(paste0("Saved rRNA End StackBarPlot to ", OUTPUT_DIR))







##### Co-Expression
## More than 4 samples are needed for spearman analysis.
if(
    length(case_samples) > 4 && length(control_samples) > 4
){

message("Starting analyze rRNA co-expression correlation...")

out <- expr %>%
  group_by(Subtype) %>%
  summarise(across(all_of(sample_names), sum, na.rm=TRUE)) %>%
  ungroup()

out <- as.data.frame(out)
subtype_order <- c("5S","5.8S","18S","28S","45S","12S","16S")
out <- out[order(subtype_order), ]
out$Subtype <- factor(out$Subtype,levels =subtype_order )
spearman_data <- as.data.frame(t(out[,-1]))
colnames(spearman_data) <- t(out[,1])

## Plot
suppressPackageStartupMessages(library(corrplot))
suppressPackageStartupMessages(library(Hmisc))
col2 = colorRampPalette(c('#213d76',  'white', '#d0343d'))

spearman_data_case <- as.matrix(spearman_data[grep("case",rownames(spearman_data)), ])
r_case <- cor(spearman_data_case, method = "spearman")
r_case[is.na(r_case)] <- 0
write.table(r_case, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_CoExpression_Spearman_case.txt"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


spearman_data_control <- as.matrix(spearman_data[grep("control",rownames(spearman_data)), ])
r_control <- cor(spearman_data_control, method = "spearman")
r_control[is.na(r_control)] <- 0

# Unify order
p5 <- corrplot(r_case, order = "AOE", tl.col = "black", tl.srt = 45)
custom_order <- colnames(p5$corr)
r_control <- r_control[
  match(custom_order, rownames(r_control)),
  match(custom_order, colnames(r_control))
]
write.table(r_control, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_CoExpression_Spearman_control.txt"), sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)


pdf(paste(OUTPUT_DIR, "/",GSE_ID,  "_rRNA_coExpression_Corrplot.pdf",sep=""), height = 5, width = 10)
par(mfrow = c(1, 2))
corrplot(r_case, order = "AOE", tl.col = "black", tl.srt = 45,col = col2(10), title = "Case", mar = c(0, 0, 1, 0))
corrplot(r_control,  tl.col = "black", tl.srt = 45,col = col2(10), title = "Control", mar = c(0, 0, 1, 0))
dev.off()

message(paste0("Saved rRNA CoExpression Plot to ", OUTPUT_DIR))

png(filename = paste0(OUTPUT_DIR, "/", GSE_ID, "_rRNA_coExpression_Corrplot.png"), width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
corrplot(r_case, order = "AOE", tl.col = "black", tl.srt = 45,col = col2(10), title = "Case", mar = c(0, 0, 1, 0))
corrplot(r_control,  tl.col = "black", tl.srt = 45,col = col2(10), title = "Control", mar = c(0, 0, 1, 0))
dev.off()

}

message("Plots finished.")
