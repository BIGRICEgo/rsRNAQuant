suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(tidyr)
    library(DESeq2)
    library(ggplot2)
    library(RColorBrewer)
})
library(jsonlite)


option_list = list(
  make_option(c("-i", "--input"), type = "character", help = "GSE id"),
  make_option(c("-o", "--output"), type = "character", help = "Output directory for plots and data"),
  make_option(c("-s", "--samples"), type = "character", help = "Comma-separated list of sample names"),
  make_option(c("-g", "--groups"), type = "character", help = "String representation of the groups dictionary")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

GSE_ID = opt$input
BASIC_OUTPUT_DIR <- opt$output
SAMPLES <- strsplit(opt$samples, ",")[[1]]
GROUPS <- fromJSON(opt$groups)
OUTPUT_DIR <- file.path(BASIC_OUTPUT_DIR, "plots")
case_samples <- GROUPS$case
control_samples <- GROUPS$control
sample_names <- SAMPLES


Class <- sapply(SAMPLES, function(s) {
  if (s %in% case_samples) {
    "case"
  } else if (s %in% control_samples) {
    "control"
  }
})

##### Differential Analysis
# Enough samples are needed
if ( length(case_samples) >= 3 && length(control_samples)  >=3) {
  print("Start analyzing ....")

  matrix_file_path <- file.path(BASIC_OUTPUT_DIR, "stats", paste0(GSE_ID, "_rRNA_count_matrix.txt"))
  expr <- read.delim(matrix_file_path, header = TRUE, sep = "\t")

  df_matrix <- expr[, c("Sequence", sample_names)]
  rownames(df_matrix) <- df_matrix[,c("Sequence")]
  df_matrix <- df_matrix[,-1]
  countData <- df_matrix[rowMeans(df_matrix)>10,] 

  condition <- factor(Class)
  condition <- factor(Class, levels = c("control", "case"))
  colData <- data.frame(row.names=colnames(countData), condition)
  dds <- DESeqDataSetFromMatrix(countData = round(countData), colData = colData, design = ~ condition)
  #head(dds)
  dds1 <- DESeq(dds) 

  res <- results(dds1)
  #View(res)
  summary(res)
  res_data <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
  res_data$color <- ifelse(res_data$padj<0.05 & abs(res_data$log2FoldChange)>= 1,ifelse(res_data$log2FoldChange > 1,'red','blue'),'gray')

  res_data$Sequence <- rownames(res_data)
  res_data <- res_data[, c(ncol(res_data), 1:(ncol(res_data)-1))]
  res_data_noNA <- res_data[!is.na(res_data$pvalue), ]
  res_data_noNA <- res_data[!is.na(res_data$padj), ]
  write.table(res_data_noNA, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_DA_result_noNA.txt"), sep="\t", quote=FALSE, row.names=FALSE)

  n_total <- nrow(res_data_noNA)
  n_up <- sum(res_data_noNA$color == "red", na.rm = TRUE)
  prop_up <- round(n_up / n_total * 100, 2)
  n_down <- sum(res_data_noNA$color == "blue", na.rm = TRUE)
  prop_down <- round(n_down / n_total * 100, 2)
  n_gray <- sum(res_data_noNA$color == "gray", na.rm = TRUE)
  prop_gray <- round(n_gray / n_total * 100, 2)

  summary_df <- data.frame(
    "Info" = c(
      "rsRNA_total",
      "adjusted p-value",
      "LFC >0 (Up,red)",
      "LFC >0 (Up,blue)",
      "Gray"
    ),
    "Value" = c(
      n_total,
      "0.05",
      paste0(n_up, " (", prop_up, "%)"),
      paste0(n_down, " (", prop_down, "%)"),
      paste0(n_gray, " (", prop_gray, "%)")
    ),
    stringsAsFactors = FALSE
  )
  write.table(summary_df, file = paste0(OUTPUT_DIR, "/", GSE_ID, "_DA_summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE, fileEncoding = "UTF-8")
  message(paste0("Saved Summary File && DA_result File  to ", OUTPUT_DIR))



  ## Plot
  color <- c(red = "red",gray = "#676662",blue = "steelblue") # "gray"
  p1<- ggplot(
    res_data_noNA, aes(log2FoldChange, -log10(padj), col = color, show.legend = F)) +  
    geom_jitter(width = 0.2, height = 0, alpha = 0.7) + # width控制横向抖动幅度
    theme_bw() +
    scale_color_manual(values = color) +
    labs(x="log2 (fold change)",y="-log10 (padj)") +
    geom_hline(yintercept = -log10(0.05), lty=4,col="gray",lwd=0.6) +
    geom_vline(xintercept = c(-1, 1), lty=4,col="gray",lwd=0.6) +
    theme(#legend.position = "none",
      panel.grid=element_blank(),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 10) )+
    theme(plot.title = element_text(size = 10,hjust = 0.5))


  ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_VolcanoPlot.pdf",sep=""), plot = p1,width = 4, height = 3, dpi = 300)
  ggsave(filename = paste(OUTPUT_DIR, "/",GSE_ID, "_VolcanoPlot.png",sep=""), plot = p1,width = 4, height = 3, dpi = 300)
  message(paste0("Saved Volcano Plot to ", OUTPUT_DIR))



} else {
  message("Samples of case and control need more than 3 or just 3, case=", length(case_samples), ", control=", length(control_samples))
}


##### Find Biomarkers
if(
    length(case_samples) > 3 && length(control_samples) > 3
){
  print("Start finding Biomarkers....")

  # Select DE rsRNAs
  res_data_noNA <- res_data_noNA[, !(colnames(res_data_noNA) %in% c("lfcSE", "stat"))]
  res_data_colored <- res_data_noNA[res_data_noNA$color != "gray", ]

  rpm_file_path <- file.path(BASIC_OUTPUT_DIR, "stats", paste0(GSE_ID, "_rRNA_expression_matrix.txt"))
  expr_rpm <- read.delim(rpm_file_path, header = TRUE, sep = "\t")
  res_df <- merge(res_data_colored, expr_rpm, by = "Sequence", all.x=TRUE)


  suppressPackageStartupMessages({
      library(caret)
      library(glmnet)
      library(MASS)
      library(pROC)
      library(randomForest)
      library(ggplot2)
      library(ggrepel)

  })


  x <- res_df[,sample_names]
  rownames(x) <- res_df$Sequence
  x <- t(x) 
  class(x) 
  y <- Class

  # Stratified sampling
  set.seed(111)
  index <- createDataPartition(y, p = 0.7, list = FALSE, times = 1)
  head(index)

  x_scaled <- x
  # Separate train and test data
  x_train <- x_scaled[index, ]
  x_test <- x_scaled[-index, ]
  dim(x_train)
  dim(x_test)

  y_train <- y[index]
  y_test <- y[-index]
  y_train <- ifelse(y_train == "case", 1, 0)
  y_test <- ifelse(y_test == "case", 1, 0)

  # lambda Path
  la <- glmnet(x_train, y_train, family="binomial", intercept = F, alpha=1)
  p1 <- plot(la, xvar = "lambda", lwd=1.5,  label = F)
  fit <- cv.glmnet(x = x_train, y = y_train, nfolds = 3, family = "binomial", alpha = 1)
  p2 <- plot(fit)
  # Check lambda.min & lambda.1se
  #print(fit$lambda.min)
  #print(fit$lambda.1se)
  best_lambda <- fit$lambda.min
  # find coefficients of best model
  best_model <- glmnet(x = x_train, y = y_train, alpha = 1, lambda = best_lambda, family = "binomial")
  coefficients <- coef(best_model)
  cf_info <- coefficients@x
  names(cf_info) <- coefficients@Dimnames[[1]][coefficients@i+1]

  sf_train <- x_train[,as.logical(as.numeric(coefficients[-1, ]))] 
  cf_length <- length(cf_info)
  # Group signature
  if (cf_length == 2) {
      group <- sf_train*cf_info[2]+cf_info[1][1]
      sf_train_group <- as.matrix(cbind(as.data.frame(sf_train),group))
      colnames(sf_train_group) <- c(names(cf_info)[2], "group")
  } else {
      group <- as.numeric(cf_info[1]) + as.matrix(sf_train[, 1:(cf_length-1)]) %*% as.numeric(cf_info[2:cf_length])
      names(group) <- 'group'
      sf_train_group <- as.matrix(cbind(as.data.frame(sf_train),group))
  }

  # Feature Importance
  rf_model <- randomForest(sf_train_group, as.factor(y_train), ntree = 500,importance = TRUE,proximity = TRUE)
  p3 <- varImpPlot(rf_model)

  pdf(paste(OUTPUT_DIR, "/",GSE_ID,  "_Lasso.pdf",sep=""), height = 5, width = 10)
  plot(la, xvar = "lambda", lwd=1.5,  label = F)
  plot(fit)
  dev.off()

  png(paste(OUTPUT_DIR, "/",GSE_ID,  "_Lasso.png",sep=""), height = 5, width = 10, units = "in", res = 300)
  par(mfrow = c(1, 2))
  plot(la, xvar = "lambda", lwd=1.5,  label = F)
  plot(fit)
  dev.off()

  pdf(paste(OUTPUT_DIR, "/",GSE_ID,  "_RF.pdf",sep=""), height = 5, width = 20)
  varImpPlot(rf_model, cex.names = 0.3, main="Variable Importance Plot")
  dev.off()
  
  png(paste(OUTPUT_DIR, "/",GSE_ID,  "_RF.png",sep=""), height = 5, width = 20, units = "in", res = 300)
  varImpPlot(rf_model, cex.names = 0.3, main="Variable Importance Plot")
  dev.off()


  # Expression of selected features 
  selected_features <- coefficients@Dimnames[[1]][coefficients@i+1]
  df_selected_features <- selected_features[-1]
  df_selected_features <- res_df %>% filter(Sequence %in% c(df_selected_features)) 
  print(df_selected_features)
  # Add cofficients 
  df_selected_features$cofficients <- cf_info[-1]
  write.table(df_selected_features, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_lasso_selected_features_rpm.txt"), sep="\t", quote=FALSE, row.names=FALSE)
  write.table(cf_info, file=paste0(OUTPUT_DIR, "/",GSE_ID, "_lasso_selected_features_cofficients.txt"), sep="\t", quote=FALSE, row.names=FALSE)

      

  # ROC
  sf_test <- x_test[,as.logical(as.numeric(coefficients[-1, ]))] 
  if (cf_length == 2) {
      group <- sf_test*cf_info[2]+cf_info[1][1]
      message(" Select only one sequence as biomarker.")
      sf_test_group <- as.matrix(cbind(as.data.frame(sf_test),group))
      colnames(sf_test_group) <- c(names(cf_info)[2], "group")
  } else {
      group <- as.numeric(cf_info[1]) + as.matrix(sf_test[, 1:(cf_length-1)]) %*% as.numeric(cf_info[2:cf_length])
      sf_test_group <- as.matrix(cbind(as.data.frame(sf_test),group))
  }


  roc_data_list <- list()
  auc_values_df <- data.frame(Feature = character(), AUC = numeric(), stringsAsFactors = FALSE)


  for (col_name in colnames(sf_test_group)) {
    roc_obj <- roc(y_test, sf_test_group[, col_name])
    auc_value <- auc(roc_obj)
    label_with_auc <- paste0(col_name, " (AUC=", round(auc_value, 3), ")")
    roc_data_list[[label_with_auc]] <- roc_obj
    auc_values_df <- rbind(auc_values_df, data.frame(Feature = col_name, AUC = auc_value))
  }

  write.table(auc_values_df, 
              file = paste0(OUTPUT_DIR, "/", GSE_ID, "_ROC_AUC_values.txt"), 
              sep = "\t", 
              quote = FALSE, 
              row.names = FALSE)

  n_colors <- length(colnames(sf_test_group))
  my_colors <- colorRampPalette(c("red", "purple", "blue"))(n_colors)

  p_ROC <- ggroc(roc_data_list, legacy.axes = TRUE) +
    scale_color_manual(values = my_colors) +
    theme_bw() +
    guides(color = guide_legend(title = "Features", ncol = 1)) + # 给legend加一个标题
    theme(legend.position = "right")+
    annotate(geom = "segment", x = 1, y = 0, xend = 0, yend = 1, col = "gray", linetype = 6)

  # Adjust the size of the legend if there are many features
  if (length(roc_data_list) > 8) {
    p_ROC <- p_ROC + theme(legend.text = element_text(size = 6), 
                          legend.key.size = unit(0.5, "cm")) 
  }

  ggsave(paste(OUTPUT_DIR, "/", GSE_ID, "_ROC.pdf", sep=""), p_ROC, width = 10, height = 8)
  ggsave(paste(OUTPUT_DIR, "/", GSE_ID, "_ROC.png", sep=""), p_ROC, width = 10, height = 8)



} else {
  message("Conditions: Samples of case and control need more than 3, case =", length(case_samples), ", control=", length(control_samples))
}

message("Plots Finished.")

