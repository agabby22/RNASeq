
##PLOTS for the DESeq analysis to look at the change in gene expression of BTN3A1/2/3 proteins 

##change the pattern in order to plot different DESeq results
#1. for controls
#2. for 103.2 interactions terms
#3. for 20.1 interactions terms 

all_datasets <- list.files("../output_updated/genset_enrichment/", pattern = "control.*.RData", full.names = TRUE)

interesting_genes <- c("BTN3A1", "BTN3A2", "BTN3A3")

plots <- NULL

for (file in all_datasets) {
  
  load(file)
}

final_results <- final_results %>%
  rename(
    "hgnc_symbol" = "symbol"
  )


long_data <- final_results %>%
  dplyr::filter(symbol %in% interesting_genes) %>% 
  pivot_longer(cols = contains("norm_counts"),
               names_to = "sample",
               values_to = "counts"
  ) %>%
  mutate(sample = gsub("norm.counts.", "", sample)) %>% 
  left_join(file_metadata)

#for te baseline and also treatment samples of interaction terms for both 20.1 and 103.2
text <- long_data %>% 
  distinct(symbol, counts, FC, padj) %>% 
  mutate(FC = FC*1,
         padj = formatC(padj, format = "e", digits = 2),
         label = glue("FC = {round(FC, digits = 2)} \n
                   padj = {padj}"))

#for the control samples of interaction terms for both 103.2 and 20.1 
text2 <- long_data %>% 
  distinct(symbol, counts, FC, padj) %>% 
  mutate(FC = FC*-1,
         padj = formatC(padj, format = "e", digits = 2),
         label = glue("FC = {round(FC, digits = 2)} \n
                   padj = {padj}"))

#adding an interaction term by joining the "treatment type" with the "length of treatment"
long_data$interaction<-interaction(long_data$treatment, long_data$treatment_length)

####subsetting data for interaction terms --------------------------------------

# treatment sample

long_data_treatment <- long_data[long_data$interaction == "20.1.12h" | long_data$interaction == "20.1.24h", ]

plots <-   
  long_data %>% 
  #mutate(dev_stage = factor(dev_stage, levels = c("D70", "D100"))) %>%  # Adding this line to reorder the bars for the baseline
  ggplot() +
  stat_summary(
    aes(symbol, counts, fill = interaction),
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    position = position_dodge(.9)
  ) +
  geom_bar(
    aes(symbol, counts, fill = interaction),
    stat = "summary",
    fun = "mean",
    position = position_dodge(.9)
  ) +
  # Adding text labels with geom_text using filtered 'text' data with FC and p values for each change in gene expression
  geom_text(data = text %>% distinct(symbol, .keep_all = TRUE), 
            aes(symbol, counts, label = label),
            position = position_dodge(width = 0.8), 
            vjust = -3.5, hjust = 0.55, size = 3) + #Adjusting the size and position 
  geom_jitter(aes(symbol, counts, fill = dev_stage), # Adding jitter plot to visualize individual data points
              stat = "identity",
              position = position_dodge(.9)) +
  ylim(0, max(long_data$counts + long_data$counts*0.1)) + #setting y-axis limits to make the plot more readable
  scale_y_log10(limits = c(1,1e6), expand = c(0, 0)) +
  ylab(expression(-log[10]*"(gene counts)")) + #plotting -10 log of the gene counts
  theme_minimal(base_size = 16) +
  theme(axis.title.x = element_blank(), #removing the title of x axis
        plot.title = element_text(hjust = 0.5)) + #to center the plot tile 
  ggtitle("Barplots for Treatment Interaction Effect for 20.1") #adding the title for the plot to indicate what it shows

print(plots)

#save plots
setwd("/Users/gabrielaarciszewska/Desktop/Dissertation/Analysis_25sept/output_updated/DESeq_results//plots/")
ggsave("barplot_interaction_20.1.png", plots, width = 12, height = 9, bg="white")


  