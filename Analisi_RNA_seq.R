library(writexl)
library(DESeq2)
library(pheatmap)



#CLEANING DATASET
count_matrix <- count_matrix[, !colnames(count_matrix) %in% c("Chr","Start","End","Strand","Length")]
espressione_regolare <- "^.*\\.BAM_Prev\\.(.*?)_R1_001\\.fastq\\.gzAligned.*$"
colnames(count_matrix) <- gsub(
  pattern = espressione_regolare, 
  replacement = "\\1",      
  x = colnames(count_matrix)
)


#METADATA CREATION
# 1.  Selecting column names without Geneid
column_names <- colnames(count_matrix)

# 2. Creation of the metadata DataFrame, using column_names as a basis
metadata <- data.frame(
# Set row names using the full ID
  row.names = column_names
)

# 3. Extraction and Creation of Columns
 
 #Condition
metadata$Condition <- sub(
  pattern = "_\\d+_S\\d+$", 
  replacement = "", 
  x = column_names
)

# CellLine
metadata$CellLine <- gsub(
  pattern = "(_.*)$",
  replacement = "",
  x = column_names
)

# Treatment
metadata$Treatment <- sub(
  pattern = ".*_(.*)$", 
  replacement = "\\1", 
  x = metadata$Condition
)

# 4. Conversion to Factors 
metadata$Condition <- factor(metadata$Condition)
metadata$Treatment <- factor(metadata$Treatment)
metadata$CellLine <- factor(metadata$CellLine)



#Set GeneID as rownames and remove the GeneID column
rownames(count_matrix)<- count_matrix$Geneid
count_matrix$Geneid<-NULL 




#DIFFERENTIAL ANALYSIS
dds<-DESeqDataSetFromMatrix(countData=count_matrix,colData=metadata,design=~Condition)
dds <- DESeq(dds)

#Graph creation: PCA and heatmap with correlation values
# Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)
# Plot PCA 
plotPCA(rld, intgroup="Condition")
# Extract the rlog matrix from the object
rld_mat <- assay(rld)
# Compute pairwise correlation values
rld_cor <- cor(rld_mat)
# Plot heatmap
pheatmap(rld_cor)


#Normalization
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<-as.data.frame(normalized_counts)


#Definition of comparisons to be made
contrasts_to_run <- list(
  c("MCF7_TamR_Tam", "MCF7_TamR_EtOH"), 
  c("MCF7_WT_Tam", "MCF7_WT_EtOH"),   
  c("MCF7_TamR_Tam", "MCF7_WT_Tam"),   
  c("MCF7_TamR_EtOH", "MCF7_WT_EtOH")    
)

# Initialising an empty list to save filtered results
filtered_results_list <- list()

#2. Cycle through comparison pairs, perform analysis, and filter
for (contrast_pair in contrasts_to_run) {
  
  # Extract the names of the treatment group and the control group
  treatment_group <- contrast_pair[1]
  control_group <- contrast_pair[2]
  
  # Create a descriptive name for the Excel file sheet
  contrast_name <- paste0(treatment_group, "_vs_", control_group)
  
  # Run the comparison and get the results
  res <- results(dds, contrast = c("Condition", treatment_group, control_group))
  
# 3. Filter results for padj < 0.05 and remove NA values  
  res_filtered <- subset(res, !is.na(padj) & padj < 0.05)
  
 # If there are no significantly expressed genes, skip saving
  if(nrow(res_filtered) == 0) {
    print(paste("Nessun gene significativamente espresso trovato per:", contrast_name))
    next
  }
  
  # Convert the results into a data.frame
  res_df <- as.data.frame(res_filtered)
  
  # Add the names of the genes
  res_df$gene_name <- rownames(res_df)
  
  # Reorder the columns to put the gene name at the beginning
  res_df <- res_df[, c("gene_name", setdiff(names(res_df), "gene_name"))]
  
  # Add the filtered dataframe to the list of results
  filtered_results_list[[contrast_name]] <- res_df
  
  # Print a message to track progress
  print(paste("Analisi e filtraggio completati per:", contrast_name))
}

# 4. Save the filtered results in a single Excel file

output_filename <- "risultati_geni_significativi_MCF7_run1.xlsx"

# If the list of results is not empty, save it.
if (length(filtered_results_list) > 0) {
  write_xlsx(filtered_results_list, path = output_filename)
  print(paste("File salvato con successo in:", output_filename))
} else {
  print("Nessun risultato da salvare. Controlla i tuoi contrasti o la soglia di significatività.")
}



























