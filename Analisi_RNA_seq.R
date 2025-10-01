library(writexl)
library(DESeq2)
library(pheatmap)
setwd("C:/Users/amydg/Desktop/Università/Tirocinio/Tesi/Deseq/Analisi_nuova/Run1_MCF7")
count_matrix<-read.table("prev_raw_count.txt",header = T,sep = "\t")



#PULIZIA DATASET
count_matrix <- count_matrix[, !colnames(count_matrix) %in% c("Chr","Start","End","Strand","Length")]
espressione_regolare <- "^.*\\.BAM_Prev\\.(.*?)_R1_001\\.fastq\\.gzAligned.*$"
colnames(count_matrix) <- gsub(
  pattern = espressione_regolare, 
  replacement = "\\1",      
  x = colnames(count_matrix)
)


#CREAZIONE METADATA
# 1. Isola i nomi delle colonne dei campioni (tutte tranne la prima "Geneid")
column_names <- colnames(count_matrix)

# 2. Creazione del DataFrame metadata, usando column_names come base
metadata <- data.frame(
  # Imposta i nomi di riga INIZIALMENTE, usando l'ID completo (cruciale per l'allineamento)
  row.names = column_names
)

# 3. Estrazione e Creazione delle Colonne

# SampleID_Clean: Rimuove la parte replica/S-number (e.g., _1_S47)
metadata$Condition <- sub(
  pattern = "_\\d+_S\\d+$", 
  replacement = "", 
  x = column_names
)

# CellLine: Estrae la linea cellulare (es. MCF7)
metadata$CellLine <- gsub(
  pattern = "(_.*)$",
  replacement = "",
  x = column_names
)

# Treatment: Estrae l'ultima parte della stringa 'Condition' (es. EtOH o Tam)
metadata$Treatment <- sub(
  pattern = ".*_(.*)$", 
  replacement = "\\1", 
  x = metadata$Condition
)

# 4. Conversione a Fattori (NECESSARIO per DESeq2)
metadata$Condition <- factor(metadata$Condition)
metadata$Treatment <- factor(metadata$Treatment)
metadata$CellLine <- factor(metadata$CellLine)

rownames(count_matrix)<- count_matrix$Geneid
count_matrix$Geneid<-NULL 




#ANALISI DIFFERENZIALE
dds<-DESeqDataSetFromMatrix(countData=count_matrix,colData=metadata,design=~Condition)
dds <- DESeq(dds)

#PCA+HEATMWAP
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


#NORMALIZZAZIONE 
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts<-as.data.frame(normalized_counts)


#CONFRONTI
contrasts_to_run <- list(
  c("MCF7_TamR_Tam", "MCF7_TamR_EtOH"), 
  c("MCF7_WT_Tam", "MCF7_WT_EtOH"),   
  c("MCF7_TamR_Tam", "MCF7_WT_Tam"),   
  c("MCF7_TamR_EtOH", "MCF7_WT_EtOH")    
)

# Inizializzazione di una lista vuota per salvare i risultati filtrati
filtered_results_list <- list()

# 2. Cicla attraverso le coppie di confronto, esegui l'analisi e filtra
for (contrast_pair in contrasts_to_run) {
  
  # Estrai i nomi del gruppo di trattamento e del gruppo di controllo
  treatment_group <- contrast_pair[1]
  control_group <- contrast_pair[2]
  
  # Crea un nome descrittivo per il foglio del file Excel
  contrast_name <- paste0(treatment_group, "_vs_", control_group)
  
  # Esegui il contrasto e ottieni i risultati
  res <- results(dds, contrast = c("Condition", treatment_group, control_group))
  
  # 3. Filtra i risultati per padj < 0.05 e rimuovi i valori NA
  
  res_filtered <- subset(res, !is.na(padj) & padj < 0.05)
  
  # Se non ci sono geni significativamente espressi, salta il salvataggio
  if(nrow(res_filtered) == 0) {
    print(paste("Nessun gene significativamente espresso trovato per:", contrast_name))
    next
  }
  
  # Converti i risultati in un data.frame per una manipolazione più semplice
  res_df <- as.data.frame(res_filtered)
  
  # Aggiungi i nomi dei geni come una nuova colonna, prendendoli dai nomi delle righe
  res_df$gene_name <- rownames(res_df)
  
  # Riordina le colonne per mettere il nome del gene all'inizio
  res_df <- res_df[, c("gene_name", setdiff(names(res_df), "gene_name"))]
  
  # Aggiungi il dataframe filtrato alla lista dei risultati
  filtered_results_list[[contrast_name]] <- res_df
  
  # Stampa un messaggio per tracciare l'avanzamento
  print(paste("Analisi e filtraggio completati per:", contrast_name))
}

# 4. Salva i risultati filtrati in un unico file Excel

output_filename <- "risultati_geni_significativi_MCF7_run1.xlsx"

# Se la lista dei risultati non è vuota, salvala
if (length(filtered_results_list) > 0) {
  write_xlsx(filtered_results_list, path = output_filename)
  print(paste("File salvato con successo in:", output_filename))
} else {
  print("Nessun risultato da salvare. Controlla i tuoi contrasti o la soglia di significatività.")
}



























