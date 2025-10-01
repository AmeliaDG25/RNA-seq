# 🧬 Analisi di Espressione Differenziale (DEG) su Linee Cellulari MCF7

Questo repository contiene lo script R completo per l'analisi dei dati di sequenziamento RNA (RNA-seq) focalizzata sulla resistenza al Tamoxifene nel cancro al seno.

L'analisi è stata eseguita utilizzando il pacchetto **DESeq2** per identificare i geni differenzialmente espressi (DEG) tra i gruppi sperimentali.

---

## 🔬 Contesto Sperimentale: Linee Cellulari e Trattamenti

L'esperimento si concentra sulla comprensione della resistenza al Tamoxifene attraverso quattro condizioni distinte:

| Linea Cellulare | Descrizione | Trattamento | Ruolo |
| :--- | :--- | :--- | :--- |
| **MCF7 WT** | Linea cellulare di **adenocarcinoma mammario umano** (sensibile al Tamoxifene). | **Tam** (Tamoxifene) | Indagine della risposta standard al farmaco. |
| **MCF7 WT** | Linea cellulare di **adenocarcinoma mammario umano** (sensibile al Tamoxifene). | **EtOH** (Etanolo) | Controllo veicolo per isolare l'effetto del farmaco. |
| **MCF7 TamR** | Sub-clone delle MCF7 WT che ha sviluppato **resistenza al Tamoxifene**. | **Tam** (Tamoxifene) | Indagine dei meccanismi di resistenza e vie alternative. |
| **MCF7 TamR** | Sub-clone delle MCF7 WT che ha sviluppato **resistenza al Tamoxifene**. | **EtOH** (Etanolo) | Controllo veicolo per la linea resistente. |

---

## 🛠 Requisiti e Pipeline di Analisi

Lo script R in questa repository esegue l'intera pipeline di analisi DEG: dalla preparazione dei dati al salvataggio dei risultati significativi.

### Pacchetti R Richiesti

Lo script si basa sui seguenti pacchetti R (CRAN e Bioconductor):

* **`DESeq2`**: Core package per l'analisi di espressione differenziale.
* **`pheatmap`**: Per la creazione di visualizzazioni di correlazione tra campioni.
* **`writexl`**: Per l'output finale dei risultati in formato `.xlsx`.

### Punti Chiave dell'Analisi (Estratto dello Script)

1.  **Pulizia Metadati**: I nomi delle colonne della matrice di conteggio sono puliti tramite espressioni regolari per creare il `colData` di DESeq2, definendo le colonne **`Condition`**, **`CellLine`** e **`Treatment`**.
2.  **Modello Statistico**: Viene utilizzato un modello di design **`~Condition`** per la costruzione dell'oggetto `DESeqDataSetFromMatrix`.
3.  **Visualizzazione Esplorativa**: Viene eseguita una trasformazione logaritmica (`rlog`) per la creazione di **Plot PCA** e **Heatmap di Correlazione**, verificando la qualità del campione.
4.  **Contrasti Differenziali**: Vengono eseguiti quattro confronti chiave:

    * `MCF7_TamR_Tam` vs `MCF7_TamR_EtOH` (Effetto Tamoxifene nella resistenza)
    * `MCF7_WT_Tam` vs `MCF7_WT_EtOH` (Effetto Tamoxifene nella sensibilità)
    * `MCF7_TamR_Tam` vs `MCF7_WT_Tam` (Resistenza vs Sensibilità, trattati)
    * `MCF7_TamR_EtOH` vs `MCF7_WT_EtOH` (Resistenza vs Sensibilità, controlli)

5.  **Filtro**: I risultati di espressione differenziale vengono filtrati in modo rigoroso per **`padj < 0.05`** e salvati.

