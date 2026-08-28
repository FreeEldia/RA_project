## Methoden

Voor deze studie werd publieke **paired-end RNA-sequencingdata** van synoviumbiopten gebruikt uit Platzer et al. (2019): vier vrouwelijke patiënten met gevestigde reumatoïde artritis (54–66 jaar) en vier vrouwelijke controles (15–42 jaar). De RA-patiënten waren ACPA-positief.
![Figuur 1](Figuren/Figuur_workflow.png)

<sub>**Figuur 1.** Overzicht van de uitgevoerde transcriptomics-analyse.</sub>

De analyse werd uitgevoerd in **R (versie [4.5.2])** volgens het stroomschema in **Figuur 1**. Reads werden met **Rsubread [2.24.0]** gemapt tegen het humane referentiegenoom **GRCh38 (GCF_000001405.26)**. De resulterende BAM-bestanden werden gesorteerd en geïndexeerd met **Rsamtools [2.26.0]**.

Met `featureCounts()` werd een count matrix gegenereerd met het bijbehorende NCBI GTF-annotatiebestand. Differentiële genexpressie tussen RA en controles werd bepaald met **DESeq2 [1.50.2]**. Genen met een aangepaste p-waarde (*padj*) < 0,05 werden als significant beschouwd.

Met **goseq [1.62.0]** werd een Gene Ontology (GO)-analyse uitgevoerd. Op basis van de GO-resultaten werd de **B-cell receptor signaling pathway (hsa04662)** geselecteerd. Log2-fold changes werden met **Pathview [1.50.0]** gekoppeld aan Entrez Gene-ID's en gevisualiseerd op de humane KEGG-pathway.
