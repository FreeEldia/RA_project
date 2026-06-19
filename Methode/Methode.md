## Methoden

Voor deze studie werd gebruikgemaakt van publieke paired-end RNA-sequencingdata afkomstig van synoviumbiopten van vier vrouwelijke patiënten met gevestigde reumatoïde artritis (54–66 jaar) en vier vrouwelijke controlepersonen zonder RA (15–42 jaar).
![Figuur 1](Figuren/Figuur_workflow.png)

<sub>**Figuur 1.** Overzicht van de uitgevoerde transcriptomics-analyse.</sub>

Zoals weergegeven in **Figuur 1** werden de analyses uitgevoerd in R (versie 4.5.2). Eerst werd het humane referentiegenoom GRCh38 (GCF_000001405.26) geïndexeerd met behulp van het package **Rsubread**. Vervolgens werden paired-end RNA-seq reads gemapt tegen het referentiegenoom met de functie `align()`. De resulterende BAM-bestanden werden gesorteerd en geïndexeerd met **Rsamtools**.

Voor het bepalen van genexpressie werd met `featureCounts()` een count matrix gegenereerd op basis van een GTF-annotatiebestand afkomstig van dezelfde referentieversie als het gebruikte genoom. De count matrix bevatte gen-specifieke readaantallen voor alle acht samples.

Differentiële genexpressie tussen de RA-groep en controlegroep werd bepaald met het package **DESeq2**. Hierbij werden de ruwe counts genormaliseerd en statistisch vergeleken. Genen met een aangepaste p-waarde (*padj*) kleiner dan 0,05 werden als significant beschouwd.

Vervolgens werd met het package **goseq** een Gene Ontology (GO)-analyse uitgevoerd om verrijkte biologische processen te identificeren. De resultaten werden gevisualiseerd in een volcano plot en een GO-verrijkingsplot. Op basis van de GO-resultaten werd een KEGG-pathway geselecteerd voor verdere analyse. Hiervoor werd met het package **Pathview** de **B-cell receptor signaling pathway (hsa04662)** gevisualiseerd.
