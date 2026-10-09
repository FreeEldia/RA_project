# RA_project
## Transcriptomics analyse van synoviumweefsel van patiënten met reumatoïde artritis

## Inhoud/structuur

In deze repository zijn het R-script, de figuren en het verslag van de transcriptomics-analyse van reumatoïde artritis (RA) opgenomen.

- [Script](Script/) – R-script voor de verwerking en analyse van RNA-seq-data.
- [Figuren](Figuren/) – Workflow, volcano plot, GO-grafiek en KEGG-visualisatie.
- [README.md](README.md) – Verslag van het onderzoek.

### Inhoudsopgave

- [Introductie](#introductie)
- [Methoden](#methoden)
- [Resultaten](#resultaten)
- [Conclusie](#conclusie)
- [Databeheer](#databeheer)
- [Referenties](#referenties)



## Introductie
Reumatoïde artritis (RA) is een chronische systemische auto-immuunziekte die voornamelijk de synoviale gewrichten aantast. De ziekte wordt gekenmerkt door ontsteking van het synovium, wat uiteindelijk kan leiden tot kraakbeenafbraak, boterosie en verlies van gewrichtsfunctie. Hoewel de exacte oorzaak van RA nog niet volledig bekend is, spelen genetische aanleg, omgevingsfactoren en ontregeling van het immuunsysteem een belangrijke rol bij het ontstaan van de ziekte (Gabriel, 2001). Een belangrijk kenmerk van RA is de aanwezigheid van autoantistoffen, waaronder anti-citrullinated protein antibodies (ACPA), die vaak al vóór het ontstaan van klinische symptomen aantoonbaar zijn (Majithia & Geraci, 2007).
Transcriptomics maakt het mogelijk om op grote schaal genexpressie te bestuderen en biedt daardoor inzicht in de moleculaire mechanismen die betrokken zijn bij ziekteprocessen. Door verschillen in genexpressie tussen patiënten en gezonde controles te analyseren, kunnen betrokken genen en biologische pathways worden geïdentificeerd (Wang et al., 2009). Eerdere studies hebben aangetoond dat immuunactivatie, B-celactiviteit, cytokinesignalering en ontstekingsprocessen een centrale rol spelen bij RA (McInnes & Schett, 2011).
In deze studie werd RNA-seq data afkomstig van synoviumbiopten van vier patiënten met vastgestelde RA en vier controlepersonen geanalyseerd. Het doel van het onderzoek was om genen en biologische processen te identificeren die differentieel tot expressie komen bij RA en om de betrokken pathways verder te onderzoeken met behulp van Gene Ontology (GO)- en KEGG-pathwayanalyses. 

## Methoden

### Onderzoeksmateriaal en databron

Voor dit onderzoek werd gebruikgemaakt van publieke RNA-sequencingdata van synoviumbiopten van patiënten met reumatoïde artritis (RA) en gezonde controlepersonen. De oorspronkelijke data zijn afkomstig uit de [NCBI Gene Expression Omnibus (GEO), dataset GSE89408](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89408), behorend bij [BioProject PRJNA352076](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA352076). De oorspronkelijke dataset is onder andere beschreven door [Guo et al. (2017)](https://doi.org/10.4049/jimmunol.1601988).

Voor de analyse werden acht vrouwelijke proefpersonen geselecteerd: vier patiënten met reumatoïde artritis (54–66 jaar) en vier gezonde controles (15–42 jaar). De RA-patiënten waren positief voor anti-citrullinated protein antibodies (ACPA). De sequencingdata waren afkomstig van het Illumina HiSeq 2000-platform.

De gebruikte SRA-runs waren:

| Onderzoeksgroep | SRA-runs | Aantal |
|---|---|---|
| RA | SRR4785979, SRR4785980, SRR4785986, SRR4785988 | 4 |
| Gezonde controles | SRR4785819, SRR4785820, SRR4785828, SRR4785831 | 4 |

### RNA-seq-verwerking en genkwantificatie

De bio-informatica-analyses werden uitgevoerd in **R (versie 4.5.2)**. Voor de verwerking van ruwe sequencingdata werden paired-end FASTQ-subsets gebruikt, aangeduid als `subset40k`.

Het humane referentiegenoom GRCh38 ([NCBI RefSeq-assembly GCF_000001405.26](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)) werd geïndexeerd met `buildindex()` uit **Rsubread (versie 2.24.0)** ([Liao et al., 2019](https://doi.org/10.1093/nar/gkz114)). Vervolgens werden de paired-end reads met `align()` tegen het referentiegenoom gemapt.

De verkregen BAM-bestanden werden gesorteerd en geïndexeerd met `sortBam()` en `indexBam()` uit **[Rsamtools](https://bioconductor.org/packages/Rsamtools/) (versie 2.26.0)**.

Voor de genkwantificatie werd `featureCounts()` uit Rsubread gebruikt. Hierbij werden de parameters `isPairedEnd = TRUE`, `isGTFAnnotationFile = TRUE`, `GTF.attrType = "gene_id"` en `useMetaFeatures = TRUE` toegepast.

Als annotatie werd het NCBI RefSeq-bestand `genomic.gtf` gebruikt. Uit de metadata van dit bestand bleek dat het afkomstig was van **GRCh38.p14**, assembly **GCF_000001405.40**, annotatierelease **RS_2025_08** van **1 augustus 2025**. Deze annotatie is gekoppeld aan de [NCBI GRCh38.p14-assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/).

De afzonderlijke verwerking van de FASTQ-subsets resulteerde in `RA_countmatrix.csv`. Voor de differentiële genexpressieanalyse en daaropvolgende functionele analyses werd de count matrix **`count_matrix_RA.txt`** gebruikt. Deze omvatte **29.407 genen en acht samples**. De definitieve DESeq2-, GO- en KEGG-resultaten zijn gebaseerd op deze TXT-count matrix en niet op de count matrix uit de FASTQ-subsets.

### Differentiële genexpressieanalyse

Differentiële genexpressie tussen RA-patiënten en gezonde controles werd onderzocht met **DESeq2 (versie 1.50.2)** ([Love et al., 2014](https://doi.org/10.1186/s13059-014-0550-8)).

Met `DESeqDataSetFromMatrix()` werd een dataset aangemaakt op basis van de count matrix en een groepstabel waarin iedere SRA-run aan de juiste onderzoeksgroep was gekoppeld. Hierbij werd het model `design = ~ treatment` gebruikt.

Vervolgens werd `DESeq()` uitgevoerd voor normalisatie, dispersieschatting en statistische toetsing. Met `results()` werden de log2-fold changes en aangepaste p-waarden bepaald. De vergelijking werd expliciet ingesteld als **RA ten opzichte van gezonde controles** met `contrast = c("treatment", "Reuma", "control")`.

Een positieve log2-fold change duidde op een hogere genexpressie bij RA en een negatieve log2-fold change op een lagere expressie bij RA.

De p-waarden werden met de Benjamini-Hochberg-methode gecorrigeerd voor multiple testing. Genen met een aangepaste p-waarde (`padj`) kleiner dan **0,05** werden als significant differentieel geëxpresseerd beschouwd. Voor de afzonderlijke aantallen op- en neergereguleerde genen werden daarnaast respectievelijk `log2FoldChange > 1` en `log2FoldChange < -1` gebruikt.

De resultaten werden gevisualiseerd in een volcano plot met **[EnhancedVolcano](https://bioconductor.org/packages/EnhancedVolcano/)**.


### Gene Ontology-verrijkingsanalyse

Om te onderzoeken welke biologische functies en processen oververtegenwoordigd waren onder de differentieel geëxpresseerde genen, werd een Gene Ontology (GO)-verrijkingsanalyse uitgevoerd met **goseq (versie 1.62.0)** ([Young et al., 2010](https://doi.org/10.1186/gb-2010-11-2-r14)).

Eerst werden alle genen geselecteerd waarvoor een aangepaste DESeq2-p-waarde beschikbaar was. Vervolgens werd een binaire genvector opgesteld waarin genen met `padj < 0,05` de waarde 1 kregen en de overige geteste genen de waarde 0.

Met `nullp()` werd een *probability weighting function* berekend om te corrigeren voor mogelijke selectiebias. Hierbij werd de som van de read counts per gen als biasvariabele (`bias.data`) gebruikt.

De GO-overrepresentatieanalyse werd uitgevoerd met `goseq()` en de annotatie-instellingen `genome = "hg19"` en `id = "geneSymbol"`.

De p-waarden voor overrepresentatie (`over_represented_pvalue`) werden met `p.adjust(method = "BH")` gecorrigeerd voor multiple testing. GO-termen met een aangepaste p-waarde kleiner dan **0,05** werden als significant beschouwd.

De tien hoogst gerangschikte significante GO-termen werden met **[dplyr](https://dplyr.tidyverse.org/)** geselecteerd en met **[ggplot2](https://ggplot2.tidyverse.org/)** gevisualiseerd. In de dotplot werd het percentage differentieel geëxpresseerde genen binnen de GO-term weergegeven op de horizontale as. De puntgrootte vertegenwoordigde het aantal differentieel geëxpresseerde genen en de kleur de aangepaste p-waarde.

### KEGG-pathwayvisualisatie

Voor de verdere interpretatie van immuungerelateerde genexpressieverschillen werd de humane **B-cell receptor signaling pathway** geselecteerd ([KEGG: hsa04662](https://www.kegg.jp/pathway/hsa04662)).

De gensymbolen uit de DESeq2-resultaten werden met `mapIds()` uit **[AnnotationDbi](https://bioconductor.org/packages/AnnotationDbi/)** en de humane annotatiedatabase **[org.Hs.eg.db](https://bioconductor.org/packages/org.Hs.eg.db/)** gekoppeld aan Entrez Gene-ID's. Bij meerdere mogelijke koppelingen werd met `multiVals = "first"` de eerste gevonden Entrez-ID gebruikt.

Genen zonder geldige Entrez-ID of met ontbrekende of niet-eindige log2-fold changes werden uitgesloten van de visualisatie.

Met **Pathview (versie 1.50.0)** ([Luo & Brouwer, 2013](https://doi.org/10.1093/bioinformatics/btt285)) werden de log2-fold changes geprojecteerd op de geselecteerde KEGG-pathway. Hierbij werden `pathway.id = "hsa04662"`, `species = "hsa"`, `gene.idtype = "ENTREZID"` en `limit = list(gene = 5)` gebruikt.

De Pathview-figuur werd gebruikt om veranderingen in genexpressie binnen de geselecteerde pathway te visualiseren. Er werd geen afzonderlijke statistische KEGG-verrijkingsanalyse uitgevoerd.

### Reproduceerbaarheid 
Het R-script met de uitgevoerde analyses is beschikbaar in de [Script-map van de GitHub-repository](https://github.com/FreeEldia/RA_project/tree/main/Script).

## Resultaten

De differentiële expressieanalyse liet verschillen in genexpressie zien tussen RA-patiënten en controles. In totaal waren 5119 genen significant differentieel geëxpresseerd (padj < 0,05). Daarvan hadden 2487 genen een log2FC > 1 en 2084 genen een log2FC < -1. Hiervan waren **2487 genen opgereguleerd** en **2084 genen neergereguleerd**. De verdeling van de differentiële genexpressie is weergegeven in **Figuur 2**. Onder de opvallende genen bevonden zich onder andere *BCL2A1*, *ADAMDEC1* en meerdere immunoglobuline-gerelateerde genen, waaronder *IGHV3-53*, *IGHV1-69*, *IGHG4*, *IGHV4-31* en *IGKV2-28*.

![Figuur 2](Figuren/VolcanoplotRA.png)


<sub>**Figuur 2.** Volcano plot van genexpressieverschillen tussen RA-patiënten en controles. De x-as geeft de log2 fold change weer en de y-as de −log10 p-waarde. Grijze punten voldoen niet aan de significantiecriteria, groene punten voldoen aan het log2-fold-changecriterium en rode punten voldoen zowel aan het p-waarde- als het log2-fold-changecriterium.. </sub>

De Gene Ontology-analyse liet zien dat verschillende immuungerelateerde processen significant verrijkt waren (Figuur 3). De meest verrijkte termen waren onder andere immunoglobulin complex, adaptive immune response, leukocyte activation, immune response en immune system process. Vooral de aanwezigheid van *immunoglobulin complex* en *adaptive immune response* wijst erop dat een relatief groot aantal differentieel geëxpresseerde genen betrokken is bij de adaptieve immuunrespons.

![Figuur 3](Figuren/GO.png)

<sub>**Figuur 3**. Top 10 verrijkte GO-termen van de differentieel geëxpresseerde genen. De grootte van de punten geeft het aantal differentieel geëxpresseerde genen per GO-term weer. De kleur geeft de statistische significantie van de verrijking weer.</sub>


Op basis van de verrijking van immuungerelateerde GO-termen werd de **B-cell receptor signaling pathway (hsa04662)** onderzocht. In deze KEGG-pathway werden zowel positieve als negatieve veranderingen in genexpressie waargenomen (**Figuur 4**). Verschillende componenten van de B-cel-signaleringsroute vertoonden veranderde expressie. In combinatie met de immunoglobuline-gerelateerde genen en GO-termen laat dit zien dat genen betrokken bij B-celfunctie verschillen in expressie tussen het RA- en controleweefsel.

![Figuur 4](Figuren/hsa04662.pathview.png)

<sub>**Figuur 4**. KEGG B-cell receptor signaling pathway (hsa04662) met genexpressieveranderingen van RA ten opzichte van controles. Rode vakken geven een positieve log2 fold change aan en groene vakken een negatieve log2 fold change.</sub>

## Conclusie

In deze studie werd transcriptomics gebruikt om verschillen in genexpressie tussen synoviumweefsel van patiënten met reumatoïde artritis en controles te onderzoeken. De differentiële expressieanalyse identificeerde een groot aantal genen waarvan de expressie significant verschilde tussen beide groepen. De Gene Ontology-analyse liet zien dat vooral immuungerelateerde processen, waaronder adaptieve immuunrespons, immuunactivatie en immunoglobuline-gerelateerde functies, sterk vertegenwoordigd waren.
De daaropvolgende KEGG-analyse van de B-cell receptor signaling pathway toonde aan dat meerdere genen binnen deze route veranderingen in expressie vertoonden. Deze bevinding sluit aan bij de bekende rol van B-cellen en autoantistofproductie in de pathogenese van reumatoïde artritis. De aanwezigheid van verrijkte immunologische processen en veranderingen binnen B-celgerelateerde signaleringsroutes ondersteunt het belang van adaptieve immuunmechanismen bij deze ziekte.
Een beperking van deze studie is het relatief kleine aantal monsters en het gebruik van subsets van de oorspronkelijke sequencingdata. Toekomstig onderzoek zou gebruik kunnen maken van grotere datasets en aanvullende pathwayanalyses om de betrokken moleculaire mechanismen verder te karakteriseren. Desondanks tonen de resultaten aan dat transcriptomics een waardevolle methode is om biologische processen en genen te identificeren die betrokken zijn bij de ontwikkeling van reumatoïde artritis.
## Databeheer
GitHub werd gebruikt voor versiebeheer en documentatie van de transcriptomics-analyse. Het gebruikte R-script is in de repository opgeslagen en per analysestap voorzien van commentaar, zodat zichtbaar is hoe de ruwe RNA-seq-data zijn verwerkt tot de uiteindelijke resultaten. Figuren en overige outputbestanden zijn in afzonderlijke mappen opgeslagen. Wijzigingen aan bestanden werden met Git-commits vastgelegd, waardoor eerdere versies behouden blijven en aanpassingen aan de analyse traceerbaar zijn. Grote ruwe sequencingbestanden zijn vanwege hun bestandsgrootte niet in de repository opgenomen; de gebruikte dataset en referentiebestanden zijn daarom beschreven met hun oorspronkelijke bron en accessienummers. Deze structuur maakt de analyse transparanter en maakt het mogelijk om de uitgevoerde stappen met dezelfde inputbestanden en software opnieuw uit te voeren

## Referenties
## AI disclaimer 
Voor het maken van dit verslag is AI gebruikt voor het controleren van spelling, grammatica en ondersteuning bij programmeervragen


### Databronnen

- NCBI Gene Expression Omnibus. *RNAseq study of synovial biopsies*. GEO Series GSE89408. https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89408

- Guo, Y., et al. (2017). CD40L-Dependent Pathway Is Active at Various Stages of Rheumatoid Arthritis Disease Progression. *The Journal of Immunology, 198*(11), 4490–4501. https://doi.org/10.4049/jimmunol.1601988

- NCBI. *Homo sapiens genome assembly GRCh38.p14*, GCF_000001405.40. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/

- NCBI. *Homo sapiens reference genome GRCh38*, GCF_000001405.26. https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

### R en softwarepackages

- R Core Team. *R: A Language and Environment for Statistical Computing*. R Foundation for Statistical Computing. https://www.r-project.org/

- Liao, Y., Smyth, G. K., & Shi, W. (2019). The R package Rsubread is easier, faster, cheaper and better for alignment and quantification of RNA sequencing reads. *Nucleic Acids Research, 47*(8), e47. https://doi.org/10.1093/nar/gkz114

- Rsamtools. *Bioconductor package*. https://doi.org/10.18129/B9.bioc.Rsamtools

- Love, M. I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology, 15*, 550. https://doi.org/10.1186/s13059-014-0550-8

- EnhancedVolcano. *Bioconductor package*. https://doi.org/10.18129/B9.bioc.EnhancedVolcano

- Young, M. D., Wakefield, M. J., Smyth, G. K., & Oshlack, A. (2010). Gene ontology analysis for RNA-seq: accounting for selection bias. *Genome Biology, 11*, R14. https://doi.org/10.1186/gb-2010-11-2-r14

- goseq. *Bioconductor package*. https://doi.org/10.18129/B9.bioc.goseq

- Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis*. Springer. https://ggplot2.tidyverse.org/

- dplyr. *R package for data manipulation*. https://dplyr.tidyverse.org/

- AnnotationDbi. *Bioconductor package*. https://doi.org/10.18129/B9.bioc.AnnotationDbi

- org.Hs.eg.db. *Bioconductor human gene annotation database*. https://doi.org/10.18129/B9.bioc.org.Hs.eg.db

- Luo, W., & Brouwer, C. (2013). Pathview: An R/Bioconductor package for pathway-based data integration and visualization. *Bioinformatics, 29*(14), 1830–1831. https://doi.org/10.1093/bioinformatics/btt285