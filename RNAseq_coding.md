# RNAseq analysis  

## to get exonic region bed file from hg38   
`wget -O - "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz" | gunzip -c | grep 'transcript_type "protein_coding"' | awk '($3=="exon") {printf("%s\t%s\t%s\n",$1,int($4)-1,$5);}' |sort -T . -k1,1 -k2,2n | bedtools merge > hg38_exonic_chr.bed`  

## paths  
`GENOME=/path/to/Homo_sapiens_assembly38.fasta`  
`ANOTATION=/path/to/hg38.exonic.chr.gtf`  
`GENOME_Index=/path/to/GENOME_Index/`  

## QC    

`for i in *.fastq.gz; do fastqc --outdir QC_raw $i; done`  
`multiqc .`  

## Trimming if necessary  
`trim_galore --paired --quality 20 --phred33 --clip_R1 3 --clip_R2 3 --stringency 5 --trim-n --output_dir cleaned R1.fastq.gz R2.fastq.gz`  

## adding chr (optional)
`awk 'OFS="\t" {if (NR > 5) $1="chr"$1; print}' Homo_sapiens.GRCh38.103.gtf >hg38.103.gtf`  

## Creating all the required index  

`STAR --runThreadN 8 --runMode genomeGenerate --genomeFastaFiles $GENOME --sjdbGTFfile $ANOTATION --sjdbOverhang 99 --genomeDir $GENOME_Index`  

## Mapping  

`cat list_files.txt | while read a b c;`  
`	do`  
`	STAR --genomeDir $GENOME_Index --readFilesIn $a $b --readFilesCommand zcat --outFileNamePrefix $c --outSAMtype BAM SortedByCoordinate --runThreadN 8 --outSAMunmapped Within --outSAMattributes Standard`  
`done`  

## index bam files    

`for i in *.bam; do samtools index -@8 $i; done`  

## convertion to sam (Optional)
`for i in *.bam; do samtools view -@8 -h *.bam > cleaned_${i}.sam; done`  

## keep only mapping qulity >=20  

`for i in *.bam;`  
`    do`  
`    samtools view --threads 8 -q 20 -h -b $i >cleaned/cleaned_${i};`  
`done`  


## bam QC    

`cat list_bams.txt | while read a b;`  
`	do`  
`	samtools flagstat -@8 $a > alignment_qc/$b`  
`done`  


`cat list_bams.txt | while read a b;`  
`	do`  
`	bam_stat.py -i $a > alignment_qc/$b`  
`done`  

`cat list_bams.txt | while read a b;`  
`	do`  
`	geneBody_coverage.py -i $a -r $ANOTATION -o alignment_qc/$b &`  
`done`  

`cat list_bams.txt | while read a b;`  
`	do`
`	read_distribution.py -r $ANOTATION -i $a > alignment_qc/$b`  
`done`  

`for i in *_alignment_STAR_Aligned.sortedByCoord.out.bam;`  
`    do`  
`    echo -e "$i \n"`  
`    echo -e "count of unmapped reads \n"`  
`    samtools view -f 4 $i | wc -l`  
`    echo -e "\ncount of mapped reads \n"`  
`    samtools view -F 4 $i | wc -l`  
`    echo -e "\n mapping qulity \n"`  
`    samtools view -q 20 $i | wc -l`  
`    echo -e "\n Total counts of reads\n"`  
`    samtools view $i | wc -l`  
`    echo -e "\n Done!! Next file\n"`  
`done`  

`multiqc .`  

## counting    

`featureCounts -p -T 8 -g gene_name -a $ANOTATION -o raw_featCounts_genes.txt *bam`  

## R DEseq   

## import packages  

`library(clusterProfiler)`  
`library(biomaRt)`  
`library(ReactomePA)`  
`library(DOSE)`  
`library(KEGG.db)`  
`library(genefilter)`  
`library(GO.db)`  
`library(topGO)`  
`library(gage)`  
`library(ggsci)`  
`library(dplyr)`  
`library(DESeq2)`  
`library(ggplot2)`  
`library(pheatmap)`  
`library(RColorBrewer)`  
`library(org.Hs.eg.db)`  

## import raw data (raw counts)  
`countdata <- read.table("raw_featCounts_genes.txt", header = TRUE, skip = 1, row.names = 1)`  

## rename header
`colnames(countdata) <- gsub(".bam", "", colnames(countdata), fixed = T)`  
`head(countdata)`  

## remove length column  
`countdata <- countdata[ ,c(-1:-5)]`  

## meta file in tab sep txt file  
`metadata <- read.delim("metadata.txt", row.names = 1)`  
`metadata$sampleid <- row.names(metadata)`  
`metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]`  

`ddsMat <- DESeqDataSetFromMatrix(countData = countdata, colData = metadata, design = ~Group)`  

`ddsMat <- DESeq(ddsMat)`  

`results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)`  

`summary(results)`  

`mcols(results, use.names = T)`  

## enrich data
`results$description <- mapIds(x = org.Hs.eg.db,keys = row.names(results), column = "GENENAME", keytype = "SYMBOL", multiVals = "first")`  

`results$symbol <- row.names(results)`  

`results$entrez <- mapIds(x = org.Hs.eg.db, keys = row.names(results), column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")`  

`results$ensembl <- mapIds(x = org.Hs.eg.db, keys = row.names(results), column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")`  

`results_sig <- subset(results, padj < 0.05)`  

## writing outputs into files

`write.table(x = as.data.frame(counts(ddsMat), normalized = T), file = 'cleaned_normalized_counts.txt', sep = '\t', quote = F, col.names = NA)`  
`write.table(x = counts(ddsMat[row.names(results_sig)], normalized = T), file = 'cleaned_normalized_counts_significant.txt', sep = '\t', quote = F, col.names = NA)`  
`write.table(x = as.data.frame(results), file = "cleaned_results_gene_annotated.txt", sep = '\t', quote = F, col.names = NA)`  
`write.table(x = as.data.frame(results_sig), file = "cleaned_results_gene_annotated_significant.txt", sep = '\t', quote = F, col.names = NA)`  


## PCA Plot

`ddsMat_rlog <- rlog(ddsMat, blind = FALSE)`  

`tiff("cleaned_PCA_top_500.tiff", units="in", width=5, height=5, res=300)`  
`plotPCA(ddsMat_rlog, intgroup = "Group", ntop = 500) + `  
`  theme_bw() + `  
`  geom_point(size = 5) + `  
`  scale_y_continuous(limits = c(-5, 5)) + `  
`  ggtitle(label = "PCA (multiple comparison)", `  
`          subtitle = "Top 500 most variable genes") `  
`dev.off()`  

## heatmap  
`ddsMat_rlog <- rlog(ddsMat, blind = FALSE)`  

## select most seg 50 genes    
`mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:50, ]`  

## add col data  
`annotation_col = data.frame(`  
`  Group = factor(colData(ddsMat_rlog)$Group), `  
`  Replicate = factor(colData(ddsMat_rlog)$Replicate),`  
`  row.names = colData(ddsMat_rlog)$sampleid`  
`)`  

## colours  
`ann_colors = list(`  
`  Group = c(Control = "red", IFRNS = "darkorange", FRNS_SDNS = "lightblue", SRNS = "green"),`  
`  Replicate = c(Rep1 = "darkred", Rep2 = "forestgreen")`  
`)`  

## plot  
`tiff("cleaned_heatmap.tiff", units="in", width=10, height=10, res=300)`  
`pheatmap(mat = mat, `  
`         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255),`  
`         scale = "row",`  
`         annotation_col = annotation_col,`  
`         annotation_colors = ann_colors,`  
`         fontsize = 6.5, `  
`         cellwidth = 55,`  
`         show_colnames = T)`  
`dev.off()`  

`tiff("cleaned_PD_plot.tiff", units="in", width=5, height=5, res=300)`  
`plotDispEsts(ddsMat)`  
`dev.off()`  

## Single gene plot
`ddsMat_rlog <- rlog(ddsMat, blind = FALSE)`  
`top_gene <- rownames(results)[which.min(results$log2FoldChange)]`  
`tiff("cleaned_top_single_plot.tiff", units="in", width=5, height=5, res=300)`  
`plotCounts(dds = ddsMat, `  
`           gene = top_gene, `  
`           intgroup = "Group", `  
`           normalized = T, `  
`           transform = T)`  
`dev.off()`  

## pathway
`results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)`  
`gene_matrix <- results_sig_entrez$log2FoldChange`  
`names(gene_matrix) <- results_sig_entrez$entrez`  

## enrich with kigg data
`kegg_enrich <- enrichKEGG(gene = names(gene_matrix),`  
`                          organism = 'human',`  
`                          pvalueCutoff = 0.05,`  
`                          qvalueCutoff = 0.10)`  

`tiff("cleaned_kegg_enrich_bar_plot.tiff", units="in", width=5, height=5, res=300)`  
`barplot(kegg_enrich, `  
`        drop = TRUE, `  
`        showCategory = 10,`  
`        title = "KEGG Enrichment Pathways",`  
`        font.size = 8)`  
`dev.off()`  

`go_enrich <- enrichGO(gene = names(gene_matrix),`  
`                      OrgDb = 'org.Hs.eg.db',`  
`                      readable = T,`  
`                      ont = "BP",`  
`                      pvalueCutoff = 0.05,`  
`                      qvalueCutoff = 0.10)`  

`tiff("cleaned_GO_enrich_bar_plot.tiff", units="in", width=5, height=5, res=300)`  
`barplot(go_enrich, `  
`        drop = TRUE, `  
`        showCategory = 10,`  
`        title = "GO Biological Pathways",`  
`        font.size = 8)`  
`dev.off()`  
