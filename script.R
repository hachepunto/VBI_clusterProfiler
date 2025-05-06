#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(tidyverse)

#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)

# Lectura de la tabla de genes diferencialemente expresados
degs = readRDS("data/degs.RDS")

# necesitamos el log2 fold change 
original_gene_list <- degs$logFC

# Nombramos el vector
names(original_gene_list) <- degs$ESGN

# eliminamos cualquier NA 
gene_list<-na.omit(original_gene_list)

# odernamos la lista en orden decreciente (requerido por clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# extraemos los genes significativos (p ajustada < 0.05)
sig_genes_df = subset(degs, adj.P.Val < 0.05)

# para los resultados significativos, queremos filtrar por log2fold change
genes <- sig_genes_df$logFC

# nombramos el vector
names(genes) <- sig_genes_df$ESGN

# omitimos posibles NAs
genes <- na.omit(genes)

# filtramos por mínimo log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]

go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

head(go_enrich)

#BiocManager::install("enrichplot")
library(enrichplot)
upsetplot(go_enrich)

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

dotplot(go_enrich)

# Convertimos los gene IDs para la función enrichKEGG
# Podría ser que perdamos algunos genes aquí dado que algunas conversiones no son compatibles
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# Removemos IDS duplicados (aquí se usa "ENSEMBL", pero debería de ser lo que hayamos usado como keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Creamos un nuevo dataframe df2 el cual solo tiene genes que se an convertido exitosamente usando la función bitr

degs2 = degs[degs$ESGN %in% dedup_ids$ENSEMBL,]

# Creamos una nueva columna en degs2 con los ENTREZ ids correspondientes
degs2$Y = dedup_ids$ENTREZID

# Creamos un vector de genes del universo
kegg_gene_list <- degs2$logFC

# Nombramos el vector con los ENTREZ ids
names(kegg_gene_list) <- degs2$Y

# Verificamos que no haya NAs 
kegg_gene_list<-na.omit(kegg_gene_list)

# Ordenamos la lista en orden decreciente (requerido por for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# Extraemos los resultados significativos del degs2
kegg_sig_genes_degs = subset(degs2, adj.P.Val < 0.05)

# También filtraremos por log2fold change
kegg_genes <- kegg_sig_genes_degs$logFC

# Nombramos el vector con los IDs convertidos
names(kegg_genes) <- kegg_sig_genes_degs$Y

# eliminanmos NAs
kegg_genes <- na.omit(kegg_genes)

# Y ahora si filtramos por log2fold change
kegg_genes <- names(kegg_genes)[abs(kegg_genes) > 2]

# Creamos el objeto enrichKEGG

kk <- enrichKEGG(gene=kegg_genes, 
                universe=names(kegg_gene_list),
                organism="hsa",
                pvalueCutoff = 0.05, 
                keyType = "ncbi-geneid")
head(kk)

# Barplot

barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

# Dotplot

dotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)



# Pathview

#BiocManager::install("pathview")
library(pathview)

# Produce una gráfica de KEGG (PNG)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3])

# Produce una gráfica diferente (PDF)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3], kegg.native = FALSE)



# GO Gene Set Enrichment Analysis con clusterProfiler

#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)


#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)

# Lectura de la tabla de genes diferencialemente expresados

degs = readRDS("data/degs.RDS")

# necesitamos el log2 fold change 
original_gene_list <- degs$logFC

# Nombramos el vector
names(original_gene_list) <- degs$ESGN

# eliminamos cualquier NA 
gene_list<-na.omit(original_gene_list)

# odernamos la lista en orden decreciente (requerido por clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH")

head(gse)


# Dotplot

#BiocManager::install("DOSE")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)


# Ridgeplot

ridgeplot(gse) + labs(x = "enrichment distribution")

# GSEA Plot

# Usamos el objeto 'geneset' para que siempre coincidan el título y el gene set correspondiente.
geneset = 1
gseaplot(gse, by = "all", title = gse$Description[geneset], geneSetID = geneset)


# KEGG Gene Set Enrichment Analysis con clusterProfiler

# Convertir genes IDs para la función gseKEGG
# Podría ser que se perdieran algunos genes por la conversión
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# se elimnan ids duplicados (aquí usamos "ENSEMBL", pero debemos usar lo que hayamos empleado en el argumento "fromType")
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Creamos un nuevo dataframe degs2 en cual solo contiene los genes que hicieron match usando ña función bitr
degs2 = degs[degs$ESGN %in% dedup_ids$ENSEMBL,]

# Creamos una nueva coñumna en degs2 con los correspondientes ENTREZ IDs
degs2$Y = dedup_ids$ENTREZID

# Creamos un vector con el universo de genes
kegg_gene_list <- degs2$logFC

# Nombramos el vector con los ENTREZ ids
names(kegg_gene_list) <- degs2$Y

# eliminamos NAs 
kegg_gene_list<-na.omit(kegg_gene_list)

# Ordenamos los datos en orden decreciente (requerido por clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

head(kk2, 10)

# Dotplot

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

# Ridgeplot

ridgeplot(kk2) + labs(x = "enrichment distribution")

# GSEA Plot

# Usamos el objeto 'geneset' para que siempre coincidan el título y el gene set correspondiente.
geneset = 1
gseaplot(kk2, by = "all", title = gse$Description[geneset], geneSetID = geneset)


# Pathview

#BiocManager::install("pathview")
library(pathview)

# Produce una gráfica de KEGG (PNG) Usaremos hsa05322 proque el primer pathway más enriquecido fue el mismo que graficamos arriba
hsa <- pathview(gene.data=gene_list, pathway.id="hsa05322", species = "hsa", gene.idtype=gene.idtype.list[3])

# Produce una gráfica diferente (PDF)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa05322", species = "hsa", gene.idtype=gene.idtype.list[3], kegg.native = FALSE)

