# Introducción a `clusterProfiler`

## Enriquecimientos 

`Un análisis de enriquecimiento funcional suele aplicarse para extraer información biológica relevante, ya sea de una lista filtrada de genes relevantes o bien de todo el genoma evaluado por un valor estadístico.` El ejemplo más conocido de esta lista podría ser la lista filtrada de genes diferencialmente expresados (por algún punto de corte de logFold Change) o bien la lista de todos los genes estimado por su valor de logFold Change.

Si bien hay varios métodos para hacer enriquecimiento[^1], los dós métodos más comunes son el Análisis de sobrerepresentación (*Over-representation analysis* o ORA) y el Análisis de Enriquecimiento de Conjuntos de Genes (*Gene Set Enrichment Analysis* o GSEA). Ambos incluidos en `clusterProfiler`.

[^1]: Para más detalles ver: [Gómez-Romero, *et al.* (**2021**) "Bioinformatics of Functional Categories Enrichment" in *Bioinformatics and Human Genomics Research* eBook ISBN: 9781003005926](https://www.taylorfrancis.com/chapters/edit/10.1201/9781003005926-14/bioinformatics-functional-categories-enrichment-laura-g%C3%B3mez-romero-hugo-tovar-enrique-hern%C3%A1ndez-lemus).

| Análisis | Lo que require de entrada | Como luce la salida |  ✅ Pros | ⚠️ Cons |
| -------- | :-----------------------: | :-----------------: | :-----: | :-----: |
| ORA (Over-representation Analysis) | Una lista de genes (no se requeiren estadísticos) | Una prueba hipergeométrica por vía | - Sencillo <br> - Barato computacionalmente para calcular *p*-values | - Ignora estadísticas asociadas con los genes y requiere puntos de corte arbitrarios <br> - Asume independencia de genes y pathways |
| GSEA (Gene Set Enrichment Analysis) | Una lista de IDs de genes con un valor estadístico por gen | Un puntaje de enriquecimiento (*enrichment score* ES) por vía | - Incluye todos los genes (Sin puntos de corte arbitrarios) <br> |  - Las permutaciones pueden ser computacionalmente costosas <br> - No considera sobrelape de vías |



## *Over-representation analysis*


![ORA](https://hbctraining.github.io/DGE_workshop_salmon_online/img/go_proportions.png "Over-representation analysis")

## *Gene Set Enrichment Analysis*

[Video sobre GSEA](https://youtu.be/bT00oJh2x_4)

![GSEA](images/gsea.png "Gene Set Enrichment Analysis")

### Enrichr

Herramienta en linea para hacer ORA
https://maayanlab.cloud/Enrichr/


### Entonces ¿Para qué queremos clusterProfiler?

# Over-representation analysis

## GO *Over-representation analysis* con `clusterProfiler`

### Instalar y cargar paquetes
```r
#BiocManager::install("clusterProfiler")
library(clusterProfiler)
library(tidyverse)
```

### Anotación
Vamos a usar los genes diferencialmente expresados del trabajo [*The regulatory landscape of retinoblastoma: a pathway analysis perspective*](https://doi.org/10.6084/m9.figshare.c.5975228), por lo que usaremos la anotación de humano "org.Hs.eg.db". Puedes ver las anotaciones disponibles [aquí](http://bioconductor.org/packages/release/BiocViews.html#___OrgDb).

```r
#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
```

### Preparar entrada

```r
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
```

### Crear objeto enrichGO

Parámetros:
  
**ont** Ontología de GO a utilizar. Opciones: "BP", "MF" o "CC"  
**keyType** Esta puede variar depepndiendo de la anotación (gene ids). Por ejemplo para *"org.Hs.eg.db"*, las opciones son:   
  
<p style='text-align: justify;'> "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS" "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL" 
"GENENAME" "GENETYPE" "GO" "GOALL"        "IPI"          "MAP"          "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"
"PATH"         "PFAM"         "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL" "UCSCKG"       "UNIPROT" </p>
  
Puedes checar las opciones de la anotación de tu organismo usando la función `keytypes`, por ejemplo `keytypes(org.Hs.eg.db)`. 

#### Creación del objeto `enrichResult`

```r
go_enrich <- enrichGO(gene = genes,
                      universe = names(gene_list),
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
```
### Output

#### Table of results
```r
head(go_enrich)
```

### Upset Plot

Enfatiza el sobrelape de genes a lo largo de los diferentes set de genes.

```r
#BiocManager::install("enrichplot")
library(enrichplot)
upsetplot(go_enrich)
```

### Barplot

```r
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)
```

### Dotplot
```r
dotplot(go_enrich)
```

## KEGG *Over-representation analysis*

Si queremos enriquecer patways de **KEGG** necesitamos convertir los ids para poder usar la función `gseKEGG()`. Para ello podemos usar la función `bitr` que está incluida el `clusterProfiler`. Es normal que algunos mensajes y warnings salgan al usar esta función. 

En la función `bitr`, el parámetro `fromType` debe ser el tipo de `keyType` del de la función `gseGO` que vimos anteriormente (la fuente original). Este parámetro es usado dos veces más la crear los `dedup_ids` y el `degs2`.  

El otro parámetro `toType` en la función `bitr` también debe ser otra de las anitaciones disponibles el `keyTypes("org.Hs.eg.db")` y debe mapear a alguna de 'kegg', 'ncbi-geneid', 'ncbi-proteinid' o 'uniprot' porque la función `gseKEGG()` solo acepta alguna de estas 4 opciones como entrara en el parámetro `keytype`. En el caso de org.Hs.eg.db usaremos 'ENTREZID' para `toType`, ya que este se corresponde con el ncbi-geneid y tiene mejor match con 'ENSEMBL', pero podríamos también usar 'UNIPROT'. 

Como nuestra entrada inicial usaremos nuestra `original_gene_list` que creamos para el enriquecimiento en GO.

### Prepare Data

```r
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

```
### Creamos el objeto `enrichKEGG`

Parámetros:

**organism** código KEGG del organismo: La lista completa está aquí: https://www.genome.jp/kegg/catalog/org_list.html (se necesita el código de 3 letras). En el caso del humano es "hsa"  
**keyType** one of 'kegg', 'ncbi-geneid', 'ncib-proteinid' or 'uniprot'.  

```r
kk <- enrichKEGG(gene=kegg_genes, 
                universe=names(kegg_gene_list),
                organism="hsa",
                pvalueCutoff = 0.05, 
                keyType = "ncbi-geneid")
head(kk)
```

### Barplot
```r
barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)
```

### Dotplot
```r
dotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)
```

## Pathview

Con este paquete se puede visualizar nuestros genes en las vías de KEGG generando una imagen PNG y un PDF *diferente*.  
  
Parámetros:

**gene.data** Esta es la lista `gene_list` creada arriba, el tipo de ids tiene que coincidir con el de `gene.idtype`
**pathway.id** Tenemos que elegir aquí nosotros alguna. Como la idea es visualizar nuestras vías significativamente enriquecidas, pordemos usar los ids de los pathways que nos sale en `head(kk)`.  
**species** El id de `organism` de la función `enrichKEGG`, en este caso "hsa"
**gene.idtype** el tipo de ids utilizados en `gene.data` pero sacado del objeto `gene.idtype.list`. En este caso tomamos el tercer elemento de esta lista.

```r
#BiocManager::install("pathview")
library(pathview)

# Produce una gráfica de KEGG (PNG)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3])

# Produce una gráfica diferente (PDF)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3], kegg.native = FALSE)
```

Las imágenes se salvan en su directorio de trabajo.


# Gene Set Enrichment Analysis

## GO *Gene Set Enrichment Analysis* con `clusterProfiler`

### Instalación y cargado de paquetes

```r
#BiocManager::install("clusterProfiler")
#BiocManager::install("pathview")
#BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
```

### Anotación
Vamos a usar, de nuevo, los genes diferencialmente expresados del este *lindo* trabajo [*The regulatory landscape of retinoblastoma: a pathway analysis perspective*](https://doi.org/10.6084/m9.figshare.c.5975228). Por lo que usaremos la anotación de humano "org.Hs.eg.db" de nuevo. Recuerden que la lista completa está [aquí](http://bioconductor.org/packages/release/BiocViews.html#___OrgDb).

```r
#BiocManager::install("org.Hs.eg.db", character.only = TRUE)
library("org.Hs.eg.db", character.only = TRUE)
```

### Prepare Input
```r
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
```

### Gene Set Enrichment

Parámetros:

**keyType** Igual que en el ORA, el tipo de ids utilizados en nuestra lista. Receurden que los tipos permitidos se enlistan con el comado `keytypes("org.Hs.eg.db")`. Recueden que si usan algo distinto a Humano deben cambiar su anotación.   
**ont** Ontología. Alguno de "BP" (procesos biológicos), "MF" (función molecular), "CC" (componente celular) o "ALL" (todas)  
**minGSSize** tamaño mínimo de geneSet para analizar.   
**maxGSSize** tamaño máximo de genes anotados para probar. 
**pvalueCutoff** punto de corte del p-value.   
**pAdjustMethod** metodo para ajustar la p, puede ser uno de estos: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" 


```r
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = "org.Hs.eg.db", 
             pAdjustMethod = "BH")
```
Noten que hemos dado toda la lista de genes con su valor de lofFC. Estás enriqueciendo sin un criterio arbitrario (elegido por nosotros), **puro dato duro y maduro**.

## Outputs

### Tabla de resultados
```r
head(gse)
```

### Dotplot
```r
#BiocManager::install("DOSE")
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
```

### Ridgeplot

Agrupados por pathways, se generan gráficos de densidad utilizando la frecuencia del logFold Change por gen dentro de cada set. Útil para interpretar vías reguladas al alza o a la baja.

```r
ridgeplot(gse) + labs(x = "enrichment distribution")
```

### GSEA Plot  

Método tradicional para visualizar resultados de GSEA.  
  
Parámetros:

**GeneSetID** Entero. Corresponde al pathway en el objeto 'gse'. El primero pathway es el 1, el segundo el 2, *etc*. 

```r
# Usamos el objeto 'geneset' para que siempre coincidan el título y el gene set correspondiente.
geneset = 1
gseaplot(gse, by = "all", title = gse$Description[geneset], geneSetID = geneset)
```


## KEGG *Gene Set Enrichment Analysis* con `clusterProfiler`

Para hacer GSEA con KEGG usaremos la función `gseKEGG()` y tambien cambiaremmos de tipos de ids con `bitr` como en ORA (arriba). 

Como input inicial usaremos `original_gene_list` que creamos para GSEA de GO.

### Prepare Input
```r
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

```
### Creación del objeto gseKEGG

Parámetros:

**organism** Cómo en el ORA (arriba) código del organismo de KEGG. En este caso también "hsa" de humano. Lista completa: https://www.genome.jp/kegg/catalog/org_list.html (need the 3 letter code).  
**ont** Ontología. Alguno de "BP" (procesos biológicos), "MF" (función molecular), "CC" (componente celular) o "ALL" (todas)  
**minGSSize** tamaño mínimo de geneSet para analizar.   
**maxGSSize** tamaño máximo de genes anotados para probar. 
**pvalueCutoff** punto de corte del p-value.   
**pAdjustMethod** metodo para ajustar la p, puede ser uno de estos: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none" 

```r
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = "hsa",
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")
```

```r
head(kk2, 10)
```

### Dotplot
```r
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
```

### Ridgeplot

Como en GSEA GO, útil para interpretar vías reguladas al alza o a la baja.

```r
ridgeplot(kk2) + labs(x = "enrichment distribution")
```

## GSEA Plot  

Método tradicional para visualizar resultados de GSEA.  
  
Parámetros:

**GeneSetID** Entero. Corresponde al pathway en el objeto 'gse'. El primero pathway es el 1, el segundo el 2, *etc*. 

```r
# Usamos el objeto 'geneset' para que siempre coincidan el título y el gene set correspondiente.
geneset = 1
gseaplot(kk2, by = "all", title = gse$Description[geneset], geneSetID = geneset)
```



## Pathview

Obviamente podemos usar Pathview para mapear los pathways que salieron significativamente representados en nuestro GSEA.

Parámetros:

**gene.data** Esta es la lista `gene_list` creada arriba, el tipo de ids tiene que coincidir con el de `gene.idtype`
**pathway.id** Tenemos que elegir aquí nosotros alguna. Como la idea es visualizar nuestras vías significativamente enriquecidas, pordemos usar los ids de los pathways que nos sale en `head(kk2, 10)`.  
**species** El id de `organism` de la función `enrichKEGG`, en este caso "hsa"
**gene.idtype** el tipo de ids utilizados en `gene.data` pero sacado del objeto `gene.idtype.list`. En este caso tomamos el tercer elemento de esta lista.

```r
#BiocManager::install("pathview")
library(pathview)

# Produce una gráfica de KEGG (PNG) Usaremos hsa05322 proque el primer pathway más enriquecido fue el mismo que graficamos arriba
hsa <- pathview(gene.data=gene_list, pathway.id="hsa05322", species = "hsa", gene.idtype=gene.idtype.list[3])

# Produce una gráfica diferente (PDF)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa05322", species = "hsa", gene.idtype=gene.idtype.list[3], kegg.native = FALSE)
```

Las imágenes se salvan en su directorio de trabajo.


### Más recursos:

### Sobre clusterProfiler 

Sitio web: https://bioconductor.org/packages/release/bioc/vignettes/clusterProfiler/inst/doc/clusterProfiler.html

Tutorial: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html

### Sobre enriquecimiento

artículo "Pathway size matters: the influence of pathway granularity on over-representation (enrichment analysis) statistics" https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-021-07502-8


Tutorial "Functional Analysis for RNA-seq" https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/10_FA_over-representation_analysis.html