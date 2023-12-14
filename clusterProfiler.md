# Introducción a `clusterProfiler`

## Enriquecimientos 

Un análisis de enriquecimiento funcional suele aplicarse para extraer información biológica relevante ya sea de una lista filtrada de genes relevantes. El ejemplo más conocido de esta lista podría ser los genes diferencialmente expresados resultado de un análisis de Expresión diferencial.

Si bien hay varios métodos para hacer enriquecimiento[^1], los dós métodos más comunes son el Análisis de sobrerepresentación (*Over-representation analysis* o ORA) y el Análisis de Enriquecimiento de Conjuntos de Genes (*Gene Set Enrichment Analysis* o GSEA). Ambos incluidos en `clusterProfiler`.


| Análisis | Lo que require de entrada | Como luce la salida |  ✅ Pros | ⚠️ Cons |
| -------- | :-----------------------: | :-----------------: | :-----: | :-----: |
| ORA (Over-representation Analysis) | Una lista de genes (no se requeiren estadísticos) | Una prueba hipergeométrica por vía | - Sencillo <br> - Barato computacionalmente para calcular p-values | - Ignora estadísticas asociadas con los genes y requiere puntos de corte arbitrarios <br> - Asume independencia de genes y pathways |
| GSEA (Gene Set Enrichment Analysis) | Una lista de IDs de genes con un valor estadístico por gen | Un puntaje de enriquecimiento (*enrichment score* ES) por vía | - Incluye todos los genes (Sin puntos de corte arbitrarios) <br> |  - Las permutaciones pueden ser computacionalmente costosas <br> - No considera sobrelape de vías |



## *Over-representation analysis*


![ORA](https://hbctraining.github.io/DGE_workshop_salmon_online/img/go_proportions.png "Over-representation analysis")


### Enrichr

Herramienta en linea para hacer OVA
https://maayanlab.cloud/Enrichr/


### Entonces ¿Para qué queremos clusterProfiler?





[^1] Para más detalles ver [Gómez-Romero, et al. (2021) Bioinformatics of Functional Categories Enrichment in *Bioinformatics and Human Genomics Research* eBook ISBN: 9781003005926]("https://www.taylorfrancis.com/chapters/edit/10.1201/9781003005926-14/bioinformatics-functional-categories-enrichment-laura-g%C3%B3mez-romero-hugo-tovar-enrique-hern%C3%A1ndez-lemus").