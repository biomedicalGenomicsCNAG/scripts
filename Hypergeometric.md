# Hypergeometric distribution in R

`cells<-table(sub$Tissue, Idents(sub))`

### head(cells)
```
                      CD14+ CD16+   DC Macrophages_1 Macrophages_2 Platelet
  blood                  3717   201   43            31             0      467
  brain                    15     8    0             6             0        0
  digestive system          0     0    3             0             0        0
  esophagus                 0     4    6             0             0        0
  eye                       7     0    0             9             0        0
  heart                   327   282  145           584           915        0
  immune system          6998  1684 1556            42            28       78
  kidney                   18    21    7             8             3        0
  liver                  1444   297   33           199             2        6
  lung                      3     4    1             7            14        0
  prostate gland           21     8   13            20            23        0
  skeletal muscle organ     5     9    1             0             1        0
  spleen                  538    57   60            89             0       74
  thymus                    3     8   35             9             0        0
  
  ```
  
`sign_table<-matrix(0, nrow = nrow(cells), ncol = ncol(cells))`
`sign_table= as.data.frame(sign_table)`

`for(c in 1: ncol(cells)){
 for (r in 1: nrow(cells)){
  x=1-phyper(cells[r,c]-1, sum(cells[,c]), sum(cells)- sum(cells[,c]), sum(cells[r,])) 
   sign_table[r,c]= paste(as.numeric(x))
 }
  sign_table[,c]= as.numeric(sign_table[,c])
}`

`row.names(sign_table) = row.names(cells)`
`colnames(sign_table) = colnames(cells)`

`sign_table[sign_table==0] <- 0.00001`
`mat= -log10(sign_table)`


`library(ComplexHeatmap)`
`library(gplots)`

`annotation = as.data.frame(colnames(cells))`
`names(annotation)= "Name"`
`annotation = as.data.frame(annotation)`
`unique(sub$Tissue)`


`colSums(cells)`
`rowSums(cells)`
`column_ha = HeatmapAnnotation(cells_subpopulation = anno_barplot(colSums(cells)), 
                              Clusters = annotation$Name,
                              col = list(Clusters = cols))`
`row_ha = rowAnnotation(cells_tissue = anno_barplot(rowSums(cells)))`


`Heatmap(as.matrix(mat), column_names_gp = gpar(fontsize = 10), name= "-log10(pvalue)",
        top_annotation = column_ha, right_annotation = row_ha, show_column_names = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
  if(as.matrix(mat)[i, j] > 1.30103) {
    grid.text("*", x, y)
  }
})`
