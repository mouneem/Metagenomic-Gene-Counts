setwd('~/Personal/Metagenomic-Gene-Counts')
library('reshape2')
library('ggplot2')
library('viridis')
library('ggrepel')

PA.matrix <- read.table('input.Rtab',header = TRUE)
PA.matrix <- PA.matrix[ , colnames(PA.matrix) != 'IDN1' ]

PA_melted <- melt(PA.matrix, id ='Gene') # transform matrix to 2d table
PA_melted <- PA_melted[PA_melted$value > 0 , ] # keep only 1s
gene_counts <- as.data.frame(table(PA_melted$Gene))
N_genome <- length(unique( PA_melted$variable  )) # number of genomes based on the number of coloumns
genes_core <- gene_counts[gene_counts$Freq == N_genome ,  ]$Var1
write.csv(as.data.frame(genes_core), 'output/core_genes.csv')

# genes by genome
ggplot(PA_melted , aes(variable, fill = variable))+
  viridis::scale_fill_viridis(discrete=TRUE, option = 'A') +
  geom_bar(width = .5) +
  theme(text = element_text(size=8), axis.text.x = element_text(angle=45, hjust=1)) 

gg <- as.data.frame(table(PA_melted[,c('variable')]) )
ggplot(gg , aes(reorder(Var1 , -Freq), Freq, fill = Var1))+
  viridis::scale_fill_viridis(discrete=TRUE, option = 'A') +
  geom_bar(width = .5, stat = 'identity') +
  ggrepel::geom_label_repel(aes( label = Var1), 
                            label.size = 0 , 
                            direction = 'y' ) +
  theme_bw()+
  theme(text = element_text(size=8), axis.text.x = element_text(angle=45, hjust=1))

write.csv(as.data.frame(gg), 'output/genes_by_genome.csv')

# specific
genes_unique <- gene_counts[gene_counts$Freq == 1 ,  ]$Var1
genes_unique_forGenomes <- PA_melted[PA_melted$Gene %in% genes_unique, ]
genes_unique_forGenomes <- as.data.frame(table(genes_unique_forGenomes$variable))
write.csv(as.data.frame(genes_uniques), 'output/genes_unique_forGenomes.csv')


# accessory
gene_accessory <- gene_counts[gene_counts$Freq > 1 & gene_counts$Freq < N_genome,  ]
gene_accessory_counts <- PA_melted[PA_melted$Gene %in% gene_accessory$Var1, ]
gene_accessory_counts <- as.data.frame(table(gene_accessory_counts$variable))
colnames(gene_accessory_counts) <- c('Var1','Accessory')
write.csv(as.data.frame(genes_uniques), 'output/genes_accessory.csv')


all.data <- merge(genes_unique_forGenomes , gene_accessory_counts)
colnames(all.data) <- c('Var1', 'Specific','Accessory')
all.data$core <- length(genes_core)
all.data_m <- melt(all.data, id = 'Var1')

write.csv(as.data.frame(genes_uniques), 'output/all_types.csv')

ggplot(data = all.data_m, aes(x = Var1, y = value, fill = variable) )+ 
  geom_bar(position="stack", stat="identity")+
  viridis::scale_fill_viridis(discrete=TRUE, option = 'E') +
  theme_bw()+
  theme(text = element_text(size=8), axis.text.x = element_text(angle=45, hjust=1))
# counts
length(genes_core) #1900
length(genes_uniques$value) #8356
length(gene_accessory$Var1) #5693


