library(dplyr)
library(ggplot2)
library(ggsci)
library(lattice)
library(micropan)
library(shipunov)
library(ape)
library(CollessLike)
library(wesanderson)
library(venn)
library(scales)
library(tidytext)
library(ggfortify)
library(ggtree)
library(hgu95av2.db)
library(package = affyLib, character.only = TRUE)
library(topGO)


genome_array_matrix <- read.table('/home/anton/serratia_pangenomics_2022/phylogeny/mash_serratia.tsv', row.names=1,header = TRUE,  sep = "\t")

genome_array_matrix <- 1-genome_array_matrix 

bootsrap_clusters <- Bclust(as.matrix(genome_array_matrix),method.d="euclidean", method.c="complete", iter=1000)
plot(bootsrap_clusters, labels=FALSE)
dend <- as.phylo(bootsrap_clusters$hclust)
dend$node.label <- bootsrap_clusters$values
dend$tip.label

write.tree(dend,file="hclust_mash_fna.nwk")

drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    # panel.text(x, y, e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet,
                   at=seq(min(df), max(df), length.out=100),
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1)),
                   scales=list(cex=0,tck = c(0,0)), xlab=list(label=''),
                   ylab=list(label=''), panel=myPanel_a))
}

mypal = colorRampPalette(c('#FAFAFA','#FF4343'))
jet = mypal(100)

mash_tree <-read.tree("/home/anton/serratia_pangenomics_2022/phylogeny/hclust_mash_fna.nwk")
mash_tree$tip.label

ggtree(mash_tree) 

drawCoolHM(as.matrix(genome_array_matrix[mash_tree$tip.label, mash_tree$tip.label]))


##Assemblies' properties
Asseblies_table = read.table('/home/anton/serratia_pangenomics_2022/assembly_data/serratia_assembls_tab_fixed_full_annotations_no_dups.csv',  
                             header=T, sep=',', stringsAsFactors = F, comment.char = "", check.names = FALSE, quote = "")

colnames(Asseblies_table) 
mean(Asseblies_table$Total_length) #5297794
mean(Asseblies_table$CDS_num) #4843
mean(Asseblies_table$Mean_CDS_size) # 362
mean(Asseblies_table$Num_hypothetical) #313
mean(Asseblies_table$GC_content) #362

mean(Asseblies_table$GC_content)# 0.595
ggplot(Asseblies_table, aes(x=reorder(Assembly,-Mean_CDS_size) , y=Mean_CDS_size)) +
  geom_point(stat='identity',alpha=0.7, shape=21, color='black', size=4, fill='darkgreen')+
  #geom_bar(stat='identity', alpha=0.9, color='black')+
  theme_bw()+
  xlab('Assembly') +
  ylab('Mean CDS size, a.a.')+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         strip.text = element_text(size = 20))+
  ylim(250, 330)+ 
  scale_fill_gradientn("GC content, %", colours=c( 'white','grey','#9C8881', '#AD3434'),na.value = "transparent",
                       breaks=c(0.58, 0.59, 0.60),labels=c(0.58,0.59, 0.60),
                       limits=c(0.58,0.602))+
  ylim(3462811, 5842069) #'#FF7C00', '#315015' # 'red', 'green' 
# 9//4

ggplot(Asseblies_table, aes(x=N50, y=GC_content, fill=GC_content, size=Contigs)) +
  geom_point(stat='identity',alpha=0.7, shape=21, color='black')+ # size=5
  theme_bw()+
  xlab('N50') +
  ylab('GC content, %')+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         strip.text = element_text(size = 20))+ 
  scale_fill_gradientn("GC content, %", colours=c( 'white','grey','#9C8881', '#AD3434'),na.value = "transparent",
                       breaks=c(0.58, 0.59, 0.60),labels=c(0.58,0.59, 0.60),
                       limits=c(0.58,0.602))+
  geom_smooth(color='black', method='lm')+ guides(size = "none")

# 8//3

all_plot <- ggplot(Asseblies_table, aes(x=reorder(Assembly, CDS_num),y=CDS_num)) + 
  geom_bar(position="stack",stat="identity",fill='grey')

all_plot + geom_bar(data=Asseblies_table, aes(x=reorder(Assembly, CDS_num), 
                                         y=Num_hypothetical), position="stack",stat="identity",fill="#083775") +
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background =  element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        axis.title.x = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12))+
  ylab('Number of CDS')+
  xlab("Assembly")+
  guides(fill= guide_legend(title="Protein sequences"))
# 9//4

#U-curves
panaroo_presence <- read.table('/home/anton/serratia_pangenomics_2022/pangenomes/panaroo_pangenome_res_num.csv',  
                               header=T, sep='\t', stringsAsFactors = F, quote = "")
roary_presence <- read.table('/home/anton/serratia_pangenomics_2022/pangenomes/roary_pangenome_res_num.csv',  
                             header=T, sep='\t', stringsAsFactors = F, quote = "")
peppan_presence <- read.table('/home/anton/serratia_pangenomics_2022/pangenomes/peppan_pangenome_res_num.csv',  
                              header=T, sep='\t', stringsAsFactors = F, quote = "")

num_stat_paparoo <- as.data.frame(panaroo_presence[,1] %>%  table())
num_stat_paparoo$instr='Panaroo'
colnames(num_stat_paparoo) <- c('num', 'Freq', 'Instr')

num_stat_roary <- as.data.frame(roary_presence[,1] %>%  table())
num_stat_roary$instr='Roary'
colnames(num_stat_roary) <- c('num', 'Freq', 'Instr')

num_stat_peppan<- as.data.frame(peppan_presence[,1] %>%  table())
num_stat_peppan$instr='PEPPAN'
colnames(num_stat_peppan) <- c('num', 'Freq', 'Instr')

all_instr_curves <- rbind(num_stat_paparoo, num_stat_roary, num_stat_peppan)
write.table(all_instr_curves, 'Data_for_U_curves.csv', sep='\t', row.names = F)

ggplot(all_instr_curves, aes(x=as.numeric(num) , y=Freq, fill=Instr)) +
  geom_bar(stat="identity",width=0.9, alpha=0.8, position="dodge")+ #, col='black'
  theme_bw()+ scale_fill_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of genes')+
  xlab("Number of genomes")+
  guides(fill= guide_legend(title="Tool"))
  #facet_wrap(~Instr)

# 12/4 - all
# default - Panaroo only

# Heaps law

panaroo_presence_matr <- panaroo_presence[,-c(1:4)] %>%  as.data.frame() #panaroo
roary_presence_matr <- roary_presence[,-c(1:15)] %>%  as.data.frame() #roary
peppan_presence_matr <- peppan_presence[,-c(1:2)] %>%  as.data.frame() #PEPPAN

peppan_presence_matr[peppan_presence_matr!='' ]<-1
peppan_presence_matr[peppan_presence_matr=='' ]<-0
pep_matrix <- matrix(as.numeric(unlist(peppan_presence_matr)),nrow = nrow(peppan_presence_matr),ncol = 73)
t_pep <- t(pep_matrix)

roary_presence_matr[roary_presence_matr!='' ]<-1
roary_presence_matr[roary_presence_matr=='' ]<-0
roary_matrix <- matrix(as.numeric(unlist(roary_presence_matr)),nrow = nrow(roary_presence_matr),ncol = 73)
t_roary<- t(roary_matrix)

panaroo_presence_matr[panaroo_presence_matr!='' ]<-1
panaroo_presence_matr[panaroo_presence_matr=='' ]<-0
pan_matrix <- matrix(as.numeric(unlist(panaroo_presence_matr)),nrow = nrow(panaroo_presence_matr),ncol = 73)
t_pan <- t(pan_matrix)

heaps_res_pan <- heaps(t_pan, n.perm = 1000) #0.6411718 
heaps_res_roary <- heaps(t_roary, n.perm = 1000) #0.5682325 
heaps_res_pep <- heaps(t_pep, n.perm = 1000) #0.6413087 

std <- function(x) sd(x,na.rm = T)/sqrt(length(x))

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = std(x[[col]]))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

make_pangenome_permutations <- function(t_pan){
  permut_df <- data.frame(genome_num=c(), gene_sum=c())
  for (perm_iter in 1:73){
    perms <- sample(c(1:73), 1)
    perm_sums <-sum(t_pan[perms,]==1)
    permut_row=data.frame(genome_num=1, gene_sum=perm_sums)
    permut_df <- rbind(permut_df, permut_row)
  }
  
  for (i in 2:73){
    print(i)
    for (perm_iter in 1:100){
      perms <- sample(c(1:73), i)
      perm_sums <-sum(colSums(t_pan[perms,]!= 0)!=0)
      permut_row=data.frame(genome_num=i, gene_sum=perm_sums)
      permut_df <- rbind(permut_df, permut_row)
    }
  }
  permut_stat <- data_summary(permut_df, varname="gene_sum", 
                              groupnames='genome_num')
  
  permut_stat$up <- permut_stat$gene_sum + permut_stat$sd
  permut_stat$down <- permut_stat$gene_sum - permut_stat$sd
  return(permut_stat)
}

panaroo_perms <- make_pangenome_permutations(t_pan)
roary_perms <- make_pangenome_permutations(t_roary)
peppan_perms <- make_pangenome_permutations(t_pep)

panaroo_perms$Tool <- 'Panaroo'
roary_perms$Tool  <- 'Roary'
peppan_perms$Tool <- 'PEPPAN'

permut_stat <- rbind(panaroo_perms, roary_perms, peppan_perms)

ggplot(permut_stat, aes(x=genome_num, y=gene_sum, col=Tool)) +
  geom_line(size = 1.3)+ scale_color_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  geom_errorbar(aes(x=genome_num, y=gene_sum, ymax=up, ymin=down),width = 0.35, size=1.1)+
  theme_bw()+
  theme(
    axis.title.x=element_text(face="bold", color="black", 
                              size=20),
    panel.background = element_blank(), 
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(color='black', 
                               size=16),
    axis.title.y = element_text(face="bold", color="black", 
                                size=18),
    legend.title=element_text(face="bold",size=14), 
    legend.text=element_text(size=12),
    axis.text.x = element_text(color='black', 
                               size=16))+
  ylab('Number of genes')+
  xlab("Number of genomes")

#9//4

#Pangenome genes


#Panaroo: 3700 - core, 1530 -shell, 15166-total
#PEPPAN:3600, 1747, 15188
#Roary:	3341,	2530, 20380


Gene_core_pangenome=data.frame(genes=c(15166,15188,20380,3700,3600,3341,1530,1747,2530),
                               gene_type =c('Total','Total','Total','Core','Core','Core',
                                            'Shell','Shell', 'Shell'),
                               Instr=c('Panaroo','PEPPAN','Roary','Panaroo','PEPPAN','Roary','Panaroo','PEPPAN','Roary'))


ggplot(Gene_core_pangenome, aes(x=gene_type , y=genes, fill=Instr)) +
  geom_bar(position="dodge", stat="identity",width=0.9, alpha=0.7, col='black')+
  theme_bw()+ scale_fill_manual(values = c("#a60b0b", '#2980b9' ,'#bdbd00'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of genes')+
  xlab("Genes' type")+
  guides(fill= guide_legend(title="Tool"))

Gene_core_panaroo <- Gene_core_pangenome[Gene_core_pangenome$Instr=='Panaroo' & Gene_core_pangenome$gene_type!='total',]
Gene_core_panaroo <- rbind(Gene_core_panaroo, data.frame(genes=15166-3700-1530, gene_type='Accessory', Instr='Panaroo'))
Gene_core_panaroo$gene_type <- c('Core','Shell', 'Accessory')
ggplot(Gene_core_panaroo, aes(x='',y=genes, fill=gene_type)) +
  geom_bar(stat="identity", width=1, alpha=0.55, col='black')+
  coord_polar(theta = "y")+
  theme_bw()+
  geom_text(aes(label = genes),
            position = position_stack(vjust = 0.5), size=9)+
  geom_text(aes(label = genes),
            position = position_stack(vjust = 0.5), size=9)+
  geom_text(aes(label = genes),
            position = position_stack(vjust = 0.5), size=9)+
  theme_bw()+ scale_fill_lancet()+
  scale_size(guide = 'none')+
  theme( axis.text.x = element_blank(),
         axis.title.x=element_blank(),
         panel.border = element_blank(),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_blank(),
         axis.title.y = element_blank(),
         legend.title=element_text(face="bold",size=14), 
         legend.text=element_text(size=12),
         axis.ticks = element_blank(),
         axis.line.x =  element_blank(),
         axis.line.y =  element_blank())+
  guides(fill= guide_legend(title="Genes' type"))



#Tress balance and suppot
tree_table <- read.table('/home/anton/serratia_pangenomics_2022/phylogeny/trees_support.tsv',  
                         header=T, sep='\t', stringsAsFactors = F)
tree_table$Colles_Like <-0
tree_table$Sackin <-0
tree_table$Cophenetic <-0

for (i in 1:nrow(tree_table)){
  tree <- read.tree(paste0('/home/anton/serratia_pangenomics_2022/phylogeny/new_topol/',tree_table[i,2]))
  balace_vec <- as.character(c(balance.indices(tree),'reported'))
  tree_table[i,]$Colles_Like <- balace_vec[1]
  tree_table[i,]$Sackin <- balace_vec[2]
  tree_table[i,]$Cophenetic <- balace_vec[3]
} 

tree_table$Colles_Like <- as.numeric(tree_table$Colles_Like)
tree_table$Sackin <- as.numeric(tree_table$Sackin)
tree_table$Cophenetic <- as.numeric(tree_table$Cophenetic)

tree_table[7,4:6] <- c(864, 1071, 15049)
tree_table[10,4:6] <- c(861, 1004, 15078)
tree_table[8,4:6] <- c(859, 1134, 15011)
tree_table[1,4:6] <- c(656, 854, 13678)

tree_table$Tree <- c('Roary_acc', 'ANI',
                     'Panaroo_core_MAFFT', 'Panaroo_core_Prank', 'Panaroo_part_MAFFT','Panaroo_part_Prank',
                     'Roary_core_MAFFT', 'Roary_core_Prank', 'Roary_part_MAFFT','Roary_part_Prank',
                     'Panaroo_acc', 'PEPPAN_acc')

tree_table$Balance_ind <- tree_table$Colles_Like/max(tree_table$Colles_Like) + tree_table$Sackin/max(tree_table$Sackin) +tree_table$Cophenetic/max(tree_table$Cophenetic)
write.table(tree_table, 'Trees_properties.csv', sep='\t', row.names = F)

tree_table_fixed <- tree_table[c(3:6, 11),]
tree_table_fixed$Balance_ind <- tree_table_fixed$Colles_Like/max(tree_table_fixed$Colles_Like) + tree_table_fixed$Sackin/max(tree_table_fixed$Sackin) +tree_table_fixed$Cophenetic/max(tree_table_fixed$Cophenetic)

ggplot(tree_table_fixed, aes(x=Balance_ind , y=Support, fill=Tree)) +
  geom_point(shape=21, size=9.5, alpha=0.95)+
  theme_bw() + 
  scale_fill_manual(values = wes_palette("Darjeeling2", 14, type = "continuous"))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x=element_text(color='black', 
                                 size=16))+
  ylab('Mean support')+
  xlab('Total balance index')
#Default


## Tree topologies


tree_dist = read.table('/home/anton/serratia_pangenomics_2022/phylogeny/tree_topologies_new.csv', row.names=1,
                       header=T, sep=',', stringsAsFactors = F)

colnames(tree_dist) <- c('Panaroo_core_MAFFT', 'Panaroo_core_Prank', 'Panaroo_part_MAFFT','Panaroo_part_Prank',
                         'Roary_core_MAFFT', 'Roary_core_Prank', 'Roary_part_MAFFT','Roary_part_Prank',
                         'Panaroo_acc','Roary_acc', 'ANI', 'PEPPAN_acc')
rownames(tree_dist) <- colnames(tree_dist) 

#order_vec <- c('ANI', 'PEPPAN_acc', 'Panaroo_acc', 'Roary_acc', 
#                     'Panaroo_core_MAFFT', 'Panaroo_core_Prank', 'Panaroo_part_MAFFT','Panaroo_part_Prank',
#                     'Roary_core_MAFFT', 'Roary_core_Prank', 'Roary_part_MAFFT','Roary_part_Prank')


get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

tree_dist <- get_upper_tri(tree_dist)
is.na(tree_dist) < 0
melted_dist <- melt(as.matrix(tree_dist), na.rm = TRUE)

melted_dist$value <- round(melted_dist$value , 3)

melted_dist<- melted_dist[order(melted_dist$Var1),]

Pan_vect <- c('Panaroo_core_MAFFT', 'Panaroo_core_Prank', 'Panaroo_part_MAFFT','Panaroo_part_Prank', 'Panaroo_acc')
melted_dist_panaroo <- melted_dist[melted_dist$Var1 %in% Pan_vect & melted_dist$Var2 %in% Pan_vect,]

ggplot(data = melted_dist, aes(Var2, Var1, fill = value))+
  geom_tile(color = 'black')+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 4.5)+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(0.7,1), space = "Lab", 
                      name="Tree simillarity") +
  theme_bw() + xlab('Trees') + ylab('Trees') +
  theme( axis.text.y = element_text(color='black', 
                                    size=22),
         axis.title.y=element_text(color="black", 
                                   size=24),
         panel.background = element_rect(fill = "#DCDCDC"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_line(size = 0.1, linetype = 'dashed',
                                         colour = "black"),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 22, hjust = 1),
         axis.title.x = element_text( color="black", 
                                      size=24))+
  coord_fixed()
#12/7
#14/9



topol_ids_doubled<-data.frame(tree_names=c(as.character(melted_dist$Var1), as.character(melted_dist$Var2)), 
                              tree_ids=c(melted_dist$value,melted_dist$value)) %>% filter(tree_ids!=1.0)


topol_ids_doubled$tree_names=as.factor(topol_ids_doubled$tree_names)

mean_id_for_trees <- topol_ids_doubled  %>% group_by(tree_names) %>% summarise(tree_ids=mean(tree_ids),num=n()) %>% as.data.frame()

median(topol_ids_doubled$tree_ids) #0.9215
mean(topol_ids_doubled$tree_ids) #0.924
max(topol_ids_doubled$tree_ids) #0.987
min(topol_ids_doubled$tree_ids) #0.831
mean_id_for_trees[3,2] <- 0.9447375

melted_dist_core <- melted_dist[!melted_dist$Var1 %in% c('ANI', 'Panaroo_acc','PEPPAN_acc','Roary_acc') & !melted_dist$Var2 %in% c('ANI', 'Panaroo_acc','PEPPAN_acc','Roary_acc'),]

topol_ids_doubled<-data.frame(tree_names=c(as.character(melted_dist[]$Var1), as.character(melted_dist$Var2)), 
                              tree_ids=c(melted_dist$value,melted_dist$value)) %>% filter(tree_ids!=1.0)

ggplot(data = mean_id_for_trees, aes(tree_names, tree_ids))+
  geom_bar(stat = 'identity', aes(x = reorder(tree_names, -tree_ids)),
           color='black',fill='#A60B0B',alpha=0.75)+
  theme_bw() + ylab('Mean similarity') + xlab('Trees') +
  geom_hline(yintercept = 0.9215,linetype="dashed", color = "#2980B9",size =2.5)+
  theme( axis.text.y = element_text(color='black', 
                                    size=24),
         axis.title.y=element_text( color="black", 
                                    size=28),
         panel.background = element_blank(),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         legend.position = "none",
         axis.text.x = element_text(color='black', 
                                    angle = 60, vjust = 1, 
                                    size = 24, hjust = 1),
         axis.title.x = element_text(color="black", 
                                     size=28)
  ) +  scale_y_continuous(limits=c(0.8,1),oob = rescale_none)

#10/7

#Host Venn diagramms
Host_distributions <- read.table('/home/anton/serratia_pangenomics_2022/host_annotations/gene_groups_distibution/genes_to_host_attributions_all.csv', 
                                 header=T, sep=',', stringsAsFactors = F)


x <- list(
  'None' = Host_distributions[Host_distributions$Host=='None',1], 
  'Human' = Host_distributions[Host_distributions$Host=='Human',1], 
  'Insect' = Host_distributions[Host_distributions$Host=='Insect',1],
  'Plant' = Host_distributions[Host_distributions$Host=='Plant',1]
)



venn(x, ilab=TRUE, zcolor = c('#53B6B6','#B17129','#3E1881','#2F690B'),ellipse = T,opacity=0.7, ilcs = 1.6, sncs = 1.5)


Host_distributions_scoary <- read.table('/home/anton/serratia_pangenomics_2022/host_annotations/gene_groups_distibution/genes_to_host_attributions_scoary_attr.csv', 
                                        header=T, sep=',', stringsAsFactors = F)

scoary_hits_all <- read.table('/home/anton/serratia_pangenomics_2022/host_annotations/scoary_results/Filtered_results/All_filt_hits_pos.csv', 
                              header=T, sep=',', stringsAsFactors = F, quote = "")

x <- list(
  'Human' = unique(Host_distributions_scoary[Host_distributions_scoary$Host=='Human',1]), 
  'Insect' = unique(Host_distributions_scoary[Host_distributions_scoary$Host=='Insect',1]),
  'Plant' = unique(Host_distributions_scoary[Host_distributions_scoary$Host=='Plant',1])
)


x <- list(
  'Human' = scoary_hits_all[scoary_hits_all$Host=='Human',1], 
  'Insect' = scoary_hits_all[scoary_hits_all$Host=='Insect',1],
  'Plant' = scoary_hits_all[scoary_hits_all$Host=='Plant',1]
)

venn(x, zcolor = c('#B17129','#3E1881','#2F690B'), opacity=0.7, ilcs = 1.6, sncs = 1.)

#Host statistics of frequency
Host_stat_tab_filt <- scoary_hits_all %>% group_by(num,Host) %>% dplyr::summarize(freq=n()) %>% as.data.frame()
Host_stat_tab$Host <- as.factor(Host_stat_tab$Host )
levels(Host_stat_tab$Host) <- c('Human','Insect','Plant')

Host_stat_all <- left_join( Host_distributions_scoary, scoary_hits_all[,c(1,7)], "Gene") %>% group_by(num,Host) %>% dplyr::summarize(freq=n()) %>% as.data.frame()
Host_stat_mean<- scoary_hits_all %>% group_by(Host) %>% dplyr::summarize(m=mean(num))
# human - 24, insect - 9, plant - 4

ggplot(scoary_hits_all, aes(x=num, fill=Host))+
  geom_density(color="black", ,alpha=0.7)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Frequency')+
  xlab("Number of genomes") +scale_fill_manual(values = c('#B17129','#3E1881','#2F690B'))+
  guides(fill= guide_legend(title="Host"))

ggplot(Host_stat_tab_filt[Host_stat_tab_filt$Host=='Human',], aes(x=num,y= freq, fill=Host)) + 
  geom_bar(stat="identity",alpha=0.5) +
  geom_bar(data=Host_stat_tab_filt[Host_stat_tab_filt$Host=='Insect',],aes(x=num,y= freq, fill=Host), stat="identity",alpha=0.7)+
  geom_bar(data=Host_stat_tab_filt[Host_stat_tab_filt$Host=='Plant',],aes(x=num,y= freq, fill=Host), stat="identity",alpha=0.7)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of genes')+
  xlab("Number of genomes") +scale_fill_manual(values = c('#B17129','#3E1881','#2F690B'))+
  guides(fill= guide_legend(title="Host"))

#Percent of hypotherical proteins in Scoary output

scoary_hits_all$hypoth_factor <- ifelse(scoary_hits_all$Annotation=='hypothetical protein', 'Hypothetical', 'Non_hypothetical')
hypothetical_scoary <- scoary_hits_all %>% group_by(hypoth_factor,Host) %>% dplyr::summarize(freq=n()) %>% as.data.frame()

hypothetical_scoary <- rbind(hypothetical_scoary, c('All','Human', sum(hypothetical_scoary[hypothetical_scoary$Host=='Human', 3])))
hypothetical_scoary$freq=as.numeric(hypothetical_scoary$freq)
hypothetical_scoary <- rbind(hypothetical_scoary, c('All','Insect',sum(hypothetical_scoary[hypothetical_scoary$Host=='Insect', 3])))
hypothetical_scoary$freq=as.numeric(hypothetical_scoary$freq)
hypothetical_scoary <- rbind(hypothetical_scoary, c('All','Plant',sum(hypothetical_scoary[hypothetical_scoary$Host=='Plant', 3])))
hypothetical_scoary$freq=as.numeric(hypothetical_scoary$freq)

hypothetical_scoary <- hypothetical_scoary[hypothetical_scoary$hypoth_factor!='Non_hypothetical',]
hypothetical_scoary <- reshape(hypothetical_scoary, idvar = "Host", timevar = "hypoth_factor", direction = "wide")

hypothetical_scoary$ratio <- hypothetical_scoary$freq.Hypothetical/hypothetical_scoary$freq.All

all_plot <- ggplot(hypothetical_scoary, aes(x=reorder(Host, freq.All),y=freq.All)) + 
  geom_bar(position="stack",stat="identity",fill='grey')

all_plot + geom_bar(data=hypothetical_scoary, aes(x=reorder(Host, freq.All), 
                                              y=freq.Hypothetical), position="stack",stat="identity",fill="#083775") +
  theme_bw()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background =  element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.text.x = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        axis.title.x = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12))+
  ylab('Number of clusters')+
  xlab("Host")+
  guides(fill= guide_legend(title="Protein sequences"))

#Number of fsctors per assembly
num_factors_per_asmbl  <- read.table('/home/anton/serratia_pangenomics_2022/host_annotations/Num_gene_per_host_per_assembly_ordered_lables.csv', 
                                     header=T, sep='\t', stringsAsFactors = F)
#Num_clusters_Human, Num_clusters_Plants, Num_clusters_Insect Percent_human Percent_insect Percent_plant

num_factors_per_asmbl %>% group_by(Host) %>% summarise()

human_df <- data.frame(Host=num_factors_per_asmbl$Host, num_clust =num_factors_per_asmbl$Num_clusters_Human, factor='Human_clusters')
insect_df <- data.frame(Host=num_factors_per_asmbl$Host, num_clust =num_factors_per_asmbl$Num_clusters_Insect, factor='Insect_clusters')
plant_df <- data.frame(Host=num_factors_per_asmbl$Host, num_clust =num_factors_per_asmbl$Num_clusters_Plants, factor='Plant_clusters')

per_cluster_nums <- rbind(human_df, insect_df, plant_df)
mean_clust_dist <- per_cluster_nums %>% group_by(Host, factor) %>% summarise(mean_clust=median(num_clust)) %>% as.data.frame()

#Boxplot mean numbers of Scoary factors
ggplot(per_cluster_nums, aes(y=num_clust, x=Host, fill=factor)) +
  stat_boxplot()+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of clusters')+
  xlab("Host")+scale_fill_manual(values = c('#B17129','#3E1881', '#2F690B'))+
  guides(fill= guide_legend(title="Host"))

#Geom point of Scoary factors' dependance
ggplot(num_factors_per_asmbl, aes(x=Percent_plant, y=Percent_human, fill=Host)) +
  geom_point(alpha=0.8, col='black', shape=21, size=5)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Insect-attributed genes')+
  xlab("Human-attributed genes")+
  scale_size(range = c(2, 10))+scale_fill_manual(values = c('#B17129','#3E1881', '#53B6B6','#2F690B'))+
  guides(fill= guide_legend(title="Host"))


#Virulence factors adundance and identity percent
VFDB_results <- read.table('/home/anton/serratia_pangenomics_2022/virulence/VirDB_top_hits_cov_70_id_70.csv', 
                                 header=T, sep='\t', stringsAsFactors = F)

mean(VFDB_results$Genome_num) #67
sum(VFDB_results_filt$Genome_num>=69) #6714
# 6714/7955 84%

ggplot(VFDB_results_filt, aes(x=Genome_num)) +
  geom_bar(position="dodge", stat="count",width=1, alpha=0.8, col='black')+
  theme_bw()+ scale_fill_manual(values = c('#bdbd00'))+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of virulence genes')+
  xlab("Number of genomes")

VFDB_results_ident <- VFDB_results_filt %>% group_by(Genome_num) %>% summarize(mean_id=mean(pident), num=n())
mean(VFDB_results_filt$pident)

ggplot(VFDB_results_ident, aes(x=Genome_num, y=mean_id, size=num, fill=mean_id)) +
  geom_point(alpha=0.9, col='black', shape=21)+
  theme_bw()+
  scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(68,100), space = "Lab", 
                                                     name="Sequence similarity")+
  geom_hline(yintercept = 81.145, color='grey', linetype='dotted',size=2)+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Identity percent')+
  xlab("Number of genomes")+
  scale_size(range = c(2, 10))+
  guides(size = guide_legend(title="Number of genes"))

#Refernce phylogeny
ref_tree <- read.tree('/home/anton/serratia_pangenomics_2022/phylogeny/best_reference_tree.nwk')

p <- ggtree(ref_tree,alpha =1) %<+% num_factors_per_asmbl
p <-p+ scale_y_reverse()

p+ 
  geom_tippoint(aes(color=Host), size=4)+
  scale_color_manual(values = c('#B17129','#3E1881', '#53B6B6','#2F690B'))+
  geom_tiplab(size=0, show.legend=T, color='black') #6//9 - no tips

ref_plot_data <- p$data %>% as.data.frame()
ref_plot_data_filtered <- ref_plot_data[c(1:73),] %>% arrange(y) %>% select(-c(parent,node,branch,branch.length,isTip,x,y,angle))
colnames(ref_plot_data_filtered)[1] <- 'Assembly'

ref_tree_order <- ref_plot_data_filtered$Assembly


write.table(ref_plot_data_filtered, 'Num_gene_per_host_per_assembly_ordered_lables.csv', sep='\t', row.names = F)

#Number of host specificity factors for heatmap

#Absolute numbers
human_points <- data.frame(Assembly=num_factors_per_asmbl$Assembly, num_clust =num_factors_per_asmbl$Num_clusters_Human, factor='Humen')
insect_points <- data.frame(Assembly=num_factors_per_asmbl$Assembly, num_clust =num_factors_per_asmbl$Num_clusters_Insect, factor='Insect')
plant_points <- data.frame(Assembly=num_factors_per_asmbl$Assembly, num_clust =num_factors_per_asmbl$Num_clusters_Plants, factor='Plant')

#Relative fold
human_points_rel <- data.frame(Assembly=num_factors_per_asmbl$Assembly, 
                           num_clust =num_factors_per_asmbl$Num_clusters_Human/(num_factors_per_asmbl$Num_clusters_Human+
                                                                                  num_factors_per_asmbl$Num_clusters_Insect+
                                                                                  num_factors_per_asmbl$Num_clusters_Plants), factor='Human')
insect_points_rel <- data.frame(Assembly=num_factors_per_asmbl$Assembly, 
                           num_clust =num_factors_per_asmbl$Num_clusters_Insect/(num_factors_per_asmbl$Num_clusters_Human+
                                                                                  num_factors_per_asmbl$Num_clusters_Insect+
                                                                                  num_factors_per_asmbl$Num_clusters_Plants), factor='Insect')
plant_points_rel <- data.frame(Assembly=num_factors_per_asmbl$Assembly, 
                           num_clust =num_factors_per_asmbl$Num_clusters_Plants/(num_factors_per_asmbl$Num_clusters_Human+
                                                                                  num_factors_per_asmbl$Num_clusters_Insect+
                                                                                  num_factors_per_asmbl$Num_clusters_Plants), factor='Plant')

per_cluster_nums_point <- rbind(human_points, insect_points, plant_points)
per_cluster_nums_point <- per_cluster_nums_point %>% arrange(factor(per_cluster_nums_point$Assembly, levels = ref_tree_order)) %>% unique()
per_cluster_nums_point$variable <- with(per_cluster_nums_point,factor(Assembly,levels = rev(ref_tree_order)))

ggplot(per_cluster_nums_point[per_cluster_nums_point$num_clust!=0,], aes(x=factor, y=variable, fill=factor, size=num_clust)) +
  geom_point(alpha=0.8, col='black', shape=21)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 70, vjust = 1, 
                                   size = 18, hjust = 1),
        axis.ticks.y= element_blank(),
        legend.position = 'none')+
  ylab('Assembly')+
  xlab("Host")+
  scale_fill_manual(values = c('#B17129','#3E1881','#2F690B'))+ 
         coord_fixed(ratio = 0.4)+scale_size(range = c(0.1, 7)) #7/9

#Insect_Plant_Filt_scoary Human_Insect_Filt_scoary Human_Plant_Filt_scoary

#Number of virulence factors per gene group
Gene_groups_df <- read.table('/home/anton/serratia_pangenomics_2022/Gene_groups_attributions_all_pos.csv', 
                                 header=T, sep='\t', stringsAsFactors = F)

vir_percent_summary <- function(group_df, names_vec){
  group_df$vir_flag <- ifelse(is.na(Gene_groups_df$Vir_ID), 'non', 'vir')
  gene_df <- group_df[group_df$Group %in% names_vec,]
  group_table <- gene_df %>% group_by(Group, vir_flag) %>% summarize(num=n()) %>% as.data.frame()
  group_table$Group = factor(group_table$Group, levels=names_vec)
  return(group_table)
}

pan_genes_tab <- vir_percent_summary(Gene_groups_df, c('Core', 'Accessory','Unique'))
scoary_host_genes_tab_filt <- vir_percent_summary(Gene_groups_df, c('Human_Filt_scoary',
                                                                   'Insect_Filt_scoary', 'Plant_Filt_scoary'))

scoary_host_genes_tab_strict_filt <- vir_percent_summary(Gene_groups_df, c('Human_StrictFilt_scoary',
                                                                    'Insect_StrictFilt_scoary', 'Plant_StrictFilt_scoary'))


#Number og genes' plots
vir_percent_plot <- ggplot(scoary_host_genes_tab_strict_filt[scoary_host_genes_tab_strict_filt$vir_flag=='non',], aes(x=Group,y=num)) + 
  geom_bar(position="stack",stat="identity",fill='grey')

vir_percent_plot + 
  geom_bar(data=scoary_host_genes_tab_strict_filt[scoary_host_genes_tab_strict_filt$vir_flag=='vir',],aes(x=Group,y=num), 
           position="stack",stat="identity",fill="#92310F") +
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of genes')+
  xlab("Gene Group")




#Percent of virulence genes' plots
percent_vec <- scoary_host_genes_tab_strict_filt[scoary_host_genes_tab_strict_filt$vir_flag=='vir',3]/scoary_host_genes_tab_strict_filt[scoary_host_genes_tab_strict_filt$vir_flag=='non',3]*100
percent_df <- data.frame(host=c('Human', 'Insect', 'Plant'),percent=percent_vec)
ggplot(percent_df, aes(x=host,y=percent)) + 
  geom_bar(position="stack",stat="identity",fill='#92310F')+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Percent of virulence factors')+
  xlab("Host")

for_tile_scoary <- Gene_groups_df[Gene_groups_df$Group %in% c('Human_Filt_scoary',
                                                              'Insect_Filt_scoary', 'Plant_Filt_scoary','Insect_Plant_Filt_scoary',
                                                              'Human_Insect_Filt_scoary','Human_Plant_Filt_scoary'),]


for_tile_scoary <- for_tile_scoary %>% arrange(factor(for_tile_scoary$Genome, levels = ref_tree_order)) %>% unique()
for_tile_scoary$variable <- with(for_tile_scoary,factor(Genome,levels = rev(ref_tree_order)))
for_tile_scoary$pseudo_num <- for_tile_scoary$Num_genomes
for_tile_scoary[for_tile_scoary$Group=='Human_Filt_scoary',9] <- for_tile_scoary[for_tile_scoary$Group=='Human_Filt_scoary',9]+1000
for_tile_scoary[for_tile_scoary$Group=='Insect_Filt_scoary',9] <- for_tile_scoary[for_tile_scoary$Group=='Insect_Filt_scoary',9]+100


for_tile_virulence <- Gene_groups_df[Gene_groups_df$Group %in% c('Vir_All'),]
for_tile_virulence <- for_tile_virulence %>% arrange(factor(for_tile_virulence$Genome, levels = ref_tree_order))
for_tile_virulence$variable <- with(for_tile_virulence,factor(Genome,levels = rev(ref_tree_order)))

tile_acc_genes <- Gene_groups_df[Gene_groups_df$Group %in% c('Accessory'),]
tile_acc_genes <- tile_acc_genes %>% arrange(factor(tile_acc_genes$Genome, levels = ref_tree_order)) %>% unique()
tile_acc_genes$pseudo_num <- tile_acc_genes$Num_genomes
tile_acc_genes$variable <- with(tile_acc_genes,factor(Genome,levels = rev(ref_tree_order)))

# #784555 - human_insect, #706D1A - human_plant, #374146 insect_plant #pink -all acc
?geom_tile
ggplot(data =tile_acc_genes[tile_acc_genes$Num_genomes>0, ], aes(reorder(Cluster, -Num_genomes), variable, fill = Group))+ #  tile_acc_genes[tile_acc_genes$Num_genomes>1, ] #Group Vir_ID
  geom_tile(alpha=0.9)+ 
  scale_fill_manual(values = c('#BE0C91'))+
  #scale_fill_manual(values = c('#B17129','#3E1881','#2F690B'),
  #                  labels = c("Human", "Insect",'Plant'))+
  #scale_fill_gradient(low = "#DCDCDC", high = "#A60B0B", limit = c(68,100), space = "Lab", 
  #                    name="Sequence similarity")+
  theme_bw() + xlab('Gene cluster') + ylab('Assembly') +
  theme( axis.text.y =  element_blank(),
         axis.title.y=element_text(color="black", 
                                   size=24),
         panel.background = element_rect(fill = "white"),
         panel.grid.minor = element_blank(), 
         panel.grid.major = element_blank(),
         #axis.text.x = element_text(size=8, angle = 60, vjust = 1, hjust = 1),
         axis.text.x = element_blank(),
         axis.title.x = element_text( color="black", 
                                      size=24),
         axis.ticks.x = element_blank(),
         axis.ticks.y = element_blank())+
  coord_fixed(ratio = 50) +
  guides(fill = guide_legend(title="Host"))#ratio = 10, 30 3.5 150 -acc genes
#12/9

Gene_groups_df <- read.table('/home/anton/serratia_pangenomics_2022/Gene_groups_attributions_all_pos.csv', 
                             header=T, sep='\t', stringsAsFactors = F)

mge_num <- read.table('/home/anton/serratia_pangenomics_2022/mobile_elements/all_elements_num_no_dups.csv', 
                      header=F, sep='\t', stringsAsFactors = F)
extr_vec_main <- c('Core', 'Accessory','Unique', 'Vir_Core','Vir_Accessory',
                   'Human_StrictFilt_scoary','Plant_StrictFilt_scoary','Insect_StrictFilt_scoary')
#phages IS GI
asmbl_host <- read.table('/home/anton/serratia_pangenomics_2022/host_annotations/serratia_biosample_host.csv', 
                         header=T, sep='\t', stringsAsFactors = F)

asmbl_host$NCBI <- paste0(asmbl_host$NCBI,'_genomic.gbff')
mge_num$Host <- asmbl_host[match(mge_num$V1, asmbl_host$NCBI),3]

ggplot(mge_num, aes(x=V2,y=V3)) +
  geom_point(shape=21, color='black', size=6, fill='#0F6FA2', alpha=0.75)+
  geom_smooth(method='lm', color='darkgrey', alpha=0.45)+
  scale_size(guide = 'none')+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background =  element_blank(), 
        axis.line = element_blank(),
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=18), 
        legend.text=element_text(size=16),
        axis.text.x = element_text(color='black',  size = 16))+
  guides(fill= guide_legend(title="Loci type"))+
  ylab('Number of insertions')+ #'Number of genetic islands' 'Number of insertions' 'Number of prophages'
  xlab("Number of phages")
 
median(mge_num$V2) 3
median(mge_num$V3) 10
median(mge_num$V4) 7
host_stat_phages <- mge_num %>% group_by(Host) %>% summarize(mean_phages = median(V2))
host_stat_IS <- mge_num %>% group_by(Host) %>% summarize(mean_phages = median(V3))
host_stat_GIs <- mge_num %>% group_by(Host) %>% summarize(mean_phages = median(V4))

anova <- aov(V4 ~ Host, data = mge_num)
tukey <- TukeyHSD(anova)
#c('#53B6B6','#B17129','#3E1881','#2F690B')
mge_split <- data.frame(Asmbl = mge_num$V1, Num=mge_num$V2, MGE_type= 'Phages', Host=mge_num$Host)
mge_split <- rbind(mge_split,
                   data.frame(Asmbl = mge_num$V1, Num=mge_num$V3, MGE_type= 'IS', Host=mge_num$Host))
mge_split <- rbind(mge_split,
                   data.frame(Asmbl = mge_num$V1, Num=mge_num$V4, MGE_type= 'GI', Host=mge_num$Host))

mge_split$Host <- droplevels(mge_split$Host)
mge_split$Host <- factor(mge_split$Host, levels=c( 'Human','Insect', 'Plant', 'None'))

ggplot(mge_split, aes(x=Host, y=Num, fill=Host)) + 
  geom_boxplot()+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of MGEs')+
  xlab("Host") +scale_fill_manual(values = c('#B17129','#3E1881','#2F690B', '#53B6B6'))+
  guides(fill= guide_legend(title="Host"))+
  facet_wrap(~MGE_type, scales = 'free_y')
#10/6

write.table(mge_split, 'MGEs_stat_sum_per_host.csv', sep='\t', row.names = F)

cor.test(mge_num$V2, mge_num$V3, method=c("pearson" )) #phages/IS  #0.1626121  corr 0.1693 p-value
cor.test(mge_num$V2, mge_num$V4, method=c("pearson" )) #phages/GI #0.5003556  corr 6.555e-06 p-value
cor.test(mge_num$V3, mge_num$V4, method=c("pearson" )) #IS/GI #0.8047942  corr 2.2e-16 p-value

IS_bed <- read.table('/home/anton/serratia_pangenomics_2022/mobile_elements/IS_intersect_no_dup.bed', 
                     header=F, sep='\t', stringsAsFactors = F)
phages_bed <- read.table('/home/anton/serratia_pangenomics_2022/mobile_elements/phages_intersect_genes_no_dup.bed', 
                         header=F, sep='\t', stringsAsFactors = F)
GI_bed <- read.table('/home/anton/serratia_pangenomics_2022/mobile_elements/GI_intersect_no_dup.bed', 
                     header=F, sep='\t', stringsAsFactors = F)

MGE_stat_df <- data.frame(num_MGE=c(), Gene_num=c(), MGE_type=c(), Gene_group=c())
Gene_groups_extr <-Gene_groups_df[Gene_groups_df$Group %in% extr_vec_main,]

for (group in levels(as.factor(Gene_groups_extr$Group))){
  
  total_gene_num <- length(Gene_groups_extr[Gene_groups_extr$Group==group,1])
  
  IS_group <- IS_bed[IS_bed$V8 %in% Gene_groups_extr[Gene_groups_extr$Group==group,1], ] %>% summarize(IS_num=n())
  IS_df <- data.frame(num_MGE=IS_group$IS_num, Gene_num=total_gene_num, MGE_type='IS', Gene_group=group)
  
  phages_group <- phages_bed[phages_bed$V8 %in% Gene_groups_extr[Gene_groups_extr$Group==group,1], ] %>% summarize(phages_num=n())
  phages_df <- data.frame(num_MGE=phages_group$phages_num, Gene_num=total_gene_num, MGE_type='Phages', Gene_group=group)
  
  GI_group <- GI_bed[GI_bed$V8 %in% Gene_groups_extr[Gene_groups_extr$Group==group,1], ] %>% summarize(GI_num=n())
  GI_df <- data.frame(num_MGE=GI_group$GI_num, Gene_num=total_gene_num, MGE_type='GI', Gene_group=group)
  
  
  MGE_stat_df <- rbind(MGE_stat_df, IS_df)
  MGE_stat_df <- rbind(MGE_stat_df, phages_df)
  MGE_stat_df <- rbind(MGE_stat_df, GI_df)
}

MGE_stat_df$MGE_percent <- MGE_stat_df$num_MGE/MGE_stat_df$Gene_num
MGE_stat_df <- MGE_stat_df %>% group_by(MGE_type) %>% summarize(num_MGE=num_MGE,
                                                                Gene_num=Gene_num,
                                                                Gene_group=Gene_group, 
                                                                MGE_percent=MGE_percent,
                                                                MGE_norm=num_MGE/max(num_MGE))

MGE_stat_df$Gene_group <- droplevels(MGE_stat_df$Gene_group)
MGE_stat_df$Gene_group <- factor(MGE_stat_df$Gene_group, levels=c( 'Core','Accessory',
                                                         'Unique','Vir_Core','Vir_Accessory', 'Human_StrictFilt_scoary', 
                                                         'Insect_StrictFilt_scoary','Plant_StrictFilt_scoary'))

levels(MGE_stat_df$Gene_group) <-c('Core','Accessory', 
                                'Unique','Vir_Core','Vir_Accessory', 
                                'Human', 'Insect','Plant')

#9F08A4, #D089D3, #A00017, #CD1313, #DA6969
ggplot(MGE_stat_df, aes(y=MGE_percent,x= Gene_group , fill=Gene_group)) + #log10(1+MGE_norm)
  geom_bar(stat="identity",alpha=0.9, position="dodge") +
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 18, hjust = 1))+
  ylab('Number of MGEs')+
  xlab("Gene group") +scale_fill_manual(values = c('#9F08A4', '#D089D3', '#F4D2F5', '#CD1313', '#DA6969',
    '#B17129','#3E1881','#2F690B'))+
  guides(fill= guide_legend(title="Host"))+
  facet_wrap(~MGE_type, scales = 'free_y')

#10/6

MGE_stat_per_protein <- data.frame(Prot_acc= c(), Assembly = c(), Host=c(), Group = c(), MGE = c(), ID=c())

for (protein_acc in IS_bed$V8){
  ID_raw = strsplit(IS_bed[IS_bed$V8==protein_acc,4], ";")[[1]][1]
  ID = strsplit(ID_raw, "=")[[1]][2]
  assembly_df <- Gene_groups_extr[Gene_groups_extr$Gene==protein_acc, ]
  if (nrow(assembly_df) !=0){
  for (row_num in (1:nrow(assembly_df))){
    df_sub <- assembly_df[row_num,]
    df_el <- data.frame(Prot_acc=protein_acc, Assembly= df_sub$Genome, 
                        Host= df_sub$Aseembly_host, Group=df_sub$Group, MGE='IS', ID=ID)
    MGE_stat_per_protein <- rbind(MGE_stat_per_protein,df_el)
  }
  }
}

for (protein_acc in phages_bed$V8){
  ID_raw = strsplit(phages_bed[phages_bed$V8==protein_acc,4], ";")[[1]][1]
  ID = strsplit(ID_raw, "=")[[1]][2]
  asmbl = phages_bed[phages_bed$V8==protein_acc,5]
  ID=paste(ID, asmbl, sep='_')
  assembly_df <- Gene_groups_extr[Gene_groups_extr$Gene==protein_acc, ]
  if (nrow(assembly_df) !=0){
    for (row_num in (1:nrow(assembly_df))){
      df_sub <- assembly_df[row_num,]
      df_el <- data.frame(Prot_acc=protein_acc, Assembly= df_sub$Genome, 
                          Host= df_sub$Aseembly_host, Group=df_sub$Group, MGE='Phages', ID=ID)
      MGE_stat_per_protein <- rbind(MGE_stat_per_protein,df_el)
    }
  }
}

for (protein_acc in GI_bed$V8){
  ID_raw = strsplit(GI_bed[GI_bed$V8==protein_acc,4], ";")[[1]][1]
  ID = strsplit(ID_raw, "=")[[1]][2]
  asmbl = GI_bed[GI_bed$V8==protein_acc,5]
  ID=paste(ID, asmbl, sep='_')
  assembly_df <- Gene_groups_extr[Gene_groups_extr$Gene==protein_acc, ]
  if (nrow(assembly_df) !=0){
    for (row_num in (1:nrow(assembly_df))){
      df_sub <- assembly_df[row_num,]
      df_el <- data.frame(Prot_acc=protein_acc, Assembly= df_sub$Genome, 
                          Host= df_sub$Aseembly_host, Group=df_sub$Group, MGE='GI', ID=ID)
      MGE_stat_per_protein <- rbind(MGE_stat_per_protein,df_el)
    }
  }
}

strict_filt_MGEs <- unique(MGE_stat_per_protein[MGE_stat_per_protein$Group %in% c('Human_StrictFilt_scoary', 
                                                                      'Insect_StrictFilt_scoary','Plant_StrictFilt_scoary'),]) %>% 
  group_by(MGE, ID, Group, Host,Assembly) %>% summarize(num=n())

write.table(strict_filt_MGEs, 'MGEs_stat_for_Strict.csv', sep='\t', row.names = F)

ggplot(test, aes(x=Group, y=num, fill=Group)) + 
  geom_boxplot()+
  theme_bw()+
  facet_wrap(~MGE, scales='free_y')

func_annot_stat <- read.table('/home/anton/serratia_pangenomics_2022/functional_annotation/Annotations_for_groups_stat.csv',header = TRUE,  sep = "\t")
func_annot_stat$Percent_ident <- round(func_annot_stat$Num_ident/func_annot_stat$Num_genes,1)

extr_vec_all <- c('Accessory', 'Core','Unique', 'Human_Filt_scoary', 'Human_StrictFilt_scoary', 'Insect_Filt_scoary','Insect_StrictFilt_scoary',
                  'Plant_Filt_scoary','Plant_StrictFilt_scoary')

extr_vec_main <- c('Core', 'Accessory','Unique', 'Vir_Core','Vir_Accessory',
                   'Human_StrictFilt_scoary','Plant_StrictFilt_scoary','Insect_StrictFilt_scoary')

func_annot_stat$ident_factor <- factor(ifelse(func_annot_stat$Percent_ident==1.0, 'Ident', 'Non-ident'))
func_annot_stat$ident_factor <- factor(func_annot_stat$ident_factor, levels=c('Non-ident','Ident'))
func_percent_ident_stat_extr <- func_annot_stat[func_annot_stat$Group %in% extr_vec_main, ]
func_percent_ident_stat_extr_for_plot <- func_percent_ident_stat_extr %>% select(Group, Ontology, ident_factor) %>% 
  group_by(Group, Ontology, ident_factor) %>% summarize(freq=n()) %>% as.data.frame()

func_percent_ident_stat_extr_for_plot$Group <- droplevels(func_percent_ident_stat_extr_for_plot$Group)
func_percent_ident_stat_extr_for_plot$Group <- factor(func_percent_ident_stat_extr_for_plot$Group, levels=c( 'Core','Accessory',
              'Unique','Vir_Core','Vir_Accessory', 'Human_StrictFilt_scoary', 'Insect_StrictFilt_scoary','Plant_StrictFilt_scoary'))

levels(func_percent_ident_stat_extr_for_plot$Group) <-c('Core','Accessory', 
                                                        'Unique','Vir_Core','Vir_Accessory', 
                                                        'Human', 'Insect','Plant')

for (group in levels(func_percent_ident_stat_extr_for_plot$Group)) {
  for (ontology in levels(func_percent_ident_stat_extr_for_plot$Ontology)){
    Non_ident_df <- func_percent_ident_stat_extr_for_plot[func_percent_ident_stat_extr_for_plot$Group==group &func_percent_ident_stat_extr_for_plot$Ontology==ontology
                                                            ,]
    if (nrow(Non_ident_df)==1){
      print(Non_ident_df)
      add_df <- data.frame(Group=group, Ontology=ontology, ident_factor='Non-ident', freq=0) 
      func_percent_ident_stat_extr_for_plot <- rbind(func_percent_ident_stat_extr_for_plot, add_df)
    }
  }
}
ident_vect <- func_percent_ident_stat_extr_for_plot[func_percent_ident_stat_extr_for_plot$Group %in% c('Accessory','Core') & func_percent_ident_stat_extr_for_plot$ident_factor=='Ident',]
non_ident_vect <- func_percent_ident_stat_extr_for_plot[func_percent_ident_stat_extr_for_plot$Group %in% c('Accessory','Core') & func_percent_ident_stat_extr_for_plot$ident_factor=='Non-ident',]
mean(ident_vect[,4]/(non_ident_vect[,4]+ident_vect[,4]))

ggplot(func_percent_ident_stat_extr_for_plot, aes(x= ident_factor, y=freq, fill=Group)) + 
  geom_bar(position='dodge', stat='identity', width=0.5)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 16, hjust = 1))+
  ylab('Frequency')+
  xlab("Type of Clusters")+
  facet_wrap(~Ontology) +scale_fill_manual(values = c('#9F08A4', '#D089D3', '#F4D2F5', '#CD1313', '#DA6969',
                                                      '#B17129','#3E1881','#2F690B'))

  #guides(fill= guide_legend(title="Host"))

non_ident_groups <- func_percent_ident_stat_extr[func_percent_ident_stat_extr$Group %in% c('Accessory','Core') & func_percent_ident_stat_extr$ident_factor=='Non-ident', ]
length(unique(non_ident_groups[non_ident_groups$Percent_annot<1,  3]))
length(unique(non_ident_groups[non_ident_groups$Percent_annot<=1,  3]))

non_annot_groups <- func_percent_ident_stat_extr[func_percent_ident_stat_extr$Group %in% c('Accessory','Core') & func_percent_ident_stat_extr$Percent_annot<1 & func_percent_ident_stat_extr$Percent_annot!=0, ]
length(unique(non_annot_groups[,  3]))
print(unique(non_annot_groups[,  3]))


#func_percent_ident_stat_extr$Group %in% c('Core','Accessory'),]
#func_percent_ident_stat_extr[func_percent_ident_stat_extr$Group %in% c('Human_StrictFilt_scoary', 'Insect_StrictFilt_scoary', 'Plant_StrictFilt_scoary'),]
#Num_terms, Num_annot
p <- ggplot(func_percent_ident_stat_extr[func_percent_ident_stat_extr$Group %in% c('Core','Accessory'),], aes(x=Num_genomes, y=Num_terms, fill = Percent_annot, size=Num_term_groups)) +
  geom_point(stat='identity',alpha=0.6, shape=21, color='black')+
  theme_bw()+
  xlab('Number of genomes') +
  #ylab('Number of annotations')+
  ylab('Number of terms')+
  theme( axis.text.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14),
         strip.text = element_text(size = 15))+
  scale_fill_gradientn(colours=c("white", "darkgreen"))+
  facet_wrap(~Ontology, scales='free')

p$labels$fill <- "Percent of annotated genes" 
p$labels$size <- "Number of annotation terms" 
#1006 - width for TIFF pic
#12/7
#13/8 for terms, 1206

Mean_anot <- func_annot_stat[func_annot_stat$Group %in% extr_vec_main,] %>% select(Group, Ontology, Percent_annot) %>% 
  group_by(Group, Ontology) %>% summarize(mean_annot=mean(Percent_annot)) %>% as.data.frame()

Mean_anot$Group <- droplevels(Mean_anot$Group)
Mean_anot$Group <- factor(Mean_anot$Group, levels=c('Core','Accessory', 
                                                     'Unique','Vir_Core','Vir_Accessory', 'Human_StrictFilt_scoary', 'Insect_StrictFilt_scoary',
                                                    'Plant_StrictFilt_scoary'))

levels(Mean_anot$Group) <-c('Core','Accessory', 
                           'Unique','Vir_Core','Vir_Accessory', 
                            'Human', 'Insect','Plant')

ggplot(Mean_anot, aes(x= Group,y= mean_annot, fill=Group)) +
  geom_bar(position='dodge', stat='identity', width=0.5)+
  theme_bw()+
  theme(axis.title.x=element_text(face="bold", color="black", 
                                  size=20),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(color='black', 
                                   size=16),
        axis.title.y = element_text(face="bold", color="black", 
                                    size=18),
        legend.title=element_text(face="bold",size=14), 
        legend.text=element_text(size=12),
        axis.text.x = element_text(color='black', 
                                   angle = 60, vjust = 1, 
                                   size = 10, hjust = 1))+
  ylab('Mean annotation percent')+
  xlab("Group") +scale_fill_manual(values = c('#9F08A4', '#D089D3', '#F4D2F5', '#CD1313', '#DA6969',
                                               '#B17129','#3E1881','#2F690B'))+
  facet_wrap(~Ontology)
#12/7

GSEA_results_all <- read.table('/home/anton/serratia_pangenomics_2022/functional_annotation/GSEA_results_COG_KEGG_reference_set.csv',
                               header = TRUE,  sep = "\t", , quote = "")

GSEA_results_all$x_1 <-GSEA_results_all$x_1+1
no_annot <- GSEA_results_all[GSEA_results_all$Group=='Non_annot',]

main_pic_vec <- c('Core', 'Accessory','Unique', 'Vir_Core','Vir_Accessory',
                  'Human_StrictFilt_scoary','Plant_StrictFilt_scoary','Insect_StrictFilt_scoary') 
GSEA_results_all <- GSEA_results_all[GSEA_results_all$Group %in% main_pic_vec,]

GSEA_results_all$Group <- droplevels(GSEA_results_all$Group)
GSEA_results_all$Group <- factor(GSEA_results_all$Group, levels=c( 'Core','Accessory',
                     'Unique','Vir_Core','Vir_Accessory', 'Human_StrictFilt_scoary', 'Insect_StrictFilt_scoary','Plant_StrictFilt_scoary'))

levels(GSEA_results_all$Group) <-c('Core','Accessory', 
                                     'Unique','Vir_Core','Vir_Accessory', 
                                      'Human', 'Insect','Plant')

Top_gene_types <- GSEA_results_all[with(GSEA_results_all, order(- x_1, pval)), ] %>% 
  group_by(Group, Ontology) %>% dplyr::slice(1:10) %>% as.data.frame()

Top_gene_types_no_annot <- no_annot[with(no_annot, order(- x_1, pval)), ] %>% 
  group_by(Group, Ontology) %>% dplyr::slice(1:10) %>% as.data.frame()


#darkgreen - COG_gene_types_enrichments
#darkred - GO_biological_process_gene_types_enrichments
#darkorange - GO_cellular_component_gene_types_enrichments
# #BBB227 - GO_molecular_function_gene_types_enrichments
#darkblue - KEGG_gene_types_enrichments

#Top_gene_types_no_annot  Top_gene_types
ggplot(Top_gene_types[ Top_gene_types$Ontology=='COG',], 
       aes(x=Group, fill = pval, y=Annotation, size=x_1/n)) +
  geom_point(shape=21, color='black')+
  theme_bw()+
  xlab('Group') +
  ylab('Annotation')+
  scale_size_continuous(range = c(1.5,6))+
  theme( axis.text.x = element_text(color='black', 
                                    angle = 65, vjust = 1, 
                                    size = 12, hjust = 1),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         #panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14,hjust = 0.8),
         strip.text = element_text(size = 18))+
  guides(size = guide_legend(title="Enrichment ratio"))+
  scale_fill_gradientn(colours=c( "darkgreen", "white")) + #darkgreen, darkred, darkorange, #BBB227, darkblue
  coord_fixed(ratio = 1.5) #ratio=1.5 - ref - non_annot 1.4 - default 0.65 -KEGG
#10/6 -COG non_annot - default #10/8 -KEGG

#/home/anton/serratia_pangenomics_2022/functional_annotation/dissimilarity_matrix/All_jac.csv
#/home/anton/serratia_pangenomics_2022/functional_annotation/dissimilarity_matrix/All_sim.csv
#/home/anton/serratia_pangenomics_2022/functional_annotation/dissimilarity_matrix/Extr_jac.csv
#/home/anton/serratia_pangenomics_2022/functional_annotation/dissimilarity_matrix/Extr_sim.csv

distance_matrix <- read.table('/home/anton/serratia_pangenomics_2022/functional_annotation/dissimilarity_matrix/Extr_sim.csv', row.names=1,header = TRUE,  sep = "\t")

sc_dist<-data.frame(t(distance_matrix))

wss <- (nrow(sc_dist)-1)*sum(apply(sc_dist,2,var))
for (i in 1:5) wss[i] <- sum(kmeans(sc_dist,
                                    centers=i)$withinss,nstart=25,iter.max=1000)
wss_gf=data.frame(clust=c(1:5), wss=wss)
ggplot(wss_gf, aes(clust,wss))+
  geom_line()+theme_bw()+ xlab('Number of clusters') + ylab('WSS')+
  geom_vline(xintercept =4, linetype="dashed")+
  theme(axis.text.y = element_text(color='black', 
                                   size=12),
        axis.title.y=element_text(color="black", 
                                  size=16),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=12),
        axis.title.x = element_text( color="black", 
                                     size=16)
  ) 

dist_kmeans <- kmeans(sc_dist,centers=4)
dist_kmeans$cluster
sc_dist$name=rownames(sc_dist)

sc_dist$clust= as.factor(dist_kmeans$cluster)

autoplot(dist_kmeans, 
         data=sc_dist, label = TRUE,frame = TRUE,
         frame.type = 'norm',label.size = 3,label.col='black',alpha=0)+
  scale_colour_manual(values = c('black','black','black','black'))+
  scale_fill_manual(values = c("#A9A9A9","#A60B0B",'blue','orange'))+
  theme_bw()+
  theme(axis.text.y = element_text(color='black', 
                                   size=24),
        axis.title.y=element_text(color="black", 
                                  size=28),
        panel.background = element_blank(), 
        axis.line = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(color='black', 
                                   size=24),
        axis.title.x = element_text(color="black", 
                                    size=28),
        legend.position = 'none'
  ) + 
  scale_x_continuous(expand = c(.1, .1)) 
scale_y_continuous(expand = c(.1, .1))

#10/6
func_terms_for_GO <- read.table('/home/anton/serratia_pangenomics_2022/functional_annotation/Annotation_results_per_ontology_cleaned.csv',header = TRUE,  sep = "\t", quote = "", stringsAsFactors = FALSE)
func_terms_GO_mearged <- data.frame(Accession = as.character(func_terms_for_GO[,1]))
func_terms_GO_mearged$GO <- paste(func_terms_for_GO$GO.biological_process,
                                  func_terms_for_GO$GO.molecular_function,
                                  func_terms_for_GO$GO.cellular_component, sep=';')

func_terms_GO_mearged$GO <- gsub(";-;", "", func_terms_GO_mearged$GO)
func_terms_GO_mearged$GO <- gsub(";-", "", func_terms_GO_mearged$GO)
func_terms_GO_mearged$GO <- gsub("-;", "", func_terms_GO_mearged$GO)
func_terms_GO_mearged$GO <- gsub("-", "", func_terms_GO_mearged$GO)
func_terms_GO_mearged$Accession <- as.character(func_terms_GO_mearged$Accession)

GO_list <- list()

#for (ind in 1:nrow(func_terms_GO_mearged)){
for (ind in 1:nrow(func_terms_GO_mearged)){
  GO_vec=c()
  GO_split <- strsplit(func_terms_GO_mearged[ind,2], ";")[[1]]
  if (length(GO_split)==0){
    GO_vec <- c(GO_vec, '')
  } else {
    for (GO in GO_split){
      GO_vec <- c(GO_vec, GO)
    }
  }
  
  GO_list[[func_terms_GO_mearged[ind,1]]] <- GO_vec
} 
  

Groups_df <- read.table('/home/anton/serratia_pangenomics_2022/Groups_all_for_GO.csv',header = TRUE,  sep = "\t", quote = "", stringsAsFactors = FALSE)


ontology_df <- data.frame(ontology=c('GO:biological_process','GO:cellular_component','GO:molecular_function'),
                          code= c('BP','CC','MF'), num_terms = c(2600, 250, 1800))

enrich_GO_df <- data.frame(GO.ID=c(), Term=c(),Annotated=c(),Significant=c(), Expected=c(), classicFisher=c(),
                           group=c(), ontology=c())

for (group in levels(as.factor(Groups_df$Group))){
  Genes_intr_vec <- Groups_df[Groups_df$Group==group,2]
  geneNames_All <- names(GO_list)
  geneList_GO  <- factor(as.integer(geneNames_All  %in% Genes_intr_vec ))
  names(geneList_GO ) <- geneNames_All 
  for (ontology in c('BP','CC','MF')){
    print(c(group, ontology))
    GOdata  <- new("topGOdata", ontology = ontology, allGenes = geneList_GO ,
                   annot = annFUN.gene2GO, gene2GO = GO_list)
    
    resultFisher_GO  <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
    top_nodes <- ontology_df[ontology_df$code==ontology,3]
    allRes_GO  <- GenTable(GOdata , classicFisher = resultFisher_GO, topNodes = top_nodes)
    ontology_name <- as.character(ontology_df[ontology_df$code==ontology,1])
    allRes_GO$group <- group
    allRes_GO$ontology <- ontology_name
    allRes_GO[allRes_GO$classicFisher=='<1e-30', 6] <- 1e-30
    allRes_GO[allRes_GO$classicFisher=='< 1e-30',6] <- 1e-30
    allRes_GO$classicFisher <- as.numeric(allRes_GO$classicFisher)
    allRes_GO <- allRes_GO[allRes_GO$classicFisher<=0.05,]
    enrich_GO_df <- rbind(enrich_GO_df, allRes_GO)
  }
}

write.table(enrich_GO_df, 'GO_topGO_enrichments_main.csv', sep='\t', row.names = F)




GO_results_all <-  read.table('/home/anton/serratia_pangenomics_2022/functional_annotation/GO_topGO_enrichments_main.csv',header = TRUE,  sep = "\t",stringsAsFactors = FALSE)

GO_results_all$group <- factor(GO_results_all$group)
GO_results_all <- GO_results_all[GO_results_all$GO.ID!='GO:0015050']


main_pic_vec <- c('Core', 'Accessory','Unique', 'Vir_Core','Vir_Accessory',
                  'Human_StrictFilt_scoary','Plant_StrictFilt_scoary','Insect_StrictFilt_scoary') 
GO_results_all <- GO_results_all[GO_results_all$group %in% main_pic_vec,]


GO_results_all$group <- droplevels(GO_results_all$group)
GO_results_all$group <- factor(GO_results_all$group, levels=c( 'Core','Accessory',
                                                                   'Unique','Vir_Core','Vir_Accessory', 'Human_StrictFilt_scoary', 'Insect_StrictFilt_scoary','Plant_StrictFilt_scoary'))

levels(GO_results_all$group) <-c('Core','Accessory', 
                                   'Unique','Vir_Core','Vir_Accessory', 
                                   'Human', 'Insect','Plant')

Top_enrich_sig <- GO_results_all[with(GO_results_all, order(-Significant, classicFisher)), ] %>% 
  group_by(group, ontology) %>% dplyr::slice(1:7) %>% as.data.frame()

Top_gene_types_no_annot <- no_annot[with(no_annot, order(-Significant, classicFisher)), ] %>% 
  group_by(group, ontology) %>% dplyr::slice(1:10) %>% as.data.frame()

#darkred - GO:biological_process
#darkorange - GO:cellular_component
# #BBB227 - GO:molecular_function

ggplot(Top_gene_types_no_annot[ Top_gene_types_no_annot$ontology=='GO:biological_process',], 
       aes(x=group, fill = classicFisher, y= Term, size=Significant/Annotated)) +
  geom_point(shape=21, color='black')+
  theme_bw()+
  xlab('Group') +
  ylab('Annotation')+
  scale_size_continuous(range = c(1.5,6))+
  theme( axis.text.x = element_text(color='black', 
                                    angle = 65, vjust = 1, 
                                    size = 12, hjust = 1),
         axis.title.x=element_text(face="bold", color="black", 
                                   size=14),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         #panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank(),
         axis.text.y = element_text(color='black', 
                                    size=12),
         axis.title.y = element_text(face="bold", color="black", 
                                     size=14,hjust = 0.8),
         strip.text = element_text(size = 18))+
  guides(size = guide_legend(title="Enrichment ratio"))+
  scale_fill_gradientn(colours=c( "darkred", "white"))  +#darkgreen, darkred, darkorange, #BBB227, darkblue
  coord_fixed(ratio = 1.3) #ratio=1.3 - MF 
#non_annot - default

#12/9 10/7 -CC
