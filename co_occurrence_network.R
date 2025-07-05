library(igraph)
library(Hmisc)
library(ggplot2)
library(ggpubr)
library(ggsignif)

ASV_table= read.table("C:/Users/RM_net.txt",fill=T,check.names =F,row.names =1,header=T, sep="\t")
ASV_table <- t(ASV_table)
tax_table= read.table("C:/Users/final/formatted_taxonomy.txt",fill=T,check.names =F,row.names =1,header=T, sep="\t")

prev_prop_df = apply(X= ASV_table,
                     MARGIN = 2,
                     FUN = function(x){sum(x>0)/dim(ASV_table)[1]})

total_relative_abundance = colSums(ASV_table)/sum(colSums(ASV_table))
is.above.prev = prev_prop_df >0.2
table(is.above.prev)

is.above.abundance =total_relative_abundance >0.0001
table(is.above.abundance)
table(is.above.prev & is.above.abundance)
target_ASV = colnames(ASV_table)[is.above.prev & is.above.abundance]

ASV_table.prop = as.data.frame(t(apply(ASV_table, 1, function(ASV) ASV/sum(ASV))))
target_ASV_table.prop = ASV_table.prop[, target_ASV]
print(rowSums(target_ASV_table.prop))

correlation_matrix<-rcorr(as.matrix(target_ASV_table.prop),type="spearman")
                          
correlation_threshold<-0.7
p_threshold<-0.001

CorrDF <- function(cormat,pmat){
  ut<- upper.tri(cormat)
  data.frame(
    source = rownames(cormat)[col(cormat)[ut]],
    target = rownames(cormat)[row(cormat)[ut]],
    corr =(cormat)[ut],
    weight = abs((cormat)[ut]),
    p = pmat[ut],
    p.adjust = p.adjust(pmat[ut], method="BH"),
    cor_type = ifelse((cormat)[ut]>0,"Positive","Negative"))}

cor_df<-CorrDF(correlation_matrix$r , correlation_matrix$P)

significant_cor_df <- cor_df[which(cor_df$weight >= correlation_threshold),]
significant_cor_df <- cor_df[which(cor_df$p.adjust <p_threshold),] 
  
g = graph_from_data_frame(significant_cor_df, directed=F)

significant_tax= tax_table[V(g)$name,]
V(g)$Domain= significant_tax$Domain
V(g)$Phylum=significant_tax$Phylum
V(g)$Class=significant_tax$Class
V(g)$Order =significant_tax$Order
V(g)$Family= significant_tax$Family
V(g)$Genus=significant_tax$Genus
V(g)$Species =significant_tax$Species

node_list = data.frame(node_id = names(V(g)),significant_tax)
edge_list = data.frame(edge_id = paste0("edge_", c(1:length(E(g)))),significant_cor_df)

table(edge_list$cor_type)

nodes_num = length(V(g))
nodes_num
edges_num = length(E(g))
edges_num
positive.corr_num =sum(E(g)$corr>0)
positive.corr_num

positive.corr_prop = positive.corr_num/edges_num 
positive.corr_prop
negative.corr_num=sum(E(g)$corr<0)
negative.corr_num

negative.corr_prop = negative.corr_num/edges_num 
negative.corr_prop
average_degree = mean(degree(g))
average_degree
#平均度
average_shortest_path_length = mean_distance(g, directed = FALSE)
average_shortest_path_length

network_diameter = diameter(g,directed =FALSE)
network_diameter
network_density = edge_density(g)
network_density

modularity = modularity(g, membership(cluster_louvain(g)),directed = F)
modularity

clustering_coefficient<-transitivity(g, type="global")
clustering_coefficient
network_parameter =data.frame(nodes_num,
                              edges_num,
                              positive.corr_num,
                              positive.corr_prop,
                              negative.corr_num,
                              negative.corr_prop,
                              average_degree,
                              average_shortest_path_length,
                              network_diameter,
                              network_density,
                              clustering_coefficient,
                              modularity)

output_dir = "C:/Users"
write.table(node_list, file = file.path(output_dir, "network.node_list.tsv"),sep="\t",quote =F,row.names =F)
write.table(edge_list, file = file.path(output_dir, "network.edge_list.tsv"),sep="\t",quote =F,row.names =F)
write.table(network_parameter, file = file.path(output_dir, "network_parameter.tsv"),sep="\t",quote =F,row.names =F)
write_graph(g,file.path(output_dir, 'network.graphml'),format = 'graphml')

