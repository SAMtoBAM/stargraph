##this scripts is a skeleton for plotting networks created by stargraph

options(warn = -1)

library(igraph)
library(ggraph)
library(tidygraph)
library(tidyverse)
library(ggforce)
#library(concaveman)
library(ggpubr)


####generate network plots using the pairwise weights (jaccard similarities based on sourmash kmer signatures)
edges_jac = read.csv(file="PATHTOOUTPUT/PREFIX.starships_SLRs.pairwise_jaccard.tsv", header=T, sep='\t')
nodes = read.csv(file="PATHTOOUTPUT/PREFIX.starships_SLRs.metadata.tsv", header=T, sep='\t')

nodes = nodes %>% distinct(name, .keep_all = TRUE)

graph_jac = tbl_graph(nodes = nodes, edges = edges_jac, directed = FALSE)

##plotting all the Starship and SLR elements
g1=ggraph(graph_jac, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), edge_colour = "black", show.legend = FALSE) +
  geom_node_point(size = 3, aes(colour = type)) +
  theme_void() +
  theme(legend.position = "top",plot.margin = unit(rep(1, 4), "cm"),legend.title=element_blank())+
  scale_color_manual(values=c("blue", "red"))


##one potential measure of false positives would the number of elements without any connections
##we could then compare the number of SLRs and Starships with no connections
##in terms of both total number and proportion

##create data frame with connectivity measured as 'degree' with 0 meaning no connections
##then classify those with 0 as unconnected in 'connection_status'
connectivity_jac=as.data.frame(graph_jac %>%
                             activate(nodes) %>%
                             mutate(degree = centrality_degree(), connection_status = if_else(degree > 0, "connected", "unconnected")))

##now we can plot the number of each groups
connectivity_prop_jac= connectivity_jac %>%
  group_by(type, connection_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  filter(connection_status == "unconnected")

g2=ggplot(connectivity_prop_jac)+geom_col(aes(x=type, y=proportion*100, fill =type))+
  theme_pubr(legend = "none", x.text.angle = 65)+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Percentage of total nodes without connections") +
  xlab(element_blank())

g3=ggplot(connectivity_prop_jac)+
  geom_col(aes(x=type, y=n, fill =type))+
  theme_pubr(legend = "none", x.text.angle = 65)+
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Number of total nodes without connections") +
  xlab(element_blank())+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))


##plot the network plus the counts of unconnected
plot_jac=ggarrange(g1,g2,g3, ncol = 3, widths = c(10,1,1))

suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/network_plotting.jaccard.png", plot = plot_jac, units = "in", height = 300, width = 400, limitsize = FALSE)))
suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/network_plotting.jaccard.svg", plot = plot_jac, units = "in", height = 300, width = 400, limitsize = FALSE)))


################################################################
###now do the same but for the containment values
edges_cont = read.csv(file="PATHTOOUTPUT/PREFIX.starships_SLRs.pairwise_containment.tsv", header=T, sep='\t')
#nodes = read.csv(file="PATHTOOUTPUT/PREFIX.starships_SLRs.metadata.tsv", header=T, sep='\t')

#nodes = nodes %>% distinct(name, .keep_all = TRUE)

graph_cont = tbl_graph(nodes = nodes, edges = edges_cont, directed = FALSE)

##plotting all the Starship and SLR elements
g4=ggraph(graph_cont, layout = "fr") +
  geom_edge_link(aes(edge_alpha = weight), edge_colour = "black", show.legend = FALSE) +
  geom_node_point(size = 3, aes(colour = type)) +
  theme_void() +
  theme(legend.position = "top",plot.margin = unit(rep(1, 4), "cm"),legend.title=element_blank())+
  scale_color_manual(values=c("blue", "red"))


##one potential measure of false positives would the number of elements without any connections
##we could then compare the number of SLRs and Starships with no connections
##in terms of both total number and proportion

##create data frame with connectivity measured as 'degree' with 0 meaning no connections
##then classify those with 0 as unconnected in 'connection_status'
connectivity_cont=as.data.frame(graph_cont %>%
                             activate(nodes) %>%
                             mutate(degree = centrality_degree(), connection_status = if_else(degree > 0, "connected", "unconnected")))

##now we can plot the number of each groups
connectivity_prop_cont= connectivity_cont %>%
  group_by(type, connection_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(
    total = sum(n),
    proportion = n / total
  ) %>%
  filter(connection_status == "unconnected")

g5=ggplot(connectivity_prop_cont)+geom_col(aes(x=type, y=proportion*100, fill =type))+
  theme_pubr(legend = "none", x.text.angle = 65)+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Percentage of total nodes without connections") +
  xlab(element_blank())

g6=ggplot(connectivity_prop_cont)+
  geom_col(aes(x=type, y=n, fill =type))+
  theme_pubr(legend = "none", x.text.angle = 65)+
  scale_fill_manual(values=c("blue", "red")) +
  ylab("Number of total nodes without connections") +
  xlab(element_blank())+
  scale_y_continuous(labels = scales::label_number(accuracy = 1))


##plot the network plus the counts of unconnected
plot_cont=ggarrange(g4,g5,g6, ncol = 3, widths = c(10,1,1))


suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/network_plotting.containment.png", plot = plot_cont, units = "in", height = 300, width = 400, limitsize = FALSE)))
suppressMessages(suppressWarnings(ggsave("PATHTOOUTPUT/network_plotting.containment.svg", plot = plot_cont, units = "in", height = 300, width = 400, limitsize = FALSE)))

