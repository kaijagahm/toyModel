# Testing out tidygraph
library(tidygraph)
library(igraph) # tidygraph can coexist happily with igraph!
library(ggraph)
library(tidyverse)

a <- create_ring(10)
plot(a) # plots it; it looks just like an igraph object.

# Can create deterministic graphs using `create` functions, and can create probabilistic graphs using `play` functions.

# set some parameters
nnodes = 40 # how many individuals to include
dens = 0.15 # density of edges

# Create a graph
g <- play_erdos_renyi(n = nnodes, p = dens, loops = FALSE, directed = FALSE) %>%
  mutate(name = 1:nnodes,
         deg = degree(.),
         Sex = sample(c("F", "M"), size = nnodes, replace = T),
         Age = rnorm(n = nnodes, mean = 30, sd = 10),
         degCen = centrality_degree())

# Plot the graph
g %>%
  ggraph()+
  geom_edge_link(edge_width = 0.1)+
  geom_node_point(aes(col = Sex, size = Age))+
  geom_node_text(aes(label = name, size = Age*0.5), 
                 colour = 'white', vjust = 0.5) + 
  theme_graph()

