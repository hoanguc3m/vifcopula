setwd("/home/hoanguc3m/Dropbox/WP2")
library("igraph")

pdf("img/onefcop.pdf", width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(3,3,3,1))

g <- graph.empty (6, directed = FALSE)
new_edges <- c(1,2, 1,3, 1,4, 1,5, 1,6)
g <- add.edges(g, new_edges)
V(g)$number <- c(1,2,3,4,5,6)
V(g)[ number == 1 ]$color <- "orange"
V(g)[ number > 1 ]$color <- "green"
E(g)$width <- 5
V(g)$name <- c(expression(v[1]), expression(u[1]), expression(u[2]),
               expression(u[3]), expression(u[4]), expression(v[5]))
V(g)$size      <- rep(50,6)
#vertex_attr(g)
plot(g, main = "Factor level 1")

g <- graph.empty (6, directed = FALSE)
new_edges <- c(1,2, 1,3, 1,4, 1,5, 1,6)
g <- add.edges(g, new_edges)
V(g)$number <- c(1,2,3,4,5,6)
V(g)[ number == 1 ]$color <- "orange3"
V(g)[ number > 1 ]$color <- "green"
E(g)$width <- 3
V(g)$name <- c(expression(v[2]), expression(u[1|v[1]]), expression(u[2|v[1]]),
               expression(u[3|v[1]]), expression(u[4|v[1]]), expression(u[5|v[1]]))
V(g)$size      <- rep(40,6)
plot(g, main = "Factor level 2")
dev.off()



pdf("img/bifcop.pdf", width = 9, height = 4.5)
par(mfrow =c(1,2))
par(mar=c(3,3,3,1))

g <- graph.empty (6, directed = FALSE)
new_edges <- c(1,2, 1,3, 1,4, 1,5, 1,6)
g <- add.edges(g, new_edges)
V(g)$number <- c(1,2,3,4,5,6)
V(g)[ number == 1 ]$color <- "orange"
V(g)[ number > 1 ]$color <- "green"
E(g)$width <- 5
V(g)$name <- c(expression(v[1]), expression(u[1]), expression(u[2]),
               expression(u[3]), expression(u[4]), expression(v[5]))
V(g)$size      <- rep(50,6)
#vertex_attr(g)
plot(g, main = "Bifactor level 1")

g <- graph.empty (7, directed = FALSE)
new_edges <- c(1,2, 1,3, 1,4, 7,5, 7,6)
g <- add.edges(g, new_edges)
V(g)$number <- c(1,2,3,4,5,6,7)
V(g)[ number > 1 ]$color <- "green"
V(g)[ number == 1 ]$color <- "orange3"
V(g)[ number == 7 ]$color <- "orange3"
E(g)$width <- 3
V(g)$name <- c(expression(v[21]), expression(u[1|v[1]]), expression(u[2|v[1]]),
               expression(u[3|v[1]]), expression(u[4|v[1]]), expression(u[5|v[1]]),
               expression(v[22]))
V(g)$size      <- rep(40,7)
plot(g, main = "Bifactor level 2")
dev.off()

pdf("img/nested.pdf", width = 4.5, height = 4.5)
par(mfrow =c(1,1))
par(mar=c(3,3,3,1))

g <- graph.empty (8, directed = FALSE)
new_edges <- c(1,2, 1,3, 2,4, 2,5, 2,6, 3,7,3,8)
g <- add.edges(g, new_edges)
V(g)$number <- c(1,2,3,4,5,6,7,8)
V(g)[ number == 1 ]$color <- "orange"
V(g)[ number > 1 ]$color <- "yellow"
V(g)[ number > 3 ]$color <- "green"
E(g)$width <- 5
V(g)$name <- c(expression(v[1]), expression(v[2]), expression(v[3]),
               expression(u[1]), expression(u[2]),
               expression(u[3]), expression(u[4]), expression(v[5]))
V(g)$size      <- rep(30,8)
#vertex_attr(g)
plot(g, main = "Nested copula")

dev.off()
