library(igraph)
source("local_hipathia.R")
# other example of fusion protein
# Q03164, P51825, or, P42568, or, Q03111
# A & (B |C | D) OR [(A&B) |C | D]
# ->[(A&B) & (C | D) =  U1 / U2 ] => A,/,B,/,C,D
########################## test
# Let's N -> F
# case 1: (P52948, and, Q96L73, or, P29375, or, Q13206, or, P31269) 
# -> (A & B) & (C | D |E) =>  0=min (min(A,B) , max(C,D,E)) X
# o-> X -> 0
sif1 <- data.frame("N-hsapath1-N","activation","N-hsapath1-FF")
att1 <- data.frame(ID=c(sif1[,1],sif1[,3]), 
                   label= c("N","FF"),
                   genesList = c("X,/,Y,/,Z,W,O","FF"))
# Case2: o -> A & B -> 0   => 0= P(min(A,B), C,D, E) 
# o -> C -> 0
# o -> D -> 0
# o -> E -> 0
sif2 <- data.frame(c("N-hsapath2-AB","N-hsapath2-C", "N-hsapath2-D","N-hsapath2-E"), rep("activation", 4), rep("N-hsapath2-F",4))
att2 <- data.frame(ID=paste0("N-hsapath2-",c("AB","C", "D","E","F")),
                   label=c("AB","C", "D","E","F"),
                   genesList = c("A,/,B","C","D","E","F"))
ig1<- creatGraph(sif1,att1)
ig2<- creatGraph(sif2, att2)
pathigraphs <- list()
pathigraphs$hsapath1<-list()
pathigraphs$hsapath2<-list()
pathigraphs <- completGraph(ig1,"hsapath1",pathigraphs)
pathigraphs <- completGraph(ig2,"hsapath2",pathigraphs)
metaginfo <- local_create_metaginfo_object(pathigraphs, "hsa", by.user = TRUE)

visNetwork::visIgraph(metaginfo$pathigraphs$hsapath1$graph)
visNetwork::visIgraph(metaginfo$pathigraphs$hsapath2$graph)
metaginfo$eff.norm
# Drug  afect (A, B, C , D) 
A<- 0.8
B<- 0.8
C<-0.8

probabilities<- rbind(A,B,C)
# if plain node
node_val<- apply(probabilities, 2, stats::quantile, 0.9, 
      na.rm = TRUE)
quantile(probabilities)
# if complex then 
node_val<- colMins(probabilities_mat, 
                                   na.rm = TRUE)

# coputing
s1 <- prettyifelse(nactivators > 0, hipathia:::colProds(1 - 
                                               activator_signals, na.rm = TRUE), rep(0, length(node_val)))
s2 <- prettyifelse(ninhibitors > 0, apply(apply(inhibitor_signals, 
                                                1, max) + apply(inhibitor_signals, 1, min) - 
                                            inhibitor_signals, 2, prod), rep(1, length(node_val)))
signal <- (1 - s1) * s2