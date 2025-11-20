library(igraph)

# Sieć Erdős-Rényi o stu wierzchołkach i prawdopodobieństwie krawędzi = 0.05.
g <- erdos.renyi.game(p.or.m=0.05, n=100)

# Podsumowanie grafu
summary(g)
#IGRAPH 6c1abfc U--- 100 232 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n)

# Czy graf jest ważony? - Nie, graf nie jest ważony, żeby był ważony musielibyśmy mieć U-W- --> wagi przypiszemy w poźniejszym kroku

# Wylistowanie wierzchołków
V(g)

# Wylistowanie krawędzi 
E(g)

# Ustawienie wag wszystkich krawędzi na losowe z zakresu 0.01 do 1
E(g)$weight <- runif(length(E(g)), 0.01, 1)

# Podsumowanie grafu
summary(g)

# IGRAPH 6c1abfc U-W- 100 232 -- Erdos-Renyi (gnp) graph
# + attr: name (g/c), type (g/c), loops (g/l), p (g/n), weight (e/n)

# Czy graf teraz jest ważony? - Tak, teraz graf jest grafem ważonym, mamy U-W- W oznacza weight

# Stopnie każdego węzła
degree(g)

# Histogram stopni węzłów
hist(degree(g))

# Connected components w grafie
cl <- clusters(g)
cl

#$membership
# [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
# [63] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1

# $csize
# [1] 100

# $no
# [1] 1

# Liczba connected components to w moim uruchomieniu 1, aby zapewnić determistyczność można by ustawić seed. 

# PageRank
pr <- page.rank(g)$vector

plot(g, vertex.size=pr*300,
     vertex.label=NA, edge.arrow.size=.2)


