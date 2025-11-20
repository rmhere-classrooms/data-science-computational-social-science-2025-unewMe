library(igraph)

# Graf wedle modelu Barabási-Albert z tysiącem węzłów
g <- barabasi.game(1000)

# Wizualizacja grafu layoutem Fruchterman & Reingold
layout <- layout.fruchterman.reingold(g)
plot(g, layout=layout, vertex.size=2,
     vertex.label=NA, edge.arrow.size=.2)

# Znajdź najbardziej centralny węzeł według miary betweenness, jaki ma numer? -> w tym uruchomieniu 21
b <- betweenness(g)
b
which.max(b)

# Jaka jest średnica grafu? W tym uruchomieniu -> 10
diameter(g)

# 6. 
# Erdős-Rényi ma losowe krawędzie między parami węzłów (stopnie ~ rozkład Poissona, brak wyraźnych hubów), 
# a Barabási–Albert dodaje nowe węzły „preferencyjnie” do już dobrze podłączonych (prawo potęgowe stopni, powstają huby i sieć scale-free).


