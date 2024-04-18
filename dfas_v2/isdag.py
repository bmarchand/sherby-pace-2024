import networkx as nx

G = nx.read_weighted_edgelist("1.sol.gr", create_using=nx.DiGraph)

print(nx.is_directed_acyclic_graph(G))