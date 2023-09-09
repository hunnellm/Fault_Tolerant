def non_adj_twins(g):
    twins=[]
    for v in g.vertices():
        for w in g.vertices():
            if w>v:
                if g.neighbors(v)==g.neighbors(w): #neighborhoods are open by default
                        twins.append([v,w])
    return(twins)
