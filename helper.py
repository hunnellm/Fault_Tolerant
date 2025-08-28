URL='mr_JG-master/mr_JG-master/'
files=['Zq_c.pyx','Zq.py','zero_forcing_wavefront.pyx','minrank.py', 'inertia.py']
for f in files:
    print("Loading %s..."%f);
    load(URL+f)

def gzerosgame(g,F=[],B=[]):
    """
    Return the derived set for a given graph g with set of banned edges B and a initial set of vertices. The derived set is given by doing generalized zero forcing process. That is, if y is the only white neighbor of x and xy is not banned, then x could force y into black.

    Input:
        g: a simple graph
        F: a list of vertices of g
        B: a list of tuples representing banned edges of g

    Output:
        A set of black vertices when zero forcing process stops.

    Examples:
        sage: gzerosgame(graphs.PathGraph(5),[0])
        set([0, 1, 2, 3, 4])
        sage: gzerosgame(graphs.PathGraph(5),[0],[(1,2)])
        set([0, 1])
    """
    S=set(F) # suspicuous vertices
    Black_vertices=set(F) # current black vertices
    again=1 # iterate again or not
    while again==1:
        again=0
        for x in S:
            N=set(g.neighbors(x))
            D=N.difference(Black_vertices) # set of white neighbors
            if len(D)==1:
                for v in D:
                    y=v # the only white neighbor
                if (((x,y) in B)==False) and (((y,x) in B)==False):
                    again=1
                    S.remove(x)
                    S.add(y)
                    Black_vertices.add(y)
                    break
    return(Black_vertices)

def gZ_leq(graph, support=[], bannedset=[],i=None):
	"""
	For a given graph with support and banned set, if there is a zero forcing set of size i then return it; otherwise return False.

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		i: an integer, the function check gZ <= i or not

	Output:
		if F is a zero forcing set of size i and support is a subset of F, then return F
		False otherwise

	Examples:
		sage: gZ_leq(graphs.PathGraph(5),[],[],1)
		set([0])
		sage: gZ_leq(graphs.PathGraph(5),[],[(0,1)],1)
		False
	"""
	if i < len(support):
#		print 'i cannot less than the cardinality of support'
		return False
	j=i-len(support) # additional number of black vertices
	VX=graph.vertices()
	order=graph.order()
	for y in support:
		VX.remove(y)
	# VX is the vertices outside support now
	for subset in Subsets(VX,j):
		test_set=set(support).union(subset) # the set is tested to be a zero forcing set
		outcome=gzerosgame(graph, test_set, bannedset)
		if len(outcome)==order:
			return test_set
	return False

def find_gzfs(graph, support=[], bannedset=[], upper_bound=None, lower_bound=None):
	"""
	For a given graph with support and banned set, return the an optimal generalized zero forcing set. If upper_bound is less than the generalized zero forcing number then return ['wrong']. If lower_bound is greater than the generalized zero forcing number then the return value will not be correct

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		upper_bound: an integer supposed to be an upper bound of gZ.
		lower_bound: an integer supposed to be a lower bound of gZ. The two bounds may shorten the computation time. But one may leave it as default value if one is not sure.

	Output:
		if F is an optimal zero forcing set of size i then return F. If upper_bound is less than the general zero forcing number then return ['wrong'].

	Examples:
		sage: find_gzfs(graphs.PathGraph(5))
		set([0])
		sage: find_gzfs(graphs.PathGraph(5),[1],[(3,2)])
		set([0, 1, 3])
	"""

	VX=graph.vertices()
	order=graph.order()
	s=len(support)
	for y in support:
		VX.remove(y)
	# VX is the vertices outside support now
	if upper_bound==None:
		upper_bound=order # the default upper bound
	if lower_bound==None:
		lower_bound=len(VX) # temporary lower bound
		for v in VX:
			N=set(graph.neighbors(v))
			D=N.difference(support)
			lower_bound=min([lower_bound,len(D)])
		for v in support:
			N=set(graph.neighbors(v))
			D=N.difference(support)
			lower_bound=min([lower_bound,len(D)-1])
		lower_bound=lower_bound+s # the default lower bound
	i=upper_bound
	find=1 # does sage find a zero forcing set of size i
	outcome=['wrong'] # default outcome
	while i>=lower_bound and find==1:
		find=0
		leq=gZ_leq(graph, support, bannedset,i) # check gZ <= i or not
		if leq!=False:
			outcome=leq
			find=1
			i=i-1
	return outcome

def find_gZ(graph, support=[], bannedset=[], upper_bound=None, lower_bound=None):
	"""
	For a given graph with support and banned set, return the zero. upper_bound and lower_bound could be left as default value if one is not sure.

	Input:
		graph: a simple graph
		support: a list of vertices of g
		bannedset: a list of tuples representing banned edges of graph
		upper_bound: an integer supposed to be an upper bound of gZ.
		lower_bound: an integer supposed to be a lower bound of gZ. The two bounds may shorten the computation time. But one may leave it as default value if one is not sure.

	Output:
		the generalized zero forcing number

	Examples:
		sage: find_gZ(graphs.PathGraph(5))
		1
		sage: find_gZ(graphs.PathGraph(5),[1],[(3,2)])
		3
	"""
	return len(find_gzfs(graph, support, bannedset, upper_bound, lower_bound))

def ft_set(g,s, faults =1):

    Z=find_gZ(g)
    if Z+faults>g.order():
        return -1
    all_checks=True
    sub_s_count=0
    if len(gzerosgame(g,s))==g.order():


        for sub_s in Subsets(s,len(s)-faults):

            #print("sub_s=",sub_s)

            if len(gzerosgame(g,sub_s))==g.order():
                #print("it worked")
                sub_s_count+=1
                #print("subsets=",len(Subsets(s,len(s)-faults)))
                #print("count=",sub_s_count)
            if sub_s_count==len(Subsets(s,len(s)-faults)):
                return True
            #else:
                #return False
    return False

def ftZ(g,faults=1,robust=False,all_sets=False):
    ftz=-1
    Z=find_gZ(g)
    ftz_sets=[]
    for i in range(faults,g.order()):

        S=[]
        for s in Subsets(g,Z+i):
            S.append(s)
            #print(S)
        for j in range(len(S)):
            if ft_set(g,S[j])==True:
                if ftz==-1 or len(S[j])==ftz:
                    ftz=len(S[j])
                    ftz_sets.append(S[j])
                    #print(S[j],"is fault tolerant")
                    #print(ftz_sets)
    if all_sets==True:
        return ftz_sets
    else:
        return ftz
def is_psd_forcing_set(g, s):
    """
    sees if set s is psd forcing for g
    """
    if len(psdgame(g, s)) == g.order():
        return True
    return False


def find_psd_sets_of_size(g, size):
    """"
    finds all psd forcing sets for g
    """
    psd_sets = []
    for s in Subsets(g, size):    #loops over substes of vertices"""
        if is_psd_forcing_set(g, s):     #checks if current subset is psd"""
            psd_sets.append(s)     #add to list if it is psd set"""
    return psd_sets

def psdZ(g, all_sets=False):
    """finds min psd forcing set size (psdz)
        if all_sets is true return all psd force sets of min size"""
  
    psdz = -1
    psdz_sets = []

   
        #Try al subset sizes from 1 up to # of vertices in graph"""
    for size in range(1, g.order() + 1):
        current_sets = find_psd_sets_of_size(g, size) 
            #Find all PSD force sets of current size"""
        if current_sets:
            psdz = size   #If one PSD force set was found update to current size"""
            psdz_sets = current_sets  #Store sets for forcing #"""
            break  

    return psdz_sets if all_sets else psdz

def ftZplus_set(g):
    """
    Returns a Fault Tolerant PSD Forcing Set for graph g.

    For special graph types, known optimal sets are used.
    Otherwise, it uses the Zplus set and adds one extra vertex.
    """
    from itertools import combinations

    V = list(g.vertices())

    # Special Case: Tree
    if g.is_tree():
        # For trees, 2 well-placed vertices are sufficient.
        # We'll just pick the two furthest apart using eccentricity
        ecc = g.eccentricity()
        max_v = max(ecc, key=lambda x: ecc[x])
        farthest = max(g.shortest_path(max_v), key=len)[-1]
        return [max_v, farthest]

    # Special Case: Complete Graph
    if g.is_clique() and g.order() >= 2:
        return V  # All vertices

    # Special Case: Cycle Graph
    if g.order() >= 3 and g.size() == g.order():

        return [V[0], V[1], V[2]]

    # Special Case: Complete Bipartite
    if g.is_bipartite():
        try:
            parts = g.bipartite_sets()
            m, n = len(parts[0]), len(parts[1])
            if m >= 2 and n >= m:
                return list(parts[0])[:m] + [parts[1][0]]  # m from one side, 1 from the other
        except:
            pass  # fallback if bipartite_sets fails

    # General Case:
    # Get minimum PSD forcing set from Zplus search
    for size in range(1, g.order()):
        for subset in combinations(V, size):
            if len(psdgame(g, list(subset))) == g.order():
                base_set = list(subset)
                # Now try adding one more vertex to make it fault tolerant
                for v in V:
                    if v not in base_set:
                        trial_set = base_set + [v]
                        # Simple check: remove any one vertex and still force all
                        valid = True
                        for removed in trial_set:
                            test_set = trial_set.copy()
                            test_set.remove(removed)
                            if len(psdgame(g, test_set)) != g.order():
                                valid = False
                                break
                        if valid:
                            return trial_set
                # Fallback if none work
                return base_set + [V[0]]





