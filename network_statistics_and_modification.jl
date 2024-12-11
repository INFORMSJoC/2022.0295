function tree_edge_weights(T, source) #computes number of vertices without power for each edge failure
    n = nv(T)
    E = edges(T)
    w = []
    W = sum(1 for v in vertices(T))
    for e in E
        rem_edge!(T,e)
        D = gdistances(T, source; sort_alg=QuickSort)
        new_weight = W
        for v in vertices(T)
            if D[v] < 100000
                new_weight = new_weight- 1
            end
        end
        append!(w, new_weight)
        add_edge!(T,e)
    end
    return w 
end

function tree_edge_weights_meta(T, source) #computes total weight of vertices without power for each edge failure
    n = nv(T)
    E = edges(T)
    w = []
    if typeof(T) == SimpleGraph{Int64}
        W = nv(T)
    else
        W = sum(get_prop(T,v,:weight) for v in vertices(T))
    end
    for e in E
        U = rem_edge!(T,e)
        D = gdistances(T, source; sort_alg=QuickSort)
        new_weight = W
        for v in vertices(T)
            if D[v] < 100000
                if typeof(T) == SimpleGraph{Int64}
                    new_weight -= 1
                else
                    new_weight = new_weight- get_prop(T,v,:weight)
                end
            end
        end
        append!(w, new_weight)
        add_edge!(T,e)
    end
    return w
end

function tree_edge_lengths(G) #computes vector of tree edge lengths
    [slen(G,src(e),dst(e)) for e in edges(G)]
end

function delete_bridge_vertex(G,v) #deletes bridge vertex v from G and splits its weight between its neighbors
    u = neighbors(G,v)[1]
    w = neighbors(G,v)[2]
    add_edge!(G,u,w)
    set_prop!(G, u, :weight, get_prop(G,u,:weight) + get_prop(G,v,:weight)/2)
    set_prop!(G, w, :weight, get_prop(G,w,:weight) + get_prop(G,v,:weight)/2)
    rem_vertex!(G,v)
    G
end

function tree_energy_cost(w,lengths) #computes the total energy loss of a tree with matching edgewise vectors of flow w and lengths
    sum(w[i]^2*lengths[i] for i=1:length(w))
end

function contract_leaf(G,e) #contracts the specified edge and assiG_ns the sum of the weights of its incident vertices to the new contracted vertex
    if typeof(e) == Graphs.SimpleGraphs.SimpleEdge{Int64}
        if has_edge(G,e)
            if degree(G,src(e))==1
                u=src(e)
                v=dst(e)
            else
                u=dst(e)
                v=src(e)
            end
            set_prop!(G, v, :long, (get_prop(G,u,:weight)*get_prop(G,u,:long) + get_prop(G,v,:weight)*get_prop(G,v,:long))/(get_prop(G,u,:weight)+get_prop(G,v,:weight)))
            set_prop!(G, v, :lat, (get_prop(G,u,:weight)*get_prop(G,u,:lat) + get_prop(G,v,:weight)*get_prop(G,v,:lat))/(get_prop(G,u,:weight)+get_prop(G,v,:weight)))
            rem_vertex!(G,u)
        end
    else
        u = e[1]
        v = e[2]
        if get_prop(G,u,:weight)+get_prop(G,v,:weight) == 0
           set_prop!(G, v, :long, (get_prop(G,u,:long) + get_prop(G,v,:long))/2)
           set_prop!(G, v, :lat, (get_prop(G,u,:lat) + get_prop(G,v,:lat))/2)
        else
            set_prop!(G, v, :long, (get_prop(G,u,:weight)*get_prop(G,u,:long) + get_prop(G,v,:weight)*get_prop(G,v,:long))/(get_prop(G,u,:weight)+get_prop(G,v,:weight)))
            set_prop!(G, v, :lat, (get_prop(G,u,:weight)*get_prop(G,u,:lat) + get_prop(G,v,:weight)*get_prop(G,v,:lat))/(get_prop(G,u,:weight)+get_prop(G,v,:weight)))
            set_prop!(G, v, :weight, get_prop(G,u,:weight) + get_prop(G,v,:weight))
        end
        rem_vertex!(G,u)
    end
    G
end

function exchange_connected(T,e,s) #perform the branch exchange step of removing edge e and adding switch s
    if has_edge(T,s) || !has_edge(T,e)
        return false
    else
        rem_edge!(T,e)
        add_edge!(T,s)
        r = is_connected(T)
        rem_edge!(T,s)
        add_edge!(T,e)
        return r
    end
end

function connectivity(G,T) #computes the connectivity binary matrix of network G with spanning tree T: each row and each column of the matrix is indexed by the edges of G, and all entries are 0 except those where the row corresponds to a non-tree edge s, the column corresponds to a tree edge f, and T-f+s remains a tree
    U = deepcopy(T)
    E = collect(edges(G))
    n = nv(U)
    M = zeros(ne(G),ne(G))
    for i = 1:length(E)
        s = E[i]
        for j = 1:length(E)
            f = E[j]
            if has_edge(U,s) || !has_edge(U,f)
                M[i,j] = 0
            else
                M[i,j] = exchange_connected(U,f,s)
            end
        end
    end
    M
end

function contract_metatree(G, weight_threshold) #contracts all bridges of network G and all leaves that are under the specified weight threshold
    H = deepcopy(G)
    n = nv(H)
    update_made = 1
    steps = 0
    while update_made == 1
        update_made = 0
        steps += 1
        for v in vertices(H)
            if has_vertex(H,v)
                if get_prop(H,v,:weight) < weight_threshold && Graphs.degree(H,v) == 1
                    H = contract_leaf(H, [v, neighbors(H,v)[1]])
                    update_made = 1
                end
            end
        end
        for v in vertices(H)
            if has_vertex(H,v)
                if Graphs.degree(H,v) == 2
                    H = delete_bridge_vertex(H,v)
                    update_made = 1
                end
            end
        end
    end
    H
end

function slen(tree,u,v) #computes the length of the edge between u and v using the geodesic
    round(haversine((get_prop(tree,u,:long),get_prop(tree,u,:lat)),(get_prop(tree,v,:long),get_prop(tree,v,:lat))),digits=3)
end


