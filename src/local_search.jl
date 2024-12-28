function update_connectivity(old_M,G,T,e,s_e) #updates the connectivity matrix old_M of G and T after removing edge e and adding switch s_e from and to T, respectively
    M = deepcopy(old_M)
    E = collect(edges(G))
    n = nv(T)
    Tnew = copy(T)
    rem_edge!(Tnew,e)
    add_edge!(Tnew,s_e)
    C_e = [] #edges covered by e in T
    Cnew_e = [] #edges covered by e in Tnew
    for f in edges(Tnew)
        if exchange_connected(T,f,e)
            push!(C_e,f)
        end
        if exchange_connected(Tnew,f,e)
            push!(Cnew_e,f)
        end
    end
    C_e_unique = C_e
    append!(C_e_unique, Cnew_e)
    C_e_unique = unique(C_e_unique)
    C_s = [] #switches covering e in T
    Cnew_s = [] #switches covering s_e in Tnew
    for s in edges(G)
        if !(s in edges(T))
            if exchange_connected(T,e,s)
                push!(C_s,s)
            end
        end
        if !(s in edges(Tnew))
            if exchange_connected(Tnew,s_e,s)
                push!(Cnew_s,s)
            end
        end
    end
    C_s_unique = C_s
    append!(C_s_unique, Cnew_s)
    C_s_unique = unique(C_s_unique)
    eindex = findlast(j->j == e, E)
    s_eindex = findlast(j->j == s_e, E)
    for f in C_e_unique
        i = findlast(j->j == f, E)
        for s in C_s_unique
            if (f in C_e && s in C_s) || (f in Cnew_e && s in Cnew_s)
                M[findlast(j->j == s, E), i] = exchange_connected(Tnew, f, s)
                #counter += 1
            end
        end
    end
    for i = 1:ne(G)
        M[eindex, i] = exchange_connected(Tnew, E[i], e)
        M[i, s_eindex] = exchange_connected(Tnew, s_e, E[i])
        M[s_eindex, i] = 0
        M[i, eindex] = 0
    end
    M
end

function local_search_weighted_product(G,start_T,start_M,source,steps,out_comp,out_iter) #Given graph G with given starting spanning tree start_T, optional connectivity matrix start_M, and source, take up to the specified number of branch exchange steps that improve the product of the three objectives. The search stops if none of the stopping_parameter many exchange attempts succeed.
    results = Array{Float64,2}(undef,steps+1,4)
    trees = Vector{Any}(undef,steps)
    T = deepcopy(start_T)
    E = deepcopy(collect(edges(G)))
    m = ne(G)
    graph_results = Array{Float64,2}(undef,steps+1,4)
    w =  tree_edge_weights_meta(T, source)
    l = tree_edge_lengths(T)
    if start_M == 0
        M = connectivity(G,T)
    else
        M = start_M
    end
    counter = 0
    local_min = false
    finalM = M
    for k = 1:steps #search iterations
        improved = false
        counter = 0 #number of exchanges attempted in this step
        row_indices = filter(j -> !(E[j] in collect(edges(T))),[i for i = 1:m])
        column_indices = filter(j -> E[j] in collect(edges(T)),[i for i = 1:m])
        trees[k] = column_indices
        product_results = greedy_product(M[row_indices, column_indices], w,l)
        incumbent_reconnect = product_results[2]
        incumbent_SAIDI = product_results[3]
        incumbent_energy = tree_energy_cost(w,l)
        incumbent_product = incumbent_reconnect*incumbent_SAIDI*incumbent_energy
        graph_results[k,1] = incumbent_reconnect
        graph_results[k,2] = incumbent_SAIDI
        graph_results[k,3] = incumbent_energy
        graph_results[k,4] = incumbent_product
        possible_exchanges = shuffle(filter(exchange -> M[exchange[1],exchange[2]] == 1, [[s,t] for s in row_indices, t in column_indices]))
        if !local_min
            while !improved && counter < length(possible_exchanges)
                counter += 1
                exchange_attempt = possible_exchanges[counter]
                edge_out = E[exchange_attempt[2]]
                edge_in = E[exchange_attempt[1]]
                newM = update_connectivity(M,G,T,edge_out,edge_in)
                add_edge!(T,edge_in)
                rem_edge!(T,edge_out)
                w =  tree_edge_weights_meta(T, source)
                l = tree_edge_lengths(T)
                E_T = collect(edges(T))
                row_indices = filter(j -> !(E[j] in E_T),[i for i = 1:m])
                column_indices = filter(j -> E[j] in E_T,[i for i = 1:m])
                newM_T = newM[row_indices, column_indices]
                new_product_results = greedy_product(newM[row_indices, column_indices],w,l)
                new_reconnect = new_product_results[2]
                new_SAIDI = new_product_results[3]
                new_energy = tree_energy_cost(w,l)
                new_product = new_reconnect*new_SAIDI*new_energy
                if new_product < incumbent_product
                    improved = true
                end
                if !improved
                    rem_edge!(T,edge_in)
                    add_edge!(T,edge_out)
                end
                if improved
                    M = newM
                    w =  tree_edge_weights_meta(T, source)
                    l = tree_edge_lengths(T)
                    println("comp $out_comp iter $out_iter step $k after $counter attempts")
                    if k == steps
                        graph_results[k+1,1] = new_reconnect
                        graph_results[k+1,2] = new_SAIDI
                        graph_results[k+1,3] = new_energy
                        graph_results[k+1,4] = new_product
                    end
                end
            end
        else
            graph_results[k+1,:] = graph_results[k,:]
        end
        if !improved
            local_min=true
        end
        if k==steps
            finalM = M
        end
    end
    best_reconnect = minimum(graph_results[:,1])
    best_SAIDI = minimum(graph_results[:,2])
    best_energy = minimum(graph_results[:,3])
    best_product = minimum(graph_results[:,4])
    results[:,1] = graph_results[:,1]/best_reconnect
    results[:,2] = graph_results[:,2]/best_SAIDI
    results[:,3] = graph_results[:,3]/best_energy
    results[:,4] = graph_results[:,4]/best_product
    gr_df = DataFrame(out_comp = [comp for i = 1:size(results)[1]], iter = [out_iter for i = 1:size(results)[1]], RTime =graph_results[:,1],SAIDI =graph_results[:,2],Energy =graph_results[:,3],Product =graph_results[:,4])
    outpath = string("output/","$out_comp","_iter_$out_iter",".csv")
    CSV.write(outpath, gr_df,writeheader=true)
    return [graph_results, results, T,finalM,trees]
end
