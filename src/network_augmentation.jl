function max_switch_list(G,d,fw) #computes a list of candidate switches to add the G that each covers a maximal set of edges over the set of all candidate switches of length at most d
    r = argmax(betweenness_centrality(G))
    edge_lookup = Dict{Array{Int},Int}()
    for i = 1:ne(G)
        e = collect(edges(G))[i]
        edge_lookup[[src(e),dst(e)]] = i
        edge_lookup[[dst(e),src(e)]] = i
    end
    path_to_source = []
    for i = 1:nv(G)
        p = enumerate_paths(fw,i,r)
        push!(path_to_source,[edge_lookup[[p[j],p[j+1]]] for j = 1:length(p)-1])
    end
    edge_lookup = Dict{Array{Int},Int}()
    maximal_switches = []
    for u = 1:nv(G)-1
        u_path = Set(path_to_source[u])
        u_switches = []
        potential_v = [v for v=u+1:nv(G) if slen(G,u,v) <= d]
        v_distances = [slen(G,u,v) for v in potential_v]
        potential_v = reverse(potential_v[sortperm(v_distances)])
        for v in potential_v
            v_path = Set(path_to_source[v])
            contained = 0
            proposed_cover = symdiff(u_path,v_path)
            for c in u_switches
                if issubset(proposed_cover,c)
                    contained = 1
                end
                if issubset(c,proposed_cover)
                    filter!(x->x!=c,u_switches)
                end
            end
            if contained == 0
                push!(maximal_switches,[u,v,[s for s in proposed_cover]])
                push!(u_switches,proposed_cover)
            end
        end
    end
    edge_exposure = tree_edge_weights_meta(G,r).*tree_edge_lengths(G)
    return [edge_exposure, maximal_switches]
end

function redundant_greedy_switch_addition(G,edge_exposure_input,maximal_switch_input,max_added) #Gives a subset of switches from the maximal switch list to add using the greedy algorithm. At each step, the switch covering the greatest SAIDI exposure is added, and the exposures of the edges it covers are halved.
    total_exposure = sum(edge_exposure_input)
    edge_exposure = edge_exposure_input/total_exposure
    edge_exposure_red = edge_exposure_input/total_exposure #copy of edge_exposure to be updated during the algorithm as coverage increases
    maximal_switch_list = []
    maximal_switches = deepcopy(maximal_switch_input)
    exposure_accumulated = [0.0]
    length_used = [0.0]
    for i = 1:min(max_added,length(maximal_switches))
        covered_exposure = [sum(edge_exposure_red[s[3]]) for s in maximal_switches]
        if maximum(covered_exposure) ==  0
            break
        end
        maximal_switch = deepcopy(maximal_switches[argmax(covered_exposure)])
        push!(maximal_switch_list,maximal_switch)
        push!(length_used, last(length_used) + slen(G,maximal_switch[1],maximal_switch[2]))
        push!(exposure_accumulated, last(exposure_accumulated) + sum(edge_exposure[i] for i in maximal_switch[3]))
        for i in maximal_switch[3]
            edge_exposure_red[i] /= 2
            edge_exposure[i] = 0
        end
    end
    return maximal_switch_list,length_used, exposure_accumulated,sum(edge_exposure)
end



