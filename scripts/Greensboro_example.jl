using CSV,DataFrames,Graphs,MetaGraphs,SimpleWeightedGraphs,Base.Threads,StatsBase,Random,JuMP,GLPK,Distances,Pkg,Tables,Statistics

include("../src/local_search.jl")
include("../src/mrt_solution_algorithms.jl") 
include("../src/network_augmentation.jl")
include("../src/network_statistics_and_modification.jl")

rng = Random.MersenneTwister(182) #reset random seed

Bus_location_list = CSV.read("../data/Greensboro_US_buses.csv",DataFrame)
Line_list = CSV.read("../data/Greensboro_US_lines.csv",DataFrame)

vertex_count = size(Bus_location_list)[1]
edge_count = size(Line_list)[1]

G_ns = MetaGraph(vertex_count) #create the base network without auxiliary switches
for i = 1:vertex_count
    set_prop!(G_ns, i, :kW, Bus_location_list[i,5])
    set_prop!(G_ns, i, :kVAr, Bus_location_list[i,6])
    set_prop!(G_ns, i, :weight, get_prop(G_ns, i, :kW)) #set the default weight used in load calculations
    set_prop!(G_ns, i, :apparent_power, Bus_location_list[i,7])
    set_prop!(G_ns, i, :name, Bus_location_list[i,1])
    set_prop!(G_ns, i, :long, Bus_location_list[i,3])
    set_prop!(G_ns, i, :lat, Bus_location_list[i,4])
end
for i = 1:edge_count
    if Line_list[i,4] == "n"
        u = Line_list[i,1]
        v = Line_list[i,2]
        add_edge!(G_ns,u,v)
        set_prop!(G_ns,u,v,:resistance,Line_list[i,5]) #resistance proportional to line length
        set_prop!(G_ns,u,v,:type,Line_list[i,3])
    end
end

C_ns = connected_components(G_ns) #Break the network into its constituent connected components. Component #19 corresponds to a portion of the transmission network, so we set it aside.
no_switch_component_index = Dict{Int,Tuple{Int64,Int64}}()
for i = 1:vertex_count
    for comp = 1:length(C_ns)
        if count(k->(k == i), C_ns[comp]) > 0
            no_switch_component_index[i] = (comp, count(k->(k<=i), C_ns[comp]))
            set_prop!(G_ns, i, :component, convert(Int,comp))
            set_prop!(G_ns, i, :originalindex, convert(Int,i))
        end
    end
end

contraction_threshold = 10
length_threshold = 1000 #maximum switch length of 1000m
contracted_graphs = Dict{Int,MetaGraph}()
contracted_trees = Dict{Tuple,MetaGraph}()
contracted_fw = Dict{Tuple,Graphs.FloydWarshallState{Float64,Int64}}()
contracted_connectivity_matrices = Dict{Int,Array}()
chosen_switches = Dict{Int,Vector}()

#indices of switch and non-switch edges in the contracted graphs
contracted_rows = Dict{Int,Vector}()
contracted_cols = Dict{Int,Vector}()
contracted_coverage = Dict{Int,Vector}()
contracted_centers = Dict{Int,Int}()
weights = Dict{Int,Vector}()
lengths = Dict{Int,Vector}()
reduced_connectivity_matrices = Dict{Int,Array}()

#solutions to MRT with both objectives
IP_SAIDI_results = Dict{Int,Tuple{Float64,Matrix{Float64},Vector{Float64}}}()
IP_recon_results = Dict{Int,Tuple{Float64,Matrix{Float64},Vector{Float64}}}()

#objective values for both SAIDI and R-Time, depending on which was selected as the MRT objective
IP_SAIDI_S_results = Dict{Int,Float64}()
IP_recon_R_results = Dict{Int,Float64}()
IP_SAIDI_R_results = Dict{Int,Float64}()
IP_recon_S_results = Dict{Int,Float64}()

max_switches = Dict{Tuple{Int, Int, Int}, Vector}() 
red_greedy_switches = Dict{Tuple{Int,Int,Int}, Tuple{Vector{Any}, Vector{Float64}, Vector{Float64}, Float64}}() 

#expected outage times for each vertex, depending on the MRT objective
SAIDI_coverage_times = Dict{Int,Vector}()
RT_coverage_times = Dict{Int,Vector}()


## There are 18 components. I am currently only running for a single component. 
for comp = 1:18 #can be multithreaded using @threads functionality as desired
    println("running for component", comp) 
    contracted_trees[(comp,contraction_threshold)] = contract_metatree(induced_subgraph(G_ns,C_ns[comp])[1],contraction_threshold)
    contracted_fw[(comp,contraction_threshold)] = floyd_warshall_shortest_paths(contracted_trees[(comp,contraction_threshold)])
    max_switches[(comp,contraction_threshold,length_threshold)] = max_switch_list(contracted_trees[(comp,contraction_threshold)],length_threshold,floyd_warshall_shortest_paths(contracted_trees[(comp,contraction_threshold)]))
    red_greedy_switches[(comp,contraction_threshold,length_threshold)] =
    redundant_greedy_switch_addition(contracted_trees[(comp,contraction_threshold)],max_switches[(comp,contraction_threshold,length_threshold)][1],max_switches[(comp,contraction_threshold,length_threshold)][2],200)
    contracted_graphs[comp] = deepcopy(contracted_trees[(comp,contraction_threshold)])
    S = red_greedy_switches[(comp,contraction_threshold,length_threshold)][1]
    for s in S[1:min(50,length(S))]
        add_edge!(contracted_graphs[comp],s[1],s[2])
    end
    for i = 1:edge_count
        if Line_list[i,4] == "y" && no_switch_component_index[Line_list[i,1]][1] == comp
            u = no_switch_component_index[Line_list[i,1]][2]
            v = no_switch_component_index[Line_list[i,2]][2]
            add_edge!(contracted_graphs[comp],u,v)
            set_prop!(contracted_graphs[comp],u,v,:resistance,Line_list[i,5]) #resistance proportional to line length
            set_prop!(contracted_graphs[comp],u,v,:type,Line_list[i,3])
        end
    end
    contracted_centers[comp] = argmax(betweenness_centrality(contracted_graphs[comp]))
    contracted_coverage = last(red_greedy_switches[(comp,contraction_threshold,length_threshold)][3])
    contracted_connectivity_matrices[comp]=connectivity(contracted_graphs[comp],contracted_trees[(comp,contraction_threshold)])
    contracted_cols[comp] = []
    contracted_rows[comp] = []
    for i = 1:ne(contracted_graphs[comp])
        e = collect(edges(contracted_graphs[comp]))[i]
        if has_edge(contracted_trees[comp,contraction_threshold],e)
            push!(contracted_cols[comp],i)
        else
            push!(contracted_rows[comp],i)
        end
    end
    reduced_connectivity_matrices[comp] = contracted_connectivity_matrices[comp][contracted_rows[comp],contracted_cols[comp]]
    weights[comp] = tree_edge_weights_meta(contracted_trees[(comp,contraction_threshold)],contracted_centers[comp])
    lengths[comp] = tree_edge_lengths(contracted_trees[(comp,contraction_threshold)])
    IP_SAIDI_results[comp] = IP_opt_SAIDI(reduced_connectivity_matrices[comp], weights[comp], lengths[comp])
    IP_recon_results[comp] = IP_opt_reconnection(reduced_connectivity_matrices[comp], lengths[comp])
    IP_SAIDI_S_results[comp] = IP_SAIDI_results[comp][1][1]
    IP_recon_R_results[comp] = IP_recon_results[comp][1][1]
    IP_SAIDI_R_results[comp] = reconnection_cost(reduced_connectivity_matrices[comp],  lengths[comp],IP_solution_to_order(IP_SAIDI_results[comp][2]))
    IP_recon_S_results[comp] = SAIDI_cost(reduced_connectivity_matrices[comp], weights[comp], lengths[comp],IP_solution_to_order(IP_recon_results[comp][2]))
    SAIDI_coverage_times[comp] = vertex_coverage_time(contracted_graphs[comp],contracted_trees[(comp,contraction_threshold)],contracted_centers[comp],lengths[comp],contracted_connectivity_matrices[comp],IP_solution_to_order(IP_SAIDI_results[comp][2]),contracted_rows[comp],contracted_cols[comp])
    RT_coverage_times[comp] = vertex_coverage_time(contracted_graphs[comp],contracted_trees[(comp,contraction_threshold)],contracted_centers[comp],lengths[comp],contracted_connectivity_matrices[comp],IP_solution_to_order(IP_recon_results[comp][2]),contracted_rows[comp],contracted_cols[comp]) 

end

# Convert and save IP_SAIDI_results
IP_SAIDI_results_df = DataFrame([(k, v[1], v[2], v[3]) for (k, v) in IP_SAIDI_results], [:ID, :ObjectiveValue, :Matrix, :Vector])
CSV.write("../output/IP_SAIDI_results.csv", IP_SAIDI_results_df)

# Convert and save IP_recon_results
IP_recon_results_df = DataFrame([(k, v[1], v[2], v[3]) for (k, v) in IP_recon_results], [:ID, :ObjectiveValue, :Matrix, :Vector])
CSV.write("../output/IP_recon_results.csv", IP_recon_results_df)

# Save the objective values
IP_SAIDI_S_results_df = DataFrame([(k, v) for (k, v) in IP_SAIDI_S_results], [:ID, :SAIDIObjective])
CSV.write("../output/IP_SAIDI_S_results.csv", IP_SAIDI_S_results_df)

IP_recon_R_results_df = DataFrame([(k, v) for (k, v) in IP_recon_R_results], [:ID, :RTimeObjective])
CSV.write("../output/IP_recon_R_results.csv", IP_recon_R_results_df)

IP_SAIDI_R_results_df = DataFrame([(k, v) for (k, v) in IP_SAIDI_R_results], [:ID, :SAIDIObjectiveWithRTime])
CSV.write("../output/IP_SAIDI_R_results.csv", IP_SAIDI_R_results_df)

IP_recon_S_results_df = DataFrame([(k, v) for (k, v) in IP_recon_S_results], [:ID, :RTimeObjectiveWithSAIDI])
CSV.write("../output/IP_recon_S_results.csv", IP_recon_S_results_df)

# Save expected outage times
SAIDI_coverage_times_df = DataFrame([(k, v) for (k, v) in SAIDI_coverage_times], [:ID, :SAIDITimes])
CSV.write("../output/SAIDI_coverage_times.csv", SAIDI_coverage_times_df)

RT_coverage_times_df = DataFrame([(k, v) for (k, v) in RT_coverage_times], [:ID, :RTimes])
CSV.write("../output/RT_coverage_times.csv", RT_coverage_times_df)


