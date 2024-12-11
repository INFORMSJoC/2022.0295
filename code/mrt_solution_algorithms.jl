function IP_opt_reconnection(M,l) #solves the IP formulation of MRT with R-Time objective for a network with connectivity matrix M and edge probability failure vector l (lengths)
    m = size(M,1)
    n = size(M,2)
    B = JuMP.Model(GLPK.Optimizer) #IP formulation

    @variable(B, x_b[i=1:n, t=1:m], Bin)
    @variable(B, y_b[j=1:m, t=1:m], Bin)

    @constraint(B, edge_order[t=1:m], sum(y_b[j,t] for j=1:m) == 1) #connect at most one edge at time t
    @constraint(B, edge_redundancy[j=1:m], sum(y_b[j,t] for t=1:m) == 1) #for efficiency
    @constraint(B, failure_cover[i=1:n], sum(x_b[i,t] for t=1:m) == min(1,sum(M[:,i]))) #each edge that can fail is covered by at least one edge
    @constraint(B, covered_feasibility[i=1:n, t=1:m], x_b[i,t] <= sum(M[j,i] * y_b[j,t] for j=1:m)) #an edge can only be covered if a corresponding switch is turned on

    @objective(B, Min, sum(t * l[i]*x_b[i, t] for i=1:n, t=1:m)/n) #average time for network to be reconnected

    optimize!(B)

    optimal_objective = objective_value(B)

    X = zeros(m,m)
    for s = 1:m
        for t = 1:m
            X[s,t] = value.(y_b[s,t])
        end
    end
    return optimal_objective, X, [sum(t *value.(x_b[i, t]) for t=1:m) for i=1:n]
end

function IP_opt_SAIDI(M,w,l) #solves the IP formulation of MRT with SAIDI objective for a network with connectivity matrix M, edge weights w, and edge probability failure vector l (lengths)
    m = size(M,1)
    n = size(M,2)
    B = JuMP.Model(GLPK.Optimizer) #IP formulation

    @variable(B, x_b[i=1:n, t=1:m], Bin)
    @variable(B, y_b[j=1:m, t=1:m], Bin)

    @constraint(B, edge_order[t=1:m], sum(y_b[j,t] for j=1:m) == 1) #connect at most one edge at time t
    @constraint(B, edge_redundancy[j=1:m], sum(y_b[j,t] for t=1:m) == 1) #for efficiency
    @constraint(B, failure_cover[i=1:n], sum(x_b[i,t] for t=1:m) == min(1,sum(M[:,i]))) #each edge that can fail is covered by at least one edge
    @constraint(B, covered_feasibility[i=1:n, t=1:m], x_b[i,t] <= sum(M[j,i] * y_b[j,t] for j=1:m)) #an edge can only be covered if a corresponding switch is turned on

    @objective(B, Min, sum(t * l[i]*w[i]*x_b[i, t] for i=1:n, t=1:m)/n) #average time for network to be reconnected

    optimize!(B)

    optimal_objective = objective_value(B)
    X = zeros(m,m)
    for s = 1:m
        for t = 1:m
            X[s,t] = value.(y_b[s,t])
        end
    end
    return optimal_objective, X, [sum(t *value.(x_b[i, t]) for t=1:m) for i=1:n]
end

function LP_opt_reconnection(M,l) #solves the LP relaxation of MRT with R-Time objective for a network with connectivity matrix M and edge probability failure vector l (lengths)
    m = size(M,1)
    n = size(M,2)
    B = JuMP.Model(GLPK.Optimizer) #IP formulation

    @variable(B, x_b[i=1:n, t=1:m] >=0)
    @variable(B, y_b[j=1:m, t=1:m] >=0)

    @constraint(B, edge_order[t=1:m], sum(y_b[j,t] for j=1:m) == 1) #connect at most one edge at time t
    @constraint(B, edge_redundancy[j=1:m], sum(y_b[j,t] for t=1:m) == 1) #for efficiency
    @constraint(B, failure_cover[i=1:n], sum(x_b[i,t] for t=1:m) == min(1,sum(M[:,i]))) #each edge that can fail is covered by at least one edge
    @constraint(B, covered_feasibility[i=1:n, t=1:m], x_b[i,t] <= sum(M[j,i] * y_b[j,t] for j=1:m)) #an edge can only be covered if a corresponding switch is turned on

    @objective(B, Min, sum(t * l[i]*x_b[i, t] for i=1:n, t=1:m)/n) #average time for network to be reconnected

    optimize!(B)

    optimal_objective = objective_value(B)

    X = zeros(m,m)
    for s = 1:m
        for t = 1:m
            X[s,t] = value.(y_b[s,t])
        end
    end
    return optimal_objective, X,[sum(t *value.(x_b[i, t]) for t=1:m) for i=1:n]
end

function LP_opt_SAIDI(M,w,l) #solves the LP relaxation of MRT with SAIDI objective for a network with connectivity matrix M, edge weights w, and edge probability failure vector l (lengths)
    m = size(M,1)
    n = size(M,2)
    B = JuMP.Model(GLPK.Optimizer) #IP formulation

    @variable(B, x_b[i=1:n, t=1:m] >=0)
    @variable(B, y_b[j=1:m, t=1:m] >=0)

    @constraint(B, edge_order[t=1:m], sum(y_b[j,t] for j=1:m) == 1) #connect at most one edge at time t
    @constraint(B, edge_redundancy[j=1:m], sum(y_b[j,t] for t=1:m) == 1) #for efficiency
    @constraint(B, failure_cover[i=1:n], sum(x_b[i,t] for t=1:m) == min(1,sum(M[:,i]))) #each edge that can fail is covered by at least one edge
    @constraint(B, covered_feasibility[i=1:n, t=1:m], x_b[i,t] <= sum(M[j,i] * y_b[j,t] for j=1:m)) #an edge can only be covered if a corresponding switch is turned on

    @objective(B, Min, sum(t * l[i]*w[i]*x_b[i, t] for i=1:n, t=1:m)/n) #average time for network to be reconnected

    optimize!(B)

    optimal_objective = objective_value(B)
    X = zeros(m,m)
    for s = 1:m
        for t = 1:m
            X[s,t] = value.(y_b[s,t])
        end
    end
    return optimal_objective, X,[sum(t *value.(x_b[i, t]) for t=1:m) for i=1:n]
end

function IP_solution_to_order(M)
    [round(argmax(M[:,t])) for t = 1:size(M)[2]]
end

function greedy_reconnection(M,l) #applies the greedy algorithm for MRT with R-Time objective to a network with connectivity matrix M and edge probability failure vector l (lengths)
    local C
    C = 1*M
    m = size(C, 1)
    n = size(C, 2)
    for i=1:m
        for j=1:n
            C[i,j]=C[i,j]*l[j]
        end
    end
    counter = 0
    cost = 0
    order = []
    while sum(sum(C, dims=1), dims=2)[1,1] > 0
        counter += 1
        D = sum(C, dims=2)
        i = argmax(D, dims=1)[1][1]
        edge_cost = D[i,1]
        for j = 1:n
            if C[i,j] > 0
                for k = 1:m
                    C[k,j] = 0
                end
            end
        end
        cost += counter * edge_cost
        push!(order,i)
    end
    return([cost/n,order])
end

function greedy_SAIDI(M,w,l) #applies the greedy algorithm for MRT with SAIDI objective to a network with connectivity matrix M, edge weight vector w, and edge probability failure vector l (lengths)
    local C
    C = 1*M
    m = size(C, 1)
    n = size(C, 2)
    for i=1:m
        for j=1:n
            C[i,j]=C[i,j]*w[j]*l[j]
        end
    end
    counter = 0
    cost = 0
    order = []
    while sum(sum(C, dims=1), dims=2)[1,1] > 0
        counter += 1
        D = sum(C, dims=2)
        i = argmax(D, dims=1)[1][1]
        edge_cost = D[i,1]
        cost += counter * edge_cost
        for j = 1:n
            if C[i,j] > 0
                for k = 1:m
                    C[k,j] = 0
                end
            end
        end
        append!(order,i)
    end
    return([cost/n,order])
end

function reconnection_cost(M,l, order) #given a reconnection order, return the R-Time metric value for a network with connectivity matrix M and failure probabilities l (proportional to edge length)
    local C
    C = 1*M
    m = size(C, 1)
    n = size(C, 2)
    for i=1:m
        for j=1:n
            C[i,j]=C[i,j]*l[j]
        end
    end
    cost = 0
    for counter = 1:length(order)
        D = sum(C, dims=2)
        i = order[counter]
        edge_cost = D[i,1]
        cost += counter * edge_cost
        for j = 1:n
            if C[i,j] > 0
                for k = 1:m
                    C[k,j] = 0
                end
            end
        end
    end
    return(cost/n)
end

function SAIDI_cost(M,w,l,order) #given a reconnection order, return the SAIDI metric value for a network with connectivity matrix M, weights w, and lengths l
    local C
    C = 1*M
    m = size(C, 1)
    n = size(C, 2)
    for i=1:m
        for j=1:n
            C[i,j]=C[i,j]*w[j]*l[j]
        end
    end
    cost = 0
    for counter = 1:length(order)
        D = sum(C, dims=2)
        i = order[counter]
        edge_cost = D[i,1]
        cost += counter * edge_cost
        for j = 1:n
            if C[i,j] > 0
                for k = 1:m
                    C[k,j] = 0
                end
            end
        end
    end
    return(cost/n)
end

function greedy_product(M,w,l) #modified greedy algorithm for MRT which at each step selects the switch with the greatest marginal increase in the product of the weights used for SAIDI and {\sc R-Time}, respectively.
    m = size(M, 1)
    n = size(M, 2)
    C = zeros(m,n)
    S = zeros(m,n)
    for i=1:m
        for j=1:n
            S[i,j]=M[i,j]*w[j]*l[j]
            C[i,j]=M[i,j]*l[j]
        end
    end
    uncovered_edges = [i for i=1:n if sum(M[:,i]) == 0]
    recon_cost = sum(l[i] for i in uncovered_edges)
    saidi_cost = sum(w[i]*l[i] for i in uncovered_edges)
    order = []
    for counter = 1:m
        D = sum(C, dims=2)
        DS = sum(S, dims=2)
        i = argmax([(recon_cost + counter * D[s,1]) * (saidi_cost + counter * DS[s,1]) for s = 1:m])
        recon_cost += counter * D[i,1]
        saidi_cost += counter * DS[i,1]
        for j = 1:n
            if C[i,j] > 0
                for k = 1:m
                    C[k,j] = 0
                    S[k,j] = 0
                end
            end
        end
        push!(order,i)
    end
    return([recon_cost*saidi_cost/n^2,recon_cost/n,saidi_cost/n,order])
end

function msalphac(M,x,w,c,iter) #given a connectivity matrix M, LP relaxation solution x, weights w, and covering parameter c, run alpha-point rounding for the specified number of iterations
    m = size(M,1)
    n = size(M,2)
    p = 2/(c-1)
    K = zeros(m,2m)
    for t = 1:2m
        Kt = c*p/(t*(t+1)^p)
        for tp = 1:m
            K[tp,t]=Kt * tp^p
        end
    end
    z = zeros(m,2m)
    for s = 1:m
        alpha = rand()
        t = 1
        maxed = 0
        while maxed == 0 && t<=2m
            z[s,t] = sum(K[tp,t]*x[s,tp] for tp =1:min(t,m))
            if z[s,t] >=1
                maxed = 1
            end
            t+=1
        end
    end
    sumz = zeros(m,2m)
    for s = 1:m
        for t = 1:2m
            sumz[s,t] = sum(z[s,tp] for tp=1:t)
        end
    end
    results = []
    for sample = 1:iter
        threshold = [rand() for i =1:m]
        tau = [minimum([2m+(t-2m)*(sumz[s,t]>=threshold[s]) for t=1:2m]) for s=1:m]
        sigma = sortperm(tau)
        objective = 0
        covered = deepcopy(M)
        for t = 1:m
            for i = 1:n
                if covered[sigma[t],i] == 1
                    objective += t*w[i]
                    for j = 1:m
                        covered[j,i] = 0
                    end
                end
            end
        end
        push!(results,deepcopy(objective)/n)
    end
    return results
end

function vertex_coverage_time(G,T,source,lengths,M,order,row_index,col_index) #returns a vector of expected reconnection times after each vertex experiences an outage, given a specified covering order
    vertex_times = [[] for i = 1:nv(G)]
    n = nv(T)
    E = collect(edges(G))
    w = []
    for i = 1:ne(G)
        if i in col_index
            e = E[i]
            e_weight = [slen(G,src(e),dst(e))]
            r_time = 0
            for j = 1:length(order)
                if M[row_index[order[j]],i] == 1
                    r_time = j
                    break
                end
            end
            if r_time > 0
                rem_edge!(T,e)
                D = gdistances(T, source; sort_alg=QuickSort)
                add_edge!(T,e)
                v_below = []
                for v in vertices(T)
                    if D[v] > 100000
                        push!(vertex_times[v],[r_time,e_weight])
                    end
                end 
            end 
        end 
    end
    waiting_times = []
    for v = 1:n 
        if length(vertex_times[v]) == 0
            push!(waiting_times,0)
        else
            push!(waiting_times,sum(vertex_times[v][i][1]*vertex_times[v][i][2] for i = 1:length(vertex_times[v]))/sum(vertex_times[v][i][2] for i = 1:length(vertex_times[v])))
        end 
    end 
    return waiting_times
end