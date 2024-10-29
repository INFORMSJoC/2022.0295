
### Running local search 50 times on each component, for 100 time steps. 

num_local_search_steps = 100
comp=1

for comp=1:18
	for out_iter=1:50
	println("component: ",comp)
	local_search_weighted_product(contracted_graphs[comp], contracted_trees[comp,10], 0, contracted_centers[comp], num_local_search_steps,comp,out_iter)
	end
end