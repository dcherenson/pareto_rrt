module ParetoRRTStar
using StaticArrays
abstract type AbstractProblem{T} end

const ParetoPath{N,F} = NamedTuple{(:path, :cost), Tuple{Vector{Int64}, SVector{N,F}}}

struct Node{T, N, F}
    state::T
    pareto_paths::Vector{ParetoPath{N,F}}
end   

# recursive definition of cost to get to a node
# function pareto_front_costs(i::Int64, nodes::Vector{Node{T,N,F}}) where {T,N,F}
#     if isempty(nodes[i].pareto_paths)
#         return [zeros(SVector{N})]
#     end

#     pareto_front = SVector{N,F}[]
#     for parent in nodes[i].pareto_paths
#         append!(pareto_front, pareto_front_costs(parent.index, nodes) .+ Ref(parent.incremental_cost))
#     end
#     return pareto_front
    
# end

function dominates(cost1, cost2, ϵ)
    return all(cost1 .<= cost2 .+ ϵ) && any(cost1 .< cost2 .+ ϵ)
    # return all((1 + ϵ[1]) .* cost1 .<= cost2)
end


function compute_pareto_front(cand_paths::Vector{ParetoPath{N,F}}, ϵ::SVector{N,F}) where {N,F}
    # get the maximum point set of a pareto front
    # this is a set of points that are not dominated by any other point in the set
    # this is a brute force method, but it works for small sets
    paths = ParetoPath{N,F}[]
    skip_list = falses(length(cand_paths))
    for i in 1:length(cand_paths)
        if skip_list[i]
            continue
        end
        # check if this point is dominated by any other point in the set
        dominated = false
        for j in 1:length(cand_paths)
            if i != j && dominates(cand_paths[j].cost, cand_paths[i].cost, ϵ)
                if !dominates(cand_paths[i].cost, cand_paths[j].cost, ϵ)
                    dominated = true
                    break
                else # mutually dominated
                    # if mutally dominated, then we can skip j
                    skip_list[j] = true
                end
            end
        end
        if !dominated
            push!(paths, cand_paths[i])
        end
        # if length(parents) >= 10
        #     break
        # end
    end
    @assert !isempty(paths)
    return paths
end

# function compute_pareto_front(cand_parents::Vector{Int64}, costs::Vector{SVector{2,F}}, ϵ::SVector{2,F}) where {F}
#     # https://en.wikipedia.org/wiki/Maxima_of_a_point_set 
#     #sort by first element
#     sorted_idxs = sortperm(costs, by = x -> x[1])

#     front = SVector{2,F}[]
#     parents = Int64[]
#     min_y = Inf
#     min_x = -Inf
    
#     for i in sorted_idxs
#         if costs[i][2] < min_y + ϵ[2] && costs[i][1] >= min_x + ϵ[1]
#             push!(front, costs[i])
#             push!(parents, cand_parents[i])
#             min_y = costs[i][2]
#             min_x = costs[i][1]
#         end
#     end
#     return parents, front
# end

# a method to keep growing a tree
function rrt_star!(problem, nodes, max_iters; do_rewire=true)

    # i -> index in nodes vector
    # x -> actual state
    # n -> full node

    sizehint!(nodes, length(nodes)+max_iters)

    for iter = 1:max_iters

        # grab a random state
        x_rand = sample_free(problem)

        # find the nearest node
        i_nearest = nearest(problem, nodes, x_rand)
        x_nearest = nodes[i_nearest].state

        # get the new point
        x_new = steer(problem, x_nearest, x_rand)

        # check that it is obstacle free
        collfree =  collision_free(problem, x_nearest, x_new)

        if collfree 
            # get the set of nearby nodes (this should return an index set)
            I_near = near(problem, nodes, x_new)

            # determine the best one to connect to
            inc_cost = path_cost(problem, x_nearest, x_new)
            if length(nodes[i_nearest].pareto_paths) == 0
                # if the nearest node has no pareto paths, then we need to add it
                p_cand = [(path=[i_nearest], cost=inc_cost)]
            else
                p_cand = [(path=[p.path; i_nearest], cost=p.cost .+ inc_cost) for p in nodes[i_nearest].pareto_paths]
            end

            for i_near in I_near
                if i_near == i_nearest
                    continue
                end
                x_near = nodes[i_near].state
                inc_cost = path_cost(problem, x_near, x_new)
                if length(nodes[i_near].pareto_paths) == 0
                    # if the nearest node has no pareto paths, then we need to add it
                    p_near = [(path=[i_near], cost=inc_cost)]
                else
                    p_near = [(path=[p.path; i_near], cost=p.cost .+ inc_cost) for p in nodes[i_near].pareto_paths]
                end
                
                if collision_free(problem, x_near, x_new)
                    p_cand = compute_pareto_front([p_cand; p_near], problem.ϵ)
                end
            end
            # if !isfinite(c_min)
            #     continue
            # end
            
            # add in the new edge
            n_new = Node(x_new, p_cand)
            push!(nodes, n_new)

            iter % 10 == 0 && println("iter: $iter, new node: $(length(nodes)), pareto parents: $(length(p_cand))")

            # rewire
            # if do_rewire
            #     i_new = length(nodes)
            #     for i_near in I_near
            #         if nodes[i_near].parent_index != 0
            #             x_near = nodes[i_near].state
            #             pc = path_cost(problem, x_new, x_near) 
            #             if (cost(i_new, nodes) + pc < cost(i_near, nodes) ) && collision_free(problem, x_new, x_near)
            #                 # # change the parent of i_near to i_new
            #                 # nodes[i_near].parent_index = i_new
            #                 # nodes[i_near].incremental_cost = pc
            #                 new_node = Node(nodes[i_near].state, i_new, pc)
            #                 nodes[i_near] = new_node
            #             end
            #         end
            #     end
            # end
        end
    end
end

function get_best_path(problem::P, nodes::Vector{Node{T}}, x_goal; rev = true) where {T, P <: AbstractProblem{T}}
    # for each node, try to connect it to the goal
    best_cost = Inf
    best_node = -1
    I_near = near(problem, nodes, x_goal)
    for i in I_near
        node = nodes[i]
        node_cost = cost(i, nodes)
        incremental_cost = path_cost(problem, node.state, x_goal)

        if  node_cost + incremental_cost < best_cost
            # check if the path was collision free
            if collision_free(problem, node.state, x_goal)
                best_cost = node_cost + incremental_cost
                best_node = i
            end
        end
    end

    if best_node == -1
        return best_cost, [x_goal]
    else
        node = nodes[best_node]
        path = [x_goal, node.state]
        while node.parent_index != 0
            node = nodes[node.parent_index]
            push!(path, node.state)
        end
        return best_cost, rev ? reverse(path) : path
    end
end

# the following functions need to be defined for your problem
function sample_free(problem::P) where {T, P <: AbstractProblem{T}}
    throw(MethodError(sample_free, (problem, )))
end

function nearest(problem::P, nodes, x_rand) where {T, P<: AbstractProblem{T}}
    throw(MethodError(nearest, (problem, nodes, x_rand)))
end

function near(problem::P, nodes,  x_new) where {T, P<: AbstractProblem{T}}
    throw(MethodError(near, (problem, nodes, x_new)))
end

function steer(problem::P, x_nearest, x_rand) where {T, P<: AbstractProblem{T}}
    throw(MethodError(steer, (problem, x_nearest, x_rand)))
end

function collision_free(problem::P, x_nearest, x_new) where {T, P<: AbstractProblem{T}}
    throw(MethodError(collision_free, (problem, x_nearest, x_new)))
end

function path_cost(problem::P, x_near, x_new) where {T, P <: AbstractProblem{T}}
    throw(MethodError(path_cost, (problem, x_near, x_new)))
end


end





