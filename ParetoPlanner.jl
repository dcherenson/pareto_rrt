using .ParetoRRTStar
using StaticArrays
using .World
using Dubins

# get just the distance between two nodes
function dubins_distance(q1, q2, turning_radius=0.1)

    e, p = Dubins.dubins_shortest_path(q1, q2, turning_radius, 1e-3)
    @assert e == Dubins.EDUBOK

    return Dubins.dubins_path_length(p)
end

@kwdef struct ParetoPlannerProblem{F, VW} <: ParetoRRTStar.AbstractProblem{SVector{3,F}}
    ϵ::SVector{2,F} = SVector(0.5, 0.5) # cost epsilon
    domain::Tuple{SVector{3,F}, SVector{3,F}} # rectangle defined by opposite corners
    near_radius::F = 30.0
    unsafe_zones::VW
    turning_radius::F = 10.0
    max_velocity::F = 10.0
    risk_zone::SVector{3,F} # x,y,radius
end

struct ParetoPlanner{TP, TN}
    prob::TP # ParetoPlanner{Float64}
    nodes::TN # Vector{ParetoRRTStar.Node{SVector{3,Float64}}}
end

function sample_domain(P::ParetoPlannerProblem)
    v = @SVector rand(3)
    return P.domain[1] + (P.domain[2] - P.domain[1]) .* v
end

function ParetoRRTStar.sample_free(problem::ParetoPlannerProblem)
    while true
        q = sample_domain(problem)
        # check if q is unsafe
        unsafe = World.is_unsafe(problem.unsafe_zones, q)
        if !unsafe
            return q
        end
    end
end

function ParetoRRTStar.nearest(P::ParetoPlannerProblem, nodes, x_rand)
    best_dist = Inf
    best_ind = -1

    for i=1:length(nodes)
        if (nodes[i].state[1] - x_rand[1])^2 + (nodes[i].state[2] - x_rand[2])^2 > best_dist^2
            continue
        end
        d = dubins_distance(nodes[i].state, x_rand, P.turning_radius)
        if d < best_dist
            best_dist = d
            best_ind = i
        end
    end
    return best_ind
end

function ParetoRRTStar.near(P::ParetoPlannerProblem, nodes, x_new)
    ind = Int[]
    # sizehint!(ind, length(nodes)/2)
    for i=1:length(nodes)
        if (nodes[i].state[1] - x_new[1])^2 + (nodes[i].state[2] - x_new[2])^2 > P.near_radius^2
            continue
        end
        d = dubins_distance(nodes[i].state, x_new, P.turning_radius)
        if d < P.near_radius
            push!(ind, i)
        end
    end
    return ind
end

function ParetoRRTStar.steer(problem::ParetoPlannerProblem, x_nearest, x_rand; max_travel_dist = 5.1)

    # first plan the full dubins path
    errcode, path = Dubins.dubins_shortest_path(x_nearest, x_rand, problem.turning_radius, 1e-3)
    @assert errcode == Dubins.EDUBOK

    L = Dubins.dubins_path_length(path) 

    # get the node 20% of the waythrough
    s = min(L, max_travel_dist)

    errcode, x_new = Dubins.dubins_path_sample(path, s)
    @assert errcode == Dubins.EDUBOK
    
    return SVector{3}(x_new)
    
end

function ParetoRRTStar.path_cost_collision_free(problem::ParetoPlannerProblem, x_nearest, x_new; step_size=0.1)

    # create the path
    errcode, path = Dubins.dubins_shortest_path(x_nearest, x_new, problem.turning_radius, 1e-3)
    @assert errcode == Dubins.EDUBOK

    L = Dubins.dubins_path_length(path)
    R = 0.0
    x = 0.0
    while x < L
        errcode, q = Dubins.dubins_path_sample(path, x)

        @assert errcode == Dubins.EDUBOK errcode
    
        # check each point for collision 
        if q[1] < problem.domain[1][1] || q[1] > problem.domain[2][1] || q[2] < problem.domain[1][2] || q[2] > problem.domain[2][2] || World.is_unsafe(problem.unsafe_zones, q)
            return SVector(0.0,0.0), false
        end
        r_sq = ((q[1] - problem.risk_zone[1])^2 + (q[2] - problem.risk_zone[2])^2)
        R += 10*exp(-r_sq / (4 * problem.risk_zone[3])) * step_size
        x += step_size
    end

    return SVector(L,R), true
end

# function ParetoRRTStar.path_cost(problem::ParetoPlannerProblem, x_near, x_new)

#     # create the path
#     errcode, path = Dubins.dubins_shortest_path(x_near, x_new, problem.turning_radius, 1e-3)
#     @assert errcode == Dubins.EDUBOK    

#     L = Dubins.dubins_path_length(path)

#     R = 0.0

#     while 

#     return SVector(L, C)
#     # return SVector(L)
# end