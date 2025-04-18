using .ParetoRRTStar
using Plots
using Dubins
using StaticArrays

function Plots.plot!(path::DubinsPath, N = 50; kwargs...)
    errcode, samples = Dubins.dubins_path_sample_many(path, Dubins.dubins_path_length(path)/N )
    @assert errcode == Dubins.EDUBOK

    # push the final state too 
    _, pt = Dubins.dubins_path_endpoint(path)
    push!(samples, pt)

    plot!([s[1] for s in samples], [s[2] for s in samples]; kwargs...)
end

function Plots.plot!(nodes::Vector{ParetoRRTStar.Node{T,N,F}}, r::F; kwargs...) where {T,N,F}
    px = [n.state[1] for n in nodes]
    py = [n.state[2] for n in nodes]

    for node in nodes
        for path in node.pareto_paths 
            parent = path.path[end]
            q0 = nodes[parent].state
            q1 = node.state
        
            err, path = dubins_shortest_path(q0, q1, r, 1e-3)
          
            plot!(path; kwargs...)
        end
    end
    # scatter!(px, py, label=false)

end

function plot_pareto_front(nodes::Vector{ParetoRRTStar.Node{T,N,F}}, index::Int64, r::F; kwargs...) where {T,N,F}
    pareto_paths = nodes[index].pareto_paths
    for pareto_path in pareto_paths
        q0 = nodes[index].state
        for i in reverse(pareto_path.path)
            q1 = nodes[i].state
            err, path = dubins_shortest_path(q1, q0, r, 1e-3)
            q0 = q1
            plot!(path; kwargs...)
        end
    end
end