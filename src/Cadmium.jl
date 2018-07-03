module Cadmium

import Plots.Shape

function closest_key(d::Dict, new_key)
    if isempty(d)
        return new_key
    else
        for old_key in keys(d)
            if isapprox(old_key, new_key)
               return old_key
            end
        end
        return new_key
    end
end

function next_key(d::Dict, half_key, spent_half_key)
    for key in keys(d)
        if any(i->i==half_key,key) && !any(i->i==spent_half_key,key)
            println("This is reached")
            return key
        end
    end
    return []
end

function assemble_shape(dcf::DualContourField)
#     p_dict = Dict{Array{Float64,1}, Array{Float64,1}}()
    p_dict=Dict()
    index_dict = Dict()
    empty_array = []#Array{Array{Float64,1},1}()
    for leaf in allleaves(dcf.root)
        if !isempty(leaf.data.residual)
            clean_keys=copy(empty_array)
            for p in leaf.data.p
                dict_key = closest_key(p_dict, p)
                push!(clean_keys, dict_key)
                p_dict[dict_key]=push!(get(p_dict, dict_key, copy(empty_array)), getindex(leaf))
            end
            index_dict[getindex(leaf)] = clean_keys
        end
    end

    previous_p = first(keys(p_dict))
    previous_idcs =pop!(p_dict, previous_p)
    previous_idx = previous_idcs[1] #hardcoding for 2D, update with more sophisticated code for 3D

    if isfinite(previous_idcs[1].data.residual)
        verts=[previous_idcs[1].data.qef_min]
    else
        verts=[]
    end

    is_closed = false

    while !isempty(p_dict) && !is_closed
        next_idx = filter(i->i!=previous_idx, previous_idcs)
        if length(next_idx)>1
            min_idx=find(i->isfinite(i.data.residual), next_idx)
            if isempty(min_idx)
                min_idx=indmax([maximum(i.boundary.widths) for i in next_idx])
            else
                min_idx = first(min_idx)
            end
        else
            min_idx = 1
        end
        #allows faces to be traversed even if they don't contain qef points, like planes
        if isfinite(next_idx[min_idx].data.residual)
            push!(verts,next_idx[min_idx].data.qef_min)
        end
        next_ps = index_dict[next_idx[min_idx]] #may be replaceable by a cell reference
        next_p = filter(i->!isapprox(previous_p,i), next_ps)
        try
            previous_idcs = pop!(p_dict, next_p...)
        catch
            is_closed = true #replace with more rigorous logic; perhaps tracking the original p and seeing when it's called again
        end
        previous_idx = next_idx[min_idx]
        previous_p=first(next_p)
    end
    return hcat(verts...)
end

function Plots.Shape(dcf::DualContourField)
    verts = assemble_shape(dcf)
    Plots.Shape(verts[1,:], verts[2,:])
end

function Plots.Shape(
        csg::Surface,
        origin::AbstractArray,
        widths::AbstractArray,
        rtol=1e-3,
        atol=1e-3,
        surfcellmax=1e-1)

    dcf=DualContourField(csg, origin, widths, rtol, atol, surfcellmax)
    Plots.Shape(dcf)
end


end # module
