using Vec
using CollisionChecking
using LinearAlgebra
using Plots


function seg_intersection(L1::LineSegment,L2::LineSegment)
    A = [Vector(L1.B-L1.A) Vector(L2.A-L2.B)]
    b = L2.A - L1.A
    if rank(A) == 2
        # (0 < x[1] < 1) => intersection occurs within L1
        # (0 < x[2] < 1) => intersection occurs within L2
        x = inv(A)*b
        return x
    end
    return [NaN,NaN]
end

function offset_polygon(p::ConvexPolygon,d_offset)
    if d_offset < 0.0
        p_edges = [get_edge(p,i) for i in 1:p.npts]
        bisectors = []
        angles = []
        dists = []
        for j in 1:p.npts
            i = (j == 1) ? p.npts : j-1
            k = (j == p.npts) ? 1 : j+1
            # i, j, k
            v1 = p.pts[i] - p.pts[j]
            v1 = v1 / norm(v1)
            v2 = p.pts[k] - p.pts[j]
            v2 = v2 / norm(v2)
            # i = (j == 1) ? p.npts : j-1
            # e1 = p_edges[i]
            # e2 = p_edges[j]
            # v1 = (e1.A - e1.B) / norm(e1.A-e1.B)
            # v2 = (e2.B - e2.A) / norm(e2.B - e2.A)
            v = (v1 + v2) / norm(v1 + v2)
            b = LineSegment(p.pts[j],p.pts[j]+v)
            push!(bisectors,b)
        end
        remaining_edges = []
        for j in 1:p.npts
            k = (j == p.npts) ? 1 : j+1
            x = seg_intersection(bisectors[j],bisectors[k])
            if maximum(x) > -d_offset
                ed = get_edge(p,j)
                v = rot(ed.B-ed.A, -π/2.0)
                v = v / norm(v)
                e = ed + (d_offset*v) # offset
                # e = LineSegment(e.A,e.A + (e.B-e.A)/norm(e.B-e.A)) # normalize
                push!(remaining_edges,e)
            end
        end
        pts = Vector{VecE2}()
        for (j,e1) in enumerate(remaining_edges)
            k = (j == length(remaining_edges)) ? 1 : j+1
            e2 = remaining_edges[k]
            x = seg_intersection(e1,e2)
            push!(pts,e1.A+(e1.B-e1.A)*x[1])
        end
        return ConvexPolygon(pts)
    else
        pts = Vector{VecE2}()
        for j in 1:p.npts
            i = (j == 1) ? p.npts : j-1
            k = (j == p.npts) ? 1 : j+1
            # i, j, k
            v1 = p.pts[i] - p.pts[j]
            v1 = v1 / norm(v1)
            v2 = p.pts[k] - p.pts[j]
            v2 = v2 / norm(v2)
            v = (v1 + v2) / norm(v1 + v2)
            # b = LineSegment(p.pts[j],p.pts[j]+v)
            # push!(bisectors,b)
            dv = rot(v1,-atan(v2.y,v2.x))
            θ = (π - atan(dv.y,dv.x))/2.0
            # push!(angles,θ)
            d = 2*d_offset*sin(θ)
            push!(pts,p.pts[j]-d*v)
        end
        return ConvexPolygon(pts)
    end
end

function plot_polygon!(plt,p::ConvexPolygon)
    if length(p.pts) > 0
        x = [v.x for v in p.pts]
        push!(x,p.pts[1].x)
        y = [v.y for v in p.pts]
        push!(y,p.pts[1].y)
        plot!(plt,x,y)
    end
end

function offset_polygon(p::ConvexPolygon,d_offset)
    p_edges = [get_edge(p,i) for i in 1:p.npts]
    n_edges = []
    for e in p_edges
        v = rot((e.B-e.A) / norm(e.B-e.A), -π/2.0)
        push!(n_edges,e+(v*d_offset))
    end
    pts = Vector{VecE2}()
    if d_offset < 0.0
        left_neighbors = [-1 for e in n_edges]
        right_neighbors = [-1 for e in n_edges]
        left_dists = [0.0 for e in n_edges]
        right_dists = [1.0 for e in n_edges]
        for (i,e1) in enumerate(n_edges)
            # for (j,e2) in enumerate(n_edges)
            for j in cyclic_shift!(collect(1:length(n_edges)),i)
                e2 = n_edges[j]
                if i != j
                    # θ = mod(atan(rot(e2.B-e2.A,-atan(e1.B-e1.A))),2π)
                    θ = mod(mod(atan(e2.B-e2.A),2π) - mod(atan(e1.B-e1.A),2π),2π)
                    if θ < π
                        x = seg_intersection(e1,e2)
                        if !any(isnan,x)
                            if (x[1] < right_dists[i]) && (θ <= π) # (x[2] < 1.0)
                                right_dists[i] = x[1]
                                right_neighbors[i] = j
                            end
                            if (x[2] > left_dists[j]) && (θ <= π) # (x[1] >= 1.0)
                                left_dists[j] = x[2]
                                left_neighbors[j] = i
                            end
                        end
                    end
                end
            end
        end
        idxs = [i for i in 1:length(n_edges) if left_dists[i] < right_dists[i]]
        for j in idxs
            e = n_edges[j]
            push!(pts,e.A+(e.B-e.A)*right_dists[j])
        end
        return ConvexPolygon(pts)
    else
        for (i,e1) in enumerate(n_edges)
            j = (i == length(p_edges)) ? 1 : i+1
            e2 = n_edges[j]
            x = seg_intersection(e1,e2)
            push!(pts,e1.A+(e1.B-e1.A)*x[1])
        end
        return ConvexPolygon(pts)
    end
end

# function offset_polygon3(p::ConvexPolygon,d_offset)
#     if d_offset < 0.0
#         p_edges = [get_edge(p,i) for i in 1:p.npts]
#         bisectors = []
#         angles = []
#         dists = []
#         for j in 1:p.npts
#             i = (j == 1) ? p.npts : j-1
#             e1 = p_edges[i]
#             e2 = p_edges[j]
#             v1 = (e1.A - e1.B) / norm(e1.A-e1.B)
#             v2 = (e2.B - e2.A) / norm(e2.B - e2.A)
#             v = (v1 + v2) / norm(v1 + v2)
#             b = LineSegment(p.pts[j],p.pts[j]+v)
#             push!(bisectors,b)
#             dv = rot(v1,-atan(v2.y,v2.x))
#             θ = (π - atan(dv.y,dv.x))/2.0
#             push!(angles,θ)
#             d = 2*d_offset*sin(θ)
#         end
#         new_p = ConvexPolygon([p.pts[j]+dists[j]*bisectors[j] for j in 1:length(dists)])
#         n_edges = [get_edge(new_p,j) for j in 1:new_p.npts]
#         remaining_edges = []
#         for (j,n_e) in enumerate(n_edges)
#             p_e = p_edges[j]
#             if dot(n_e.B-n_e.A,p_e.B-p_e.A) > 0
#                 push!(remaining_edges,j)
#             end
#         end
#         n_bisectors = []
#         for (idx,j) in enumerate(remaining_edges)
#             i_idx = (idx == 1) ? length(remaining_edges) : idx-1
#             i = remaining_edges[i_idx]
#             if (j-i) == 1 || (i == length(p_edges) && j == 1)
#                 push!(n_bisectors,bisectors[j])
#             else
#                 e1 = p_edges[i]
#                 e2 = p_edges[j]
#                 v1 = e1.A - e1.B
#                 v1 = v1 / norm(v1)
#                 v2 = e2.B - e2.A
#                 v2 = v2 / norm(v2)
#                 v = (v1 + v2) / norm(v1 + v2)
#                 b = LineSegment(p.pts[j],p.pts[j]+v)
#                 push!(n_bisectors,b)
#             end
#         end
#         for j in 1:p.npts
#             k = (j == p.npts) ? 1 : j+1
#             x = seg_intersection(bisectors[j],bisectors[k])
#             if maximum(x) > -d_offset
#                 ed = get_edge(p,j)
#                 v = rot(ed.B-ed.A, -π/2.0)
#                 v = v / norm(v)
#                 e = ed + (d_offset*v) # offset
#                 # e = LineSegment(e.A,e.A + (e.B-e.A)/norm(e.B-e.A)) # normalize
#                 push!(remaining_edges,e)
#             end
#         end
#         pts = Vector{VecE2}()
#         for (j,e1) in enumerate(remaining_edges)
#             k = (j == length(remaining_edges)) ? 1 : j+1
#             e2 = remaining_edges[k]
#             x = seg_intersection(e1,e2)
#             push!(pts,e1.A+(e1.B-e1.A)*x[1])
#         end
#         return ConvexPolygon(pts)
#     else
#         pts = Vector{VecE2}()
#         for j in 1:p.npts
#             i = (j == 1) ? p.npts : j-1
#             k = (j == p.npts) ? 1 : j+1
#             # i, j, k
#             v1 = p.pts[i] - p.pts[j]
#             v1 = v1 / norm(v1)
#             v2 = p.pts[k] - p.pts[j]
#             v2 = v2 / norm(v2)
#             v = (v1 + v2) / norm(v1 + v2)
#             # b = LineSegment(p.pts[j],p.pts[j]+v)
#             # push!(bisectors,b)
#             dv = rot(v1,-atan(v2.y,v2.x))
#             θ = (π - atan(dv.y,dv.x))/2.0
#             # push!(angles,θ)
#             d = 2*d_offset*sin(θ)
#             push!(pts,p.pts[j]-d*v)
#         end
#         return ConvexPolygon(pts)
#     end
# end
#
# function offset_polygon2(p::ConvexPolygon,d_offset)
#     p_edges = [get_edge(p,i) for i in 1:p.npts]
#     p_bisectors = []
#     for j in 1:p.npts
#         i = (j == 1) ? p.npts : j-1
#         e1 = p_edges[i]
#         e2 = p_edges[j]
#         # k = (j == p.npts) ? 1 : j+1
#         # i, j, k
#         v1 = e1.A - e1.B
#         v1 = v1 / norm(v1)
#         v2 = e2.B - e2.A
#         v2 = v2 / norm(v2)
#         v = (v1 + v2) / norm(v1 + v2)
#         b = LineSegment(p.pts[j],p.pts[j]+v)
#         push!(p_bisectors,b)
#     end
#     L = length(p_edges)
#     n_edges = []
#     dropped_edges = []
#     for j in 1:length(p_edges)
#         k = (j == p.npts) ? 1 : j+1
#         x = seg_intersection(bisectors[j],bisectors[k])
#         if maximum(x) > -d_offset
#             ed = get_edge(p,j)
#             v = rot(ed.B-ed.A, -π/2.0)
#             v = v / norm(v)
#             e = ed + (d_offset*v) # offset
#             # e = LineSegment(e.A,e.A + (e.B-e.A)/norm(e.B-e.A)) # normalize
#             push!(n_edges,e)
#         else
#             push!(dropped_edges,j)
#         end
#     end
#     n_bisectors = []
#     for j in 1:length(p_bisectors)
#
#     end
#     if length(dropped_edges) > 0
#
#     end
#
#     if d_offset < 0.0
#         intersections = []
#         pts = Vector{VecE2}()
#         for (j,e1) in enumerate(p_edges)
#             xmin = Inf
#             idx = 0
#             for k in cyclic_shift!(collect(1:length(p_edges)),j)
#                 e2 = p_edges[k]
#                 x = seg_intersection(e1,e2)
#                 if x[1] < xmin
#                     xmin = x[1]
#                     idx = k
#                 else
#                     break
#                 end
#             end
#             push!(pts,e1.A+(e1.B-e1.A)*x[1])
#         end
#         return ConvexPolygon(pts)
#     else
#         pts = Vector{VecE2}()
#         for j in 1:p.npts
#             i = (j == 1) ? p.npts : j-1
#             k = (j == p.npts) ? 1 : j+1
#             # i, j, k
#             v1 = p.pts[i] - p.pts[j]
#             v1 = v1 / norm(v1)
#             v2 = p.pts[k] - p.pts[j]
#             v2 = v2 / norm(v2)
#             v = (v1 + v2) / norm(v1 + v2)
#             # b = LineSegment(p.pts[j],p.pts[j]+v)
#             # push!(bisectors,b)
#             dv = rot(v1,-atan(v2.y,v2.x))
#             θ = (π - atan(dv.y,dv.x))/2.0
#             # push!(angles,θ)
#             d = 2*d_offset*sin(θ)
#             push!(pts,p.pts[j]-d*v)
#         end
#         return ConvexPolygon(pts)
#     end
# end

# hexagon
# r = 1.0
# α = collect(range(0,stop=2π,length=7))[1:end-1]
# p = ConvexPolygon([VecE2(r*cos(a),r*sin(a)) for a in α])
# sort_pts!(p)
# irregular hexagon
r = 1.0
p = ConvexPolygon([VecE2(r*cos(a),r*sin(a)) for a in mod.(randn(6)*2π,2π)])
sort_pts!(p)
d_offset = 0.5
# pr = offset_polygon2(p,d_offset)
plt = plot(aspect_ratio=:equal)
plot_polygon!(plt,p)
pr = offset_polygon(p,d_offset)
sort_pts!(pr)
# plot!(plt)
plot_polygon!(plt,pr)
