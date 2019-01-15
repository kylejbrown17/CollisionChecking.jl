module CollisionChecking

# Collision Detection routines based on Minkowski techniques.
# Borrows heavily from AutomotiveDrivingModels: see https://github.com/sisl/AutomotiveDrivingModels.jl/blob/master/src/2d/utils/minkowski.jl

using Parameters
using Printf
using Colors
using LinearAlgebra
using StaticArrays
using Vec

export
    Shape,
    ConvexPolygon,
    Circle,
    Rectangle,

    # has_intersection,
    check_collision,
    nearest_free_space,
    nearest_in_bounds_space,

    cyclic_shift_left!,
    cyclic_shift!,
    get_signed_area,
    get_edge,
    get_center,
    sort_pts!,
    ensure_pts_sorted_by_min_polar_angle!,
    shift!,
    rotate!,
    mirror!,
    minkowski_sum!,
    minkowski_difference!,
    is_colliding,
    get_collision_time

abstract type Shape end
@with_kw mutable struct ConvexPolygon{T} <: Shape
    pts::Vector{T} = Vector{VecE2{Float64}}() # ordered counterclockwise along polygon boundary s.t. first edge has minimum polar angle in [0,2π]
    npts::Int = length(pts) # number of pts that are currently used (since we preallocate a longer array)
end
ConvexPolygon(npts::Int) = ConvexPolygon(Array{VecE2}(undef, npts), 0)
ConvexPolygon(pts::Vector{T} where T <: VecE2) = ConvexPolygon(pts, length(pts))
@with_kw struct Circle{T} <: Shape
    x::T = NaN
    y::T = NaN
    r::T = NaN
end
@with_kw struct Rectangle{T} <: Shape
    pt1::VecE2{T}
    pt2::VecE2{T}
end
function ConvexPolygon(rect::Rectangle)
    ConvexPolygon([
        VecE2(min(rect.pt1.x,rect.pt2.x),min(rect.pt1.y,rect.pt2.y)),
        VecE2(max(rect.pt1.x,rect.pt2.x),min(rect.pt1.y,rect.pt2.y)),
        VecE2(max(rect.pt1.x,rect.pt2.x),max(rect.pt1.y,rect.pt2.y)),
        VecE2(min(rect.pt1.x,rect.pt2.x),max(rect.pt1.y,rect.pt2.y)),
    ])
end
function Rectangle(polygon::ConvexPolygon)
    Rectangle(
        pt1=VecE2(minimum(pt.x for pt in polygon.pts),minimum(pt.y for pt in polygon.pts)),
        pt2=VecE2(maximum(pt.x for pt in polygon.pts),maximum(pt.y for pt in polygon.pts))
    )
end

function Vec.intersects(seg::LineSegment,poly::ConvexPolygon)
    for i in 1:length(poly.pts)
        seg2 = get_edge(poly,i)
        if intersects(seg,seg2)
            return true
        end
    end
    return false
end
function Vec.intersects(seg::LineSegment,objects::Vector{ConvexPolygon})
    for o in objects
        if intersects(seg,o)
            return true
        end
    end
    return false
end
# function has_intersection(p1::T,p2::T,q1::T,q2::T) where T <: AbstractVector
#     """
#     Returns true if there is an intersection
#     between seg1 and seg2
#     """
#     # L1: (p1) -- (p2)
#     # L2: (q1) -- (q2)
#     A = [Vector(p2-p1) Vector(q1-q2)]
#     b = q1 - p1
#     if rank(A) == 2
#         x = inv(A)*b
#         s = x[1]
#         t = x[2]
#         if (0.0 <= s <= 1.0) && (0.0 <= t <= 1.0)
#             return true # intersection
#         end
#     else
#         # is q1 in p1 -- p2?
#         a = q1-p1
#         b = p2-p1
#         c = dot(a,b) / (norm(b)^2)
#         if (0.0 <= c <= 1.0) && (norm(c*b - a) == 0)
#             return true # colinear overlap
#         end
#         # is q2 in p1 -- p2?
#         a = q2-p1
#         b = p2-p1
#         c = dot(a,b) / (norm(b)^2)
#         if (0.0 <= c <= 1.0) && (norm(c*b - a) == 0)
#             return true # colinear overlap
#         end
#         # is p1 in q1 -- q2?
#         a = p1-q1
#         b = q2-q1
#         c = dot(a,b) / (norm(b)^2)
#         if (0.0 <= c <= 1.0) && (norm(c*b - a) == 0)
#             return true # colinear overlap
#         end
#     end
#     return false
# end
# function has_intersection(seg1::Tuple{T,T},seg2::Tuple{T,T}) where T <: AbstractVector
#     """
#     Returns true if there is an intersection
#     between seg1 and seg2
#     """
#     has_intersection(seg1[1],seg1[2],seg2[1],seg2[2])
# end
# function has_intersection(seg::Tuple{T,T},pts::Vector{S}) where {T <: AbstractVector,S <: AbstractVector}
#     """
#     Returns true if there is an intersection
#     between the line segment and the polygon
#     """
#     p1,p2 = seg
#     for (i,q1) in enumerate(pts)
#         if i < length(pts)
#             q2 = pts[i+1]
#             if has_intersection(p1,p2,q1,q2)
#                 return true
#             end
#         else
#             q2 = pts[1]
#             if has_intersection(p1,p2,q1,q2)
#                 return true
#             end
#         end
#     end
#     return false
# end
# function has_intersection(pt1::T,pt2::T,objects::Vector{ConvexPolygon}) where T <: AbstractVector
#     """
#     checks to see if there is a clear line of sight between two points through
#     an environment described by a collection of obstacles
#     """
#     for o in objects
#         if has_intersection((pt1,pt2),o.pts)
#             return true
#         end
#     end
#     return false
# end

function check_collision(pt::VecE2,polygon::ConvexPolygon)
    n = length(polygon.pts)
    for i in 1:n
        pt1 = polygon.pts[i]
        if i < n
            pt2 = polygon.pts[i+1]
        else
            pt2 = polygon.pts[1]
        end
        q = pt - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        d = sign(cross([q;0],[v;0])[end])
        if d > 0.0
            return false
        end
    end
    return true
end
check_collision(polygon::ConvexPolygon,pt::VecE2) = check_collision(pt,polygon)
function check_collision(circle::Circle,polygon::ConvexPolygon)
    n = length(polygon.pts)
    for i in 1:n
        pt1 = polygon.pts[i]
        if i < n
            pt2 = polygon.pts[i+1]
        else
            pt2 = polygon.pts[1]
        end
        q = VecE2(circle.x,circle.y) - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        d = sign(cross([q;0],[v;0])[end])*norm(q - v*dot(q,v))
        if d > circle.r
            return false
        end
    end
    return true
end
function check_collision(polygon1::ConvexPolygon,polygon2::ConvexPolygon)
    n1 = length(polygon1.pts)
    n2 = length(polygon2.pts)
    for pt in polygon1.pts
        if check_collision(pt,polygon2)
            return true
        end
    end
    for pt in polygon2.pts
        if check_collision(pt,polygon1)
            return true
        end
    end
    for i in 1:n1
        pt1 = polygon1.pts[i]
        if i < n1
            pt2 = polygon1.pts[i+1]
        else
            pt2 = polygon1.pts[1]
        end
        for j in 1:n2
            pt3 = polygon2.pts[j]
            if j < n2
                pt4 = polygon2.pts[j+1]
            else
                pt4 = polygon2.pts[1]
            end
            if has_intersection((pt1,pt2),(pt3,pt4))
                return true
            end
        end
    end
    return false
end
check_collision(rect::Rectangle,polygon::ConvexPolygon) = check_collision(ConvexPolygon(rect),polygon)
function check_collision(pt::VecE2,circle::Circle)
    norm(VecE2(circle.x,circle.y)-pt) < circle.r
end
function check_collision(circle1::Circle,circle2::Circle)
    """ Between two circles """
    if norm([circle1.x-circle2.x;circle1.y-circle2.y]) < (circle1.r + circle2.r)/2.0
        return true
    else
        return false
    end
end
check_collision(polygon::ConvexPolygon,circle::Circle) = check_collision(circle,polygon)
function check_collision(rect::Rectangle,circle::Circle)
    """ Between a circle and a rectangular obstacle """
    if (circle.x + circle.r) > min(rect.pt1.x,rect.pt2.x) && (circle.x - circle.r) < max(rect.pt1.x,rect.pt2.x)
        if (circle.y + circle.r) > min(rect.pt1.y,rect.pt2.y) && (circle.y - circle.r) < max(rect.pt1.y,rect.pt2.y)
            return true
        end
    end
    return false
end
function check_collision(pt::VecE2,rect::Rectangle)
    """ Between a circle and a rectangular obstacle """
    if  min(rect.pt1.x,rect.pt2.x) < pt.x < max(rect.pt1.x,rect.pt2.x)
        if min(rect.pt1.y,rect.pt2.y) < pt.y < max(rect.pt1.y,rect.pt2.y)
            return true
        end
    end
    return false
end
check_collision(circle::Circle,rect::Rectangle) = check_collision(rect,circle)
function check_collision(rect1::Rectangle,rect2::Rectangle)
    if max(rect1.x1,rect1.x2) > min(rect2.x1,rect2.x2) && min(rect1.x1,rect1.x2) < max(rect2.x1,rect2.x2)
        if max(rect1.y1,rect1.y2) > min(rect2.y1,rect2.y2) && min(rect1.y1,rect1.y2) < max(rect2.y1,rect2.y2)
            return true
        end
    end
    return false
end

function Base.in(pt::VecE2,object::T where T <: Shape)
    check_collision(pt,object)
end
function Base.in(circle::Circle,polygon::ConvexPolygon)
    n = length(polygon.pts)
    for i in 1:n
        pt1 = polygon.pts[i]
        if i < n
            pt2 = polygon.pts[i+1]
        else
            pt2 = polygon.pts[1]
        end
        q = VecE2(circle.x,circle.y) - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        d = sign(cross([q;0],[v;0])[end])*norm(q - v*dot(q,v))
        if d > -circle.r
            return false
        end
    end
    return true
end
function Base.in(polygon1::ConvexPolygon,polygon2::ConvexPolygon)
    for pt in polygon1.pts
        if !check_collision(pt,polygon2)
            return false
        end
    end
    return true
end
function Base.in(rect::Rectangle,polygon::ConvexPolygon)
    Base.in(ConvexPolygon(rect),polygon)
end
function Base.in(polygon::ConvexPolygon,circle::Circle)
    cpt = VecE2(circle.x,circle.y)
    for pt in polygon.pts
        if norm(cpt-pt) > circle.r
            return false
        end
    end
    return true
end
function Base.in(rect::Rectangle,circle::Circle)
    Base.in(ConvexPolygon(rect),circle)
end
function Base.in(c1::Circle,c2::Circle)
    return norm([c1.x,c1.y]-[c2.x,c2.y]) + c1.r < c2.r
end
function Base.in(circle::Circle,rect::Rectangle)
    if min(rect.pt1.x,rect.pt2.x) + circle.r <= circle.x <= max(rect.pt1.x,rect.pt2.x) - circle.r
        if min(rect.pt1.y,rect.pt2.y) + circle.r <= circle.y <= max(rect.pt1.y,rect.pt2.y) - circle.r
            return true
        end
    end
    return false
end
function Base.in(polygon::ConvexPolygon,rect::Rectangle)
    Base.in(polygon,ConvexPolygon(rect))
end
function Base.in(rect1::Rectangle,rect2::Rectangle)
    Base.in(ConvexPolygon(rect1),ConvexPolygon(rect2))
end

function nearest_free_space(circle::Circle,rect::Rectangle)
    """
        get vector to nearest free space between circular agent and rectangular obstacle
    """
    Δx = 0.0
    Δy = 0.0
    rcenter = VecE2((rect.pt1.x+rect.pt2.x)/2,(rect.pt1.y+rect.pt2.y)/2)
    rwidth = VecE2(abs(rect.pt2.x-rect.pt1.x),abs(rect.pt2.y-rect.pt1.y))
    if abs(circle.x - rcenter[1]) <= circle.r+rwidth[1]/2
        if abs(circle.y - rcenter[2]) <= circle.r+rwidth[2]/2
            targetX = rcenter[1]
            targetY = rcenter[2]
            if circle.x < rcenter[1]
                targetX = (rcenter[1]-(rwidth[1]/2+circle.r))
            else
                targetX = (rcenter[1]+(rwidth[1]/2+circle.r))
            end
            if circle.y < rcenter[2]
                targetY = (rcenter[2]-(rwidth[2]/2+circle.r))
            else
                targetY = (rcenter[2]+(rwidth[2]/2+circle.r))
            end
            Δx = targetX - circle.x
            Δy = targetY - circle.y
        end
    end
    if abs(Δx) < abs(Δy)
        return VecE2(Δx, 0.0)
    else
        return VecE2(0.0, Δy)
    end
end
function nearest_free_space(circle::Circle,polygon::ConvexPolygon)
    """
        get vector to nearest free space between circular agent and rectangular obstacle
    """
    Δmin = VecE2(Inf,Inf)
    n = length(polygon.pts)
    for i in 1:n
        pt1 = polygon.pts[i]
        if i < n
            pt2 = polygon.pts[i+1]
        else
            pt2 = polygon.pts[1]
        end
        # projection
        q = VecE2(circle.x,circle.y) - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        b = v*dot(q,v) # footpoint
        d = sign(cross(q,v)) # direction (outward is positive)
        dv = q - b
        # dv = sign(cross([q;0],[v;0])[end])*(q - b) # normal vector pointing outward from polygon face
        if d > 0 && norm(dv) > circle.r
            return VecE2(0.0,0.0)
        else
            target = circle.r * d*dv/norm(dv)
            Δ = target - dv
            if norm(Δ) < norm(Δmin)
                Δmin = VecE2(Δ)
            end
        end
    end
    return Δmin
end
function nearest_in_bounds_space(circle::Circle,polygon::ConvexPolygon)
    """
        get vector to nearest free space between circular agent and rectangular obstacle
    """
    Δmin = VecE2(Inf,Inf)
    n = length(polygon.pts)
    for i in 1:n
        pt1 = polygon.pts[i]
        if i < n
            pt2 = polygon.pts[i+1]
        else
            pt2 = polygon.pts[1]
        end
        # projection
        q = VecE2(circle.x,circle.y) - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        b = v*dot(q,v) # footpoint
        d = sign(cross(q,v)) # direction (outward is positive)
        dv = q - b
        if d < 0 && norm(dv) > circle.r
            return VecE2(0.0,0.0)
        else
            target = circle.r * -d*dv/norm(dv) # pointing inside
            Δ = target - dv
            if norm(Δ) < norm(Δmin)
                Δmin = VecE2(Δ)
            end
        end
    end
    return Δmin
end

######################################

function cyclic_shift_left!(arr::AbstractVector, d::Int, n::Int=length(arr))
    #=
    Perform a cyclic rotation of the elements in the array, in place
        d = amount to shift
        n = length of array (if we only want to work with first n elements)
    =#
    for i in 1 : gcd(d, n)
    # move i-th values of blocks

        temp = arr[i]
        j = i
        while true
            k = j + d
            if k > n
                k = k - n
            end
            if k == i
                break
            end
            arr[j] = arr[k]
            j = k
        end
        arr[j] = temp
    end
    arr
end
function cyclic_shift!(arr::AbstractVector, d::Int)
    idxs = collect(1:length(arr))
    i = (d == length(arr)) ? 1 : d + 1
    idxs = [idxs[i:end];idxs[1:i-1]]
    arr = arr[idxs]
end

######################################

function get_signed_area(pts::Vector{T} where T <: VecE2, npts::Int=length(pts))

    # https://en.wikipedia.org/wiki/Shoelace_formula
    # sign of -1 means clockwise, sign of 1 means counterclockwise

    retval = pts[npts].x*pts[1].y - pts[1].x*pts[npts].y
    for i in 1 : npts-1
        retval += pts[i].x * pts[i+1].y
        retval -= pts[i+1].x * pts[i].y
    end

    retval / 2
end
function get_edge(pts::Vector{T} where T <: VecE2, i::Int, npts::Int=length(pts))
    a = pts[i]
    b = i+1 ≤ npts ? pts[i+1] : pts[1]
    LineSegment(a,b)
end

######################################

function Base.iterate(poly::ConvexPolygon, i::Int=1)
    if i > length(poly)
        return nothing
    end
    return (poly.pts[i], i+1)
end

Base.length(poly::ConvexPolygon) = poly.npts
Base.isempty(poly::ConvexPolygon) = poly.npts == 0
Base.copy(p::ConvexPolygon) = ConvexPolygon(copy(p.pts))
get_edge(P::ConvexPolygon, i::Int) = get_edge(P.pts, i, P.npts)
get_signed_area(poly::ConvexPolygon) = get_signed_area(poly.pts, poly.npts)
function Base.empty!(poly::ConvexPolygon)
    poly.npts = 0
    poly
end
function Base.copyto!(dest::ConvexPolygon, src::ConvexPolygon)
    dest.npts = src.npts
    copyto!(dest.pts, 1, src.pts, 1, src.npts)
    dest
end
function Base.push!(poly::ConvexPolygon, v::VecE2)
    push!(poly.pts, v)
    poly.npts = length(poly.pts)
    poly
end
function Base.in(v::VecE2, poly::ConvexPolygon)
    # DOES include pts on the physical boundary
    previous_side = 0
    for i in 1 : length(poly)
        seg = get_edge(poly, i)
        affine_segment = seg.B - seg.A
        affine_point = v - seg.A
        current_side = sign(cross(affine_point,affine_segment)) # sign indicates side
        if current_side > 0 # outside
            return false
        end
    end
    true
end
# Base.in(v::VecSE2{Float64}, poly::ConvexPolygon) = in(convert(VecE2, v), poly)
function Base.show(io::IO, poly::ConvexPolygon)
    @printf(io, "ConvexPolygon: len %d (max %d pts)\n", length(poly), length(poly.pts))
    for i in 1 : length(poly)
        print(io, "\t"); show(io, poly.pts[i])
        print(io, "\n")
    end
end

get_center(poly::ConvexPolygon) = sum(poly.pts) / poly.npts
function Vec.get_distance(poly::ConvexPolygon, v::VecE2; solid::Bool=true)
    if solid && in(v, poly)
        0.0
    else
        min_dist = Inf
        for i in 1 : length(poly)
            seg = get_edge(poly, i)
            min_dist = min(min_dist, get_distance(seg, v))
        end
        min_dist
    end
end
function Vec.get_distance(v::VecE2, poly::ConvexPolygon; solid::Bool=true)
    get_distance(poly,v;solid=solid)
end
function sort_pts!(poly::ConvexPolygon, npts::Int=poly.npts)
    """
    sorts points by polar angle
    """
    c = get_center(poly)
    sort!(poly.pts,by=pt->mod(atan(pt-c)+π/2,2π))
    poly
end
# function ensure_pts_sorted_by_min_polar_angle!(poly::ConvexPolygon, npts::Int=poly.npts)
#
#     @assert(npts ≥ 3)
#     sort_pts!(poly)
#     @assert(sign(get_signed_area(poly)) == 1) # must be counter-clockwise
#
#     # ensure that edges are sorted by minimum polar angle in [0,2π]
#     angle_start = Inf
#     index_start = -1
#     for i in 1 : npts
#         seg = get_edge(poly.pts, i, npts)
#
#         θ = mod(atan(seg.B - seg.A),2π)
#         if θ < angle_start
#             angle_start = θ
#             index_start = i
#         end
#     end
#
#     if index_start != 1
#         cyclic_shift_left!(poly.pts, index_start-1, npts)
#     end
#     poly
# end
function shift!(poly::ConvexPolygon, v::VecE2)
    for i in 1 : length(poly)
        poly.pts[i] += v
    end
    sort_pts!(poly)
    poly
end
function rotate!(poly::ConvexPolygon, θ::Float64)
    for i in 1 : length(poly)
        poly.pts[i] = Vec.rot(poly.pts[i], θ)
    end
    sort_pts!(poly)
    poly
end
function mirror!(poly::ConvexPolygon)
    for i in 1 : length(poly)
        poly.pts[i] = -poly.pts[i]
    end
    sort_pts!(poly)
    poly
end

function minkowski_sum!(P::ConvexPolygon, Q::ConvexPolygon, retval::ConvexPolygon=ConvexPolygon())
    #=
    For two convex polygons P and Q in the plane with m and n vertices, their Minkowski sum is a
    convex polygon with at most m + n vertices and may be computed in time O (m + n) by a very simple procedure,
    which may be informally described as follows.

    Assume that the edges of a polygon are given and the direction, say, counterclockwise, along the polygon boundary.
    Then it is easily seen that these edges of the convex polygon are ordered by polar angle.
    Let us merge the ordered sequences of the directed edges from P and Q into a single ordered sequence S.
    Imagine that these edges are solid arrows which can be moved freely while keeping them parallel to their original direction.
    Assemble these arrows in the order of the sequence S by attaching the tail of the next arrow to the head of the previous arrow.
    It turns out that the resulting polygonal chain will in fact be a convex polygon which is the Minkowski sum of P and Q.
    =#
    empty!(retval)
    index_P = 1
    index_Q = 1
    θp = get_polar_angle(get_edge(P, index_P))
    θq = get_polar_angle(get_edge(Q, index_Q))
    while index_P ≤ length(P) || index_Q ≤ length(Q)
        # select next edge with minimum polar angle
        if θp == θq
            seg_p = get_edge(P, index_P)
            seg_q = get_edge(Q, index_Q)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg_p.B - seg_p.A + seg_q.B - seg_q.A)
            index_P += 1
            θp = index_P ≤ length(P) ? get_polar_angle(get_edge(P, index_P)) : Inf
            index_Q += 1
            θq = index_Q ≤ length(Q) ? get_polar_angle(get_edge(Q, index_Q)) : Inf
        elseif θp ≤ θq
            seg = get_edge(P, index_P)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg.B - seg.A)
            index_P += 1
            θp = index_P ≤ length(P) ? get_polar_angle(get_edge(P, index_P)) : Inf
        else
            seg = get_edge(Q, index_Q)
            O = isempty(retval) ? P.pts[1] + Q.pts[1] : retval.pts[retval.npts]
            push!(retval, O + seg.B - seg.A)
            index_Q += 1
            θq = index_Q ≤ length(Q) ? get_polar_angle(get_edge(Q, index_Q)) : Inf
        end
    end
    sort_pts!(retval)
    retval
end
function minkowski_sum!(P::ConvexPolygon, Q::Circle, retval::ConvexPolygon=ConvexPolygon())
    """
    Not a true minkowski_sum, since it leaves the corners sharp
    """
    for i in 1:length(P)
        seg1 = get_edge(P,i)
        seg2 = i > 1 ? get_edge(P,i-1) : get_edge(P,length(P))
        seg2 = LineSegment(seg2.B,seg2.A) # reverse
        v1 = rot(normalize(seg1.B - seg1.A),-π/2)
        v2 = rot(normalize(seg2.B - seg2.A),π/2)
        v = (v1 + v2)/2
        a = norm(v)
        b = norm(v2-v1)/2
        c = b^2 / a
        v = (Q.r * (a+c)) * normalize(v)
        push!(retval,P.pts[i] + (v + VecE2(Q.x,Q.y)))
    end
    # sort_pts!(retval)
    retval
end
function minkowski_difference!(P::ConvexPolygon, Q::ConvexPolygon, retval::ConvexPolygon=ConvexPolygon())
    #=
    The minkowski difference is what you get by taking the minkowski sum of a shape and the mirror of another shape.
    So, your second shape gets flipped about the origin (all of its points are negated).

    The idea is that you do a binary operation on two shapes to get a new shape,
    and if the origin (the zero vector) is inside that shape, then they are colliding.
    =#
    minkowski_sum!(P, mirror!(Q),retval)
end
function minkowski_difference!(P::ConvexPolygon, Q::Circle, retval::ConvexPolygon=ConvexPolygon())
    minkowski_sum!(P,Q)
end

function is_colliding(P::ConvexPolygon, Q::ConvexPolygon)
    poly = minkowski_difference!(P, Q)
    in(VecE2(0.0,0.0), poly)
end
function Vec.get_distance(P::ConvexPolygon, Q::ConvexPolygon)
    poly = minkowski_difference!(P, Q)
    get_distance(VecE2(0,0), poly)
end
function Vec.get_distance(P::ConvexPolygon, C::Circle)
    poly = minkowski_difference!(P, C)
    get_distance(VecE2(0,0), poly)
end

# function nearest_in_bounds_space(circle::Circle,polygon::ConvexPolygon)
#     """
#         get vector to nearest free space between circular agent and rectangular obstacle
#     """
#     poly = copy(polygon)
#     i = 1
#     while i <= length(poly)
#         npts = length(poly)
#         seg = get_edge(poly,i)
#         v = rot(normalize(seg.B - seg.A),π/2)
#         plane = LineSegment(plane.A+v*circle.r,plane.B+v*circle.r)
#         idxs = collect(1:length(poly))
#         idxs = [idxs[mod(i+2,npts):end];idxs[1:mod(i+1,npts)]]
#         for j in idxs
#
#         end
#         pt_sides = [sign(cross(pt,seg.B-seg.A)) for pt in poly.pts]
#         valid_idxs = collect(1:length(poly.pts))[pt_sides .<= 0]
#         new_pt = intersection(plane,get_edge(poly,valid_idxs[1]-1))
#         prev_side = pt_sides[i]
#     end
# # end

######################################

end # module
