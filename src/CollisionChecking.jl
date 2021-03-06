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
    get_width,
    get_height,
    get_origin,

    cyclic_shift_left!,
    cyclic_shift!,
    get_signed_area,
    get_edge,
    get_center,
    shift!,
    rotate!,
    mirror!,
    sort_pts!,
    remove_colinear_edges,

    penetration,
    check_collision,
    get_displacement,
    get_closest_point,
    seg_intersection,
    minkowski_sum!,
    minkowski_difference!,
    offset_polygon,
    nearest_free_space,
    nearest_in_bounds_space,

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
    c::VecE2{T} = VecE2(NaN,NaN)
    r::T        = NaN
end
# Circle(c::VecE2,r) = Circle(c.x,c.y,r)
Circle(x,y,r) = Circle(VecE2(x,y),r)
# get_center(c::Circle) = VecE2(c.x,c.y)
get_center(c::Circle) = c.c
@with_kw struct Rectangle{T} <: Shape
    pt1::VecE2{T} = VecE2(NaN,NaN)
    pt2::VecE2{T} = VecE2(NaN,NaN)
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
        VecE2(minimum(pt.x for pt in polygon.pts),minimum(pt.y for pt in polygon.pts)),
        VecE2(maximum(pt.x for pt in polygon.pts),maximum(pt.y for pt in polygon.pts))
    )
end
function Rectangle(circle::Circle)
    Rectangle(
        VecE2(get_center(circle).x - circle.r,get_center(circle).x + circle.r),
        VecE2(get_center(circle).y - circle.r,get_center(circle).y + circle.r)
    )
end

get_width(rect::Rectangle) = rect.pt2.x - rect.pt1.x
get_height(rect::Rectangle) = rect.pt2.y - rect.pt1.y
get_origin(rect::Rectangle) = rect.pt1
get_width(shape::Union{Circle,ConvexPolygon}) = get_width(Rectangle(shape))
get_height(shape::Union{Circle,ConvexPolygon}) = get_height(Rectangle(shape))
# Rectangle(rect::Rectangle) = Rectangle(rect.pt1,rect.pt2)
ConvexPolygon(poly::ConvexPolygon) = ConvexPolygon(poly.pts)

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
get_center(poly::ConvexPolygon) = sum(poly.pts) / poly.npts
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
function sort_pts!(poly::ConvexPolygon, npts::Int=poly.npts)
    """
    sorts points by polar angle
    """
    c = get_center(poly)
    sort!(poly.pts,by=pt->mod(atan(pt-c)+π/2,2π))
    poly
end
function remove_colinear_edges(poly::ConvexPolygon)
   to_remove = Set{Int}()
   for i in 1:poly.npts
       j = (i == poly.npts) ? 1 : i+1
       e1 = get_edge(poly,i)
       e2 = get_edge(poly,j)
       if parallel(e1,e2)
           push!(to_remove,j)
       end
   end
   idxs = setdiff(Set{Int}(collect(1:poly.npts)),to_remove)
   ConvexPolygon(poly.pts[sort(collect(idxs))])
end
function Base.empty!(poly::ConvexPolygon)
    poly.npts = 0
    poly.pts = Vector{VecE2}()
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
function Base.in(v::VecE2, polygons::Vector{P} where {P <: ConvexPolygon})
    for p in polygons
        if v in p
            return true
        end
    end
    return false
end

function Base.show(io::IO, poly::ConvexPolygon)
    @printf(io, "ConvexPolygon: len %d (max %d pts)\n", length(poly), length(poly.pts))
    for i in 1 : length(poly)
        print(io, "\t"); show(io, poly.pts[i])
        print(io, "\n")
    end
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
function Vec.intersects(seg::LineSegment,objects::Vector{P} where {P <: ConvexPolygon})
    for o in objects
        if intersects(seg,o)
            return true
        end
    end
    return false
end
Vec.intersects(v1::VecE2,v2::VecE2,shape) = Vec.intersects(LineSegment(v1,v2),shape)

function penetration(v::VecE2,polygon::ConvexPolygon)
    previous_side = 0
    d_max = Inf
    for i in 1 : length(polygon)
        seg = get_edge(polygon, i)
        affine_segment = seg.B - seg.A
        affine_point = v - seg.A
        current_side = sign(cross(affine_point,affine_segment)) # sign indicates side
        if current_side > 0 # outside
            return 0.0
        else
            d_max = min(get_distance(seg, v),d_max)
        end
    end
    return d_max
end
function penetration(pt::VecE2,polygons::Vector{P} where {P <: ConvexPolygon})
    for p in polygons
        if pt in p
            return penetration(pt,p)
        end
    end
    return 0.0
end
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
        if d > 0.0 # Changed this from >= 0.0. Not sure if it will cause major trouble...
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
        q = get_center(circle) - pt1
        v = (pt2-pt1)/norm(pt2-pt1)
        d = sign(cross([q;0],[v;0])[end])*norm(q - v*dot(q,v))
        if d >= circle.r
            return false
        end
    end
    return true
end
function check_collision(polygon1::ConvexPolygon,polygon2::ConvexPolygon)
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
    for i in 1:polygon1.npts
        e1 = get_edge(polygon1, i)
        for j in 1:polygon2.npts
            e2 = get_edge(polygon2,j)
            if intersects(e1,e2)
                return true
            end
        end
    end
    return false
end
check_collision(rect::Rectangle,polygon::ConvexPolygon) = check_collision(ConvexPolygon(rect),polygon)
function check_collision(pt::VecE2,circle::Circle)
    norm(get_center(circle)-pt) < circle.r
end
function check_collision(circle1::Circle,circle2::Circle)
    """ Between two circles """
    if norm(get_center(circle1)-get_center(circle2)) < (circle1.r + circle2.r)/2.0
        return true
    else
        return false
    end
end
check_collision(polygon::ConvexPolygon,circle::Circle) = check_collision(circle,polygon)
function check_collision(rect::Rectangle,circle::Circle)
    """ Between a circle and a rectangular obstacle """
    if (get_center(circle).x + circle.r) > min(rect.pt1.x,rect.pt2.x) && (get_center(circle).x - circle.r) < max(rect.pt1.x,rect.pt2.x)
        if (get_center(circle).y + circle.r) > min(rect.pt1.y,rect.pt2.y) && (get_center(circle).y - circle.r) < max(rect.pt1.y,rect.pt2.y)
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
function check_collision(shape1::T where {T<:Union{Shape,VecE2}},shapes::Vector{S} where {S<:Shape})
    for s in shapes
        if check_collision(shape1,s)
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
        q = get_center(circle) - pt1
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
    cpt = get_center(circle)
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
    return norm(get_center(c1)-get_center(c2)) + c1.r < c2.r
end
function Base.in(circle::Circle,rect::Rectangle)
    if min(rect.pt1.x,rect.pt2.x) + circle.r <= get_center(circle).x <= max(rect.pt1.x,rect.pt2.x) - circle.r
        if min(rect.pt1.y,rect.pt2.y) + circle.r <= get_center(circle).y <= max(rect.pt1.y,rect.pt2.y) - circle.r
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

function get_closest_point(P::VecE2, seg::LineSegment)
    """
    returns closest point on `seg` to `P`
    """
    ab = seg.B - seg.A
    pb = P - seg.A

    denom = normsquared(ab)
    if denom == 0.0
        return VecE2(0.0,0.0)
    end

    r = (ab⋅pb)/denom

    if r ≤ 0.0
        return seg.A
    elseif r ≥ 1.0
        return seg.B
    else
        return (seg.A + r*ab)
    end
end
function get_displacement(P::VecE2, seg::LineSegment)
    """
    returns vector from `P` to the closest point on `seg`
    """
    return get_closest_point(P,seg) - P
end
function get_displacement(P::VecE2, poly::ConvexPolygon;solid::Bool=true)
    if solid && P ∈ poly
        return VecE2(0.0,0.0)
    else
        Δmin = Inf
        Δv = VecE2(NaN,NaN)
        for i in 1:poly.npts
            seg = get_edge(poly,i)
            v = get_displacement(P,seg)
            if norm(v) < Δmin
                Δv = v
                Δmin = norm(v)
            end
        end
        return Δv
    end
end
function get_displacement(P::VecE2, circle::Circle;solid::Bool=true)
    if P ∈ circle
        return VecE2(0.0,0.0)
    else
        c = get_center(circle) - P
        return c * (1.0 - circle.r/norm(c))
    end
end
function get_displacement(c1::Circle, c2::Circle;solid::Bool=true)
    get_displacement(get_center(c1),Circle(get_center(c2),c1.r+c2.r))
end
function get_displacement(poly::ConvexPolygon,P::VecE2;solid::Bool=true)
    return -get_displacement(P,poly)
end

function Vec.get_distance(v1::T where T <: Union{VecE,VecSE}, v2::S where S <: Union{VecE,VecSE})
    norm(VecE2(v2.x-v1.x,v2.y-v1.y))
end
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
        push!(retval,P.pts[i] + (v + get_center(Q)))
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

function nearest_free_space(circle::Circle,rect::Rectangle)
    """
        get vector to nearest free space between circular agent and rectangular obstacle
    """
    Δx = 0.0
    Δy = 0.0
    rcenter = VecE2((rect.pt1.x+rect.pt2.x)/2,(rect.pt1.y+rect.pt2.y)/2)
    rwidth = VecE2(abs(rect.pt2.x-rect.pt1.x),abs(rect.pt2.y-rect.pt1.y))
    if abs(get_center(circle).x - rcenter[1]) <= circle.r+rwidth[1]/2
        if abs(get_center(circle).y - rcenter[2]) <= circle.r+rwidth[2]/2
            targetX = rcenter[1]
            targetY = rcenter[2]
            if get_center(circle).x < rcenter[1]
                targetX = (rcenter[1]-(rwidth[1]/2+circle.r))
            else
                targetX = (rcenter[1]+(rwidth[1]/2+circle.r))
            end
            if get_center(circle).y < rcenter[2]
                targetY = (rcenter[2]-(rwidth[2]/2+circle.r))
            else
                targetY = (rcenter[2]+(rwidth[2]/2+circle.r))
            end
            Δx = targetX - get_center(circle).x
            Δy = targetY - get_center(circle).y
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
        q = get_center(circle) - pt1
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
    p = offset_polygon(polygon,-circle.r)
    get_displacement(get_center(circle),p)
    # if c ∈ p
    #     return VecE2(0.0,0.0)
    # end
    # Δmin = VecE2(Inf,Inf)
    # n = length(polygon.pts)
    # for i in 1:n
    #     pt1 = polygon.pts[i]
    #     if i < n
    #         pt2 = polygon.pts[i+1]
    #     else
    #         pt2 = polygon.pts[1]
    #     end
    #     # projection
    #     q = VecE2(circle.x,circle.y) - pt1
    #     v = (pt2-pt1)/norm(pt2-pt1)
    #     b = v*dot(q,v) # footpoint
    #     d = sign(cross(q,v)) # direction (outward is positive)
    #     dv = q - b
    #     if d < 0 && norm(dv) > circle.r
    #         return VecE2(0.0,0.0)
    #     else
    #         target = circle.r * -d*dv/norm(dv) # pointing inside
    #         Δ = target - dv
    #         if norm(Δ) < norm(Δmin)
    #             Δmin = VecE2(Δ)
    #         end
    #     end
    # end
    # return Δmin
end
function nearest_in_bounds_space(circle::Circle,rect::Rectangle)
    nearest_in_bounds_space(circle, ConvexPolygon(rect))
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


######################################

end # module
