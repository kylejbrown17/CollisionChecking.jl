using CollisionChecking
using LinearAlgebra, Vec
using Test

################################################################################
############################## Collision Checking ##############################
################################################################################
# Circle on Circle
@test check_collision(Circle(0.0,0.0,1.0),Circle(0.5,0.0,1.0)) == true
@test check_collision(Circle(0.0,0.0,1.0),Circle(2.5,0.0,1.0)) == false
# ConvexPolygon and VecE2
p = ConvexPolygon([
    VecE2(0.0,0.0),
    VecE2(1.0,0.0),
    VecE2(1.0,1.0),
    VecE2(0.0,1.0),
    ])
@test check_collision(VecE2(0.5,0.5),p) == true
@test check_collision(p,VecE2(1.5,1.5)) == false
# ConvexPolygon and Circle
@test check_collision(Circle(1.25,0.5,1.0),p) == true
@test check_collision(p,Circle(2.25,0.5,1.0)) == false
@test check_collision(Circle(0.5,0.5,1.0),p) == true
# ConvexPolygon and ConvexPolygon
p1 = ConvexPolygon([
    VecE2(-2.0,0.0),
    VecE2(2.0,0.0),
    VecE2(2.0,1.0),
    VecE2(-2.0,1.0),
    ])
p2 = ConvexPolygon([
    VecE2(0.0,-2.0),
    VecE2(1.0,-2.0),
    VecE2(1.0,2.0),
    VecE2(0.0,2.0),
    ])
p3 = ConvexPolygon([
    VecE2(1.0,0.5),
    VecE2(3.0,0.5),
    VecE2(3.0,1.5),
    VecE2(1.0,1.5),
    ])
p4 = ConvexPolygon([
    VecE2(1.0,2.5),
    VecE2(3.0,2.5),
    VecE2(3.0,3.5),
    VecE2(1.0,3.5),
    ])
@test check_collision(p1,p2) == true
@test check_collision(p1,p3) == true
@test check_collision(p3,p4) == false

################################################################################
############################## Is Contained Within #############################
################################################################################
p = ConvexPolygon([
    VecE2(-2.0,-2.0),
    VecE2(2.0,-2.0),
    VecE2(2.0,2.0),
    VecE2(-2.0,2.0),
    ])
@test Base.in(Circle(0.0,0.0,1.0),p) == true
@test Base.in(Circle(0.0,0.0,10.0),p) == false
@test Base.in(Circle(0.0,0.0,1.0),Rectangle(p)) == true
@test Base.in(Circle(0.0,0.0,10.0),Rectangle(p)) == false
@test Base.in(p,Circle(0.0,0.0,10.0)) == true
@test Base.in(Rectangle(p),Circle(0.0,0.0,10.0)) == true
@test Base.in(p,Circle(0.0,0.0,1.0)) == false
@test Base.in(Rectangle(p),Circle(0.0,0.0,1.0)) == false

################################################################################
############################### Nearest Free Space #############################
################################################################################
p = ConvexPolygon([
    VecE2(-2.0,-2.0),
    VecE2(2.0,-2.0),
    VecE2(2.0,2.0),
    VecE2(-2.0,2.0),
    ])
@test norm(nearest_free_space(Circle(4.0,0.0,1.0),p) - [0.0;0.0]) < 0.000001
@test norm(nearest_free_space(Circle(4.0,0.0,1.0),Rectangle(p)) - [0.0;0.0]) < 0.000001
@test norm(nearest_free_space(Circle(2.1,0.0,1.0),p) - [0.9;0.0]) < 0.000001
@test norm(nearest_free_space(Circle(2.1,0.0,1.0),Rectangle(p)) - [0.9;0.0]) < 0.000001
@test norm(nearest_free_space(Circle(1.9,0.0,1.0),p) - [1.1;0.0]) < 0.000001
@test norm(nearest_free_space(Circle(1.9,0.0,1.0),Rectangle(p)) - [1.1;0.0]) < 0.000001


# @test norm(nearest_in_bounds_space(Circle(4.0,0.0,1.0),p) - [-0.0;0.0]) < 0.000001
# @test norm(nearest_in_bounds_space(Circle(4.0,0.0,1.0),Rectangle(p)) - [0.0;0.0]) < 0.000001
# @test norm(nearest_in_bounds_space(Circle(2.1,0.0,1.0),p) - [-1.1;0.0]) < 0.000001
# @test norm(nearest_in_bounds_space(Circle(2.1,0.0,1.0),Rectangle(p)) - [-1.1;0.0]) < 0.000001
# @test norm(nearest_in_bounds_space(Circle(1.9,0.0,1.0),p) - [-0.9;0.0]) < 0.000001
# @test norm(nearest_in_bounds_space(Circle(1.9,0.0,1.0),Rectangle(p)) - [-.9;0.0]) < 0.000001

α = collect(range(0,stop=.9*2π,length=10))
pts = [VecE2(cos(a),sin(a)) for a in α]
poly1 = ConvexPolygon(pts)
poly2 = ConvexPolygon(pts)
mirror!(poly2)
@test length(poly1) <= length(minkowski_sum!(poly1,poly2)) <= length(poly2)
@test length(poly1) <= length(minkowski_difference!(poly1,poly2)) <= length(poly2)
shift!(poly2,VecE2(2.0,0.0))
@test is_colliding(poly1,poly2) == false

p = ConvexPolygon([
    VecE2(-2.0,-2.0),
    VecE2(2.0,-2.0),
    VecE2(2.0,2.0),
    VecE2(-2.0,2.0),
    ])
sort_pts!(p)
@test get_distance(p,Circle(0.0,0.0,1.0)) == 0.0
@test round(get_distance(p,Circle(4.0,0.0,1.0)),digits=4) == 1.0

@test norm(nearest_in_bounds_space(Circle(4.0,0.0,1.0),p) - VecE2(-3.0,0.0)) < 0.000001
