module Raygun

using Colors
using LinearAlgebra
using ProgressBars
using StaticArrays

Vec3 = SVector{3, Float64}
Point3 = SVector{3, Float64}

norm2(vec::Vec3) = dot(vec, vec)
unit_vector(vec::Vec3) = vec/sqrt(norm2(vec))

struct Ray
    origin::Point3
    direction::Vec3
end

position(ray::Ray, t::Float64) = origin + t*direction

function hit_sphere(center::Point3, radius::Float64, ray::Ray)
    oc = ray.origin - center
    a = dot(ray.direction, ray.direction)
    b = 2*dot(oc, ray.direction)
    c = dot(oc,oc) - radius*radius
    Δ = b*b - 4*a*c
    return (Δ > 0)
end

function color(ray::Ray)
    if hit_sphere(Point3(0,0,-1), 0.5, ray)
        return RGB(1,0,0)
    end
    unit_direction = unit_vector(ray.direction)
    t = 0.5*(unit_direction[2] + 1.0)
    return weighted_color_mean(t, RGB(0.5,0.7,1.0), RGB(1.0, 1.0, 1.0))
end

function render()
    aspect_ratio = 16/9
    width = 400
    height = convert(Int, width/aspect_ratio)

    viewport_height = 2.0
    viewport_width = viewport_height*aspect_ratio
    focal_length = 1.0

    origin = Point3(0,0,0)
    horizontal = Vec3(viewport_width,0,0)
    vertical = Vec3(0,viewport_height,0)
    lower_left = origin - horizontal/2 - vertical/2 - Vec3(0,0,focal_length)

    img = Array{RGB, 2}(undef, height, width)
    for j in ProgressBar(1:height)
        for i in 1:width
            u,v = i/width, (height-j)/height
            ray = Ray(origin, lower_left + u*horizontal + v*vertical - origin)
            img[j,i] = color(ray)
        end
    end
    return img
end
end # module
