module JuliaJuliaSet

using LinearAlgebra
using Plots, Colors

# ----- Functions for generating Julia and Mandelbrot sets -----

function julia(f::Function, c::Complex{Float64}, z::Complex{Float64}, threshold::Float64)
    for i = 1:255
        z = f(z, c)
        if abs(z) >= threshold
            return i
        end
    end
    return -255
end

function mandelbrot(f::Function, c::Complex{Float64}, threshold::Float64)
    z = 0
    for i = 1:255
        z = f(z, c)
        if abs(z) >= threshold
            return i
        end
    end
    return -0
end

function setdims(max_coord, min_coord, resolution)
    dim = (max_coord - min_coord) * resolution
    @assert dim != 0
    return dim
end

struct SetParams
    min_coord::Complex{Float64}
    max_coord::Complex{Float64}
    resolution::Int64
    width::Int64
    height::Int64
    threshold::Float64
    nr_frames::Int64
    function set_p(min_coord, max_coord, resolution, width, height, threshold, nr_frames) 
        return new(min_coord, 
            max_coord, 
            resolution, 
            setdims(max_coord.re, min_coord.re, resolution), 
            setdims(max_coord.im, min_coord.im, resolution), 
            threshold,
            nr_frames) 
    end
end 

function genplane(set_p::SetParams)
    real = range(set_p.min_coord.re, set_p.max_coord.re,length=set_p.width)
    imag = range(set_p.min_coord.im, set_p.max_coord.im,length=set_p.height)
    complexplane = zeros(Complex{Float64},(set_p.width,set_p.height))
    for (i,x) ∈ enumerate(real)
        complexplane[:,i] .+= x
    end
    for (i,y) ∈ enumerate(imag)
        complexplane[i,:] .+= (y * 1im)
    end
    return reverse!(complexplane, dims=1)
end

function genmandelbrotset(set_p::SetParams, f::Function, plane::Matrix{ComplexF64})
    return mandelbrot.(f, plane, set_p.threshold)
end

function genjuliaset(set_p::SetParams, c, f::Function, plane::Matrix{ComplexF64})
    return julia.(f, c, plane, set_p.threshold)
end

function genjuliaprogression(set_p::SetParams, P::Path, f::Function, plane::Matrix{ComplexF64})
    c_vec = pointsoncruve(P,sep_p.nr_frames)
    return [genjuliaset(set_p, c, f, plane) for c ∈ c_vec]
end

# ----- Functionality for defining paths through the complex plane -----

struct Path
    parameterization::Function
    start::ComplexF64
    ending::ComplexF64
end

Circle_Path(r::Real)::Path = Path(t-> (r*ℯ^(t*1im)),0,2π)

pointsoncurve(P::Path, n::Integer) = [P.parameterization(t) for t in range(P.start, P.ending,length=n)]

# ----- Functions for animating series of sets -----

function animateprogression(sets::Vector{Matrix{Int64}}, set_p::SetParams, file_name::String ="~/GIFs/julia_set.gif")
    anim = @animate for set ∈ sets
        heatmap(set, size=(set_p.width,set_p.height), color=:terrain, leg=false)
    end
    gif(anim, file_name, fps=30)
end