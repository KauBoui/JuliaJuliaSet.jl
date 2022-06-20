module JuliaJuliaSet

using LinearAlgebra
using Distributed
using Plots

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
            (max_coord.re - min_coord.re) * resolution ; @assert width != 0,
            (max_coord.im - min_coord.im) * resolution ; @assert height != 0, 
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

function genmandelbrotset(set_p::SetParams, f::Function, plane::Matrix{Complex64})
    return mandelbrot.(f, plane, set_p.threshold)
end

function genjuliaset(set_p::SetParams, c, f::Function, plane::Matrix{Complex64})
    return julia.(f, c, plane, set_p)
end


