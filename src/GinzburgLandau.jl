module GinzburgLandau

# package code goes here

dct_nodes(N, L) = (0:N-1)*L/N + L/(2N)


AdamsMoulton = [1          0          0         0         0          0;
                1/2        1/2        0         0         0          0;
                5/12       8/12      -1/12      0         0          0;
                9/24       19/24     -5/24      1/24      0          0;
                251/720    646/720   -264/720   106/720  -19/720     0;
                475/1440   1427/1440 -798/1440  482/1440 -173/1440   27/1440]

AdamsBashforth=[1          0          0          0        0          0;
                3/2        -1/2       0          0        0          0;
                23/12      -16/12     5/12       0        0          0;
                55/24      -59/24     37/24     -9/24     0          0;  
                1901/720   -2774/720  2616/720  -1274/720  251/720   0;
                4277/1440  -7923/1440 9982/1440 -7298/1440 2877/1440 -475/1440]

type GLECoeffs
    σ::Float64
    γ::Float64
    c₁::Float64
    c₃::Float64
    k::Vector{Complex{Float64}}
end

function GLECoeffs(σ, γc, c₁, c₃)
    k = zeros(Complex{Float64}, 3)
    γ = γc/c₁
    k[1] = σ
    k[2] = γ*(1 + im*c₁)
    k[3] = -σ*(1 - im*c₃)
    GLECoeffs(σ, γ, c₁, c₃, k)
end

import Base.getindex

getindex(coeffs::GLECoeffs, i) = coeffs.k[i]


type GLESolver1
    L::Float64
    N::Int
    Δt::Float64
    δ::Vector{Complex{Float64}}
    g::Vector{Complex{Float64}}
end

function GLESolver1(L, N, Δt, k::GLECoeffs)
    
    δ = zeros(Complex{Float64}, N)
    
    for i = 1:N
        w = π/L*(i-1)
        δ[i] = 1.0 - Δt*(k[1] - k[2]*w^2)
    end

    g = zeros(Complex{Float64}, N)
    GLESolver1(L, N, Δt, δ, g)
end

nmodes(slv::GLESolver1) = slv.N
nodes(slv::GLESolver1) = dct_nodes(slv.N, slv.L)
tstep(slv::GLESolver1, i) = (i-1)*slv.Δt


function gle_step!(slv::GLESolver1, k::GLECoeffs, a₀::AbstractVector{Complex{Float64}},
                   a₁::AbstractVector{Complex{Float64}})
    # Calcularo f
    g = slv.g
    dt = slv.Δt
    N = slv.N
    δ = slv.δ
    
    copy!(a₁, a₀)
    idct!(a₁)
    
    for i = 1:N
        g[i] = k[3] * abs(a₁[i])^2 * a₁[i]
    end

    dct!(g)
    for i = 1:N
        a₁[i] = ( a₀[i] + dt*g[i] ) / δ[i]
    end


    idct!(a₁)
    for i = 1:N
        g[i] = k[3] * abs(a₁[i])^2 * a₁[i]
    end
    dct!(g)
    
    for i = 1:N
        a₁[i] = ( a₀[i] + dt*g[i] ) / δ[i]
    end

    return nothing
    
end



type GLEAdamsBashforth
    L::Float64
    N::Int
    s::Int
    Δt::Float64
    f::Vector{Vector{Complex{Float64}}
    a₀::Vector{Complex{Float64}}
    idx::Vector{Integer}
    step::Int
    function GLEAdamsBashforth(L, N, s, Δt)
        f = Vector{Vector{Complex{Float64}}}(s)
        for i = 1:s
            f[i] = zeros(Complex{Float64}, N)
        end
        a₀ = zeros(Complex{Float64}, N)
        new(L, N, s, Δt, f, a₀, zeros(Int, s), 0)
    end

    
end
nmodes(slv::GLEAdamsBashforth) = slv.N
nodes(slv::GLEAdamsBashforth) = dct_nodes(slv.N, slv.L)
tstep(slv::GLEAdamsBashforth, i) = (i-1)*slv.Δt
nstep(slv::GLEAdamsBashforth) = slv.s


function init_gle!(slv::GLEAdamsBashforth, a::AbstractVector)
    s = nstep(slv)
    slv.step = 0
    fill!(slv.idx, 0)
    for i = 1:s
        fill!(slv.f[i], 0)
    end
    fill!(slv.a₀, a)
    
end


function simple_step!(slv::GLEAdamsBashforth, a::AbstractVector{Complex{Float64}})
    

                    
end # module
