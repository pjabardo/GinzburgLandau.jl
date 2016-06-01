module GinzburgLandau

# package code goes here

dct_nodes(N, L) = (0:N-1)*L/N + L/(2N)


AdamsMoulton = [1     0      0     0    0    0;
                1     1      0     0    0    0;
                5     8     -1     0    0    0;
                9     19    -5     1    0    0;
                251   646  -264  106  -19    0;
                475   1427 -798  482 -173   27]

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
    a1::Vector{Complex{Float64}}
end

function GLESolver1(L, N, Δt, k::GLECoeffs)
    
    δ = zeros(Complex{Float64}, N)
    
    for i = 1:N
        w = π/L*(i-1)
        δ[i] = 1.0 - Δt*(k[1] - k[2]*w^2)
    end

    g = zeros(Complex{Float64}, N)
    a1 = zeros(Complex{Float64}, N)
    GLESolver1(L, N, Δt, δ, g, a1)
end

nmodes(slv::GLESolver1) = slv.N
nodes(slv::GLESolver1) = dct_nodes(slv.N, slv.L)


function gle_init!(slv::GLESolver1, a::Vector{Complex{Float64}})

    copy!(slv.a1, a)
    return nothing
end

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



    

                    
end # module
