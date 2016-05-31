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
    α::Vector{Complex{Float64}}
    w::Vector{Float64}
    f::Vector{Complex{Float64}}
    a1::Vector{Complex{Float64}}
end

function GLESolver1(L, N, Δt, k::GLECoeffs)
    
    α = zeros(Complex{Float64}, N)
    
    for i = 1:N
        w = π/L*(i-1)
        α[i] = 1.0 - Δt(k[1] - k[2]*w^2)
    end

    f = zeros(Complex{Float64}, N)
    a1 = zeros(Complex{Float64}, N)
    GLESover1(L, N, Δt, α, w, f, a1)
end

function gle_step!(slv::GLESolver1, k::GLECoeffs, a::Vector{Complex{Float64}})
    # Calcularo f
    f = slv.f
    a1 = slv.a1
    dt = slv.Δt
    N = slv.N
    copy!(a1, a)
    idct!(a1)
    
    for i = 1:N
        f[i] = k[3] * abs(a1[i])^2 * a1[i]
    end

    dct!(f)
    α = slv.α
    for i = 1:N
        a1[i] = ( a[i] + dt*f[i] ) / α[i]
    end


    idct!(a1)
    for i = 1:N
        f[i] = k[3] * abs(a1[i])^2 * a1[i]
    end
    dct!(f)
    
    for i = 1:N
        a[i] = ( a[i] + dt*f[i] ) / α[i]
    end

    return a
    
end



    

                    
end # module
