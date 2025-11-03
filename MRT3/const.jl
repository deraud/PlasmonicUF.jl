module Constants_MRT3

using StaticArrays

export Q, ex, ey, w, M, Minv, ρ0, G, τ, s

const Q::Int = 9
const ex = [0, 1, 0, -1, 0, 1, -1, -1, 1]
const ey = [0, 0, 1, 0, -1, 1, 1, -1, -1]
const w = [4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36]

const M = [
        1  1  1  1  1  1  1  1  1;
        -4 -1 -1 -1 -1  2  2  2  2;
        4 -2 -2 -2 -2  1  1  1  1;
        0  1  0 -1  0  1 -1 -1  1;
        0 -2  0  2  0  1 -1 -1  1;
        0  0  1  0 -1  1  1 -1 -1;
        0  0 -2  0  2  1  1 -1 -1;
        0  1 -1  1 -1  0  0  0  0;
        0  0  0  0  0  1 -1  1 -1
    ]

const Minv = inv(M)

const ρ0::Float64 = 1.0
const G::Float64 = -5.0

const τ = 1.0
const s = [0.0, 1.4, 1.4, 0.0, 1.2, 0.0, 1.2, 1/τ, 1/τ]


end