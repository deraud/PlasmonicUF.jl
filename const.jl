module Constants

using StaticArrays

export Nx, Ny, Nt, NtScale, CellType

const Nx::Int = 128
const Ny::Int = 128
const Nt::Int = 1000
const NtScale::Int = 1

module CellType
    const SOLID = 2^0         
    const FLUID = 2^1         
    const GAS = 2^2
end

end