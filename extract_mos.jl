using LinearAlgebra

#EXTRACT DATA

function extract_N(file_name)
    open(file_name) do f
    g=[parse(Int,x) for x in split(readline(f))]
end
end
g=extract_N("data_dir/nums_orbitals.dat")
Nb=g[1]
Nd=g[2]
Ns=g[3]

function extract_matrix(file_name,Nb)
    @info "Extract matrix :", file_name
    open(file_name) do f
        lines = readlines(f)
        g= zeros(Nb,Nb)
        for i=1:Nb
            g[i,:]=[parse(Float64,x) for x in split(lines[i]) ]
        end
        g
    end
end

# Extract data
Pd = Symmetric(extract_matrix("data_dir/Pd.dat",Nb))
Ps = Symmetric(extract_matrix("data_dir/Ps.dat",Nb))
S = Symmetric(extract_matrix("data_dir/overlap_matrix.dat",Nb))
S12=real(sqrt(S))
Sm12=inv(S12);

# Test validity of given matrices in the orthonormal basis.
@assert round(tr(Pd)) == Nd
@assert round(tr(Ps)) == Ns
@assert round.(S*(Sm12*Sm12))==I(Nb)


#Extract occupied Mos
@info "Orthonormalizing"
Nd = Int(round(tr(Pd)))
Ns = Int(round(tr(Ps)))
Cd  = eigen(-Pd).vectors[:,1:Nd]
Cs  = eigen(-Ps).vectors[:,1:Ns]

# Gram-Schmidt
"""
    Orthonormalize vector b with the n first columns of A
"""
function gram_schmidt(A::AbstractArray{T}, b::AbstractVector{T},n::Integer) where T
    k = size(A,2)  #number of orthogonal vectors
    dot_prods = zeros(T,k)  
    
    for (i,a_i) in enumerate(eachcol(A[:,1:n]))
        dot_prods[i] = dot(a_i,b)
        axpy!(-dot_prods[i],a_i,b)
    end
    nrm = norm(b)
    b *= one(T)/nrm
    b
end

"""
    Uses the preceding routine to obtain a set of orthonormalized orbitals
"""
function generate_new_orbs(Cd::AbstractArray{T},Cs::AbstractArray{T},
                           Nb::Integer,Nd::Integer,Ns::Integer) where T
    #Initialize new orbs, fill virtuals randomly
    N_occupied = Nd+Ns
    new_orbitals = zeros(T,Nb,Nb)
    new_orbitals[:,1:Nd+Ns] = hcat(Cd,Cs)
    virtuals = rand(T,Nb,Nb-(Nd+Ns))

    for (i,v) in enumerate(eachcol(virtuals))
        ortho_v = gram_schmidt(new_orbitals,v,N_occupied+(i-1))
        new_orbitals[:,N_occupied+i] = ortho_v
    end

    new_orbitals
end

# orthonormalization and "de"-orthogonalize basis
new_orbitals = Sm12*generate_new_orbs(Cd,Cs,Nb,Nd,Ns)

# test output
Cd_new = S12*new_orbitals[:,1:Nd]
Cs_new = S12*new_orbitals[:,Nd+1:Nd+Ns]
@assert round(tr(Cd_new*Cd_new'))==Nd
@assert round(tr(Cs_new*Cs_new'))==Ns

#write orbitals
open("data_dir/new_orbitals.dat", "w") do out
    @info "Writing new orbitals",out
    for i in 1:Nb
        line = " "
        for x in new_orbitals[i,:]
            line *=string(x)*" "
        end
       println(out, line)
    end
end
