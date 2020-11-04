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

Pd = Symmetric(extract_matrix("data_dir/Pd.dat",Nb))
Ps = Symmetric(extract_matrix("data_dir/Ps.dat",Nb))



S = Symmetric(extract_matrix("data_dir/overlap_matrix.dat",Nb))
S12=real(sqrt(S))
Sm12=inv(S12);

#Extract occupied Mos
@info "Orthonormalizing"
Nd = Int(round(tr(Pd)))
Ns = Int(round(tr(Ps)))

Cd  = eigen(-Pd).vectors[:,1:Nd]
Cs  = eigen(-Ps).vectors[:,1:Ns]

orbitals = diagm(ones(Nb))
orbitals[:,1:Nd+Ns] = hcat(Cd,Cs)



# Generate orthonormal virtual orbitals from AO basis.

function gm_ortho_bloc(A,B)
    # A has to be already orthonormal orbitals
    # B is a guess for virtual to be orthonormalized with A
    B = Matrix(qr(B).Q) #Orthogonalize B if not orthogonal with QR ortho
    B = B - A * A'B #Gram-Schmidt per Bloc
    @assert norm(A'B) < 1e-10
    A,B
end

@views function generate_virtuals(orbitals, n_occupied)
    occ_orb = orbitals[:,1:n_occupied]
    virt_orb = orbitals[:,n_occupied+1:end]
    occ_orb,virt_orb = gm_ortho_bloc(occ_orb,virt_orb)
    hcat(occ_orb,virt_orb)
end

new_orbitals = Sm12*generate_virtuals(orbitals,Nd+Ns)

# # Test properties
# Cd_new = new_orbitals[:,1:Nd]
# Cs_new = new_orbitals[:,Nd+1:Nd+Ns]
# Pd_new = Cd_new*Cd_new'
# Ps_new = Cs_new*Cs_new'


# S12=real(sqrt(S))
# Sm12=inv(S12);





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
