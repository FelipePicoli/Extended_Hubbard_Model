
############################################################################################################################
####################################################--- Functions ---#######################################################
############################################################################################################################

####################################################--- MPO construction ---################################################

function MPO_construction(N,sites,U,V,ts)

    ost=OpSum()   # chain object
    osU=OpSum()   # resultion object
    osV=OpSum()   # intersite repulsion object

    for i=1:1:(N-1)  

        ost+=ts[i],"Cdagup",i,"Cup",i+1
        ost+=ts[i],"Cdagup",i+1,"Cup",i
        ost+=ts[i],"Cdagdn",i,"Cdn",i+1
        ost+=ts[i],"Cdagdn",i+1,"Cdn",i
            
        osV+=V[i],"Ntot",i,"Ntot",i+1
        
    end

    for i=1:1:N 

        osU+=U[i],"Nupdn",i 

    end

    Hchain=MPO(ost,sites)
    HU=MPO(osU,sites)
    HV=MPO(osV,sites)

    return Hchain+HU+HV

end

####################################################--- DMRG optimization ---################################################

function DMRG_optm(N,Nup,Ndn,H,sites,Nsweep,maxD,cutoff=1e-8,psi0_D=10)

    initial_state=get_state(N,Nup,Ndn)

    psi0=randomMPS(sites,initial_state; linkdims=psi0_D)

    energy,psi=dmrg(H,psi0;nsweeps=Nsweep,maxdim=maxD,cutoff=cutoff)

    return energy,psi
end

####################################################--- Initial state construction ---################################################

function get_state(L,Nup,Ndn)

    state=fill("Emp",L)

    p=Nup+Ndn

    for i=1:1:L 

        if i%2==0
            state[i]="Up"
        else
            state[i]="Dn"
        end

    end

    return state
end

####################################################--- One particle Reduced Density Matrix ---################################################

#=
function One_particle_RDM(Psi)

    L=length(Psi)

    rho_1=zeros(ComplexF64,2*L,2*L)

    rho_upup=correlation_matrix(Psi,"Cdagup","Cup")
    rho_dndn=correlation_matrix(Psi,"Cdagdn","Cdn")
    rho_updn=correlation_matrix(Psi,"Cdagup","Cdn")

    for i=1:1:L 
        for j=1:1:L 

            rho_1[i,j]=rho_upup[i,j]
            rho_1[i,j+L]=rho_updn[i,j]
            rho_1[i+L,j+L]=rho_dndn[i,j]
            rho_1[i+L,j]=rho_updn[i,j]'

        end
    end

    return rho_1./L
end
=#

function One_particle_RDM(Psi,i)

    orthogonalize!(Psi,i)

    Ai = Psi[i]
    Ai_dag = dag(prime(Ai,"Site"))

    rho=Ai*Ai_dag

    return rho, Array(rho.tensor)
end

####################################################--- Two particle Reduced Density Matrix ---################################################

function Two_sites_RDM(Psi,i,j)

    if i>j
        i,j = j,i
    end

    orthogonalize!(Psi,i)  # Move o centro de ortogonalidade para i 
    
    Psi_bra=dag(Psi) # Calcula o dag
    prime!(Psi_bra,"Link")  # prima o dag

    if i!=1
            
        li_1=linkind(Psi,i-1) # pega o nome da perda a esquerda de i

        rho=prime(Psi[i],li_1)*prime(Psi_bra[i],"Site")

    else

        rho=Psi[i]*prime(Psi_bra[i],"Site")

    end

    # as contrações entre os indices 

    for k=(i+1):1:(j-1)

        rho*=Psi[k]
        rho*=Psi_bra[k]

    end

    if j!=length(Psi)

        lj=linkind(Psi,j)

        rho*=prime(Psi[j],lj)
        rho*=prime(Psi_bra[j],"Site")

    else

        rho*=Psi[j]
        rho*=prime(Psi_bra[j],"Site")

    end

    return rho
end

function Two_particle_cor_matrix_sites(Psi,site_i,site_j)

    basis=[
    (1,4), # (0,duplo)
    (2,2), # (up,up)
    (3,3), # (dn,dn)
    (2,3), # (up,dn)
    (3,2), # (dn,up)
    (4,1)  # (duplo,0)
    ]

    rho_ij=Two_sites_RDM(Psi,site_i,site_j)

    Matrix=zeros(ComplexF64,6,6)

    si_out,si_in,sj_out,sj_in=inds(rho_ij)

    for i=1:1:6
        for j=1:1:6

            sio,sjo=basis[i]
            sii,sji=basis[j]

            Matrix[i,j]+=rho_ij[si_out=>sio,sj_out=>sjo,si_in=>sii,sj_in=>sji]

        end
    end

    return Matrix
end
 
function Two_particle_RDM(Psi)

    L=length(Psi)

    matrix=zeros(ComplexF64,6*L,6*L)

    for site_i=1:1:L
        for site_j=1:1:L

            if site_i<site_j

                println(site_i,"   ",site_j)

                m=Two_particle_cor_matrix_sites(Psi,site_i,site_j)

                i1=6*(site_i-1)+1
                i2=6*site_i
                j1=6*(site_j-1)+1
                j2=6*site_j

                matrix[i1:i2,j1:j2].+=m
                matrix[j1:j2,i1:i2].+=m'

            end
        end
    end

    return matrix
end

#=
function Two_particle_RDM(Psi,sites)
    
    L=length(sites)
    siz=(size(Psi[1])[1])^2*L
    size_site=(size(Psi[1])[1])^2
    
    rho2=zeros(ComplexF64,siz,siz)
    
    idxi=0
    
    for i=1:1:L
        idxj=0

        for j=(i+1):1:L
            
            _,rho_matrix=RDM_contraction_two_sites(Psi,i,j)
            
            println(size(rho2[(i+idxi):(i+idxi+size_site-1),(i+idxj):(i+idxj+size_site-1)]))
            println(size(rho_matrix))

            rho2[(1+idxi):(1+idxi+size_site-1),(1+idxj):(1+idxj+size_site-1)].+=rho_matrix
        
            idxj+=size_site 

        end 
        
        idxi+=size_site

    end 
    
    return (rho2.+rho2').*2(L*(L-1)) 
end

function Two_particle_RDM(Psi,sites)

    L=length(Psi)
    spins=["up","dn"]

    rho_2=zeros(ComplexF64,L,2,L,2,L,2,L,2)

    for i=1:1:L, j=1:1:L, k=1:1:L, l=1:1:L
        for (is,si) in enumerate(spins), (js,sj) in enumerate(spins), (ks,sk) in enumerate(spins), (ls,sl) in enumerate(spins)

            if yes_or_not(si,sj,sk,sl,i,j,k,l)
                rho_2[i,is,j,js,k,ks,l,ls]=Four_op_correlation(Psi,i,j,k,l,si,sj,sk,sl,sites)
            end

        end
    end
 
    return 2/(L*(L-1)).*reshape(rho_2,(2*L)^2,(2*L)^2)
end

function Two_particle_RDM_block(Psi,sites)

    L=length(Psi)
    spins=["up","dn"]

    rho_2=zeros(ComplexF64,L,2,L,2,L,2,L,2)

    rho_upupupup=zeros(ComplexF64,L,L,L,L)
    rho_dndndndn=zeros(ComplexF64,L,L,L,L)
    rho_dnupdnup=zeros(ComplexF64,L,L,L,L)

    for i=1:1:L, j=1:1:L, k=1:1:L, l=1:1:L
       
        if yes_or_not("up","up","up","up",i,j,k,l)

            rho_upupupup[i,j,k,l]=Four_op_correlation(Psi,i,j,k,l,"up","up","up","up",sites)

        end

        if yes_or_not("dn","dn","dn","dn",i,j,k,l)

            rho_dndndndn[i,j,k,l]=Four_op_correlation(Psi,i,j,k,l,"dn","dn","dn","dn",sites)

        end
        
        if yes_or_not("dn","up","dn","up",i,j,k,l)

            rho_dnupdnup[i,j,k,l]=Four_op_correlation(Psi,i,j,k,l,"dn","up","dn","up",sites)

        end

    end

    rho_upupupup=reshape(rho_upupupup,L^2,L^2)
    rho_dndndndn=reshape(rho_dndndndn,L^2,L^2)
    rho_dnupdnup=reshape(rho_dnupdnup,L^2,L^2)

    mat=2/(L*(L-1)).*cat(rho_upupupup,rho_dndndndn,rho_dnupdnup;dims=(1,2))

    return (mat.+mat')./2
end

function Four_op_correlation(Psi,i,j,k,l,si,sj,sk,sl,sites)

    O = OpSum()
    O += 1.0,"Cdag"*si, i, "Cdag"*sj, j,"C"*sk,  k,"C"*sl,  l

    H = MPO(O, sites)

    return inner(Psi,Apply(H,Psi))
end

function map_spin(s)
    if s=="up"
        return +1.0
    else
        return -1.0
    end
end

function spin_conserve(si,sj,sk,sl)

    s=0.0

    s=s+map_spin(si)+map_spin(sj)-map_spin(sk)-map_spin(sl)

    if s==0.0
        return true
    else
        return false
    end   
end

function yes_or_not(si,sj,sk,sl,i,j,k,l)

    yes_or_not_spin=spin_conserve(si,sj,sk,sl)

    yes_or_not_paulli=true

    if i==j && si==sj

        yes_or_not_paulli=false

    elseif k==l && sk==sl

        yes_or_not_paulli=false

    end

    return yes_or_not_paulli*yes_or_not_spin
end
=#
####################################################--- von Neumann Entropy ---################################################

function S_vNeumann(rho,N)

    lambs=max.(0.0,eigvals(rho))
    
    return -1.0*sum(lambs.*log.(lambs))-log(N)
end

####################################################--- Main Function ---################################################

function Main_func_cluster(U,V,N,Nsweep,maxD,cutoff)

    Npart=floor(Int,N/2)
    Nup=Npart+N%2
    Ndn=N-Nup
    
    sites=siteinds("Electron", N; conserve_qns = true)

    H=MPO_construction(N,sites,U.*ones(N),V.*ones(N),-1.0.*ones(N))

    E0,Psi=DMRG_optm(N,Nup,Ndn,H,sites,Nsweep,maxD,cutoff)

    _,rho1=One_particle_RDM(Psi,Int(N/2))

    C1_l=sum(abs.(rho1))-sum(diag(abs.(rho1)))

    Entro_1=S_vNeumann(rho1,N)

    _,rho2=Two_particle_RDM(Psi,Int(N/2),Int(N/2+1))

    C2_llp1=sum(abs.(rho2))-sum(diag(abs.(rho2)))

    Entro_2llp1=S_vNeumann(rho2,N)

    _,rho2=Two_particle_RDM(Psi,Int(N/2),Int(N/2+2))

    C2_llp2=sum(abs.(rho2))-sum(diag(abs.(rho2)))

    Entro_2llp2=S_vNeumann(rho2,N)

    return E0, Entro_1, C1_l, Entro_2llp1, C2_llp1, Entro_2llp2, C2_llp2
end
