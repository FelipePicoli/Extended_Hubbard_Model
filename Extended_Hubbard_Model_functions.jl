
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

    psi0=random_mps(sites,initial_state; linkdims=psi0_D)

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

####################################################--- Two particle Reduced Density Matrix ---################################################

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
 
    return reshape(rho_2,(2*L)^2,(2*L)^2)
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

####################################################--- von Neumann Entropy ---################################################

function S_vNeumann(Psi,N)

    rho_1=One_particle_RDM(Psi)

    lambs=max.(0.0,eigvals(rho_1))
    
    return -1.0*sum(lambs.*log.(lambs))-log(N)
end