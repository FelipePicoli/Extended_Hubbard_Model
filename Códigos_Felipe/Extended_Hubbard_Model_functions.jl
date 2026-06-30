
############################################################################################################################
####################################################--- Functions ---#######################################################
############################################################################################################################

####################################################--- MPO construction ---################################################

function MPO_construction(N,sites,U,V)

    os=OpSum()   # chain object

    for i=1:1:(N-1)  

        os+=-1.0,"Cdagup",i,"Cup",i+1
        os+=-1.0,"Cdagup",i+1,"Cup",i
        os+=-1.0,"Cdagdn",i,"Cdn",i+1
        os+=-1.0,"Cdagdn",i+1,"Cdn",i
            
        os+=V,"Ntot",i,"Ntot",i+1
        
    end
#=
    for i in [1, N]

    os += V, "Ntot", i

    os += -2*V, "Nupdn", i

    end

    os+=-V,"Ntot",1,"Ntot",1
    os+=2*V,"Ntot",1
    os+=-V,"Id",1

    os+=-V,"Ntot",N,"Ntot",N
    os+=2*V,"Ntot",N
    os+=-V,"Id",N
  

    os+=1e-2,"Ntot",1
    os+=-1e-2,"Ntot",N
=#
    for i=1:1:N 

        os+=U,"Nupdn",i 

    end

    return MPO(os,sites)

end

####################################################--- Getting Initial State ---################################################

function get_region_phase_diagram(U,V)

    alpha=0.2

    get_metallic_lower_boundary(U) = U <= 0 ? -exp(alpha * U) : -alpha * U - 1.0
    get_metallic_upper_boundary_h(U) = -0.1 * U
    get_metallic_upper_boundary_negative(U) = -2.5 * alpha * U

    lower = get_metallic_lower_boundary(U)
    upper_h = get_metallic_upper_boundary_h(U)
    upper_neg = get_metallic_upper_boundary_negative(U)

    if V <= lower
        return "PS"
    elseif (U < 0 && V < 0 && V > lower) ||
           (U > 0 && V > 0 && V < upper_h) ||
           (U > 0 && V < 0 && V > lower && V < upper_neg)
        return "METALLIC"
    elseif (U < 0 && V > 0) ||
           (U > 0 && V > 0 && V >= 2 * U)
        return "CDW"
    elseif (U > 0 && V < 0 && V > lower) ||
           (U > 0 && V > 0 && V >= upper_h && V < 2 * U)
        return "SDW"
    end
    return "PS"
end

function state_ehm_diagram(L, Nup, Ndn, U, V)
    state = fill("Emp", L)
    region = get_region_phase_diagram(U, V)

    if region == "METALLIC"
        println("METALLIC")
        state = product_metallic_state(L, Nup, Ndn)
    elseif region == "CDW"
        println("CDW")
        state = product_cdw_state(L, Nup, Ndn)
    elseif region == "SDW"
        println("SDW")
        state = product_sdw_state(L, Nup, Ndn)
    else
        println("PS")
        state = product_ps_state(L, Nup, Ndn)
    end
    return state
end

function product_ps_state(N, Nup, Ndn)
   
    Nocc=Int(N/2)

    state = fill("Emp", N)

    center = div(N,2)
           
    start = center - div(Nocc,2) + 1
    stop  = start + Nocc - 1

    for i in start:stop
        state[i] = "UpDn"
    end

    return state
end

function product_sdw_state(L, Nup, Ndn)

    state = fill("Up", L)

    for i=1:1:L 

        if isodd(i)
            state[i]="Dn"

        end
    end
    return state
end


function product_cdw_state(L, Nup, Ndn)

    state = fill("Emp", L)
   
    for i=1:1:L 

        if isodd(i)
            state[i]="UpDn"
        end
    end
    return state
end


function product_metallic_state(L, Nup, Ndn)
    state = fill("Up", L)
    for i=1:1:L 

        if isodd(i)
            state[i]="Dn"

        end
    end

    return state
end

####################################################--- DMRG optimization ---################################################

function DMRG_optm(U,V,N,Nup,Ndn,H,sites,Nsweep;maxD=100,cutoff=1e-8)

  # initial_state=get_state(N,Nup,Ndn) #state_ehm_diagram(N, Nup, Ndn, 0,0) 

    #initial_state = [i <= N÷2 ? "UpDn" : "Emp" for i in 1:N]
    if V>=0.0
        initial_state=product_cdw_state(N, Nup, Ndn)
    else
        initial_state=product_ps_state(N, Nup, Ndn)
    end
    #initial_state=state_ehm_diagram(N,Nup,Ndn,U,V)

   # initial_state[Int(N/2)+1:L].="UpDn"

    psi0=productMPS(sites,initial_state)
  #  psi0=randomMPS(sites,initial_state)#; linkdims=psi0_D)

    energy,psi=dmrg(H,psi0,nsweeps=Nsweep,maxdim=maxD,cutoff=cutoff)

    psi_=copy(apply(H,psi))
    E2=inner(psi_,psi_)

    println("Ground state variance: ",E2-energy^2)

    return energy,psi
end

function centered_density_state(N, Nocc)

    state = fill("Emp", N)

    center = div(N,2)
           
    start = center - div(Nocc,2) + 1
    stop  = start + Nocc - 1

    for i in start:stop
        state[i] = "UpDn"
    end

    return state
end

function DMRG_optm_convergence(U,V,N,Nup,Ndn,H,sites,Nsweep,maxD;cutoff=1e-8,psi0_D=10)

    results=zeros(Float64,Nsweep,8)
 
    if V>=0.0
        initial_state=product_cdw_state(N, Nup, Ndn)
    else
        initial_state=product_ps_state(N, Nup, Ndn)
    end

    println(initial_state)
    psi=productMPS(sites,initial_state)

    for is=1:1:Nsweep

        if is>length(maxD)
            mD=maxD[end]
        else
            mD=maxD[is]
        end

        energy,psi=dmrg(H,psi;nsweeps=1,maxdim=[mD],cutoff=cutoff)

        # Variância 

        psi_=copy(apply(H,psi))
        var=Float64(inner(psi_,psi_))-Float64(energy)^2

        # Ocupação do sítio central

        nc=expect(psi,"Ntot")

        # entropia do sítio central

        _,rho_c=One_site_RDM(psi,Int(N/2))

        entropy_c=S_vNeumann(rho_c)/2

        _,rho_c=One_site_RDM(psi,Int(N/2)+1)

        entropy_c+=S_vNeumann(rho_c)/2

        # entropia do sítio central e dos dois laterais

        _,rho_l=Two_sites_RDM(psi,sites,Int(N/2),Int(N/2)-1)
        _,rho_r=Two_sites_RDM(psi,sites,Int(N/2),Int(N/2)+1)

        entropy_l=S_vNeumann(rho_l)/2
        entropy_r=S_vNeumann(rho_r)/2

        _,rho_l=Two_sites_RDM(psi,sites,Int(N/2)+1,Int(N/2)+2)
        _,rho_r=Two_sites_RDM(psi,sites,Int(N/2)+1,Int(N/2))

        entropy_l+=S_vNeumann(rho_l)/2
        entropy_r+=S_vNeumann(rho_r)/2

        Rug,Entro_2=Rugosity_2sites2(psi,sites,Int(N/2),1;vN=true)

        maximumD=maximum(i -> linkdim(psi, i), 1:length(psi)-1)

        println("Sweep: ",is,"/",Nsweep, " Energy: ", energy, " Variance: ", var, " Max. Bound Dim.: ",maximumD, " Central site occu: ",nc[Int(N/2)], " One site entro.: ",entropy_c, " Two site entro. left & right: ", entropy_l," ",entropy_r," Two site Rug. left & right: ",Rug[1]," ",Rug[2])
        println(nc)
        println("##########################################################################################################################################################################")

        results[is,1]+=is 
        results[is,2]+=energy
        results[is,3]+=var
        results[is,4]+=maximumD
        results[is,5]+=nc[Int(N/2)]
        results[is,6]+=entropy_c
        results[is,7]+=entropy_l
        results[is,8]+=entropy_r

    end

    return energy,psi

end


####################################################--- Initial state construction ---################################################

function get_state(L, Nup, Ndn)
    state = fill("Emp", L)

    # coloca ups
    for i in 1:Nup
        state[i] = "Up"
    end

    # coloca downs
    for i in (Nup+1):(Nup+Ndn)
        state[i] = "Dn"
    end

    return state
end

#=
function random_metallic_state(L, Nup, Ndn)
    state = fill("Emp", L)
    p = Nup + Ndn
    for i in 1:L
        j = L - i
        if(p > j)
            println(p," ",j)
            state[j] = "UpDn"
            p -= 2
        elseif (p > 0)
            println(p," ",j)
            state[j] = j % 2 == 1 ? "Up" : "Dn"
            p -= 1
        end
    end
    return state
end

function random_cdw_state(L, Nup, Ndn)
    state = fill("Emp", L)
    Nup_extra = Int.(Nup % Ndn)
    for i in 1:2:(L - Nup_extra)
        state[i] = "UpDn"
        Nup-=1
    end
    if Nup != 0
        state[L] = "Up"
    end
    return state
end
function random_sdw_state(L, Nupd::Int, Ndnd::Int)
    state = fill("Emp", L)
    Nup_extra = Int.(Nupd % Ndnd)
    for i=1:2:(L-Nup_extra)
        state[i] = "Up"
        state[i+1] = "Dn"
    end
    if Nupd != 0
        state[L] = "Up"
    end
    return state
end


function random_ps_state(L, Nup, Ndn)
    state = fill("Emp", L)

    Ndbl = min(Nup, Ndn)
    Nup -= Ndbl
    Ndn -= Ndbl

    for i in 1:Ndbl
        state[i] = "UpDn"
    end

    next_site = Ndbl + 1

    if Nup > 0
        state[next_site] = "Up"
        next_site += 1
    elseif Ndn > 0
        state[next_site] = "Dn"
        next_site += 1
    end
    return state
end
#=
    Returns a state for the point (U, V) of the Extended Hubbard Model.
=#
function get_state(L, Nup::Int, Ndn::Int, U, V)
    state = fill("Emp", L)
    region = get_region_phase_diagram(U, V)
    if region == "METALLIC"
        state = random_metallic_state(L, Nup, Ndn)
    elseif region == "CDW"
        state = random_cdw_state(L, Nup, Ndn)
    elseif region == "SDW"
        state = random_sdw_state(L, Nup, Ndn)
    else
        state = random_ps_state(L, Nup, Ndn)
    end
    return state
end
#=
    Trying to separate the phase-diagram states not only
    in big squared blocks.
    See notebook ploting_ehm_diagram_selection_of_states.ipynb for an example of plot
    using this function. I try to made it look like the plot schematic
    diagram for the EHM.
=#
function get_region_phase_diagram(u, v)
    alpha = 0.5

    # Boundary functions
    get_metallic_lower_boundary(u) = u <= 0 ? -exp(alpha * u) : -alpha * u - 1.0
    get_metallic_upper_boundary_h(u) = -0.1 * u
    get_metallic_upper_boundary_negative(u) = -2.5 * alpha * u

    lower = get_metallic_lower_boundary(u)
    upper_h = get_metallic_upper_boundary_h(u)
    upper_neg = get_metallic_upper_boundary_negative(u)

    if v <= lower
        return "PS"
    elseif (u < 0 && v < 0 && v > lower) ||
           (u > 0 && v > 0 && v < upper_h) ||
           (u > 0 && v < 0 && v > lower && v < upper_neg)
        return "METALLIC"
    elseif (u < 0 && v > 0) ||
           (u > 0 && v > 0 && v >= 2 * u)
        return "CDW"
    elseif (u > 0 && v < 0 && v > lower) ||
           (u > 0 && v > 0 && v >= upper_h && v < 2 * u)
        return "SDW"
    end
    return "PS"
end
=#

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

function One_site_RDM(Psi,i)

    orthogonalize!(Psi,i)

    Ai = Psi[i]
    Ai_dag = dag(prime(Ai,"Site"))

    rho=Ai*Ai_dag

    return rho, reshape(Array(rho.tensor),4,4)
end

####################################################--- Two particle Reduced Density Matrix ---################################################

function Two_sites_RDM(Psi,sites,i,j)

    if i>j
        i,j = j,i
    end

    orthogonalize!(Psi,i)  # Move o centro de ortogonalidade para i 
    
    Psi_bra=prime(dag(Psi),linkinds(Psi)) # Calcula o dag

    rho=prime(Psi[i],linkinds(Psi,i-1))*prime(Psi_bra[i],sites[i])

    for k=(i+1):1:(j-1)

        rho*=Psi[k]
        rho*=Psi_bra[k]

    end

    rho*=prime(Psi[j],linkinds(Psi,j))*prime(Psi_bra[j],sites[j])

    return rho, RDM_to_matrix(rho)
end

function RDM_to_matrix(rho)
    inds_rho = inds(rho)

    s_ket = Index[]
    s_bra = Index[]

    for i in inds_rho
        if i == noprime(i)
            push!(s_ket, i)
        else
            push!(s_bra, i)
        end
    end

    Cket = combiner(s_ket...)
    Cbra = combiner(s_bra...)

    M = rho * Cket * Cbra

    # pega os índices combinados
    i_ket = inds(M)[1]
    i_bra = inds(M)[2]

    return Array(M, i_ket, i_bra)
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

function Two_particle_RDM2(Psi,sites)

    L = length(sites)
    spins = ["up","dn"]
    spin_val = Dict("up"=>1,"dn"=>-1)

    rho = zeros(ComplexF64,2L,2L,2L,2L)

    for i in 1:L-1, j in i+1:L
    for k in 1:L-1, l in k+1:L

        # Hermiticidade: só metade

        if (i,j) > (k,l)
            continue
        end

        for (is,si) in enumerate(spins), (js,sj) in enumerate(spins)
        for (ks,sk) in enumerate(spins), (ls,sl) in enumerate(spins)

            # conservação de spin

            if spin_val[si] + spin_val[sj] != spin_val[sk] + spin_val[sl]
                continue
            end

            val = Four_op_correlation(Psi,i,j,k,l,si,sj,sk,sl,sites)

            a = (i-1)*2 + is
            b = (j-1)*2 + js
            c = (k-1)*2 + ks
            d = (l-1)*2 + ls

            rho[a,b,c,d] = val

            # antissimetria

            rho[b,a,c,d] = -val
            rho[a,b,d,c] = -val
            rho[b,a,d,c] = val

            # hermiticidade

            rho[c,d,a,b] = conj(val)
            rho[d,c,a,b] = -conj(val)
            rho[c,d,b,a] = -conj(val)
            rho[d,c,b,a] = conj(val)

        end
        end

    end
    end

    return 2/(L*(L-1)).*reshape(rho,(2*L)^2,(2*L)^2)
end

function Two_particle_RDM3(Psi,sites)

    L = length(sites)
    spins = ["up","dn"]
    spin_val = Dict("up"=>1,"dn"=>-1)

    rho = zeros(ComplexF64,(2L)^2,(2L)^2)

    for i in 1:L-1, j in i+1:L
    for k in 1:L-1, l in k+1:L

        # Hermiticidade: só metade

        if (i,j) > (k,l)
            continue
        end

       #println(i," ",j," ",k," ",l)

        for (is,si) in enumerate(spins), (js,sj) in enumerate(spins)
        for (ks,sk) in enumerate(spins), (ls,sl) in enumerate(spins)
            # conservação de spin

            if spin_val[si] + spin_val[sj] != spin_val[sk] + spin_val[sl]
                continue
            end

            val = Four_op_correlation(Psi,i,j,k,l,si,sj,sk,sl,sites)

            a=pair_index(i,si,j,sj,L)
            b=pair_index(k,sk,l,sl,L)
            
            a_=pair_index(j,sj,i,si,L)
            b_=pair_index(l,sl,k,sk,L)

            rho[a,b] = val
            rho[b,a] = conj(val)

            rho[a_,b] = -val
            rho[b,a_] = -conj(val)

            rho[a,b_] = -val
            rho[b_,a] = -conj(val)

            rho[a_,b_] = val
            rho[b_,a_] = conj(val)

        end
        end

    end
    end

    for i in 1:L

    for (si_idx,si) in enumerate(spins), (sj_idx,sj) in enumerate(spins)
    for (sl_idx,sl) in enumerate(spins), (sk_idx,sk) in enumerate(spins)
        # só spins diferentes
        

        if spin_val[si] + spin_val[sj] != spin_val[sk] + spin_val[sl]
            continue
        end

        if si == sj || sk == sl
            continue
        end

        # conservação de spin automática aqui
        val = Four_op_correlation(Psi,i,i,i,i,si,sj,sk,sl,sites)

        a = pair_index(i,si,i,sj,L)
        b = pair_index(i,sk,i,sl,L)

        rho[a,b] = val
        rho[b,a] = conj(val)
    
        #println(a," ",b," ",val)
    end
    end
end

    return 2/(L*(L-1)).*rho
end

function pair_index(i, si, j, sj, L)

    # índice spin-orbital
    a = 2*(i-1) + (si == "up" ? 1 : 2)
    b = 2*(j-1) + (sj == "up" ? 1 : 2)

    N = 2*L

    return (a-1)*N + b
end

function Nup_Ndn(Psi, sites,i,j)

    L = length(sites)
    val = 0.0

        val = inner(Psi,
            apply(op(sites,"Nup",i),
            apply(op(sites,"Ndn",j), Psi)))


    return val
end

function Ndn_Nup(Psi, sites,i,j)

    L = length(sites)
    val = 0.0

        val = inner(Psi,
            apply(op(sites,"Ndn",i),
            apply(op(sites,"Nup",j), Psi)))


    return val
end

function Four_op_correlation(Psi,i,j,k,l,si,sj,sk,sl,sites)

    O = OpSum()
    O += 1.0,"Cdag"*si, i, "Cdag"*sj, j,"C"*sl,  l,"C"*sk,  k # Se trocar a ordem os autovaores fica positivo

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

function ni_nj(psi,sites,i,j)

    ampo = AutoMPO()
    add!(ampo, 1.0, "Nup", i, "Ndn", j)

    Op = MPO(ampo, sites)

    val = inner(psi, Op, psi)

    return val
end


function Si_Sj(psi, sites, i, j)
    ampo = AutoMPO()
    add!(ampo, 1.0, "Sz", i, "Sz", j)
    Op = MPO(ampo, sites)

    return inner(psi, Op, psi)
end

####################################################--- von Neumann Entropy ---################################################

function S_vNeumann(rho)

    lambs_pre=real.(eigvals(rho))

    lambs=lambs_pre[lambs_pre.>0]

    return -1.0*sum(lambs.*log.(lambs))
end


####################################################--- Rugosity ---################################################

function Rugosity_2sites(Psi,sites,i_center,range_,fac;vN=true)

    N=length(Psi)

    Rug=Float64[]
    Entro=Float64[]

    for d=(i_center-range_):1:(i_center+range_)

        if d!=i_center

            _,rho=Two_sites_RDM(Psi,sites,i_center,d)

            push!(Rug,sum(rho)/fac)

            if vN 
                
                push!(Entro,S_vNeumann(rho,N))

            end
            
        end
    end

    return Rug, Entro
end



function Rugosity_2sites2(Psi,sites,i_center,range_;vN=true)

    N=length(Psi)

    Rug1=ComplexF64[]
    Entro1=ComplexF64[]
    Rug2=ComplexF64[]
    Entro2=ComplexF64[]

    for d=(i_center-range_):1:(i_center+range_)

        if d!=i_center

            println(i_center, " ",d)

            rho,rho_mat=Two_sites_RDM(Psi,sites,i_center,d)

            psi=T_less_2sites(sites,i_center,d)

            psiT=psi[1]*psi[2]

            psiT_ket = dag(psiT) 
            psiT_bra = prime(psiT)

            push!(Rug1,(psiT_bra * (rho * psiT_ket))[])


            if vN 
                
                push!(Entro1,S_vNeumann(rho_mat))

            end
            
        end
    end

    i_center+=1
    
    for d=(i_center-range_):1:(i_center+range_)

        if d!=i_center

            println(i_center, " ",d)

            rho,rho_mat=Two_sites_RDM(Psi,sites,i_center,d)

            psi=T_less_2sites(sites,i_center,d)

            psiT=psi[1]*psi[2]

            psiT_ket = dag(psiT) 
            psiT_bra = prime(psiT)

            push!(Rug2,(psiT_bra * (rho * psiT_ket))[])


            if vN 
                
                push!(Entro2,S_vNeumann(rho_mat))

            end
            
        end
    end


    return (Rug1.+reverse(Rug2))./2, (Entro1.+reverse(Entro2))./2
end

####################################################--- Testure less state ---################################################

function T_less_2sites(sites,i,j)

    if i>j
        i,j = j,i
    end

    sites_ij=[sites[i],sites[j]]

    # Duble occupation states

    psi_vac_dub=productMPS(sites_ij,["UpDn","Emp"])
    psi_dub_vac=productMPS(sites_ij,["Emp","UpDn"])

    # Singlet and triplet states

    psi_up_dn=productMPS(sites_ij,["Up","Dn"])
    psi_dn_up=productMPS(sites_ij,["Dn","Up"])
    psi_up_up=productMPS(sites_ij,["Up","Up"])
    psi_dn_dn=productMPS(sites_ij,["Dn","Dn"])

    psi_S=(psi_up_dn-psi_dn_up)/sqrt(2)
    psi_T0=(psi_up_dn+psi_dn_up)/sqrt(2)
    psi_Tp=psi_up_up
    psi_Tm=psi_dn_dn

    psi=(psi_dub_vac+psi_vac_dub+psi_Tp+psi_Tm+psi_S+psi_T0)/sqrt(6)

    return psi
end

####################################################--- Gate Construction ---################################################

function Gates_list(dt,L,Ut,Vt,sites)

    gates = ITensor[]

    # -------- odd bonds --------
    for i=1:2:(L-1)
        h=-1.0*op("Cdagup", sites[i])*op("Cup",sites[i+1])-1.0*op("Cdagup", sites[i+1])*op("Cup",sites[i])-1.0*op("Cdagdn", sites[i])*op("Cdn",sites[i+1])-1.0*op("Cdagdn", sites[i+1])*op("Cdn",sites[i])+Vt*op("Ntot", sites[i])*op("Ntot",sites[i+1])
        push!(gates,exp(-dt/2*h))
    end

    # -------- even bonds --------
    for i=2:2:(L-1)
         h=-1.0*op("Cdagup", sites[i])*op("Cup",sites[i+1])-1.0*op("Cdagup", sites[i+1])*op("Cup",sites[i])-1.0*op("Cdagdn", sites[i])*op("Cdn",sites[i+1])-1.0*op("Cdagdn", sites[i+1])*op("Cdn",sites[i])+Vt*op("Ntot", sites[i])*op("Ntot",sites[i+1])
        push!(gates,exp(-dt*h))
    end

    # -------- odd bonds again --------
    for i=1:2:(L-1)
         h=-1.0*op("Cdagup", sites[i])*op("Cup",sites[i+1])-1.0*op("Cdagup", sites[i+1])*op("Cup",sites[i])-1.0*op("Cdagdn", sites[i])*op("Cdn",sites[i+1])-1.0*op("Cdagdn", sites[i+1])*op("Cdn",sites[i])+Vt*op("Ntot", sites[i])*op("Ntot",sites[i+1])
        push!(gates,exp(-dt/2*h))
    end

    # -------- onsite --------

    for i=1:(L-1)
        h=Ut*op("Nupdn", sites[i])*op("Id",sites[i+1])
        push!(gates,exp(-dt*h))
    end

    h=Ut*op("Id", sites[end-1])*op("Nupdn",sites[end])
    push!(gates,exp(-dt*h))

    return gates
end

####################################################--- Time evolution ---################################################

function time_evolve_TEDB(psi,dt,U0,V0,tauU,tauV,sites,tf,Nup,Ndn,Nsweep,maxD_dmrg,cutoff_drmg; cutoff=1e-8, maxdim=500)
    
    L=length(sites)
    Nsteps=Int(tf/dt)
    Tau=zeros(Float64,Nsteps)
    Utau=zeros(Float64,Nsteps)
    Vtau=zeros(Float64,Nsteps)
    fidelity=zeros(Float64,Nsteps)
    fidelity_loc=zeros(Float64,Nsteps)
    Rug=zeros(Float64,Nsteps)
    Entro=zeros(Float64,Nsteps)
    psi0=copy(psi)

    println("Time Evolution: Starting")

    for step=1:1:Nsteps

        println("Evolution step: ",step,"/",Nsteps," (t=",step*dt,")")

        Ut=ramp(U0,1,tauU,step*dt)
        Vt=ramp(V0,1,tauV,step*dt)

        Vtau[step]+=Vt 
        Utau[step]+=Ut

        gates=Gates_list(-1im*dt,L,Ut,Vt,sites)

        @time psi=apply(gates,psi;cutoff,maxdim)#, returninfo=true)

        println("Maximum Bound Dimension: ",maximum(i -> linkdim(psi, i), 1:length(psi)-1))

        normalize!(psi)

        Ht=MPO_construction(L,sites,Ut,Vt)

        _,psi_loc=DMRG_optm(Ut,Vt,L,Nup,Ndn,Ht,sites,Nsweep;maxD=maxD_dmrg,cutoff=cutoff_drmg)

        Tau[step]+=step*dt
        fidelity[step]+=abs2(inner(psi0,psi))
        fidelity_loc[step]+=abs2(inner(psi,psi_loc))

        r,e=Rugosity_2sites2(psi,sites,Int(L/2),1;vN=true)

        Rug[step]+=real(r[1]) 
        Entro[step]+=real(e[1])

        println("-------------------")
        println("Fidelity: ",fidelity[step])
        println("Fidelity inst.: ",fidelity_loc[step])
        println("Rugosity: ",r[1] )
        println("Entropy: ",e[1])
        println("###############################################################")

    end 

    println("Time Evolution: Finished")

    return Tau, fidelity,fidelity_loc, Rug, Entro
end


function time_evolve_TDVP(psi,dt,U0,V0,tauU,tauV,sites,tf; cutoff=1e-8, maxdim=500)
    
    L=length(sites)
    Nsteps=Int(tf/dt)
    Tau=Float64[]
    fidelity=Float64[]
    psi0=copy(psi)

    println("Time Evolution: Starting")

    for step=1:1:Nsteps

        println("Evolution step: ",step,"/",Nsteps," (t=",step*dt,")")

        println("Hamiltonian construction: Starting")

        Ut=ramp(U0,1,tauU,step*dt)
        Vt=ramp(V0,1,tauV,step*dt)

        @time Ht=MPO_construction(L,sites,Ut,Vt)

        println("Hamiltonian construction: Finished")
        println("------------------------------------")

        println("TDVP step: Starting")

        @time psi=tdvp(H,-1im*dt,psi;maxdim,cutoff)

        println("TDVP step: Finished")
        
        println("Maximum Bound Dimension: ",maximum(i -> linkdim(psi, i), 1:length(psi)-1))

        normalize!(psi)

        push!(Tau,step*dt)
        push!(fidelity,abs2(inner(psi0,psi)))

        println("###############################################################")

    end 

    println("Time Evolution: Finished")

    return Tau, fidelity
end


function ramp(a,b,tau,t)

    if tau==0

        return a 
    else

        return a+b*t/tau
    end
end

####################################################--- Main Function ---################################################

function Main_func_cluster(U,V,N,Nsweep,maxD_dmrg,cutoff_dmrg,Spintot)

    Npart=floor(Int,N/2)

    if Spintot=="0"
        Nup=Npart+N%2
        Ndn=N-Nup
    elseif Spintot=="1"
        Nup=Npart+N%2+1
        Ndn=N-Nup
    end

    sites=siteinds("Electron", N; conserve_qns = true)

    H=MPO_construction(N,sites,U,V)

    E0,Psi=DMRG_optm(U,V,N,Nup,Ndn,H,sites,Nsweep;maxD=maxD_dmrg,cutoff=cutoff_dmrg)

    _,rho1=One_site_RDM(Psi,Int(N/2))

    C1_l=sum(abs.(rho1))-tr(abs.(rho1))

    Entro_1=S_vNeumann(rho1)

    Rug,Entro_2=Rugosity_2sites2(Psi,sites,Int(N/2),3;vN=true)

    return E0, Entro_1, C1_l, Entro_2, Rug, expect(Psi,"Ntot")
end

####################################################--- Main Function Time ---################################################

function Main_func_time(U0,V0,N,Nsweep,Spintot,tauU,tauV,tf,dt,cutoff_time,maxD_time,cutoff_drmg,maxD_dmrg)

    Npart=floor(Int,N/2)

    if Spintot=="0"
        Nup=Npart+N%2
        Ndn=N-Nup
    elseif Spintot=="1"
        Nup=Npart+N%2+1
        Ndn=N-Nup
    end

    sites=siteinds("Electron", N; conserve_qns = true)

    H=MPO_construction(N,sites,U0,V0)

    E0,Psi=DMRG_optm(U0,V0,N,Nup,Ndn,H,sites,Nsweep;maxD=maxD_dmrg,cutoff=cutoff_drmg)

    Tau,fidelity,fidelity_loc,Rug,Entro=time_evolve_TEDB(Psi,dt,U0,V0,tauU,tauV,sites,tf,Nup,Ndn,Nsweep,maxD_dmrg,cutoff_drmg;cutoff=cutoff_time,maxdim=maxD_time)
    
    return Tau, fidelity,fidelity_loc, Rug, Entro
end
####################################################--- Saving Function ---################################################

function saving_file(file_name,dados;file_=true)

    if file_

        file1=open(file_name,"w")  # file's name -> W value, r means random, 1 is the delta's value and 100 is the number of points of time, essemble's size
            
        for j=1:1:length(dados[:,1])

            for i=1:1:(length(dados[j,:])-1)

                write(file1,string(dados[j,i])," ")

            end

            write(file1,string(dados[j,end]),"\n")

        end

        close(file1)

    end
    
end
