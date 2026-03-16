##########################################################################################################################

using ITensors, ITensorMPS, LinearAlgebra, Base.Threads,ITensors.SiteTypes

include("Extended_Hubbard_Model_functions.jl")

############################################################################################################################

const N=10        # Number of sites
const Nsweep=20
const cutoff=1e-8
const maxD=[2,2,2,2,2,8,8,8,8,8,20,20,20,20,20,50,50,50,50,50,100,100,200,200,400]

U=1.0
V=1.0

Npart=floor(Int,N/2)
Nup=Npart+N%2
Ndn=N-Nup
    
sites=siteinds("Electron", N; conserve_qns = true)

H=MPO_construction(N,sites,U.*ones(N),V.*ones(N),-1.0.*ones(N))

E0,Psi=DMRG_optm(N,Nup,Ndn,H,sites,Nsweep,maxD,cutoff)