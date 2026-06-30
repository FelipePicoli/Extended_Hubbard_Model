############################################################################################################################
###########################################--- Extended Hubbard Model ---###################################################
############################################################################################################################


############################################################################################################################
###################################################--- Parameters ---#######################################################
############################################################################################################################

using ITensors, ITensorMPS, LinearAlgebra, Base.Threads,ITensors.SiteTypes

include("Extended_Hubbard_Model_functions.jl")
-
############################################################################################################################
file_=true
N=64        # Number of sites
Spintot="0"
Nsweep=50
cutoff=1e-30
maxD=[20,20,20,20,20,30,30,30,30,30,30,40,40,40,40,40,50,50,50,50,50,55,55,55,55,55,55,55,60,65,65,65,70,70,70,70,75,80,85,85,90,90,100,150,200,400,400,300]
 #maxD=[2,2,2,2,2,8,8,8,8,8,20,20,20,20,20,50,50,50,50,50,100,100,200,200,400,600,800,800,800,1200,1200,1200,1600,1600,1600,2000,2000,2000][20,20,50,50,100,100,200,200,400,400,800,800,1200,1200,1600]#

Ul=[1.6]#LinRange(-6,6,30)
Vl=LinRange(-3.0,+3.0,50)

resultsS=zeros(Float64,length(Ul)*length(Vl),5+6)
resultsR=zeros(Float64,length(Ul)*length(Vl),3+6)

for iu=1:1:length(Ul)
    for iv=1:1:length(Vl)

        println(iu," => ",Ul[iu],", ",iv," => ",Vl[iv])

        idx = (iu-1)*length(Vl) + iv

        @time E0,s1,c1,s2,rug,n=Main_func_cluster(Ul[iu],Vl[iv],N,Nsweep,maxD,cutoff,Spintot)

        resultsS[idx,1]=Ul[iu]
        resultsS[idx,2]=Vl[iv]
        resultsS[idx,3]=E0
        resultsS[idx,4]=s1
        resultsS[idx,5]=c1
        resultsS[idx,6:end]=s2

        resultsR[idx,1]=Ul[iu]
        resultsR[idx,2]=Vl[iv]
        resultsR[idx,3]=E0
        resultsR[idx,4:end]=rug
        

    end
end

saving_file("3teste_Entro_EHM_N"*string(N)*"_lengU"*string(length(Ul))*"_lengV"*string(length(Vl))*"_Nsweep"*string(Nsweep)*"_Spintot"*string(Spintot)*".txt",resultsS;file_) 

saving_file("3teste_Rug_EHM_N"*string(N)*"_lengU"*string(length(Ul))*"_lengV"*string(length(Vl))*"_Nsweep"*string(Nsweep)*"_Spintot"*string(Spintot)*".txt",resultsR;file_) 
