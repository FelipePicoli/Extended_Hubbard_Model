############################################################################################################################
###########################################--- Extended Hubbard Model ---###################################################
############################################################################################################################


############################################################################################################################
###################################################--- Parameters ---#######################################################
############################################################################################################################

using ITensors, ITensorMPS, LinearAlgebra, Base.Threads,ITensors.SiteTypes

include("Extended_Hubbard_Model_functions.jl")

############################################################################################################################

const N=64        # Number of sites
const Nsweep=50
const cutoff=1e-8
const maxD=[2,2,2,2,2,8,8,8,8,8,20,20,20,20,20,50,50,50,50,50,100,100,200,200,400]

const Ul=LinRange(-6,6,30)
const Vl=LinRange(-3,3,30)

results=zeros(Float64,length(Ul)*length(Vl),9)

Threads.@threads  for iu=1:1:length(Ul)
    for iv=1:1:length(Vl)

        idx = (iu-1)*length(Vl) + iv

        E0,s1,c1,s21,c21,s22,c22=Main_func_cluster(Ul[iu],Vl[iv],N,Nsweep,maxD,cutoff)

        results[idx,1]+=Ul[iu]
        results[idx,2]+=Vl[iv]
        results[idx,3]+=E0
        results[idx,4]+=s1
        results[idx,5]+=c1
        results[idx,6]+=s21
        results[idx,7]+=c21
        results[idx,8]+=s22
        results[idx,9]+=c22

    end
end

file=open("EHM_N"*string(N)*"_lengU"*string(length(Ul))*"_lengV"*string(length(Vl))*"_Nsweep"*string(Nsweep)*".txt","w")  # file's name -> W value, r means random, 1 is the delta's value and 100 is the number of points of time, essemble's size
    
for i=1:1:size(results,1)
        
    write(file,string(results[i,1])," ",string(results[i,2])," ",string(results[i,3])," ",string(results[i,4])," ",string(results[i,5])," ",string(results[i,6])," ",string(results[i,7])," ",string(results[i,8])," ",string(results[i,9]),"\n")
                 
end

close(file)
