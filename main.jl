#Ferrari Leon
#M1 ORO
#Version de test
#Julia JuMP
#DM1 - Metaheuristiques
using JuMP, GLPKMathProgInterface,PyPlot
#tHIS STRUCT DESCRIBE THE COMMON CHARACTERISTICS OF THE PROBLEM
# IN ORDER TO USE SMALLER INDIVIDUAL I USE PROBLEM AND INDIVIDUAL
type Problem
   NBvariables::Int32
   NBconstraints::Int32
   Variables::Vector{Int32}
   LeftMembers_Constraints::SparseMatrixCSC{Int32,Int32}
   Utility::Array{Float64,2}
   Bonus::Int64
end
#Individual describe a possible solution, it will be used as an individu for the population of the genetic algorithm
#Solution is the variable containing the solution
#It will be used during the crossover the mutation etc
type Individual
   CurrentObjectiveValue::Int32
   Solution::Vector{Int32}
   CurrentVarUsed::Vector{Int32}
   LastRightMemberValue_Constraint::Vector{Int32} #Useless with the freedom vector
   Freedom::Vector{Int32}
end
#Individual two --> Containing only the solution
#We will have to make our solution feasible at each move
type Genome
   CurrentObjectiveValue::Int32
   Solution::Vector{Int32}
end
type Grasp
   Values::Vector{Float64}
   Probability::Vector{Float64}
end
include("heuristics_spp1.jl")
include("ToolsAndPlot.jl")
include("genetic1.jl")
function main()
   N::Int32             = 100
   Ngen::Int32          = 30
   Pb::Problem          = ReadFile("./Data/pb_100rnd0100.dat")
   Grasp1::Grasp        = Grasp([0.50,0.65,0.75,0.80,0.90],[0.20,0.20,0.20,0.20,0.20])
   Pb,Population::Vector{Genome},SumObj::Int32 = InitPopulation(Pb,N,Grasp1)
   #Evolution(Population,BPP,100,Grasp1,SumObj)
   #@time Evolution(Population,Pb,N,Ngen,Grasp1,SumObj)
   #println(Population[1].CurrentObjectiveValue)
end

main()
