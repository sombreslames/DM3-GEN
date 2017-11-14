#Ferrari Leon
#M1 ORO
#Version de test
#Julia JuMP
#DM1 - Metaheuristiques
using JuMP, GLPKMathProgInterface,Gadfly,DataFrames
#tHIS STRUCT DESCRIBE THE COMMON CHARACTERISTICS OF THE PROBLEM
# IN ORDER TO USE SMALLER INDIVIDUAL I USE PROBLEM AND INDIVIDUAL
type Problem
   NBvariables::Int32
   NBconstraints::Int32
   Variables::Vector{Int32}
   LeftMembers_Constraints::SparseMatrixCSC{Int32,Int32}
   Utility::Array{Float64,2}
   Bonus::Int64
   SumObj::Int64
   MinObj::Int64
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
RdSeed = RandomDevice()
include("heuristics_spp.jl")
include("ToolsAndPlot.jl")
include("Genetic.jl")
function main()
   Time::Float64        = 0.0
   AVG                  = 0
   N::Int32             = 100
   Ngen::Int32          = 30
   Grasp1::Grasp        = Grasp([0.50,0.65,0.75,0.80,0.90],[0.20,0.20,0.20,0.20,0.20])
   filename = string("pb_100rnd0100.dat")
   Pb::Problem          = ReadFile(string("./Data/",filename))
   println("Beginning evolution #IDONTBELIEVEINIT")
   max = 0
   min = typemax(Int32)
   sv = Vector{Vector{Genome}}(Ngen)
   ngen = 0
   maxit = 5
   for i = 1:1:maxit
      tic()
      Pb,Population::Vector{Genome} = InitPopulation(Pb,N,Grasp1)
      Population = Evolution(Population,Pb,N,Ngen,0.01)
      Time+= toc()
      AVG+=Population[1].CurrentObjectiveValue
      if Population[1].CurrentObjectiveValue < min
         min = Population[1].CurrentObjectiveValue
      end
      if Population[1].CurrentObjectiveValue > max
         max = Population[1].CurrentObjectiveValue
      end

   end
   println("Max : ",max ," | Min : ",min)
   println("Average time : ", round(Time/maxit,2))
   println("Average known value is ",round(AVG/maxit,2))
   #PlotGeneticAlgorithm(sv,N,ngen,filename)
end

main()
