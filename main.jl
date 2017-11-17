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
   IndexRow::Vector{Int32}
   IndexColumn::Vector{Int32}
   IndexRowOR::Vector{Int32}
   IndexColumnOR::Vector{Int32}
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
   TimeG::Float64       = 0.0
   TimeCons::Float64    = 0.0
   AVG                  = 0
   N::Int32             = 100
   Ngen::Int32          = 30
   Grasp1::Grasp        = Grasp([0.50,0.65,0.75,0.80,0.90],[0.20,0.20,0.20,0.20,0.20])
   filename = string("pb_2000rnd0100.dat")
   Pb::Problem          = ReadFile(string("./Data/",filename))
   println("Beginning evolution")
   println("Problem is ",filename)
   max = 0
   min = typemax(Int32)
   sv = Vector{Vector{Genome}}(Ngen)
   ngen = 0
   maxit = 5
   println("\\begin{tabular}{|c|c|c|c|c|c|c|c|c|}")
   println("Problem & \\multicolumn{7}{|c|}{Genetic algorithm} \\\\")
   println(" & \\multicolumn{2}{|c|}{Parameter} & \$ \\tilde{Z_{min }} \$ & \$ \\tilde{Z_{avg}} \$ & \$ \\tilde{Z_{max}} \$ & Time GRASP & Time Genetic \\\\")
   for g::Int32 in [50,100]
      for l in [0.01,0.05]
         max = 0
         min = typemax(Int32)
         AVG                  = 0
         TimeCons = 0.0
         TimeG = 0.0
         for i = 1:1:maxit
            TimeCons += @elapsed Pb,Population::Vector{Genome} = InitPopulation(Pb,N,Grasp1)
            TimeG += @elapsed Population = Evolution(Population,Pb,g,Ngen,l)

            AVG+=Population[1].CurrentObjectiveValue
            if Population[1].CurrentObjectiveValue < min
               min = Population[1].CurrentObjectiveValue
            end
            if Population[1].CurrentObjectiveValue > max
               max = Population[1].CurrentObjectiveValue
            end

         end
         println("\\hline")
         println(filename[4:length(filename)]," & ",g," & ",l," & ",min," & ",round(AVG/maxit,2)," & ",max," & ", round(TimeCons/maxit,2)," & ",round(TimeG/maxit,2),"\\\\")
      end
   end
   println("\\hline")
   println("\\end{tabular}")
   #PlotGeneticAlgorithm(sv,N,ngen,filename)
end

main()
