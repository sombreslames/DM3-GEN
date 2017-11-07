function InitPopulation(Pb::Problem,N::Int32,Grasp::Grasp)
   println("####################### Problem Info #######################")
   println("# Var : ",Pb.NBvariables, " | Constraints : ", Pb.NBconstraints)
   Pb.Bonus,NB,Pb = ComputeUtilityPlus(Pb)
   println("# f(x): ", Pb.Bonus, " | Variable considered : ",NB)
   println("############################################################")
   Max::Int32      = 0
   Min::Int32      = typemax(Int32)
   G::Int32        = length(Grasp.Values)
   Average     = zeros(Float64,G)
   NBruns      = zeros(Float64,G)
   Population  = Array{Genome}(N)
   for i = 1:1:N
      IndivTemp         = Individual(Pb.Bonus,zeros(Int64,Pb.NBvariables),zeros(Int64,0),zeros(Int64,Pb.NBconstraints),zeros(Int64,Pb.NBvariables))
      IndexAlpha,Alpha  = ReactiveGrasp(Grasp)
      IndivTemp         = GraspConstruction(Pb,IndivTemp,Alpha)
      Population[i]     = Genome(IndivTemp.CurrentObjectiveValue,IndivTemp.solution)
      if Population[i].CurrentObjectiveValue > Max
         Max = Population[i].CurrentObjectiveValue
      elseif  Population[i].CurrentObjectiveValue < Min
         Min = Population[i].CurrentObjectiveValue
      end
      Average[IndexAlpha] += Population[i].CurrentObjectiveValue
      NBruns[IndexAlpha] += 1
   end
   TotSum = sum(Average)
   println("###################### Population Info #####################")
   println("# Moyenne : ",TotSum/sum(NBruns))
   println("# MIN : ", Min, " | MAX : ",Max)
   println("############################################################")
   for j in 1:1:G
      Average[j] = Average[j] / NBruns[j]
   end
   Population = sort(Population,by=a->a.CurrentObjectiveValue,rev=true)
   return Pb,Population,TotSum
end

function Evolution(Population::Vector{Individual},Pb::Problem,NPop::Int32,Ngen::Int32,Grasp::Grasp,SumObj::Int32)
   for i  = 1:1:Ngen
      for j = 1:1:NPop
         p1,p2 = BinaryTourmanent(Population,NPop,false,SumObj)
         child = CrossoverMethod1(Pb,Population[p1],Population[p2])
         InsertAndReplace(Population,child)
         ProbMut = rand()
         if ProbMut > 0.8
            Population[j] = Mutation(Pb,Population[j])
         end
      end
   end


end
function InsertAndReplace(Population::Vector{Individual},Indi::Individual)
   index = searchsortedfirst(Population,Indi,by=x->x.CurrentObjectiveValue,rev=true)
   insert!(Population,index,Indi)
   PopSize = length(Population)
   deleteat!(Population,PopSize)
end

function Mutation(Pb::Problem,Indi::Genome)
   VarToZero      = rand(Indi.CurrentVarUsed)
   #println(Indi.CurrentVarUsed)
   Ok,Indi        = SetToZero(Pb,Indi,VarToZero)
   for i in 1:1:Pb.NBvariables
      Index = convert(Int,Pb.Utility[1,i])
      if Indi.Freedom[Index] == 0 && Indi.Solution[Index] == 0
            Ok,Indi = SetToOne(Pb,Indi,Index)
            if Ok
               #println("Mutated successfully")
            end
      end
   end
   return Indi
end
function RouletteSelection(Population::Vector{Genome},N::Int32,SumObj::Int32)
   Value::Float64 =  0
   p1::Int32 = 0 ; p2::Int32=0
   Rand1 = rand()
   Rand2 = rand()
   for i = 1:1:N
      Value += Population[i].CurrentObjectiveValue/SumObj
      if Rand1 < Value
         p1 = i
      end
      if Rand2 < Value
         p2 = i
      end
      if p1 != 0 && p2 != 0
         return p1,p2
      end
   end
end
function ParentSelection(Population::Vector{Genome},N::Int32)
   Parent1::Int32     = rand(1:100)
   Parent2::Int32     = rand(1:100)
   different::Bool         = true
   while different
      if Parent1 != Parent2
         different   = false
      else
         Parent1     = rand(1:100)
         Parent2     = rand(1:100)
      end
   end
   return Parent1,Parent2
end
function BinaryTourmanent(Population::Vector{Genome},N::Int32,Mode::Bool,SumObj::Int32)
   #Mode = true ==> Same Probability for each individual
   p1::Int32=0;  p2::Int32=0
   p3::Int32=0;  p4::Int32=0
   if Mode
      p1,p2 = ParentSelection(Population,N)
      p3,p4 = ParentSelection(Population,N)
   else
      p1,p2 = RouletteSelection(Population,N,SumObj)
      p3,p4 = RouletteSelection(Population,N,SumObj)
   end
   if Population[p1].CurrentObjectiveValue > Population[p2].CurrentObjectiveValue
      if Population[p3].CurrentObjectiveValue > Population[p4].CurrentObjectiveValue
         return p1,p3
      else
         return p1,p4
      end
   end
   if Population[p3].CurrentObjectiveValue > Population[p4].CurrentObjectiveValue
      return p2,p3
   end
   return p2,p4
end
#  Child[i] = parent1[i] if mask[i]=1
#  Child[i] = parent2[i] if mask[i]=0
function CrossoverMethod1(Pb::Problem,Parent1::Genome,Parent2::Genome)
   Child = Individual(Pb.Bonus,zeros(Int64,Pb.NBvariables),zeros(Int64,0),zeros(Int64,Pb.NBconstraints),zeros(Int64,Pb.NBvariables))
   Mask = rand(0:1,Pb.NBvariables)
   for i in 1:1:Pb.NBvariables
      if Mask[i] == 1 && Parent1.Solution[i] == 1
         ChildTemp.Solution[i] = 1
         ChildTemp.CurrentObjectiveValue += Pb.Variables[i]
      elseif Mask[i] == 0 && Parent2.Solution[i] == 1
         ChildTemp.Solution[i] = 1
         ChildTemp.CurrentObjectiveValue += Pb.Variables[i]
      end
   end
   #=println("###################### Crossover Info ######################")
   println("# Parent 1 : ",Parent1.CurrentObjectiveValue," | Parent 2 :",Parent2.CurrentObjectiveValue)
   println("# Child : ",Child.CurrentObjectiveValue)
   println("############################################################")
   =#return Child
end
#   println("# Parent 2 : ",Parent2.CurrentObjectiveValue,"\n","# ",Parent2.Solution)
