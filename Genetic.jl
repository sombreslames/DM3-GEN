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
      Population[i]     = Genome(IndivTemp.CurrentObjectiveValue,IndivTemp.Solution)
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
   #GRASP AVERGARDE VALUE
   for j in 1:1:G
      Average[j] = Average[j] / NBruns[j]
   end
   Population = sort(Population,by=a->a.CurrentObjectiveValue,rev=true)
   Pb.MinObj = Min
   Pb.SumObj = TotSum
   return Pb,Population
end

function Evolution(Population::Vector{Genome},Pb::Problem,NPop::Int32,Ngen::Int32,Grasp::Grasp)
   rng = MersenneTwister(1234);
   for i  = 1:1:Ngen
      for j = 1:1:NPop
         p1,p2 = BinaryTourmanent(Population,NPop,false,Pb)
         #here repair and mutate
         child1,child2 = CrossoverMethod1(Pb,Population[p1],Population[p2],Pb.SumObj)
         InsertAndReplace(Population,child)
         ProbMut = rand()
         if ProbMut > 0.8
            Population[j] = Mutation(Pb,Population[j])
         end
      end
   end


end


function RepairAndMutation(Pb::Problem,Indi::Genome,RandSeed::MersenneTwister)
   RanPermut = randperm(RandSeed, Pb.NBconstraints)
   for i in eachindex(RanPermut)
      for j in 1:1:Pb.NBconstraints
         if Pb.LeftMembers_Constraints[RanPermut[i],j] ==1

         end
      end
   end
   return Indi
end
function RouletteSelection(Population::Vector{Genome},N::Int32,Pb::Problem)
   Value::Float64 =  0
   p1::Int32 = 0 ; p2::Int32=0
   Rand1 = rand()
   Rand2 = rand()
   TotDec = Pb.SumObj - (N*Pb.MinObj)
   for i = 1:1:N
      Value += (Population[i].CurrentObjectiveValue-Pb.MinObj)/TotDec
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
function BinaryTourmanent(Population::Vector{Genome},N::Int32,Pb::Problem)
   #Mode = true ==> Same Probability for each individual
   p1::Int32=0;  p2::Int32=0
   p3::Int32=0;  p4::Int32=0

   p1,p2 = RouletteSelection(Population,N,Pb)
   p3,p4 = RouletteSelection(Population,N,Pb)

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
function CrossoverMethod(Pb::Problem,Parent1::Genome,Parent2::Genome)
   Child1 = Individual(Pb.Bonus,zeros(Int64,Pb.NBvariables))
   Child2 = Individual(Pb.Bonus,zeros(Int64,Pb.NBvariables))
   #Soyons malin utilisons les probabilites
   for i in 1:1:Pb.NBvariables
      if Parent1.Solution[i]  == Parent1.Solution[2]
         Child1.Solution[i] = Parent1.Solution[i]
         Child2.Solution[i] = Parent1.Solution[i]
         if Parent1.Solution[i] == 1
            Child1.CurrentObjectiveValue += Pb.Variables[i]
            Child2.CurrentObjectiveValue += Pb.Variables[i]
         end
      else
         Probk1 = rand(1:2)
         if Probk1 == 1
            Probk2 = 2
            pk = Parent1.CurrentObjectiveValue / ( Parent1.CurrentObjectiveValue + Parent2.CurrentObjectiveValue)
            if rand > pk
               Child1.Solution[i] = Parent1.Solution[i]
               if Parent1.Solution[i] == 1
                  Child1.CurrentObjectiveValue += Pb.Variables[i]
               end
            else
               Child1.Solution[i] = Parent2.Solution[i]
               if Parent2.Solution[i] == 1
                  Child1.CurrentObjectiveValue += Pb.Variables[i]
               end
            end
         else
            Probk2 = 1
            pk = Parent2.CurrentObjectiveValue / ( Parent1.CurrentObjectiveValue + Parent2.CurrentObjectiveValue)
            if rand > pk
               Child1.Solution[i] = Parent2.Solution[i]
               if Parent2.Solution[i] == 1
                  Child1.CurrentObjectiveValue += Pb.Variables[i]
               end
            else
               Child1.Solution[i] = Parent1.Solution[i]
               if Parent1.Solution[i] == 1
                  Child1.CurrentObjectiveValue += Pb.Variables[i]
               end
            end
         end
      end
   end
   println("###################### Crossover Info ######################")
   println("# Parent 1 : ",Parent1.CurrentObjectiveValue," | Parent 2 :",Parent2.CurrentObjectiveValue)
   println("# Child 1 : ",Child1.CurrentObjectiveValue,"# Child 2 : ",Child2.CurrentObjectiveValue)
   println("############################################################")
   return Child1,Child2
end
function InsertAndReplace(Population::Vector{Genome},Indi::Genome)
   index = searchsortedfirst(Population,Indi,by=x->x.CurrentObjectiveValue,rev=true)
   insert!(Population,index,Indi)
   PopSize = length(Population)
   deleteat!(Population,PopSize)
end
