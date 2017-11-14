function InitPopulation(Pb::Problem,N::Int32,Grasp::Grasp)
   #=println("####################### Problem Info #######################")
   println("# Var : ",Pb.NBvariables, " | Constraints : ", Pb.NBconstraints)=#
   Pb.Bonus,NB,Pb = ComputeUtilityPlus(Pb)
   #=println("# f(x): ", Pb.Bonus, " | Variable considered : ",NB)
   println("############################################################")=#
   Population  = Array{Genome}(N)
   TotSum      = 0
   for i = 1:1:N
      IndivTemp         = Individual(Pb.Bonus,zeros(Int64,Pb.NBvariables),zeros(Int64,0),zeros(Int64,Pb.NBconstraints),zeros(Int64,Pb.NBvariables))
      IndexAlpha,Alpha  = ReactiveGrasp(Grasp)
      IndivTemp         = GraspConstruction(Pb,IndivTemp,Alpha)
      Population[i]     = Genome(IndivTemp.CurrentObjectiveValue,IndivTemp.Solution)
      TotSum += Population[i].CurrentObjectiveValue
   end
   Population = sort(Population,by=a->a.CurrentObjectiveValue,rev=true)
   #=println("###################### Population Info #####################")
   println("# Moyenne : ",TotSum/N)
   println("# MIN : ", Population[N].CurrentObjectiveValue, " | MAX : ",Population[1].CurrentObjectiveValue)
   println("############################################################")=#

   Pb.MinObj = Population[N].CurrentObjectiveValue
   Pb.SumObj = TotSum
   return Pb,Population
end

function Evolution(Population::Vector{Genome},Pb::Problem,NPop::Int32,Ngen::Int32,stoplimit::Float64)
   it = 0
   for i  = 1:1:Ngen
       for j = 1:2:NPop
         #here repair and mutate
         ProbMut = rand(RdSeed)
         if ProbMut > 0.3
            p1,p2          = BinaryTourmanent(Population,NPop,Pb)
            child1,child2  = CrossoverMethod(Pb,[Population[p1],Population[p2]])

            #PLus divers en pop mais converge moins vite
            #=if child1.Solution == child2.Solution
               if child1.Solution != Population[p1].Solution
                  if child1.Solution != Population[p2].Solution
                     child1 = RepairAndMutationSparse(Pb,child1)
                     child1 = AugmentIndividual(Pb,child1)
                  end
               end

            else
               if child1.Solution != Population[p1].Solution
                  if child1.Solution != Population[p2].Solution
                     child1 = RepairAndMutationSparse(Pb,child1)
                     child1 = AugmentIndividual(Pb,child1)
                  end
               end
               if child2.Solution != Population[p1].Solution
                  if child2.Solution != Population[p2].Solution
                     child2 = RepairAndMutationSparse(Pb,child2)
                     child2 = AugmentIndividual(Pb,child2)
                     InsertAndReplace(Pb,Population,child2)
                  end
               end
            end
            InsertAndReplace(Pb,Population,child1)=#
            #Converge + vite mais moins divers en pop
            if child1.Solution != Population[p1].Solution && child1.Solution != Population[p2].Solution
               child1 = RepairAndMutationSparse(Pb,child1)
               child1 = AugmentIndividual(Pb,child1)
            end
            if child2.Solution != Population[p1].Solution && child2.Solution != Population[p2].Solution
               child2 = RepairAndMutationSparse(Pb,child2)
               child2 = AugmentIndividual(Pb,child2)
            end
            InsertAndReplace(Pb,Population,child2)
            InsertAndReplace(Pb,Population,child1)

         end
         it += 2
      end
      #=println("###################### Generation Info #####################")
      println("# Generation : ",i)
      println("# Moyenne : ",round(Pb.SumObj/NPop,2))
      println("# MIN : ", Pb.MinObj, " | MAX : ",Population[1].CurrentObjectiveValue)
      println("# Time spend : ",time, "s | ",round((it/(Ngen*NPop))*100,2),"% done")
      println("############################################################")=#
      if ((Population[1].CurrentObjectiveValue - Population[NPop].CurrentObjectiveValue) / (Pb.SumObj/NPop)) < stoplimit
         break
      end
   end
   return Population
end
function RepairAndMutationSparse(Pb::Problem,Indi::Genome)

   #=println("###################### Repair Info ######################")
   println("# Before repair : ",Indi.CurrentObjectiveValue)=#
   IndexRow,IndexColumn,Value = findnz(Pb.LeftMembers_Constraints)
   size = length(IndexRow)
   SaturatedCtrt = 0
   for i in 1:1:size
      if Indi.Solution[IndexColumn[i]] == 1

         Pos = findfirst(IndexRow,IndexRow[i])
         SaturatedCtrt = 0
         while Pos != 0
            if IndexColumn[Pos] != IndexColumn[i]
               if Indi.Solution[IndexColumn[Pos]] == 1
                  Indi.Solution[IndexColumn[Pos]] = 0
                  Indi.CurrentObjectiveValue -= Pb.Variables[IndexColumn[Pos]]
               end
            end
            Pos = findnext(IndexRow,IndexRow[i],Pos+1)
         end
      end
   end

   #=println("# After repair : ",Indi.CurrentObjectiveValue)
   println("############################################################")=#

   return Indi
end

function AugmentIndividual(Pb::Problem,Indi::Genome)
   IndexRow,IndexColumn,Value = findnz(Pb.LeftMembers_Constraints)
   #=println("###################### Mutation Info ######################")
   println("# Before mutation : ",Indi.CurrentObjectiveValue)=#
   #Freedom = zeros(Int32,Pb.NBvariables)
   size = length(IndexColumn)
   for i = 1:1:Pb.NBvariables
      if Indi.Solution[i] == 0 #&& Freedom[i] == 0
         Ok = false
         PosV = findfirst(IndexColumn,i)
         while Ok  == false && PosV != 0 && IndexColumn[PosV] == i
            PosC = findfirst(IndexRow,IndexRow[PosV])
            while PosC != 0
               if Indi.Solution[IndexColumn[PosC]] == 1
                  Ok = true
                  break
               end
               #Freedom[IndexColumn[PosC]] -= 1
               PosC = findnext(IndexRow,IndexRow[PosV],PosC+1)
            end
            PosV = findnext(IndexColumn,i,PosV+1)
         end
         if Ok == false
            Indi.Solution[i] = 1
            Indi.CurrentObjectiveValue += Pb.Variables[i]
         end
      end
   end
   #=println("# After mutation : ",Indi.CurrentObjectiveValue)
   println("############################################################")=#
   return Indi
end

function RouletteSelection(Population::Vector{Genome},N::Int32,Pb::Problem)
   Value::Float64 = 0.0
   p1::Int32 = 0
   p2::Int32 = 0
   Rand1::Float64 = rand(RdSeed)
   Rand2::Float64 = rand(RdSeed)
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
   return rand(RdSeed,1:N),rand(RdSeed,1:N)
end
function BinaryTourmanent(Population::Vector{Genome},N::Int32,Pb::Problem)
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
#  FUsion operator yeah it rock
function CrossoverMethod(Pb::Problem,Parents::Vector{Genome})
   Child1 = Genome(Pb.Bonus,zeros(Int32,Pb.NBvariables))
   Child2 = Genome(Pb.Bonus,zeros(Int32,Pb.NBvariables))
   #Soyons malin utilisons les probabilites
   for i in 1:1:Pb.NBvariables
      if Parents[1].Solution[i]  == Parents[2].Solution[i]
         Child1.Solution[i] = Parents[1].Solution[i]
         Child2.Solution[i] = Parents[1].Solution[i]
      else
         Probk1 = rand(RdSeed,1:2)
         Probk2 = 0
         if Probk1 == 1
            Probk2 = 2
         else
            Probk2 = 1
         end
         pk = Parents[Probk1].CurrentObjectiveValue / ( Parents[1].CurrentObjectiveValue + Parents[2].CurrentObjectiveValue)
         pkPrim = Parents[Probk2].CurrentObjectiveValue / ( Parents[1].CurrentObjectiveValue + Parents[2].CurrentObjectiveValue)
         if rand(RdSeed) <= pk
            Child1.Solution[i] = Parents[Probk1].Solution[i]
         else
            Child1.Solution[i] = Parents[Probk2].Solution[i]
         end
         if rand(RdSeed) <= pkPrim
            Child2.Solution[i] = Parents[Probk2].Solution[i]
         else
            Child2.Solution[i] = Parents[Probk1].Solution[i]
         end
      end
      if Child1.Solution[i] == 1
         Child1.CurrentObjectiveValue += Pb.Variables[i]
      end
      if Child2.Solution[i] == 1
         Child2.CurrentObjectiveValue += Pb.Variables[i]
      end
   end
   #println("###################### Crossover Info ######################")
   #println("# Parent 1 : ",Parents[1].CurrentObjectiveValue," | Parent 2 :",Parents[2].CurrentObjectiveValue)
   #println("# Child 1 : ",Child1.CurrentObjectiveValue," | Child 2 : ",Child2.CurrentObjectiveValue)
   #println("############################################################")
   return Child1,Child2
end

function InsertAndReplace(Pb::Problem,Population::Vector{Genome},Indi::Genome)
   index = searchsortedfirst(Population,Indi,by=x->x.CurrentObjectiveValue,rev=true)
   insert!(Population,index,Indi)
   Pb.SumObj += Indi.CurrentObjectiveValue
   PopSize = length(Population)
   Pb.SumObj -= Population[PopSize].CurrentObjectiveValue
   deleteat!(Population,PopSize)
   Pb.MinObj = Population[PopSize-1].CurrentObjectiveValue
end
