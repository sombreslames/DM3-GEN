function GraspConstruction(Pb::Problem,Indi::Individual, Alpha::Float64)
   Available::Int32            = Pb.NBvariables
   IndiTemp::Individual             = deepcopy(Indi)
   RandomCandidateList         = Vector{Int32}(1)
   Utility = deepcopy(Pb.Utility)
   AboveTheLimit::Int32        = 0
   while Available >= 1
      AboveTheLimit        = 0
      RandomCandidateList  = empty!(RandomCandidateList)
      LimitSelect = (minimum(Utility[2,:]) + (Alpha * (maximum(Utility[2,:])-minimum(Utility[2,:]))))
      for i = 1:1:Available
         if Utility[2,i] >= LimitSelect
            RandomCandidateList = push!(RandomCandidateList,Utility[1,i])
            AboveTheLimit       += 1
         else
            break
         end
      end
      if AboveTheLimit == 0
         break
      else
         RandomPickedCandidate    = rand(RandomCandidateList)
         answer,Indi              = SetToOne(Pb,Indi,convert(Int32,RandomPickedCandidate))
			if answer == false
				println("Failed")
			end
      end
      Available,Utility    = UpdateUtility(Pb,Indi)
   end
   return Indi
end
function ComputeUtilityPlus(Pb::Problem)
   BonusObj::Int32   = 0
   UtilitiesIndex    = Float64[]
   UtilitiesValues   = Float64[]
   VariablesCoef     = Int64[]
   Constraint_Matrix = Array{Bool}(Pb.NBconstraints,0)
   Inc::Int32        = 1
   for i = 1:1:Pb.NBvariables
      nb::Int32            = sum(Pb.LeftMembers_Constraints[:,i])
      if nb != 0
         UtilitiesIndex    = push!(UtilitiesIndex,Inc)
         UtilitiesValues   = push!(UtilitiesValues,(Pb.Variables[i]/nb))
         Constraint_Matrix = hcat(Constraint_Matrix,Pb.LeftMembers_Constraints[:,i])
         VariablesCoef     = push!(VariablesCoef,Pb.Variables[i])
         Inc += 1
      else
         Pb.NBvariables          -= 1
         BonusObj += Pb.Variables[i]
      end
   end
   Utilities         = Matrix(0,Inc-1)
   Utilities         = vcat(Utilities,UtilitiesIndex')
   Utilities         = vcat(Utilities,UtilitiesValues')
   Pb.Variables      = deepcopy(VariablesCoef)
   Pb.LeftMembers_Constraints = deepcopy(Constraint_Matrix)
   Pb.Utility        = deepcopy(Utilities)
   Pb.Utility        = sortcols(Pb.Utility, rev=true, by = x -> x[2])
   return BonusObj,Inc-1,Pb
end

function UpdateUtility(Pb::Problem,Unit::Individual)
   UtilitiesIndex    = Float64[]
   UtilitiesValues   = Float64[]
   Inc::Int32        = 1
   size::Int32       = length(Pb.Utility[1,:])
   for i = 1:1:size
		index = convert(Int,Pb.Utility[1,i])
      if Unit.Freedom[index] == 0 && Unit.Solution[index] == 0
         UtilitiesIndex    = push!(UtilitiesIndex,Pb.Utility[1,i])
         UtilitiesValues   = push!(UtilitiesValues,Pb.Utility[2,i])

         Inc += 1
      end
   end
   Utilities         = Matrix(0,Inc-1)
   Utilities         = vcat(Utilities,UtilitiesIndex')
   Utilities         = vcat(Utilities,UtilitiesValues')
   return Inc-1,Utilities
end

function SetToZero(Pb::Problem,Unit::Individual, x::Int32)
   if Unit.Solution[x] == 1
      for j in 1:1:Pb.NBconstraints
            if Pb.LeftMembers_Constraints[j,x] == 1
               Unit.LastRightMemberValue_Constraint[j] = 0
               for i in 1:1:Pb.NBvariables
                  if Pb.LeftMembers_Constraints[j,i] == 1
                     Unit.Freedom[i]+=1
                  end
               end
            end
      end
   else
      return false,Unit
   end
   Unit.CurrentVarUsed          = deleteat!(Unit.CurrentVarUsed,findin(Unit.CurrentVarUsed,x))
   Unit.Solution[x]             = 0
   Unit.CurrentObjectiveValue   -=Pb.Variables[x]
   return true,Unit
end

function SetToOne(Pb::Problem,Unit::Individual, x::Int32)
   if Unit.Freedom[x] == 0 && Unit.Solution[x] == 0
      for j in 1:1:Pb.NBconstraints
            if Pb.LeftMembers_Constraints[j,x] == 1
               if Unit.LastRightMemberValue_Constraint[j] ==  0
                  Unit.LastRightMemberValue_Constraint[j] = 1
                  for i in 1:1:Pb.NBvariables
                     if Pb.LeftMembers_Constraints[j,i] == 1
                        Unit.Freedom[i]-=1
                     end
                  end
               else
                  return false,Unit
               end
            end
      end
   else
      return false,Unit
   end
   Unit.CurrentVarUsed       = push!(Unit.CurrentVarUsed,x)
   Unit.Solution[x]  = 1
   Unit.CurrentObjectiveValue+=Pb.Variables[x]
   return true,Unit
end

function UpdateReactiveGrasp(AlphaProba::Vector{Float64},Average::Vector{Float64},Worst::Float64,Max::Float64)
   NewValue = Vector{Float64}(length(AlphaProba))
   for i in eachindex(AlphaProba)
      NewValue[i] = ( Average[i] - Worst ) / ( Max - Worst )
   end
   SumOfNew = sum(NewValue)
   for i in eachindex(AlphaProba)
      AlphaProba[i] = NewValue[i] / SumOfNew
   end
   return AlphaProba
end

function ReactiveGrasp(Grasp1::Grasp)
   Proba       = rand()
   Val  = 0
   for i in eachindex(Grasp1.Probability)
      Val += Grasp1.Probability[i]
      if Proba <= Val
         return i, Grasp1.Values[i]
      end
   end
   return length(Grasp1.Probability),Grasp1.Values[length(Grasp1.Values)]
end
#=
#Rules 1 : You will now call GRASP the cycle of foundation
#Rules 2 : No computer sciences student can control others people feelings
#Rules 3 : Waiting for climate change, after all URSS invented it.
function SimulatedAnnealing(CS::CurrentSolution,InitTemperature::Float64,CoolingCoef::Float64,StepSize::Int,MinTemp::Float64)
   CSTemp      = deepcopy(CS)
   CSBest      = deepcopy(CS)
   Temperature = InitTemperature
   ClimateChange         = true
	Historic 	= Int[]
   nbRun       = 0
   while ClimateChange
      for i in 1:1:StepSize
         LocalCS        = AddOrElseDrop(CSTemp)
         #LocalCS 			= GetRandomNeighbour(CSTemp)
         DeltaObj       = LocalCS.CurrentObjectiveValue - CSTemp.CurrentObjectiveValue
			ValueOf			= exp(DeltaObj/Temperature)
			RandValue 		= rand()
         if DeltaObj > 0 || ValueOf > RandValue
				#println("Solution accepted : f(x) ",CSTemp.CurrentObjectiveValue, " --> "," f'(x) : ",CSTemp.CurrentObjectiveValue)
            CSTemp      = deepcopy(LocalCS)
				Historic		= push!(Historic,CSTemp.CurrentObjectiveValue)
            if CSTemp.CurrentObjectiveValue > CSBest.CurrentObjectiveValue
               CSBest   =  deepcopy(CSTemp)
					println("Improved ! We got : ",CSBest.CurrentObjectiveValue)
            end
         end
      end
      nbRun += 1
      #println("Nb of run : ",StepSize * nbRun)
      Temperature *= CoolingCoef
      if Temperature < MinTemp
         ClimateChange = false
      end
   end
   return Historic,CSBest
end

function AddOrElseDrop(CS::CurrentSolution)
   nb,FreeVar     = UpdateUtility(CS)
   CSTemp         = deepcopy(CS)
   if nb > 0
      RandomVar       = convert(Int,FreeVar[1,rand(1:end)])
      answer,CSTemp   =  SetToOne(CSTemp,convert(Int,RandomVar))
      if answer
         return CSTemp
      else
         println("Duck Duck Duck")
      end
   else
      RandomlyPickedUsedVar = rand(CS.CurrentVarUsed)
      answerz,CSTemp        = SetToZero(CSTemp,RandomlyPickedUsedVar)
      if answerz
         return CSTemp
      else
         println("Damn damn damn")
      end
   end
   return nothing
end
#Un petit N puissance 4 au calme coder en 5 min because no need to opti bro
=#
