function ReadFile(FileName::String)
   workingfile    = open(FileName)
   NBcons,NBvar   = parse.(split(readline(workingfile)))
   Coef           = parse.(split(readline(workingfile)))
   LeftMembers_Constraints    = spzeros(NBcons,NBvar)
   RightMembers_Constraints   = Vector(NBcons)
   for i = 1:1:NBcons
         readline(workingfile)
         RightMembers_Constraints[i]=1
         for val in split(readline(workingfile))
            LeftMembers_Constraints[i, parse(val)]=1
         end
   end
   close(workingfile)
   IndexRow,IndexColumn,Value = findnz(LeftMembers_Constraints)
   OrderedByCons = vcat(IndexColumn',IndexRow')
   OrderedByCons = sortcols(OrderedByCons,by=a->a[2])
   return Problem(NBvar, NBcons, Coef, IndexRow,IndexColumn,OrderedByCons[2,:],OrderedByCons[1,:],LeftMembers_Constraints,Array{Float64}(2,NBvar),0,0,0)
end
function gen_colors(n)
  cs = distinguishable_colors(n,
      [colorant"#FE4365", colorant"#eca25c"], # seed colors
      lchoices=Float64[58, 45, 72.5, 90],     # lightness choices
      transform=c -> deuteranopic(c, 0.1),    # color transform
      cchoices=Float64[20,40],                # chroma choices
      hchoices=[75,51,35,120,180,210,270,310] # hue choices
  )

  convert(Vector{Color}, cs)
end
function PlotGeneticAlgorithm(Population::Array{Vector{Genome}},N::Int32,Ngen::Int32,filename::String)
   nbpoint = N*Ngen
   x = collect(1:N)
   y = Array{Int32}(N*Ngen,3)
   Plot = Vector{Gadfly.Plot}(Ngen)
   for i = 1:1:N
      for j = 1:1:Ngen
         y[(i-1)*Ngen + j,1] = Population[j][i].CurrentObjectiveValue
         y[(i-1)*Ngen + j,2] = j
         y[(i-1)*Ngen + j,3]  = i
         #=if i == 1
            rename!(y,:j,string("individu ",j))
         end=#
      end

   end
   y = convert(DataFrame, y)
   names!(y, [Symbol("f(x)"),Symbol("Generation"),Symbol("Numero")])
   set_default_plot_size(35cm, 25cm)
   plot(y, x="Numero",y="f(x)",
         color="Generation",
         Scale.color_discrete(),
         Guide.title(string(Ngen," Generation pour le probleme ",filename)),
         Theme(default_color=color("red"))
   )
   #plot(y, x="Generation 6",y="Generation 7" , Geom.point)
   #plot(y,x="SepalLength", y="SepalWidth", Geom.point)
end
