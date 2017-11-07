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
   return Problem(NBvar, NBcons, Coef, LeftMembers_Constraints,Array{Float64}(2,NBvar),0)
end
function plotRunGrasp(iname,zinit, zls, zbest)
    figure("Examen d'un run",figsize=(6,6)) # Create a new figure
    title("GRASP-SPP | zConst/zLS/zBest | "iname)
    xlabel("Itérations")
    ylabel("valeurs de z(x)")
    ylim(0, maximum(zbest)+2)

    nPoint = length(zinit)
    x=collect(1:nPoint)
    xticks([1,convert(Int64,ceil(nPoint/4)),convert(Int64,ceil(nPoint/2)), convert(Int64,ceil(nPoint/4*3)),nPoint])
    plot(x,zbest, linewidth=2.0, color="green", label="meilleures solutions")
    plot(x,zls,ls="",marker="^",ms=2,color="green",label="toutes solutions améliorées")
    plot(x,zinit,ls="",marker=".",ms=2,color="red",label="toutes solutions construites")
    vlines(x, zinit, zls, linewidth=0.5)
    legend(loc=4, fontsize ="small")
end
function plotAnalyseGrasp(iname, x, zmoy, zmin, zmax)
    figure("bilan tous runs",figsize=(6,6)) # Create a new figure
    title("GRASP-SPP | zMin/zMoy/zMax | "iname)
    xlabel("Itérations (pour nbRunGrasp)")
    ylabel("valeurs de z(x)")
    ylim(0, zmax[end]+2)

    nPoint = length(x)
    intervalle = [reshape(zmoy,(1,nPoint)) - reshape(zmin,(1,nPoint)) ; reshape(zmax,(1, nPoint))-reshape(zmoy,(1,nPoint))]
    xticks(x)
    errorbar(x,zmoy,intervalle,lw=1, color="black", label="zMin zMax")
    plot(x,zmoy,linestyle="-", marker="o", ms=4, color="green", label="zMoy")
    legend(loc=4, fontsize ="small")
end
