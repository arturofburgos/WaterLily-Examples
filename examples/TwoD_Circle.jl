using WaterLily
include("../src/TwoD_plots.jl")

function circle(n,m;Re=250,U=1)
    radius, center = m/8, m/2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body)
end

# Initialize and run
sim = circle(3*2^6,2^7)
sim_gif!(sim,duration=10,clims=(-5,5),plotbody=true)
