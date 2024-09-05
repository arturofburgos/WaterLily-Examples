using WaterLily
using StaticArrays
include("../src/ThreeD_plots.jl")
norm(x::StaticArray) = √(x'*x)

function make_sim(;n=2^6, U=1, Re=1000, mem=Array)
    # define the body
    R = n/2
    function sdf(xyz, t)
        norm(xyz) - R # sphere
    end

    map(xyz, t) = xyz - SA[2n/3,0,0] # place/move the center

    # Return Simulation
    return Simulation((2n,n,n), (U,0,0), R;
        ν=U*R/Re, body=AutoBody(sdf,map), mem)
end

function mirrorto!(a, b)
    # Notice, this function could be wrong

    n = size(b,1)
    p = size(b,2)
    # q = size(b,3)
    a[1:n,reverse(1:p),:].=b #apparently good, in accordance with the video tutorial
    # a[reverse(1:n),reverse(1:p),:].=b #apparently wrong, not in accordance to the video tutorial
    a[:, reverse(1:p), :].=a[:, 1:p, :]

    # Apparently all wrong - Dont forget I am not trying to mirror the body, but the domain that encloses it, to gain more potential fluid cells to be altered by the fluid flow
    # a[reverse(1:n),:,:].=a[:,:,:]
    # a[:,reverse(p+1:2p),:].=a[:,1:p,:]
    # a[:,reverse(p+1:2p),:].=a[:,reverse(1:p),:]
    return a
end

# sim = make_sim();

using CUDA
@assert CUDA.functional()

begin
    sim = make_sim(;mem = CuArray);
    sim_step!(sim, 0.1)
    
    # Create CPU buffer arrays for geometry flow viz 
    a = sim.flow.σ
    d = similar(a,size(inside(a))) |> Array; # one quadrant
    md = similar(d,(2,2,1).*size(d))  # hold mirrored data

    # Set up geometry viz
    geom = geom!(md,d,sim) |> Observable;
    fig, _, _ = GLMakie.mesh(geom, alpha=1, color=:lightblue)

    
    volume!(d) #comment if necessary
    volume!(md) #comment if necessary
    fig
end