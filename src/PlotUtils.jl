module PlotUtils

import GLMakie: mesh!, Axis3, on;
import StatsBase: Histogram, fit, normalize;

export createAxis3D, bar!, link_cameras, link_cameras2, histbar!;

const rect_faces = [ 
    1 2 3; 3 4 1; # bottom
    1 6 2; 1 6 5; # side y-
    2 7 3; 2 7 6; # side x+
    3 8 4; 3 8 7; # side y+
    4 5 1; 4 5 8; # side x-
    5 6 7; 7 8 5; # top
];

function bar!(axis, x, y, z1, z2; width=1, color)
    h = 0.5*width;
    vertices = [
        x-h  y-h  z1;
        x+h  y-h  z1;
        x+h  y+h  z1;
        x-h  y+h  z1;
        x-h  y-h  z2;
        x+h  y-h  z2;
        x+h  y+h  z2;
        x-h  y+h  z2;
    ];
    mesh!(axis, vertices, rect_faces, shading=false, color=color);
end

function createAxis3D(pos; title = "", aspect=(1,1,1));
    Axis3(pos, title=title, xlabel="rows", ylabel="columns", zlabel="", aspect=aspect, azimuth=-0.25*π, elevation=π/7, perspectiveness=0.3)
end

function link_cameras(axes::Axis3...; step=0.01)
    len = length(axes);
    for i in 1:len
        on(axes[i].azimuth) do x
            for j in 1:(i-1); if abs(x - axes[j].azimuth[]) > step; axes[j].azimuth[] = x; end end
            for j in (i+1):len; if abs(x - axes[j].azimuth[]) > step; axes[j].azimuth[] = x; end end
        end
        on(axes[i].elevation) do x
            for j in 1:(i-1); if abs(x - axes[j].elevation[]) > step; axes[j].elevation[] = x; end end
            for j in (i+1):len; if abs(x - axes[j].elevation[]) > step; axes[j].elevation[] = x; end end
        end
    end
end

function link_cameras2(axes::Axis3...; step=0.01)
    len = length(axes);
    az = axes[1].azimuth[];
    el = axes[1].elevation[];
    for i in 1:len
        on(axes[i].azimuth) do x
            if abs(x - az) > step; 
                az = x;
                for j in 1:(i-1); axes[j].azimuth[] = x; end
                for j in (i+1):len; axes[j].azimuth[] = x; end
            end
        end
        on(axes[i].elevation) do x
            if abs(x-el) > step
                el = x;
                for j in 1:(i-1); axes[j].elevation[] = x; end
                for j in (i+1):len; axes[j].elevation[] = x; end
            end
        end
    end
end

function histbar!(axis, x, y, values, color; onone=1)
    h = fit(Histogram, values);
    h = normalize(h, mode=:probability);
    if length(h.weights) == 1
        # all the values are the same
        bar!(axis, x, y, values[1]-onone*0.5, values[1]+onone*0.5, color=color(values[1]));
    else 
        edges = h.edges[1]; # it's a tuple, idk why
        for (i, w) in enumerate(h.weights)
            bar!(axis, x, y, edges[i], edges[i+1], width=sqrt(w), color=color((edges[i]+edges[i+1])*0.5));
        end
    end
end

# struct InverseRange
#     min::Float64;
#     max::Float64;
#     bins::Int64;
#     α::Float64;
#     InverseRange(min::Float64, max::Float64, bins::Int64) = new(min, max, bins, bins/(max-min));
# end

# function getbin(ir::InverseRange{T}, value::{<:T})::Int64 where T;
#     if value == ir.max;
#         ir.bins
#     else
#         floor(Int64, (value-ir.min)*ir.α)
#     end
# end

# struct BinsIter
#     ir::InverseRange
#     bins::Vector{UInt64}
#     BinsIter(ir::InverseRange, values) = begin
#         bins = zeros(UInt64, ir.bins);
#         for v ∈ values;
#             bins[getbin(ir, v)] += 1;
#         end
#         return new(ir, bins);
#     end
# end

# Base.step(r::InverseRange) = (r.max-r.min)/r.bins;
# Base.iterate(bi::BinsIter) = ()

# function histbar2!(axis, x, y, values, color; onone=1)
#     h = fit(Histogram, values);
#     h = normalize(h, mode=:probability);
#     if length(h.weights) == 1
#         # all the values are the same
#         bar!(axis, x, y, values[1]-onone*0.5, values[1]+onone*0.5, color=color(values[1]));
#     else 
#         edges = h.edges[1]; # it's a tuple, idk why
#         for (i, w) in enumerate(h.weights)
#             bar!(axis, x, y, edges[i], edges[i+1], width=sqrt(w), color=color((edges[i]+edges[i+1])*0.5));
#         end
#     end
# end

end