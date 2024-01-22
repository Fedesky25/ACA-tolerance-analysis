using GLMakie;
using ColorSchemes;
import StatsBase: Histogram, fit, normalize;

include("./src/Common.jl");
using .Common;

include("./src/PlotUtils.jl");
using .PlotUtils;

function histbar!(axis, x, y, values; color)
    h = fit(Histogram, values);
    h = normalize(h, mode=:probability);
    edges = h.edges[1]; # it's a tuple, idk why
    for i in 1:length(h.weights)
        bar!(axis, x, y, edges[i], edges[i+1], width=h.weights[i], color=color((edges[i]+edges[i+1])*0.5));
    end
end

struct TolGroup
    rows::UInt;
    cols::UInt;
    tols1::Vector{Float32};
    tols2::Vector{Float32};
    count::UInt;
end

infos = parseBlocksInfos("./input/2500_1e-4/Ts-3.txt");
sel = filter(i -> i.cols > 15 && i.rows < 15, infos);
data = PartitionData(sel);

mat = map(d -> d === nothing ? nothing : TolGroup(
    d[1].rows, d[1].cols,
    map(i -> i.diff1 == 0.0 ? -20 : 0.5*log10(i.diff1/i.exact), d),
    map(i -> i.diff2 == 0.0 ? -20 : 0.5*log10(i.diff2/i.exact), d),
    length(d)
), group(data));

min_tol = typemax(Float32);
max_tol = typemin(Float32);
min_num = typemax(UInt);
max_num = typemin(UInt);
for d in mat
    if d === nothing continue end
    (min1,max1) = extrema(d.tols1);
    (min2,max2) = extrema(d.tols2);
    min_tol = min(min_tol, min1, min2);
    max_tol = max(max_tol, max1, max2);
    (minc,maxc) = extrema(d.count);
    if minc < min_num; min_num = minc; end
    if maxc > max_num; max_num = maxc; end
end

clrmap = cgrad(:thermal, rev=true);
getclr = t -> get(clrmap, (t-min_tol)/(max_tol-min_tol));
cmap2 = cgrad(:matter, rev=true);

aspect = begin (rows, cols) = Common.size(data); rows > cols ? (1.0, cols/rows, 1.0) : (rows/cols, 1.0, 1.0) end;
fig = Figure(; resolution=(1600,500));
ax1 = createAxis3D(fig[1,1], L"\log_{10}(\varepsilon)", title="Before recompression", aspect=aspect);
ax2 = createAxis3D(fig[1,2], L"\log_{10}(\varepsilon)", title="After recompression", aspect=aspect);
ax3 = createAxis3D(fig[1,4], "Frequency", title="Size count", aspect=aspect);
Δ = 0.1*(max_tol-min_tol);
zlims!(ax1, min_tol-Δ, max_tol+Δ);
zlims!(ax2, min_tol-Δ, max_tol+Δ);
Colorbar(fig[1,3], colormap=clrmap, limits=(min_tol, max_tol));
Colorbar(fig[1,5], colormap=cmap2, limits=(min_num, max_num));
link_cameras(ax1, ax2, ax3);

for d in mat
    if d === nothing continue end
    histbar!(ax1, d.rows, d.cols, d.tols1; color=getclr);
    histbar!(ax2, d.rows, d.cols, d.tols2; color=getclr);
    bar!(ax3, d.rows, d.cols, 0, d.count; color=get(cmap2, (d.count-min_num)/(max_num-min_num)));
end

fig