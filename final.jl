using GLMakie;
using ColorSchemes;
import StatsBase: Histogram, fit, normalize;

include("./src/PlotUtils.jl");
using .PlotUtils;

# 1in = 72pt
# good dpi = 300 px/in = 4.1667 px/pt
# width = 7.16in = 515.52pt = 2148px
# 515.5 is too small => take directly 2148px but multiply fontsize by 4.1667 (9*4.1667 ≈ 37)
update_theme!(fontsize=22);
PLOT_WIDTH = 2148;
LABEL_SIZE = 35;

struct MinMax{T}
    min::T;
    max::T;
end

Base.show(io::IO, mm::MinMax) = print(io, "(min $(mm.min), max $(mm.max))");
Base.size(mm::MinMax) = mm.max - mm.min;

struct NormData 
    exact::Float64
    diff1::Float64
    diff2::Float64
end

struct BlockData 
    l2::NormData
    op::NormData
    sigma::MinMax{Float64}
    k1::UInt
    k2::UInt
    rows::UInt
    cols::UInt
end

struct SizedBlocksData
    rows::UInt
    cols::UInt
    count::UInt
    surface::UInt
    ϵF1::Vector{Float64} # frobenius errors before
    ϵF2::Vector{Float64} # frobenius error after
    ϵσ1::Vector{Float64}
    ϵσ2::Vector{Float64}
    κ1::Vector{Float64}
    κ2::Vector{Float64}
    ρ::Vector{Float64}
end

struct BinEdges 
    min::Float64;
    max::Float64;
    points::Int64;
end




function getBlocks(filename)
    f = open(filename, "r")
    infos = map(readlines(f)) do line
        c = split(line, " ");
        rows = parse(UInt, c[11]);
        cols = parse(UInt, c[12]);
        if cols < rows; (rows,cols) = (cols,rows); end
        return BlockData(
            NormData(
                parse(Float64, c[1]),
                parse(Float64, c[2]),
                parse(Float64, c[3])
            ),
            NormData(
                parse(Float64, c[4]),
                parse(Float64, c[5]),
                parse(Float64, c[6])
            ),
            MinMax(
                parse(Float64, c[7]),
                parse(Float64, c[8]),
            ),
            parse(UInt, c[9]),
            parse(UInt, c[10]),
            rows, cols
        );
    end;
    close(f);
    return infos;
end

function group(blocks::Vector{BlockData})
    g = Dict{Tuple{UInt, UInt}, Vector{BlockData}}();
    for b in blocks
        key = (b.rows,b.cols);
        if haskey(g, key)
            push!(g[key], b);
        else
            g[key] = [b];
        end
    end
    return g;
end

function dimension(gb)
    min_cols = min_rows = typemax(UInt);
    max_cols = max_rows = typemin(UInt);
    for (pos, _) in gb
        if pos[1] < min_rows; min_rows = pos[1]; 
        elseif pos[1] > max_rows; max_rows = pos[1];
        end
        if pos[2] < min_cols; min_cols = pos[2]; 
        elseif pos[2] > max_cols; max_cols = pos[2];
        end
    end
    return (MinMax(min_rows,max_rows), MinMax(min_cols,max_cols));
end

function dimension(dv::Vector{SizedBlocksData})
    min_cols = min_rows = typemax(UInt);
    max_cols = max_rows = typemin(UInt);
    for d ∈ dv
        if d.rows < min_rows; min_rows = d.rows; 
        elseif d.rows > max_rows; max_rows = d.rows;
        end
        if d.cols < min_cols; min_cols = d.cols; 
        elseif d.cols > max_cols; max_cols = d.cols;
        end
    end
    return (MinMax(min_rows,max_rows), MinMax(min_cols,max_cols));
end

function getlims(min, max, margin) 
    Δ = (max-min)*margin;
    return (min-Δ, max+Δ);
end

function createAxesCouple(fig, col::Int, label, aspect, min::Float64, max::Float64, color; rev=true, clrlims::Union{Nothing,Tuple{Number,Number}}=nothing)
    clr = cgrad(color, rev=rev);
    if clrlims === nothing; clrlims = (min,max); end
    getclr = t -> get(clr, (t-clrlims[1])/(clrlims[2]-clrlims[1]));
    ax1 = createAxis3D(fig[1, col], aspect=aspect);
    ax2 = createAxis3D(fig[2, col], aspect=aspect);
    Colorbar(fig[3, col], colormap=clr, limits=clrlims, vertical=false, flipaxis=false, label=label, labelsize=LABEL_SIZE);
    Δ = (max-min)*0.05;
    zlims!(ax1, min-Δ, max+Δ);
    zlims!(ax2, min-Δ, max+Δ);
    return (ax1, ax2, getclr);
end

function plotBlocks(blocks::Vector{BlockData}; filter::Function = (_) -> true)
    grouped = group(blocks);
    data = Vector{SizedBlocksData}(undef,0);
    ϵF_min = ϵσ_min = κ_min = ρ_min = typemax(Float64);
    ϵF_max = ϵσ_max = κ_max = ρ_max = typemin(Float64);
    surface_min = typemax(UInt);
    surface_max = typemin(UInt);

    for (pos, blocks) in grouped
        len = length(blocks);
        ϵF1 = Vector{Float64}(undef, len);
        ϵF2 = Vector{Float64}(undef, len);
        ϵσ1 = Vector{Float64}(undef, len);
        ϵσ2 = Vector{Float64}(undef, len);
        κ1 = Vector{Float64}(undef, len);
        κ2 = Vector{Float64}(undef, len);
        α = (pos[1]+pos[2])/(pos[1]*pos[2]);
        ρ = Vector{Float64}(undef,len);
        count = 0;

        for b in blocks
            # frobenius relative error
            ϵF = b.l2.diff1 == 0.0 ? -20 : 0.5*log10(b.l2.diff1/b.l2.exact);
            if !filter(ϵF); continue; end;
            count += 1;
            ϵF1[count] = ϵF;
            ϵF2[count] = b.l2.diff2 == 0.0 ? -20 : 0.5*log10(b.l2.diff2/b.l2.exact);
            ϵF_min = min(ϵF_min, ϵF1[count], ϵF2[count]);
            ϵF_max = max(ϵF_max, ϵF1[count], ϵF2[count]);
            # spectral relative error
            ϵσ1[count] = b.op.diff1 == 0.0 ? -20 : log10(b.op.diff1/b.op.exact);
            ϵσ2[count] = b.op.diff2 == 0.0 ? -20 : log10(b.op.diff2/b.op.exact);
            ϵσ_min = min(ϵσ_min, ϵσ1[count], ϵσ2[count]);
            ϵσ_max = max(ϵσ_max, ϵσ1[count], ϵσ2[count]);
            # compression ratio
            κ1[count] = b.k1*α;
            κ2[count] = b.k2*α;
            κ_min = min(κ_min, κ1[count], κ2[count]);
            κ_max = max(κ_max, κ1[count], κ2[count]);
            # singular value ratio
            ρ[count] = log10(b.sigma.max/b.sigma.min);
            if ρ[count] > ρ_max; ρ_max = ρ[count] elseif ρ[count] < ρ_min; ρ_min = ρ[count] end
        end

        if count > 0
            surface = count*pos[1]*pos[2];
            if surface > surface_max; surface_max = surface; elseif surface < surface_min; surface_min = surface; end;
            if count != len
                resize!(ϵF1, count);
                resize!(ϵF2, count);
                resize!(ϵσ1, count);
                resize!(ϵσ2, count);
                resize!(κ1, count);
                resize!(κ2, count);
                resize!(ρ, count);
            end
            push!(data,SizedBlocksData(
                pos[1], pos[2], count, surface, 
                ϵF1, ϵF2, ϵσ1, ϵσ2, 
                κ1, κ2, ρ
            ));
        end

    end

    
    blockDims = dimension(data);
    aspect = begin 
        rows = size(blockDims[1]);
        cols = size(blockDims[2]);
        rows > cols ? (1.0, cols/rows, 1.0) : (rows/cols, 1.0, 1.0) 
    end;

    println("Calculation completed -> plotting");

    height = floor(Int, PLOT_WIDTH * 4/7);
    fig = Figure(; resolution=(PLOT_WIDTH,height));
    frobenius = createAxesCouple(fig, 1, L"\log(\varepsilon_{F})", aspect, ϵF_min, ϵF_max, :thermal);
    spectral = createAxesCouple(fig, 2, L"\log(\varepsilon_{\sigma})", aspect, ϵσ_min, ϵσ_max, :seaborn_rocket_gradient);
    compression = createAxesCouple(fig, 3, L"\text{Compression ratio } \kappa", aspect, κ_min, κ_max, :berlin, rev=false, clrlims=(0,2));
    # σ_ratio = createAxesCouple(fig, 4, L"\log(\rho)", aspect, sr_min, sr_max, :linear_ternary_red_0_50_c52_n256);

    clr_ρ = cgrad(:linear_ternary_red_0_50_c52_n256, rev=true);
    ax_ρ = createAxis3D(fig[1,4], aspect=aspect);
    Colorbar(fig[1,5], colormap=clr_ρ, limits=(ρ_min,ρ_max), label=L"\log(\rho)", labelsize=LABEL_SIZE);
    zlims!(ax_ρ, ρ_min, ρ_max);
    get_ρ_clr = t -> get(clr_ρ, (t-ρ_min)/(ρ_max-ρ_min));

    clr_surface = cgrad(:matter, rev=true);
    ax_surface = createAxis3D(fig[2,4], aspect=aspect);
    Colorbar(fig[2,5], colormap=clr_surface, limits=(surface_min,surface_max), label="Surface", labelsize=LABEL_SIZE);
    zlims!(ax_surface, 0, nothing);

    link_cameras2(frobenius[1], frobenius[2], spectral[1], spectral[2], compression[1], compression[2], ax_ρ, ax_surface);
    for ax in [frobenius[1], frobenius[2], spectral[1], spectral[2], compression[1], compression[2], ax_ρ, ax_surface]
        xlims!(ax, blockDims[1].min-1, blockDims[1].max+1);
        ylims!(ax, blockDims[2].min-1, blockDims[2].max+1);
    end

    for sbd in data
        histbar!(frobenius[1], sbd.rows, sbd.cols, sbd.ϵF1, frobenius[3]);
        histbar!(frobenius[2], sbd.rows, sbd.cols, sbd.ϵF2, frobenius[3]);
        histbar!(spectral[1], sbd.rows, sbd.cols, sbd.ϵσ1, spectral[3]);
        histbar!(spectral[2], sbd.rows, sbd.cols, sbd.ϵσ2, spectral[3]);
        histbar!(compression[1], sbd.rows, sbd.cols, sbd.κ1, compression[3], onone=0.05);
        histbar!(compression[2], sbd.rows, sbd.cols, sbd.κ2, compression[3], onone=0.05);
        histbar!(ax_ρ, sbd.rows, sbd.cols, sbd.ρ, get_ρ_clr, onone=0.1);
        bar!(ax_surface, sbd.rows, sbd.cols, 0, sbd.surface, color=get(clr_surface, (sbd.surface-surface_min)/(surface_max-surface_min)));
    end
    fig
end

blocks = getBlocks("input/2500_1e-4_bis/Ts.txt");
sel = filter(b -> b.cols < 50, blocks);
f = plotBlocks(sel, filter=ϵ -> ϵ >= -7)