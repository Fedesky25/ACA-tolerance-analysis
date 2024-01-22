using GLMakie
using ColorSchemes
using StatsBase

rect_faces = [
    # bottom
    1 2 3;
    3 4 1;
    # side y-
    1 6 2;
    1 6 5;
    # side x+
    2 7 3;
    2 7 6;
    # side y+
    3 8 4;
    3 8 7;
    # side x-
    4 5 1;
    4 5 8;
    # top
    5 6 7;
    7 8 5;
];

function bar(x, y, from, to, axis, clr)
    vertices = [
        x-0.5  y-0.5  from;
        x+0.5  y-0.5  from;
        x+0.5  y+0.5  from;
        x-0.5  y+0.5  from;
        x-0.5  y-0.5  to;
        x+0.5  y-0.5  to;
        x+0.5  y+0.5  to;
        x-0.5  y+0.5  to;
    ];
    mesh!(axis, vertices, rect_faces, shading=false, color=clr);
    # lines!(axis, [
    #     x-0.5  y-0.5  to;
    #     x+0.5  y-0.5  to;
    #     x+0.5  y+0.5  to;
    #     x-0.5  y+0.5  to;
    #     x-0.5  y-0.5  to;
    # ], color=:black);
    # lines!(axis, [x-0.5  y-0.5  from; x-0.5  y-0.5  to;]);
    # lines!(axis, [x+0.5  y-0.5  from; x+0.5  y-0.5  to;]);
    # lines!(axis, [x-0.5  y+0.5  from; x-0.5  y+0.5  to;]);
    # lines!(axis, [x+0.5  y+0.5  from; x+0.5  y+0.5  to;]);
end

function minmax_noNaN(arr) 
    min = typemax(eltype(arr));
    max = typemin(eltype(arr));
    for v in arr
        if isnan(v)
            continue
        end
        if v < min
            min = v;
        elseif v > max
            max = v;
        end 
    end
    return min, max
end

function bars3D(data, xshift, yshift, colormap, axis, cmap_pos)
    cmap = colorschemes[colormap];
    min, max = minmax_noNaN(data);
    rows, cols = size(data);
    for i in 1:rows, j in 1:cols
        z = data[i,j];
        if isnan(z)
            continue;
        end
        bar(i+xshift, j+yshift, min, z, axis, get(cmap, (z-min)/(max-min)));
    end
    Colorbar(cmap_pos, colormap=colormap, limits=(min, max));
end

function create_axis3D(zlabel, pos; title = "");
    Axis3(pos, title=title, xlabel="rows", ylabel="columns", zlabel=zlabel, aspect=(1,1,1), elevation=π/6, perspectiveness=0.5)
end

struct BlockInfo
    exact::Float32
    diff1::Float32
    diff2::Float32
    k1::UInt
    k2::UInt
    rows::UInt
    cols::UInt
end

function parseBlocksInfo(filename)
    f = open(filename, "r")
    infos = map(readlines(f)) do line
        c = split(line, " ");
        return BlockInfo(
            parse(Float32, c[1]),
            parse(Float32, c[2]),
            parse(Float32, c[3]),
            parse(UInt, c[4]),
            parse(UInt, c[5]),
            parse(UInt, c[6]),
            parse(UInt, c[7]),
        );
    end;
    close(f);
    return infos;
end

function calcTols(diff, exact) 
    map((d,e) -> e == 0.0 ? NaN : 0.5*log10(d/e), diff, exact)
end

function averageRatios(sum_ratios, count)
    map((s,c) -> c == 0.0 ? NaN32 : convert(Float32, s/c), sum_ratios, count)
end

function showTotalTolerances(infos::Vector{BlockInfo})
    tot = sum(i -> i.exact, infos);
    tol1 = sqrt(sum(i -> i.diff1, infos)/tot);
    tol2 = sqrt(sum(i -> i.diff2, infos)/tot);
    println("Tolerance: $tol1 -> $tol2");
end

function toleranceHistogram(data::Vector{BlockInfo})
    tols1 = map(i -> 0.5*log10(i.diff1/i.exact), data);
    tols2 = map(i -> 0.5*log10(i.diff2/i.exact), data);
    fig = Figure(; resolution=(700,1000));
    ax1 = Axis(fig[1,1], title="Before recompression");
    ax2 = Axis(fig[2,1], title="After recompression")
    linkxaxes!(ax1,ax2);
    hist!(ax1, tols1, normalization = :probability,
        bar_labels = :values, label_formatter = x-> round(x, digits=2), label_size = 15,
        strokewidth = 0.5, strokecolor = (:black, 0.5),
        color=:values, colormap=:deep
    );
    hist!(ax2, tols2, normalization = :probability, 
        bar_labels = :values, label_formatter = x-> round(x, digits=2), label_size = 15,
        strokewidth = 0.5, strokecolor = (:black, 0.5),
        color=:values, colormap=:deep
    );
    display(fig);
    return fig;
end

infos = parseBlocksInfo("input/2500_1e-4/Ts-2.txt");

# selection = filter(i -> (i.rows<25 && i.cols>75) || (i.cols<25 && i.rows>75), infos);
selection = filter(i -> (i.rows>75 && i.cols>75), infos);
toleranceHistogram(selection);
showTotalTolerances(selection);
println("range: $(minimum(i -> i.rows, selection)) - $(maximum(i -> i.rows, selection))");


# print overall tolerance
# showTotalTolerances(infos);

# find cols, rows bounds
col_min = typemax(UInt);
row_min = typemax(UInt);
col_max = 0;
row_max = 0;
for bi in infos
    if bi.rows < row_min
        global row_min = bi.rows
    elseif bi.rows > row_max
        global row_max = bi.rows
    end
    if bi.cols < col_min
        global col_min = bi.cols
    elseif bi.cols > col_max
        global col_max = bi.cols
    end
end

# computes average tolerances and compression ratio for each row and column
rows = row_max-row_min+1;
cols = col_max-col_min+1;
exact = zeros(Float32, rows, cols);
sum_diff1 = zeros(Float32, rows, cols);
sum_diff2 = zeros(Float32, rows, cols);
sum_ratio1 = zeros(Float32, rows, cols);
sum_ratio2 = zeros(Float32, rows, cols);
count = zeros(UInt, rows, cols);
for i in infos 
    exact[i.rows-row_min+1, i.cols-col_min+1] += i.exact;
    sum_diff1[i.rows-row_min+1, i.cols-col_min+1] += i.diff1;
    sum_diff2[i.rows-row_min+1, i.cols-col_min+1] += i.diff2;
    sum_ratio1[i.rows-row_min+1, i.cols-col_min+1] += i.k1*(i.rows+i.cols)/(i.rows*i.cols);
    sum_ratio2[i.rows-row_min+1, i.cols-col_min+1] += i.k2*(i.rows+i.cols)/(i.rows*i.cols);
    count[i.rows-row_min+1, i.cols-col_min+1] += 1;
end

if false
    figc = Figure(; resolution=(800,700));
    bars3D(count, row_min-1, col_min-1, :deep, create_axis3D("Frequency", figc[1,1], title="Partition count"), figc[1,2]);
    display(figc)
end

# tols = map((diff, real) -> real == 0.0 ? NaN : 0.5*log10(diff/real), sum_diff1, exact);
# ratios = map((s,c) -> convert(Float32, c == 0 ? NaN : s/c), sum_ratio, count);

# deep curl CMRmap
fig = Figure(; resolution=(1200, 1000));

tol_title = L"Tolerance  $\epsilon^2 \;=\; \frac{\sum ||\tilde{A} - A||^2}{\sum ||A||^2}$";
ratio_title = L"Compression ratio  $\nu \;=\; \frac{k(n+m)}{nm}$";
Label(fig[1, 2], tol_title, justification=:center, fontsize=20, tellwidth=false);
Label(fig[1, 4], ratio_title, fontsize=20, tellwidth=false);
Label(fig[2, 1], "Before SVD", rotation = pi/2, halign=:center, tellheight=false);
Label(fig[3, 1], "After SVD", rotation = pi/2, halign=:center, tellheight=false);

tols1 = calcTols(sum_diff1, exact);
tols2 = calcTols(sum_diff1, exact);
ratio1 = averageRatios(sum_ratio1, count);
ratio2 = averageRatios(sum_ratio2, count);

bars3D(tols1, row_min-1, col_min-1, :thermal, create_axis3D(L"\log_{10}(\epsilon)", fig[2,2]), fig[2,3]);
bars3D(tols2, row_min-1, col_min-1, :thermal, create_axis3D(L"\log_{10}(\epsilon)", fig[3,2]), fig[3,3]);
bars3D(ratio1, row_min-1, col_min-1, :thermal, create_axis3D(L"\nu", fig[2,4]), fig[2,5]);
bars3D(ratio2, row_min-1, col_min-1, :thermal, create_axis3D(L"\nu", fig[3,4]), fig[3,5]);

# tol_title = L"\text{Tolerance} \;=\; \sqrt{\left(\sum ||\tilde{A} - A||^2\right)/\left(\sum ||A||^2\right)} \;=\; \epsilon";
# bars3D(tols, row_min-1, col_min-1, :thermal, create_axis3D(tol_title, L"\log_{10}(\epsilon)", fig[1,1]), fig[1,2]);
# bars3D(ratios, row_min-1, col_min-1, :thermal, create_axis3D(L"\text{Compression ratio} \;=\; \frac{k(n+m)}{nm} \;=\; \nu", L"\nu", fig[1,3]), fig[1,4]); 
# bars3D(count, row_min-1, col_min-1, :thermal, create_axis3D("Size count", "Frequency", fig[1,5]), fig[1,6]); 


fig