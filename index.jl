using Plots

struct BlockInfo
    norm_diff::Float32
    norm_real::Float32
    rows::UInt
    cols::UInt
    k::UInt
end

# parse the file
f = open("input/Ts.txt", "r")
infos = map(readlines(f)) do line
    c = split(line, " ");
    return BlockInfo(
        parse(Float32, c[1]),
        parse(Float32, c[2]),
        parse(UInt, c[3]),
        parse(UInt, c[4]),
        parse(UInt, c[5]),
    );
end;
close(f)

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
sum_diff = zeros(Float32, row_max-row_min+1, col_max-col_min+1);
sum_real = zeros(Float32, row_max-row_min+1, col_max-col_min+1);
sum_ratio = zeros(Float32, row_max-row_min+1, col_max-col_min+1);
count = zeros(UInt, row_max-row_min+1, col_max-col_min+1);
for i in infos 
    sum_diff[i.rows-row_min+1, i.cols-col_min+1] += i.norm_diff;
    sum_real[i.rows-row_min+1, i.cols-col_min+1] += i.norm_real;
    sum_ratio[i.rows-row_min+1, i.cols-col_min+1] += i.k*(i.rows+i.cols)/(i.rows*i.cols);
    count[i.rows-row_min+1, i.cols-col_min+1] += 1;
end

# surface plot
x = collect(col_min:1:col_max);
y = collect(row_min:1:row_max);

tols = map((diff, real) -> real == 0.0 ? NaN : 0.5*log10(diff/real), sum_diff, sum_real);
ratios = map((s,c) -> c == 0 ? NaN : s/c, sum_ratio, count);
p1 = heatmap(x, y, tols, title = "Tolerance (log10)");
p2 = heatmap(x, y, ratios, title = "Compression");
plot(p1, p2, xlabel = "Columns", ylabel = "Rows", layout=(2,1), size=(500,1000), aspect_ratio=1)
# plot(p1, p2, xlabel = "Columns", ylabel = "Rows", camera=(20,60), xflip=true, layout=(2,1), size=(500,1000))

# Plots.svg("export.svg");