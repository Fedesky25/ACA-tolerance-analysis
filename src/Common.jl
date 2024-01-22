module Common

export BlockInfo, parseBlocksInfos, tolerance, PartitionData, group;

struct BlockInfo
    exact::Float32
    diff1::Float32
    diff2::Float32
    k1::UInt
    k2::UInt
    rows::UInt
    cols::UInt
end

function parseBlocksInfos(filename)
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

function tolerance(data::Vector{BlockInfo}; type::Symbol = :total)
    if type == :total
        tot = sum(i -> i.exact, data);
        s1 = s2 = 0;
        for i in data 
            s1 += i.diff1;
            s2 += i.diff2;
        end
        return (sqrt(s1/tot), sqrt(s2/tot));
    elseif type == :avg
        n1 = n2 = 0;
        for i in data
            n1 += sqrt(i.diff1/i.exact);
            n2 += sqrt(i.diff2/i.exact);
        end
        len = length(data);
        return (n1/len, n2/len);
    elseif type == :rms 
        area = n1 = n2 = 0;
        for i in data
            area += i.rows*i.cols;
            n1 += i.diff1/i.exact*i.rows*i.cols;
            n2 += i.diff2/i.exact*i.rows*i.cols;
        end
        return (sqrt(n1/area), sqrt(n2/area));
    elseif type == :sqrt
        d = n1 = n2 = 0;
        for i in data
            d += sqrt(i.rows*i.cols);
            n1 += sqrt(i.diff1/i.exact*i.rows*i.cols);
            n2 += sqrt(i.diff2/i.exact*i.rows*i.cols);
        end
        return (n1/d, n2/d);
    else
        error("Invalid tolerance selection")
    end
end

struct MinMax
    min::UInt;
    max::UInt;
end

Base.show(io::IO, mm::MinMax) = print(io, "(min $(mm.min), max $(mm.max))");
Base.size(mm::MinMax) = mm.max - mm.min;

struct PartitionData
    rows::MinMax;
    cols::MinMax;
    data::Vector{BlockInfo};

    PartitionData(data::Vector{BlockInfo}) = begin
        if length(data) == 0
            error("Empty vector of block info");
        end
        r_min = c_min = typemax(UInt);
        r_max = c_max = 0;
        for i in data
            if i.rows < r_min
                r_min = i.rows;
            elseif i.rows > r_max
                r_max = i.rows;
            end
            if i.cols < c_min
                c_min = i.cols;
            elseif i.cols > c_max
                c_max = i.cols;
            end
        end
        return new(MinMax(r_min, r_max), MinMax(c_min, c_max), data);
    end
end

Base.show(io::IO, d::PartitionData) = print(io, "PartitionData{rows: ", d.rows, ", cols: ", d.cols, ", ", length(d.data), " BlockInfo}");
Base.size(d::PartitionData) = (d.rows.max-d.rows.min+1, d.cols.max-d.cols.min+1);

function group(data::PartitionData)
    dims = size(data);
    mat = Matrix{Union{Nothing, Vector{BlockInfo}}}(nothing, dims);
    for d in data.data
        i = d.rows - data.rows.min + 1;
        j = d.cols - data.cols.min + 1;
        if mat[i,j] === nothing
            mat[i,j] = [ d ];
        else
            push!(mat[i,j], d);
        end
    end
    return mat;
end

end