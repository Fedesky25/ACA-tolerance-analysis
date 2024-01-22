module GaussNewton

using LinearAlgebra;

mutable struct FnOut{T}
    y::T
    d::Vector{T}
    FnOut{T}(num::Int) where T = begin
        x = new{T}();
        x.d = Vector{T}(undef, num);
        return x;
    end
end

struct GNData{T}
    x::Vector
    y::Vector
    len::Int
    r::Vector{T}
    J::Matrix{T}
    fn::Function
    betas::Vector{T}
    num::Int
    out::FnOut{T}
    TSS::Float64
end

function createGNData(fn::Function, x::Vector, y::Vector, guess::Vector; T::Type = Real)
    len = length(x);
    if length(y) != len throw("x and y do not have same lenght") end
    num = length(guess);
    return GNData(
        x, y, len,
        Vector{T}(undef, len),
        Matrix{T}(undef, num, len),
        fn, copy!(Vector{T}(undef,num), guess), num,
        FnOut{T}(num),
        totalSumSquare(y)
    );
end

function totalSumSquare(v::Vector) 
    mean = sum(v) / length(v);
    return sum(x -> abs2(x - mean), v);
end

function newguess(data::GNData{T}, guess::Vector{T}) where T
    if length(guess) != data.num throw("new guess has wrong lenght") end
    copy!(data.betas, guess);
end

"""
Calculates the residue of the Gauss-Newton data.

If force is false, the residue is calculated using the current value, 
otherwise the function will be called to compute it.
"""
function residue(data::GNData, force::Bool = false)
    if force 
        map!(
            (y,x) -> begin
                data.fn(x, data.betas, data.out);
                return y - data.out.y;
            end, 
            data.r, 
            data.y, 
            data.x
        );
    end
    return sum(abs2, data.r);
end

function step(data::GNData, α = 1) 
    for i in 1:data.len
        data.fn(data.x[i], data.betas, data.out);
        data.r[i] = data.y[i] - data.out.y;
        data.J[(data.num*(i-1)+1):(data.num*i)] = data.out.d;
    end
    map!((β, Δ) -> β+α*Δ, data.betas, data.betas, qr(transpose(data.J)) \ data.r);
end

function solve(data::GNData; iterations::Int = 5, α = 1)
    for _ ∈ 1:iterations step(data, α); end
    return (data.betas, 1-residue(data, true)/data.TSS);
end

function track(callback::Function, data::GNData{T}; iterations::Int = 5, α=1) where T
    pass = Vector{T}(undef,data.num);
    i = 0;
    while i < iterations
        i += 1;
        for i in 1:data.len
            data.fn(data.x[i], data.betas, data.out);
            data.r[i] = data.y[i] - data.out.y;
            data.J[(data.num*(i-1)+1):(data.num*i)] = data.out.d;
        end
        callback(1-residue(data)/data.TSS, copy!(pass, data.betas));
        map!((β, Δ) -> β+α*Δ, data.betas, data.betas, qr(transpose(data.J)) \ data.r);
    end
    callback(1 - residue(data,true)/data.TSS, copy!(pass,data.betas));
end

export createGNData, step, residue, newguess, solve;

end