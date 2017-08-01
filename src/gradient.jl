export gradient

function gradient(f::Function, x::Array, tol::Float64)
    Df  = ForwardDiff.gradient(f, x)
    j = 1
    while dot(Df, Df) > tol
        lambda = get_lambda(f, -Df, x)
        x -= lambda*Df
        Df = ForwardDiff.gradient(f, x)
        j+=1
        if j > 1e3
            return x
        end
    end
    return x
end
