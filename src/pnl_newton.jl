export pnl_newton, pnl_quasinewton

function pnl_newton(f::Function, x::Array, tol::Float64)
    Df = ForwardDiff.gradient(f, x)
    fx = f(x)
    j = 1
    while dot(Df, Df) > tol
        Hf = ForwardDiff.hessian(f, x)
        #if !isposdef(Hf)
        #    return "fail"
        #end
        d = -Hf\Df
        lambda = 1.0
        fxnew = f(x + lambda*d)
        while fxnew > fx + 0.02*lambda*dot(Df, d)
            lambda *= 0.25
            fxnew = f(x + lambda*d)
        end
        x += lambda*d
        fx = fxnew
        Df = ForwardDiff.gradient(f, x)
        j += 1
        if j > 1e3
            return x
        end
    end
    return x
end

function pnl_quasinewton(f::Function, x::Array, tol::Float64)
    Df = ForwardDiff.gradient(f, x)
    fx = f(x)
    H = speye(length(x))
    j = 1
    while dot(Df, Df) > tol
        d = -H*Df
        lambda = get_lambda(f, d, x)
        fxnew = f(x + lambda*d)

        x += lambda*d
        p = lambda*d
        Dfnew = ForwardDiff.gradient(f, x)
        q = Dfnew - Df
        H += (p*p')/(p'*q)[1] - (H*q)*(q'*H)/(q'*H*q)[1]
        Df = Dfnew
        fx = fxnew
        j += 1
        if j > 1e3
            return x
        end
    end
    return x
end