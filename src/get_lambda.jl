export get_lambda

function get_lambda(f::Function, d::Array, x::Array)
    const tol = 1e-3 #deveria depender da norma do gradiente?#
    lambda = 1.0
    g(lambda) = f(x + lambda*d)
    gl = g(lambda)
    Dg = ForwardDiff.derivative(g, lambda)
    j = 1
    while abs(Dg) > tol
        lambda -= gl/Dg
        Dg = ForwardDiff.derivative(g, lambda)
        gl = g(lambda)
        j += 1
        if j > 1e3
            return lambda
        end
    end
    return lambda
end
