module pnl
    Pkg.add("ForwardDiff")
    using ForwardDiff

    include("get_lambda.jl")
    include("pnl_gradient.jl")
    include("pnl_newton.jl")
end