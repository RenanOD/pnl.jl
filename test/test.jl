using Base.Test

function test_gradient(f, x0, exact, tol)
    approx = gradient(f, x0, tol)
    @test_approx_eq_eps(dot((approx - exact), (approx - exact)), 0.0, tol*1e3)
end

function test_newnewton(f, x0, exact, tol)
    approx = newnewton(f, x0, tol)
    @test_approx_eq_eps(dot((approx - exact), (approx - exact)), 0.0, tol*1e3)
end

function test_quasinewton(f, x0, exact, tol)
    approx = quasinewton(f, x0, tol)
    @test_approx_eq_eps(dot((approx - exact), (approx - exact)), 0.0, tol*1e3)
end

function newtest()

    F = Function[]
    x0 = Array[]
    exact = Array[]
    tol = 1e-6

    push!(F, x-> x[1]^2 + x[2]^2 + x[3]^2)
    push!(x0, [1.0, -2.0, 0.5])
    push!(exact, [0.0, 0.0, 0.0])

    push!(F, x-> 100*(x[2]-x[1]^2)^2 + (1-x[1])^2)
    push!(x0, [1.5, 0.7])
    push!(exact, [1.0, 1.0])


    push!(F, x-> (x[1]^3 + x[2]^2)^2 + 2*(x[2] - x[1] - 4)^4)
    push!(x0, [-1.0, 0.5])
    push!(exact, [-1.7281632, 2.27183680])

    for i = 1:3
        test_gradient(F[i], x0[i], exact[i], tol)
        test_newnewton(F[i], x0[i], exact[i], tol)
        test_quasinewton(F[i], x0[i], exact[i], tol)
    end
end
newtest()
