using LinearAlgebra

xi = [0.235, 0.240, 0.250, 0.255, 0.260, 0.280, 0.295, 0.300, 0.305]
yi = [1.2080, 1.2126, 1.2217, 1.2263, 1.2355, 1.2493, 1.2633, 1.2680, 1.2726]

h = diff(xi)
b = zeros(length(xi))

for i in 2:length(xi)-1
    b[i] = 3 * ((yi[i+1] - yi[i]) / h[i] - (yi[i] - yi[i-1]) / h[i-1])
end

A = zeros(length(xi), length(xi))
A[1,1] = 1
A[end, end] = 1

for i in 2:length(xi)-1
    A[i, i-1] = h[i-1]
    A[i, i] = 2 * (h[i-1] + h[i])
    A[i, i+1] = h[i]
end

M = A \ b

function cubic_spline(x, i)
    h_i = xi[i+1] - xi[i]
    a = (M[i+1] - M[i]) / (6 * h_i)
    b = M[i] / 2
    c = (yi[i+1] - yi[i])/h_i - h_i*(M[i+1] + 2*M[i])/6
    d = yi[i]
    return a * (x - xi[i])^3 + b * (x - xi[i])^2 + c * (x - xi[i]) + d
end

println("Значения сплайна в узловых точках:")
for i in 1:length(xi)-1
    println("S($(xi[i])) = ", cubic_spline(xi[i], i))
end

println("\nЗначения сплайна в произвольных точках: ")
test_points = [0.237, 0.248, 0.265, 0.290, 0.303]
for x in test_points
    for i in 1:length(xi)-1
        if xi[i] <= x <= xi[i+1]
            println("S($x) = ", cubic_spline(x, i))
            break
        end
    end
end