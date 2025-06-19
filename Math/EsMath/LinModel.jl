import Pkg; Pkg.add ("LinearAlgebra")
using Plots
import Pkg; Pkg.add ("Statistics")
function my_mean(data:: Vector(T}) where T
	Ss = sum(data)
	return ss / length(data)
end

mу_х = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
my_y = [2.20, 2.18, 1.87, 1.85, 1.77, 1.62, 1.57, 1.27, 1.05, -2.08, -2.55, -0.10, -0.41, -1.00, -1.19, -1.56, -2.08, -2.61, -3.37, -3.86] 
my_x = collect(my_x)
my_y = collect(my_y)

lin_x = hcat(ones (length(my_x)), my_x)
lin_b = lin_x \ my_y
lin_y = lin_x * lin_b

quad_x = hcat(ones(length(my_x)), my_х, my_х.^2)
quad_b = quad_x \ my_y
quad_y = quad_x * quad_b
lin_mse = my_mean((my_y .- lin_y).^2）
quad_mse = my_mean((my_y .- quad_y).^2)

println("коэффициенты линейной модели: ", lin_b) 
println("Ско линейной модели: " , lin_mse)
println("Коэффициенты квадратичной модели: ", quad_b)
println("ско квадратичной модели: ", quad mse)

plot(my_x, my_y, seriestype= scatter, label="анные", xabel="x"， ylabel="y"）
plot! (my_x, lin_y, label="Линейная модель", Iw=2)
plot! (my_x, quad_y, label="квадратичная модель", lw=2）