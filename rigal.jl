using GLM
using GLMakie

function classify_trend(x_init, y;
    significant=0.05,
)
    # Standardise timeline between zero and one
    if length(x_init) <= 2
        return (; class=:none, r2=missing, α0=missing, α1=missing, α2=missing, p_α0=missing, p_α1=missing, p_α2=missing, model=missing)
    end
    xmin, xmax = extrema(x_init)
    x = (x_init .- xmin)# ./ xmax
    # x = Float64.(x_init)
    # ymin, ymax = extrema(y_init)
    # y = (y_init .- ymin) ./ ymax

    model = lm(@formula(y ~ x + x^2), (; x, y))
    ct = coeftable(model)
    α0, α1, α2 = ct.cols[1]
    p_α0, p_α1, p_α2 = ct.cols[ct.pvalcol]

    Xm = mean(x)
    δ = (maximum(x) - minimum(x)) * 0.25

    Y(x) = α0 + x * α1 + α2 * x^2
    Ẏ(x) = α1 + 2 * α2 * x
    γ̇(x) = (-12α2^2 * (2α2 * x + α1)) / (1 + (2α2 * x + α1)^2)^2

    velocity = Ẏ(Xm)
    # Is the squared term sigificant
    constant_threashold = 0.1
    class = if (p_α1 < significant) && (p_α2 < significant) # Nonlinear
        if (Ẏ(Xm - δ) > 0 && Ẏ(Xm + δ) > 0)
            acceleration = sign(γ̇(Xm)) * sign(α2)
            acceleration < 0 ? :accelerated_increase : :decelerated_increase
        elseif (Ẏ(Xm - δ) < 0 && Ẏ(Xm + δ) < 0)
            acceleration = sign(γ̇(Xm)) * sign(α2)
            acceleration < 0 ? :accelerated_decline : :decelerated_decline
        else
            acceleration = sign(α2)
            class = α2 > 0 ? :convex : :concave
        end
    elseif (p_α2 < significant)
        acceleration = sign(α2)
        class = α2 > 0 ? :convex : :concave
    else # Linear
        model = lm(@formula(y ~ x), (; x, y))
        ct = coeftable(model)
        α0, α1 = ct.cols[1]
        p_α0, p_α1 = ct.cols[ct.pvalcol]
        α2 = missing
        p_α2 = missing
        class = if (p_α1 < significant)
            velocity = abs(α1)
            α1 > 0 ? :constant_increase : :constant_decline
        else
            :stable
        end
    end
    (; class, r2=r2(model), α0, α1, α2, p_α0, p_α1, p_α2, model)
end


# @testset "test trends" begin
#     n = 20
#     x = range(0, 1, n)
#     y = rand(n) .+ range(0, 1, n)
#     trends = (; 
#         stable = range(0.5, 0.5, n),
#         constant_increase = range(0, 1, n),
#         constant_decline = range(1, 0, n),
#         concave = sin.(range(0, π, n)),
#         convex = 1 .- sin.(range(0, π, n)),
#         accelerated_increase = range(0, 1, n) .^ 2,
#         accelerated_decline = 1 .- range(0, 1, n) .^ 2,
#         decelerated_increase = 1 .- range(1, 0, n) .^ 2,
#         decelerated_decline = range(1, 0, n) .^ 2,
#     )
#     classifications = map(y -> classify_trend(x, y), trends)
#     classifications.constant_decline
#     keys(trends) == map(x -> x.class, Tuple(classifications))

#     # fig = Figure()
#     # ax = Axis(fig[1, 1])
#     # foreach(trends, keys(trends)) do trend, title
#     #     text = replace(string(title), '_' => ' ')
#     #     Makie.lines!(ax, x, trend; label=text)
#     #     Makie.text!(ax, 0.57, trend[57]; text, align=(:center, :center))
#     # end
# end

