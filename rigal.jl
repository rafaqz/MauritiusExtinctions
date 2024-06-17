using GLM
using GLMakie

function classify_trend(x, y)
    model = lm(@formula(y ~ x + x^2), (; x, y))
    ct = coeftable(model)
    α0, α1, α2 = ct.cols[1]
    _, p_α1, p_α2 = ct.cols[ct.pvalcol]

    Xm = mean(x)
    δ = (maximum(x) - minimum(x)) * 0.25

    Y(x) = α0 + x * α1 + α2 * x^2
    Ẏ(x) = α1 + 2 * α2 * x
    # Tx_m(x) = fy_deriv1(Xm) .* (x .- x_m) .+ fy(Xm)  

    velocity = Ẏ(Xm)
    # Is the squared term sigificant
    constant_threashold = 0.1
    if (p_α1 < 0.05)
        if p_α2 < 0.05 # Nonlinear
            direction = if (Ẏ(Xm - δ) > 0 && Ẏ(Xm + δ) > 0)
                :increase
            elseif (Ẏ(Xm - δ) < 0 && Ẏ(Xm + δ) < 0)
                :decline
            else
                acceleration = sign(α2)
                if isapprox(α2, 0; atol=1e-10)
                    :stable
                elseif α2 > 0 
                    :convex 
                else
                    @show α2
                    :concave
                end
                :nil
            end
        else # Linear
            velocity = abs(α1)
            if isapprox(α1, 0; atol=1e-10)
                :stable
            elseif α1 > 0
                :constant_increase
            else
                :constant_decline
            end
        end
    else
        acceleration = 0.0
        velocity = 0.0
        @show p_α1 p_α2
        :stable
    end
end

x = range(0, 1, 100)
y = rand(100) .+ range(0, 1, 100)

trends = (; 
    stable = range(0.5, 0.5, 100),
    constant_increase = range(0, 1, 100),
    constant_decline = range(1, 0, 100),
    concave = sin.(range(0, π, 100)),
    convex = 1 .- sin.(range(0, π, 100)),
    accelerated_increase = range(0, 1, 100) .^ 2,
    accelerated_decline = 1 .- range(0, 1, 100) .^ 2,
    decelerated_increase = 1 .- range(1, 0, 100) .^ 2,
    decelerated_decline = range(1, 0, 100) .^ 2,
)
map(y -> classify_trend(x, y), trends) |> pairs

fig = Figure()
ax = Axis(fig[1, 1])
foreach(trends, keys(trends)) do trend, title
    text = replace(string(title), '_' => ' ')
    Makie.lines!(ax, x, trend; label=text)
    Makie.text!(ax, 0.57, trend[57]; text, align=(:center, :center))
end

