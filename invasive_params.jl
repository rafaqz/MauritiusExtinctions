using Graphs, GLMakie, GraphMakie, NetworkLayout, InvertedIndices
using CSV, DataFrames

# convenience functions
# picks n random numbers between low and high, that are no more than a factor 2 apart
function rejection_range(n, low = 0, high = 1, max_multiplier = 2)
    while true
        cont = false
        r = rand(n) .* (high - low) .+ low
        for i in 1:n, j in 1:n
            if r[i] / r[j] > max_multiplier || r[j] / r[i] > max_multiplier
                cont = true
            end
        end
        cont || return r
    end
end

# visualize a single interaction matrix
function plot_interactions(mat, nodes)
    vals = [mat[j, i] for i in axes(mat, 1), j in axes(mat, 2) if i != j]
    cols = ifelse.(vals .>= 0, :blue, :red)
    g = complete_digraph(5)
    f = Figure(size = (800, 800))
    graphplot(f[1,1], g, layout = Shell(), nlabels = nodes, edge_width = 3abs.(vals), edge_color = cols)

    return f
end

# visualize an aggregate from a vector of interaction matrices
function plot_densities(mat, nodes)
    mm = stack(mat)
    f = Figure(size = (900, 1200))
    for i in 1:5, j in 1:5
        if i != j
            col = sum(mm[i,j,:]) < 0 ? :red : :blue
            ax = Axis(f[i,j], title = "$(nodes[i]) on $(nodes[j])", limits = (extrema(mm)..., nothing, nothing))
            GLMakie.density!(ax, mm[i,j,:], color = (col, 0.3), strokearound = true, strokewidth = 1, strokecolor = col)
        end
    end
    f
end

# Metadata
# const nodes = ["cat", "black", "norway", "mouse", "pig"]
# const bodysizes = [4000, 184, 346, 35, 169000]

function interaction_matrix(bodysizes)
    # Create interaction matrix
    cat = 1
    black = 2
    norway = 3
    mouse = 4
    pig = 5

    ints = zeros(5,5)

    # cats
    cat_order = [norway, black, mouse] # these are the relative rankings of the effect of cats on that animal
    ints[cat, black:mouse] = -sort(rand(3))[sortperm(cat_order)]
    ints[black:mouse, cat] = [-ints[cat,i] * sqrt(bodysizes[i]) / sqrt(bodysizes[cat]) for i in black:mouse]
    ints[cat, pig] = -rand() * 0.5minimum(abs.(ints[cat, black:mouse]))
    ints[pig, cat] = -rand() * 0.5minimum(abs.(ints[cat, black:mouse]))

    # rodents
    competition_to_predation_max = 0.3
    randvals = -rejection_range(6, 0, competition_to_predation_max*abs(ints[cat, black]), 2)
    for i in black:mouse, j in black:mouse
        if i != j
            ints[i, j] = pop!(randvals) * sqrt(bodysizes[i]) / sqrt(bodysizes[j])
        end
    end

    # pigs
    randpigs = -rejection_range(6, 0, minimum(abs(x) for x in ints if x!= 0), 2)
    for i in black:mouse
        ints[i, 5] = pop!(randpigs) 
        ints[5, i] = pop!(randpigs) 
    end

    return ints
end

pred_df = CSV.read("tables/animals.csv", DataFrame)
nodes = pred_df.name[1:5]
bodymass = pred_df.mass[1:5]

# one matrix
ints = interaction_matrix(bodymass)
plot_interactions(ints, nodes)

# multiple
mat = map(_-> interaction_matrix(bodymass), 1:100)
plot_densities(mat, nodes)

