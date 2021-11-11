
function random_constraint_match(categories, init, target, result)
    adjust(cat) = count(==(cat), init) - count(==(cat), target)
    catlist(cat) = make_list_of_all_cells(init, cat)

    result .= init

    for cat in categories
        if adjust(cat) > 0
            for sample in 1:adjust(cat)
                cell = rand(CartesianIndices(cellist))
                addtolist(celllist, cell)
                deleteat!(catlist(cat), cell)
                result[cell] = cat
            end
        end
    end
    for cat in categories
        if adjust(cat) < 0
            for sample in 1:-adjust(cat)
                cell = rand(CartesianIndices(cellist))
                result[cell] = cat
                deleteat!(celllist, cell)
            end
        end
    end
end

function growing_clusters(categories, init, target, result)
    adjust(cat) = count(==(cat), init) - count(==(cat), target)
    result .= init
    edgelist = NamedTuple{(:cell, :cat, :neighborhood_cat),Tuple{Int,Int,Int}[]

    while !adjust(cat) == 0
        for cell in CartesianIndices(init) 
            for nb in ((-1, 0), (0, 1), (0, -1), (1, 0))  # nb is short for neighbour
                cat = result[cell]
                nb_cell = Tuple(cell) .+ nb
                neighborhood_cat = result[nb_cell...]
                catright = adjust(cat) > 0
                nbright = adjust(neighborhood_cat) < 0
                if catright && nbright 
                    edge = (; cell, cat, neighborhood_cat)
                    push!(edgelist, edge)
                end 
            end
        end
    end

    if !isempty(edgelist)
        for cat in categories 
            if adjust(cat) < 0
                randomcat = 0
                while adjust(randomcat) <= 0 
                    cell = rand(init)
                    randomcat = result[cell]
                end
                result(cell) = cat
                # adjust(randomcat)-- # decrement 
                # adjust(cat)++ - # increment 
            end
        end 
    end

    while !isempty(edgelist)
        edge = rand(edgelist)
        deleteat!(edgelist, edge) 
        stillover = adjust(edge.cat) < 0
        stillunder = adjust(edge.neighborhood_cat) < 0
        stillsame = result[edge.cell] == edge.cat
        if stillover && stillunder && stillsame 
            # adjust(edge.cat)-- 
            # adjust(edge.neighborhood_cat)++
            result[edge.cell] = edge.neighborhood_cat 
        end
    end 
end

