
struct LandUseChange{R,W,PC,PF,PDU,PR,CP,UP} <: NeighborhoodRule{R,W}
    prob_cleared::PC
    prob_urbanised::PF
    prob_deurbanised::PDU
    prob_deregrown::PR
    clearing_pressure::CP
    urbanisation_pressure::UP
end

function DynamicGrids.modifyrule(data, rule::LandUseChange)
    ( forrested, cleared, urban) = rule.land_use_category
    human_pop = aux(data)[:human_pop][currentframe(data)]
    nurban = count(==(urban), data[:LU])
    ncleared = count(==(cleared), data[:LU])
    @set! rule.urbanisation_pressure = population / nurban
    @set! rule.clearing_pressure = population / ncleared
end

function DynamicGrids.applyrule(data, rule::LandUseChange, lu, I)
    ( forrested, cleared, urban) = rule.land_use_category
    hood = neighbourhood(rule)
    if lu == forrested
        if count(==(cleared), hood) > (length(hood) ÷ 2)
            rand() < rule.prob_cleared ? cleared : lu
        else
            lu
        end
    elseif lu == cleared
        if count(==(urban), hood) > (length(hood) ÷ 2)
            rand() < rule.prob_urbanised : urban : lu
        else
            lu
        end
    else
        lu
    end
end

