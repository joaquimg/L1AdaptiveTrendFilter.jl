compute_components! = Vector{Any}(TOTALCOMPONENTS)

compute_components![STEP] = compute_step!
compute_components![SLOPE] = compute_slope!
compute_components![SPIKE] = compute_spike!
compute_components![SIN] = compute_sin!
compute_components![COS] = compute_cos!



function compute_estimate!(IT, β)

    β_new, β0 = stdβ2usualβ(β,IT,d)

    y = zeros(IT.obs)
    for i in IT.components
        compute_components![i](y, β_new, IT)
    end
    for i in IT.components
        for j in 1:IT.obs
            y[j] += β0[i]
        end
    end

    return y
end

function compute_step!(y, β, IT)

    for j in IT.elements[STEP]
        for i in (j+1):IT.obs
            y[i] = y[i] + β[j]
        end
    end
end

function compute_slope!(y, β, IT)

    for j in IT.elements[SLOPE]
        for i in (j+1):IT.obs
            y[i] = y[i] + (i-j)*β[j]
        end
    end
end

function compute_spike!(y, β, IT)

    for i in IT.elements[SPIKE]
        y[i] = y[i] + β[i]
    end
end


function compute_sin!(y, β, IT, d::CDdata) = compute_sin!(y, β[SIN], IT, d.fs)
function compute_sin!(y, β, IT, f)

    for j in IT.elements[SIN]
        for i in 1:IT.obs
            y[i] = y[i]+sin(f[j]*i)*β[j]
        end
    end
end

function compute_cos!(y, β, IT, d::CDdata) = compute_cos!(y, β[COS], IT, d.fc)
function compute_cos!(y, β, IT, f)

    for j in IT.elements[COS]
        for i in 1:IT.obs
            y[i] = y[i]+cos(f[j]*i)*β[j]
        end
    end
end



function stdβ2usualβ(β,IT,d)

    β_new = copy(β)

    β0 = zeros(TOTALCOMPONENTS)


    for i in IT.components
        for j in IT.elements[i]
            β_new[i][j] = β[i][j]/d.σ[i][j]
        end
    end
    for i in IT.components
        for j in IT.elements[i]
            β_0[i] -= β[i][j]*d.μ[i][j]/d.σ[i][j]
        end
    end

    return β_new, β0
end



