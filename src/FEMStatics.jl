# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMStatics.jl/blob/master/LICENSE

""" Nonlinear quasistatic analysis for JuliaFEM. """
module FEMStatics

using Reexport
@reexport using FEMBase

mutable struct Static <: AbstractAnalysis
    t0 :: Float64
    t1 :: Float64
    initial_increment_time :: Float64
    min_increment_time :: Float64
    max_increment_time :: Float64
    max_iterations :: Int64
end

function Static()
    return Static(0.0, 1.0, Inf, 1.0e-5, Inf, 20)
end

function run!(analysis::Analysis{Static})

    p = analysis.properties
    dt = p.initial_increment_time
    time = p.t0
    n_increment = 0
    max_iterations = 2
    done = false
    while !done
        # determine dt
        dt = clamp(dt, p.min_increment_time, p.max_increment_time)
        dt = min(dt, p.t1-time)
        n_increment += 1
        time = min(time+dt, p.t1)
        converged = false
        info("Starting increment $n_increment at time $time.")
        local i
        for i=1:p.max_iterations
            info("Increment # $n_increment. Iteration # $i. Time increment $dt.")
            converged = true
            if converged
                info("Increment # $n_increment converged in $i iterations.")
                break
            end
        end
        if i == p.max_iterations && !converged
            error("Did not converge in $(p.max_iterations) iterations.")
        end
        if isapprox(time, p.t1)
            info("Step completed!")
            done = true
        end
    end

    return nothing
end

export Static

end
