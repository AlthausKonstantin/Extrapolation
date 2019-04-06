
"""
    AdaptiveExtrapolationD

This module is  implementing an
adaptive extrapolation of the explicit midpoint rule according to Deuflhard.
"""
module AdaptiveExtrapolationD



"""
    (Δ,Δx,statisitic) = (mySolver::solver)(f::Function, x₀::Vector{T}, t₀::S, tEnd::S; <options>)
        where {T<:Number,S<:AbstractFloat}

Computes the grid function `xΔ` on the grid `Δ` approximating the solution of the initial value problem

        x′ = f(t,x),  x(t₀) = x₀

on the interval `[t₀,tEnd]` (or `[tEnd,t₀]` if `tEnd < t₀`). The integrator is the extrapolated explicit midpoint rule.
The algorithm adaptively controls stepsize and extrapolation order in the line of thought of Deuflhard.
Note that `solver` is a callable instance of the `struct solver`.

***
## Initialize `solver`


    solver(N::Integer, sequence::String)

Initializes and returns the an instance of `solver`.
`N ≧ 1` is the maximal order of extrapolation. `sequence` specifies the subdividing sequence
used. The options are
* `Harmonic`, that is `1, 2, 3, 4, 5, 6,...`
* `Romberg`, that is `1, 2, 4, 8, 16, 32,...`
* `Bulirsch`, that is `1, 2, 3, 4, 6, 8,... `
The structure contains all quantities for the adaptive extrapolation of the the explicit
midpoint rule for all orders in `1:N` that can be tabulated before the actual computation .
These are the *subdividing sequence* and the *weights* for the extrapolation based on the first barycentric formula.

***
## Options

The following options can be passed as keyword arguments.
* `tol::S > 0` is tolerance used for accuracy check. Default is `1e-3`.
* `relativeScaling::Vector{S}` and `absoluteScaling::Vector{S}` are vectors of the same length as `x₀`
    containing no negative elements.
    Default (in every component) is `relativeScaling[i] = 1.0` and `absoluteScaling[i] = 1e-3`.
    The accuracy check is passed iff the estimated error satisfies in every component:

        |relativeError[i]| ≦ relativeScaling[i]⋅tol & |absoluteError[i]| ≦ absoluteScaling[i]⋅tol

* `firstGuessStepsize::S` is the length of the first step. Its modulus must be greater than
    `eps(S)` Default is `1e-3`.
* `maximalSteplength::S > 0` is the maximal absolute stepsize used for the integration. Default is `Inf`.
* `minimalStepsizeScaling::S` and `maximalStepsizeScaling::S` are safety factors for the stepsize
    selection. Default are `0.25` and `4.0`. The new stepsize satisfies:

         0 ≦ minimalStepsizeScaling ≦ newStepsize / oldStepsize ≦ maximalStepsizeScaling

* `safetyStepsizeSelection::S>0`. Use `safetyStepsizeSelection⋅tol` instead of `tol`. Default is `0.25`.
* `minimumOrder::Int64 ≧ 1` and `maximumOrder::Int64 ≦ N` are the minimal and maximal order of extrapolation for the integration (`N` is the number that has been passed to the constructor). Default are `1` and `N`.
* `firstGuessOrder::Int64` is the extrapolation order for the first step. It must satisfy
    `minimumOrder ≦ firstGuessOrder ≦ maximumOrder`. Default is `round(Int64,(maximumOrder-minimumOrder)/2)`.
* `maximumReduction::Int64` is the upper bound of the number of successive reductions allowed per step.
    Exceeding this boundary leads to ending the integration prematurely (`Δ[end] ≠ tEnd`). Default is `10`.
* `maximumSteps::Int64 ≧ 0` is the maximal length of the array `Δ`. Exceeding this boundary leads to
    ending the integration prematurely (`Δ[end] ≠ tEnd`). Default is `10 000`.
* `rescaleWeights::Bool` If `true` Lagrange interpolation is used otherwise the frist barycentric fromula. Default is `false`
* `detailedStatistics::Bool` flag for the content of `statisitic`. Default is `false`
    + `true`: `statisitic[:,i]` contains [extrapolation order; stepsize; number of `f`-evaluations;
        number of reductions] for the ith step.
    + `false`: `statisitic` contains [total number of evaluations of the right side `f`, total number of reductions].
"""
struct solver
    """
    The structure contains all quantities for the adaptive extrapolation of the the explicit
    midpoint rule that can be tabulated before the acutal computation.
    These are for all orders in `1:N` (`N` is provided the constructor):
    The subdividing sequence, the weights and scaling used for
    the extrapolation based on the first barycentric formula.
    """
   subdividingSequence::Array{BigInt,1}

   # weights and scaling factors for extrapolation operators
   ω::Array{Rational{BigInt},2}
   ρ::Array{Rational{BigInt},1}

   # weights and scaling factors for internal extrapolation operators (used for error estimate)
   ω2::Array{Rational{BigInt},2}
   ρ2::Array{Rational{BigInt},1}

   # constructor
   function solver(N::Integer, sequence::String)
       # check input and initialize subdividing sequence subdividingSequence for barycentric weights
       if N < 0
           error("Order of extrapolation must not be negative. But I got N = $N")
       end
       if sequence == "Harmonic"
           subdividingSequence = [BigInt(n+1) for n = 0:N]
       elseif sequence == "Romberg"
           subdividingSequence = [BigInt(2)^n for n = 0:N]
       elseif sequence == "Bulirsch"
           subdividingSequence = [n==0 ? BigInt(1) : (isodd(n) ? BigInt(2)^Int64(n/2+0.5) : 3*BigInt(2^Int64(n/2-1))) for n = 0:N]
       else
           error("Name of subdividing sequence must be ''Harmonic'', ''Romberg'' or ''Bulirsch''. But I got ''$sequence''")
       end

       # compute nodes corresponding to the subdividing sequence subdividingSequence
       nodes = BigInt(1).// subdividingSequence.^2

       # compute barycentric weights for internal extrapolation operators
       ω2 = zeros(Rational{BigInt},N,N)
       ω2[1,:] = ones(Rational{BigInt},1,N)
       for n = 2:N
           distance = nodes[2:n] .- nodes[n+1]
           ω2[1:(n-1),n] = ω2[1:n-1,n-1] .// distance
           ω2[n,n] = 1 // prod(-distance)
       end

       # compute barycentric weights for extrapolation operators
       ω = zeros(Rational{BigInt},N+1,N+1)
       for n = 1:N
           ω[n+1,(n+1):(N+1)] = ω2[n,n:N] // (nodes[n+1]-nodes[1])
           ω[1,n] = 1 // prod(nodes[1].-nodes[2:n])
       end
       ω[1,N+1] = 1 // prod(nodes[1].-nodes[2:N+1])

       #rescale barycentric weights to obtain weights of 1. barycentric formula
       for m = 1:(N+1)
           ω[1:m,m] = - ω[1:m,m] .// nodes[1:m]
           if 2 <= m
               ω2[1:m-1,m-1] = - ω2[1:m-1,m-1] .// nodes[2:m]
           end
       end

       # compute scaling factors for internal extrapolation operators
       ρ2 = ones(Rational{BigInt},N)
       ρ2[1] = -nodes[2]
       for n = 1:(N-1)
           ρ2[n+1] = -ρ2[n]*nodes[n+2]
       end

       # compute scaling factors for extrapolation operators
       ρ = -nodes[1]*[BigInt(1); ρ2]

       # initialize structure
       new(subdividingSequence,ω,ρ,ω2,ρ2)
   end
end

function (mySolver::solver)(
            f::Function, x₀::Vector{T}, t₀::S, tEnd::S;
            tol::S = parse(S,"1e-3"),
            relativeScaling::Vector{S} =ones(S,size(x₀)),
            absoluteScaling::Vector{S} = fill(parse(S,"1e-3"),size(x₀)),
            firstGuessStepsize::S = parse(S,"1e-3"),
            maximalSteplength::S = parse(S,"Inf"),
            minimalStepsizeScaling::S = parse(S,"0.02"),
            maximalStepsizeScaling::S = parse(S,"4"),
            safetyStepsizeSelection::S =  parse(S,"0.25"),
            minimumOrder::Int64 = 1,
            maximumOrder::Int64 = length(mySolver.subdividingSequence)-1,
            firstGuessOrder = round(Int64,(maximumOrder-minimumOrder)/2),
            maximumReduction::Int64 = 10,
            maximumSteps::Int64 = 10000,
            rescaleWeights::Bool = false,
            detailedStatistics::Bool = false,
            kwargs...
            )where {T<:Number,S<:AbstractFloat}
    # initialize
    # 1. constants
    d = length(x₀) # systems dimension
    if 1<= minimumOrder <= maximumOrder <= length(mySolver.ρ)-1
        n_ex, N_ex = minimumOrder, maximumOrder
    else
        error("Minimal and maximal order must satisfy: 1 ≦ minimumOrder ≦ maximumOrder ≦ $(length(mySolver.ρ)-1)!")
    end

    if length(relativeScaling) == length(absoluteScaling) == d
        if all(el->el>=zero(S),relativeScaling ) && all(el->el>=zero(S),absoluteScaling)
            σᵣ, σₐ = relativeScaling,  absoluteScaling # scaling factors for the  error estimate
        else
            error("All elements of relativeScaling and absoluteScaling must be non-negative!")
        end
    else
        error("relativeScaling and absoluteScaling must be vectors of length $(d)!")
    end

    if zero(S) <= minimalStepsizeScaling <= maximalStepsizeScaling
        p1, p2 = minimalStepsizeScaling, maximalStepsizeScaling
    else
        error("Minimal and maximal stepsize scaling must satisfy: 0.0 <= minimalStepsizeScaling <= maximalStepsizeScaling!")
    end


    # 2. arrays
    if maximumSteps >= 0
        Δ = zeros(S, maximumSteps + 1) # time grid
        xΔ = zeros(T, d, maximumSteps + 1) # solution
    else
        error("maximumSteps must not be negative!")
    end
    X = zeros(T, d, N_ex + 1) # storage for the internal discretisations obtained by the explicit midpoint rule
    λ = zeros(S, N_ex - n_ex +1) # storage for scaling factors of stepsize. λ[k] contains the scalar for extr. order (k + n_ex - 1)
    τ_opt =zeros(S, N_ex - n_ex +1) # storage for optimal stepsize. τ_opt[k] contains the stepsize for extr. order (k + n_ex - 1)
    s = [2*sum(Int64.(mySolver.subdividingSequence[1:n+1])) - n for n in n_ex:N_ex] # s[k] is the number of stages for  extrapolation order (k + n_ex - 1)
    statistics = detailedStatistics ? zeros(S,4,maximumSteps) : zeros(Int64,2)


    # 3. initialization for integration loop
    Δ[1], xΔ[:,1] = t₀, x₀ # store inital values
    tRest = tEnd - Δ[1] # remaining size of the interval
    if maximalSteplength > zero(S)
        maximalSteplength = min(abs(tRest),maximalSteplength)
    else
        error("maximalSteplength must be positive!")
    end

    if abs(firstGuessStepsize) > eps(S)
        τ = flipsign(min(abs(firstGuessStepsize),maximalSteplength),tRest) # stepsize for first step, accomodate backward integration
    else
        error("Absolute value of firstGuessStepsize must be greater than eps($S) = $(eps(S))!")
    end
    if n_ex <= firstGuessOrder <= N_ex
        m = firstGuessOrder # extrapolation order for first step
    else
        error("firstGuessOrder must be between $n_ex and $(N_ex)!")
    end

    counterStep = 1
    counterReduction = 0
    counterEvaluation = 0
    converged = false

    finished = maximumSteps == 0

    x_n = ones(T,d) # storage for the latest solution
    ϵ_n = zero(S) # storage for the latest error estimate
    if safetyStepsizeSelection > 0
        tol = safetyStepsizeSelection*tol # the tolerance for accuracy checks
    else
        error("safetyStepsizeSelection must be positive!")
    end

    # integration loop
    while !finished
        # the nth step: the approximation for the time t₀ + τ is saved in xΔ[n+1].

        # check if the computation must be aborted
        if counterStep > maximumSteps
            println("Not finished. Exceeded maximal number of steps ($maximumSteps)")
            return (Δ[1:counterStep], xΔ[:,1:counterStep], detailedStatistics ? statistics[:,1:counterStep-1] : statistics)
        end
        if counterReduction > maximumReduction
            println("Not finished. Too many reductions in step number $(counterStep)!")
            return (Δ[1:counterStep], xΔ[:,1:counterStep], detailedStatistics ? statistics[:,1:counterStep-1] : statistics)
        end
        if abs(τ) <= eps(S)
            println("Stepsize is becoming too small in step $(counterStep)!")
            return (Δ[1:counterStep], xΔ[:,1:counterStep], detailedStatistics ? statistics[:,1:counterStep-1] : statistics)
        end

        # preperations for current step:
        n_win, N_win = max(n_ex, m - 1), min(N_ex, m + 1) # order window
        n = n_win # start with smalles order in the order window


        # compute all necessary internal approximations:
        f₀ = f(t₀,x₀) # the information for the Euler step
        counterEvaluation = counterEvaluation + 1

        for i = 0:n
            X[:,i+1] = explicitMidpointRule(t₀,x₀,f₀,f,τ,2*Int64(mySolver.subdividingSequence[i+1]))
            counterEvaluation = counterEvaluation + 2*Int64(mySolver.subdividingSequence[i+1]) - 1
        end

        #  compute all information relating to an extrapolation order ≦ n_win:
        for i = n_ex : n
            if rescaleWeights
                # rescaling weights is equivalent to Lagrange interpolation
                x_i= X[:, 1:(i+1)]*T.(broadcast(*,mySolver.ω,mySolver.ρ')[1:(i+1),(i+1)]) # discretisation of order i
                xx_i = X[:, 2:(i+1)]*T.(broadcast(*,mySolver.ω2,mySolver.ρ2')[1:i, i]) # its internal counterpart
            else
                x_i = T.(mySolver.ρ[i+1])*(X[:, 1:(i+1)]*T.(mySolver.ω[1:(i+1), (i+1)]))
                xx_i = T.(mySolver.ρ2[i])*(X[:, 2:(i+1)]*T.(mySolver.ω2[1:i, i]))
            end
            σ = max.(σₐ, σᵣ.*abs.(x_i)) # scaling for error estimate
            ϵ_i = sqrt(sum((abs.(x_i-xx_i)./σ).^2)/d) # scaled mean square norm of the estimated error
            λ[i-n_ex+1] = min(p2, max(p1, (tol/ϵ_i)^(1/(2i+1)))) # factor for optimal stepsize, including safety scaling
            if i == n
                x_n, ϵ_n = x_i, ϵ_i # save solution and error for the current order n
            end
        end

        # check if soltution with for n in the order window can be accepted
        while n <= N_win
            if ϵ_n <= tol
                # accept order n
                converged = true
                break
            elseif ϵ_n <= tol^(s[n-n_ex+1]/s[N_win-n_ex+1])
                # reject order n but pass convergence monitor
                n = n + 1
                # compute x_n, λ_n and ϵ_n for new n:
                X[:,n+1] = explicitMidpointRule(t₀,x₀,f₀,f,τ,2*Int64.(mySolver.subdividingSequence[n+1]))
                counterEvaluation = counterEvaluation + 2Int64(mySolver.subdividingSequence[n+1]) - 1
                if rescaleWeights
                    # rescaling weights is equivalent to Lagrange interpolation
                    x_n= X[:, 1:(n+1)]*T.(broadcast(*,mySolver.ω,mySolver.ρ')[1:(n+1),(n+1)])
                    xx_n = X[:, 2:(n+1)]*T.(broadcast(*,mySolver.ω2,mySolver.ρ2')[1:n, n])
                else
                    x_n = T.(mySolver.ρ[n+1])*(X[:, 1:(n+1)]*T.(mySolver.ω[1:(n+1), (n+1)])) # approximation of extrapolation order n
                    xx_n = T.(mySolver.ρ2[n])*(X[:, 2:(n+1)]*T.(mySolver.ω2[1:n, n])) # and its internal counterpart
                end
                σ = max.(σₐ, σᵣ.*abs.(x_n)) # scaling for error estimate
                ϵ_n = sqrt(sum((abs.(x_n-xx_n)./σ).^2)/d) # scaled mean square norm of the estimated error
                λ[n-n_ex+1] = min(p2,max(p1,(tol/ϵ_n)^(1/(2n+1)))) # factor for optimal stepsize, including safety scaling

            else
                # reject order n and not pass convergence monitor
                break
            end
        end

        if converged
            #  step is accepted, save x_n for the time t₀ + τ and update statistics
            Δ[counterStep+1] = t₀ + τ
            xΔ[:, counterStep+1] = x_n
            tRest = tEnd - Δ[counterStep+1]
            maximalSteplength = min(abs(tRest),maximalSteplength)
            if detailedStatistics
                statistics[:,counterStep] = [n, τ, counterEvaluation, counterReduction]
            else
                statistics = statistics + [counterEvaluation,counterReduction]
            end

            #  check if tEnd is already reached
            if abs(tRest) < eps(S)
                finished = true
            else
                # preperations for the next step:
                t₀, x₀  = Δ[counterStep+1], xΔ[:, counterStep+1] # new inital data

                ## compute optimal order m and stepsize τ for the next step acc. Deuflhard
                # first compute the optimal extrapolation order
                temp = (n_ex:n) .- n_ex .+ 1 # index range of computed quantities
                τ_opt[temp] = min.(abs(τ)*λ[temp], maximalSteplength) # safety boundary for τ_opt
                ω = s[temp]./(τ_opt[temp]) # work per step

                m = argmin(ω) + n_ex - 1 # optimal order

                # check if we may increase m
                if m == n < N_win
                    # compute scaling for order m+1, including safety scaling
                    λ[n-n_ex+2] = min(p2, max(p1, (tol^(s[n-n_ex+1]/s[n-n_ex+2])/ϵ_n)^(1/(2n+1))))
                    # check if work dereases from order m to  m+1
                    if ω[end] > s[n-n_ex+2]/(τ_opt[n-n_ex+2] = min(abs(τ)*λ[n-n_ex+2], maximalSteplength))
                        m = m + 1
                    end
                end
                τ = flipsign(τ_opt[m-n_ex+1],τ) # accomodate backward integration

                counterStep = counterStep + 1
                converged = false
                counterEvaluation = 0
                counterReduction = 0
                fill!(λ,zero(S))
                fill!(τ_opt,zero(S))
            end
        else
            counterReduction = counterReduction + 1
            # compute reduced stepsize τ for order m. Use latest error estimate ϵ_n if estimate of order m is not available
            if n < m
                 λ[m-n_ex+1] = min(p2,max(p1,(tol^(s[n-n_ex+1]/s[m-n_ex+1])/ϵ_n)^(1/(2n+1))))
            end
            τ_opt[m-n_ex+1] = min(abs(τ)*λ[m-n_ex+1],maximalSteplength) # safety boundary for τ_opt
            τ = flipsign(τ_opt[m-n_ex+1],τ) # accomodate backward integration
        end
    end # integration loop
    return (Δ[1:counterStep+1], xΔ[:,1:counterStep+1], detailedStatistics ? statistics[:,1:counterStep] : statistics)
end

"""
    explicitMidpointRule(t₀,x₀,f₀,f,τ,N)

Computes the value ``x_\\Delta(t_0 + \\tau)`` recursively. Here ``x_\\Delta`` is the grid-function generated
by the explicit midpoint rule on the equidistant grid `[t₀ + n/N*τ for n = 0:N]`
with the initial data `(t₀,x₀)` and the right side `f`.
In the starting step `f₀` is used instead of `f(t₀,x₀)`.
"""
function explicitMidpointRule(t₀::S,x₀::Vector{T},f₀::Vector{T},f::Function,τ::S,N) where{S<:AbstractFloat, T<:Number}
    τ_n = τ/N
    x = x₀ + τ_n*f₀
    for i = 0:N-2
        temp = x
        x = x₀+2*τ_n*f(t₀+(i+1)*τ_n,x)
        x₀ = temp
    end
    return x
end

export solver
end # module AdaptiveExtrapolationD
