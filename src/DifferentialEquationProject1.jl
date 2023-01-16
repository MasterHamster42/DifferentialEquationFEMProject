module DifferentialEquationProject1

using GLMakie

function makeMainMatrix(n::Integer, leftBoundary::Real, rightBoundary::Real)
    numIntervals = n-1
    shortIntervalLen = (rightBoundary-leftBoundary)//numIntervals
    matrix = zeros(Float64, n, n)

    for i in 1:n-1, j in i-1:i+1
        if j == 0 continue end

        shortLeftBoundary = xiVal(max(i, j)-1, shortIntervalLen)
        shortRightBoundary = xiVal(min(i, j)+1, shortIntervalLen)
        
        # calculating integrals
        if shortLeftBoundary < 1 && shortRightBoundary > 1
            matrix[i,j] = gauseLegendre(x -> eFunction(x, 1)*testFunctionDerivative(i, x, shortIntervalLen, n)*testFunctionDerivative(j, x, shortIntervalLen, n), 
                        shortLeftBoundary, 1)
            matrix[i,j] += gauseLegendre(x -> eFunction(x, 1)*testFunctionDerivative(i, x, shortIntervalLen, n)*testFunctionDerivative(j, x, shortIntervalLen, n), 
                        1, shortRightBoundary)
        elseif j == i
            #dividing test function into two linear parts
            shortMiddleBoundary = xiVal(i, shortIntervalLen)

            matrix[i,j] = gauseLegendre(x -> eFunction(x, 1)*testFunctionDerivative(i, x, shortIntervalLen, n)*testFunctionDerivative(j, x, shortIntervalLen, n), 
            shortLeftBoundary, shortMiddleBoundary)
            matrix[i,j] += gauseLegendre(x -> eFunction(x, 1)*testFunctionDerivative(i, x, shortIntervalLen, n)*testFunctionDerivative(j, x, shortIntervalLen, n), 
            shortMiddleBoundary, shortRightBoundary)
        else
            matrix[i,j] = gauseLegendre(x -> eFunction(x, 1)*testFunctionDerivative(i, x, shortIntervalLen, n)*testFunctionDerivative(j, x, shortIntervalLen, n), 
            shortLeftBoundary, shortRightBoundary)
        end

        matrix[i,j] -=eFunction(0, 1)*testFunction(i, 0, shortIntervalLen, n)*testFunction(j, 0, shortIntervalLen, n)
    end

    # setting dirichlet's condition
    matrix[n, n] = 1
    matrix[n-1,n] = 0

    return matrix
end

# i -> 1:n
function xiVal(i::Integer, intervalLen::Real)
    return min(max((i-1)*intervalLen, 0), 2)
end
    
function gauseLegendre(integrateOver::Function, leftBoundary::Real, rightBoundary::Real)
    argument = 1/sqrt(3)
    halfLength = (rightBoundary-leftBoundary)/2
    argumentElongated = argument*halfLength

    return integrateOver(argumentElongated + (rightBoundary+leftBoundary)/2)*halfLength + integrateOver(-argumentElongated + (rightBoundary+leftBoundary)/2)*halfLength
end
    
    
function testFunctionDerivative(i::Integer, x::Real, intervalLen::Real, n::Integer)
    if i > n || i < 0 return 0 end

    if xiVal(i, intervalLen) > x
        return 1/intervalLen
    end

    return -1/intervalLen
end

function testFunction(i::Integer, x::Real, intervalLen::Real, n::Integer)
    if i > n || i < 0 return 0 end
    if x < xiVal(i-1, intervalLen) || x > xiVal(i+1, intervalLen) return 0 end
    
    if xiVal(i, intervalLen) > x
        return 1/intervalLen*(x-xiVal(i-1,intervalLen))
    end
    
    return -1/intervalLen*(x-xiVal(i+1,intervalLen))
end

function eFunction(x::Real, mid::Real) x > mid ? 5 : 3 end


function calculateFunctionValues(weights::Array, leftBoundary::Real, rightBoundary::Real)
    n = length(weights)
    intervalLen = (rightBoundary-leftBoundary)/(n-1)
    Y = Array{Float64}(undef, n)
    X = Array{Float64}(undef, n)
    for i in 1:n
        X[i] = (i-1)*intervalLen
        Y[i] = weights[i]*testFunction(i, X[i], intervalLen, n)
    end
    return X, Y
end 

function solveEquation(n::Integer)
    A = makeMainMatrix(n, 0, 2)
    B = [-10*eFunction(0, 1)*testFunction(i, 0.0, 2/(n-1), n) for i in 1:n]
    B[end] = 0
    weights = A\B
    return calculateFunctionValues(weights, 0, 2)
end


function main(n=10)
    X, Y = solveEquation(n)
    plotResult(X, Y)
end


function plotResult(X, Y)

    f = Figure(resolution = (1290,800), fonts = (; regular= "CMU Serif"))
    ax = Axis(f[1, 1], xlabel = L"x", ylabel = ylabel = L"u(x)",
    xlabelsize = 30, ylabelsize = 30
    )
    Label(f[0,:], L"Elastic\  ddeformation\ ggraph", color = "#F79D1EFF", fontsize = 32, tellwidth=false)
    scatterlines!(ax, X, Y)
    settings_grid = GridLayout(f[2,1])

    input_field = Textbox(settings_grid[1,2], placeholder = "Number of elements", validator = x->(typeof(tryparse(Int64, x)) != Nothing && parse(Int64, x) > 1), tellwidth = false, )
    clear_button = Button(settings_grid[1,3], label="Clear plot", tellwidth=false)
    completiontime_field = Label(settings_grid[1,1], tellwidth=false, text="Execution time:")
    
    
    on(input_field.stored_string) do n
        n = parse(Int64, n)
        execution_time = @elapsed X, Y = solveEquation(n)
        completiontime_field.text[] = "Execution time: $(execution_time)s"
        if n > 50
            lines!(ax, X, Y)
        else
            scatterlines!(ax, X, Y)
        end
    end
    on(clear_button.clicks) do n
        empty!(ax.scene)
    end
    f
end

end # module DifferentialEquationProject1
