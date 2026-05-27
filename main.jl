# ==============================================================================
# Minimum Cost Domination Problem
# ==============================================================================

using JuMP
using CPLEX

const MOI = JuMP.MOI

# Direction vectors and movement configurations
const ROOK_DIRS = [(-1, 0), (1, 0), (0, -1), (0, 1)]
const BISHOP_DIRS = [(-1, -1), (-1, 1), (1, -1), (1, 1)]
const QUEEN_DIRS = vcat(ROOK_DIRS, BISHOP_DIRS)
const KNIGHT_MOVES = [(-2, -1), (-2, 1), (-1, -2), (-1, 2), (1, -2), (1, 2), (2, -1), (2, 1)]

"""
    verify_board(board::Matrix{Char})

Verifies if the given n x n board configuration is fully dominated using an 
efficient early-exit scan. Pawns ('P') are assumed to attack upwards.
"""
function verify_board(board::Matrix{Char})
    n = size(board, 1)
    
    for i in 1:n, j in 1:n
        # Pieces dominate their own occupied squares
        if board[i, j] != '.'
            continue
        end

        attacked = false

        # 1. Pawns attack upwards (a pawn at row i+1 attacks row i)
        for (di, dj) in [(1, -1), (1, 1)]
            ni, nj = i + di, j + dj
            if 1 <= ni <= n && 1 <= nj <= n && board[ni, nj] == 'P'
                attacked = true
                break
            end
        end
        attacked && continue

        # 2. Knight attacks
        for (di, dj) in KNIGHT_MOVES
            ni, nj = i + di, j + dj
            if 1 <= ni <= n && 1 <= nj <= n && board[ni, nj] == 'K'
                attacked = true
                break
            end
        end
        attacked && continue

        # 3. Diagonals (Queen & Bishop)
        for (di, dj) in BISHOP_DIRS
            ni, nj = i + di, j + dj
            while 1 <= ni <= n && 1 <= nj <= n
                if board[ni, nj] != '.'
                    if board[ni, nj] == 'Q' || board[ni, nj] == 'B'
                        attacked = true
                    end
                    break # Blocked by any piece
                end
                ni += di; nj += dj
            end
            if attacked
                break
            end
        end
        attacked && continue

        # 4. Straights (Queen & Rook)
        for (di, dj) in ROOK_DIRS
            ni, nj = i + di, j + dj
            while 1 <= ni <= n && 1 <= nj <= n
                if board[ni, nj] != '.'
                    if board[ni, nj] == 'Q' || board[ni, nj] == 'R'
                        attacked = true
                    end
                    break # Blocked by any piece
                end
                ni += di; nj += dj
            end
            if attacked
                break
            end
        end

        if !attacked
            return false # Found an undefended square
        end
    end
    return true
end

"""
    solve_analytical(n::Int)

Generates the constructive Bishop-Pawn configuration for board size n >= 7.

Note: This solution always achieves n - 2 for even and n - 3 for odd, which matches
the true optimum at least up to n = 23.
"""
function solve_analytical(n::Int)
    if n < 7
        error("The analytical solution is known to be non-optimal for n < 7. Please use an ILP method.")
    end
    board = fill('.', n, n)
    cost = 0
    
    if iseven(n)
        c = div(n, 2)
        # Bishop column
        for r in 2:(n - 1)
            board[r, c] = 'B'
        end
        cost += (n - 2) * 3

        # 4 extra Pawns
        board[n, c] = 'P'
        board[1, c] = 'P'
        board[n + 1 - c, n] = 'P'
        board[n - c, n] = 'P'
        cost += 4
    else # odd n >= 7
        c_py = div(n, 2)
        # Bishops columns
        for r in 2:2:(n - 1)
            board[n + 1 - r, c_py + 1] = 'B'
        end
        for r in 2:2:(n - 3)
            board[n + 1 - r, c_py] = 'B'
        end
        cost += (n - 2) * 3

        # 3 extra Pawns
        board[n + 1 - c_py, n - 1] = 'P'
        board[2, c_py] = 'P'
        board[n + 2 - c_py, n] = 'P'
        cost += 3
    end
    return board, cost
end

# Helper to find potential line-of-sight attackers without blocking
function get_unblocked_attackers(i, j, n)
    rooks = Tuple{Int,Int}[]
    bishops = Tuple{Int,Int}[]
    for k in 1:n
        k != j && push!(rooks, (i, k))
        k != i && push!(rooks, (k, j))
    end
    for d in 1:n
        for (di, dj) in [(-d, -d), (-d, d), (d, -d), (d, d)]
            ni, nj = i + di, j + dj
            if 1 <= ni <= n && 1 <= nj <= n
                push!(bishops, (ni, nj))
            end
        end
    end
    return unique(rooks), unique(bishops)
end

# Helper to add ray propagation constraints for the full ILP formulation
function add_ray_constraints!(model, x, o, a, t, dirs, n)
    for (d, (di, dj)) in enumerate(dirs)
        for i in 1:n, j in 1:n
            pi, pj = i - di, j - dj
            if 1 <= pi <= n && 1 <= pj <= n
                @constraint(model, a[i, j, d] == t[pi, pj, d])
            else
                @constraint(model, a[i, j, d] == 0)
            end
            @constraint(model, t[i, j, d] <= a[i, j, d] + x[i, j])
            @constraint(model, t[i, j, d] >= x[i, j])
            @constraint(model, t[i, j, d] <= 1 - o[i, j] + x[i, j])
            @constraint(model, t[i, j, d] >= a[i, j, d] - o[i, j])
        end
    end
end

"""
    solve_chessboard_domination(n::Int; method=:analytical)

Solves the chessboard domination problem using CPLEX.
"""
function solve_chessboard_domination(n::Int; method::Symbol=:analytical)
    if method == :analytical
        board, cost = solve_analytical(n)
        return board, cost, verify_board(board)
    end

    model = Model(CPLEX.Optimizer)
    set_silent(model)

    # Standard piece binary variables
    @variable(model, q[1:n, 1:n], Bin) # Queen (9)
    @variable(model, r[1:n, 1:n], Bin) # Rook (5)
    @variable(model, b[1:n, 1:n], Bin) # Bishop (3)
    @variable(model, k[1:n, 1:n], Bin) # Knight (3)
    @variable(model, p[1:n, 1:n], Bin) # Pawn (1)

    # Minimize total standard weights
    @objective(model, Min, sum(p) + 3*sum(k) + 3*sum(b) + 5*sum(r) + 9*sum(q))

    for i in 1:n, j in 1:n
        @constraint(model, q[i, j] + r[i, j] + b[i, j] + k[i, j] + p[i, j] <= 1)
    end

    if method == :lazy
        for i in 1:n, j in 1:n
            rooks, bishops = get_unblocked_attackers(i, j, n)
            queen_attack = q[i, j] + sum(q[row, col] for (row, col) in rooks; init=0) + sum(q[row, col] for (row, col) in bishops; init=0)
            rook_attack = r[i, j] + sum(r[row, col] for (row, col) in rooks; init=0)
            bishop_attack = b[i, j] + sum(b[row, col] for (row, col) in bishops; init=0)
            
            knight_attack = k[i, j]
            for (di, dj) in KNIGHT_MOVES
                ni, nj = i + di, j + dj
                if 1 <= ni <= n && 1 <= nj <= n
                    knight_attack += k[ni, nj]
                end
            end
            
            # Pawn attacks (pawns at row i+1 attack upwards to row i)
            pawn_attack = p[i, j]
            for (di, dj) in [(1, -1), (1, 1)]
                ni, nj = i + di, j + dj
                if 1 <= ni <= n && 1 <= nj <= n
                    pawn_attack += p[ni, nj]
                end
            end

            @constraint(model, queen_attack + rook_attack + bishop_attack + knight_attack + pawn_attack >= 1)
        end

        function lazy_callback(cb_data)
            status = callback_node_status(cb_data, model)
            status != MOI.CALLBACK_NODE_STATUS_INTEGER && return

            # Optimized bulk query to reduce Julia-C boundary crossings
            q_vals = callback_value.(cb_data, q)
            r_vals = callback_value.(cb_data, r)
            b_vals = callback_value.(cb_data, b)
            k_vals = callback_value.(cb_data, k)
            p_vals = callback_value.(cb_data, p)

            board = fill('.', n, n)
            for r_idx in 1:n, c_idx in 1:n
                if q_vals[r_idx, c_idx] > 0.5
                    board[r_idx, c_idx] = 'Q'
                elseif r_vals[r_idx, c_idx] > 0.5
                    board[r_idx, c_idx] = 'R'
                elseif b_vals[r_idx, c_idx] > 0.5
                    board[r_idx, c_idx] = 'B'
                elseif k_vals[r_idx, c_idx] > 0.5
                    board[r_idx, c_idx] = 'K'
                elseif p_vals[r_idx, c_idx] > 0.5
                    board[r_idx, c_idx] = 'P'
                end
            end

            if !verify_board(board)
                cut_expr = AffExpr(0.0)
                for r_idx in 1:n, c_idx in 1:n
                    # Add terms for each variable to form an integer exclusion cut
                    cut_expr += (q_vals[r_idx, c_idx] > 0.5) ? (1.0 - q[r_idx, c_idx]) : q[r_idx, c_idx]
                    cut_expr += (r_vals[r_idx, c_idx] > 0.5) ? (1.0 - r[r_idx, c_idx]) : r[r_idx, c_idx]
                    cut_expr += (b_vals[r_idx, c_idx] > 0.5) ? (1.0 - b[r_idx, c_idx]) : b[r_idx, c_idx]
                    cut_expr += (k_vals[r_idx, c_idx] > 0.5) ? (1.0 - k[r_idx, c_idx]) : k[r_idx, c_idx]
                    cut_expr += (p_vals[r_idx, c_idx] > 0.5) ? (1.0 - p[r_idx, c_idx]) : p[r_idx, c_idx]
                end
                con = @build_constraint(cut_expr >= 1.0)
                MOI.submit(model, MOI.LazyConstraint(cb_data), con)
            end
        end

        set_attribute(model, MOI.LazyConstraintCallback(), lazy_callback)

    elseif method == :full
        @variable(model, o[1:n, 1:n], Bin) # Occupied matrix
        for i in 1:n, j in 1:n
            @constraint(model, o[i, j] == q[i, j] + r[i, j] + b[i, j] + k[i, j] + p[i, j])
        end

        @variable(model, ar[1:n, 1:n, 1:4], Bin)
        @variable(model, tr[1:n, 1:n, 1:4], Bin)
        @variable(model, ab[1:n, 1:n, 1:4], Bin)
        @variable(model, tb[1:n, 1:n, 1:4], Bin)
        @variable(model, aq[1:n, 1:n, 1:8], Bin)
        @variable(model, tq[1:n, 1:n, 1:8], Bin)

        add_ray_constraints!(model, r, o, ar, tr, ROOK_DIRS, n)
        add_ray_constraints!(model, b, o, ab, tb, BISHOP_DIRS, n)
        add_ray_constraints!(model, q, o, aq, tq, QUEEN_DIRS, n)

        for i in 1:n, j in 1:n
            knight_attack = AffExpr(0.0)
            for (di, dj) in KNIGHT_MOVES
                ni, nj = i + di, j + dj
                if 1 <= ni <= n && 1 <= nj <= n
                    knight_attack += k[ni, nj]
                end
            end

            # Pawn attacks (pawns at row i+1 attack upwards to row i)
            pawn_attack = AffExpr(0.0)
            for (di, dj) in [(1, -1), (1, 1)]
                ni, nj = i + di, j + dj
                if 1 <= ni <= n && 1 <= nj <= n
                    pawn_attack += p[ni, nj]
                end
            end

            sliding_attack = sum(ar[i, j, d] for d in 1:4) +
                             sum(ab[i, j, d] for d in 1:4) +
                             sum(aq[i, j, d] for d in 1:8)

            @constraint(model, o[i, j] + knight_attack + pawn_attack + sliding_attack >= 1)
        end
    end

    optimize!(model)

    board = fill('.', n, n)
    vars = [q, r, b, k, p]
    chars = ['Q', 'R', 'B', 'K', 'P']
    for i in 1:n, j in 1:n
        for (v, char) in zip(vars, chars)
            if value(v[i, j]) > 0.5
                board[i, j] = char
            end
        end
    end

    return board, Int(round(objective_value(model))), verify_board(board)
end

"""
    print_board(board::Matrix{Char})

Renders the board
"""
function print_board(board::Matrix{Char})
    n = size(board, 1)
    for r in 1:n
        println(join(board[r, :], " "))
    end
end

function get_method_name(method::Symbol)
    if method == :analytical
        return "General Solution"
    elseif method == :lazy
        return "ILP with Lazy Constraints"
    elseif method == :full
        return "ILP with Full Encoding"
    else
        return string(method)
    end
end

"""
    main()

Quick demonstration on boards n = 8 to n = 12.
"""
function main()
    methods = [:analytical, :lazy, :full]

    for m in methods
        solve_chessboard_domination(8; method=m)  # Pre-compilation step
    end

    for n in 8:12
        println("----------------------------------------------------------------------")
        println("Board Size: $(n)x$(n)")
        println("----------------------------------------------------------------------")
        for m in methods
            println("Running $(get_method_name(m))")
            
            start_time = time()
            board, weight, verified = solve_chessboard_domination(n; method=m)
            elapsed = time() - start_time
            
            println("Total cost: ", rpad(string(weight), 6), 
                    " | Verified: ", rpad(string(verified), 5), 
                    " | Time: ", round(elapsed, digits=4), "s")
            
            print_board(board)
            println()
        end
        println()
    end
end

main()

# (c) Mia Muessig