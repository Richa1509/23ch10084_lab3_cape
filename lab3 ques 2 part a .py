# --- Physical Constants ---
L = 10.0
T0 = 300.0
TL = 400.0
T_inf = 200.0
alpha = 0.05
beta = 2.7e-9

def solve_linear_system(A, b):
    """Solves Ax = b using Gaussian Elimination"""
    n = len(b)
    # Forward elimination
    for i in range(n):
        pivot = A[i][i]
        for j in range(i + 1, n):
            factor = A[j][i] / pivot
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            b[j] -= factor * b[i]
    # Back substitution
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        sum_ax = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (b[i] - sum_ax) / A[i][i]
    return x

def part_a_fdm(nodes=11):
    """Solves the nonlinear BVP using Finite Difference Method"""
    dx = L / (nodes - 1)
    
    # Initial guess for T: Linear interpolation between T0 and TL
    T = [T0 + (TL - T0) * i / (nodes - 1) for i in range(nodes)]
    
    # Newton-Raphson iterations
    max_iter = 25
    tolerance = 1e-6
    
    for iteration in range(max_iter):
        # F is the residual vector, J is the Jacobian matrix
        F = [0.0] * nodes
        J = [[0.0] * nodes for _ in range(nodes)]
        
        # Boundary Condition at x = 0
        F[0] = T[0] - T0
        J[0][0] = 1.0
        
        # Interior points (1 to nodes-2)
        for i in range(1, nodes - 1):
            # The discretized equation at node i
            # (T[i-1] - 2*T[i] + T[i+1]) / dx^2 - alpha*(T[i] - T_inf) - beta*(T[i]^4 - T_inf^4) = 0
            F[i] = (T[i-1] - 2*T[i] + T[i+1]) / (dx**2) - \
                   alpha * (T[i] - T_inf) - beta * (T[i]**4 - T_inf**4)
            
            # Derivatives for Jacobian Matrix (partial derivatives of F[i] wrt T[i-1], T[i], T[i+1])
            J[i][i-1] = 1.0 / (dx**2)
            J[i][i] = -2.0 / (dx**2) - alpha - 4 * beta * (T[i]**3)
            J[i][i+1] = 1.0 / (dx**2)
            
        # Boundary Condition at x = L
        F[nodes - 1] = T[nodes - 1] - TL
        J[nodes - 1][nodes - 1] = 1.0
        
        # Calculate update step (delta = -J^-1 * F)
        neg_F = [-val for val in F]
        delta = solve_linear_system(J, neg_F)
        
        # Update Temperature values
        for i in range(nodes):
            T[i] += delta[i]
            
        # Check for convergence
        if sum(abs(d) for d in delta) < tolerance:
            break
            
    return T

# --- Execute and Print Results ---
nodes_count = 11
results = part_a_fdm(nodes=nodes_count)

print("-" * 30)
print(f"{'x (Length)':<12} | {'Temp (T)':<10}")
print("-" * 30)
for i in range(nodes_count):
    x_val = (L / (nodes_count - 1)) * i
    print(f"{x_val:<12.1f} | {results[i]:<10.4f}")
print("-" * 30)