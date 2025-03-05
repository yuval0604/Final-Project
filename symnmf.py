import sys
import numpy as np
import symnmf

def main():
    # Check correct number of arguments
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)

    # Parse arguments
    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    # Validate k and goal
    valid_goals = ['symnmf', 'sym', 'ddg', 'norm']
    if goal not in valid_goals:
        print("An Error Has Occurred")
        sys.exit(1)

    # Read input data
    try:
        X = np.loadtxt(file_name, delimiter=',')
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    # Set random seed as specified
    np.random.seed(1234)

    # Process based on goal
    if goal == 'symnmf':
        # Compute W matrix first
        W = symnmf.norm(X)
        
        # Initialize H randomly
        m = np.mean(W)
        H = np.random.uniform(0, 2 * np.sqrt(m/k), (X.shape[0], k))
        
        # Call C extension to perform symNMF
        H_final = symnmf.symnmf(H, W)
        
        # Print results formatted to 4 decimal places
        for row in H_final:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    elif goal == 'sym':
        result = symnmf.sym(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    elif goal == 'ddg':
        result = symnmf.ddg(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    elif goal == 'norm':
        result = symnmf.norm(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))

if __name__ == "__main__":
    main()
