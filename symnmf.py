import sys
import numpy as np
import symnmf

def main():
    # command line arguments
    if len(sys.argv) != 4:
        print("An Error Has Occurred")
        sys.exit(1)

    try:
        k = int(sys.argv[1])
        goal = sys.argv[2]
        file_name = sys.argv[3]
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    valid_goals = ['symnmf', 'sym', 'ddg', 'norm']
    if goal not in valid_goals:
        print("An Error Has Occurred")
        sys.exit(1)

    # read file
    try:
        X = np.loadtxt(file_name, delimiter=',')
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    np.random.seed(1234)

    # Calculate full symNMF
    if goal == 'symnmf':
        # Compute W
        W = symnmf.norm(X)
        
        # Initialize H
        m = np.mean(W)
        H = np.random.uniform(0, 2 * np.sqrt(m/k), (X.shape[0], k))
        
        # Calculate H
        H_final = symnmf.symnmf(W,H)
        
        # Print H
        for row in H_final:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    # Calculate similarity matrix
    elif goal == 'sym':
        result = symnmf.sym(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    # Calculate diagonal degree matrix
    elif goal == 'ddg':
        result = symnmf.ddg(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))
    
    # Calculate normalized similarity matrix
    elif goal == 'norm':
        result = symnmf.norm(X)
        for row in result:
            print(','.join(['{:.4f}'.format(val) for val in row]))

if __name__ == "__main__":
    main()
