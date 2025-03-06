import sys
import math
import numpy as np
from sklearn.metrics import silhouette_score
import symnmf

# kmeans algorithm from HW1
def run_kmeans(k_str, filename):

    eps = 1e-4

    f2 = filename[-4:]
    if f2 != ".txt":
        raise ValueError()

    with open(filename, "r") as f1:
        N = sum(1 for _ in f1)

    with open(filename, "r") as f1:
        d = len(f1.readline().split(","))

    if not k_str.isdigit() or int(k_str) <= 1 or int(k_str) >= N or k_str[0] == "0":
        raise ValueError()

    k_val = int(k_str)

    points = []
    with open(filename, "r") as f1:
        for line in f1:
            points.append(line.strip().split(","))

    centroids = points[:k_val]

    clusters = [[] for _ in range(k_val)]

    def dist(point1, point2):
        s = 0
        for i in range(d):
            s += (float(point1[i]) - float(point2[i]))**2
        return math.sqrt(s)

    def point_to_cluster():
        for cl in clusters:
            cl.clear()
        for point in points:
            nearest = float('inf')
            c_idx = 0
            idx = 0
            for centroid in centroids:
                dis = dist(point, centroid)
                if dis < nearest:
                    nearest = dis
                    c_idx = idx
                idx += 1
            clusters[c_idx].append(point)

    def cluster_avg():
        idx = 0
        flag = 0
        for cl in clusters:
            if not cl:
                idx += 1
                continue
            tmp = []
            l = len(cl)
            for j in range(d):
                s = 0.0
                for i in range(l):
                    s += float(cl[i][j])
                s /= l
                tmp.append(s)
            tmp1 = centroids[idx]
            centroids[idx] = tmp

            if len(tmp1) > 0:
                if dist(tmp1, tmp) > eps:
                    flag = 1
            idx += 1
        if flag == 0:
            return False

    iteration = 300
    for _ in range(iteration):
        point_to_cluster()
        if cluster_avg() == False:
            break

    final_centroids = []
    for c in centroids:
        final_centroids.append([float(x) for x in c])
    return final_centroids


# filter points based on the centroids gotten from Kmeans
def assign_points_to_centroids(data, centroids):
    labels = []
    d = len(data[0])
    for p in data:
        best_dist = float('inf')
        best_idx = 0
        for i, c in enumerate(centroids):
            s = 0.0
            for j in range(d):
                diff = float(p[j]) - float(c[j])
                s += diff*diff
            dval = math.sqrt(s)
            if dval < best_dist:
                best_dist = dval
                best_idx = i
        labels.append(best_idx)
    return labels


def main():
    if len(sys.argv) < 3:
        print("An Error Has Occurred")
        sys.exit(1)

    k_str = sys.argv[1]
    filename = sys.argv[2]

    try:
        final_centroids = run_kmeans(k_str, filename)
    except:
        print("An Error Has Occurred")
        sys.exit(1)

    data = np.loadtxt(filename, delimiter=",", dtype=float)
    labels_km = assign_points_to_centroids(data.tolist(), final_centroids)
    score_km = silhouette_score(data, labels_km)

    W_list = symnmf.norm(data.tolist())
    W = np.array(W_list, dtype=float)

    np.random.seed(1234)
    N = data.shape[0]
    k_val = int(k_str)
    m = float(np.mean(W))
    init_scale = 2.0*math.sqrt(m/k_val)
    H_init = np.random.uniform(0.0, init_scale, size=(N, k_val))

    finalH_list = symnmf.symnmf(W.tolist(), H_init.tolist())
    finalH = np.array(finalH_list, dtype=float)

    labels_nmf = np.argmax(finalH, axis=1)
    score_nmf = silhouette_score(data, labels_nmf)

    print(f"nmf: {score_nmf:.4f}")
    print(f"kmeans: {score_km:.4f}")


if __name__ == "__main__":
    main()
