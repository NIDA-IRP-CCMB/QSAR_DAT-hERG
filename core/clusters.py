import operator as op
from functools import reduce
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from yellowbrick.cluster import KElbowVisualizer


# Given input dataframe matrix containing data and range of # of clusters
# Return optimal number of clusters to use
def get_kmeans_k(df, method='silhouette', n_min=2, n_max=20, n_init=50, max_iter=300, random_state=42):
    if n_max > len(df):
        n_max = len(df)
    if method == 'silhouette':
        k_list = []
        for i in range(n_min,n_max):
            kmeans = KMeans(init='random', n_clusters=i, n_init=n_init, max_iter=max_iter, random_state=random_state, n_jobs=-1)
            k_list.append(silhouette_score(df,kmeans.fit_predict(df)))
        n_clusters = k_list.index(max(k_list))+n_min
        return n_clusters
    elif method.lower() == 'sse':
        model = KMeans()
        visualizer = KElbowVisualizer(model, k=(n_min,n_max), timings=False)

        visualizer.fit(df)        # Fit the data to the visualizer
        visualizer.show()        # Finalize and render the figure

        
def kmeans(in_df, n_clusters, n_init=50, max_iter=300, random_state=42):
    df = in_df.copy()
    kmeans = KMeans(init='random', n_clusters=n_clusters, n_init=n_init, max_iter=max_iter, random_state=random_state)
    kmeans.fit(df)
    df["cluster_id"] = kmeans.labels_
    return df.copy()


def ncr(n, k):
    k = min(k, n-k)
    numer = reduce(op.mul, range(n, n-k, -1), 1)
    denom = reduce(op.mul, range(1, k+1), 1)
    return numer // denom