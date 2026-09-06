import sys
import io
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt

# Force UTF-8 encoding for console output
if sys.stdout.encoding.lower() != 'utf-8':
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def run_fel(mode, md1_path, md2_path, md3_path, eps, min_pts, temp):
    # 1. Loading and joining the input data
    df1 = pd.read_csv(md1_path, sep=r'\s+')
    df2 = pd.read_csv(md2_path, sep=r'\s+')
    df3 = pd.read_csv(md3_path, sep=r'\s+')

    df1['Replica'] = 'MD1'
    df2['Replica'] = 'MD2'
    df3['Replica'] = 'MD3'

    data = pd.concat([df1, df2, df3], ignore_index=True)
    coords = data[['RMSD', 'Rg']].values

    # 1. GENERATING k-NN DISTANCE GRAPH (kNNdistplot)
 
    nbrs = NearestNeighbors(n_neighbors=min_pts).fit(coords)
    distances, _ = nbrs.kneighbors(coords)

    k_distances = np.sort(distances[:, -1])

    fig_knn, ax_knn = plt.subplots(figsize=(7, 5), dpi=300)
    ax_knn.plot(k_distances, color='blue', linewidth=1.5, label=f'{min_pts}-NN distance')
    ax_knn.axhline(y=eps, color='red', linestyle='--', linewidth=1.2, label=f'chosen eps = {eps}')

    ax_knn.set_xlabel("Ponits / Frames (sorted by distance)")
    ax_knn.set_ylabel(f"{min_pts}-NN Distance [A]")
    ax_knn.set_title(f"k-NN Distance Plot (minPts = {min_pts})")
    ax_knn.legend(loc='upper left')
    ax_knn.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.savefig("kNN_distance_plot.png", dpi=300)
    plt.close(fig_knn)
    print(f"[Python] Generated k-NN graph: kNN_distance_plot.png (eps={eps}, minPts={min_pts})")

    # If only KNN mode is selected, stop execution here
    if mode == "knn":
        print("[Python] Mode 'knn': Script completed after generating the k-NN graph.")
        return

    # 2. FULL FEL AND DBSCAN ANALYSIS ('full' mode)
    
    print("[Python] Starting full FEL analysis...")

    # 2D Kernel Density Estimation (KDE)
    x = data['RMSD'].values
    y = data['Rg'].values

    k = gaussian_kde([x, y])
    xi_line = np.linspace(x.min(), x.max(), 200)
    yi_line = np.linspace(y.min(), y.max(), 200)
    xi, yi = np.meshgrid(xi_line, yi_line)

    positions = np.vstack([xi.ravel(), yi.ravel()])
    density = k(positions).reshape(xi.shape)

    # Free energy (FE)
    kB = 0.008314  # kJ/mol·K units
    min_nonzero = density[density > 0].min()
    density[density == 0] = min_nonzero

    fe = -kB * temp * np.log(density)
    fe -= fe.min()
    
    f_max = np.quantile(fe, 0.99)
    fe = np.clip(fe, a_min=None, a_max=f_max)

    # DBSCAN clustering
    db = DBSCAN(eps=eps, min_samples=min_pts).fit(coords)
    data['cluster'] = db.labels_

    cluster_stats = []
    rep_structs = []

    for cluster_id in set(data['cluster']):
        if cluster_id == -1:
            continue
            
        c_data = data[data['cluster'] == cluster_id]
        pop = len(c_data)
        rmsd_mean = c_data['RMSD'].mean()
        rg_mean = c_data['Rg'].mean()
        
        prob = pop / len(data)
        fe_cluster = -kB * temp * np.log(prob)
        
        cluster_stats.append({
            'cluster': cluster_id,
            'population': pop,
            'RMSD_mean': rmsd_mean,
            'Rg_mean': rg_mean,
            'FE': fe_cluster
        })

        center = np.array([[rmsd_mean, rg_mean]])
        nn = NearestNeighbors(n_neighbors=1).fit(c_data[['RMSD', 'Rg']].values)
        idx = nn.kneighbors(center, return_distance=False)[0][0]
        rep_structs.append(c_data.iloc[idx])

    df_stats = pd.DataFrame(cluster_stats)
    if not df_stats.empty:
        df_stats['FE'] -= df_stats['FE'].min()
        df_stats = df_stats.sort_values(by='population', ascending=False)
        df_stats.to_csv("Cluster_stats_monomer.txt", sep="\t", index=False)

    if rep_structs:
        df_rep = pd.DataFrame(rep_structs)
        df_rep.to_csv("Rep_frame_monomer.txt", sep="\t", index=False)

    # Visualization FEL plot
    fig, ax = plt.subplots(figsize=(7, 6), dpi=300)
    
    cf = ax.contourf(xi, yi, fe, levels=15, cmap='jet')
    fig.colorbar(cf, ax=ax, label="Free Energy (kJ/mol)")
    ax.contour(xi, yi, fe, levels=15, colors='black', linewidths=0.4, alpha=0.6)

    colors = {'MD1': 'yellow', 'MD2': 'orange', 'MD3': 'magenta'}
    for replica_name, group in data.groupby('Replica'):
        ax.scatter(group['RMSD'], group['Rg'], label=replica_name, 
                   color=colors.get(replica_name, 'white'), s=5, alpha=0.5)

    ax.set_xlabel("RMSD [A]")
    ax.set_ylabel("Radius of Gyration [A]")
    ax.set_title("Monomer Global FEL")
    ax.legend(loc='upper right')

    plt.tight_layout()
    plt.savefig("Global_FEL_RMSD_Rg_monomer.png", dpi=400)
    plt.close(fig)
    print("[Python] Full analysis completed successfully!")

if __name__ == "__main__":
    mode_arg = sys.argv[1]
    md1, md2, md3 = sys.argv[2], sys.argv[3], sys.argv[4]
    eps_val = float(sys.argv[5])
    pts_val = int(sys.argv[6])
    t_val = float(sys.argv[7])

    run_fel(mode_arg, md1, md2, md3, eps_val, pts_val, t_val)