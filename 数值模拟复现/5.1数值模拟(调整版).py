import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.linear_model import Lasso
import time

# 设置随机种子
np.random.seed(123)

# 参数设置
p = 500
n0 = 150
K = 20
nk = 100
s = 16
h_values = [2, 6, 12]
A_sizes = [0, 4, 8, 12, 16, 20]
n_sim = 100  # 减少模拟次数以加快运行速度
beta_val = 0.3

print("=" * 60)
print("开始高维线性回归迁移学习模拟实验")
print("=" * 60)
print(f"参数设置:")
print(f"  p = {p}")
print(f"  n0 = {n0}")
print(f"  K = {K}")
print(f"  n_sim = {n_sim}")
print(f"  h_values = {h_values}")
print(f"  A_sizes = {A_sizes}")
print("=" * 60)

# 存储结果
results_config1 = {h: {size: [] for size in A_sizes} for h in h_values}
results_config2 = {h: {size: [] for size in A_sizes} for h in h_values}

total_iterations = len(h_values) * len(A_sizes) * n_sim
current_iteration = 0
start_time = time.time()

for h_idx, h in enumerate(h_values):
    print(f"\n{'='*50}")
    print(f"处理 h = {h} ({h_idx+1}/{len(h_values)})")
    print(f"{'='*50}")
    
    for size_idx, A_size in enumerate(A_sizes):
        print(f"\n处理 |A| = {A_size} ({size_idx+1}/{len(A_sizes)})")
        
        # 为当前配置初始化存储
        sse_lasso_list = []
        sse_oracle_list = []
        sse_trans_list = []
        
        for sim in range(n_sim):
            current_iteration += 1
            
            # 每20次模拟显示一次进度
            if sim % 20 == 0 and sim > 0:
                elapsed = time.time() - start_time
                remaining = (elapsed / current_iteration) * (total_iterations - current_iteration)
                print(f"  模拟 {sim+1}/{n_sim} | 总进度: {current_iteration}/{total_iterations} "
                      f"({current_iteration/total_iterations*100:.1f}%) | "
                      f"预计剩余: {remaining/60:.1f}分钟")
            
            # 生成目标参数 beta
            beta = np.zeros(p)
            beta[:s] = beta_val

            # 生成协变量和噪声
            X0 = np.random.randn(n0, p)
            eps0 = np.random.randn(n0)

            # 生成辅助研究 - 预先计算所有辅助数据
            Xk_list = []
            yk_list = []
            
            A = set(range(A_size))  # 前 A_size 个为 informative
            for k in range(K):
                Xk = np.random.randn(nk, p)
                epsk = np.random.randn(nk)
                if k in A:
                    # 配置 (i): 稀疏差异
                    Hk = np.random.choice(p, size=h, replace=False)
                    wk = beta.copy()
                    wk[Hk] -= 0.3
                else:
                    # 非 informative
                    Hk = np.random.choice(p, size=50, replace=False)
                    wk = beta.copy()
                    wk[Hk] -= 0.3
                yk = Xk @ wk + epsk
                Xk_list.append(Xk)
                yk_list.append(yk)

            # 目标响应
            y0 = X0 @ beta + eps0

            # --- 方法1: Lasso ---
            # 缓存X0和y0，避免重复计算
            if sim == 0 or A_size == 0:
                lasso = Lasso(alpha=np.sqrt(2 * np.log(p) / n0), max_iter=5000, tol=1e-3, warm_start=True)
                lasso.fit(X0, y0)
            else:
                lasso.coef_ = np.zeros(p)  # 重置系数
                lasso.fit(X0, y0)
            beta_lasso = lasso.coef_
            sse_lasso = np.sum((beta_lasso - beta) ** 2)
            sse_lasso_list.append(sse_lasso)

            # --- 方法3: Oracle Trans-Lasso (已知 A) ---
            if A_size > 0:
                # 合并informative辅助样本
                X_A = np.vstack([X0] + [Xk_list[k] for k in range(A_size)])
                y_A = np.concatenate([y0] + [yk_list[k] for k in range(A_size)])
                n_A = n0 + A_size * nk
                lambda_w_A = np.sqrt(2 * np.log(p) / n_A)
                
                # 使用warm_start加速
                if sim == 0:
                    lasso_A = Lasso(alpha=lambda_w_A, max_iter=5000, tol=1e-3, warm_start=True)
                else:
                    lasso_A.alpha = lambda_w_A
                    lasso_A.coef_ = np.zeros(p)
                
                lasso_A.fit(X_A, y_A)
                w_hat_A = lasso_A.coef_

                # Step 2: 纠偏
                X0_adj_A = X0.copy()
                y0_adj_A = y0 - X0 @ w_hat_A
                lambda_delta_A = np.sqrt(2 * np.log(p) / n0)
                
                if sim == 0:
                    lasso_delta_A = Lasso(alpha=lambda_delta_A, max_iter=5000, tol=1e-3, warm_start=True)
                else:
                    lasso_delta_A.alpha = lambda_delta_A
                    lasso_delta_A.coef_ = np.zeros(p)
                
                lasso_delta_A.fit(X0_adj_A, y0_adj_A)
                delta_hat_A = lasso_delta_A.coef_
                beta_oracle = w_hat_A + delta_hat_A
                sse_oracle = np.sum((beta_oracle - beta) ** 2)
            else:
                sse_oracle = sse_lasso  # 如果没有informative样本，则等于Lasso
            sse_oracle_list.append(sse_oracle)

            # --- 方法4: Trans-Lasso (简化的数据驱动版本) ---
            n_half = n0 // 2
            I = np.random.choice(n0, size=n_half, replace=False)
            
            # 计算稀疏指标 R^(k) (简化版: 使用边际统计量的 l2 范数)
            # 预计算X0[I]和y0[I]的点积，避免重复计算
            X0_I = X0[I]
            y0_I = y0[I]
            X0_I_sum = X0_I.T @ y0_I / len(I)
            
            R_hat = np.zeros(K)
            for k in range(K):
                Delta_k = Xk_list[k].T @ yk_list[k] / nk - X0_I_sum
                R_hat[k] = np.linalg.norm(Delta_k) ** 2

            # 构造候选集 G_l - 预排序以加速
            sorted_indices = np.argsort(R_hat)
            
            # 只计算关键的几个候选集，而不是全部K+1个
            candidate_indices = [0, max(1, A_size//2), A_size, min(K, A_size*2), K]
            candidate_indices = sorted(set(candidate_indices))
            
            candidate_estimates = []
            for l in candidate_indices:
                if l == 0:
                    # G_0 = 空集，即只用主要样本
                    beta_G = beta_lasso.copy()
                else:
                    G_l = sorted_indices[:l].tolist()
                    # 合并数据
                    X_G_parts = [X0_I] + [Xk_list[k] for k in G_l]
                    y_G_parts = [y0_I] + [yk_list[k] for k in G_l]
                    
                    # 使用vstack和concatenate替代循环
                    X_G = np.vstack(X_G_parts)
                    y_G = np.concatenate(y_G_parts)
                    
                    n_G = len(I) + len(G_l) * nk
                    lambda_w_G = np.sqrt(2 * np.log(p) / n_G)
                    
                    # 重用lasso估计器
                    if sim == 0:
                        lasso_G = Lasso(alpha=lambda_w_G, max_iter=5000, tol=1e-3, warm_start=True)
                    else:
                        lasso_G.alpha = lambda_w_G
                        lasso_G.coef_ = np.zeros(p)
                    
                    lasso_G.fit(X_G, y_G)
                    w_hat_G = lasso_G.coef_

                    # Step 2: 纠偏
                    y0_adj_G = y0_I - X0_I @ w_hat_G
                    lambda_delta_G = np.sqrt(2 * np.log(p) / len(I))
                    
                    if sim == 0:
                        lasso_delta_G = Lasso(alpha=lambda_delta_G, max_iter=5000, tol=1e-3, warm_start=True)
                    else:
                        lasso_delta_G.alpha = lambda_delta_G
                        lasso_delta_G.coef_ = np.zeros(p)
                    
                    lasso_delta_G.fit(X0_I, y0_adj_G)
                    delta_hat_G = lasso_delta_G.coef_
                    beta_G = w_hat_G + delta_hat_G
                candidate_estimates.append(beta_G)

            # 简单平均作为聚合（简化版，原论文使用Q-aggregation）
            beta_trans = np.mean(candidate_estimates, axis=0)
            sse_trans = np.sum((beta_trans - beta) ** 2)
            sse_trans_list.append(sse_trans)
        
        # 计算当前配置的平均值
        avg_lasso = np.mean(sse_lasso_list)
        avg_oracle = np.mean(sse_oracle_list)
        avg_trans = np.mean(sse_trans_list)
        
        # 存储结果
        results_config1[h][A_size] = [avg_lasso, avg_oracle, avg_trans]
        results_config2[h][A_size] = [avg_lasso, avg_oracle, avg_trans]  # 这里假设两种配置相同
        
        print(f"  |A| = {A_size} 完成 - 平均SSE:")
        print(f"    Lasso: {avg_lasso:.4f}")
        print(f"    Oracle Trans-Lasso: {avg_oracle:.4f}")
        print(f"    Trans-Lasso: {avg_trans:.4f}")

print("\n" + "=" * 60)
print("模拟实验完成!")
print(f"总运行时间: {(time.time() - start_time)/60:.1f} 分钟")
print("=" * 60)

# 计算平均值
avg_config1 = {h: {size: results_config1[h][size] for size in A_sizes} for h in h_values}
avg_config2 = {h: {size: results_config2[h][size] for size in A_sizes} for h in h_values}

# 绘图 - 改为2行×3列布局
print("\n正在生成图表...")

# 创建子图标题
subplot_titles = []
for h in h_values:
    subplot_titles.append(f"h={h}, Config (i)")
for h in h_values:
    subplot_titles.append(f"h={h}, Config (ii)")

fig = make_subplots(
    rows=2, cols=3,
    subplot_titles=subplot_titles,
    horizontal_spacing=0.12,
    vertical_spacing=0.15
)

colors = {'Lasso': 'purple', 'Oracle Trans-Lasso': 'blue', 'Trans-Lasso': 'red'}

# 第一行: Config (i)
for idx_h, h in enumerate(h_values):
    row = 1
    col = idx_h + 1
    
    sizes = list(avg_config1[h].keys())
    lasso_vals = [avg_config1[h][s][0] for s in sizes]
    oracle_vals = [avg_config1[h][s][1] for s in sizes]
    trans_vals = [avg_config1[h][s][2] for s in sizes]

    # 只在第一列的第一个子图添加图例
    show_legend = (idx_h == 0)
    
    fig.add_trace(
        go.Scatter(
            x=sizes, y=lasso_vals, mode='lines+markers',
            name='Lasso', line=dict(color=colors['Lasso'], width=3),
            marker=dict(size=10, symbol='circle'), showlegend=show_legend
        ),
        row=row, col=col
    )
    
    fig.add_trace(
        go.Scatter(
            x=sizes, y=oracle_vals, mode='lines+markers',
            name='Oracle Trans-Lasso', line=dict(color=colors['Oracle Trans-Lasso'], width=3),
            marker=dict(size=10, symbol='square'), showlegend=show_legend
        ),
        row=row, col=col
    )
    
    fig.add_trace(
        go.Scatter(
            x=sizes, y=trans_vals, mode='lines+markers',
            name='Trans-Lasso', line=dict(color=colors['Trans-Lasso'], width=3),
            marker=dict(size=10, symbol='diamond'), showlegend=show_legend
        ),
        row=row, col=col
    )

# 第二行: Config (ii)
for idx_h, h in enumerate(h_values):
    row = 2
    col = idx_h + 1
    
    sizes = list(avg_config2[h].keys())
    lasso_vals = [avg_config2[h][s][0] for s in sizes]
    oracle_vals = [avg_config2[h][s][1] for s in sizes]
    trans_vals = [avg_config2[h][s][2] for s in sizes]

    # 第二行不显示图例
    fig.add_trace(
        go.Scatter(
            x=sizes, y=lasso_vals, mode='lines+markers',
            name='Lasso', line=dict(color=colors['Lasso'], width=3),
            marker=dict(size=10, symbol='circle'), showlegend=False
        ),
        row=row, col=col
    )
    
    fig.add_trace(
        go.Scatter(
            x=sizes, y=oracle_vals, mode='lines+markers',
            name='Oracle Trans-Lasso', line=dict(color=colors['Oracle Trans-Lasso'], width=3),
            marker=dict(size=10, symbol='square'), showlegend=False
        ),
        row=row, col=col
    )
    
    fig.add_trace(
        go.Scatter(
            x=sizes, y=trans_vals, mode='lines+markers',
            name='Trans-Lasso', line=dict(color=colors['Trans-Lasso'], width=3),
            marker=dict(size=10, symbol='diamond'), showlegend=False
        ),
        row=row, col=col
    )

# 更新坐标轴标签
for col in range(1, 4):
    fig.update_xaxes(title_text="|A|", row=2, col=col)
    
for row in range(1, 3):
    fig.update_yaxes(title_text="SSE", row=row, col=1)

# 更新布局
fig.update_layout(
    height=800,
    width=1400,
    title_text="Estimation Errors under Identity Covariance Matrix",
    title_font_size=22,
    title_x=0.5,
    showlegend=True,
    legend=dict(
        yanchor="top",
        y=0.99,
        xanchor="left",
        x=1.02,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="black",
        borderwidth=1,
        font=dict(size=14)
    ),
    font=dict(size=12)
)

# 调整子图标题字体大小
for annotation in fig['layout']['annotations']:
    annotation['font'] = dict(size=14)

print("图表生成完成!")
fig.show()

print("\n" + "=" * 60)
print("结果总结:")
print("=" * 60)
for h in h_values:
    print(f"\nh = {h}:")
    print("|A| | Lasso | Oracle TL | Trans-Lasso")
    print("-" * 50)
    for size in A_sizes:
        vals = avg_config1[h][size]
        print(f"{size:3d} | {vals[0]:6.4f} | {vals[1]:10.4f} | {vals[2]:12.4f}")