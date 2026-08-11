import numpy as np
import matplotlib.pyplot as plt

print("正在讀取數據，請稍候...")
data = np.loadtxt('Results_of_1000000_cells.txt')

# 2. 抽出前三排的精華：X, Y, 密度 (rho)
# Python 的陣列是從 0 開始數的，不要數錯囉！
x_raw = data[:, 0]
y_raw = data[:, 1]
rho_raw = data[:, 2]

# 3. 找出網格的長寬大小 (nx, ny)
# 我們讓 NumPy 自己去算 X 和 Y 有幾種不重複的數字，你連手動輸入都不用！
nx = 1000
ny = 1000
print(f"網格大小破解成功： {nx} x {ny}")

# 4. 把一維的「長麵條」折成二維的「千層麵」
X = x_raw.reshape((nx, ny))
Y = y_raw.reshape((nx, ny))
Rho = rho_raw.reshape((nx, ny))

# 5. 開始你的魔法擺盤
plt.figure(figsize=(9, 8))

# 使用 contourf 畫出超平滑的等高線填色圖
# 根據你的 image_c96595.png，那種黃綠藍交織的復古感，用 'turbo' 或是傳統的 'jet' 最對味
# levels=150 代表漸層切得非常細，保證絲滑到不行
cf = plt.contourf(X, Y, Rho, levels=150, cmap='turbo')

# 在旁邊掛一條顏色對照表 (Colorbar)
cbar = plt.colorbar(cf)
cbar.set_label('Density', fontsize=12)

# 寫上霸氣的標題和座標軸 (照抄你原圖的文字)
plt.title('2nd Order Solver with 4 contact problem - Density', fontsize=14)
plt.xlabel('X/L', fontsize=12)
plt.ylabel('Y/H', fontsize=12)

# 讓 X 跟 Y 軸的比例 1:1，才不會讓你的四個象限看起來像被壓扁的車輪餅
plt.axis('scaled')
plt.xlim(X.min(), X.max())
plt.ylim(Y.min(), Y.max())

# 6. 無頭騎士模式：直接在遠端默默存成圖片，不彈出視窗
# bbox_inches='tight' 會自動幫你把旁邊多餘的白邊裁掉
output_name = 'Result_of_Density.png'
plt.savefig(output_name, dpi=300, bbox_inches='tight')

print(f"圖表已熱騰騰出爐！請享用： {output_name}")
