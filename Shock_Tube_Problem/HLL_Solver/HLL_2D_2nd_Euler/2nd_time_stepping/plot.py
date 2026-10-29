import numpy as np
import matplotlib.pyplot as plt

# 1. 載入資料 (對應 MATLAB 的 load)
# 提醒：把你 C 語言輸出的檔名替換掉這裡的 'Result.txt'
# 如果你的檔案第一行有英文字母標頭，可以加上 skiprows=1 來跳過
file_name = 'Results_of_5000_cells_X_direction'
try:
    data = np.loadtxt(file_name)
    print("✅ 資料載入成功！準備開獎...")
except FileNotFoundError:
    print(f"❌ 找不到檔案 {file_name}，是不是檔名打錯，或是大雨把它沖走了？")
    exit()

# 2. 提取欄位 (⚠️ 記得 Python 是從 0 開始數的！)
# 列對應: [0: X, 1: Y, 2: rho, 3: u, 4: v, 5: T, 6: P]
X = data[:, 0]       # 取出所有的 X 座標
rho = data[:, 2]     # 取出所有的密度 (rho)

# 3. 設定畫布並畫出密度圖 (對應 MATLAB 的 plot)
plt.figure(figsize=(10, 6)) # 設定一下畫布大小，看起來比較霸氣
plt.plot(X, rho, label='2nd-Order HLL', color='blue', linewidth=1.5)

# 4. 加上 Legend, title 這些你在 MATLAB 會做的事
plt.title("1D Shock Tube - Density Profile", fontsize=16, fontweight='bold')
plt.xlabel("Location (X)", fontsize=14)
plt.ylabel("Density (rho)", fontsize=14)

# 加個網格，看起來專業度瞬間 +10
plt.grid(True, linestyle='--', alpha=0.7) 

# 顯示圖例
plt.legend(fontsize=12)

# 5. 存成高畫質 PNG 檔 (把原本的 plt.show() 刪掉換成這行)
# dpi=300 可以讓你的圖擁有可以放上 Paper 的高清畫質！
plt.savefig('Result_of_Density.png', dpi=300, bbox_inches='tight')

print("✅ 圖檔已經幫你偷偷存進同一個資料夾啦，快去看看！")
