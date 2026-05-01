# Spectrum & DWT Analysis 

頻譜特徵分析與離散小波轉換（DWT）分析的 Jupyter Notebooks。

## 資料夾結構

```
spectrum_dwt_analysis/
├── spectrum.ipynb        # 頻譜分析
├── dwt.ipynb             # 離散小波轉換分析
├── requirements.txt      # 相依套件
└── README.md
```

## 輸入資料

### spectrum_features_all_parts.csv

兩個 notebook 共用同一份特徵 CSV，請將此檔案放在與 notebook 相同的目錄下再執行。

此 CSV 由 **spectrum.ipynb Section 1** 自動產生，來源為 MongoDB 匯出的 8 個 JSON 分割檔（`ORMongodb.GTECHDaqRawData.json.part1` … `part8`）。每一列代表一筆合格的 DAQ 量測紀錄（雙 channel 均有完整 51200 點）。

#### 欄位一覽（共 82 欄）

##### 1. 識別與品質欄位（13 欄，分析時會排除）

| 欄位名稱 | 型別 | 說明 |
|---|---|---|
| `row_index` | int | 依讀取順序編號（從 1 開始），代表時間序列的順序 |
| `_id` | str | MongoDB ObjectId |
| `ch1_exists` | int (0/1) | CH1 資料是否存在 |
| `ch2_exists` | int (0/1) | CH2 資料是否存在 |
| `ch1_len_original` | int | CH1 原始樣本數 |
| `ch2_len_original` | int | CH2 原始樣本數 |
| `expected_len` | int | 期望樣本數（預設 51200） |
| `ch1_len_ok` | int (0/1) | CH1 樣本數是否等於 expected_len |
| `ch2_len_ok` | int (0/1) | CH2 樣本數是否等於 expected_len |
| `fs_ch1` | float | CH1 取樣頻率（Hz，預設 51200） |
| `fs_ch2` | float | CH2 取樣頻率（Hz，預設 51200） |
| `ch1_spec_len` | int | CH1 單邊頻譜長度（= N/2 + 1） |
| `ch2_spec_len` | int | CH2 單邊頻譜長度 |

##### 2. 頻譜特徵欄位（各 channel 34 欄 × 2 = 68 欄）

以下欄位對 CH1（前綴 `ch1_`）與 CH2（前綴 `ch2_`）各出現一次：

**基本統計（5 欄）**

| 欄位 | 說明 |
|---|---|
| `{ch}_spec_mean` | 頻譜振幅均值 |
| `{ch}_spec_std` | 頻譜振幅標準差 |
| `{ch}_spec_max` | 頻譜振幅最大值 |
| `{ch}_spec_rms` | 頻譜振幅 RMS |

**主頻特徵（3 欄）**

| 欄位 | 說明 |
|---|---|
| `{ch}_dom_freq` | 主頻（Hz），功率最大的頻率 |
| `{ch}_dom_mag` | 主頻振幅 |
| `{ch}_dom_power_ratio` | 主頻功率佔總功率比 |

**頻譜形狀（3 欄）**

| 欄位 | 說明 |
|---|---|
| `{ch}_spec_centroid` | 頻譜重心（Hz） |
| `{ch}_spec_bandwidth` | 頻譜帶寬（加權標準差，Hz） |
| `{ch}_spec_entropy` | 頻譜熵（以 bit 為單位） |

**頻帶能量（12 欄，6 個頻帶 × 2）**

| 頻帶 | `_energy` 欄 | `_ratio` 欄 |
|---|---|---|
| 0–20 Hz | `{ch}_band_0_20_energy` | `{ch}_band_0_20_ratio` |
| 20–100 Hz | `{ch}_band_20_100_energy` | `{ch}_band_20_100_ratio` |
| 100–1000 Hz | `{ch}_band_100_1000_energy` | `{ch}_band_100_1000_ratio` |
| 1000–5000 Hz | `{ch}_band_1000_5000_energy` | `{ch}_band_1000_5000_ratio` |
| 5000–10000 Hz | `{ch}_band_5000_10000_energy` | `{ch}_band_5000_10000_ratio` |
| 10000–20000 Hz | `{ch}_band_10000_20000_energy` | `{ch}_band_10000_20000_ratio` |

> `_energy`：該頻帶內功率之和；`_ratio`：除以總功率

**低高頻比（1 欄）**

| 欄位 | 說明 |
|---|---|
| `{ch}_low_high_energy_ratio` | 0–1000 Hz 能量 / 1000 Hz 以上能量 |

**Top-5 峰值（10 欄）**

| 欄位 | 說明 |
|---|---|
| `{ch}_peak1_freq` … `{ch}_peak5_freq` | 第 1–5 大振幅峰值的頻率（Hz） |
| `{ch}_peak1_mag` … `{ch}_peak5_mag` | 第 1–5 大振幅峰值的振幅 |

> 若峰值不足 5 個，剩餘欄位填 `NaN`

##### 3. 雙 Channel 交叉特徵（3 欄）

| 欄位 | 說明 |
|---|---|
| `ch12_spec_corr` | CH1 與 CH2 頻譜振幅的 Pearson 相關係數 |
| `ch12_spec_abs_mean_diff` | 頻譜振幅絕對差的均值 |
| `ch12_spec_power_ratio` | CH1 總功率 / CH2 總功率 |

#### 欄位數量總結

| 類別 | 欄數 |
|---|---|
| 識別 / 品質（排除用） | 13 |
| CH1 頻譜特徵 | 34 |
| CH2 頻譜特徵 | 34 |
| 雙 Channel 交叉特徵 | 3 |
| **合計** | **84** |

> DWT / 相關分析時會排除 13 個品質欄位，實際用於分析的特徵為 **69 欄**（注意：`ch1_spec_len`、`ch2_spec_len` 也被排除）。

---

## spectrum.ipynb

### 功能說明

| Section | 內容 |
|---------|------|
| 1. 整合 part1–part8 + 畫頻譜圖 | 讀取多個 JSON 資料片段，合併後畫出各 channel 的頻譜圖 |
| 2. 頻譜圖統計分析 | 使用 StandardScaler + IsolationForest 做離群值偵測 |
| 3. 頻譜變數萃取畫時間序列 | 從 CSV 萃取頻譜特徵並畫各變數的時間趨勢 |
| 4. 熱圖分析 | Pearson/Spearman 相關矩陣（greedy 排序 + scipy hierarchical clustering） |

### 輸出

- `correlation_plots/` — 相關矩陣 heatmap PNG 及 CSV
- `correlation_clustered/` — hierarchical clustering 排序後的 heatmap

---

## dwt.ipynb

### 功能說明

| Section | 內容 |
|---------|------|
| 1. 變數相關圖 | Pearson 相關 heatmap，greedy 排序使高相關變數相鄰 |
| 2. Clustered 相關圖 | 以絕對相關值重新排序後的 heatmap |
| 3. DWT 分析 | 對 69 個特徵做 db4 小波分解（level=3），計算 A、D1、D2、D3 間的相關係數 |
| 4. 加輔助線（v1） | Guided 圖：原始訊號 envelope + A 慢變趨勢 + D1/D2 強度；Cycle-level 統計圖 |
| 5. 加輔助線（v2） | 在 v1 基礎上額外加入 rolling max(abs(D1/D2)) 輔助線 |

### 小波設定

| 參數 | 預設值 |
|------|--------|
| `WAVELET` | `db4` |
| `LEVEL` | `3` |
| `USE_ZSCORE_BEFORE_DWT` | `True` |
| `ANCHOR_FEATURE` | `ch2_dom_power_ratio` |

### 輸出

- `all_variables_A_D_relationship/` — 每個變數的 A/D 關係圖 + summary CSV
- `all69_wavelet_guided_and_cycle_level/` — guided 圖 + cycle-level 指標圖（v1）
- `all69_wavelet_guided_and_cycle_level_v2/` — guided 圖 + cycle-level 指標圖（v2）

---

## 安裝相依套件

```bash
pip install -r requirements.txt
```

## 執行環境

- Python 3.8+
- Jupyter Notebook / JupyterLab
