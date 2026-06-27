# diaPASEF best deconvolution package

このフォルダは、diaPASEF Fig.6 デモデータで最も良かった2つの方法を
見返すための最小パッケージです。

## 入っているもの

### Data

```text
data/matrices/fig6_ms2_chromatogram_matrix_rt_by_fragment.tsv
data/matrices/fig6_ms2_mobilogram_matrix_im_by_fragment.tsv
data/fragments/fig6_site_localizing_fragments.tsv
```

### Scripts

```text
scripts/test_diapasef_deconv.R
scripts/run_best_msdial_like.R
scripts/run_best_mcr_als.R
```

`test_diapasef_deconv.R` は検証全体を含むスクリプトです。
`run_best_msdial_like.R` と `run_best_mcr_als.R` は、元の検証スクリプトを実行し、
見るべきベスト画像を表示するためのラッパーです。

### Docs

```text
docs/diapasef_deconv_investigation_summary.md
docs/diapasef_commonA_deconv_notes.md
```

### Results

ベスト2手法の画像:

```text
results/images/best_msdial_like_window_baseline_summary.png
results/images/best_mcr_als_peak_selected_Aonly_summary.png
```

対応するCSV:

```text
results/tables/ms2decr_rt070_im060_window_baseline_commonA_purity.csv
results/tables/ms2decr_rt070_im060_window_baseline_commonA_top_fragments.csv
results/tables/mcr_als_strategy_comparison.csv
results/tables/mcr_als_overcomplete_peak_selected_Aonly_purity.csv
results/tables/mcr_als_overcomplete_peak_selected_Aonly_top_fragments.csv
results/tables/mcr_als_overcomplete_component_table.csv
```

## ベスト結果

### 1. MS-DIAL-like

方法:

```text
RT model: MS2DecR msdial_model, ideal_min = 0.70
IM model: MS2DecR-like mobilogram model, ideal_min = 0.60, detect_fwhm = 3
Fit: peak window + RT/IM baseline terms + shared A
```

結果:

```text
M1 -> S8_phospho purity 0.968
M2 -> S1_phospho purity 0.771
```

画像:

```text
results/images/best_msdial_like_window_baseline_summary.png
```

### 2. MCR-ALS

方法:

```text
fastICA initialization
overcomplete joint MCR-ALS, k = 5, maxiter = 500
A 正規化 (MS2DecR optimizeALS と同じ)
最小誤差 iteration 選択
RT peak picking: wid = 5, max_local_peaks = 2
IM peak picking: wid = 7, max_local_peaks = 3
A は再推定しない (5成分 ALS の結果をそのまま使用)
```

ピークピッキング結果:

```text
成分    RT通過  IM通過  両方通過
ALS1    ○       ×       ×
ALS2    ×       ×       ×
ALS3    ×       ○       ×
ALS4    ○       ○       ○  → S8_phospho
ALS5    ○       ○       ○  → S1_phospho
```

結果:

```text
ALS4 -> S8_phospho purity 0.814
ALS5 -> S1_phospho purity 0.719
```

## 結論

今回の Fig.6 データでは、MS-DIAL-like のモデル抽出と window/baseline つき
共通 `A` 推定が最も良い結果でした。

MCR-ALS は、A 正規化・最小誤差選択・maxiter=500 に修正し、
IM 側の max_local_peaks を 3 に緩めることで S8=0.814 まで改善しましたが、
S1=0.719 は MS-DIAL-like (0.771) に届きませんでした。

なお、ICA 単独（ALS 繰り返しなし）では S1 purity が 0.96〜1.00 と
非常に高くなりますが、この評価は fragment の site ラベル割合のみで、
正解 MS/MS スペクトルとの比較は未実施です。
