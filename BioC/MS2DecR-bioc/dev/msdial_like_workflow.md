# MS-DIAL Like Deconvolution Workflow

## Overview

この文書は、[dev/msdial_like.R](c:/Users/hyama/Documents/R/MS2DecR/dev/msdial_like.R) がどのような手順で
MS-DIAL 論文に近いデコンボリューションを行っているかを整理したものです。

この実装では、`MS1 precursor` は既知であると仮定しています。
そのため、論文中の `peak spotting` のうち `MS1` で precursor spot を見つける工程は省略し、
ユーザーが指定した `premz` と `rt_range` を出発点にしています。

## Input

入力として使う主な情報は次のとおりです。

- `premz`
  対象 precursor の m/z
- `rt_range`
  対象とする保持時間範囲
- `dia_data`
  DIA の元データ

スクリプトではまず `generate_ms2_matrix()` を使って、
指定 precursor / RT 範囲に対応する `MS/MS` 行列を作成します。

- 行方向: 時間
- 列方向: fragment m/z
- 値: 強度

## Analysis Flow

### 1. Focused Peak の決定

まず、MS2 行列 `Y` を時間方向に足し合わせて `total MS2 chromatogram` を作ります。
このクロマトグラムに対して以下を行います。

- 線形重み付き平滑化
- ベースライン補正
- クロマトグラムピーク検出

ここで最も強いピークを `focused peak` とし、
その `left edge`, `apex`, `right edge` を以降の基準にします。

## 2. 各 m/z クロマトグラムの前処理

各 fragment m/z について、対応するクロマトグラムを 1 本ずつ処理します。

各クロマトグラムに対して行う処理は次のとおりです。

- 平滑化
- ベースライン補正
- ピーク検出

ここで得られたピークのうち、
`focused peak` の RT 範囲と重なるものを候補として扱います。

## 3. Model Peak Candidate の選抜

重なったピークそれぞれについて、以下の指標を計算します。

- `ideal slope`
  apex に向かって単調増加し、apex 後に単調減少しているかを評価
- `sharpness`
  ピーク幅に対する鋭さの指標
- `apex_int`
  apex 強度

この実装では、次の条件を満たすものを `model peak candidate` としています。

- `ideal_slope >= 0.95`
- `sharpness >= sharp_min`
- `apex_int >= int_min`

## 4. Sharpness Trace と Matched Filter

候補ピークについて、各 apex scan に `sharpness` を配置し、
`sharpness trace` を作ります。

次に、この配列へ `second Gaussian derivative filter`
（実装上は `matchedfilter()`）を適用し、
`matched wave` を求めます。

この matched wave の極大位置を使って、
デコンボリューションに使うモデルピークの中心候補を決めます。

## 5. M1 / M2 / M3 の決定

MS-DIAL 論文の考え方に合わせて、
以下の 3 本をモデルピークとして扱います。

- `M2`
  focused peak に最も近い中心モデル
- `M1`
  その左隣のモデル
- `M3`
  その右隣のモデル

今回のデータでは、必ずしも 3 本そろうとは限りません。
たとえば左側候補がなければ `M1` は欠損します。

## 6. 最小二乗による分解

各 fragment m/z のクロマトグラムを、
モデルピークとベースライン項の線形結合で表します。

概念的には次の形です。

```text
C(n) = a*M1(n) + b*M2(n) + c*M3(n) + d*n + e
```

ここで、

- `C(n)` は各 fragment m/z の実測クロマトグラム
- `M1(n), M2(n), M3(n)` はモデルクロマトグラム
- `d*n + e` はベースライン補正用の線形項

です。

各 m/z について最小二乗で係数を推定し、
各モデル成分への寄与を求めます。

## 7. 出力される結果

このスクリプトでは、主に次の結果を出力します。

- `model_matrix`
  モデルクロマトグラム `M1/M2/M3`
- `coeff_matrix`
  各 m/z に対する回帰係数
- `spectra`
  各成分ごとに再構成したデコンボリューション後スペクトル
- `reconstructed_chrom`
  各成分ごとに再構成したクロマトグラム

## Output Figures

スクリプト実行後、次の図が保存されます。

- [msdial_like_result.png](c:/Users/hyama/Documents/R/MS2DecR/dev/msdial_like_result.png)
  全体の流れをまとめた図
- [msdial_like_chromatograms.png](c:/Users/hyama/Documents/R/MS2DecR/dev/msdial_like_chromatograms.png)
  モデルクロマトと再構成クロマトの比較
- [msdial_like_components.png](c:/Users/hyama/Documents/R/MS2DecR/dev/msdial_like_components.png)
  各成分について
  左に再構成クロマト、右に対応スペクトルを並べた図

## Important Notes

- この実装は、MS-DIAL 論文の後段を参考にした `MS2Dec-like` な近似実装です。
- `MS1 peak spotting` は省略しています。
- `ideal slope`, `sharpness`, baseline correction は、
  論文本文をもとに再構成した近似であり、Supplementary の式そのものを完全再現したものではありません。
- そのため、MS-DIAL 本体の結果と完全一致することを保証するものではありません。

## Recommended Reading Order

コードを読むときは、次の順に追うと分かりやすいです。

1. パラメータ設定
2. `select_focused_peak()`
3. `extract_model_peak_candidates()`
4. `choose_model_triplet()`
5. `deconvolve_ms2dec()`
6. `result` の中身と出力図
