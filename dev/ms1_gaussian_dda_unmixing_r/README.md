# MS1ガウスピーク分解を使ったDDA-MS2スペクトルアンミキシング

このフォルダは、以下のコンセプトを確認するためのRシミュレーションです。

```text
MS1 EICをガウス関数でピーク分解する
→ MS1由来の混合行列Aを作る
→ chimeric DDA-MS2スペクトルを非負制約付きで分解する
```

ロイシン/イソロイシンのような異性体をイメージしたtoy simulationです。実際のLeu/Ileライブラリスペクトルを使った例ではありません。

## コンセプト

ロイシンとイソロイシンのような異性体では、precursor m/zが同じになります。

そのため、MS1では同じm/zのEICとして観測されます。

```text
EIC_total(t) = Leu(t) + Ile(t)
```

もしLC上で2つのピークが少しずれていれば、この合計EICを2本のガウス関数で近似できます。

```text
EIC_total(t) ~= G_Leu(t) + G_Ile(t)
```

DDA-MS2スキャン時刻 `t_j` における2本のガウス関数の高さを、混合行列 `A` として使います。

```text
A[j, Leu] = G_Leu(t_j)
A[j, Ile] = G_Ile(t_j)
```

観測されたDDA-MS2スペクトルは、次のように表します。

```text
Y ~= A S
```

ここで、

- `Y`: 観測されたchimeric DDA-MS2スペクトル行列
- `A`: MS1ガウスピーク分解から得た混合行列
- `S`: 推定したいLeu/IleそれぞれのMS2スペクトル

です。

このスクリプトでは、`A` を固定し、`S` を非負制約付き最小二乗で推定します。

これはMCR-ALSではありません。  
MS1から混合行列 `A` を先に決めて、MS2スペクトル `S` だけを推定する方法です。

## ファイル構成

```text
simulate_leu_ile_unmixing.R
README.md
output/
```

`output/` には、スクリプト実行後のCSVファイルと図が出力されます。

## 必要な環境

base Rのみで動きます。

CRANパッケージは不要です。

## 実行方法

このフォルダ内で実行する場合:

```bash
Rscript simulate_leu_ile_unmixing.R
```

親フォルダから実行する場合:

```bash
Rscript ms1_gaussian_dda_unmixing_r/simulate_leu_ile_unmixing.R
```

出力は、スクリプトの場所を基準にして以下に保存されます。

```text
ms1_gaussian_dda_unmixing_r/output/
```

## シミュレーション設定

同じprecursor m/zを持つ2つのMS1成分を作ります。

```text
Leu apex: 5.08 min
Ile apex: 5.32 min
```

DDA-MS2スキャンは、この2つのピークが重なるRT範囲でスパースに取得される設定です。

各DDA-MS2スペクトルは、以下の式で生成します。

```text
Y_j = G_Leu(t_j) S_Leu + G_Ile(t_j) S_Ile + noise
```

toy fragment m/zは以下です。

```text
30, 41, 44, 69, 86, 132
```

toy spectrumは説明用に単純化しています。

```text
Leu-like: 44と86が強い
Ile-like: 41, 69, 86が強い
```

これらは実際のライブラリスペクトルではありません。

## 出力ファイル

スクリプトを実行すると、以下のCSVが出力されます。

```text
output/ms1_eic.csv
output/design_matrix_A.csv
output/observed_chimeric_ms2.csv
output/estimated_spectra.csv
output/summary.txt
```

また、以下の図が出力されます。

```text
output/01_ms1_eic_and_dda_times.png
output/02_design_matrix_A.png
output/03_true_vs_estimated_spectra.png
output/04_wrong_merged_feature.png
```

## 図の意味

### 01_ms1_eic_and_dda_times.png

Leu/Ileに相当する2本のMS1ガウスピーク、合計EIC、DDA-MS2スキャン時刻を示します。

### 02_design_matrix_A.png

MS1ガウスピークから作った混合行列 `A` を示します。

各DDA-MS2スキャンにおいて、Leu成分とIle成分がどれくらい寄与しているかを表します。

### 03_true_vs_estimated_spectra.png

真のtoy Leu/Ileスペクトルと、chimeric DDA-MS2から推定したスペクトルを比較します。

MS1由来の `A` がうまく作れていれば、真値に近いスペクトルが復元されます。

### 04_wrong_merged_feature.png

対照計算です。

MS1でLeu/Ileを分けず、1本のsame-m/z EICとして扱った場合を示します。

この場合、Leu/Ileを別々には分けられず、混合スペクトルだけが得られます。

## 期待される結果

デフォルトの乱数シードでは、推定スペクトルは真のtoy spectrumにかなり近くなります。

典型的には以下のような結果になります。

```text
Correlation true vs estimated Leu: close to 1
Correlation true vs estimated Ile: close to 1
```

厳密な数値は、Rのバージョンや乱数生成の違いで少し変わる可能性があります。

## 解釈

このシミュレーションで確認したいことは、以下です。

```text
MS1ガウスピーク分解によって、
DDA-MS2スキャン時刻ごとのLeu/Ile混合比を推定できれば、
chimeric DDA-MS2スペクトルを線形アンミキシングできる。
```

一方で、限界もあります。

```text
同じm/zのMS1信号を1本のEICとして扱うと、
MS2側でもLeu/Ileを分離できない。
```

したがって、この方法の重要点はALSではありません。

重要なのは、

```text
MS1から、分解に使える混合行列Aを作れるか
```

です。

## 実データで確認すべきこと

実データに適用する場合は、少なくとも以下を確認します。

- same-m/z EICが2本のガウスピークで説明できるか
- DDA-MS2がピークの異なる位置で取得されているか
- 混合行列 `A` が成分数に対してfull rankか
- `A` のcondition numberが大きすぎないか
- NNLS後の再構成残差が小さいか
- 推定スペクトルが標準品やライブラリと一致するか

2つの異性体が完全共溶出し、MS1由来の混合比がMS2スキャン間で変化しない場合、この方法では分離できません。
