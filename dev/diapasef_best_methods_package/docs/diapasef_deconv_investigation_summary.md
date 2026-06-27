# diaPASEF デコンボリューション検証まとめ

## 目的

diaPASEF の MS/MS データでは、RT 方向のクロマトグラムと IM 方向の
モビログラムがあります。今回の基本モデルは次の形です。

```text
X = C A
Y = M A
```

ここで、

- `X`: RT x fragment のクロマトグラム行列
- `Y`: IM x fragment のモビログラム行列
- `C`: モデルクロマトグラム
- `M`: モデルモビログラム
- `A`: 共通の推定 MS/MS スペクトル

`C` と `M` が固定されているなら、基本的な共通スペクトル `A` は
最小二乗法で次のように推定できます。

```text
A = solve(C'C + M'M)(C'X + M'Y)
```

検証には、ローカルにある diaTracer Fig.6 のリン酸化ペプチド例を
使いました。fragment 列には `S1_phospho` と `S8_phospho` のラベルが
付いているので、これを実用上の正解ラベルとして評価しました。

## 関連ファイル

メインの検証スクリプト:

```text
C:\Users\hyama\Documents\diaPASEF\MS2DecR-main\dev\test_diapasef_deconv.R
```

サマリー出力フォルダ:

```text
C:\Users\hyama\Documents\diaPASEF\MS2DecR-main\dev\diapasef_deconv_summary
```

最終的に一番重要な画像:

```text
C:\Users\hyama\Documents\diaPASEF\MS2DecR-main\dev\diapasef_deconv_summary\ms2decr_rt070_im060_window_baseline_deconvolution_summary.png
```

## 試したこと

### 1. site ごとの平均モデル

最初の MS-DIAL-like 試作では、既知の site ラベルごとに強い fragment
trace を選び、それらを平均してモデルクロマトグラムとモデルモビログラムを
作りました。

これは共通 `A` の式が動くかを見るための簡易テストで、MS2DecR 本体の
`msdial_model()` と同じ方法ではありません。

結果:

```text
S1 purity: 0.797
S8 purity: 0.810
```

数値だけ見ると悪くありませんでしたが、RT 側のモデルクロマトグラムが
きれいな単峰になっていませんでした。特に S8 側は複数ピークを含んでいる
ように見えました。

### 2. MS2DecR の `msdial_model()` を厳しめの条件で使う

次に、MS2DecR 本体の `msdial_model()` を使ってモデルクロマトグラムを
抽出しました。最初はデフォルトに近い、やや厳しめの条件です。

```text
ideal_min = 0.90-0.95
detect_fwhm = 5
sharp_min = 0.02
```

この条件では、モデルは1本しか選ばれませんでした。

```text
M2: mz 1233.555, apex RT 18.479
```

MS2DecR 本体の流れで、

```text
msdial_model() -> msdial_deconv()
```

を実行した結果:

```text
M2 -> S1_phospho purity 0.774
```

つまり、MS2DecR の使い方自体は間違っていませんでしたが、
`ideal_min` が厳しすぎて、S8 側のモデルクロマトグラムが選ばれませんでした。

### 3. RT 側の `ideal_min` を緩める

RT 側のモデルクロマトグラム抽出で `ideal_min` を下げました。

```text
RT ideal_min = 0.70
detect_fwhm = 5
```

この条件では、2本のモデルクロマトグラムが選ばれました。

```text
M1: mz 788.3589, apex RT 18.27225  S8-like
M2: mz 654.3333, apex RT 18.45602  S1-like
```

MS2DecR 本体の `msdial_deconv()` での結果:

```text
M1 -> S8_phospho purity 0.931
M2 -> S1_phospho purity 0.772
```

この結果から、RT 側のモデルクロマトグラム抽出は、パラメータを調整すれば
かなりうまくいくことが分かりました。

### 4. モデルモビログラム

モデルモビログラムも、モデルクロマトグラムと同じ考え方で抽出できるはずです。
つまり、各 fragment の mobilogram を 1D trace として扱い、
ピーク候補を検出し、apex 位置でグルーピングし、代表モデルを選びます。

今回の IM 側で有効だった設定:

```text
IM ideal_min = 0.60
detect_fwhm = 3
sharp_min = 0.02
```

この条件では、次のモデルモビログラムが選ばれました。

```text
M2: mz 654.3333, IM 0.9410  S1-like
M3: mz 694.8160, IM 0.9635  S8-like
```

RT 側と IM 側のモデルは、代表 fragment の site ラベルを使って対応づけました。

## うまくいかなかったこと

RT と IM のモデルを用意したあと、単純に全範囲で次の式を解くと、
スペクトルが大きく混ざりました。

```text
A = solve(C'C + M'M)(C'X + M'Y)
```

結果:

```text
M1 -> S1 0.515 / S8 0.485
M2 -> S1 0.720 / S8 0.280
```

このことから、基本式は数学的には正しいものの、実データに対して
全 RT 範囲・全 IM 範囲をそのまま使うのは粗すぎることが分かりました。

全範囲を使うと、目的ピーク以外のシグナル、ピーク外のノイズ、
ベースライン、モデルでは説明できない小さなズレまで、すべて `A` で
説明しようとしてしまいます。その結果、推定スペクトルが混ざります。

## 最終的に良かった方法

最も良かったのは、共通 `A` の考え方を保ちつつ、MS2DecR の
`msdial_deconv()` に近い形で、ピーク window と baseline 項を入れる方法です。

各 fragment について、次のような同時回帰を解きました。

```text
X_window = C_window A + RT baseline
Y_window = M_window A + IM baseline
```

ここで、モデル成分の係数 `A` は RT と IM で共通です。
一方、RT baseline と IM baseline は別々の nuisance parameter として扱います。

最終設定:

```text
RT model:
  ideal_min = 0.70
  detect_fwhm = 5

IM model:
  ideal_min = 0.60
  detect_fwhm = 3

Spectrum fit:
  モデル C または M が非ゼロの peak window だけを使う
  RT slope + RT intercept を入れる
  IM slope + IM intercept を入れる
```

最終結果:

```text
M1 -> S8_phospho purity 0.968
M2 -> S1_phospho purity 0.771
```

これは今回の検証で最も良い結果でした。

## 結論

全体の方針は良いです。

```text
1. モデルクロマトグラム C を抽出する
2. モデルモビログラム M を抽出する
3. RT モデルと IM モデルを対応づける
4. 共通スペクトル A を最小二乗法で推定する
```

ただし、最後の最小二乗推定では、全 RT 範囲・全 IM 範囲をそのまま使うのではなく、
ピーク window と baseline 項を入れた方が良いです。

今回の Fig.6 データでの実用的な結論:

```text
MS2DecR-style の RT モデル抽出は、ideal_min を緩めるとうまくいく。
モデルモビログラムも、同じ 1D peak picking の考え方で抽出できる。
共通 A 推定は、peak window + baseline 項ありで解くとよく分離する。
```

## 残タスク

- RT と IM のパラメータ選択を手動ではなく自動化する。
- 現在は代表 fragment の site ラベルで RT/IM モデルを対応づけているが、
  より一般的な対応づけ方法を考える。
- window + baseline つきの共通 `A` 推定を、きれいな関数として切り出す。
- Fig.6 デモデータに対して、S8/S1 purity を確認するテストを追加する。

## MCR-ALS 側の追加検証

MCR-ALS についても、MS2DecR 本体の考え方に近づけて追加検証しました。

MS2DecR 本体の `deconvICA()` は、概念的には次の流れです。

```text
1. ICA で多めの成分を出す
2. ALS で C と A を最適化する
3. detectPeaks(C, A, wid, gap) でピークらしい成分を残す
4. 選ばれた C_fin, A_fin を返す
```

diaPASEF では `C` だけでなく `M` もあるので、今回の MCR-ALS 試作では
次のようにしました。

```text
1. fastICA 初期値で k = 5 の overcomplete 成分を作る
2. joint MCR-ALS で X = C A, Y = M A を同時に最適化する
3. C と M の両方でピークらしい成分を選ぶ
4. 選択後の扱いを3通りで比較する
```

比較した3通りは次の通りです。

```text
A. peak selection のみで終了
B. 選択した C/M を固定し、window + baseline つきで A だけ再推定
C. 選択した C/M からもう一度 ALS
```

結果:

```text
method                              S1_best  S8_best  mean
A_peak_select_only                  0.718    0.773    0.745
B_fixed_CM_window_baseline_A_only   0.742    0.769    0.756
C_reALS_after_peak_select           0.745    0.633    0.689
MS-DIAL-like reference              0.771    0.968    0.870
```

MCR-ALS 系では、B の「選択した `C/M` を固定して `A` だけ再推定」が最も良い
結果でした。一方、選択後に再度 ALS を行う C は、S8 側が悪化しました。
これは、ALS が site purity ではなく再構成誤差を下げる方向に動くため、
再最適化で成分が再び混ざった可能性があります。

MCR-ALS 系のベスト画像:

```text
C:\Users\hyama\Documents\diaPASEF\MS2DecR-main\dev\diapasef_deconv_summary\mcr_als_overcomplete_peak_selected_Aonly_deconvolution_summary.png
```

ただし、今回の Fig.6 データでは、最終的な性能は MS-DIAL-like の
window + baseline 共通 `A` 推定が最も良好でした。

## ICA 関数の修正

MS2DecR 本体の `performICA()` は、以前は存在しない `ica()` 関数を直接呼んでいました。
`loadings` パッケージには `ica()` が存在しないため、このままでは
`deconvICA()` が動きません。

現在の作業ツリーでは、`performICA()` を `fastICA::fastICA()` 必須に修正しました。
SVD fallback は使わない方針です。

修正後の方針:

```text
performICA() は fastICA::fastICA() を必ず使う
fastICA が無ければ明示的にエラーにする
```

GitHub に反映すべきファイルは次の通りです。

通常パッケージ側:

```text
MS2DecR-main/R/performICA.R
MS2DecR-main/DESCRIPTION
MS2DecR-main/NAMESPACE
MS2DecR-main/man/performICA.Rd
```

Bioconductor 用コピーも管理対象にするなら、こちらも同時に反映します。

```text
MS2DecR-main/BioC/MS2DecR-bioc/R/performICA.R
MS2DecR-main/BioC/MS2DecR-bioc/DESCRIPTION
MS2DecR-main/BioC/MS2DecR-bioc/NAMESPACE
MS2DecR-main/BioC/MS2DecR-bioc/man/performICA.Rd
```

アップロードしなくてよいもの:

```text
*.tar.gz
*.Rcheck/
dev/diapasef_deconv_test/
dev/diapasef_deconv_summary/
```

`dev/` 以下の diaPASEF 検証ファイルは、論文・検証用に残すなら別途コミットしてもよいですが、
MS2DecR 本体の ICA 修正として必須なのは上記の `performICA` 関連ファイルです。

## MCR-ALS のモビログラム判定を緩めた追加検証

5成分 overcomplete MCR-ALS の後にピークピッキングする際、当初は RT と IM の両方に
ほぼ同じ判定を使っていました。

具体的には、`matchedfilter` 後の頂点が元の頂点に近いことと、局所ピーク数が少ないことを
同時に要求していました。しかし、IM 側のモビログラムは点数が少なく、細かい揺れが
局所ピークとして数えられやすいため、見た目には比較的きれいな成分でも落ちることがありました。

直近の実行では、`ALS2` はモビログラムの頂点位置自体は合っていましたが、
通常判定では `m_n_peaks = 7` となり `m_pass = FALSE` でした。
これは「ピーク位置が悪い」のではなく、「局所ピーク数の基準がモビログラムには厳しい」
という挙動です。

そこで、IM 側だけ次のような `mobility-shape` 判定を追加しました。

```text
平滑幅を 7 にする
局所ピーク数の閾値を最大強度の 50% に上げる
局所ピーク数は 6 個まで許容する
matchedfilter 頂点一致は診断値として残すが、合否条件にはしない
```

この条件では、`ALS2` は `m_pass = TRUE` になり、候補として残りました。
結果は次の通りです。

```text
method                                  S1_best  S8_best  mean
B4_mobility_shape_peak_select_only       0.756    0.767   0.761
B5_mobility_shape_fixed_CM_A_only        0.810    0.787   0.799
MS-DIAL-like window baseline reference   0.771    0.968   0.870
```

この結果から、MCR-ALS 側では IM ピッキング基準はまだ厳しすぎたと考えられます。
特に、`mobility-shape` で選んだ C/M を固定して window + baseline 付きで `A` だけ再推定すると、
S1 側は 0.810 まで改善しました。

ただし、S8 側は MS-DIAL-like の 0.968 には届いていません。
したがって、現時点の結論は次の通りです。

```text
MCR-ALS は、モビログラム判定を緩めると改善する。
成分2相当を落としていた元の基準は厳しすぎた。
ただし、Fig.6 データでの最良結果はまだ MS-DIAL-like window + baseline 共通A推定。
```

追加で出力した主なファイル:

```text
MS2DecR-main/dev/diapasef_deconv_summary/mcr_als_overcomplete_mobility_shape_component_table.csv
MS2DecR-main/dev/diapasef_deconv_summary/mcr_als_overcomplete_mobility_shape_peak_selected_deconvolution_summary.png
MS2DecR-main/dev/diapasef_deconv_summary/mcr_als_overcomplete_mobility_shape_Aonly_deconvolution_summary.png
```

## `joint_mcr_als` の MS2DecR 本体との整合性修正

`test_diapasef_deconv.R` の `joint_mcr_als` を MS2DecR 本体の `optimizeALS` と
揃えるために、以下の3点を修正しました。

### 修正1: A 正規化

MS2DecR 本体の `optimizeALS` は A の各行を正規化しています。
今回の `joint_mcr_als` は C/M を正規化していましたが、A 正規化に変更しました。

### 修正2: 最小誤差 iteration 選択

MS2DecR 本体の `optimizeALS` は全 iteration の中から最小誤差の結果を返します。
今回の `joint_mcr_als` も同様に変更しました。
ALS が後半で発散するケースに対して安全です。

### 修正3: maxiter デフォルト 500

maxiter=100 では収束が不十分でした。テストの結果:

```text
maxiter   error     S8 purity  S1 purity
100       18596.0   0.793      0.719
200       18557.0   0.807      0.719
500       18505.2   0.814      0.719
1000      18505.0   0.814      0.719
```

500 回でほぼ収束するため、デフォルトを 500 に変更しました。

## IM 側のピークピッキング基準の調整

モビログラムのピークはクロマトグラムに比べてブロードなため、
`max_local_peaks` を 2 → 3 に変更しました。

5成分 MCR-ALS (k=5) のピークピッキング結果:

```text
成分    RT c_pass (n_peaks)  IM m_pass (n_peaks)  両方通過
ALS1    TRUE   (1)           FALSE  (7)           ×
ALS2    FALSE  (3)           FALSE  (14)          ×
ALS3    FALSE  (13)          TRUE   (1)           ×
ALS4    TRUE   (1)           TRUE   (3)           ○
ALS5    TRUE   (2)           TRUE   (1)           ○
```

- RT 通過: ALS1, ALS4, ALS5
- IM 通過: ALS3, ALS4, ALS5
- 両方通過: ALS4, ALS5

ALS4 は IM の `max_local_peaks=2` では落ちていましたが（n_peaks=3）、
`max_local_peaks=3` にすることで通過するようになりました。
matchedfilter の wid を変えても n_peaks は変わりませんでした。

## 5成分 MCR-ALS の最終結果

修正後（A正規化, 最小誤差選択, maxiter=500, IM max_local_peaks=3）:

```text
成分   site          purity
ALS4   S8_phospho    0.814
ALS5   S1_phospho    0.719
```

各成分の物理的位置:

```text
ALS4: RT 17.9 min, IM 0.96 → S8 の期待位置（RT 18.27, IM 0.9635）に近い
ALS5: RT 18.3 min, IM 0.94 → S1 の期待位置（RT 18.46, IM 0.9410）に近い
```

## Lambda 調整の結果

lambda を 1e-8 から 10.0 まで変えてテストしました。

```text
lambda    selected  S1      S8      mean
1e-8      ALS4,5    0.719   0.814   0.766
1e-6      ALS4,5    0.719   0.814   0.766
1e-4      ALS4,5    0.718   0.814   0.766
1e-3      ALS4,5    0.717   0.814   0.766
1e-2      ALS4,5    0.719   0.781   0.750  ← 悪化
0.1       ALS5のみ  0.717   ---     ---    ← S8成分が落ちる
1.0+                                       ← 正則化が強すぎ
```

lambda=1e-6（現状）がほぼ最適で、変えても改善しませんでした。

## Seed 依存性の確認

fastICA の初期化が結果に影響するかを 20 seed でテストしました。
成分の番号は seed ごとに変わりますが、最終的な purity はほぼ同じ値
（S1≈0.719, S8≈0.814）に収束しました。
一部の seed は局所解に落ちますが、最良でも改善しません。

## 10成分 (k=10) の結果

成分数を増やしてテストしました。

```text
成分    RT c_pass  IM m_pass  purity       site
ALS5    TRUE       FALSE(5)   0.875        S1  ← IM で落ちる
ALS6    TRUE       TRUE       0.757        S8  ← 選ばれた
ALS7    TRUE       TRUE       0.787        S1  ← 選ばれた
ALS10   TRUE       FALSE(8)   0.866        S8  ← IM で落ちる
```

S1=0.787, S8=0.757 で、5成分の S8=0.814 より S8 側が悪化しました。
IM のピークピッキングで落ちた ALS5 (S1=0.875) や ALS10 (S8=0.866) を
拾えれば改善しますが、IM 側の基準をさらに緩めることになります。

## ICA 単独の評価

ALS の繰り返し最適化前の ICA 結果を評価しました。

クロマトグラム ICA (X のみ) から A を推定:

```text
成分    site    purity
ICA1    S1      0.960
ICA5    S1      1.000
ICA4    S8      0.718
ICA2    S8      0.653
```

ICA C + ICA M をペアリングして joint A を推定:

```text
成分    site    purity
ICA1    S1      0.834
ICA5    S1      0.994
ICA2    S8      0.646
ICA4    S8      0.580
```

ICA 段階では S1 の purity が非常に高い（0.96〜1.00）のに対し、
ALS の繰り返し後は 0.719 まで下がっています。
ALS は再構成誤差を最小化する方向に C/M を動かすため、
ICA で得られていた成分の独立性が崩れている可能性があります。

ただし、この purity 評価は fragment の site ラベル割合を見ているだけで、
正解の MS/MS スペクトルとの cosine similarity は計算していません。
ICA5 (purity 1.000) は少数の S1 fragment だけで構成されているため、
実際の S1 MS/MS を正しく再現しているかは別途検証が必要です。

## 現時点の全手法比較

```text
method                                    S1_best  S8_best  mean
ICA only (C only, no ALS)                 1.000    0.718    0.859  (*)
ICA only (joint C+M, no ALS)              0.994    0.646    0.820  (*)
MCR-ALS k=5 (修正版, maxiter=500)         0.719    0.814    0.766
MCR-ALS k=10                              0.787    0.757    0.772
MS-DIAL-like window + baseline 共通A推定   0.771    0.968    0.870

(*) purity は site ラベル割合のみの評価。正解MS/MSとの比較は未実施。
```
