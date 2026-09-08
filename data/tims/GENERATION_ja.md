# 生データからMS2DecR用デモを生成するプロセス

## 処理の流れ

```text
PRIDE PXD033904 / 60min.zip
  → 対象runのBruker .d
  → opentimsr + MS-DIAL同梱Bruker timsdata.dll
  → プリカーサーm/zを含むDIAウィンドウ・RT・IM・m/zで抽出
  → *_raw_peaks.tsv
  → RT × IM × m/zビニング → *_cube_binned.tsv
  → raw_peaksをframeごとにIM積算 → 2次元MS2行列 → RDS
```

DLLは `opentimsr::setup_bruker_so()` に渡し、TOF→m/zとscan→逆イオン移動度（1/K0）の変換に使います。検証に使用したのはMS-DIAL 5.5の `MSDIAL.console.v5.5.260319-windows-net48/lib/Bruker/timsdata.dll` です。本データではopen-source変換でIM座標のずれが見られたため、同梱スクリプトはDLLが見つからない場合に停止します。

## 抽出条件

対象ペプチドはSPSPPDGSPAATPEIRのS1/S8リン酸化位置異性体、電荷2です。

| 設定 | 既定値 |
|---|---|
| FIG6_PRECURSOR_MZ | 829.8746 |
| FIG6_FRAGMENT_MZ_LOWER / UPPER | 150 / 1500 |
| FIG6_RT_BIN_WIDTH_MIN | 0.01分（0.6秒） |
| FIG6_IM_BIN_WIDTH | 0.001 |
| FIG6_MZ_BIN_WIDTH | 0.01 Da |
| FIG6_FRAME_CHUNK_SIZE | 80 |

| 領域 | RT（分） | IM（1/K0） |
|---|---|---|
| both_deconv | 17.60–18.70 | 0.930–0.970 |
| S1_box | 18.40–18.50 | 0.932–0.948 |
| S8_box | 17.75–17.85 | 0.953–0.968 |

S1_box/S8_boxは比較用の範囲抽出です。純粋な単一成分や正解スペクトルであることは保証しません。

.d内のFrames、DiaFrameMsMsInfo、DiaFrameMsMsWindowsを読み、対象m/zを含むDIAウィンドウを選びます。本データはWindowGroup=12、scan=299–918、プリカーサー分離範囲824.35–847.29です。フラグメントはm/z 150–1500の正の信号を保持し、参照ライブラリのフラグメントによる事前選択は行いません。

現実装は該当WindowGroupのframeを選び、該当scan範囲の和集合で抽出します。本runでは該当ウィンドウは1つです。複数WindowGroupが該当する別runへの適用時は、frameごとのscan範囲の対応を追加確認してください。

## 再生成

[READMEのPowerShellコマンド](README.md#生データから再生成)をdata/timsから実行します。入力.dとDLLの場所は利用者が指定します。outputには検証済みデータがあるため、再生成先の既定の例はregeneratedです。

抽出スクリプトは各領域のraw_peaks・cube_binned・summaryと共通のcube_windows・precursor_containing_windows・cube_summaryを出力します。その後build_demo.RがRDS、検証表、確認図とgzip圧縮TSVを作ります。圧縮前TSVも残すため、再生成後の容量は配布版より大きくなります。

領域を変更する場合は、抽出時に次の環境変数が利用できます。

```powershell
$env:FIG6_CUBE_WINDOW = 'custom,17.6,18.7,0.930,0.970'
```

build_demo.Rは同梱の3領域・既定の範囲用です。領域名や条件を変更する場合は、入力名・検証条件・metadataも対応して変更してください。

## 出力の意味

- raw_peaks: frame、scan、tof、intensity、mz、inv_ion_mobility、retention_time（秒）、rt_min（分）と各ビン座標。
- cube_binned: rt_bin_min、im_bin、mz_bin、合計intensity、n_peaks、強度加重平均のmz/RT/IM。
- ビン座標: floor(値 / ビン幅) × ビン幅で求めた下端。信号のあるビンのみ記録するlong table。
- RDS: list(mz, Y, rt, frame, metadata)。Yの行は測定frame、列はm/zビン。rtは秒。
- 2次元化: raw_peaksから同一frame・m/zビン内のIM方向の強度を合計。RTは元のframe時刻を保ち、追加の時間ビニングはしない。
- Yは信号のあるframeとm/zビンを軸にし、組み合わせに信号がなければ0。全く信号のないframeは復元しない。
- 正規化・平滑化は実施していない。

既存MS2DecRのgenerate_ms2_matrix()が返すmz/Y/rtと同じ軸構成ですが、MSnbaseのビニングと完全に同じ前処理ではありません。3次元の解析にはcubeまたはraw_peaksを利用し、IMを積算したYとは区別してください。

## 検証

build_demo.Rは全3領域についてピーク数、RT/IM/m/z範囲、強度の有限性と正値、raw→cubeの総強度保存を確認します。さらに2次元Yの次元、時間軸、非負性、総強度保存を確認します。結果はoutput/validation.tsvに保存します。

MS2DecRによる成分分離の最適化や性能評価は、この生成処理の検証対象には含めません。
