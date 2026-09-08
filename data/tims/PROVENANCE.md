# Generation provenance

## 入力と環境

- 生成・検証日: 2026-09-08
- 生データ: PRIDE PXD033904（配布元と対象runは[SOURCE_DATA.md](SOURCE_DATA.md)）
- OS: Windows
- R: 4.5.1 (ucrt)
- opentimsr: 1.2.0
- data.table: 1.18.4
- 変換: opentimsr::setup_bruker_so()
- DLLの出所: MS-DIAL 5.5、配布物 `MSDIAL.console.v5.5.260319-windows-net48` の `lib/Bruker/timsdata.dll`
- DLL SHA-256: `8b533c66a5c2961dfc2a4a91a589d83f025fe99e086e7cd22a820da74a4c50cf`

DLLそのもの、元.d、個人環境の絶対パスを含むログは同梱していません。公開用には本ファイルと数値検証表を保存しています。

## 実行結果

同梱scripts/extract_fig6_3d_cube_opentimsr.Rを既定の抽出条件と上記DLLで実行し、scripts/build_demo.Rで2次元化と検証を実施しました。

| 領域 | ピーク行数 | cube行数 | 総強度 |
|---|---:|---:|---:|
| both_deconv | 197659 | 197202 | 14903345 |
| S1_box | 8948 | 8905 | 736252 |
| S8_box | 6366 | 6366 | 516612 |

2次元デモは **48 frames × 36488 m/z bins** です。全領域で抽出範囲の検証とraw→cubeの総強度保存が通過し、both_deconvについてraw→Yの総強度保存も通過しました。詳細は[output/validation.tsv](output/validation.tsv)を参照してください。

圧縮TSVは元のTSVと展開後のバイト列が同じであることも確認しています。RDSと圧縮TSVの読込み、圧縮TSVからの再構築を検証しました。

この結果はデータ抽出・変換の検証です。MS2DecRによる成分分離性能や、別runへの一般化は検証していません。
