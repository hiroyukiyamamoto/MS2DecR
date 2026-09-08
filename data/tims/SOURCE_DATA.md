# Source data

## 元の生データ

- データセット: **PXD033904**
- PRIDE: https://www.ebi.ac.uk/pride/archive/projects/PXD033904
- 生データ一覧: https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/10/PXD033904/
- **60分グラジエントのダウンロードURL:** https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/10/PXD033904/60min.zip
- FTPで取得する場合: ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/10/PXD033904/60min.zip
- 配布場所の確認元: https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD033904-3&test=no

公開ディレクトリに60min.zipがあることを2026-09-08に確認しました。この公開ZIPの再ダウンロードは実施せず、既に展開されていた次のrunを再抽出に使用しています。

```text
20220302_tims1_nElute_8cm_DOl_Phospho_60min_rep1_Slot1-95_1_1835.d
```

展開後、このrun名の.dフォルダを探して `FIG6_D_PATH` に指定してください。`fig6_run.d` はコマンド例の説明用名称です。

.dはフォルダ全体を扱い、少なくともanalysis.tdfとanalysis.tdf_binを保持します。本デモ作成に使用したファイルサイズは、それぞれ101,457,920 bytes、3,391,397,888 bytesです。生データ本体はGitHubへ同梱していません。

## 論文・関連データ

- 元のリン酸化プロテオミクス研究: https://doi.org/10.1002/pmic.202200032
- diaTracer論文（Fig.6の例）: https://doi.org/10.1038/s41467-024-55448-8
- diaTracer論文のData availability: https://pmc.ncbi.nlm.nih.gov/articles/PMC11696033/
- diaTracer変換mzML・FragPipe結果: https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?accession=MSV000094803
- 関連データFTP: ftp://massive.ucsd.edu/v07/MSV000094803/

論文のData availabilityでは、リン酸化プロテオミクスの生データをPXD033904、変換mzML・FragPipe結果をMSV000094803として区別しています。本手順は.dを直接読むため、変換mzMLは必須ではありません。
