# 調査資料の索引

調査本文・解析の保存仕様と索引を `docs/research/` にまとめています。
調査コードは `investigations/<テーマ>/`、図・ROOT・集計表・実行ログは `investigations/result/<テーマ>/` に分けています。
通常の解析結果とホーム配下の既存結果は従来の場所にあり、下表から参照できます。
今後の保存先は [AGENTS.md](../../AGENTS.md) に従います。

## 調査資料と既存結果

| 調査本文 | 図・データ・実行記録 |
|---|---|
| [data/MC 再評価](reassessment.md) | [既存の保存先](/home/murata_t/scecal_reassessment_20260929/) |
| [raw → calibration の調査](raw-calibration.md) | [既存の保存先](../../investigations/result/raw_calibration/100GeV/) |
| [impact study 後の追加調査](response-diagnostics.md) | [既存の保存先](../../investigations/result/response_diagnostics/threshold/production_v1/) |
| [e_ssa/threshold の data/MC 比較](e-ssa-threshold.md) | [既存の保存先](../../investigations/result/e_ssa_threshold_review_20260930/) |
| [HG–LG tail 生成結果の確認](hg-lg-tail.md) | [既存の保存先](../../investigations/result/hg_lg_tail_review_20260930/) |
| [全 muon run の MIP 規格化・校正](muon-summary.md) | [既存の保存先](../../analysis/result/channel_response/mip_closure/all_muon_runs_20260930/) |
| [ビーム入射分布探索の保存仕様](beam-input-scan.md) | 条件別ROOT・summary・図・実行記録の予定配置（未実装） |
| channel response の補助調査 | [図・集計・実行記録](../../investigations/result/channel_response_review_20260930/) |

## 実行方法

- [raw → calibration の調査](../../investigations/raw_calibration/README.md)
- [impact study 後の追加調査](../../investigations/response_diagnostics/README.md)
- [e_ssa/threshold の調査](e-ssa-threshold.md)
- [muon/MIP 解析コード](../../analysis/channel_response/mip_closure/)
- [muon 詳細報告書の生成コード](/home/murata_t/scecal_muon_summary_20261001/build_summary.py)

通常の解析コードは `analysis/`、再利用する調査コードは `investigations/` にあります。
調査段階では代表 energy を使用し、全 energy の解析は明示依頼された場合だけ行います。
`response_diagnostics/run.sh` は過去の複数 energy 調査を再現する旧ドライバであり、新しい調査の既定ドライバとして使いません。

## 調査コードと結果の分離（2026-10-02）

5テーマの図・ROOT・集計表・実行記録を `investigations/result/<テーマ>/` に移動しました。
コード・設定・実行方法の説明は `investigations/<テーマ>/` に残し、出力先・読込先と資料のリンクを更新しています。
結果内の `inputs/`、再現用ソースの控え、チェックサム、過去のログは内容を変更せず保持しました。記録中の旧パスは実行当時のものです。
通常の解析結果 `analysis/result/` とホーム配下の結果は移動していません。

| 移動元（investigations/ 基準） | 現在の結果保存先 |
|---|---|
| `channel_response_review_20260930/` の結果ファイルと `figures/` | [result/channel_response_review_20260930/](../../investigations/result/channel_response_review_20260930/) |
| `e_ssa_threshold_review_20260930/result/` と直下の集計TSV | [result/e_ssa_threshold_review_20260930/](../../investigations/result/e_ssa_threshold_review_20260930/) |
| `hg_lg_tail_review_20260930/completed_summary.tsv` | [result/hg_lg_tail_review_20260930/](../../investigations/result/hg_lg_tail_review_20260930/) |
| `raw_calibration/result/` | [result/raw_calibration/](../../investigations/result/raw_calibration/) |
| `response_diagnostics/result/` | [result/response_diagnostics/](../../investigations/result/response_diagnostics/) |

raw_calibration と response_diagnostics の結果内にあった `investigation.md` の互換リンクも、新しい結果保存先へ移動しています。

## 資料の集約（2026-10-02）

この段階では本文と索引だけを移動しました。数値・結論は保持し、資料内の参照先を更新しています。
資料集約時の旧ファイル名には移動先への互換リンクを設けました。結果内の2件は、上記の結果分離でリンク自体も移動しています。資料内の相対リンクは移動先を基準にしているため、閲覧時はこの索引から本文を開いてください。
過去のログやチェックサム記録は実行当時のままです。
HG–LG tail 資料にある旧 `threshold/default/50GeV/` の4つの参照先は、整理前から見つからないため、本文に注記しています。

| 資料集約時の旧ファイル | 現在の本文 |
|---|---|
| `investigations/README.md` | [README.md](README.md) |
| `investigations/e_ssa_threshold_review_20260930/README.md` | [e-ssa-threshold.md](e-ssa-threshold.md) |
| `investigations/hg_lg_tail_review_20260930/README.md` | [hg-lg-tail.md](hg-lg-tail.md) |
| `investigations/raw_calibration/result/100GeV/investigation.md` | [raw-calibration.md](raw-calibration.md) |
| `investigations/response_diagnostics/result/threshold/production_v1/investigation.md` | [response-diagnostics.md](response-diagnostics.md) |
| `/home/murata_t/scecal_muon_summary_20261001/detailed_summary_ja_20261001.txt` | [muon-summary.md](muon-summary.md) |
| `/home/murata_t/scecal_reassessment_20260929/reassessment_ja.md` | [reassessment.md](reassessment.md) |

## 以前の整理履歴

2026-09-30にコードと結果を `investigations/` へ集約。既存ROOT・PNG・TSVを再生成せず移動した。
inputs内のchecksumや過去のログに記録されたパスは実行時の記録として保持している。
旧配置からの対応は次のとおり。

| 旧パス（ScECAL_BeamTestからの相対パス） | 2026-09-30時点の移動先 |
|---|---|
| response_diagnostics/ | investigations/response_diagnostics/ |
| result/response_diagnostics/ | investigations/response_diagnostics/result/ |
| raw_calibration_investigation/ | investigations/raw_calibration/ |
| result/raw_calibration_investigation/ | investigations/raw_calibration/result/ |
| channel_response_review_20260930/ | investigations/channel_response_review_20260930/ |
