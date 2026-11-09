# Response diagnostics

`hl_intercept/threshold/production_v1` に残るdata/MC差を調べる読み取り専用の診断。
C++/ROOT + shell、図はPNGのみ。元のROOT・校正・digitizationは変更しない。

```bash
cd /home/murata_t/ScECAL_BeamTest/investigations/response_diagnostics
bash run.sh
# 任意の入出力ディレクトリ
bash run.sh INPUT_LABEL_DIRECTORY OUTPUT_DIRECTORY
```

直接実行する。ジョブ投入は行わない。既定出力は
`investigations/result/response_diagnostics/threshold/production_v1`。
再実行では診断出力だけを更新する。入力は完了済みimpact studyを使用し、
`event_histograms.root` のdataイベント総energyキャッシュを必要とする。

- `SurveyImpact`: 全19energyの全イベントを使い、低energy割合・run別集計・主ピークを調べる。
  低energyは `signed E < 0.1 * median(MC nominal signed E)`。
  core統計は `0.5*median < signed E < 1.5*median`。
  Gaussianは主ピーク付近±1.5 sigmaの反復fit。非Gaussian性があるためfit誤差だけを
  系統誤差や正式な分解能として扱わない。fit status・chi2/NDF・core medianも保存する。
- `AuditMIP`: 現在のExtract.cxx/Calibration.cxxの採用・代替規則で定数を比較する。
  元ファイルは同じ。calibration側の既存挙動は`../../analysis/impact_studies/StudyCommon.hh`で再現。
- `InspectResponse`: 0.5,1,5,30,100,200,250 GeVを詳しく調査する。
  PS dataは各ファイル最大100000イベント（今回の3energyでは全イベント）、
  SPS dataは各ファイル2000イベント。MCは全入力ファイル合計約20000イベント。
  各ファイルを10区間に分け、その区間中央の連続イベントを取る決定論的なサンプル。
  MCの総energy・切片0 energyを元のpaired_eventsと照合する。
  HG/LGの使用分岐は保存EH/EL/Eの厳密な一致から推定し、曖昧hit数と元の切片0総和との差を検証する。
- `InspectTopology`: 1,5,100,200 GeVでCoGと低energyイベントのhit数・layer数を調べる。
  dataは各ファイル2000、MCは合計約20000。ヒストグラムは選別イベント数で正規化して共通色軸で描画。
- `CheckTruthResponse`: 1 GeVの先頭10000イベント、100 GeVの先頭ファイル1000イベントについて
  同じbasenameのMC_Truthとrawを比較。entry数と全raw hitのCellID対応を確認。
  truthへの対応は生成・digitizationがentry順を保つことを前提とする。
  `truth_response.tsv` category = 2*inner10um + (digi/calib MIP不一致)。
  truth_bandはMeVで `[0.5,1),[1,2),[2,5),[5,10),[10,30),[30,100)`。
  HGはraw HG<2400だけを集計、LGはpedestal差引き後にHL傾きでHG相当へ変換。
  切片を追加する前の旧MCを対象とする検査。新しいMCではLGの切片分を分離して解釈すること。
- `PlotResponse`: 保存済み診断ROOT/TSVから比較PNGを描画する。
- `verify_results.sh`: 元MCとの一致、signed=positive+negative、channel集計との一致を検証。

詳細比較の`shower`は `0.5*MC signed median < sum(positive hit E) < 1.5*median`。
粒子同定ではなく、低energy成分を分離する診断用選択である。
`shower_cog20`はさらに、負hitとbad channelを除いた全hitから求めたlayer9・10のCoGについて
それぞれ `|x|<20 mm && |y|<20 mm` を要求する。
エネルギー和の`good`は負hit・comparisonのbad-channel maskを除外。
`central`はさらにchannel中心の `|x|,|y|<22.5 mm`。peripheralはgood-central。
CoGはcentral cutの前に計算する。MC0でもイベント・hit集合は元MCで固定する。
maskは物理CellIDに適用。今回の全sampleではmemory番号付きIDと物理IDに差はなかった。
条件付きenergyの図はdata/MC両方30イベント以上のNhit binだけを表示し、誤差はSEM。
同じNhit・同じCoG条件でも、layer分布・厳密な位置分布まで一致した比較ではない。

`event_metrics.tsv`は各選択のイベント平均。`channels.tsv`、`layers.tsv`、`hit_energy_bands.tsv`は
shower_cog20の結果。hit energy bandは `[0,.1),[.1,.3),[.3,1),[1,3),[3,10),[10,30),[30,100),[100,inf)` MeV。
`good_energy/`、`central_energy/`の図も同じ選択後の比較。
`sampled_files.tsv`に実際のファイルとsample数、`validation.tsv`に元MCとの一致精度を保存する。
今回の結論・数値・次の切り分けは [調査本文](../../docs/research/response-diagnostics.md) を参照。
