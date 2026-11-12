# Beam quantitative evaluation

保存済み `Calib_Hit` と `SSA_Hit` から、入口位置とhit levelのdata–MC一致度を定量評価する。
コードは独立したディレクトリに置き、`../ssa_shower_shape/Common.hh` の入力検証・配置・入口軸・同一位置統合を再利用する。
生成、digitization、再校正、SSA再生成、ビーム条件の自動探索は行わない。保存済みの条件を指定して評価・順位付けする。

## 実行

C++17、ROOT、Python 3とNumPy/uproot/awkwardが必要。SciPy・matplotlib・PyROOTは不要。
既定は100 GeV、各sample最大4ファイル・1500イベント、500 bootstrap。
比較全体の既定名は `--label beam_scan`、現在のMC生成条件の既定名は `--condition default`。
今後のMC候補は同じ `beam_scan` の下でcondition名を分け、候補ごとのROOT入力先を指定する。
label/conditionは保存・比較用の名前で、MCの生成条件や入力先を変更する指定ではない。
イベント数制限はphase選別後、入口重心が定義できるかの選別前に適用する。

```bash
cd ~/ScECAL_BeamTest/analysis/beam_quantitative
# 入力・保存予定の確認。ビルド・出力作成・ジョブ投入なし
bash run.sh --label beam_scan --condition default --dry-run
# ローカル実行。MCの生成条件は注記として指定する
bash run.sh --label beam_scan --condition default --material-label unspecified
# Slurm（energyごとに1ジョブ、2 CPU・8 GB・4時間）
bash submit.sh --energies 100 --label beam_scan --condition default --dry-run
bash submit.sh --energies 100 --label beam_scan --condition default
# 候補MCの別ディレクトリを指定。同じlabelではdata参照を再利用する
bash run.sh --label beam_scan --condition sx0p50_sy0p75 \
  --mc-calib-dir /path/to/calib/100GeV \
  --mc-ssa-dir /path/to/ssaCalib/100GeV \
  --mc-truth-dir /path/to/generate/100GeV
# 全入力・統計増加。同じlabelの参照条件は変えない
bash run.sh --label full_statistics --max-files 0 --max-events 0 --bootstrap 1000
# 保存済み評価から再描画
bash run.sh --label beam_scan --condition default --plots-only
```

`--data-calib-dir/--data-ssa-dir` も指定可能。各指定先はenergy末端ディレクトリ。
`--truth auto|required|off` は既存解析と同じ。明示したtruthの欠損・部分欠損は停止する。
既定dataは `ECAL_data/analysed/2023/{sps,ps}/{calib,ssaCalib}/e-/<energy>GeV`。
既定MCは `Result_MC/{calib,ssaCalib}/e-/e_ssa/threshold/<energy>GeV`。
生成truthは `Result_MC/generate/e-/e_ssa/<energy>GeV`。
5 GeV以下のdataはps、それより上はsps。入力を別productionへ自動置換しない。
`--material-label` と `--production-note` は由来の記録のみで、実際のgeometryを変更・証明しない。

## 選別と整合性

- 接続されたlayer 0–29のみ。chip5のchannel30–35、その他不正なCellIDは診断に数えて除外。
- calib/SSAの同名ファイル・entry数・Run/Event/Time/Detectorを確認。
- 共通の元CellIDの符号付きエネルギー保存を `1e-7 + 1e-7*sum(abs(pre))` MeVで確認。SSAだけの親IDは停止。
- 欠落親は件数と符号付きエネルギーを記録し、理由を自動断定しない。
- calib座標を既存EBUdecodeと照合。有限値・vector長・layer Z・ssaTagを検証。
- truthは同名・同entry数、MCのRun=0, Event_Time=Event_Num=entry、電子PDG・energy・方向を確認。
- 入口重心はSSA前の正の信号から、Xを奇数層、Yを偶数層で計算。既定layer0–3。
- X/Y双方が定義できるイベントで全hit指標を比較。軸未定義・非接続hit・欠落親はJSON cutflowに保存。
- SSAは既定 `--ssa-selection split`（tag1）。`all` はtag0/1両方。
- 同じイベントの完全一致 `(layer,x,y)` は1点に統合し、エネルギーは符号付き和を取る。
- countはゼロ・負信号を含む固有位置数、positive_energyは統合後の `max(0,sum)`。
- 追加のenergy窓・CoG cut・dead-channel mask・平滑化は適用しない。

## 定量指標と暫定許容差

探索開始前に固定する100 GeV向けの初期案で、既存測定から導いた許容差ではない。
必要なら `--tolerances-json` に同じキーを指定する。未知キー・非有限・ゼロ/負許容値は停止。

| key | 指標 | 許容差 |
|---|---|---|
| mean | mean_MC − mean_data | ±0.5 mm |
| rms_relative | RMS_MC/RMS_data − 1 | ±0.05 |
| q90_relative | q90_MC/q90_data − 1 | ±0.10 |
| tail_fraction | tail_MC − tail_data | ±0.01（1 percentage point） |
| correlation | rho_MC − rho_data | ±0.05 |
| tv | 1/2 sum(abs(p_MC − p_data)) | 0.05 |
| hits_relative | mean_hits_MC/mean_hits_data − 1 | ±0.10 |

- RMSは平均まわりの標準偏差。q90は**各sample自身の平均**からの絶対距離の、重み付き離散90%分位点。
- tailは**dataの平均位置**からの絶対距離が30 mm以上の割合。MCにも同じdata中心を用いる。
- 入口重心ではX/Yを各イベント1回数える。X–Y相関、2D TVも計算する。
- 層別hitの1Dは奇数層X、偶数層Yの細かく測れる方向。count/positive_energyを別々に評価。
- hitのfine軸はストリップ中心の実座標を1e-6 mmで丸めて集計し、図の5.3 mm表示binとは独立に評価。
- 2Dは5.3 mm幅、44×44 bin、±116.6 mm。範囲外は1つのflowカテゴリーとしてTVに含める。
- TVのshapeは全カテゴリーで合計1に規格化。hit数は採用イベント数で割る。
- TSVのTV行はdifferenceが距離本体で、data_value/mc_valueは個別sampleの物理量ではない。
- 図の1D表示・2D色表示は範囲内。flowはROOTに保持し、TV評価には含める。図の平均値で合否を決めない。
- 空分布・ゼロ幅の相対値・定義不能な相関は未判定とする。

## Bootstrapと判定

data・MCを独立に、**イベント全体**を復元抽出するordinary bootstrap。
全指標・全層で同じ再抽出重みを共有し、イベント内相関を保つ。
再抽出ごとに平均・包含幅・data中心のtail・TVを計算し直す。
95%区間はpercentile法の2.5/97.5%点による近似区間。厳密な保証や全項目同時95%保証ではない。
TVのような非負指標は有限統計による偏りがあり、低統計では合格しにくい。
bootstrapの信頼区間は統計的不確かさのみで、geometry・校正などの系統誤差を含まない。

- PASS: 区間全体が許容範囲内。
- FAIL: 区間全体が許容範囲外。
- UNDETERMINED: 境界をまたぐ、定義不能な再抽出がある、または採用イベント数が `--min-events` 未満。
- 既定のmin-eventsは各sample200。bootstrapは最低50、通常は500以上を使用する。
- scoreは `max(abs(CI_low),abs(CI_high))/tolerance`。条件全体は合格判定項目の最大score。
- 合格項目は入口重心と `--gate-layers 0,1,2,3` の検出器座標でのhit指標。
- 全30層のlab評価と入口重心補正後のmean/RMSは別の診断。alignedをビーム幅の最適化項目に含めない。
- 条件全体はFAILが1つでもあればFAIL、それ以外で未判定があればUNDETERMINED、全項目PASSでPASS。
- 同じ評価定義・上流物質注記の条件のみrankingに含め、PASS→UNDETERMINED→FAIL、各群内でscore順に並べる。
- 入力が点入射でも有限の再構成幅がある。合格は固定MC条件下での有効なビーム候補で、真のビーム分布を証明しない。

## 探索と独立検証

既定 `--phase scan` はfilename:entryの安定hashによる半分のイベントを使用する。
`--phase validation` は残り半分で、同じ元データの両phaseは重複しない。
ファイル制限は全一覧を等分した中央ファイル、イベントはeligible poolを等分したstratumから乱数選択。
同じlabel/phaseではdata入力・サンプリング・軸・抽出ソースを固定し、候補ごとに同じdataを使う。
現在のMC結果は `beam_scan` の `default` に保存し、別候補も同じlabelのrankingに追加する。
名前を変更しても既存の評価値・判定は変わらず、defaultは基準条件を表し、合格を意味しない。
`--phase all` は分割せず全eligibleイベントから抽出するため、独立検証とは呼ばない。
**phase変更はMCを別seedで生成する処理ではない。** 最終検証には独立seedで生成したMCを明示して入力する。

## 保存先

```text
analysis/result/beam_quantitative/<label>/
  root/data/<phase>/reference/<energy>GeV/events.root
  root/mc/<phase>/<condition>/<energy>GeV/{events,histograms}.root
  root/comparison/<phase>/<condition>/<energy>GeV/
    data_histograms.root, metrics.root, comparison.root
  plots/<phase>/<condition>/<energy>GeV/
    entrance_positions[_log].png, entrance_xy.png, assessment.png
    profiles_<lab|aligned>_<count|positive_energy>[_log].png
    hit_map_<data|mc>_<even|odd>_<count|positive_energy>.png
  summary/<phase>/<energy>GeV/
    <condition>.metrics.tsv, <condition>.json, ranking.tsv
  inputs/<phase>/<energy>GeV/
    reference/{configuration.json,data.tsv,complete.json}
    <condition>/{configuration.json,mc.tsv,run.log,sources/}
analysis/log/beam_quantitative/<label>/<phase>/<condition>/<energy>GeV/
  build.log, data.log, mc.log, plot.log
```

ROOT領域のcondition直下にはGeVディレクトリのみ、各GeV内にはROOTのみ。
summary・図・設定・SlurmログをROOT領域に混在させない。CoG cutを使わないので形式上のCoG階層は追加しない。
metrics.rootは全定量値・区間・判定のTTree、comparison.rootは再編集可能なcanvas。
metrics TSV/ROOTの `value_unit` と `difference_unit` はsample値と差の単位を明記する。幅そのものはmm、幅の相対差はrelative。
events.rootはbootstrapの由来を示すイベント別入力であり、保持する。
SSAの `profiles_<lab|aligned>_<count|positive_energy>[_log]` は、線形・対数表示とも5.3 mm幅で表示する。
labはストリップ中心間隔に合わせた42 bin、alignedは入口重心を引いた連続座標を同じ幅の85 binに集計する。
これにより、1 mm binがストリップ中心の間に作っていた周期的な空binを解消する。
数値評価は従来どおり実座標から計算し、TVの2D binは変更しない。
保存済み結果は通常の `--plots-only` でイベントROOTから表示binを更新する。
表示binの設定・ソースchecksumは `inputs/.../<condition>/plot_sources/` とrun.logに保存する。
凡例はcanvas上部の専用欄に置く。対応する2D mapは層ごとにdata/MCで共通の色範囲。
condition名 `reference` は共有data用の予約名。
同じlabel/phase/condition/energyは上書き拒否。失敗途中はcompleteとせず、別condition/labelを使う。
shared data抽出と同一energyの候補評価はlockで競合を拒否する。別energyは並列実行できる。
ソース控え・checksum、入力path/size/mtime、コマンドを保存する。
上流の現在のgeometry/generatorの控えは、入力MCの作成時ソースを証明しない。
truth位置は初期一次電子の直線外挿で、物質通過後の実際のECAL入口位置ではない。

## Threshold・低信号hitの追加診断

通常の新規評価では、従来のbeam評価に加えて次の診断を自動実行する。
`--skip-threshold-diagnostics` は従来の評価だけを実行する。
従来のbeam判定・rankingと、追加診断の数値・判定は分けて保存する。

1. **SSA前の物理ストリップ**を、層別hit数・信号スペクトル・位置分布・channel占有率で比較する。
   同一イベントの物理CellID（memory-cell番号を除く）を統合し、符号付きエネルギーを合算する。
   `uncut` は追加の信号カットなしで、ゼロ・負信号のストリップも数える。
   SSA後の固有分割位置数を横に表示し、分割位置数／入力ストリップ数も保存する。
2. **共通MIPカット**をSSA前の合算信号に適用する。既定は0.5、1、2 MIP以上。
   MIPは保存MeVを既定0.305 MeV/MIPで割った値。`--mip-energy-mev` で変更できる。
   `--mip-cuts 0,0.5,1,2` のように追加のカットを指定できる（最大6個）。0 MIP以上とuncutは異なる。
   各sampleのイベント集合・元の入口軸選別をカット間で固定し、ゼロhitのイベントも分母に含める。
   hit数・平均位置・fine軸RMSの差と95%イベントbootstrap区間を保存する。
   これは保存済みhitに対する診断であり、取得時に失ったhitは復元せず、SSAも再計算しない。
3. **threshold違いの保存済みMCを比較する。** 同じ生成ROOT/truth、同じsampled entry、同じdata参照、
   同じMIPカット・bootstrap・上流物質注記を要求する。MCは全候補に共通するaxis-validイベントの共通集合で
   再集計し、候補固有の軸有効率も別途記録する。generation、digitization、calibration、SSAの自動再実行は行わない。

### 既存default・pointへの追加

元のbeam評価を上書きせず、保存済みの入力一覧・イベント集合を使って追加する。
入力ROOTのsize/mtime、Run/Event/Time、入口軸との一致を確認する。元入力が変わっていれば停止する。

```bash
cd ~/ScECAL_BeamTest/analysis/beam_quantitative
bash run.sh --energy 100 --label beam_scan --condition point \
  --threshold-diagnostics-only --dry-run
bash submit.sh --energies 100 --label beam_scan --condition point \
  --threshold-diagnostics-only
# defaultについても同様
bash submit.sh --energies 100 --label beam_scan --condition default \
  --threshold-diagnostics-only
# 保存済み追加診断だけを再描画
bash run.sh --energy 100 --label beam_scan --condition point \
  --threshold-diagnostics-only --plots-only
```

通常の`--plots-only`は従来のbeam図、`--threshold-diagnostics-only --plots-only`は追加診断図を再描画する。
追加診断の再計算による上書きは行わない。別のMIP単位・カット等を使う場合は新しいlabel/conditionで評価する。

### Threshold読み込みの監査

現在のMC配置ではdigitizationログを`Analysis_edit/job/scan/<particle/condition>/<energy>GeV/out/`から推定する。
各入力名に対応し、calibration更新時刻以前で最も新しいstdoutを調べる。
別配置は`--mc-digitization-log-dir /path/to/out`で指定する。
`GET Threshold FILE FAILED`、fallbackのADC値、parameter ROOTのパス、ログchecksumをJSONに保存し、
`threshold_input_audit.png`にも状態を表示する。対応ログの控えもinputsに保持する。
これはログとの対応に基づく証拠で、producerとの厳密な紐付けを保証するものではない。
読み込み失敗が確認された候補はthreshold parameter scanの比較から除外する。
fallback ADC=0でもpedestal・calibration側の選別は残る。

### Threshold候補の評価・比較

同じ生成ROOTから作成したthreshold 0.8/1.0/1.2倍のMCを、別conditionで通常評価する。
`--mc-threshold-scale`は**既に実行したdigitizationの由来注記**で、thresholdを変更する命令ではない。
使用したparameter ROOTが残っていれば`--mc-threshold-file`も指定する。
各候補を評価するときの形は次のとおり。3つのMC入力ディレクトリは実在するものを指定する。

```bash
bash submit.sh --energies 100 --label threshold_check --condition th0p8 \
  --mc-calib-dir /path/to/calib/th0p8/100GeV \
  --mc-ssa-dir /path/to/ssaCalib/th0p8/100GeV \
  --mc-truth-dir /path/to/original/generate/100GeV \
  --truth required --mc-threshold-scale 0.8 \
  --mc-threshold-file /path/to/threshold_0p8.root
# th1p0とth1p2も、同じ元generationを指定して評価する
bash run.sh --energy 100 --label threshold_check \
  --compare-thresholds th0p8,th1p0,th1p2 \
  --threshold-comparison-name threshold_scan --dry-run
bash run.sh --energy 100 --label threshold_check \
  --compare-thresholds th0p8,th1p0,th1p2 \
  --threshold-comparison-name threshold_scan
```

全候補でparameter ROOTを指定すると、CellIDの一致、指定倍率を除いたthreshold平均の一致、sigmaの不変を
確認する。parameter ROOTが揃わない場合は倍率を注記として扱い、数値的な検証をしていない旨を保存する。
parameter ROOTが同じだから実際のproducerがそれを使用したとは断定しない。
現在のdigitizerでは`expErfThre`の読み込みに成功した場合、末尾の整数引数ではなくROOT内の値を使用する。
したがって、末尾の整数だけ変える方法で平均thresholdの倍率scanを行わない。
今回の比較機能は平均thresholdだけを変え、sigmaを固定するscan用。

### 追加される図と数値

凡例は既存と同じcanvas上部の専用領域に置く。層別の並びは従来と同じ0–29。
SSA前の物理ストリップの図名は `strip_*` とする。`strip` はSSAで分割する前の物理ストリップを指す。
1次元位置分布は実座標の中心間隔に合わせた5.3 mm幅・42 binで表示する。
従来の1 mm binではストリップ中心の間に空binが並んでいたが、その表示上の隙間を解消する。
2次元占有率は1ストリップを1 binとし、偶数層はX 5 bin × Y 42 bin、奇数層はX 42 bin × Y 5 bin。
ストリップ本体は5 × 45 mm、bin幅はgapを含む配置pitchの5.3 × 45.4 mmとする。
各binの値はその物理channelの占有率で、ストリップ内の位置分布を測定したものではない。
`--threshold-diagnostics-only --plots-only`で、保存済みストリップイベントからこれらの表示binを更新して再描画する。
元のイベント選別、bootstrap、数値評価は変更しない。更新した図に対応する旧`pre_ssa_*` PNGは削除する。
表示binの仕様と変更ソースはthreshold JSONとinputsの`threshold_plot_sources`に記録する。

| PNG | 確認する内容 |
|---|---|
| `threshold_input_audit` | digitizationでthreshold ROOTが読めたか、fallbackの記録があるか |
| `strip_hit_counts` | 左：物理ストリップ数、右：SSA後の固有位置数、層別・95%区間 |
| `strip_signal_spectrum` | 全ストリップの符号付きMIP分布、左：全体、右：低信号域 |
| `strip_signal_by_layer` | 低信号の分布を層ごとに比較 |
| `strip_occupancy_{data,mc}_{even,odd}` | SSA前ストリップの位置別占有率、層ごとに共通のdata/MC色範囲 |
| `strip_ssa_multiplicity` | 平均SSA位置数／平均物理ストリップ数 |
| `strip_mip_cut_count` | uncutと各MIPカットでの層別ストリップ数 |
| `strip_mip_cut_difference` | カット別のMC/data−1、95%区間、設定したhit数許容差 |
| `strip_mip_cut_rms` | カット別のSSA前fine軸位置RMS |
| `strip_profiles_<cut>` | カット別のSSA前fine軸位置分布、層ごとに規格化 |
| `threshold_scan_hit_counts` | 各threshold候補のSSA前後hit数を同時表示 |
| `threshold_scan_signal_spectrum` | 各threshold候補の信号分布と低信号域 |
| `threshold_scan_mip_cut_difference` | 共通MIPカットで候補間のhit数差がどう変わるか |
| `threshold_scan_response` | 横：threshold平均の倍率、縦：全層の物理ストリップ数のMC/data−1 |

channel占有率は各イベントで0/1のBernoulli値とし、接続6300 channelすべてについてWilson 95%区間を保存する。
channelごとの一括有意差判定はしない。hit数等の許容差は従来の設定を再利用する診断用で、beam rankingに追加しない。
空分布・ゼロ幅の相対値・統計不足は未判定とする。bootstrap区間は系統誤差を含まない。

追加ROOTもcondition直下はGeVディレクトリのみ、各GeV内はROOTのみ。

```text
analysis/result/beam_quantitative/<label>/
  root/data/<phase>/reference/<energy>GeV/pre_ssa_events.root
  root/mc/<phase>/<condition>/<energy>GeV/{pre_ssa_events,pre_ssa_histograms}.root
  root/comparison/<phase>/<condition>/<energy>GeV/
    pre_ssa_data_histograms.root, threshold_metrics.root, threshold_comparison.root
  plots/<phase>/<condition>/<energy>GeV/<追加図>.png
  summary/<phase>/<energy>GeV/
    <condition>.threshold.json, <condition>.threshold.metrics.tsv, <condition>.channel_occupancy.tsv
  inputs/<phase>/<energy>GeV/
    reference/threshold_configuration.json
    <condition>/threshold_sources/
  root/threshold_comparison/<phase>/<comparison>/<energy>GeV/
    common_selection.root, comparison.root, metrics.root, <各候補の集計ROOT>
  plots/threshold_comparison/<phase>/<comparison>/<energy>GeV/<比較図>.png
  summary/threshold_comparison/<phase>/<energy>GeV/<comparison>.*
  inputs/threshold_comparison/<phase>/<energy>GeV/<comparison>/
    configuration.json, conditions.tsv, run.log, sources/
analysis/log/beam_quantitative/<label>/<phase>/<condition>/<energy>GeV/threshold_*.log
analysis/log/beam_quantitative/<label>/<phase>/threshold_comparison/<comparison>/<energy>GeV/
```
