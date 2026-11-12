# Shower leakage signal proxies (calibrated hits only)

C++/ROOT共通処理。初期比較は5 GeV (PS)、20/100 GeV (SPS)。現在は0.5–250 GeVの19 energyに対応。
保存済みCalib_Hitと既存HG–LG tail判定を使用し、MC truth・再校正は使わない。

```bash
cd ~/ScECAL_BeamTest/analysis/shower_leak
bash run.sh --energy 5 --label initial --max-events 6000 --max-files 4
bash run.sh --energy 20 --label initial --max-events 6000 --max-files 4
bash run.sh --energy 100 --label initial --max-events 6000 --max-files 4
# 保存済みevent集計から再描画（入力calibration ROOTへアクセスしない）
bash run.sh --energy 100 --label initial --plots-only
# 既存環境のSlurm wrapperへ投入。既定は3 energyのみ
bash submit.sh --energies 5,20,100 --label larger --max-events 30000 --max-files 10 --dry-run
```

`--max-events`はsampleごとの選別前上限、`--max-files`はsampleごとの上限。0で制限なし。
ファイルは全一覧から等間隔で選び、その連結entryを等分した各区間からseed付き乱数で1 eventを抽出する。
ファイル先頭には偏らず、選んだファイル内ではevent数に比例する標本となる。
ただし少数ファイルの標本は全runの代表性を保証しない。`--seed`で再現、entryも保存する。
`--edge-mm`は既定10 mm。変更時は別labelで再集計する。
`--data-dir/--mc-dir/--tail-results-base/--result-dir`で入力と保存先を変更可能。
同じ出力を使う同時実行とevents.rootの上書きを禁止する。

## 共通定義

- Calib_HitのCellID (vector<int>)、Hit_Energy/Hit_X/Hit_Y/Hit_Z (vector<double>)、Run_Num/Event_Time/Event_Num (int)を型確認して読む。各fileのschemaを保存する。
- Hit_EnergyはMeV。1000で割ったGeVを保存する。この量は検出器の**calibrated signal sum**であり、beam energyや漏出energyではない。
- `StudyCommon.hh::physical/connected`を共用。layer 0–29、chip 0–5、channel 0–35（chip5は0–29）。memory cell桁をmask照合時に落とす。同一stripの複数memory hitはsignal・hit数にそのまま足す。
- resolutionのbad_channels集合を既存hg_lg_tailと同じ方法で抽出し、保存済みtail解析のmaskとも集合一致を要求する。bad channelはdata/MC共通。非接続channel、負energy hitを除外。0 energy hitはhit数に含む。NaN/Inf・vector長不一致・geometry不一致は停止。
- 両sampleで既存tail_richを除外する。`HLTailSelection.hh`がcanonical入力path、全entry数、全entryのRun/Event_Time/Event_Numを照合する。判定のselectedフラグは使用しない。既存の3 channel以上・z<-5という判定を再定義しない。
- tail参照の全calibration入力mtimeがtail_study.root以下であること、保存された校正パラメータchecksumが一致することを確認する。欠落/不一致/更新された入力では停止し、keep条件へfallbackしない。元判定が入力calibrationのhashを持たないため、mtimeを保持した書換えまで完全には検出できない。
- 正のtotal signalを持つeventがfraction診断の母集団。total=0のeventは保存し、cutflowで別計数する。狭いenergy窓・peak cutは適用しない。
- `200mm`表示はnocutで、位置cutを本当に適用しない。layer9/10のCoG未定義eventも保持する。`20mm`は両layerでenergy>0かつabs(x),abs(y)<20 mm。未定義は不合格。
- layer番号はEBUdecodeでz増加方向。upstream layer0 (z=1 mm)からdownstream layer29 (z=290.8 mm)。最後3層=27–29、最後5層=25–29。実入力Hit_Zと全採用前の接続hitで照合する。
- energy重み付き深さはmmとlayer番号の両方。最大応答layerは同値なら小さい番号。layer CoG・広がりはstrip中心の重み付き平均・標準偏差。空layerはNaN。

## 物理外縁

`analysis/calibration_parameters/EBUdecode.cxx`を直接リンク。strip実寸5×45 mm、pitch 5.3×45.4 mm。
偶数layerはx方向が長軸、奇数layerはy方向が長軸。各層の外周は全接続stripの外側面：
偶数layerで|x|=113.3、|y|=111.15 mm、奇数layerで逆。
外側面から幅wの帯と**正の面積で交差するstrip全体**を採用する：
`abs(center_x)+half_x > boundary_x-w`（yも同様）。
edge_x/edge_yを別保存し、edgeは和集合（corner二重計数なし）。mask境界は外周にしない。
stripのenergyを面積比で分配しない。長軸方向では45 mmのどこに信号が生じたか不明で、10 mm帯でも端の45 mm strip全体を数える。
これは物理的に10 mm帯へdepositしたenergyそのものでも、検出器外に漏れたenergyでもない。
mapは赤=x外縁、青=y外縁、紫=両方、灰=masked、白=内側。strip_map.tsvにはmask適用前の外縁フラグも残す。

## 保存先

`analysis/result/shower_leak/<label>/` 以下：

- `inputs/<energy>GeV/`：実行コマンド・条件・ソース・checksum・tail入力一覧の控え。
- `{data,mc}/exclude/200mm/<energy>GeV/`：events.root（除外前全sampled event）、cutflow.tsv、validation.tsv、strip_map.tsv、inputs/files.tsvとbranches.tsv。
- `{data,mc}/exclude/{200,20}mm/<energy>GeV/histograms.root`：規格化前histogram/profile、選別数、event_source。
- `comparison/exclude/{200,20}mm/<energy>GeV/`：PNG、comparison.root、summary.tsv、longitudinal.tsv、histogram_flows.tsv、indicator_intervals.tsv、impact.tsv、relative_variance.tsv、tail_association.tsv、run_dependence.tsv。
- strip mapはcomparisonの200mmだけにPNG2枚・ROOTを保存する。channelごとの図は作らない。

実行・ビルド・Slurmログは`analysis/log/shower_leak/<label>/`。途中失敗の詳細もここに記録する。
Aggregateに1 sampleの入力一覧（1行1path）とsample種別を渡し、data/MC同一処理を使用する。
Plotは両sampleのevents.rootだけを読んで全plot・tableを生成する。CoG20では集計treeを複製しない。

## 統計と解釈

1D分布はdata/MC共通軸・binning、sample内で確率規格化。2Dは各sampleの確率/bin、同じ色範囲。
ROOTには規格化前件数とTProfileを保存。横軸の最大値は共通sample全体を含める。
fractionの0は専用intervalとzero件数で明示し、1は最後のbinへ含める。flow件数も記録する。
区間別指標は=0と(0,.001],(.001,.003],(.003,.01],(.01,.03],(.03,.05],(.05,.1],(.1,.2],(.2,1]。

- MeanとRMSはbinに依存しない全event値から計算。RMSは母標準偏差。
- 中央68%は`(q84-q16)/(2*median)`（半幅）を使用。線形補間の分位点。
- Gaussianは100 binの分布の[q16,q84]にfitする固定手順。これはevent選別ではない。
  N>=100、fit成功・covariance品質>=2、正のmean/sigmaを要求。失敗はNaN、statusとfit範囲を保存。
  fit幅はcoreの記述量であり、RMSとは異なる。少数binしかないcoreでは不安定になる。
- N<100の区間はLOW_STAT、図では白抜き。Gaussianは省略する。
- 閾値は事前固定の0,.001,.003,.01,.03,.05,.1,.2。全5指標で同じ閾値を両sampleへ適用し、`fraction>threshold`をtag、`<=`を保持。
- above_fractionにはbinomial SEと68%相当Wilson区間(z=1)を併記。0件でも区間上限が残る。
- impact.tsvはbaselineとcut後のN・保持率・mean/各幅の差、低energy裾割合を記録する。
- relative_variance.tsvの差は`R_data^2-R_MC^2`。RはGaussian、中央68%、RMSの各相対幅。
  厳密な相対分散`Var(E)/mean(E)^2`に対応するのはRMS行。他2行は幅の二乗の比較である。
- 裾診断はsampleごとのcut前medianを固定し、低側E<0.8median、高側E>1.2medianと定義する。
  event除外には使用しない。tail_association.tsvは各領域で指標閾値を超える割合を示す。
- 比率の分母による相関を見るため、同じ図に後方/外縁の絶対signalとの相関も保存する。
- 統計誤差は独立eventを仮定。run間変動・校正systematic・sampling/file選択の不確かさは含まない。
  幅の変化量や相対分散差の不確かさはこの初期解析では推定しない。

**後方/外縁signal割合は漏れの指標であり、実際の漏出energy割合ではない。**
cut後の幅の縮小だけではleakage寄与は確定できない。shower形状・位置・既存選別の効果が混ざる。
特定の改善率へ閾値を最適化しない。絶対漏出energyと原因の分離はこのhit-only解析では判断できない。

## 初期実行での確認（initial）

各energyでdata 6000 event（5/20/100 GeVで4/2/3 file）、MC 4000 event（各4 file）。
MC入力は1 file 1000 eventのため、max-events=6000でも4000になる。
両sampleのmaskとtailの全entry対応が一致し、座標もEBUdecodeと一致した。
既存beam設定`CEPCScECAL/generate/mac/e-/5GeV/5GeV_0.mac`の`/gps/direction 0 0 1`と、
既存解析のfirstLayer=0、z増加の配置を照合して後方を高layer番号とした。primary truthは読まない。

±20 mm・低LG tail除外・energy窓なしの結果：

| GeV | N data / MC | Mean signal [GeV] data / MC | Gaussian σ/mean data / MC | 中央68%半幅/median data / MC | RMS/mean data / MC |
|---|---:|---:|---:|---:|---:|
|5|2444 / 1940|0.15294 / 0.21344|0.09588 / 0.06800|0.11225 / 0.07712|0.28467 / 0.07840|
|20|2635 / 1808|0.84014 / 0.85408|0.06639 / 0.04126|0.07762 / 0.04363|0.12794 / 0.04427|
|100|5081 / 3036|4.26423 / 4.03151|0.04603 / 0.03259|0.04698 / 0.02851|0.10891 / 0.02971|

最後5層signal割合>0.03のevent率（binomial SE、改善率で選んだ閾値ではない）：
5 GeVでdata 10.27±0.61%、MC 5.26±0.51%、20 GeVで6.03±0.46%、11.17±0.74%、
100 GeVで31.14±0.65%、56.36±0.90%。「dataの後方指標が常に大きい」とはならない。
100 GeVの同cutでは保持率data 68.86%、MC 43.64%。dataのRMS/meanは0.10891→0.07447、
MCは0.02971→0.02759、平均signalは4.26423→4.31205と4.03151→4.05435 GeV。
この差には保持率と形状・位置選別の違いが含まれる。

100 GeVで外縁割合>0.03はdata 19/5081（0.374±0.086%）、MC 0/3036。
MCの0件は真の率0を意味せず、Wilson上限(z=1)は0.0329%。dataの同cut後RMS/meanは0.09385。
外縁はmaskで観測されないstripが多く、かつ45 mm長軸の局在は分からないため、側方漏れ全体の頻度には換算できない。

5 GeV dataでは±20 mmでも低signalの集団が残り、RMSとGaussian coreの幅が大きく異なる。
nocutではCoG未定義も含む広い集団になり、5 GeV dataのGaussian fitは失敗（NaN）。
この集団をhit-only情報だけから漏れeventと同定することはできない。

合成12 eventで後方3/5層、corner和集合、memory-cell mask、負/zero energy、未定義CoG、tail除外、
外縁幅変更を検証。実sampleの保存treeとhistogramのN・energy/hit和・fractionも照合する。

追加の異常系検証では、entryのEvent_Num不一致と、tail結果より新しいcalibration入力を
いずれもエラーで停止することを確認した。再描画前後の全数値表のchecksum一致も確認した。
Slurm投入はdry-runを確認し、初期比較はローカル実行した。全energy実行は行っていない。
CoG依存図はfine segmentation側のlayer9 CoG xとlayer10 CoG yを表示する。
全layerのx/y CoG・幅はevent treeとhistograms.rootに保存している。

## interval_metricsの表示と表示だけの更新

横軸は0と(0,0.001]だけを左側のカテゴリ領域に残し、0.001–1を対数軸で表示する。
点線がカテゴリ領域と対数領域の境界。正の各区間のmarkerは幾何平均位置に置き、
横エラーバーは区間の下限・上限を表す（統計誤差ではない）。data/MCで同じx位置を使う。
0は幅のない点、(0,0.001]の横幅はカテゴリ表示上の幅で、対数距離ではない。
実装上は右側をu=log10(f/0.001)で描画し、TGaxisで実際のfraction値の対数目盛を付ける。
ROOT graphのx座標もこの描画座標になる。物理的な区間境界は従来どおりindicator_intervals.tsvに保存する。
各panelのplot領域の上に独立した凡例領域を設け、青=data、赤=MC、白抜き=N<100を示す。
選別、区間境界、fit値、縦の誤差は変更しない。

```bash
# interval_metricsのPNGとcomparison.root内の該当canvasだけを更新する。
# 保存済みcanvasから元の精度で値と縦誤差を読むため、再集計・再fitしない。
for e in 5 20 100; do
  bash run.sh --energy "$e" --label initial --intervals-only
done
```

`IntervalPlots.hh`は通常描画と表示だけの更新で共用する。
旧カテゴリ軸canvasと更新後のcanvasの両方から繰り返し再描画できる。

strip mapの右下の凡例は、実際のstripと同じ色のboxにx edge／y edge／x and y edges／
Masked strips／Interior stripsのラベルを添える。白いboxにも枠線を付けて区別する。
`bash run.sh --energy 100 --label initial --maps-only`で、保存済み外縁幅・maskから
strip mapのPNG2枚とstrip_maps.rootだけを更新できる。eventの再走査・fitは行わない。

## Optional sparse through-going track veto

既定は `--track-cut off`（従来通り）。`on` は低信号・少hitで前後に続く細い直線track候補を除外する。
**muon/pionのPIDではない**。非相互作用pionやmuonと整合する形状を選ぶ操作的なtagであり、
電子の低応答eventの混入・候補の取り逃しもあり得る。hadronic showerを作るpion全般を除く選択ではない。
resolution改善率やdata/MC一致を目標に閾値を最適化していない。

次の全条件のANDを候補とする。同じ閾値を5/20/100 GeVのdata/MCへ適用する。

| 条件 | 初期値 / option | 意図 |
|---|---|---|
| 正のtotal信号が上限以下 | 0.03 GeV / `--track-max-energy-gev` | 主showerより十分低い信号領域に限定。beam energyやMIP単位ではない |
| 採用hit数が上限以下 | 60 / `--track-max-hits` | 少hit |
| 正の信号があるlayer数が下限以上 | 15 / `--track-min-layers` | 検出器を縦に通る候補 |
| 最初5層・最後5層のactive layer数がそれぞれ下限以上 | 2 / `--track-min-end-layers` | 前方だけ/後方だけのeventを区別 |
| 採用hit数 / active layer数が上限以下 | 2 / `--track-max-hits-per-layer` | 多hitのshowerを保護 |
| x,yの直線からの広がりがそれぞれ上限以下 | 5.3 mm / `--track-max-rms-mm` | 短辺pitch程度の細いtrack |

active layerは既存maskとhit条件を適用した後の `layer_energy>0`。
hit数は元解析と同じで、0信号hitや同じphysical stripの複数memory hitを含む。
前後5層は0–4、25–29。欠損layerを許すため連続30層は要求しないが、maskによる取り逃しは残る。
このANDなので、低energyだけ、少hitだけ、後方hitだけでは除外しない。
0.03 GeVは5 GeVの主集団（平均信号約0.155 GeV）より低い初期診断用の固定値で、
MIP校正から導いたPID境界ではない。異なる校正設定への無検証な流用はしない。

直線性は保存済みlayer momentだけで評価する。奇数layerの短辺方向xと偶数layerの短辺方向yを、
各々 `q=a+b*z` に等layer重みでfitする。各投影で4層以上必要。
`RMS_q = sqrt(mean_layers[(CoG_q-a-b*z)^2 + width_q^2])` を使う。
layer内の広がりも含むので、広いshowerが直線状のCoGを持つだけでは候補にならない。
長辺45 mm方向のstrip中心を精密なhit位置としてfitしない。
RMSは幾何学的な広がりであり、fitの統計誤差やchi-squareではない。
曲がったtrackや途中から生じたtrackは取り逃す保守的な選択である。

```bash
# 既存initial集計を再利用。元のcalibrated ROOTは再走査しない。
bash analysis/shower_leak/run.sh --energy 5 --label initial --plots-only --track-cut on
bash analysis/shower_leak/run.sh --energy 20 --label initial --plots-only --track-cut on
bash analysis/shower_leak/run.sh --energy 100 --label initial --plots-only --track-cut on
# cutなし（元のexcludeへ再描画）
bash analysis/shower_leak/run.sh --energy 100 --label initial --plots-only --track-cut off
# 新規入力集計にも同じoptionを使用可能
bash analysis/shower_leak/run.sh --energy 100 --label track_sample --max-events 6000 --max-files 4 --track-cut on
# Slurm：投入前にコマンド確認
bash analysis/shower_leak/submit.sh --energies 5,20,100 --label initial --plots-only --track-cut on --dry-run
```

集計events.rootは選別前全eventなので、`{data,mc}/exclude/200mm/<energy>GeV/events.root`を共用する。
`on`の結果は `{data,mc,comparison}/exclude_track/{200,20}mm/<energy>GeV/` に保存し、
既存`exclude`の図・数値を上書きしない。保存先の200mmは位置nocutを意味し、track vetoの有無は別条件。
処理順はlow-LG tail除外→正のtotal→track veto→任意のCoG20。
CoG未定義でも直線fitに十分な他layerがあれば判定可能。nocutではtagされなければ保持する。

追加保存物：

- sample別 `track_selection.root`：全sampled eventのfile_index/entry/run/time/event、各条件の判定、
  layer数、hit密度、x/y RMS、候補tag、同じCoGでのcut前/後の採否。
  元入力fileは元の`exclude/200mm/<energy>GeV/inputs/files.tsv`で照合する。
- sample別 `track_cutflow.tsv`：各条件までの累積候補数、実際の除外数、CoG除外数、残存数。
- comparisonの `track_cut_impact.tsv`：同じCoGのbefore/after/removedの件数・保持率・平均・各幅。
  offではtagを診断するが実際のremovedは0。幅の定義は従来と同じ。
- `track_configuration.txt`とROOTの`track_definition`：具体的な条件式。
- 再描画のソース・入力event checksum・コマンド：`inputs/<energy>GeV/track_cut_on/`または`track_cut_off/`。
- ログ：`analysis/log/shower_leak/<label>/replot_<energy>GeV_track_on.log`など。

閾値変更は保存eventsから再選別する`--plots-only`で行う。同じon出力先を更新するので、
閾値違いを同時保持したい場合は別labelを使う。保存済みeventの再利用は後述の`--source-label`を指定する。
`--intervals-only --track-cut on`は保存時と同じ閾値が必要で、違えば停止する。
`--maps-only --track-cut on`は同じgeometryのmapをon側へ描画するだけでevent選択を変更しない。
track veto後にも残る低信号eventをmuon/pionではないと断定しない。
また、このvetoによる幅の縮小をleakage寄与の測定と解釈しない。

## Readable axes and all-energy execution

全plotの軸目盛り・軸titleを拡大し、余白とcolor bar位置を共通設定する。
`PlotStyle.hh`が描画時に適用され、intervalの特殊な横軸やstrip mapにも適用する。
保存済みhistogram・profileとinterval canvasの点から表示だけを再構築する場合は次を使う。event集計・選別・fitは実行しない。

```bash
bash analysis/shower_leak/run.sh --energy 100 --label initial --style-only --track-cut both
```

`--track-cut both`は同じ集計からoff→onを順に描画する。low-LG tail除外は両条件で有効。
PSは0.5,1,2,3,4,5 GeV、SPSは10,20,30,40,50,60,70,80,100,120,150,200,250 GeVに対応する。
`submit.sh --energies all`でこの19点を指定できる。既定の対象は従来の5,20,100のまま。
energyディレクトリ・既存tail判定が欠ける/一致しない/古い場合は停止し、選択を省略しない。
0.5 GeVも浮動小数energyとして扱い、保存先は`0.5GeV`。
track vetoの閾値は全energyで共通だが、特に低beam energyでの電子効率は個別のcutflowで確認する。

```bash
cd ~/ScECAL_BeamTest
# 全19 energy、全入力file・全event、track cutのoff/on両方。まず投入コマンドを確認。
bash analysis/shower_leak/submit.sh --energies all --label all_energy \
  --max-events 0 --max-files 0 --track-cut both --dry-run
# 上記と同じ設定で実際に投入：1 energy = 1 job、計19 jobs。
bash analysis/shower_leak/submit.sh --energies all --label all_energy \
  --max-events 0 --max-files 0 --track-cut both
# 小規模に全energyを確認する場合は別labelと制限を指定する。
bash analysis/shower_leak/submit.sh --energies all --label all_energy_sample \
  --max-events 6000 --max-files 4 --track-cut both
# 完了した全energyの集計から再選別・再描画する場合
bash analysis/shower_leak/submit.sh --energies all --label all_energy \
  --plots-only --track-cut both
```

`--max-events 0 --max-files 0`は全件処理。省略すると各sample最大6000 events・4 files。
新規集計時は未使用labelを指定する。既存events.rootの上書きは行わない。
出力は`analysis/result/shower_leak/<label>/comparison/{exclude,exclude_track}/{200,20}mm/<energy>GeV/`、
ログは`analysis/log/shower_leak/<label>/`。
`exclude`はtrack cutなし、`exclude_track`はtrack cutあり。両方とも既存low-LG tail除外後。

## v2: 低信号ピークを除く診断用の選択

`--track-cut v2` は従来の `on` のtag **または**
`beam energy >= 10 GeV && 0 < total calibrated energy <= 0.1 GeV`
を満たすeventを除外する。10 GeV未満では従来のonと同じ。
0.1 GeVはhitの校正信号の和であり、beam energyに換算した再構成energyではない。
低信号条件はhit数、貫通、直線性を要求しない。data/MCに同じ式を適用する。
閾値は `--low-signal-max-energy-gev` で変更できる。
これは0付近のピークを除いた診断用のenergy選別であり、muon/pionの同定ではない。
実際のleakageや検出器応答異常による低信号も除き得るため、cutなしとの比較を保持する。
幅の改善率を目標に閾値を最適化したものではない。

既存 `off`、`on`、`both` の定義は変更しない。v2は
`<label>/{data,mc,comparison}/exclude_track_v2/{200,20}mm/<energy>GeV/`
に出力する。200mmは位置nocut。既存のlow-LG tail除外は維持する。
図のcaptionにv2を表示し、`track_configuration.txt` とROOTに条件式を保存する。
`track_selection.root` の `old_track_tag` は旧tag、`low_signal_tag` は上記低信号条件、
`track_tag` は選択modeの候補（v2ではOR）、`selected` は実際の採否。
cutflowの旧累積条件は旧tagのまま、`candidate_low_signal` と
`candidate_combined_veto` を追加する。実際の除外数・保持率・幅は
`track_cut_impact.tsv` のbefore/after/removedを確認する。

全energyの保存済みeventを再利用し、新しいlabelへv2だけを出力する例：

```bash
cd ~/ScECAL_BeamTest
bash analysis/shower_leak/submit.sh --energies all \
  --source-label all_energy_track_off --label all_energy_track_v2 \
  --plots-only --track-cut v2 --low-signal-max-energy-gev 0.1
```

1 energy = 1 Slurm job。投入確認だけなら `--dry-run` を追加する。
`--source-label` は初回のみ、未使用の出力label（energyごと）に指定する。
入力のevents.rootを読み、calibration ROOTは再走査しない。
max-events/max-filesは保存eventの再抽出には使わず、元の集計全件を処理する。
元の入力記録を `inputs/<energy>GeV/source_inputs/` に控え、sample別file一覧、
参照先 `input_base.txt`、実行コマンド、入力ROOTのchecksum、現在のソースを保存する。
元labelのevents.rootは引き続き再解析の入力として必要なので保持する。
再描画・再選別時は `--source-label` を省略すれば記録した入力先を使う。

```bash
bash analysis/shower_leak/run.sh --energy 200 --label all_energy_track_v2 \
  --plots-only --track-cut v2 --low-signal-max-energy-gev 0.1
```

calibration ROOTから新規集計する場合も `--track-cut v2` を使える。
全件処理には従来どおり `--max-events 0 --max-files 0` を指定する。
