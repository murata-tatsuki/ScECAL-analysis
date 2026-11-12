# SSA shower shape: beam-position and reconstruction checks

`beamsize_ssa` のdata/MC差を、入射位置分布・SSA分割・ヒット数/エネルギー重みの違いに分けて調べる解析。
保存済み `Calib_Hit` と `SSA_Hit` を同一イベントで比較する。再校正・SSA再生成・点入射MC生成は行わない。
既定は **100 GeV、入口layer 0–3の重心による平行移動**。全30層を出力する。

## 実行

ROOTとC++17、Python 3（標準ライブラリのみ）が必要。既存 `channel_response/root.mk` を共用する。

```bash
cd ~/ScECAL_BeamTest/analysis/ssa_shower_shape
# 入力と出力予定の確認。ディレクトリ作成・ビルド・ジョブ投入なし
bash run.sh --dry-run
# 既定: 各sampleで最大4ファイル、選別前6000イベント
bash run.sh --label baseline
# 全ファイル・全イベント。選択の異なる実行は別labelを使用
bash run.sh --label all_events --max-files 0 --max-events 0
# 再描画。元のcalib/SSA/truth ROOTは読み直さない
bash run.sh --label baseline --plots-only
# 複数energyはenergy単位で同じコードを実行
bash submit.sh --energies 20,100 --label baseline --dry-run
bash submit.sh --energies 20,100 --label baseline
```

`run.sh --help` で全引数を表示。`submit.sh --energies all` は既存解析と同じ19 energy。
利用可能性は実行時に検証する。欠落入力を別MC条件や別SSAへ自動置換しない。
`--data-calib-dir/--data-ssa-dir/--mc-calib-dir/--mc-ssa-dir/--mc-truth-dir` はenergy末端ディレクトリを指定。
`--beam-map` は入射位置分布のROOT、`--result-dir` はlabelより上の出力基点を指定する。
Slurmでは既存の `channel_response/execute.sl` を使い、標準出力・エラーの保存先を上書き指定する。

## 入力と対応検証

既定入力（`<beam>` は5 GeV以下がps、それより上がsps）：

- data: `ECAL_data/analysed/2023/<beam>/{calib,ssaCalib}/e-/<energy>GeV/*.root`
- MC: `simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/{calib,ssaCalib}/e-/e_ssa/threshold/<energy>GeV/*.root`
- truth: 同じ `Result_MC/generate/e-/e_ssa/<energy>GeV/*.root` の `MC_Truth`
- beam template: `analysis/result/trackFit/beamsize_ssa/<beam>/e-/<energy>GeV.root` の `h_hitmap_smooth`

1. 同名ファイルを対応づける。calib全ファイルに対応するSSAが必要。
2. 既定のファイル制限は、ソートした一覧を等分し各区間の中央ファイルを選ぶ。全runの代表性は保証しない。
3. 選択ファイル内の全entryを連結して等分し、各区間からseed付き乱数で1 entryを抽出。先頭だけを使わない。0指定は制限なし。
4. calib/SSAの全entry数が等しいこと、採用候補entryの `(Run_Num, Event_Time, Event_Num, DetectorID)` が一致することを要求。並べ替えを推測して結合しない。
5. 両方に存在する**元のCellID**ごとに符号付きエネルギー和を比較。許容差は `1e-7 + 1e-7 * sum(abs(pre energy))` MeV。不一致・SSA側だけの親CellIDは停止。
6. SSAに存在しない親ストリップは除外理由を推定せず診断する。前後比較の対照条件に `pre_present` と `pre_split_parents` を用意する。
7. finite値・vector長・ssaTag=0/1・layer Zを確認。calibのX/Yは既存EBUdecodeのストリップ中心と照合する。

MC truthは独立イベントIDを持たない。現在のMCDigiがentryをEvent_Time/TriggerIDに保存する規則に従い、同名ファイル・同entry数・`Run_Num=0, Event_Time=Event_Num=entry` を全採用候補で要求して結合する。
一次粒子PDG=11、PrimaryEnergyが要求energy（GeV→MeV）と一致、前向きの有限directionも検証する。
**この照合は別productionの同名truthファイルを識別するものではない**。入力の由来をそろえること。
`--truth auto` はtruthディレクトリ全体が存在しない場合のみ無効化し、理由を表示・保存する。部分欠損・不一致は停止。
`--truth required` は存在を必須にし、`--truth off` は使用しない。
beam templateは現在の参照候補であり、そのROOTがMC生成当時と同じだったことは自動証明できない。

## 共通母集団・座標・重み

追加のCoG cut、energy窓、tail cut、不良チャンネルmaskは適用しない。元のbeamsize_ssaの差を調べるbaselineとして扱う。
物理的に接続されたlayer 0–29、chip 0–5、channel 0–35（chip5は0–29）を共通に採用。その他は数えて除外する。
既存SSAの非接続channelなどが含まれる場合、元の図とはその部分が異なる。

- **軸**: SSA前の入口4層で、正のHit_Energyを重みにXを奇数層1/3のX、Yを偶数層0/2のYから計算。長手方向の粗いストリップ中心を軸推定に使わない。
- `--axis-layers 4` で入口層数を変更可能。層間の等重み平均ではなく、各方向の採用ヒットのenergy重心。X/Y両方の正の信号和が必要。
- 共通比較は軸が定義できるイベントだけを使用。lab/alignedおよびSSA前後で**同じイベント**を使う。軸未定義イベントはevents treeとcutflowに残す。
- `lab`: 検出器の元の座標。
- `aligned`: 上記の固定X/Yを全層から引く。入射角や層ごとのシャワー揺らぎは補正しない。
- `truth`: MC一次電子の直線を各hit Zへ外挿して引く。MCだけの補助比較。dataとの比較には共通のaligned軸を用いる。
- 1D座標 `u`: 奇数層X、偶数層Y。SSA前の細かく測れる方向のみを用い、前後で同じ方向を比較する。
- `count`: 同一イベント内の完全に同じ `(layer,x,y)` を1点にまとめた位置数。0/負energyも数える。preでは同じ位置のストリップ信号をまとめる。
- `positive_energy`: 同一位置の**符号付きenergyを合算した後**、`max(0,sum)` を重みにする。MeV。入射ビームenergyではなく校正済み信号。
- 負energyのヒット数と符号付き総和、正成分和を別途保存する。SSA分割でenergyの符号が変わる場合、正成分だけの和は保存されないことに注意。

入口層自身には軸推定による自己相関がある。layer 0–3のaligned幅の縮小をそのままビーム幅の改善と解釈しない。
シャワーの層ごとの重心を引く解析ではない。

## 比較段階

| stage | 使用する信号 |
|---|---|
| pre | 接続されたSSA前の全ストリップ |
| pre_present | SSA後に元CellIDが残るSSA前ストリップ |
| pre_split_parents | ssaTag=1の出力がある元CellIDのSSA前ストリップ |
| post_all | SSA後のtag 0/1を両方含める |
| post_split | SSA後のtag 1のみ。元のbeamsize_ssaの選別に対応 |

`pre → pre_present` で親ストリップ欠落、`pre_present → pre_split_parents` で分割成功の選別、`pre_split_parents → post_split` で分割後の位置/重み変化を調べる。
SSA前後のenergy照合はこの選別より前に行う。異なる校正入力を混ぜた比較を許さない。

調査時の現行StripSplit.cxxには `chipIDs = (CellID/100000)%10` の後に `chipIDs==5 && channelIDs>29` を除外する記述があった。
これは通常のchip桁ではなくlayerの下一桁を参照する。入力例でもlayer 5/15/25の該当親ストリップ欠落を確認した。
本解析では元の処理や既存ROOTを修正せず、parent_audit.tsvで実際の欠落を記録する。履歴不明の全入力にこの理由を一律適用しない。

## 出力配置

既定の基点は `analysis/result/ssa_shower_shape/<label>/`。
CoG cutを使わないため、形式だけの `200mm` 階層は追加しない。

```text
inputs/<energy>GeV/
  configuration.json, data.tsv, mc.tsv, run.log, complete.json
  sources/, geometry/, current_upstream_sources/, plot_sources/
{data,mc}/<stage>/<lab|aligned|truth>/<count|positive_energy>/<energy>GeV/
  histograms.root, summary.tsv
{data,mc}/diagnostics/<energy>GeV/
  diagnostics.root, validation.tsv, parent_audit.tsv, cutflow.tsv
comparison/<stage>/<frame>/<weight>/<energy>GeV/
  comparison.root, fine_profiles.png, fine_profiles_log.png
  # post_splitだけ: data/mc_map_even/odd[_log].png
comparison/overview/<energy>GeV/
  comparison.root, summary.tsv, beam_closure.tsv
  rms_overview_count.png, rms_overview_positive_energy.png
  beam_closure.png, data_diagnostics.png, mc_diagnostics.png
```

ビルド・実行・Slurmログ: `analysis/log/ssa_shower_shape/<label>/`。ファイルごとの実行状況は `<energy>GeV/` 以下。
ROOTと図は既存と同様に保持し、canvasはcomparison.rootへ保存する。
同じlabel/energyの再集計は拒否、同時書込みはlockで拒否。失敗途中の結果にcomplete.jsonは付かない。
`--plots-only` は成功した集計を再利用し、図と比較ROOT・比較TSVだけを更新する。
inputsにはコマンド・入力絶対path/サイズ/mtime、ROOT版、現行ソースの控えとSHA256を記録。巨大な入力ROOT自体の全ハッシュは計算しない。
上流ソースの控えは**現在のコード**で、保存済みROOTの作成時ソースを証明するものではない。

## 図と統計

図の見出しにはenergy、SSA前後の処理内容、座標の基準、重みを表示する。
内部のstage/frame名は保存先のパスに残し、図中では次の説明に置き換える。

| パス内の名前 | 図中の表示 | 意味 |
|---|---|---|
| pre | Before SSA: all strips | SSA前の全ストリップ |
| pre_present | Before SSA: strips retained after SSA | SSA後に元CellIDが残るストリップ |
| pre_split_parents | Before SSA: strips that were split | SSAで分割された元ストリップ |
| post_all | After SSA: all hits (split + unsplit) | 分割・未分割を含むSSA後の全ヒット |
| post_split | After SSA: successfully split hits | SSA後の分割成功ヒットのみ |
| lab | Detector coordinates | 検出器の元の座標 |
| aligned | Centered on entrance centroid | イベントごとに入口4層の重心を引いた座標（既定） |
| truth | Relative to simulated electron track | MCの一次電子軌道からの距離 |
| count | Hit-position counts | 重複位置をまとめたヒット位置数による重み |
| positive_energy | Positive signal-energy weights | 同一位置の符号付き和の正成分による重み |

比較図の凡例は **Experimental data / Simulation** と表示し、各パネルの外に設けた図上部の共通欄に置く。
入射位置の比較図も5系列の凡例を図上部へ置く。2D図は図上部にsample・処理条件を表示し、色軸用の余白を確保する。
層別1D図は各パネルにlayer番号と測定方向（奇数層X／偶数層Y）を表示する。
幅の一覧図は左列が元の座標、右列が入口重心による位置補正後で、各行の見出しがSSA前後の段階を示す。
この表示変更はヒストグラムのbin、規格化、集計値、採用イベントを変更しない。

- 層別fine_profilesは5×6パネル。data黒、MC赤。全sampleで同じbin・軸、flowを含む全重みで確率規格化。表示範囲は±120 mm、保存範囲は±250 mm。linear/logを保存。
- post_splitの2D mapは既存beamsize_ssaと同じ偶数/奇数層4×4パネル。**イベントあたり**の位置数または正信号量。対応するdata/MCの同じ層で共通色軸。表示±115 mm、保存±226.6 mm。真のdeposit位置の測定とは解釈しない。
- histograms.rootは規格化前の1D/2D、絶対距離、post_splitの2D半径、採用イベント数を保持。
- summary.tsvのmean/RMSは全ヒット位置から直接加算した値（bin依存なし、flowも含む）。RMSは重み付き標準偏差。
- `abs_q50/90` は基準軸からの1D絶対距離の包含幅。labは原点基準。`radial_q50/90` はpost_splitのみの2D包含半径。1 mm bin内で線形補間し、範囲外の分位点はNaN。平均位置まわりのRMSと意味が異なる。
- `abs_tail30_fraction` は|u|>=30 mmの全重みに対する割合。flowも含む。empty layerはNaN。
- countはヒットが多いイベントの寄与が大きい。positive_energyは信号が大きいイベントの寄与が大きい。イベントごとに面積1へ規格化する解析ではない。
- ROOTのSumw2はヒット重みの二乗和。ヒット間のイベント内相関を考慮した誤差ではないため、図は線で表示し、有意差やfitの判断には使わない。イベントbootstrapや系統誤差評価は未実装。

## 診断ROOT

`events` は全sampled eventのfileIndex/entry/Run/Event/軸valid、層別点数、負energy数、符号付き/正成分energy、欠落親数・energyを保存する。
fileIndexはinputsのsample.tsvの行番号（0始まり）。npost/nsplit/nunsplitはそれぞれ固有位置数なので、tag 0/1が同じ位置に重なる場合はnpost=nsplit+nunsplitにはならない。
エネルギー配列は位置統合前のhit単位。summaryのpositive_energyとは負成分との位置統合により差が出る場合がある。

- entrance_axis: 共通母集団の再構成入口重心。
- primary_xy / primary_at_layer0: 全sampled MCの一次電子を1event 1回数える。primary_at_layer0_selectedは共通母集団のみ。
- post_layer0_all/smooth: 軸未定義も含む全sampled eventのSSA tag1位置。smoothはseed固定の±2.7 mm一様乱数付き。元の図と同じ数え方の診断で、中心軸の比較とは別母集団。
- beam_input_template: 指定した入力ビームmapの控え（MC側のみ）。beam_closureはこれらを個別に確率規格化して比較する。ヒット重みとevent重みの違い自体が検証対象。
- axis_residual: 共通入口重心とtruth軸のlayer 0位置の差（MC）。平行移動軸と傾いたtruth軸の差を含む。
- strips_by_absu/split_by_absu/unsplit_by_absu: 共通母集団での層・軸からの距離別点数。後二つの比から位置ごとのtag内訳を確認可能。
- pre/post_energy_vs_absu: 位置ごとの符号付き信号と軸からの距離。範囲外はROOT flowに保持。
- parent_audit.tsv: 物理CellIDごとの親出現数・欠落数・tag1出力のある回数・符号付きenergy和。memory桁をまとめるが、対応照合自体は元CellIDで行う。

## 解釈の順番

1. parent_audit/validationでSSA入力の一致と欠落を確認する。
2. beam_closureで入力template→一次電子→MC layer 0の変化を見る。入力templateとdata layer0は統計範囲・平滑化seedが異なり完全一致を要求しない。
3. 同じstage/weightのlab→aligned、MCのaligned→truthを比較する。
4. 軸と母集団を固定して、stageとcount→positive_energyの変化を見る。
5. ビーム幅の二重計上を確定するには点入射MCなどの別検証が必要。本解析で物理モデル差や原因を自動断定しない。
