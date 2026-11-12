# Resolution の sample 比較

現在は SingleEnergy / MultiEnergy の energy–Nhit 楕円 event cut と `afterEventCut` の生成・比較・専用計算を一時停止し、`beforeEventCut` のみを処理する。既存の after-cut 結果は保持する。

今後の標準実行ではHG–LG tail-rich eventを除外する。`--keep-tail-events` を明示した場合だけ保持する。
対応するtail判定結果がない入力では、その判定結果を先に用意し、keepへ自動で切り替えない。

SingleEnergyのsummaryでは、energy deposit・number of hitsの表示範囲をそれぞれのGaussian fitの
平均±5σに合わせる。nhits vs edepは両軸の下限を0に保ち、上限だけ各平均＋5σの1.5倍に合わせる。
範囲は元の軸内に収め、fitが無効な場合は元の表示範囲を使う。PNGと保存summary canvasの描画コピーだけを変更し、
元のヒストグラム・fit・event選択は変更しない。

Slurmの標準出力・標準エラーは `analysis/log/resolution/` に保存します。`singleEnergy/`、`compare/`、`hl_tail/<条件>/<工程>/` などの区分はその配下で維持します。実行中・待機中のジョブが使う旧 `jobs/` のログは、ジョブ終了後に移動します。

## 結果のディレクトリ構造

`analysis/result/resolution/` を基点とし、SingleEnergyは
**single → data / mc → sample → selection → CoG**、MultiEnergyは
**multi → 比較種別 → 比較名 → selection** の順で保存する。MCの名称は `mc` に統一する。
図は `figures/` 以下に同じ相対階層を使う。

```text
analysis/result/resolution/
  single/data/<sample>/<selection>/<CoG>mm/<energy>GeV_<CoG>mm.root
  single/mc/<sample>/<selection>/<CoG>mm/<energy>GeV_<CoG>mm.root
  multi/<data_vs_mc|mc_vs_mc|data_vs_data>/<比較名>/<selection>/comparison.root
  figures/single/data/<sample>/<selection>/<CoG>mm/<energy>GeV/...
  figures/single/mc/<sample>/<selection>/<CoG>mm/<energy>GeV/...
  figures/multi/<比較種別>/<比較名>/<selection>/...
  unknown/<旧相対パス>
  unknown/inputs/layout_migration_20261112/
```

現役設定のsampleはdataが `beam2023_reprocessed`、MCが `e_ssa_threshold` と
`spread_ps0p3_sps2`。選別名はdata/MCで共通の `tail_exclude_front_off`、
`tail_keep_front_off`、`tail_exclude_front20mm` 等とし、入力sampleとevent選別を分ける。
CoG 200 mmは保存先では `200mm`、図では `nocut` と表示する。

由来を確定できない旧data・`simulation/threshold`・strip結果・`comparsion`等は、
`unknown/` 以下に旧相対パスを保って移した。旧図も `unknown/figures/<旧相対パス>` に保持し、
新しいsampleと統合しない。既定のdata/MC比較は、移行前と同じdataを
`unknown/data/` から読み、比較結果も `unknown/custom/momentum_spread/comparison/` 以下へ保存する。
再処理dataを比較に使う場合は、後述の `MOMENTUM_SPREAD_DATA_RESULTS` で明示する。

移動対応は [moves.tsv](../result/resolution/unknown/inputs/layout_migration_20261112/moves.tsv)、
全ファイル一覧は [files.tsv.gz](../result/resolution/unknown/inputs/layout_migration_20261112/files.tsv.gz)、
実施記録は [run.log](../result/resolution/unknown/inputs/layout_migration_20261112/run.log) に保存する。
既存の `run.log` とROOT内部の入力パスは、当時の実行記録として書き換えない。

新規のKEKCC実行には `momentum_spread_samples.sh`、任意条件には更新済みの
`samples.example.sh` / `strip_samples.example.sh` を使う。
`analysis_paralell*.sh`、`energy_resolution*.sh`、`compare_ds*.sh`、
`hl_intercept_samples.sh`、`hl_tail_samples.sh`、`sim_sps_e_ssa_samples.sh`、
`data_sim_e_ssa_samples.sh` は過去の解析用入口として保持し、結果参照先は `unknown/` 以下とする。
これらには旧サーバーの入力パスやSlurm設定が残り、必要な入力・保存済み結果の存在も保証しない。
新規解析へそのまま流用せず、以下の歴史的な設定説明と現役のmomentum設定を区別する。

## 前方2層の位置によるevent cut

SingleEnergyAnalysisに `--front-position-range MM` を追加すると、中心(0,0)からの半幅を
Rx=Ry=MMとして、`abs(x_front) < MM && abs(y_front) < MM` を要求する。
`x_front` はlayer 1のHit_X、`y_front` はlayer 0のHit_Yをenergyで重み付けした平均。
各層の短辺方向だけを使い、正で有限のenergyと有限の座標を持つhitを、既存の `only_best` 設定に従うchannel mask適用後に集計する。
両層の位置が有効でないeventは除外する。境界は含まず、各層1hitの条件は設けない。
この位置は異なるzにある2層の重心による推定値で、trackの外挿位置ではない。

省略時は無効で、現行CoG 200/20/10/5 mmの4条件と既存の結果を維持する。
指定時はHG–LG tail除外後・CoG cut前に追加のevent選択を行い、通過したeventのenergy/Nhitを従来どおり集計する。
`--strip-range` と併用する場合も、前方位置はstripによるhit選択前の全good hitから求める。
既存の `position_cut_mode`、CoG/strip設定はそのまま保存し、前方cutは別metadataへ記録する。

```bash
cd /gpfs/group/ilc/users/murata/ScECAL_BeamTest/analysis/resolution
make SingleEnergyAnalysis
# dataだけ：19 energy × 4 CoGの76ジョブをLSFへ投入
MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh --front-position-range 20 --exclude-tail-events single momentum_spread_samples.sh
```

10 mm等に変更する場合は20の値を変更する。環境変数 `FRONT_POSITION_RANGE_MM=20` でも指定できる。
C++を直接実行する場合は、既存の位置引数の後ろに `--front-position-range 20` を付ける。
同じ設定ファイルを使う場合、前方cut付きROOTは
`analysis/result/resolution/single/data/beam2023_reprocessed/tail_exclude_front20mm/<CoG>mm/`、
図は対応する `analysis/result/resolution/figures/single/data/beam2023_reprocessed/tail_exclude_front20mm/<CoG>mm/` へ保存する。
実行記録はROOT基点の `single/data/beam2023_reprocessed/tail_exclude_front20mm/inputs/single.run.log`、LSFログは
`analysis/log/resolution/momentum_spread/data/exclude_front20mm/single/%J.out` / `.err`。
前方cutなしは同じsampleの `tail_exclude_front_off/` に保存する。MCモードでも
`single/mc/{e_ssa_threshold,spread_ps0p3_sps2}/tail_exclude_front20mm/` に分ける。
Multiでは同じオプションで対応する保存済みROOTを選ぶ。data/MC比較には両方のfront cut付きSingle結果が必要で、Multiで新たにeventを選別することはない。
前方cut指定時のdata/MC比較は、この設定で生成した同じ半幅のdata保存先を使用する。
`MOMENTUM_SPREAD_DATA_RESULTS` が別の保存先を指す場合は停止する。未指定なら自動で対応する保存先を選ぶ。
別のsample設定へこのオプションを使う場合は、その設定でもcut条件別の出力先を指定する。

dataのfront cutなしとfront20ありを比較する場合は、両方のSingleEnergy完了後に次を実行する。
`data_front` はMulti専用で、`single/data/beam2023_reprocessed/` 以下の
`{tail_exclude_front_off,tail_exclude_front20mm}/<CoG>mm/` の保存済みROOTを読み、4つのCoGそれぞれで2条件を比較する。

```bash
MOMENTUM_SPREAD_COMPARE=data_front \
bash run_samples.sh --exclude-tail-events --front-position-range 20 multi momentum_spread_samples.sh
```

LSFへ1ジョブを投入し、比較ROOTは `multi/data_vs_data/front_off_vs_front20/tail_exclude/comparison.root`、
図は `figures/multi/data_vs_data/front_off_vs_front20/tail_exclude/`（いずれも `analysis/result/resolution/` 基準）へ保存する。
全条件の図は `summary_resolution.png`、CoG別の図は `summary_resolution_{200,20,10,5}mm.png`。
SingleEnergyの再実行は不要。
実行記録はROOT基点の `multi/data_vs_data/front_off_vs_front20/tail_exclude/inputs/multi.run.log`、
LSFログは `analysis/log/resolution/momentum_spread/data_front/exclude_front20mm/multi/%J.out` / `.err`。
他の半幅を比較する場合は20を変更し、その半幅のSingleEnergy結果を先に生成する。
`MOMENTUM_SPREAD_COMPARE=data` は従来どおりdata対momentum-spread MCの比較に使う。

ROOTの `front_position_cut_enabled`、`front_position_range_mm`、`front_position_definition` に条件を保存する。
`front_position_checked_events` / `front_position_invalid_events` / `front_position_rejected_events` /
`front_position_passed_events` はtail除外後・CoG前の件数で、rejectedはinvalidを含み、checked=rejected+passed。
disabled時は半幅−1、各件数0。`hl_tail_events_after_cog` は追加した前方位置cutも通過した後の件数になる。
既存の `firstLayerHit_vs_cog.png` は前方位置cut・CoG cut前の診断図として維持する。

## hit strip 座標で集計範囲を選ぶ

SingleEnergyAnalysis に `--strip-range MM` を付けると、各 hit の校正済み strip 中心座標で
`abs(Hit_X) < MM && abs(Hit_Y) < MM` を要求する。MM は x/y 方向の半幅 [mm] で、
22.5 なら中央45×45 mmの範囲。境界は含まず、layer 0–29 の範囲内の hit だけを集計する。
この指定時は layer 9/10 の CoG によるイベント選択を行わない。
総energy、Nhit、層別energy/Nhit、1hit energy、channel別分布は同じ strip 条件を使う。
CoG図は選択前の全good hitから計算した重心を診断用に表示する。
選択範囲にhitがないイベントも `beforeEventCut` ではenergy/Nhitが0のイベントとして保持する。
その範囲の集計値とfitを使う従来の `afterEventCut` のenergy–Nhit選別は、現在停止している。
bad-channelとHG–LG tailの指定は従来どおり適用する。

shellから選ぶ場合は、[strip_samples.example.sh](strip_samples.example.sh) をコピーする。
`strip_ranges_mm=(22.5 68.2)` を編集すると、各領域のSingleEnergyと領域間比較の対象を選べる。

```bash
cd analysis/resolution
make SingleEnergyAnalysis
cp strip_samples.example.sh strip_samples.sh
# strip_samples.sh の strip_ranges_mm、入力、保存先、runnerを編集
bash run_samples.sh --dry-run single strip_samples.sh
bash run_samples.sh single strip_samples.sh
# 全SingleEnergy jobの正常終了後:
bash run_samples.sh multi strip_samples.sh

# 1領域のSingleEnergyだけを実行する場合:
bash run_samples.sh --strip-range 22.5 single strip_samples.example.sh
```

設定ファイル内では `strip_range_mm=22.5` を `single_sample` の前に指定してもよい。
stripモードでは `single_sample` のCoG引数を200とし、保存先に
`tail_exclude_front_off_strip22.5mm/200mm/` のように領域条件を含める。`*_200mm.root` の末尾は既存MultiEnergyとの読込互換のため保持する。
このモードでのCoG event cutは無効であり、ROOTの `cog_range_mm` は0となる。
MultiEnergyでは `multi_sample 'Data, strip |X|,|Y| < 22.5 mm' ...` のように凡例を指定する。
設定例は `single/data/beam2023_reprocessed/tail_exclude_front_off_strip<半幅>mm/200mm/` と、
それに対応する `figures/single/` に保存する。領域間比較は
`multi/data_vs_data/strip22.5mm_vs_strip68.2mm/tail_exclude_front_off/comparison.root` と、対応する図の保存先を使う。
22.5と68.2の範囲は入れ子で、周辺だけを分離した領域ではない。

過去のSlurm用shellにも次の指定がある。CoGの4条件ループは1回にし、ROOTは
`unknown/data/channel22.5mm/200mm/` または `unknown/mc/<threshold>/channel22.5mm/200mm/`、
図は対応する `unknown/figures/` 以下に分ける。新規実行は上記の設定例を使う。

```bash
STRIP_RANGE_MM=22.5 bash analysis_paralell_allCondition.sh
```

C++を直接呼ぶ場合は位置引数の後ろに指定する。

```bash
./SingleEnergyAnalysis output.root 1 100GeV 1 calibration.root 200 1 figures \
  --strip-range 22.5
```

`--skip-channel-plots` はSingleEnergyでも使用でき、channel別の大量のPNG/canvasだけを省略する。
`run_samples.sh` では設定の `single_options=(--skip-channel-plots)` で渡せる。
stripの設定例ではこの指定を有効にしている。総energy、Nhit、層別分布、fit、summaryは保存する。
stripモードではfit初期値を領域内の分布から求め、全検出器向け200 GeVの固定初期値は使わない。
energy/Nhitのfit下限は0にし、狭い領域の小さい信号もfit範囲に含める。
ROOTには `position_cut_mode` (`strip` / `cog`)、`strip_range_mm`（CoGモードでは−1）、
`cog_range_mm`、`channel_plots_saved` を保存する。
`--strip-range` を省略した場合は従来のCoGモードとなる。

## 任意の calibration / 条件を shell で選ぶ

`calibration ROOT → SingleEnergyAnalysis → 条件別の解析 ROOT → MultiEnergyAnalysis`
の流れは従来と同じ。SingleEnergy の出力は cut 済み event tree ではなく、分布・fit・canvas を保存した解析 ROOT。
SingleEnergy の CoG cut は layer 9 と 10 の両方で `abs(x) < range && abs(y) < range` を要求する。
threshold を変更する引数はなく、その条件は入力 calibration sample 側で選ぶ。
bad-channel 除外とGaussian fit の処理は維持する。energy–Nhit の楕円 event cut は現在適用しない。
これらに先立ち、既定では HG–LG の tail-rich event を除外する（下記）。

```bash
cd analysis/resolution
make SingleEnergyAnalysis MultiEnergyAnalysis
cp samples.example.sh my_correction_samples.sh
# my_correction_samples.sh の calibration 入力、保存先、比較するファイル、凡例を編集
bash run_samples.sh --dry-run single my_correction_samples.sh
bash run_samples.sh single my_correction_samples.sh
# SingleEnergy の完了後（Slurm 使用時は全 job の完了を待つ）:
bash run_samples.sh --dry-run multi my_correction_samples.sh
bash run_samples.sh multi my_correction_samples.sh
```

`my_correction_samples.sh` は新しい補正を比較するときの設定名の例。内容が分かる名前に変更する。
相対パスの基準は常に `analysis/resolution/`。
HL intercept用の `hl_intercept_samples.sh`、tail比較の `hl_tail_samples.sh` と
`analysis_paralell*.sh` / `compare_ds*.sh` は、冒頭に示した過去の解析用入口として保持する。

### SingleEnergy の選択

```bash
single_sample ROOT保存先 図の保存先 CoG範囲mm bad_channel除外 FILE_GLOB calibration_dir...

# 入力は calibration_dir/10GeV/*.root など。最後に複数の入力ディレクトリも指定可能。
single_sample '../result/resolution/single/data/calibration_A/tail_exclude_front_off/20mm' \
  '../result/resolution/figures/single/data/calibration_A/tail_exclude_front_off/20mm' \
  20 1 'ECAL*.root' '/path/to/calib_A/ps' '/path/to/calib_A/sps'
```

- `bad_channel除外` は `1` で既存の bad-channel 一覧を除外、`0` で除外しない。
- 同じ energy の calibration ファイルをまとめて1つの SingleEnergy に渡す。
  ファイル数は `FILE_GLOB` に実際に一致したファイルだけを数える。
- ROOT は `ROOT保存先/10GeV_20mm.root`、図は `図の保存先/10GeV/` に保存する。
- calibration / CoG 条件ごとに別の保存先を指定する。再実行時は同じ出力を更新する。
  同一実行内で ROOT 保存先が重複した場合や同じ入力が重複した場合はエラー。
- 既存 SingleEnergy に合わせ、energy ディレクトリは正の整数 GeV または `0.5GeV`。
- 標準はローカル実行。`single_runner=(sbatch ... execute_paralell.sl)` で既存 Slurm wrapper も使える。

### MultiEnergy の選択

```bash
multi_sample 'Calibration A, CoG 20 mm' '../result/resolution/single/data/calibration_A/tail_exclude_front_off/20mm/'*.root
multi_sample 'Calibration B, CoG 20 mm' '../result/resolution/single/data/calibration_B/tail_exclude_front_off/20mm/'*.root
multi_sample 'MC alternative, CoG 10 mm' '../result/resolution/single/mc/calibration_C/tail_exclude_front_off/10mm/'*.root

multi_output='../result/resolution/multi/data_vs_mc/calibration_A_vs_B_vs_MC_C/tail_exclude_front_off/comparison.root'
multi_figures='../result/resolution/figures/multi/data_vs_mc/calibration_A_vs_B_vs_MC_C/tail_exclude_front_off'
# channel 別の図は標準で省略。必要な場合だけ、ほかのオプション追加前に指定:
# multi_options=()
```

- 2系列以上を指定できる。各系列のファイルは glob または明示したリストで選択する。
- 各系列は同じエネルギー集合を指定する。ファイル名先頭の数値で並べ替えて対応付けるため、入力順は自由。
  エネルギーの欠落・重複、入力ファイルや必要な histogram / fit の欠落はエラー。
- 凡例には指定ラベルとCoGを表示する。ラベルは `MC after` / `MC before` など条件名だけでよい。
  CoGは入力名 `*_20mm.root` 等から読む。200 mmは表示上 `nocut` に統一し、保存ディレクトリ名は `200mm` のまま。
- CoG別に色・marker形状を固定する。nocut=青/丸、20 mm=紫/菱形、10 mm=緑/三角、5 mm=赤/四角。
  同じ条件のCoG違いは同じ線種。最初の条件は実線/塗りつぶし、次は点線/白抜き。
  条件は入力の `<条件>/<CoG>mm/` の `<条件>` までのパスで識別する。系列の交互登録・条件ごとの登録の両方に対応する。
  3条件以上は追加の線種を使う。CoGを含まない任意名の入力は従来の系列別配色を維持する。
- resolutionは全系列を重ねる。6系列以上では図の右側に凡例専用の領域を設ける。既存の `resolution.png` / `resolution_with_noise.png` に加え、
  `resolution_loglog.png` / `resolution_with_noise_loglog.png` を保存する。ROOTにも同名（拡張子なし）のcanvasを保存。
  linear版は常にx=0–130 GeV、y=0–0.3。log–log版は全入力の正の点を表示できる軸範囲を使う。
  4分割の従来表示でも、すべての描画panelをlog–logにする。表示だけを変更し、点・fit値・fit範囲は変えない。
  保存用のfit曲線は求めた係数から解析式として作成し、除外energyで描画上の落ち込みが生じることを防ぐ。
  energy summaryとchannel分布にも同じ色・線種・凡例を使用する。
- MultiEnergyの全canvasは `new TCanvas(..., 2560, 1440)` で作成し、`SetCanvasSize(2560, 1440)` で描画領域も固定する。PNG・ROOTともWQHDに統一し、summaryや多系列比較でも縦横のサイズを変えない。
- resolution単独4枚のmarkerはサイズ3.0とする。
- resolution summaryは全条件をまとめたものと、CoG cutごとに分けたものの両方を保存する。
  全条件は `summary_resolution.png` とROOTの同名canvasに保存する。
  CoG別は同じCoGの比較条件だけを重ね、
  `summary_resolution_{200,20,10,5}mm.png` とROOTの同名canvasに保存し、200 mmは図・凡例では `nocut` と表示する。
  `Divide(2,2)` で、
  上段linear・下段log–log、左列without noise・右列with noise。各パネルのタイトルは左列を `without noise term`、右列を `with noise term` とし、軸形式はタイトルに付けない。凡例は上段だけに表示する。
  各panelは対応する単独図の点・fit係数をそのまま使用する。
  CoGを識別できない任意名の入力も全条件の `summary_resolution.png` に含める。
  全系列を重ねるresolution単独4枚も維持する。
- CoG別の比較図 `resolution_cog_comparison.png` と `mean_linearity_cog_comparison.png` を追加保存し、
  比較ROOTにも拡張子なしの同名canvasを保存する。いずれも2560×1440の `Divide(2,2)` で、
  左上nocut（200 mm）、右上20 mm、左下10 mm、右下5 mm。各セルを上下に分ける。
- `resolution_cog_comparison` の上段は既存のwithout noise図と同じ分解能の点・fitを使用し、
  x=0–130 GeV、y=0–0.3。下段は `100 × (系列 / 参照系列 − 1)` [%]。
  同じenergyの点だけを対応付け、参照点がないenergyは省略する。参照系列自身は0の基準線とする。
  入力パスの要素が `mc` または `simulation` の系列が同じCoG内に1つだけある場合は、そのMCを参照する。
  それ以外（dataの前後比較やMCが2条件ある場合など）は、そのCoGで最初に登録した条件を参照する。
  比の誤差は系列同士を独立とした近似で、共有eventによる共分散は含めない。
- `mean_linearity_cog_comparison` の上段はSingleEnergyの `beforeEventCut/fit_gaus` のmeanと誤差を
  MeVからGeVへ換算したもの。分解能と同じ有効点判定を使い、dataの200 GeVは除外する。
  各系列を切片自由の `pol1`（μ=a+bE）でfitし、150/250 GeVや30/80 GeVも含む全有効energyを使用する。
  分解能fitの30/80 GeV除外設定はこの直線fitには適用せず、横軸は全energyを表示する。
  下段は `100 × (mean − fit) / mean` [%]、誤差はfitを固定した
  `100 × abs(fit) × errMean / mean²`。fit係数の不確かさや共分散は誤差に含めない。
  有効点が3点未満、またはfit失敗の場合もmeanの点は残し、直線と残差は描かない。
- 下段の縦軸は差と誤差から自動設定するが、大きい低energyの誤差で通常の差が潰れないよう、
  resolution比較は最大±100%、meanのlinearityは最大±10%とする。
  点の中心値が範囲外なら矢印と件数を表示し、全点・全誤差はROOTのgraphに保持する。
- 新しい2図にCLIの追加はなく、同じMultiEnergyコマンドで生成する。保存済みSingleEnergy ROOTを使えるため、
  SingleEnergyの再実行は不要。既存のresolution単独図・全条件summary・CoG別summaryも保存する。
  CoGを識別できない入力だけの場合、新しい図では全系列を最初のセルにまとめる。
  識別できるCoGと混在する場合、新しい4 CoG図には識別できる系列だけを載せる。
- `summary_<energy>GeV.png` のenergy deposit / number of hitsは、元の描画領域の広さを維持し、凡例をヒストグラム内の右寄りに重ねる。
  凡例の線見本を長くし、これらの分布の線幅を2にして条件を読み取りやすくする。
- channel 分布を有効にした場合は、選択した全系列を重ね、`multi_figures/samples/<energy>GeV/raw/` に保存する。
  比較 ROOT 内は `samples/<energy>GeV/raw/`。共通 rebin と正規化は比較する全系列から決める。
- `run_samples.sh` は標準で `--skip-channel-plots` を付け、channel 別の PNG と ROOT canvas を省略する。resolution と energy summary は通常どおり出力する。
  channel 比較も必要な場合だけ、設定の系列登録・`multi_run` より前で `multi_options=()` とする。
  その他のオプションは `multi_options+=(...)` で追加すれば、標準の省略設定を維持できる。
  C++ 実行ファイルを直接呼ぶ場合は、省略するには `--skip-channel-plots` を明示する。
- resolution は従来と同じ `beforeEventCut/fit_gaus` の **σ/μ**。CoG cut と指定した bad-channel 除外は適用済み。
  エネルギー依存 fit で 25–35 GeV / 75–85 GeV を除外する既存設定も維持している（点自体は表示する）。

### 2条件 × CoG 4条件を同時に比較する

`samples.example.sh` は `(200 20 10 5)` × after/beforeの8系列を1つの比較に登録する。
SingleEnergyは `single/<対象>/<sample>/<selection>/<CoG>mm/<energy>GeV_<CoG>mm.root`、
図も `figures/single/` 以下に同じ階層で保存する。
`sample_type` は入力に合わせて `data` または `mc` に設定し、after/beforeは異なるsample名で区別する。

```bash
for cog in "${cog_ranges[@]}"; do
  multi_sample "After" "$result_base/$sample_type/my_correction_after/$selection/${cog}mm/"*.root
  multi_sample "Before" "$result_base/$sample_type/my_correction_before/$selection/${cog}mm/"*.root
done
multi_run "../result/resolution/multi/$comparison_type/my_correction_after_vs_before/$selection/comparison.root" \
  "../result/resolution/figures/multi/$comparison_type/my_correction_after_vs_before/$selection"
```

`multi_run` は全CoGの登録後に1回呼ぶ。第1引数が比較ROOT、第2引数が図の保存先。
CoG別の個別比較も必要なら、各CoGの2系列を登録するループ内で `multi_run` を呼ぶ。
実行後は系列リストを空にする。出力の重複はエラー。`multi_output` / `multi_figures` による暗黙の1比較も利用可能。
過去の `hl_intercept_samples.sh` はSingle/MultiともSlurmを使い、multiは8系列をまとめて1ジョブ投入する。
参照先は `unknown/custom/hl_intercept/mc/{after,before}/<CoG>mm/`、比較ROOTは
`unknown/custom/hl_intercept/mc/comparison/all.root`、図は `unknown/figures/custom/hl_intercept/` 以下。
この設定の旧 `after/<CoG>mm/`・`before/<CoG>mm/` へのfallbackも、同じunknown基点内だけで行う。

```bash
# 既存のHL intercept比較はtail除外前のROOTなのでkeepを明示する。
bash run_samples.sh --keep-tail-events --dry-run multi hl_intercept_samples.sh
bash run_samples.sh --keep-tail-events multi hl_intercept_samples.sh
```

### Momentum spread あり / なしのMC比較（KEKCC）

`momentum_spread_samples.sh` は校正済みの旧SSA MC `calib/e-/e_ssa/threshold` と
新MC `calib/e-/momentum_spread/default` をSingleEnergyから解析する。
入力は `data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/` 以下。
19 energy × 2 MC × 4 CoG（200/20/10/5 mm）でSingleEnergyは152ジョブ、
MultiEnergyは全8系列をまとめる1ジョブをLSFへ投入する。キューは既定 `s`、`LSF_QUEUE` で変更できる。
既定では両MCのtail-rich eventを除外する。旧MCの判定元は
`analysis/result/channel_response/hg_lg_tail/threshold/e_ssa_kekcc/`、新MCは
`analysis/result/channel_response/hg_lg_tail/momentum_spread/default/`。
各energyの `tail_study.root` と `inputs/{data,mc}.tsv` を全eventについて先に生成する。
旧 `threshold/e_ssa/` のmanifestは旧サーバーの絶対パスを記録しているため、そのまま再利用せず、
KEKCCの入力を使い別labelへ再生成する。判定結果が不足している場合は、先に以下の38ジョブをLSFへ投入する。

```bash
cd /gpfs/group/ilc/users/murata/ScECAL_BeamTest/analysis/channel_response/hg_lg_tail
bash submit.sh --queue l --energy all --mc-sample e_ssa/threshold \
  --mc-tag threshold --label e_ssa_kekcc
bash submit.sh --queue l --energy all --mc-sample momentum_spread/default \
  --mc-tag momentum_spread --label default
bjobs
```

`run.sh`はKEKCCの入力・定数パスを既定で使い、`--mc-sample`からenergyごとのdecode/calibを選ぶ。
全入力・全eventでの正常終了が必要なため、サンプリングオプションは付けない。
不足時はSingleEnergyの投入前に停止する。
判定元は設定の `tail_results_old` / `tail_results_spread` で指定する。dry-runはコマンドの確認だけで、判定結果の存在・完全性は検証しない。
bad-channel除外は1、channel別の大量の図は省略する。HG/LG補正は入力calibrationに含まれる。

```bash
cd /gpfs/group/ilc/users/murata/ScECAL_BeamTest/analysis/resolution
make SingleEnergyAnalysis MultiEnergyAnalysis
bash run_samples.sh --dry-run single momentum_spread_samples.sh
bash run_samples.sh single momentum_spread_samples.sh
bjobs
# 全SingleEnergyジョブの正常終了後:
bash run_samples.sh --dry-run multi momentum_spread_samples.sh
bash run_samples.sh multi momentum_spread_samples.sh
```

ROOTは `single/mc/{e_ssa_threshold,spread_ps0p3_sps2}/tail_exclude_front_off/<CoG>mm/`、
比較ROOTは `multi/mc_vs_mc/e_ssa_threshold_vs_spread_ps0p3_sps2/tail_exclude_front_off/comparison.root`
（`analysis/result/resolution/` 基準）。図は `figures/single/` / `figures/multi/` 以下の対応する階層に保存する。
Singleの共同実行記録は `single/mc/inputs/momentum_spread/tail_exclude_front_off/single.run.log`、
Multiの記録は比較ROOTと同じディレクトリの `inputs/multi.run.log`。
LSFログは `analysis/log/resolution/momentum_spread/exclude/{single,multi}/%J.out` / `.err`。
`--keep-tail-events` を明示した場合のみselectionを `tail_keep_front_off` とし、LSFログの条件名も `keep` に分ける。
linear / log–logのresolution、summaryを保存し、既存のfit・event選択規則を使う。

dataとmomentum spread MCの比較には、同じ設定を `MOMENTUM_SPREAD_COMPARE=data` で実行する。
既定では、移行前と同じtail除外済みdataの `analysis/result/resolution/unknown/data/<CoG>mm/` と
`single/mc/spread_ps0p3_sps2/tail_exclude_front_off/<CoG>mm/` を再利用する。MultiEnergyは
4 CoG × 2条件の8系列を1ジョブにまとめる。現在は両方のSingleEnergy結果が全19 energy分あるため、SingleEnergyの再実行は不要。
dataの200 GeVはsummaryだけを描画し、分解能点・fitには使わない。その他の無効なenergy fitも現在のMultiEnergyの規則に従いsummaryのみになる。

```bash
cd /gpfs/group/ilc/users/murata/ScECAL_BeamTest/analysis/resolution
make MultiEnergyAnalysis
MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh --exclude-tail-events --dry-run multi momentum_spread_samples.sh
MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh --exclude-tail-events multi momentum_spread_samples.sh
bjobs
```

既定の比較ROOTは `analysis/result/resolution/unknown/custom/momentum_spread/comparison/data_exclude.root`、
図は `analysis/result/resolution/unknown/figures/custom/momentum_spread/comparison/data_exclude/`。
実行設定記録は `unknown/custom/momentum_spread/comparison/inputs/exclude/multi.run.log`（ROOT基点）、
LSFログは `analysis/log/resolution/momentum_spread/data/exclude/multi/%J.out` / `.err` に保存する。

dataのSingleEnergy図を現在の表示設定で再出力する場合も、同じ設定でdataだけを投入できる。
PS/SPSのKEKCC上のcalibをenergyごとにまとめ、19 energy × 4 CoGの76ジョブを投入する。
HG–LG tailは `threshold/e_ssa_kekcc` の生成済み判定で除外し、bad-channel除外は1、channel別の大量の図は省略する。

```bash
cd /gpfs/group/ilc/users/murata/ScECAL_BeamTest/analysis/resolution
make SingleEnergyAnalysis
MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh --exclude-tail-events --dry-run single momentum_spread_samples.sh
MOMENTUM_SPREAD_COMPARE=data bash run_samples.sh --exclude-tail-events single momentum_spread_samples.sh
bjobs
```

ROOTは `analysis/result/resolution/single/data/beam2023_reprocessed/tail_exclude_front_off/<CoG>mm/`、
図は `analysis/result/resolution/figures/single/data/beam2023_reprocessed/tail_exclude_front_off/<CoG>mm/` に保存する。
実行設定記録はROOT基点の `single/data/beam2023_reprocessed/tail_exclude_front_off/inputs/single.run.log`、LSFログは
`analysis/log/resolution/momentum_spread/data/exclude/single/%J.out` / `.err`。
この再生成したdataをmomentum spread MCとのMulti比較に使う場合は、全Singleジョブの正常終了後に次を実行する。
指定しない場合、Multiは移行前と同じdataを `analysis/result/resolution/unknown/data/` から読む。

```bash
MOMENTUM_SPREAD_COMPARE=data \
MOMENTUM_SPREAD_DATA_RESULTS=../result/resolution/single/data/beam2023_reprocessed/tail_exclude_front_off \
bash run_samples.sh --exclude-tail-events multi momentum_spread_samples.sh
```

この明示指定では比較ROOTを
`multi/data_vs_mc/beam2023_reprocessed_vs_spread_ps0p3_sps2/tail_exclude_front_off/comparison.root`、
図を対応する `figures/multi/` 以下、記録を比較ROOTと同じディレクトリの `inputs/multi.run.log` に保存する。

### Gaussian beam / SSA hitmap の既存MC比較（過去の入口）

`sim_sps_e_ssa_samples.sh` はMultiEnergy専用で、保存済みの
`unknown/default/simulation/threshold/<CoG>mm/`（凡例 `oval`）と
`unknown/custom/hl_intercept/mc/before/<CoG>mm/`（凡例 `SSA hitmap`）を比較する。
CoGは200/20/10/5 mmの8系列。200 mmは凡例で `nocut` と表示する。
両方ともtail除外導入前の結果なので、`--keep-tail-events` が必須。
この指定は保存済みの選択条件を検査するもので、SingleEnergyのevent選択を変更しない。

```bash
bash run_samples.sh --keep-tail-events --dry-run multi sim_sps_e_ssa_samples.sh
bash run_samples.sh --keep-tail-events multi sim_sps_e_ssa_samples.sh
```

ローカル実行。ROOTと入力記録は
`analysis/result/resolution/unknown/default/comparison/sim_sps_e-_ssa/{all.root,run.log}`、
図は `analysis/result/resolution/unknown/figures/default/comparison/sim_sps_e-_ssa/` に保存する。
`oval` は実線・塗りつぶし、`SSA hitmap` は点線・白抜き。CoGごとの色・marker形状は共通。

### Data / SSA hitmap MC before の既存結果比較（過去の入口）

`data_sim_e_ssa_samples.sh` は `unknown/default/data/<CoG>mm/`（凡例 `data`）と
`unknown/custom/hl_intercept/mc/before/<CoG>mm/`（凡例 `sim`）の生成済みSingleEnergy結果を比較する。
4 CoG × 2系列をローカル実行し、SingleEnergyは再計算しない。
0.5 GeVはdataの20 mmの保存fitが無効なため、全8系列から除外する（共通の18 energy、1〜250 GeV）。

```bash
bash run_samples.sh --keep-tail-events multi data_sim_e_ssa_samples.sh
```

ROOTと入力記録は `analysis/result/resolution/unknown/default/comparison/data_sim_e-_ssa/`、
図は `analysis/result/resolution/unknown/figures/default/comparison/data_sim_e-_ssa/` に保存する。
`data` は実線・塗りつぶし、`sim` は点線・白抜き。CoGの色と表記は共通。

C++ を直接呼ぶ場合も、従来の位置引数の後ろに凡例を追加できる。

```bash
./MultiEnergyAnalysis comparison.root 3 2 \
  A/10GeV_20mm.root A/20GeV_20mm.root \
  B/10GeV_20mm.root B/20GeV_20mm.root \
  C/10GeV_10mm.root C/20GeV_10mm.root \
  figures --labels 'Calibration A' 'Calibration B' 'Alternative MC'
```

## 動作確認

ROOT と同じ Python バージョンの PyROOT を使用する。

```bash
python3 tests/test_sample_comparison.py
# 全360枚の channel 図の保存まで検証する場合:
python3 tests/test_sample_comparison.py --channels
```

合成 ROOT を用いて3系列の凡例・色・resolution 値、並べ替え、小数エネルギー、
不整合の検出、従来 CLI、SingleEnergy の shell 引数を検証する。

## 従来の data / simulation 比較

```bash
make MultiEnergyAnalysis
bash compare_ds.sh               # 現在は CoG 200 mm の data / simulation
bash compare_ds_allCondition.sh   # 現在は CoG 200 / 20 / 10 / 5 mm
```

`--labels` を指定しない場合の位置引数・4分割構成は維持する。CoG配色、nocut表示、軸範囲とlog–log追加保存は共通。
過去の tail 除去なし ROOT を直接読む場合は `--keep-tail-events` を付ける。

```text
MultiEnergyAnalysis output.root N_datasim N_files input_1.root ... figure_path
```

## 従来モードの channel ごとの edep

SingleEnergy の ROOT に保存された `<energy>GeV/edep_Layer<L>_Chip<C>`
canvas から `edep_channel_<L>_<C>_<channel>` を読み、**同じ CoG 条件の data と simulation の2系列だけ**を重ねる。
入力順に依存せず条件を照合する。simulation の threshold が複数ある場合は、それぞれ別の図にする。
calibration の再読込や SingleEnergy の再実行は不要。
小数エネルギーでは、旧 SingleEnergy が使う整数のディレクトリ名（例: `0GeV`）も読む。

- 1枚の canvas は Single の raw と同じ6×6配置（channel 0〜35）。
- layer 0〜29 × chip 0〜5、現在は event cut 前のみを保存する（1 energy・1条件あたり180枚）。
  CoG が4条件なら、各 chip につき4枚を別ディレクトリに保存する。
- 各 channel の表示範囲内（0〜2 MeV）の積分を1にして分布の形を比較する。
- rebin factor は channel ごとに統計量と四分位範囲から自動選択する。
  同じ条件・channel の data / simulation に共通の factor を使い、表示は40〜200 bins に収める。
  別条件の統計量はその図の rebin factor に影響しない。
  20 entries 未満の系列は粗い bin 幅を選び、空の系列は幅の推定に使わない。
  元の bin 数を割り切る factor のみを選ぶため、余りを overflow に移さない。
  各パネルに factor を表示し、rebin の後に正規化する。
- 空のヒストグラムは割り算せずに残し、すべて空の channel は `no entries` と表示する。
- 縦軸の上限は同条件の2系列から決める。凡例は canvas 上部に置き、data / simulation・CoG・threshold を示す。
- 入力 ROOT は変更しない。元の raw と同じく、bad-channel 除外前に集計された channel 分布を比較する。
- 必要な raw canvas / histogram がない場合は、対象ファイルと名前を表示して終了する。
- 対になる data / simulation がない条件や、同じ条件の重複入力は、異なる条件を混ぜずにエラーとして知らせる。

PNG は指定した `figure_path/cog<cog_range>mm/<threshold>/<energy>GeV/raw/` に保存する。
**PNG のファイル名は従来のまま**で、条件による区別はディレクトリだけで行う。
1条件だけを比較する場合も同じ構造を使う。

```text
cog200mm/threshold/30GeV/raw/30GeV_layer0chip0.png
cog20mm/threshold/30GeV/raw/30GeV_layer0chip0.png
cog10mm/threshold/30GeV/raw/30GeV_layer0chip0.png
cog5mm/threshold/30GeV/raw/30GeV_layer0chip0.png
```

既存の `cut_30GeV_layer0chip0.png` は従来の event cut 後の図であり、現在は生成しない。
既存 shell の `figure_path` は `../result/resolution/unknown/figures/comparsion/200mm/` または
`../result/resolution/unknown/figures/comparsion/dataThre/`（既存の `comparsion` 表記を維持）。
比較 ROOT にも `cog<cog_range>mm/<threshold>/<energy>GeV/raw/edep_Layer<L>_Chip<C>`
として canvas を保存する。0.5 GeV なども出力では正確な energy 名を使う。
従来の energy summary・resolution の全条件比較図はそのまま出力する。


## HG–LG tail event を除いた resolution

SingleEnergy は既定で tail-rich event を除外してから、従来の CoG cut、
エネルギー分布作成、Gaussian fit を行う。MultiEnergy の resolution は、その
`beforeEventCut/fit_gaus` の **σ/μ**。hg_lg_tail の診断用 RMS/mean とは異なる。
tail 除外は `beforeEventCut` に適用する。従来は `afterEventCut` にも同じ除外を適用していたが、現在は生成しない。

既存の `hg_lg_tail` が保存した `tail_rich` を再利用する。
対象は **低側 z < −5 が異なる物理 channel に3個以上ある event**。
`tail_study.root` の `selected` は使わず、resolution 独自の CoG 選択を維持する。
既定の判定元は `../result/channel_response/hg_lg_tail/threshold/e_ssa/<energy>GeV/`。
`inputs/data.tsv` / `inputs/mc.tsv` の calibration ファイルの実パスで sample と file_index を決め、
entry ごとに Run_Num・Event_Time・Event_Num（tail 側 TriggerID）を照合する。
判定元にない入力、全 event を覆わないサンプル解析、ID 不一致は出力作成前に停止する。
別 calibration / MC を使う場合は、対応する全 event の tail 解析結果を指定する。

過去の設定 `hl_tail_samples.sh` は SingleEnergy でdataのPS/SPSとMC `e_ssa/threshold`を解析する。
旧サーバーの入力パスとSlurmを使うため、新規のKEKCC実行には上記momentum設定を使う。
MultiEnergy は既定で data の8系列を比較する。`HL_TAIL_MULTI_TARGET=mc` を指定すると、MC の8系列を比較する。
CoG は200 / 20 / 10 / 5 mm。設定ファイルで変更できる。

```bash
cd analysis/resolution
make SingleEnergyAnalysis MultiEnergyAnalysis

# 既定: tail-rich を除く。dry-run は投入しない。
bash run_samples.sh --dry-run single hl_tail_samples.sh
bash run_samples.sh single hl_tail_samples.sh
# 8系列比較には tail-rich を残す SingleEnergy 結果も必要（下記）。

# tail-rich を残す SingleEnergy も作成する場合:
bash run_samples.sh --keep-tail-events single hl_tail_samples.sh
# 両方の SingleEnergy job の完了後、data の8系列をまとめて比較:
bash run_samples.sh multi hl_tail_samples.sh
# MC の8系列だけを比較（出力名は mc_all_cog_tail）:
HL_TAIL_MULTI_TARGET=mc bash run_samples.sh multi hl_tail_samples.sh

# 別の判定元を指定する例。DIR の下に <energy>GeV/ が必要。
bash run_samples.sh --tail-results ../result/channel_response/hg_lg_tail/threshold/e_ssa_noHLintercept \
  single matching_samples.sh
```

`--exclude-tail-events` も明示できる。shell のオプションは `single|multi` より前に置く。
`--tail-results` は SingleEnergy に渡す判定元であり、calibration 入力は設定ファイルで選ぶ。
複数の判定元を1設定で使う場合は、それぞれの `single_sample` の前で `tail_results_dir` を設定する。
過去のHL interceptは `hl_intercept_samples.sh`、tail有無の比較は `hl_tail_samples.sh` に記録されている。

出力は次のように分け、除去あり・なしで上書きしない。

```text
analysis/result/resolution/unknown/custom/hl_tail/{data,mc}/
  exclude/<CoG>mm/<energy>GeV_<CoG>mm.root
  keep/<CoG>mm/<energy>GeV_<CoG>mm.root
analysis/result/resolution/unknown/figures/custom/hl_tail/{data,mc}/
  exclude/<CoG>mm/<energy>GeV/...
  keep/<CoG>mm/<energy>GeV/...
```

C++ を直接呼ぶ場合は、既存の位置引数の**後ろ**に
`--exclude-tail-events` / `--keep-tail-events` を追加する。
SingleEnergy のみ `--tail-results DIR` を指定できる。省略時の判定元は実行ファイルからの相対位置。
MultiEnergy は event tree を持たないため、その場で event 除外は行わない。
SingleEnergy の除外条件メタデータを検査し、指定と異なる ROOT や未完了の ROOT は拒否する。
メタデータのない旧 ROOT は `--keep-tail-events` のときだけ利用可能。
除去あり・なしを1枚に重ねる場合は、以下の明示的な比較モードを使用する。

SingleEnergy ROOT の `hl_tail_excluded` に除外有無、`hl_tail_checked_events` と
`hl_tail_rejected_events` に CoG cut 前の照合数・除外数、`hl_tail_events_after_cog` に
tail、指定時の前方位置cut、CoG cut 後（`beforeEventCut`）の件数を保存する。
除去なしでは照合を省略するため checked/rejected は0。
判定元・定義・完了状態は `hl_tail_source` / `hl_tail_definition` / `hl_tail_status` に保存する。

```bash
# ROOT と一致する Python を使う。この環境では python3.12。
python3.12 tests/test_hl_tail.py
python3.12 tests/test_sample_comparison.py
```

合成 Calib_Hit / tail tree を用い、既定の除去と除去なし、境界の2/3 channel、
重複 TriggerID、tail 側 selected に依存しない除去、全 event coverage と ID の照合、
MultiEnergy の条件不一致拒否、SingleEnergy の σ/μ の引継ぎを検証する。


### data の CoG 4条件 × tail cut 有無の8系列を重ねる

`hl_tail_samples.sh` の MultiEnergy 設定は、data の200 / 20 / 10 / 5 mmそれぞれについて
`data/exclude` と `data/keep` を登録し、8系列を1枚に重ねる。凡例は `Data tail exclude/keep` とCoG。
hl_tail設定のROOTは `unknown/custom/hl_tail/{data,mc}/<条件>/<CoG>mm/`、
図は対応する `unknown/figures/custom/hl_tail/` 以下を参照する。配置や入力の由来を推定して新しいsampleへ統合しない。
200 mmはCoG cutなし相当の条件で、図では `nocut` と表示する。

```bash
# 未作成のSingleEnergy結果がある場合に実行し、全ジョブの正常終了を待つ。
bash run_samples.sh single hl_tail_samples.sh
bash run_samples.sh --keep-tail-events single hl_tail_samples.sh

# MultiEnergyは両方の結果を読むため、keep/excludeオプションは付けない。
bash run_samples.sh --dry-run multi hl_tail_samples.sh
bash run_samples.sh multi hl_tail_samples.sh
```

図は `../result/resolution/unknown/figures/custom/hl_tail/comparison/data_all_cog_tail/` の
`resolution.png` / `resolution_with_noise.png` と、それぞれのlog–log版
`resolution_loglog.png` / `resolution_with_noise_loglog.png`。
ROOTは `../result/resolution/unknown/custom/hl_tail/comparison/data_all_cog_tail.root`。

設定の `multi_options+=(--compare-tail-selections)` により条件の混在を明示的に許可する。
各系列の全エネルギーでは同じtail条件であることと、SingleEnergyが正常終了したことを検査する。
このモードは除去有無のメタデータを必須とし、旧ROOTから条件を推測しない。
C++を直接呼ぶ場合は `--labels ... --compare-tail-selections` を指定する。
比較ROOTの `hl_tail_excluded=-1` は混在比較を表し、各系列の条件は
`hl_tail_excluded_sample_0` 以降に0（keep）または1（exclude）で保存する。

8系列比較の検証: `python3.12 tests/test_tail_comparison.py`。

### compare の一時除外: data 200 GeV

MultiEnergyAnalysis は、入力パスに `data/` を含む200 GeVを分解能の値（σ/μ・誤差）と分解能fitから除外します。MCの200 GeVは残します。ファイル一覧からは削除せず、`beforeEventCut` の分布を `summary_200GeV.png`、layer/channel比較と保存canvasに表示します。このdataの描画ではエネルギーfit値をrebinや軸範囲にも使いません。

この対象だけは `hl_tail_status=incomplete` でもsummary用として読めますが、tail選別条件の一致と描画用ヒストグラムの存在は検査します。他の入力の完了検査は維持します。除外方針は出力ROOTの `resolution_value_exclusions` に記録します。任意名のパスや表示ラベルからdataを推測しません。

### 不正なfit値を含む比較

欠落したenergy fit、非正・非有限の平均/σ、非有限・負のfit誤差は、summary分布の表示を続けながら分解能の値取得から除外します。除外した系列・energy・理由・入力パスを標準エラーと比較ROOTの `resolution_skipped_points`（TSV形式のTNamed）に記録します。hit数fitが不正な場合はfit由来の軸範囲を使用しません。data 200 GeVの一時除外も維持します。

after-cut channel 比較は現在停止しており、未完了のdata 200 GeVもbefore-cutのsummary・layer/channel分布のみを表示します。他の入力の完了状態・tail条件・必須分布の検査は維持します。
