# CoG / channel scan

layer 9・10 の重心でイベントを選び、そのイベントの指定範囲内 channel のみを集計する。
どちらの範囲も x・y 方向の半幅 [mm]。境界は含まない。

- `cog_range > 0`: layer 9・10 の両方で `abs(CoG_X) < cog_range` かつ
  `abs(CoG_Y) < cog_range` を要求する。重心を定義できない層（エネルギー和が0）が
  あるイベントは除外する。`cog_range=0` は CoG カットなし。
- CoG は **channel カット前**の全有効 hit のエネルギー重み付き重心。
  従来どおり負の hit energy と、設定に応じた bad channel を除く。
- 集計する hit は `abs(Hit_X) < channel_range` かつ
  `abs(Hit_Y) < channel_range` を満たす layer 0〜29 の channel。
  `channel_range` は正の値を指定する。200 mm は検出器全体を含む。

CoG の図も channel カット前の重心を表示する。
`beforeEventCut` / `afterEventCut` は **どちらも CoG と channel カット適用後**で、
従来の energy–Nhit イベントカットの前後を表す。
energy–Nhit カットは範囲内 channel の energy と Nhit で計算する。
この解析では従来どおり `beforeEventCut` を比較する。
同じ Nhit での energy は各 ROOT の `beforeEventCut/hit_vs_e_<energy>GeV`
（x: energy、y: Nhit）から確認できる。

## 入力と tail 除外

MCは `Result_MC/calib/e-/e_ssa/threshold/<energy>GeV/*.root` を使用する。
旧 `e-/sps/threshold` は使用しない。e_ssaのファイル名には `e-` prefixがないため、`*.root` を読む。
dataは従来と同じ2023 PS/SPSのcalibration入力を使用し、**HG–LG低側tail-richイベントを除外する**。
低側 `z < -5` が異なる物理channelに3個以上あるイベントを丸ごと除外し、MCはtailイベントを保持する。

判定は既存 `../result/channel_response/hg_lg_tail/threshold/e_ssa/<energy>GeV/` の
`tail_study.root` と `inputs/{data,mc}.tsv` を、resolution共通の `HLTailSelection.hh` で読む。
`selected` は使わず、本解析独自のCoG/event cutを維持する。
calibration実パス・全entry・Run_Num/Event_Time/Event_Numを照合し、未登録入力、
判定の欠落、ID不一致、tail判定より新しいcalibrationは出力作成前にエラーにする。
`TAIL_RESULTS` で別の判定元を指定できるが、同じ入力を全event処理した結果が必要。

## 実行

ROOTとg++を使える環境で、このディレクトリから実行する。

```bash
make -j2
# コマンドと入力を確認するだけ。ファイル作成・ジョブ投入なし。
DRY_RUN=1 bash analysis_paralell_allCondition.sh
# dataとMCを再生成。表示される配列job IDを控える。
MAX_PARALLEL=24 bash analysis_paralell_allCondition.sh
```

19 energy（0.5, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 250 GeV）
× CoGなし/20/10/5 mm × channel半幅68.2/22.5 mm × data/MCの**304 task**をSlurm配列jobで投入する。
既定は最大24 task同時、各taskは1 CPU・4 GiB。配列jobを二重投入しないよう、同名jobが残っている場合は停止する。
`CoG=0` はカットなし。200 mmは重心を定義できないイベントを除外するため0とは異なる。

比較は、全SingleEnergy taskの正常終了後に実行する。Slurmの依存関係を使う例：

```bash
# ARRAY_JOB_IDを、上で表示された配列job IDに置き換える。
ARRAY_JOB_ID=1234567
sbatch --dependency="afterok:${ARRAY_JOB_ID}" \
  --output="$PWD/../log/cog_channel_scan/compare-%j.out" \
  --error="$PWD/../log/cog_channel_scan/compare-%j.err" \
  execute_paralell.sl bash "$PWD/compare_ds_allCondition.sh"
squeue -u "$USER"
```

完了済みのSingleEnergy結果からローカルで比較する場合は `bash compare_ds_allCondition.sh`。
`dataThre_cog_channel`（8条件×data/MC）と `dataThre`（CoGなし、hit範囲2条件×data/MC）の両方を生成する。
比較前に、全入力が完了済み・dataはtail除外済み・MCはe_ssa/thresholdでtail保持、
かつROOTのCoG/channelメタデータがファイル名と一致することを検証する。旧SPS結果へのfallbackはしない。

対象を絞る場合は、**生成・比較で同じ環境変数**を指定する。

```bash
# 例: 100 GeVだけ、CoGなしと5 mm、hit範囲2条件（8 task）
ENERGIES=100 COG_RANGES="0 5" DRY_RUN=1 bash analysis_paralell_allCondition.sh
# 同じ条件で比較する場合:
ENERGIES=100 COG_RANGES="0 5" bash compare_ds_allCondition.sh
```

`CHANNEL_RANGES`、`RESULT_BASE`も両script共通。生成側の `SAMPLES=data` / `SAMPLES=mc` で片方のみ再生成可能。
比較側には両sampleの結果が必要。CoGの指定に0がない場合は `dataThre` を生成しない。
小数の表記はファイル名に使うため、両scriptで一致させる（5と5.0は別パス）。

## 出力

基点は `../result/cog_channel_scan/`。退避済みの `cog_channel_scan.bak/` は変更しない。
ROOTファイル名は `condition=cog<cog_range>mm_channel<channel_range>mm` を維持する。

- data ROOT: `data/exclude/channel<channel_range>mm/<cog_range>mm/<energy>GeV_<condition>.root`
- MC ROOT: `mc/threshold/e_ssa/channel<channel_range>mm/<cog_range>mm/<energy>GeV_<condition>.root`
- 図: ROOTと同じ階層を `figures/` の下に置く。energy別PNGはその下。
- 比較ROOT: `comparison_dataThre_cog_channel.root` と `comparison_dataThre.root`
- 比較図: `figures/comparison/dataThre_cog_channel/` と `figures/comparison/dataThre/`
- 再現用記録: `run.log`、`inputs/source_checksums.txt`。各ROOTにも全calibration実パスとsample名を保存。
- Slurmログ: `../log/cog_channel_scan/single-<array ID>_<task ID>.out/.err`
- 投入記録・全引数のtask manifest: 同ログディレクトリの `submission_<run ID>.tsv`、`tasks_<run ID>.tsv`

各単一energy ROOTには `cog_range_mm`、`channel_range_mm`、`analysis_sample`、`calibration_inputs`、
`input_events`、`hl_tail_excluded`、`hl_tail_checked_events`、`hl_tail_rejected_events`、
`hl_tail_events_after_cog`、`hl_tail_source`、`hl_tail_definition`、`hl_tail_status` を保存する。
`complete` のstatusは全処理終了時だけ書く。tail除外件数はCoG適用前で、2回のeventループを重複計数しない。

生成scriptは `--skip-channel-plots` を渡し、360枚/energyのchannel別canvasを省略する。
energy・Nhit・layer別・CoG・SiPM別の集計と比較図は通常どおり生成する。
channel別canvasも必要なら、SingleEnergyを直接呼び出してこのオプションを省く。

比較図は CoG/channel の組合せごとに色・形を割り当て、同条件の data は
塗りつぶし＋実線、simulation は白抜き＋破線で表示する。
凡例の値は `5 / 22.5` のように **CoG / channel** の順とし、単位と順序は
凡例の見出しにまとめる。CoG の0は `no cut / 22.5` のように `no cut` と表示する。
number of hits の表示範囲は全条件の分布に合わせる。

resolution 図の凡例は右上パネルだけに表示する。fit の式は共通の見出しにまとめ、
各行は条件と係数 `(a, c)` または `(a, b, c)` を表示する。
条件が多い場合は2列にし、行数に応じて文字サイズを調整する。

- `resolution.png`: noise 項なし、linear–linear
- `resolution_loglog.png`: noise 項なし、log–log
- `resolution_with_noise.png`: noise 項あり、linear–linear
- `resolution_with_noise_loglog.png`: noise 項あり、log–log

log–log は同じデータ点・fit 結果を使い、再fit は行わない。
energy と resolution がともに正かつ有限な点を表示する。
y 軸の上限は2に固定し、下限は正の値・誤差から設定する。
fit の除外範囲（25 < E < 35、75 < E < 85 GeV）でも曲線は同じ係数で連続して描画する。
非正値の点は log 図だけで非表示とし、元データ・linear 図は維持する。
linear–linear は従来の表示範囲を維持する。
比較 ROOT にも同名（拡張子なし）の4つの canvas を保存する。

さらに、同じ CoG / channel 条件の data と simulation だけを重ねた図を追加する。
条件ごとに noise 項なし・ありを1枚の PNG にまとめる。data は青、simulation は橙。

| | 左: linear–linear | 右: log–log |
| --- | --- | --- |
| 上段 | noise 項なし | noise 項なし |
| 下段 | noise 項あり | noise 項あり |

各パネルの右上に data / simulation の2本の fit 関数を表示する。
係数名 `(a, b, c)` の一覧ではなく、`σ(E)/E = 0.221/√E ⊕ 0.279/E ⊕ 0.010`
のように数値を代入した式とする（`⊕` は二乗和の平方根）。
小さな係数は指数表記で表示する。軸範囲は全条件図と揃え
（linear は x: 0〜130 GeV、y: 0〜0.3、log の y 上限は2）、点・fit 係数も同じ値を使う。
全条件を重ねた既存4枚も引き続き出力する。

ファイル名は `resolution_pair_<condition>_<threshold>.png`。
例: `resolution_pair_cog0mm_channel22.5mm_threshold.png`。
現在の8条件では合計8枚を、既存図と同じ比較図ディレクトリに保存する。
simulation の threshold が複数ある場合は、それぞれ data と1対1で比較する。
比較 ROOT にも PNG と同名（拡張子なし）の4パネル canvas を追加する。
`bash compare_ds_allCondition.sh` で既存図と一緒に生成する。


## 単独実行・従来形式

```text
SingleEnergyAnalysis output.root Nenergy E1...En E1Nfiles...EnNfiles input1.root ... cog_range_mm channel_range_mm b_only_best figure_path [options]
```

末尾の範囲引数は CoG、channel の順。
オプションは `--exclude-tail-events` / `--keep-tail-events`、`--tail-results DIR`、
`--sample-label data|e_ssa`、`--skip-channel-plots`。直接実行時は従来互換のためtail保持が既定であり、
生成scriptはdataの除外とMCの保持を必ず明示する。
比較プログラムは末尾の `--require-e-ssa-data-tail` で本番条件と来歴を検証する（比較scriptが自動指定）。
以前の channel のみの呼出し（`cog_range_mm` 省略）も CoG カットなしとして利用可能。
`MultiEnergyAnalysis` は旧 `100GeV_22.5mm.root` 形式も読み込める。

Gaussian fit の方法・下限15 MeVは従来どおり。
狭い範囲・低エネルギー・厳しい CoG 条件では、統計と fit の成立を確認すること。

## 検証

```bash
python3.12 test_channel_cut.py
```

人工の calibration ROOT を実際の `SingleEnergyAnalysis` で処理する。
CoG なしの従来形式と新形式、CoG 2.5/5 mm × channel 22.5 mm を検証する。
channel 範囲外の hit が CoG に効くこと、エネルギー重み、矩形条件、両層の条件、
境界、重心を定義できない層、bad channel、負のエネルギーを含める。
イベントカット前後の総エネルギー・Nhit・30層の分布・SiPM別分布を確認する。
テストの ROOT・図・ログは `../result/cog_channel_scan/test_channel_cut/` に保存する。
