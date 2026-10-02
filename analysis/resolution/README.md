# Resolution の sample 比較

Slurmの標準出力・標準エラーは `analysis/log/resolution/` に保存します。`singleEnergy/`、`compare/`、`hl_tail/<条件>/<工程>/` などの区分はその配下で維持します。実行中・待機中のジョブが使う旧 `jobs/` のログは、ジョブ終了後に移動します。

## 新規出力の共通ディレクトリ規則

今後の条件比較は、**対象（`data` / `mc`）→ 比較条件 → CoG** の順で保存する。
`exclude/keep`、`before/after`、その他の比較条件でも同じ規則を使い、
dataのみ・MCのみの解析でも対象の階層は省略しない。

```text
<基点>/data/exclude/<CoG>mm/
<基点>/data/keep/<CoG>mm/
<基点>/mc/exclude/<CoG>mm/
<基点>/mc/keep/<CoG>mm/

<基点>/data/before/<CoG>mm/
<基点>/data/after/<CoG>mm/
<基点>/mc/before/<CoG>mm/
<基点>/mc/after/<CoG>mm/
```

resolutionの基点は、ROOTが `analysis/result/resolution/custom/<テーマ>/`、
図が `analysis/result/resolution/figures/custom/<テーマ>/`。
ROOTはCoGディレクトリ内の `<energy>GeV_<CoG>mm.root`、
図はCoGディレクトリ内の `<energy>GeV/` 以下に保存する。
複数条件をまとめる比較ROOT・図は、従来どおり各基点の `comparison/` 以下に保存する。

新しい設定の保存先・読込先はこの階層にそろえる。
`exclude/data/` や `keep/data/` 等の旧配置は、既存結果の読込互換のために参照する場合がある。
新規出力には使用せず、既存結果の移動は明示的に依頼された場合だけ行う。

## 任意の calibration / 条件を shell で選ぶ

`calibration ROOT → SingleEnergyAnalysis → 条件別の解析 ROOT → MultiEnergyAnalysis`
の流れは従来と同じ。SingleEnergy の出力は cut 済み event tree ではなく、分布・fit・canvas を保存した解析 ROOT。
SingleEnergy の CoG cut は layer 9 と 10 の両方で `abs(x) < range && abs(y) < range` を要求する。
threshold を変更する引数はなく、その条件は入力 calibration sample 側で選ぶ。
bad-channel 除外と既存 event cut、Gaussian fit の処理は維持する。
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
HL intercept用の既存設定は `hl_intercept_samples.sh`、tail比較は `hl_tail_samples.sh`。相対パスの基準は常に `analysis/resolution/`。
既存の `analysis_paralell*.sh` / `compare_ds*.sh` も従来どおり使用できる。

### SingleEnergy の選択

```bash
single_sample ROOT保存先 図の保存先 CoG範囲mm bad_channel除外 FILE_GLOB calibration_dir...

# 入力は calibration_dir/10GeV/*.root など。最後に複数の入力ディレクトリも指定可能。
single_sample '../result/resolution/custom/my_correction/data/A/20mm' \
  '../result/resolution/figures/custom/my_correction/data/A/20mm' \
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
multi_sample 'Calibration A, CoG 20 mm' '../result/resolution/custom/my_correction/data/A/20mm/'*.root
multi_sample 'Calibration B, CoG 20 mm' '../result/resolution/custom/my_correction/data/B/20mm/'*.root
multi_sample 'MC alternative, CoG 10 mm' '../result/resolution/custom/my_correction/mc/C/10mm/'*.root

multi_output='../result/resolution/custom/comparison.root'
multi_figures='../result/resolution/figures/custom/comparison'
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
- `summary_resolution.png` とROOTの同名canvasにも4条件を保存する。`Divide(2,2)` で、
  上段linear・下段log–log、左列without noise・右列with noise。各パネルのタイトルは左列を `without noise term`、右列を `with noise term` とし、軸形式はタイトルに付けない。凡例は上段だけに表示する。
  各panelは対応する単独図の点・fit係数をそのまま使用する。
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

`hl_intercept_samples.sh` / `samples.example.sh` は `(200 20 10 5)` × after/beforeの8系列を1つの比較に登録する。
SingleEnergyは `<対象>/<条件>/<CoG>mm/<energy>GeV_<CoG>mm.root`、図も同じ対象/条件/CoG順で保存する。
HL interceptでは `mc/after/<CoG>mm/` と `mc/before/<CoG>mm/` を使う。
`samples.example.sh` の `sample_type` は入力に合わせて `data` または `mc` に設定する。
HL interceptの比較は新配置を優先し、その条件/CoGの新ディレクトリが存在しない場合だけ
旧 `after/<CoG>mm/`・`before/<CoG>mm/` を読む。
既存のHL intercept結果（before/after/comparison）はROOT・図ともに `mc/` 配下へ移動済み。
MCのみの比較なのでdata用ディレクトリは作成しない。
比較ROOTは `custom/hl_intercept/mc/comparison/all.root`、
図は `figures/custom/hl_intercept/mc/comparison/all/`（ともに `analysis/result/resolution/` 基準）。

```bash
for cog in "${cog_ranges[@]}"; do
  multi_sample "MC after" "$result_base/mc/after/${cog}mm/"*.root
  multi_sample "MC before" "$result_base/mc/before/${cog}mm/"*.root
done
multi_run "$result_base/mc/comparison/all.root" "$figure_base/mc/comparison/all"
```

`multi_run` は全CoGの登録後に1回呼ぶ。第1引数が比較ROOT、第2引数が図の保存先。
CoG別の個別比較も必要なら、各CoGの2系列を登録するループ内で `multi_run` を呼ぶ。
実行後は系列リストを空にする。出力の重複はエラー。`multi_output` / `multi_figures` による暗黙の1比較も利用可能。
現在の `hl_intercept_samples.sh` はSingle/MultiともSlurmを使い、multiは8系列をまとめて1ジョブ投入する。

```bash
# 既存のHL intercept比較はtail除外前のROOTなのでkeepを明示する。
bash run_samples.sh --keep-tail-events --dry-run multi hl_intercept_samples.sh
bash run_samples.sh --keep-tail-events multi hl_intercept_samples.sh
```

### Gaussian beam / SSA hitmap の既存MC比較

`sim_sps_e_ssa_samples.sh` はMultiEnergy専用で、保存済みの
`default/simulation/threshold/<CoG>mm/`（凡例 `oval`）と
`custom/hl_intercept/mc/before/<CoG>mm/`（凡例 `SSA hitmap`）を比較する。
CoGは200/20/10/5 mmの8系列。200 mmは凡例で `nocut` と表示する。
両方ともtail除外導入前の結果なので、`--keep-tail-events` が必須。
この指定は保存済みの選択条件を検査するもので、SingleEnergyのevent選択を変更しない。

```bash
bash run_samples.sh --keep-tail-events --dry-run multi sim_sps_e_ssa_samples.sh
bash run_samples.sh --keep-tail-events multi sim_sps_e_ssa_samples.sh
```

ローカル実行。ROOTと入力記録は
`analysis/result/resolution/default/comparison/sim_sps_e-_ssa/{all.root,run.log}`、
図は `analysis/result/resolution/figures/default/comparison/sim_sps_e-_ssa/` に保存する。
`oval` は実線・塗りつぶし、`SSA hitmap` は点線・白抜き。CoGごとの色・marker形状は共通。

### Data / SSA hitmap MC before の既存結果比較

`data_sim_e_ssa_samples.sh` は `default/data/<CoG>mm/`（凡例 `data`）と
`custom/hl_intercept/mc/before/<CoG>mm/`（凡例 `sim`）の生成済みSingleEnergy結果を比較する。
4 CoG × 2系列をローカル実行し、SingleEnergyは再計算しない。
0.5 GeVはdataの20 mmの保存fitが無効なため、全8系列から除外する（共通の18 energy、1〜250 GeV）。

```bash
bash run_samples.sh --keep-tail-events multi data_sim_e_ssa_samples.sh
```

ROOTと入力記録は `analysis/result/resolution/default/comparison/data_sim_e-_ssa/`、
図は `analysis/result/resolution/figures/default/comparison/data_sim_e-_ssa/` に保存する。
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

SingleEnergy の ROOT に保存された `<energy>GeV/(cut_)edep_Layer<L>_Chip<C>`
canvas から `edep_channel_<L>_<C>_<channel>` を読み、**同じ CoG 条件の data と simulation の2系列だけ**を重ねる。
入力順に依存せず条件を照合する。simulation の threshold が複数ある場合は、それぞれ別の図にする。
calibration の再読込や SingleEnergy の再実行は不要。
小数エネルギーでは、旧 SingleEnergy が使う整数のディレクトリ名（例: `0GeV`）も読む。

- 1枚の canvas は Single の raw と同じ6×6配置（channel 0〜35）。
- layer 0〜29 × chip 0〜5、event cut 前後の両方を保存する（1 energy・1条件あたり360枚）。
  CoG が4条件なら、各 chip・event cut 段階につき4枚を別ディレクトリに保存する。
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

各ディレクトリの `cut_30GeV_layer0chip0.png` は event cut 後。
既存 shell の `figure_path` は `../result/resolution/figures/comparsion/200mm/` または
`../result/resolution/figures/comparsion/dataThre/`（既存の `comparsion` 表記を維持）。
比較 ROOT にも `cog<cog_range>mm/<threshold>/<energy>GeV/raw/(cut_)edep_Layer<L>_Chip<C>`
として canvas を保存する。0.5 GeV なども出力では正確な energy 名を使う。
従来の energy summary・resolution の全条件比較図はそのまま出力する。


## HG–LG tail event を除いた resolution

SingleEnergy は既定で tail-rich event を除外してから、従来の CoG cut、
エネルギー分布作成、Gaussian fit を行う。MultiEnergy の resolution は、その
`beforeEventCut/fit_gaus` の **σ/μ**。hg_lg_tail の診断用 RMS/mean とは異なる。
`beforeEventCut` / `afterEventCut` の両方に同じ tail 除外を適用する。

既存の `hg_lg_tail` が保存した `tail_rich` を再利用する。
対象は **低側 z < −5 が異なる物理 channel に3個以上ある event**。
`tail_study.root` の `selected` は使わず、resolution 独自の CoG/event cut を維持する。
既定の判定元は `../result/channel_response/hg_lg_tail/threshold/e_ssa/<energy>GeV/`。
`inputs/data.tsv` / `inputs/mc.tsv` の calibration ファイルの実パスで sample と file_index を決め、
entry ごとに Run_Num・Event_Time・Event_Num（tail 側 TriggerID）を照合する。
判定元にない入力、全 event を覆わないサンプル解析、ID 不一致は出力作成前に停止する。
別 calibration / MC を使う場合は、対応する全 event の tail 解析結果を指定する。

用意した `hl_tail_samples.sh` は SingleEnergy で既存 data の PS/SPS と MC `e_ssa/threshold` を解析する。
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
HL interceptは `hl_intercept_samples.sh`、tail有無の比較は `hl_tail_samples.sh` を使用する。

出力は次のように分け、除去あり・なしで上書きしない。

```text
analysis/result/resolution/custom/hl_tail/{data,mc}/
  exclude/<CoG>mm/<energy>GeV_<CoG>mm.root
  keep/<CoG>mm/<energy>GeV_<CoG>mm.root
analysis/result/resolution/figures/custom/hl_tail/{data,mc}/
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
tail/CoG cut 後かつ従来 event cut 前の件数を保存する。
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
hl_tailの既存結果はROOT・図とも `data/<条件>/<CoG>mm/`、`mc/<条件>/<CoG>mm/` に整理済み。旧配置の重複を除去し、MultiEnergyもこの配置だけを参照する。
200 mmはCoG cutなし相当の条件で、図では `nocut` と表示する。

```bash
# 未作成のSingleEnergy結果がある場合に実行し、全ジョブの正常終了を待つ。
bash run_samples.sh single hl_tail_samples.sh
bash run_samples.sh --keep-tail-events single hl_tail_samples.sh

# MultiEnergyは両方の結果を読むため、keep/excludeオプションは付けない。
bash run_samples.sh --dry-run multi hl_tail_samples.sh
bash run_samples.sh multi hl_tail_samples.sh
```

図は `../result/resolution/figures/custom/hl_tail/comparison/data_all_cog_tail/` の
`resolution.png` / `resolution_with_noise.png` と、それぞれのlog–log版
`resolution_loglog.png` / `resolution_with_noise_loglog.png`。
ROOTは `../result/resolution/custom/hl_tail/comparison/data_all_cog_tail.root`。

設定の `multi_options+=(--compare-tail-selections)` により条件の混在を明示的に許可する。
各系列の全エネルギーでは同じtail条件であることと、SingleEnergyが正常終了したことを検査する。
このモードは除去有無のメタデータを必須とし、旧ROOTから条件を推測しない。
C++を直接呼ぶ場合は `--labels ... --compare-tail-selections` を指定する。
比較ROOTの `hl_tail_excluded=-1` は混在比較を表し、各系列の条件は
`hl_tail_excluded_sample_0` 以降に0（keep）または1（exclude）で保存する。

8系列比較の検証: `python3.12 tests/test_tail_comparison.py`。
