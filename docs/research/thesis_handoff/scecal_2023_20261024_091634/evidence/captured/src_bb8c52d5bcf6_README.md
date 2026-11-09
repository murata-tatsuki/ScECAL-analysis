# SiPM gain の採用値と calibration 互換 ROOT

既存の `sipm_gain` の比較結果から、条件ごとに再現する LED gain を選ぶ C++/ROOT コードです。bootstrap、元ADCの再読込み、再fit、温度補正は行いません。既存 calibration の置換も行いません。

入力は `analysis/result/sipm_gain/<label>/<comparison-name>/layers/layer_XX/numeric.root` の `Gain` tree です。元の compare の全30 layer の数値処理と merge が完了している必要があります。compare の plots 完了は不要です。

既存 `sipm_gain` は配下のソース全体を照合するため、完成済み比較結果をそのまま使えるよう、この後処理を隣の独立ディレクトリに置いています。

## 採用方法

PS/SPS、LED測定日時（manifest の `period`）、温度bin、hit/gain tag、memory cell、宣言された動作設定が同じものを一条件とします。DACを跨いで比較しますが、異なる日時・温度binを混ぜません。対象は memory cell 0、HitTag 1、GainTag 1 です。

既定の判定は次のとおりです。

1. LED の multi-Gaussian と FFT が、元解析でともに `candidate` になっていること。
2. その2手法の gain の対称相対差 `2*abs(g1-g2)/(g1+g2)` が5%以下であること。ほかの手法も `candidate` なら同じ一致条件を要求します。`review` の手法は採用の根拠にも拒否の根拠にも使いません。
3. 有効な温度、fit誤差、fit品質、安定性の記録があること。
4. 同一条件で3種類以上のDACが上記を満たすこと。同一DACの繰り返しをまず中央値にし、そのDAC別中央値の中央値を採用gainにします。イベント数やfit誤差による重み付けはしません。
5. 使用する個々のfitのgainが採用gainから5%以内にあり、温度の平均値の最大–最小が0.5℃以下であること。

5%・3DACは保守的な品質判定の既定値で、統計的信頼水準ではありません。複数手法は同じイベントに依存するため、独立測定として誤差を縮めません。必要なら別 `--output-label` で基準への依存性を確認できます。

electron/muon は、同じbeam・tag・memory cell、温度差0.5℃以内で、primary fit と FFT がともに candidate かつ5%以内に一致するものを `BeamChecks` に保存します。異なるrunの動作設定が未確認のため、LEDとの一致・不一致を補助的に表示し、採用gainの平均や採否判定には使いません。muonのMIPピーク位置そのものを1 p.e. gainとして扱うこともありません。

既定の primary は `multi_gaussian` です。`--primary-method envelope_comb` に変更する場合も、FFTとの一致と3DAC以上の再現性を要求します。LEDで決められないchannelをelectron/muon・隣接channel・旧gainで埋めません。

既知のdead channelを明示的に除外する場合は、物理CellID（`layer*100000+chip*10000+channel`）を1行1個にした既存のリストを `--dead-channels PATH` で指定してください。指定なしでも、根拠の不足したchannelは除外します。ただし「gain未決定」と「dead」の断定は区別します。

## 出力と branch

基点は `analysis/result/sipm_gain_selection/<output-label>/` です。`--output-label` を省略すると入力と同じlabelになります。

| 保存先 | 内容 |
|---|---|
| `inputs/configuration.json`, `inputs/source/` | 入力チェックサム、基準ROOTの形式、採用条件、ソースの記録 |
| `layers/layer_XX/selection.root` | `Decisions`、`Evidence`、`BeamChecks`。対応するTSVも保存 |
| `selection.root` | 全layerの診断treeを統合。異なる条件を混ぜた `T_Event` は作らない |
| `conditions.tsv` | 条件ごとの採用数、未採用数、calibration互換ROOTへの相対パス |
| `data/<beam>/<period>/temperature_<lower>-<upper>C/hit_1_gain_1_memory_0/settings_<hash>/gain.root` | **その条件だけの calibration 互換ROOT** |
| `plots/data/<同じ条件階層>/layer_XX/` | 条件別の `overview.png` と `channel_<physicalID>.png` |
| `plots/layers/layer_XX/plots.root`, `plot_index.tsv` | 保存canvasとPNG一覧 |

基準ROOTは既存解析の次のファイルです。

```text
/megraid01/users/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share/20230504_1425_Calib_LED_all.root
```

条件別 `gain.root` の `T_Event` には以下の7branchを同じ名前・型で必ず作ります。採用channelが0なら空のtreeになります。

| branch | ROOT型 | 値の意味 |
|---|---|---|
| `cellIDs` | `Int_t` | **旧calibration表記** `layer*1000000+chip*10000+channel` |
| `gain` | `Double_t` | 採用した1 p.e. gain [ADC/p.e.] |
| `gainError` | `Double_t` | 下記の保守的な品質幅。統計的な標準誤差ではない |
| `chiSquare` | `Double_t` | 採用値に最も近い代表LED fitの**Poisson deviance**。中央値自体のχ²ではない |
| `NDF` | `Double_t` | 同じ代表fitの自由度 |
| `temp` | `Double_t` | 使用したLED spectrumの平均温度の中央値 [℃] |
| `tempError` | `Double_t` | 温度の最大偏差と入力温度RMSの最大値。平均温度の標準誤差ではない |

`gainError` は、使用したfitの誤差の中央値、個々のgainと採用値の最大偏差、元解析の `stability_ADC` の最大値、candidate手法間のgain差の最大値、のうち最大のものです。bootstrapなしのため、統計的な1σ・信頼区間とは解釈しないでください。

`gainErrorKind` と `qualityStatistic` を追加branchとして保存し、誤差とfit品質の定義を明示します。`metadata_json` にも定義を保存します。したがって名前・型の読込互換性は保ちますが、既存の `chiSquare/NDF` のχ²カットや `gainError` による統計的重み付けを、そのまま意味が同じとして適用しないでください。

同じファイルの `Selection` tree に採用理由の詳細、`Rejected` tree に未採用・データ不足・明示deadの理由を残します。各条件で全6300物理channelをこの2treeに分け、`T_Event` には採用分だけを保存します。診断treeの `cellid` は物理CellIDであり、`T_Event.cellIDs` と桁が違うことに注意してください。既存calibrationの読込側が欠損を平均gainで補完する場合、その挙動は別途読込側で決まります。このコードは補完値を出力しません。

入力にGainTag 0など対象外の条件が含まれる場合も、その条件の除外記録を保存します。そのROOTの `T_Event` は空です。使用可能な条件は `conditions.tsv` の `selected_channels` と `status` で確認してください。

## 実行順序

ROOTとC++コンパイラが使える環境で実行してください。以下の各段階は、前のジョブが正常終了してから**都度**実行します。全段階を自動投入する処理はありません。

### 0. compile

```bash
cd /home/murata_t/ScECAL_BeamTest
make -C analysis/sipm_gain_selection -j2
mkdir -p analysis/log/sipm_gain_selection/full_grouped/{select,merge,plots}
```

### 1. layer別の採用判定（最大32並列）

```bash
sbatch --job-name=sipm_select --partition=all \
  --array=0-29%32 --cpus-per-task=1 --mem=8G --time=08:00:00 \
  --chdir=/home/murata_t/ScECAL_BeamTest \
  --output=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/select/%A_%a.out \
  --error=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/select/%A_%a.err \
  --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0; ./analysis/sipm_gain_selection/Select --label full_grouped --stage select --layer "$SLURM_ARRAY_TASK_ID"'
```

### 2. 条件別 calibration ROOT の生成（1ジョブ）

```bash
sbatch --job-name=sipm_select_merge --partition=all \
  --cpus-per-task=1 --mem=8G --time=08:00:00 \
  --chdir=/home/murata_t/ScECAL_BeamTest \
  --output=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/merge/%j.out \
  --error=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/merge/%j.err \
  --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0; ./analysis/sipm_gain_selection/Select --label full_grouped --stage merge'
```

全30layerの `SELECTION_COMPLETE` がなければ停止します。未完了分を省いて完成ROOTを作ることはありません。

### 3. plot（最大32並列）

```bash
sbatch --job-name=sipm_select_plots --partition=all \
  --array=0-29%32 --cpus-per-task=1 --mem=8G --time=08:00:00 \
  --chdir=/home/murata_t/ScECAL_BeamTest \
  --output=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/plots/%A_%a.out \
  --error=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain_selection/full_grouped/plots/%A_%a.err \
  --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0; ./analysis/sipm_gain_selection/Select --label full_grouped --stage plots --layer "$SLURM_ARRAY_TASK_ID"'
```

元の `Analyze compare --stage plots` の代わりに、この新しいコマンドで採用値の図を作ります。元の比較plotを残したまま併用できます。

既定では採用channelの個別図と全条件のoverviewを作り、未採用channelの大量の個別PNGは省きます。未採用の根拠もROOT/TSVには全件保存します。未採用channelの個別図も必要なら `--plot-scope all` を指定します。

少数の確認では、異なる出力labelで `--layer 0` だけをselectし、plotsに `--plot-limit 4` を加えられます。plotsはそのlayerの判定結果があれば動作しますが、mergeは全30layerが必要です。`--plot-limit` はlayerあたりのchannel図数で、overviewは全条件を出します。既定0は選んだscopeの全channel図です。

完了済みの段階はチェックサムを照合してスキップします。採用判定・統合のソース、採用条件、入力を変更した再解析は別の `--output-label` を指定し、全段階で同じlabelを使ってください。途中終了した未完了段階は同じコマンドで再実行できます。

表示だけを変更した場合は、compile後に**段階3だけ**を再実行します。保存済みの判定結果と元ソースの記録を検証して読み込み、描画ソースが変わっていれば同じ保存先のPNGとROOT canvasを再作成します。gainの再判定・calibration ROOTの再生成は行いません。描画用のソースと温度の換算根拠は、各layerのplot出力内の `inputs/source/` と `plot_configuration.json` に記録します。

## 図の読み方

- `overview.png`：上段は採用gainと品質幅、中段は品質幅/gain、下段は0=LEDなし・1=未決定・2=採用・3=明示deadです。横軸は `chip*36+channel`。
- 図上部の温度は `T: 23.5-24.0 ℃` のように表示します。ADC分布作成時の保存設定から区間幅を読み、下端以上・上端未満の温度範囲に換算しています。温度情報のない条件は `Temperature unavailable` と表示します。保存先も `temperature_23.5-24.0C` の形式です。診断tree内のbin番号は条件識別用として保持します。
- `channel_*.png`：横軸DAC、縦軸gain。色は4手法で固定し、塗りつぶしがcandidate、薄い白抜きがreviewです。黒い大きな丸は採用値の算出に使ったprimary fit。青線と帯は採用gainと品質幅です。
- 未採用の図には理由を表示します。`too_few_qualified_dacs` は3DACに届かない、`dac_or_repeat_inconsistent` はDAC間・繰返し間の不一致、個々のspectrumの `method_conflict_*` は手法間の不一致です。詳細は `decisions.tsv`・`evidence.tsv`。
- `beam=supports/disagrees/mixed/not_available` は温度の近いelectron/muonの補助比較です。設定未確認なので、この表示だけで採用gainを変更しません。

採用数は条件ごとに数えてください。異なる条件の採用数の合計は、重複のないchannel数ではありません。

## 実データでの小規模確認

`full_grouped` のlayer 0と4で採用判定とPNG/ROOT描画を確認しました。保存先は `analysis/result/sipm_gain_selection/full_grouped_selection_validation/` です。このlabelは検証用で、全layerの統合はしていません。

| layer | 処理時間（select） | 少なくとも1条件で採用できたchannel | 採用channel×条件の件数 |
|---|---:|---:|---:|
| 0 | 32.50秒 | 75 / 210 | 158 |
| 4 | 34.09秒 | 3 / 210 | 3 |

これはmegofl04での測定値で、Slurm待ち時間・同時I/O負荷は含みません。採用数が少ないlayerもありますが、未決定の値を補完していません。実行コマンド、条件、入力チェックサムは各layerのsummaryとinputsに保存しています。

## 温度のdirectory名への移行

旧出力の `temperature_bin_47` を `temperature_23.5-24.0C` に変更し、conditions.tsv・plot_index.tsv・完了印の参照先を更新します。温度幅は元のADC収集設定を読みます。ROOT・PNGのバイト内容は変更せず、移行前後のchecksumを検証します。途中で停止した場合も同じコマンドで再開できます。

```bash
./analysis/sipm_gain_selection/Select --label full_grouped --stage temperature-paths
```

`conditions.tsv` には `temperature_lower_C`、`temperature_upper_C` も保存します。温度不明は `temperature_unavailable` です。処理記録は `inputs/temperature_path_migration.json` に保存します。

## 必要DAC数を2にした比較

既定値は3DACのままです。`--min-dacs 2` で、手法間の一致5%・DAC間の変動5%を維持した比較を行えます。異なる条件なので出力labelを分けます。下のselectを全30layerで完了してから、merge、plotsをそれぞれ実行してください。既存のADC作成・fit・compareは再実行しません。

```bash
./analysis/sipm_gain_selection/Select --label full_grouped --output-label full_grouped_relaxed_2dac --stage select --min-dacs 2 --layer 0
# 全30layerのselect終了後
./analysis/sipm_gain_selection/Select --label full_grouped --output-label full_grouped_relaxed_2dac --stage merge
./analysis/sipm_gain_selection/Select --label full_grouped --output-label full_grouped_relaxed_2dac --stage plots --layer 0
```

Slurmでは従来の30layer配列コマンドに同じ `--output-label full_grouped_relaxed_2dac` を付け、selectに `--min-dacs 2` を追加します。ログもこのlabelの下に分けます。2DACでは再現性確認の点数が減るため、3DAC基準の結果と区別して扱います。electron/muonの温度別run統合は [sipm_gain_pool](../sipm_gain_pool/README.md) にあります。
