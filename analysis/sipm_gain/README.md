# 2023 PS/SPS の1 p.e. gain：手法別・bootstrapなし

C++17/ROOTをコンパイルして使う。Pythonは実行に不要。raw HG ADCからgain候補を推定し、LED・electron・muonおよび手法間で比較する。calibration定数の更新・resolution再計算は行わない。

比較結果から条件別の採用gain・calibration互換ROOTを作る後処理は、[sipm_gain_selection](../sipm_gain_selection/README.md) にある。既存のfitやcompareを再実行せず、採用判定・ROOT統合・plotを別段階で実行する。

## コードの配置

```text
analysis/sipm_gain/
  collect/             raw読み込み、memory cell 0のADC分布作成・LED group統合
  multi_gaussian/      自由なピーク面積を持つ等間隔Gaussian群
  envelope_comb/       滑らかな全体分布に従うGaussian群
  fft/                 残差のFFT、周期候補の検索
  autocorrelation/     残差の線形自己相関
  common/              共通モデル・ROOT/TSV IO・品質判定
  compare/             各手法のROOTを読み、差・比・図を作成
  workflow/            CLI、channel群分割、Slurm投入
  Makefile, run.sh, slurm_array.sh
```

Gaussian群の数値計算は `common/CombFit.cc` を両fitで共有する。全体分布のモデルは `common/Numerics.cc` にあり、LEDは平滑化、muonはLandau⊗Gaussian×threshold、electronは指数関数＋定数×thresholdを使用する。滑らかな分布のfitだけでは1 p.e. gainは決まらず、等間隔のピーク・周期が必要。

## ビルド

プロジェクト直下から、ROOTの環境が有効な端末で実行する。

```bash
cd /home/murata_t/ScECAL_BeamTest
make -C analysis/sipm_gain -j2
./analysis/sipm_gain/Analyze --help
```

`bash analysis/sipm_gain/run.sh ...` でも差分ビルド後にAnalyzeを実行できる。

## 出力の配置

```text
analysis/result/sipm_gain/<label>/
  inputs/manifest.tsv, source/, analysis_plan.json
  data/
    led/<sample>/
      collection/spectra.root
      multi_gaussian/part_0000_of_0004/gain.root
      envelope_comb/part_0000_of_0004/gain.root
      fft/part_0000_of_0004/gain.root
      autocorrelation/part_0000_of_0004/gain.root
      ... part_0001～0003
    electron/<sample>/... 同じ構造
    muon/<sample>/... 同じ構造
  comparison/
    all_estimates.root, all_estimates.tsv
    method_comparison.tsv, cross_sample.tsv
    summary.json, PREPARE_COMPLETE, NUMERIC_COMPLETE
    inputs/configuration.json, source/
    prepared/layer_00.root ... layer_29.root
    layers/layer_00/numeric.root, *.tsv, NUMERIC_COMPLETE
    plots/layer_00/comparison.root, channel_*.png, plot_spectra.tsv, PLOTS_COMPLETE
    ... layer_01～29
```

LEDのsampleは `sps_led_20230428_0829_Calib_dac2100` のように **PS/SPS × 測定日時 × DAC** とする。electron・muonは従来通りrun別。各partは独立した出力・ロック・完了印を持つ。共有ROOTへの同時書込みは行わない。

各partには `gain.root`（数値TTree `Gain`）、`gain_estimates.tsv`、詳細候補の`fits.json`、分布・fit・残差・FFT・自己相関の`diagnostics.root`、`figures/`、条件記録と`FIT_COMPLETE`がある。plot-limitはpartごとの上限、0は全スペクトルをPNG化する。

統合は各partの **gain.root** を読み込む。完了印・checksum・設定・入力ROOT・sourceの一致を検証し、重複した結果を拒否する。`all_estimates.root`には `Gain`・`MethodComparison`・`CrossSample` の3個のTTreeを保存する。巨大な比較表を1個のJSON文字列にはしない。ROOT・TSVとも逐次出力し、channelごとの比較時に保持するのはそのchannelのレコードだけとする。準備時のROOT索引はlayer単位、重複検査は入力単位で行う。元の手法別ROOTも保持する。

`Gain`にはreference関連の数値欄も保存し、指定がない値はNaNとなる。出力の行順はlayer・channel順へ変わるが、gain値・手法間の差/比・温度が同距離の場合の対応の選び方は維持する。`metadata_json`は条件などの小さな情報だけを含む。

比較コードは `inputs/source/` に別途記録し、fit時のソース記録は更新しない。完了印とchecksumで中間結果を再利用する。条件・元入力・比較コードを変更した場合や旧形式の比較結果がある場合は、`--comparison-name NAME`で同じlabel内の別の比較先を指定する。fitの再計算は不要。

## 入力と解析条件

LEDのdac/adacはLED印加電圧の内部設定値であり、Vへの換算は行わない。group→channelの対応は`common/LEDGroups.hh`で共有し、`analysis/led/between_files_LED.cc`と全14 groupで一致する。decodeに値があっても点灯対象外のchannelは収集時に除外する。

2023年PS/SPSのdecode ROOTを使う。LEDはSPSの`LEDCali`・PSの`lightCalibData`、electronは`e-`、muonは`mu-`からmanifestを作る。manifestはタブ区切りの `sample kind beam energy_GeV period led_group dac settings raw`。

- 入力manifestは元ファイルごとのまま保持し、`inputs/analysis_units.json`に日時・DAC別の入力一覧を保存する。`--row N`とSlurm taskのrowは、この解析単位の0始まりの番号。group不足は一覧にそのまま残し、重複group・異なるsettingsの混在は拒否する。
- LEDの14 groupは対応するchannelが重ならないため、各channelの分布を同じROOTへ収める。異なるchannelのヒストグラムを加算しない。`sid`はROOT内で振り直し、`source_sid`・`source_index`・`led_group`をSpectraに保持する。`sources_json`・`source_collections_json`に元ファイル・sampling・checksumを残す。fitの乱数seedには元のsidを用い、再配置でnull診断が変わらないようにする。
- `gain.root`にも`source_sid`・`led_group`・`source_sample`・`source_raw`を保存する。
- memory cellは **raw CellIDの百・千の位が0**、すなわち `(CellID % 10000) / 100 == 0` のhitだけを使用する。
- 物理IDは `layer*100000 + chip*10000 + channel`。接続された30 layer×210 channelを扱う。
- run/channel/HitTag/GainTag/0.5℃温度区間を分離する。既定fitはHitTag=1。温度欠損は区別して保存し、温度補正はしない。
- 既定は全channel・全event。`--channels`、`--events-per-file`で制限できる。event数制限時はファイル全体から等間隔に選ぶ。
- ADC 250～1400内の0.5～99.5%点から最大400 ADC幅を選択し、6～40 ADC/p.e.を検索する。最低120 entry。channelごとに範囲の妥当性を図で確認する。
- channel群は物理channelの通し番号 `(layer*210 + chip*36 + channel) % groups` で決める。同じchannelの温度・タグ違いは同じ群に入る。group番号は0始まり。

## bootstrapなしの意味

bootstrapは既定0で、正の指定は拒否する。event block配列の作成・100回の再fit・bootstrap共分散は実行しない。従来のblock-size、memory-mode、min-blocks、max-rel-errorオプションはこの構成では使用しない。

- `stat_error = NaN`、`uncertainty_status = not_evaluated`を明示する。誤差0という意味ではない。
- Gaussian fitの`fit_error`はMinuitの条件付き誤差として別欄に保存し、方法間の統計的な一致判定には使わない。
- `candidate`は周期・fit形状・ADC範囲変更などの品質条件を満たした候補。calibration採用値ではない。`review`は理由付きの要確認、`insufficient_data`は統計・範囲不足。
- 方法間の差・比は出力するが、pull・共分散・統計的整合性は未評価。比較図に統計誤差棒は付けない。塗りつぶしはcandidate、白抜きはreview。
- **滑らかな帰無分布から作る199組の擬似データは残す**。これはFFT/自己相関の周期検出用で、bootstrapやGaussian再fitは行わない。`null_p`は推定した全体分布に条件付けた診断で、多channelの多重検定補正は含まない。
- ADC範囲・bin幅を変える少数回の再fitも残し、変化を`stability_ADC`に保存する。完全な系統誤差ではない。
- settings未確認・温度の違いは別欄で示す。gainの単純平均や自動校正は行わない。bootstrapなしでは温度傾きの誤差付きfitも行わない。

## 作成済みのADC分布をまとめ直す

旧配置の`full_split`を残し、新しい`full_grouped`へ日時・DAC単位で統合する。rawイベントの再読み込みやfitは実行しない。electron・muonのADC ROOTは内容を変えずコピーする。

```bash
./analysis/sipm_gain/Analyze regroup \
  --from-label full_split --label full_grouped
```

実行前の件数確認は`--dry-run`、代表1組だけの処理は`--row 0`を追加する。完成した統合単位は入力条件・checksumを確認して再利用する。単純な`hadd`はhistogram名の衝突があるため使用しない。全体はLED 128組＋electron 43 run＋muon 49 run＝220単位となる。

## 段階別のSlurm実行

**既定stageはcollectのみ**。各段階の結果を確認してから次を実行する。現在の`full_split`から上記regroupを完了した場合は、①をやり直さず②へ進む。

① 新規にrawからADC分布を作る場合（labelは新しいものを指定する）：

```bash
./analysis/sipm_gain/Analyze submit \
  --manifest analysis/result/sipm_gain/full/inputs/manifest.tsv \
  --label example_grouped --stage collect \
  --channels all --events-per-file 0 \
  --max-concurrent 16 --memory 32G --walltime 24:00:00 --execute
```

LEDは日時・DACごとに1taskが対応するgroupファイルを順に読み、対象channelだけを統合する。収集用のgroup別中間ファイルは成功後に除去し、入力条件・sampling記録は最終ROOTと`inputs/`に残す。

② 全ADC分布の完成後に、手法別gainを算出する：

```bash
./analysis/sipm_gain/Analyze submit \
  --manifest analysis/result/sipm_gain/full_grouped/inputs/manifest.tsv \
  --label full_grouped --stage fit \
  --channel-groups 4 --max-concurrent 16 \
  --memory 32G --walltime 24:00:00 --execute
```

③ 全fit終了後に数値を統合・比較する。既定の `numeric` はprepare→analyze→mergeを逐次実行し、描画は行わない：

```bash
./analysis/sipm_gain/Analyze compare --label full_grouped
```

④ 数値比較を確認した後に描画する：

```bash
./analysis/sipm_gain/Analyze compare --label full_grouped --stage plots
```

### 比較をさらに分割して実行する場合

各コマンドの完了を確認してから次へ進む。bootstrapは引き続き実行しない。

```bash
# ③-a: 全fitの検証、30 layer別の中間ROOTとchannel索引の作成
./analysis/sipm_gain/Analyze compare --label full_grouped --stage prepare

# ③-b: 数値比較。--layerを省くと全30 layerを逐次処理する
./analysis/sipm_gain/Analyze compare --label full_grouped --stage analyze --layer 0

# ③-c: 全30 layerのNUMERIC_COMPLETEが揃った後に集約する
./analysis/sipm_gain/Analyze compare --label full_grouped --stage merge

# ④: 描画だけ。--layerを省くと全30 layerを逐次処理する
./analysis/sipm_gain/Analyze compare --label full_grouped --stage plots --layer 0
```

`prepare`だけが元のfitを読み、全4成果物のchecksum・設定・collection・sourceを照合する。以後はその時点の検証済み中間ROOTを使用する。元のfitを更新して比較し直す場合は新しい `--comparison-name` でprepareから実行する。referenceや温度許容幅はprepare時に指定し、後段は保存済みの条件を使う。途中状態を許す場合はprepare時に `--allow-partial` を指定し、完成後の全体比較とは比較名を分ける。

別layerのanalyze/plotsは並列実行できる。同じlayer・同じ工程への同時書込みはロックで拒否する。layerのanalyzeが完了すれば、そのlayerのplotsを実行できる。全体のmerge完了は描画の必須条件ではない。

### layer比較のSlurm投入例（最大32並列）

ROOT環境が有効な端末で実行する。ビルド後、prepareをSlurmの単一ジョブで行う場合は次を使う（同時にanalyzeは投入しない）。

```bash
cd /home/murata_t/ScECAL_BeamTest
mkdir -p analysis/log/sipm_gain/full_grouped/compare
sbatch --job-name=sipm_compare_prepare --partition=all \
  --cpus-per-task=1 --mem=8G --time=08:00:00 \
  --chdir=/home/murata_t/ScECAL_BeamTest \
  --output=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain/full_grouped/compare/prepare-%j.out \
  --error=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain/full_grouped/compare/prepare-%j.err \
  --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0; ./analysis/sipm_gain/Analyze compare --label full_grouped --stage prepare'
```

prepareの正常終了と `comparison/PREPARE_COMPLETE` を確認後、次で数値比較を投入する。30 layerなので最大同時実行数は30となる。数値比較の成果物は各 `layers/layer_XX/` に分かれる。

```bash
cd /home/murata_t/ScECAL_BeamTest
mkdir -p analysis/log/sipm_gain/full_grouped/compare
sbatch --job-name=sipm_compare --partition=all \
  --array=0-29%32 --cpus-per-task=1 --mem=8G --time=08:00:00 \
  --chdir=/home/murata_t/ScECAL_BeamTest \
  --output=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain/full_grouped/compare/analyze-%A_%a.out \
  --error=/home/murata_t/ScECAL_BeamTest/analysis/log/sipm_gain/full_grouped/compare/analyze-%A_%a.err \
  --wrap='export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ROOT_HIST=0; ./analysis/sipm_gain/Analyze compare --label full_grouped --stage analyze --layer "$SLURM_ARRAY_TASK_ID"'
```

partition名は利用環境に合わせる。8G・8時間は初期設定例であり、全規模での必要量は未測定。全taskの完了印を確認したら `--stage merge` を別途実行する。描画も同様に `--stage plots` を使う独立した配列ジョブにできる（ログ名は `plots-%A_%a.out/.err` にする）。自動的な次工程の投入は行わない。

`Analyze submit` の `--execute` を省くと準備だけで投入しない。上記の直接の `sbatch` コマンドは実行時にジョブを投入する。①から新しく収集した場合、②・③も①と同じlabelを使用する。

- fitは **解析単位×4手法×channel群**。全体220単位・4群なら3520 task（LED 2048、electron/muon 1472）。各channelで必要なfit回数は変わらず、起動・ファイル管理の負担を減らす。
- Slurmの配列は1000 task以下のwaveに分け、`afterok`で順につなぐ。同時に走るのは最大`--max-concurrent` task。group統合で1taskの担当channelは増えるため、タスク数の減少率がそのまま高速化率にはならない。
- 同じlabelへのcollect→fitの追加投入が可能。fitを投入済みのlabelへの重複投入は拒否する。旧配置のlabelは`regroup`で新しいlabelへ移してから新コードでfitする。
- ログは `analysis/log/sipm_gain/<label>/<stage>/wave_*/slurm-%A_%a.out/.err`、投入記録は`jobs.tsv`。`runtime/`に固定した実行ファイル・ソースを使う。コードを変えた場合は別labelを使う。
- `--stage all`を明示した場合だけ収集からfitまで連続実行する。比較は常に別コマンド。
- 途中比較は`compare --allow-partial --comparison-name partial`で不足partを明示する。

## 小規模なローカル実行

```bash
./analysis/sipm_gain/Analyze discover --profile pilot \
  --output analysis/result/sipm_gain/example/inputs/manifest.tsv
./analysis/sipm_gain/Analyze run \
  --manifest analysis/result/sipm_gain/example/inputs/manifest.tsv \
  --label example --stage collect \
  --channels 20010,20011,420010,420011 --events-per-file 3000
./analysis/sipm_gain/Analyze run \
  --manifest analysis/result/sipm_gain/example/inputs/manifest.tsv \
  --label example --stage fit --channel-groups 2 --plot-limit 1
./analysis/sipm_gain/Analyze compare --label example
./analysis/sipm_gain/Analyze compare --label example --stage plots
```

`run`はローカルで逐次実行する。`run --stage fit --row 0 --method fft --channel-group 0 --channel-groups 2`で1単位・1手法・1channel群だけを指定できる。Slurmのfit投入は4手法まとめて行う。

完成partは条件とchecksumが同じなら再利用する。`--refit`で指定partを更新、`--recollect`でrawから再収集する。異なる条件を並存させる場合は別labelを使う。collectionの更新は古いfitとのchecksum照合で検出する。

旧`pilot`、`pilot_cpp`、`full`、`full_retry`は旧条件の記録として保持する。`regroup`で流用できるのは、memory cell 0のみ・bootstrap保存なしを明示する`full_split`形式の完成したcollection。

## 実装確認の結果

`analysis/result/sipm_gain/split_validation/`に小規模検証結果を保存した。PS/SPS LED各1ファイル、SPS electron・muon各1ファイル、各3000 event、8 channel、2 channel群で、収集4 taskとfit32 taskのSlurm実行が完了。60行の手法別結果をROOTから統合し、欠落・重複がないことを確認した。旧`pilot_cpp`と比較できる24件のgain中心値は全件一致した。

既知のgain 15.2 ADC/p.e.を持つ合成分布では、3種の分布×4手法で差が0.03 ADC以内だった。memory cell 1の除外、bootstrap指定の拒否、ROOTの読み書き、入力変更時の検出も確認した。これは処理の動作確認であり、実データの全channelでgainが決定できたことを意味しない。

- 検証記録：`split_validation/inputs/verification.json`、`synthetic_validation.json`、`negative_checks.json`
- 比較結果：`split_validation/comparison/`（旧結果との差は`migration_comparison.tsv`）
- 検証ログ：`analysis/log/sipm_gain/split_validation/`
- 旧group別収集の投入記録：`analysis/log/sipm_gain/full_split/`（元データの来歴として保持）

`inputs/analysis_plan.json`はfit開始時に作成される。投入前の設定は上記`submission.json`と各waveの`tasks.json`で確認できる。

### LEDの日時・DAC統合の確認

`analysis/result/sipm_gain/grouping_validation/inputs/`に検証記録を保存した。

- `merged_histogram_check.json`：14 group・6425スペクトルの全ADC bin、bin誤差、entry数、温度・タグ情報が統合前と一致し、sidの重複がない。
- `gain_comparison.json`：小規模4入力・60手法レコードで、24件の有限gainとnull診断・安定性・品質判定が旧結果と一致。
- `direct_collection_check.json`：新規raw収集でも14 group・4選択channel・各25 eventのgroup/memory選別とキャッシュ再利用を確認。
- `negative_checks.json`：重複group、settings混在、元labelへの上書きを拒否し、stage省略時はcollectのみであることを確認。

今回の全体統合先は`analysis/result/sipm_gain/full_grouped/`。元の`full_split/`は保持する。全体統合の記録は`full_grouped/inputs/regroup_verification.json`、並列統合ログは`analysis/log/sipm_gain/full_grouped/regroup/`。全体fitの投入前に作成した準備記録は `full_grouped/inputs/fit_prepare.log` に残している。


### 比較処理の省メモリ化の確認

`grouping_validation/inputs/comparison_streaming/verification.json` に検証結果と測定ログを保存した。

- 既存実データ60 gain・90手法間比較と、18個の描画グラフが旧版と一致。
- 合成入力1,200 gain・1,800手法間比較・1,312 run間比較と、96個の描画グラフが旧版と一致。温度の同距離、NaN/非有限値、reference、複数ページを含む。
- 30 layerを4並列で処理し、逐次実行と同じ数値を得た。空layer・空比較表、キャッシュ再利用、ロック、checksum不一致、未完成の統合、条件変更の拒否を確認。
- 21,000件と210,000件のROOT入出力で、比較時に保持するレコード数が1 channel分の1,000件に収まることを確認。全規模の処理時間・最大メモリは未測定。
- 実データの新比較結果は `grouping_validation/comparison_streaming/`、合成入力と旧・新比較結果は `comparison_fixture/` に保存。検証時の一時プロジェクトから保存先へコピーし、出力checksumを再確認した。ログ内の `.scratch/` は検証時のパスで、合成入力の再現用ソースは `comparison_fixture/inputs/fixture_generator.cc` に保持する。

`full_grouped` の全体比較ジョブはこの検証では投入していない。既存fitを再実行せず、上記prepareから開始できる。
