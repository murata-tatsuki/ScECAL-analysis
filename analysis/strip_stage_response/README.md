# MC truth → digitization → calibration の物理ストリップ比較

threshold fit品質による共通channel maskを指定できる。`beam_quantitative` で選んだ**同じイベント**について、generationの `MC_Truth`、digitization後の `Raw_Hit`、calibration後の `Calib_Hit` を比較する。SSAで分割された小さいhitの数ではなく、各イベントで異なる**物理ストリップの数**を数える。MCで深い層のストリップがdataより多い問題を、工程ごとの削除・追加・信号応答に分けて調べる解析。

## 実行

先に同じlabel・condition・energy・phaseの `beam_quantitative` とstrip診断を完了させる。既定は `beam_scan` / `point` / `allChannel` / `100 GeV` / `scan`。beam解析の入口4層の重心による位置補正とaxis-valid選択をそのまま使い、工程ごとにイベントを選び直さない。CoG cutは追加しない。beam解析に含まれていないenergyを調べる場合は先にそちらを実行する。

```bash
cd ~/ScECAL_BeamTest/analysis/strip_stage_response

# ローカル実行。入力ROOTの内容は変更しない。
bash run.sh --condition point --energy 100 --threshold-cut allChannel
bash run.sh --condition point_noTriggerSci --energy 100 --threshold-cut allChannel

# thresholdのfit品質でchannelを選ぶ。MC再生成は不要。
bash run.sh --condition point --energy 100 --threshold-cut precision
bash run.sh --condition point --energy 100 --threshold-cut chi2_10

# 読み取りのみの事前確認
bash run.sh --condition point --energy 100 --dry-run

# Slurmで各条件を投入
bash submit.sh --condition point --energies 100
bash submit.sh --condition point_noTriggerSci --energies 100

# 複数cut × energyを一括投入。各cutは別ジョブで同時実行できる。
bash submit.sh --condition point --energies 100 \
  --threshold-cuts allChannel,fitValid,precision,chi2_10,chi2_5,chi2_3
bash submit.sh --condition point_noTriggerSci --energies 100 \
  --threshold-cuts allChannel,fitValid,precision,chi2_10,chi2_5,chi2_3

# 任意の数値をshellから指定。出力先はcut名で区別する。
for chi in 20 10 8 5 3; do
  bash submit.sh --condition point --energies 100 \
    --threshold-cut "custom_chi2_${chi}" --threshold-max-chi2ndf "$chi"
done

# 組合せcutの例
bash submit.sh --condition point --energies 100 --threshold-cut my_quality \
  --threshold-max-chi2ndf 8 --threshold-max-rel-threshold-error 0.05 \
  --threshold-max-rel-sigma-error 0.3 --threshold-min-sigma-adc 0.6 \
  --threshold-max-sigma-adc 29.5

# 複数energy。全energyを調べる場合は解析済みenergyをカンマ区切りで指定する。
bash submit.sh --condition point_noTriggerSci --energies 40,100

# 保存済み集計ROOTから再描画
bash run.sh --condition point --energy 100 --threshold-cut precision --plots-only
```

異なるcondition/cut/energyは同時に実行できる。同じcondition/cut/energyは結果側のlockで重複書込みを防ぐ。ビルドは共通lockを使う。照合済みMCイベントは `allChannel` のROOTを共有し、cacheの作成・更新のみ排他lock、解析中は共有読み取りlockを使う。各cutの集計・描画は別ファイルなので並行実行できる。同じ条件を通常実行し直すと同じ出力先を更新する。特別なoverwriteオプションは不要。`--plots-only` は既存の層別・深さ別数値を変更しない。SiPM別histogramのない旧結果では、保存済みイベントcacheと保存済みchannel maskからSiPM別の表示用分布だけを追加する。threshold fitを読み直さず、イベント・入口軸・maskを維持する。cacheのSHA256が保存記録と違う場合は停止する。

generation/calibrationのパスと選択entryは既存beam解析の `configuration.json` から読み取る。digitizationはgenerationのパスの `generate` を `decode` に置き換えて同じファイル名を探す。異なる配置の場合は `--mc-digitized-dir /path/to/100GeV` を指定する。beam解析の入力が変更されていたら、古い選択結果との混在を避けるため停止する。その場合はbeam解析を再実行してから使う。

`--label`、`--phase`、`--beam-result-dir`、`--result-dir`、`--bootstrap`、`--seed`、`--mip-energy-mev` も指定できる。既定のbootstrapは500回、seedは20261011、1 MIPは0.305 MeV。truthのMIP表示はdigitizerと同じ名目値でEdepを割った値であり、測定したADCではない。

## Threshold channel cut

これは**fit品質によるchannelの除外**であり、hit signalにthresholdをかけ直す処理ではない。MC truth・digitization・calibrationとdataに同じ物理CellID集合を適用する。入口軸・axis-valid判定・イベント数は元のbeam解析から固定し、cut後に0 hitになったイベントも分母から落とさない。元の11種類の図・MIP cut・定義を維持する。

既定のfit入力は `analysis/result/threshold/threshold.root` の `expErfThre`。現在のdigitizerが読み込むモデルと一致する。`--threshold-root FILE --threshold-tree TREE` で外から指定できる。他のfitモデルを選ぶ場合、それは同じMCに別の品質maskを適用した比較であり、MCのthresholdパラメータを変更した比較ではない。

`fitValid` の共通条件は、接続済み6300 channel中、fit treeに個別channelの行があり、threshold/sigma/誤差/χ²ndfが有限、threshold>0、sigma>0、誤差≥0、χ²ndf≥0であること。ROOT生成コードは成功fitだけを保存しているが、fit status・共分散・NDFそのものはtreeにないため独立に再確認できない。個別fitがないchannelは除外する（MC側でchip平均等の補完値を使用していても含めない）。

| cut名 | fitValidに追加する条件 | 現在のROOTで選択される本数 / 6300 |
|---|---|---:|
| `allChannel` | channel除外なし。fit ROOTを参照しない | 6300 |
| `fitValid` | 共通条件のみ | 5959 |
| `precision` | threshold相対誤差≤5%、sigma相対誤差≤30%、0.6≤sigma≤29.5 ADC | 5636 |
| `chi2_20` | χ²ndf≤20 | 4733 |
| `chi2_10` | χ²ndf≤10 | 3441 |
| `chi2_5` | χ²ndf≤5 | 1767 |
| `chi2_3` | χ²ndf≤3 | 704 |
| `chi2_2` | χ²ndf≤2 | 209 |
| `quality_loose` | χ²ndf≤10、threshold相対誤差≤10%、sigma相対誤差≤50%、0.6≤sigma≤29.5 ADC | 3234 |
| `quality_medium` | χ²ndf≤5、threshold相対誤差≤5%、sigma相対誤差≤30%、0.6≤sigma≤29.5 ADC | 1628 |
| `quality_strict` | χ²ndf≤3、threshold相対誤差≤2%、sigma相対誤差≤20%、0.6≤sigma≤29.5 ADC | 608 |

まず `allChannel → fitValid → precision` で「個別fitがない・不正」「誤差が大きい・sigmaがfit境界に張り付く」channelの影響を分け、その後 `chi2_10 → chi2_5 → chi2_3` でfit形状の影響を見る。σのfit探索範囲は0.5–30 ADCなので、precisionの0.6–29.5は境界付近を除く診断用条件。これらは検証用の候補で、品質を保証する最終基準ではない。

現在のχ²ndfの中央値は約8.3。`chi2_2` はlayer 3/17/21が空になるため最初の標準比較には勧めない。厳しいcutでMC/data差が減っても、中心のchannelを多く除外した結果かもしれない。channel本数だけでなく以下のmask図を一緒に確認する。

- `threshold_channel_coverage.png`：各層の接続channel数210と選択channel数。
- `threshold_channel_mask_even.png` / `threshold_channel_mask_odd.png`：30層の物理strip位置。0=除外、1=採用。strip pitchのbinを使用。event occupancyの図ではない。

上記3図は `allChannel` 以外で追加する。元の11図は同じ名前で出力する。選択結果・層別本数・空の層・各条件の定義はsummary/configurationに、6300 channelそれぞれのfit値と除外理由は `inputs/.../threshold_channels.tsv` に保存する。fit tableのSHA256を記録する。

任意のcut名を `--threshold-cut NAME` で指定し、数値条件を追加できる。未登録名には少なくとも1つの数値条件が必要。presetに指定した数値は上書きされるので、条件を変える場合は区別できる名前を付ける。`allChannel` に品質条件を指定することは禁止。

| 数値オプション | 条件 |
|---|---|
| `--threshold-max-chi2ndf` | χ²ndfの上限 |
| `--threshold-max-rel-threshold-error` | thresholdErr/thresholdの上限 |
| `--threshold-max-rel-sigma-error` | sigmaErr/sigmaの上限 |
| `--threshold-min-threshold-adc` / `--threshold-max-threshold-adc` | threshold値の下限/上限 [ADC] |
| `--threshold-min-sigma-adc` / `--threshold-max-sigma-adc` | sigma値の下限/上限 [ADC] |
| `--threshold-max-sigma-over-threshold` | sigma/thresholdの上限 |

数値条件は全てAND、境界値を含む。最後のthreshold値やsigma/thresholdの条件はfit品質以外に実効thresholdの分布を選ぶため、まず品質cutを単独で比較してから使う。presetの本数は現在の入力ROOTでの値で、入力変更後はdry-runの `channel_selection` を確認する。`--plots-only` は保存済みmaskの集計を描画し、fit ROOTを読み直してmaskを変更しない。

## 入力と照合

- truth：`generate/.../<energy>GeV/*.root` の `MC_Truth/CellID,Hit_Energy`。初期electronのPDG/energyも確認。
- digitization：`decode/.../<energy>GeV/*.root` の `Raw_Hit/CellID,HitTag,HG_Charge,LG_Charge`。
- calibration：`calib/.../<energy>GeV/*.root` の `Calib_Hit/CellID,Hit_Energy`。
- data：既存beam解析の `root/data/<phase>/reference/<energy>GeV/pre_ssa_events.root`。SSA直前の校正済み物理ストリップを比較基準とする。

ファイルごとのentry数を一致させ、MCの `Run_Num,Event_Time,TriggerID,CycleID` とcalibrationのevent情報を照合する。現在のMCでは時刻・trigger・cycleがtruthのentry番号。異なるID方式には自動で対応せず停止する。重心のvalid判定・位置と校正後ストリップの集合・signedエネルギーも既存beam cacheに対して確認する。

CellIDのmemory-cell桁を取り除いて物理stripを特定し、30層×210の接続済みstripを対象とする。同じイベント・同じ物理stripが複数行にある場合、truth/calibrationは**signed energyを合算**して1本とする。digitizationは存在を1本とし、HitTagの有無を記録する。ADCはphysical stripごとの最大値を保存するが、ADCの図にはrawが1行のstripだけを使う。重複行・未接続ID・非正信号・HitTag=0・前段に存在しないstripをsummaryに記録する。校正後の負信号を黙って除いて本数を減らさない。

図中の3段階の本数は以下。

- **MC truth**：接続済みstripに合算したdeposited energyが正のもの。
- **MC after digitization**：保存済み `Raw_Hit` に存在する接続済みstrip。現在のMCではHitTag=1のみ。
- **MC after calibration**：保存済み `Calib_Hit` に存在する接続済みstrip。追加MIP cutなし、非正信号も本数には含む。

残存率は前後両方に存在するstrip数／前段strip数。前段にないstripが追加されても残存率の分子に入れず、追加本数として別に記録する。平均本数・energyは全acceptedイベント（その層が0本のイベントも含む）で割る。95%区間はイベント単位のbootstrap。MCの3段階には同じbootstrap indexを使い、dataは独立に再標本化する。

## 図の見方

凡例は図上部の専用領域に配置する。通常のstage比較では灰色=truth、青=digitization、赤=calibration、黒=data。digitizationとcalibrationが一致すると曲線は重なる。

| PNG | 何を見るか |
|---|---|
| `stage_strip_counts.png` | 左：全30層のstage別物理strip数/event。右：問題になっているlayer 20–28を拡大。どの工程で本数が減るか確認。 |
| `calibrated_strip_excess.png` | 校正後MC/data−1。追加cutなし／0.5／1／2 MIPを重ねる。高い信号のstripにも過剰が残るか確認。 |
| `stage_strip_retention.png` | 左：truth→digitization、digitization→calibration、truth→calibrationの残存率。右：各工程で消えた本数と、前段にないのに現れた本数。 |
| `stage_energy_response.png` | 左：全truth deposit、digitizationを通過したstripに元々あったtruth deposit、MC/dataの再構成positive energy。右：同じtruth-positive strip上のsigned再構成energy合計／truth deposit合計。選択によるenergy lossと応答変化を分ける。 |
| `stage_strip_width.png` | layerごとの細かい方向（奇数層X、偶数層Y）の本数重みRMS。共通の入口重心を差し引いてから分布自身の平均回りに計算。幅がどの工程で変わるか確認。 |
| `truth_signal_survival.png` | 4列×2行。上段10 µm（全4–27／4–9／10–19／20–27）、下段15 µm（全0–3・28–29／0–3／28–29／29）。横軸truth deposit/名目MIP、縦軸その信号binにいたtruth stripが各工程を通過した割合。95%区間は元の設定でevent bootstrap。 |
| `stage_signal_spectra.png` | 2列×2行。上段10 µm、下段15 µm。左：全truth／digitizationを通過／calibrationを通過したstripの元のtruth deposit。右：校正後MCとdataのsignal。各SiPM種類・各曲線を保存histogram内のstrip総数で個別に規格化。可変幅binのfractionであり密度ではない。正領域のみ表示し、非正信号・範囲外はsummaryに残す。上下段は同じ軸範囲。 |
| `digitized_adc.png` | 2列×2行。上段10 µm、下段15 µm、左HG、右LG。各SiPM種類で規格化した保存済みADC分布。pedestal込み。Raw_Hitが1行のstripのみ。上下段は同じ軸範囲。飽和や出力域の異常を確認。 |
| `truth_calibration_correlation.png` | 左10 µm、右15 µm。同じstripのtruth depositと校正後signalの2D分布。色はstripエントリ数で、左右共通の色スケール。赤線は1:1。表示範囲外の数はJSONの`sipm.flows`に記録。 |
| `stage_strip_profiles_lab.png` | 30層の検出器座標の1D strip分布。各stage/layerで規格化。5.3 mm pitchで隙間のない表示binを使う。 |
| `stage_strip_profiles_aligned.png` | 同じ30層を入口重心で位置補正。MCの全stageに同じ補正値を使う。表示は約±111 mm、全範囲のhistogramはROOTに保存。 |

層をまとめる上記4図はSiPM cell sizeで分ける。15 µmはlayer 0・1・2・3・28・29・30・31、10 µmは残りの4–27。**現在の解析は接続channel定義と既存cacheに合わせてlayer 0–29を対象にするため、図の15 µmは0–3・28–29のみで、30・31は入らない。** 層別の図・元の全層/深さ別summary値・MIP cut・イベント選択は維持する。SiPM別分布は`histograms.root`の`sipm/10um/`、`sipm/15um/`、分類・範囲外件数・残存率は`metadata/sipm`と`summary.json`の`sipm`に保存する。正規化した分布は種類ごとの形状比較であり、10 µmと15 µmのhit総数の比較には使わない。SiPM種類とdepthも異なるため、差をSiPM応答だけの違いとは解釈しない。

**calibrationで本数が増えていなくても、energy応答が大きくなれば固定MIP cutを通る本数は変わる。** 本数・残存率とenergy応答・MIP cut比較を一緒に読む。

この解析は保存済みROOTの工程間の変化を特定する。`Raw_Hit` は既にthreshold選択後なので、棄却stripのADC、threshold乱数、threshold前のADC分布は復元できない。実効残存率はthresholdに加えpedestal/noise/SiPM応答を含む。dataにはtruthがないため、truthの本数がdataより多いという比較だけでGeant4やbeam入力が原因と断定しない。誤差はevent統計で、calibration等の系統誤差・fileごとの相関は含まない。

## 保存先

```text
analysis/strip_stage_response/                  コード・実行スクリプト
analysis/result/strip_stage_response/<label>/
  root/mc/<phase>/<condition>/<threshold-cut>/
    100GeV/
      events.root                              allChannelのみ：全channelの照合済みイベント
      histograms.root                          層別平均・95%区間・分布
      plots.root                               描画canvas
  plots/mc/<phase>/<condition>/<threshold-cut>/100GeV/*.png
  summary/mc/<phase>/<condition>/<threshold-cut>/100GeV/
    summary.json                               定義・注意点・audit・層別/深さ別数値
    metrics.tsv                                scope,metric,mean,low_95,high_95
  inputs/<phase>/<condition>/<threshold-cut>/100GeV/
    configuration.json, beam_configuration.json, threshold_channels.tsv
    cache_configuration.json, mc.tsv           allChannelのみ：共有cacheと抽出manifest
    complete.json, run.log, sources/
    beam_input_sources/                        既存allChannelの入力関連ソース控え
    plot_sources/                              再描画時のコードとchecksums
analysis/log/strip_stage_response/<label>/<phase>/<condition>/<threshold-cut>/
  100GeV-<jobID>.out, 100GeV-<jobID>.err          Slurm stdout/stderr
  submissions_<日時>_<PID>.tsv                  condition親にenergy/cut/jobIDを記録
```

`root/.../<condition>/<threshold-cut>/` の下はenergyディレクトリだけ、各energyディレクトリの中はROOTだけ。既存の無cut結果はユーザー指定に従って `allChannel` へ移した。maskを適用したイベントROOTのコピーは作らず、各cutのROOTはhistograms/plotsのみで、元の照合イベントROOTはallChannelに共有する。data cacheは既存beam解析を読み取り、新たなdata ROOTを重複作成しない。summary・PNG・設定・ログはROOTとは別系統に置く。元データのpath/size/mtime、beam cacheのSHA256、実行コマンド、解析コードを保存する。`beam_input_sources` は既存beam解析が保存した入力関連コードの控えであり、それだけで実際のproducerに使われたversionの証明にはならない。`complete=true` は抽出・集計・描画・入力更新検査が全て終了した時点で付ける。
