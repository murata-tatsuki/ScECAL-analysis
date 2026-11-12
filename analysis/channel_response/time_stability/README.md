# 時間・run・温度に対する応答安定性

C++/ROOTの `Collect` → `Analyze` → `Draw` を分離した診断です。既存muon追跡、現行校正定数reader、data温度補間、electron低LG tail判定を読み取りで再利用します。本番入力・校正・他解析は変更しません。

## 実行

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/channel_response/time_stability
make
bash run.sh --dry-run
bash run.sh --label pilot --energy 100 --muon-energy 100 \
  --muon-runs 29,30 --max-files 2 --max-events 250 \
  --event-bin 100 --block-events 100 --channels 20010,20013,420013
bash run.sh --label pilot --plots-only --plot-label redraw --channels 20013,420013
bash submit.sh --label batch_pilot --max-events 250 --dry-run
```

既定は各sample/particleの1ファイル・先頭1000 raw event。`--start-entry` と `--max-events` は連続raw entryを選びます。`--max-events 0 --max-files 0` は全量です。5/20 GeVには `--energy 5` 等を使い、muon energyは実在入力を確認して別指定します。不在入力をelectron MCで代用しません。全optionは `--help`。

`--collect-only` 後は同じlabelの `--analyze-only` で再開。保存configuration/manifest/flatのSHA256とschemaを確認し、入力・選別変更を拒否します。別解析のflatは受け付けません。plots-onlyは保存TSVのみを読み、数値ROOT/TSVの前後checksumを確認します。元のdecode/calib/truthを開きません。既存label・statistics・PNG先は上書きせず、plots-onlyにも新しいplot-labelが必要です。排他的directory作成とflockで出力を保護します。

Slurmは現行 `../execute.sl`、all、既定1 task/1 CPU、8 GB、2時間。`--partition/--mem/--time/--cpus` で変更可。`--cpus 4` は親1 process＋fit worker3 processを使います。dry-runは入力確認とsbatchコマンド表示だけ。logは `analysis/log/channel_response/time_stability/slurm-%j.{out,err}`、投入記録は同directoryのsubmissions.tsv。

必要branchのみ、scalar ID監査後のvectorは各入力を一度走査し全channelを同時集計。fit/描画は保存集計を使用。-O2、ROOT cache、低圧縮flat、FastPngが既定です。ROOT/TSVは親1 processだけが保存し、fit workerは分布の十分統計量だけを受け取るため入力の再走査・同一ROOT並列書込みはありません。channelsは数値を限定せず、そのchannelと対応chip/layerの図だけを限定します。

### メモリ上限と収集済みデータからの再開

`Analyze` は全hit/eventや区間の生ADC配列を保持しません。分布は固定binの度数、平均・幅・回帰は十分統計量で集計し、`--memory-mb 512`（既定）の集計バッファ上限で圧縮した一時ファイルへ退避します。最大16本ずつmergeし、同じ区間の値を統合してから従来と同じfitを行います。回帰はchannel/run切片を除いた最小二乗とblock score、相関は1 eventとblockごとの和、補正診断は検証blockのmomentを使います。選別・bin幅・分母・不確かさの定義は変更しません。浮動小数点の加算順は変わります。回帰の二次モーメントは中心化してmergeし、一定位置・角度やゼロ幅のrunで大きな数同士の引算による桁落ちを避けます。ROOT出力に8 MiB、数値表に各1 MiBの書込みbufferを使います。

`--memory-mb` はプロセスRSSそのものの上限ではありません。ROOT・branch cache・固定channel数の比較表・相関表・merge用bufferが別途必要です。Slurmの `--mem` は余裕を持って指定してください。退避先は出力内の固有 `.work-*` で、正常終了・通常の例外では自動削除します。強制終了で残った一時領域や不完全な成果物は自動削除・上書きしません。共有ディスクには一時領域も必要です。

既定の `--histogram-storage compact` は全分布の軸・entries・ROOT統計量と、非ゼロbinの番号/度数（flowを含む）を `statistics.root` の `histograms` TTreeへ保存します。記録のないbinは厳密に0です。分布や統計不足の区間は削除しません。`fits.tsv` の `histogram` は `histograms/entry_N` という参照です。`HistogramStore.hh` の `ts::readHistogram(file, reference)` でTH1Dに復元できます。旧形式は `--histogram-storage legacy` で `hist_0/`, `hist_1/`, …へ1024個ずつ保存できます。PNGはTSVだけを使い、trendsをメモリ制限付きの外部sort（1 process）で図単位に読み出します。図の点数がバッファを超える場合は間引かず、明示したエラーで停止します。

OOMなどで `statistics.root` が途中まで生成されている場合、元のlabelには上書きしません。`COLLECTED` と保存checksumを確認して、新しいlabelへ集計・描画だけを実行できます。

```bash
bash submit.sh --label muon100_run29_57_bounded \
  --reuse-collection muon100_run29_57 \
  --memory-mb 512 --channels all --mem 16G --time 24:00:00 --dry-run
bash submit.sh --label electron100_run105_107_bounded \
  --reuse-collection electron100_run105_107 \
  --memory-mb 512 --channels all --mem 16G --time 24:00:00 --dry-run
```

`--dry-run` を外すと投入します。元の入力を再収集せず、flat.rootは元結果へのsymlinkで参照します。元のlabelを保持してください。入力・選別・温度補正・bin定義は元のconfigurationを継承し、再利用時の変更は拒否します。由来は `inputs/reused_collection.tsv`、元のchecksumは `inputs/source_collection.sha256` に保存します。新label側にもconfiguration/manifest/flatのchecksumを作り、後続の再描画で検証します。

## 入力・baseline

electron dataは `ECAL_data/analysed/2023/sps/{decode,calib}/e-/<energy>GeV/`。MCは指定の `Result_MC/{decode,calib}/e-/e_ssa/threshold/<energy>GeV/`、truthは `Result_MC/generate/e-/e_ssa/<energy>GeV/`。

muonは既存muon_mip_responseの規則：data `{decode,calib}/mu-/<muon-energy>GeV/`、MC `{decode,calib}/mu-/mu_track/threshold/<muon-energy>GeV/` と `generate/mu-/mu_track/<muon-energy>GeV/`。両particleのenergyは独立です。

任意入力はinput-manifestの7列TAB区切り・絶対パス。header任意。同一raw重複拒否。manifestでは全行を使用し、max-filesは自動探索だけ、max-eventsは常に適用します。

```text
sample  particle  run  decode  calib  truth  evidence
0 muon 29 /absolute/raw.root /absolute/calib.root - -
1 electron 0 /absolute/raw.root /absolute/calib.root /absolute/truth.root /absolute/production_record.txt
0 pedestal 123 /absolute/dedicated_raw.root - - -
```

列間をTABにしてください。sample=0:data/1:MC。evidenceは生成記録など、未確認は `-`。

全raw/calibの `(Run_Num,Event_Time,TriggerID/Event_Num)` を監査し一意な保存IDでjoin。重複IDの全出現は応答から除外しentry順で補いません。MC manifest runとROOT内Run_Num=0を別保存。MC species/primary energyとcurrent MCDigiの `TriggerID=Event_Time=generated entry` 規則を照合します。これはID・species・energyの互換性確認で、歴史的production chainの証明ではありません。

production_audit.tsvは確認限界を残します。設定埋込みがなければCT適用・生成定数同一性はunconfirmed。path/mtimeでは断定しません。現行sourceも過去の実行sourceの証明ではありません。ROOT metadata、evidenceのstat、現行source、実使用定数snapshot/hashを保持します。

electron既定は `--data-tail 1 --mc-tail 0`。既存hg_lg_tail/threshold/e_ssa/<energy>GeVの完全判定をHLTailSelectionで読み、canonical calib pathと全保存IDを照合。tail-resultsで指定可。sample別checked/rejected/selectedを記録。muonにtail・electron energy・CoG cutは適用しません。decode→calibration hit残存率は調べません。

## 選別・時刻・温度

- channelはlayer0–29、chip0–4は36、chip5は30。maskの物理CellID一覧を両sampleに適用。既定は接続channelのみ、追加maskなし。
- full CellIDはmemoryを含む。muon追跡座標は物理stripにまとめ、予測stripに複数HitTag=1があればobservedは1回、ADC応答は除外。calib hit joinはfull ID。
- electron記録hit数はmask後のcalib record（負/ゼロ含む）。physical stripsはmemoryをまとめたdistinct数。rawはHitTag=1。SSA virtual hit数は未測定(-1)。e_ssaというpathだけでCalib_Hitをvirtual hitと解釈しません。
- electron主energyは非負Hit_Energyの和。signed sum、負/ゼロ数も保存。重複recordのenergyは合計。非有限energyのeventは応答未測定として選択から除外。追加energy windowなし。
- electron CoGはstrip中心と非負energyで求めるlayer9/10の各X/Y正方形cut。200 mmはnocut。shower重心は独立入射位置ではなく外部trackerはありません。
- timestampの最小正差・差のGCD・同値・後退・run/file境界を全scalarで監査。後退でsegment分割。観測刻みはclock精度を証明しません。decodeは4 byte整数保存ですが単位/epochは未確認、既定はtick/event順。
- 秒表示はseconds-per-tickとclock-evidence指定時のみ、segment開始からの経過秒。epochは推測しません。time-binは秒、tick-binはtick、event-binはraw entry幅。区間範囲/母数はintervals.tsv。空区間は行なしであり0測定ではありません。
- MC時刻はgenerated index。event順・production partitionとしてのrun・pooled比較だけ。MC実時間/実測温度trendは生成せず、20℃はnominal。
- 現行data EBUdecode.cxxをリンク。16センサーをstrip中心への逆距離で補間、距離0はセンサー値、重みを最小距離でrescale。decodeは生温度/128で℃。現行calibrationは空配列なら20℃。本診断は16個すべて有限かつ0<T<60℃を実測条件とし、欠損/代入を実測と扱いません。
- 20℃基準とF(T)=1-alpha*(T-20)を維持。alphaはlayer4–27:1.6/135、その他:3.5/230。補正は保存NewTemperature、独立sensor補間値もflatに保持。bias一定でもgainからCTを導きません。

## muon応答・統計

既存TrackGeometryで対象layer全体を除外してから非加重track推定。各投影7層以上、角度0.05 rad以下、RMS3 mm以下、45×5×2 mm strip全厚を長手1 mm・幅0.25 mm余白付きで通過。長さは2*sqrt(1+sx²+sy²) mm。分母は全予測、observedは0/1。他layerでtrack取得できた集団での条件付き割合でhardware単独効率ではありません。

ADCは現行CurrentCalibrationのpedestalを引き2/lengthを掛けます。before/afterはpedestalあり、raw HG<2200、full ID一意、係数有限・正、dataでは実測温度ありの同じhit集合。HG=4000をLGで代替しません。負/ゼロADCのみを理由に除外しません。MIP/HL採用・fallbackは現行規則で旧intercept条件をbaselineに戻しません。

ADCは共通5 ADC bin、[-100,2600]。threshold+3sigma、通過長、温度係数、HG上限から全hit共通fit区間を取り、Landau–GaussianをPoisson likelihoodでfit。MPVはLandau成分のmode補正済みMPV、Landau widthとGaussian sigmaは別保存。100 plateau hit未満、fit/covariance/quality/boundary不成立はNaNとstatus。全数・fit内数・flowを保持。chip/layerのADC fitはchannel混合で、gain差を除いた幅ではありません。

mpv_slopes.tsvは成功温度binのMPVからpooled/within-run傾き・誤差を算出。3点未満や幅0.5℃未満は未判定。位置/角度や時間相関を調整した係数ではありません。

slopes.tsvは同一hit ADCにchannel固定切片を除くpooled傾き、channel×run切片を除くrun内傾きを求め、長手/幅位置、sx/syを同時調整。温度幅不足・退化・自由度不足では確定値にしません。run_means.tsvにchannel平均を除いたrun別温度・応答。chip/layerでもchannel平均を除きます。交絡が完全に消えたとは断定しません。

残差はsample別channel平均を引き、平均ADCで割ったfractionも保存。before/afterはpairedで独立誤差から差の有意度を計算しません。平均推定に使った集団内の診断です。correlationは同じeventの観測pairのみから分母・covariance・rhoを算出。muonは1 event/layerに最大1予測stripなので同一chip/layer内の複数channel共通成分は識別不能です。異なるlayer/chipの残差・event平均残差にはdE/dx・追跡揺らぎも含まれます。electron layer energy相関はshowerを含み、回路noiseとは呼びません。

通常fit誤差とtrend SEMは独立近似。ADC温度傾きにはrun/file/raw-entry blockのcluster SE、相関にはblock jackknifeを併記。4 block未満はblock誤差NaN。時間不明時のblockはentry順で相関時間を測定したことにはなりません。

## 別系統と補正候補

beam rawのHitTag=0はmemory別HG/LG pedestal差分・幅として別集計。非hit条件付きであり無biasと仮定しません。専用rawはmanifestのpedestal行で追加可。既存ped/fit・ped/histを名称だけで専用runと認定しません。

`--gain-table /absolute/gain_estimates.tsv` で既存sipm_gain数値を再利用。値・誤差・method・period・DAC・memory・status・review理由をgain.tsvに保持し、温度に対する候補値をperiod/DAC/method別に描きます。review/unresolved/failedを確定値にせずgain slopeやbeam期間の変化を無理に求めません。time_bin/periodから実時間を作らず入力なしは未測定。gain→CT推定モデルはありません。

`--diagnose-correction 1` は時刻定義のあるelectronの経験的診断。split=blocksは連続blockを丸ごと交互に学習/検証、split=runsは偶数/奇数run。学習側だけでenergy対経過秒の直線を決め、固定モデルを未使用eventに適用。平均energy、RMS/mean、hit数、paired差を保存。hit再選択なし。誤差は検証block bootstrapで固定学習モデル条件付き（学習側の不確かさは含まない）。同一eventの平坦化を独立校正効果とは呼びません。時間・学習/検証統計不足・不正予測は未評価、既定は無効です。

## 出力と図

基点は `analysis/result/channel_response/time_stability/<label>/`。比較用ROOT/TSVは基点に置き、各行にsample/particle/group/axis/条件を保持。PNGはplots/<render>/<data|mc>/<muon|electron_exclude|electron_keep>/<CoG>mm/<energy>GeV/、muonはCoGなし。channel図はlayerN/chipN/channelN/、chip図はlayerN/chipN/、layer図はlayerN/。LEDはbeam energyなしでdata/LED_gain/period/DAC/layerN/chipN/。

| 出力 | 意味 |
|---|---|
| input_audit.tsv, run_boundaries.tsv | tick刻み/reset/run境界/ID重複/時計の確認状態 |
| counters.tsv, production_audit.tsv, root_metadata.tsv | 選別母数、対応、production/CTの確認限界 |
| flat.root | Events/Crossings/Pedestal。選択前event、同一hit、未測定を保持 |
| statistics.root, fits.tsv | 絶対数のADC/energy分布、MPV/幅/誤差/status |
| intervals.tsv, trends.tsv | 区間範囲、分母、MPV/幅/取得割合/pedestal/electron応答 |
| comparison.tsv | 同定義pooled data/MC母数・絶対値・差・誤差/status |
| slopes.tsv, mpv_slopes.tsv, run_means.tsv | run内/間傾き、交絡、温度幅/統計不足 |
| residuals.tsv, common_components.tsv, correlations.tsv | 平均差除去残差、観測可能な共通変動、pair数/block誤差 |
| gain.tsv, availability.tsv | LED fit候補と未測定項目 |
| correction_model.tsv, impact.tsv | 学習モデルと未使用event/run診断。未評価はNaN |
| plots/render/plot_index.tsv | 各PNGのpath・意味・点数 |

MPV_ADC/Landau_width_ADC/Gaussian_sigma_ADC/mean_ADC_vs_*は同一hit before/after（stage0/1）。hit_fraction_vs_*は予測通過分母。pedestal_*はmemory別HitTag0。channel_centered_fraction_vs_*は平均差除去ドリフト。electron energy/hit数/RMS_meanは絶対応答でshape規格化しません。RMS/meanはGaussian-core fit resolutionではなく通常誤差は近似。

correlation図はindex行列で物理配置図ではありません。色はdata/MC共通[-1,1]、灰色は未測定。physical_slope図は有効傾きがある場合だけ実strip座標・45×5 mm矩形、X/Y同縮尺、data/MC共通色範囲、未測定は灰色。slopes図は傾きと誤差または不足理由。held_out_impactは検証側before/afterまたは未評価理由。

小規模検証はsmoke_100GeVとvalidation_100GeVのrun.log/inputs/verification.tsvに記録。全量・全energy、Slurm実投入、本番校正更新は行いません。短区間の温度幅・channel統計不足から原因や温度係数は確定できません。

確認済みの図は `validation_100GeV/plots/checked/`（400 PNG）です。initial/redrawは初回検証の保存物です。そこでの相関mapはROOTがNaNを最小色で描くため未測定セルの判別に使わず、灰色の未測定表示を修正したcheckedを参照してください。数値ROOT/TSVは変更していません。人工入力ではcalib逆順ID join、時刻reset、明示時間単位、未使用blockでの補正診断も確認しました。詳細はinputs/additional_verification.tsvです。

5 GeV electronはSPS既定pathにありません。PS electronとSPS 100 GeV muonを併用するmanifest例をvalidation_100GeV/inputs/5GeV_manifest.tsvに保存し、pathのdry-runのみ確認しました。PSの校正定数の適合は別途確認が必要です。20 GeVは既定探索dry-runまでで実解析していません。

最終確認では、欠測thresholdがstd::maxのNaN比較で数値化されないguardも追加し、人工入力600 hitでfit範囲NaNを確認しました。実入力の定数12600行には欠測thresholdがなく、保存済み実入力結果の数値には影響しません。

## メモリ修正の検証

`memory_stable_100GeV/` は既存 `validation_100GeV` の収集済み実入力2,000 event（予測通過22,647件）を、既定512 MiBのバッファで集計・fit・400 PNGまで実行した結果です。最大RSSは588,368 KiB（約575 MiB）、終了コード0。`plots/checked/` は16 MiB設定のplots-only確認で、元の数値ROOT/TSVのchecksumを維持しています。

13種類の数値表を旧結果と比較し、108万行のtrend、67,364 fit、run内/間傾き・block誤差・相関を含めて一致を確認しました（数値の比較許容差は絶対2e-9＋相対2e-5）。68,232ヒストグラムの36,860,224 binはunderflow/overflowを含め完全一致です。人工入力では温度・位置・角度の回帰、一定幾何条件、block jackknife、未使用blockでの補正診断、16/512 MiBの退避条件を検証しました。400万eventのevent保持負荷試験は最大357,608 KiBでしたが、これはhit treeが空の人工入力であり全量の物理解析ではありません。

検証記録は同結果の `inputs/` に保存しています。OOMで止まったmuon/electron全量ジョブの再投入・完走は今回未検証です。元の `muon100_run29_57`、`electron100_run105_107` と収集済みflatは変更していません。


## 全量ジョブの進捗監査と時間制限対策（2026-10-14）

旧muonジョブ1328587は24時間の上限で停止しました。収集と77,059,121予測通過の集計は完了し、応答区間は43,797,948件中42,951,445件（98.0673%）まで数値出力されています。この件数はsample・channel/chip/layer・軸・補正前後ごとの分布数で、独立event数や計算時間の割合ではありません。実際のminimizer呼出しは414,729件、成功は198,970件でした。未処理のMC layerには少なくとも `n_all>=100` の673,708件があり、fit内数・共通plateau等の条件判定が残っています。回帰・残差/共通変動・pedestal・描画も未完了です。

監査は `analysis/result/channel_response/time_stability/timeout_audit_100GeV/` のprogress.tsvとrun.log。元の圧縮集計を全件走査し、途中TSVの改行で確定した行だけを数えています。98%を全解析の完了率や残り時間に換算しません。channel/event-orderの区間はdata 4,492,192件、MC 31,607,350件すべて `n_all<100` であり、記録数に比べてfit候補が少ないことも確認できます。

対策は以下を組み合わせます。

- `--reuse-fits SOURCE_LABEL`：収集済み入力と、元のfits.tsvの完了済みfitを新labelで再利用。入力checksum、保存Common/StreamingStats、ADC集計・fit定義を照合し、各再利用行でもキー・entries・flow・plateau範囲/数を検証します。欠測/統計不足は再判定し、末尾の未完行は無視します。元のROOTを修復・上書きしません。
- fit中は60秒または10万区間ごとにTSVをflushし、`inputs/analysis_progress.tsv` にphase・処理済み数・分母・新fit/再利用fit数を記録。時間切れ後も次の新labelから確定済みfitを再利用できます。収集後の集計・mergeと、回帰/残差/pedestalは再実行します。全工程の任意位置を復元するcheckpointではありません。
- `--cpus 4`：3 fit worker＋親1。workerは最大256分布のbatchから仕事を受け取り、raw/flatの再読込みとROOT/TSVへの書込みを行いません。`--fit-workers N` で指定可能で、N>1では親を加えたCPU数が割当を超える指定を拒否。1 worker指定は親内で逐次fitです。
- 既定はcompact ROOT。百万個の同じ疎分布の保存ループは旧214.453秒/358,856,120 byte、新6.57066秒/8,866,600 byteでした。初期化・終了を含む実時間は共有I/O待ちを含めて旧4:15.63、新1:18.68です。この単体比較にはAnalyzeの8 MiB書込みcacheは使っていません。これは保存形式の比較で、全解析の速度比ではありません。
- `--trend-axes standard`（新既定）：event順、確認できた経過秒、run、温度、pooledを集計。timestampの全scalar監査は維持し、raw timestamp ticksに対する応答trendは要求時だけ生成します。`--trend-axes all` は従来と同じtick区間も含みます。新規収集で `--tick-bin` を明示し、trend-axesを省略した場合はallです。基準温度・選別・各軸のbin幅は変えません。実行時設定は入力configurationを上書きせず `inputs/statistics_execution.tsv` とrun.logに保存します。

まず数値解析だけを、4 CPU・48時間の上限で再投入する例です。上限内での全量完走は未検証です。

```bash
cd /home/murata_t/ScECAL_BeamTest
bash analysis/channel_response/time_stability/submit.sh \
  --label muon100_run29_57_fast \
  --reuse-fits muon100_run29_57_bounded \
  --trend-axes standard --memory-mb 512 --channels all \
  --cpus 4 --mem 16G --time 2-00:00:00 --no-plots
```

`--dry-run` を追加すると投入せず設定・sbatch引数を確認します。数値解析のANALYZEDができた後、全channel描画を別jobで実行します。

```bash
bash analysis/channel_response/time_stability/submit.sh \
  --label muon100_run29_57_fast --plots-only --plot-label all_channels \
  --channels all --memory-mb 512 --cpus 1 --mem 16G --time 24:00:00
```

再びfit途中で止まった場合は、新しいlabelで `--reuse-fits muon100_run29_57_fast` を指定します。元labelと元のflat.rootは保持してください。plots-onlyには集計時のtrend-axesやfit-workersを指定せず、保存済みの数値を使います。

実入力2,000 event・22,647予測通過で、compact/legacyの13数値表（fitsの参照列を除く）と68,232分布の36,860,224binを比較しました。3 workerで新たに426 fitを計算し、同じ結果を確認。16 MiBでの旧fit再利用も確認しました。standardはtick区間だけを除いた元結果と一致します。bin度数・entriesは完全一致、ROOT統計量にはmerge順による最大約4.6e-15の相対丸め差があります。数値表の許容差は絶対2e-9＋相対2e-5です。切断したfitsの100件の再利用・残り326件の未再利用、不一致拒否、CPU超過拒否も検証しています。

各検証のログ・数値比較はtimeout_audit_100GeV、実入力結果はtimeout_compact_100GeV、timeout_reuse_100GeV、timeout_standard_100GeV、timeout_parallel_100GeVに保存。最新コードでの再利用はtimeout_resume_parallel_100GeV（新fit0・再利用389・標準318 PNG）まで確認し、plots/checkedは元ROOT/TSVのchecksumを維持しています。全量の再投入は自動では行っていません。
