# Strip内通過位置によるmuon ADC応答

同じchannel内の長手方向u・幅方向vのADC応答をdata/MCで比較するC++17/ROOT解析。共通geometry・温度関数・定数loaderとLandau–Gaussian密度を読み取りで再利用する。正式MIP定数生成は `analysis/mip_calibration/` の担当。本解析は定数の更新、本番digitizerへの適用、electron解析を行わない。

## 実行

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/channel_response/strip_position_response
bash run.sh --dry-run
# 既定は100 GeV、run44、各入力から1000 event。全量を既定にしない。
bash run.sh --label pilot
# 本依頼の少量確認例
bash run.sh --label smoke_100GeV --runs 40,44 --max-events 1000 \
  --channels 10011,20013,400000,420013 --plot-channels 10011,20013,400000,420013
bash submit.sh --label scheduled_pilot --runs 44 --max-events 1000 --dry-run
# 保存分布だけからPNG再描画。既存PNGは保存して新しいplot-labelへ出力。
bash run.sh --label verified_100GeV --plots-only --plot-label next_redraw
# 自分の完了済み収集を別条件で解析。入力ROOTは再走査しない。
bash run.sh --label coarse_bins --reuse-collection \
  ../../result/channel_response/strip_position_response/verified_100GeV --u-bins 1 --v-bins 1
# 保存histogramからfit条件だけを変える場合（元入力・Crossingsの再集計なし）
bash run.sh --label diagnostic_refit --refit-from \
  ../../result/channel_response/strip_position_response/verified_100GeV --min-fit-hits 150
```

`--energy 5`、`--energy 20`なども指定できるが、実在するmuon入力だけを選ぶ。electronと同じenergy/runが存在するとは仮定せず、なければ出力作成前に停止する。明示入力は `--input-manifest FILE`、6列のTSV: `sample(0=data,1=mc), run, decode, calib, -, truth`。`--sample`、`--runs`もmanifestに適用する。同一sample/run・同じcanonical pathの重複を拒否。既定のdata/MC組以外も明示可能。

`--max-events N`は各ファイルから等間隔のN event、0は全event。`--channels`は対象channelだけを制限し、trackは全unmasked channelで作る。`--mask`はdata/MC共通、trackにも対象にも適用する。全オプションは`--help`。既存labelの上書きを拒否し、結果の移動・削除は行わない。`--collect-only`、`--no-plots`も使用可能。

Slurmは現行`all` partition、既存`../execute.sl`を再利用する。既定2 CPU / 8 GB / 4時間、`--workers / --mem / --time / --partition`で変更する。dry-runは入力存在確認と完全なsbatchコマンド表示のみで、ビルド・保存・投入をしない。ROOT内容監査は実行時に行う。

## 入力とproductionの証拠

- data: `ECAL_data/analysed/2023/sps/{decode,calib}/mu-/<energy>GeV/`。
- MC: `CEPCScECAL_SML_Portable_update_new/Result_MC/{decode,calib}/mu-/mu_track/<mc-tag>/<energy>GeV/`。
- muon truth: 同じResult_MCの`generate/mu-/mu_track/<energy>GeV/`。既存`muon_mip_response`のmanifest規則を使用し、electron MCは流用しない。
- electron候補`Result_MC/calib/e-/e_ssa/threshold/`は本解析では未使用。electronの低LG tail除外・HG–LG intercept対策を取り消さず、古い条件をbaselineとして再実装しない。muonにはelectron tail/energy/CoG cutを適用しない。Eventsに適用件数0を明示する。
- raw/calibrationは保存`Run_Num, Event_Time, TriggerID/Event_Num`で一意性を検査してjoin。rawのCycleIDとfull CellIDも保存する。calibration側にCycleIDがないため、より粗い共通identityが重複したraw eventは全出現を除外する。entry順を一致と仮定しない。calibrationとの不一致はraw ADC選別を変更しない。
- truthは必要なprimary branchだけを読み、最大100箇所のPDGがmuonであることと保存energyを監査する。truthに保存event IDがない場合、本解析はtruthとdecodeをevent joinしない。現在のMCDigi sourceにはentryをTriggerID/CycleID/Event_Timeに保存する実装があるが、歴史的入力とそのsourceの同一性を証明できないためである。
- `production_audit.tsv`はmanifest上の対応と生成履歴の証明を区別する。CT補正や生成時MIP設定はパス名・mtimeから断定しない。現在のproducer source、定数snapshot/checksum、入力statを保存するが、それだけで過去の生成設定を確定しない。CT stateが`unverified_for_historical_input`の比較は、その制約を伴う。

## Crossings、track、ADC

既存muon Crossingsの位置・長さ・ADC定義と`TrackGeometry.hh`を再利用する。専用schemaでは、hitがない通過にもdecode温度センサーを使う。旧flatは校正後hitの温度のみを保持しており、missing hitの温度matchingと独立したADC選別を再現できないため、無条件には再利用しない。

- layer全体を除外した非加重track。奇数layerはX、偶数layerはYのfine座標。各投影7layer以上、投影角0.05 rad以内、RMS 3 mm以内が既定。対象layerのhit・ADCはtrack推定に使わない。重複fine座標をまとめてgeometry残差から各layerの点を選ぶ。
- geometryは現行data readerの`EBUdecode.cxx`を直接リンクする。strip中心、pitch45.4/5.3 mm、active45×5×2 mm。偶数layerではu=X−Xcenter、v=Y−Ycenter、奇数layerではu=Y−Ycenter、v=X−Xcenter。uの正側をSiPM側とは仮定しない。
- 全厚を通過する直線を要求し、長手1 mm、幅0.25 mmのmarginを設ける。長さは`2*sqrt(1+sx²+sy²)` mm。角度・散乱・位置分解能systematicはfit誤差に含まれない。外部trackerはなく、shower CoGも独立な入射測定として使わない。
- track作成はRaw HitTag=1、接続6300channelを対象としADC重み・energy cutなし。hitがなくても予測stripは分母に残す。observed/predictedは他layerのtrack条件下での取得割合であり、hardware単独efficiencyとは異なる。
- 同じ物理stripの複数memory-cell hitはまとめてmultiplicityを記録し、応答fitから除外する。分母からは除外しない。Eventsのrecorded_raw、HitTag1件数、物理strip数、校正record数を区別。SSAは本解析に不要で読まず、SSA仮想hit数は−1/not_read。0の測定値にはしない。
- 実測HG<2600、HG/LG finite、pedestalあり、温度正常、multiplicity=1を共通hit集合とする。HG=4000、HG>=4095、HG>=2600は別計数。4000だけでdataの強制置換だったとは断定しない。LGでHGを置換しない。
- `A=HG-float(pedestal)`、`A_path=A*2/length`、`A_current=A_path*C`の3段階を同一hit集合で保存。負/ゼロのpedestal差引きADCも保存し、histogramのflow、件数・平均にも残す。保存energyの負/ゼロ値でhit選別しない。
- 温度はrawのlayer別16sensorと現行`tempReconstruction`でstrip中心に補間。欠損・非finite・範囲外sensor（−20〜60℃外）は状態を記録し、3段階すべての応答fitから除外する。予測分母は維持。`C=1−(T−20)*alpha`、10µmは1.6/135、15µmは3.5/230。これは現行補正式の再現で、物理的20℃換算の独立な確証ではない。校正保存温度もCrossingsへ診断用に残す。

## bin、fit、基準位置

既定`--u-bins 3 --v-bins 1`は長手15 mm・幅5 mm。幅方向は分解能の未確立を考慮して保守的に1binとする。幅の細分化は例として`--v-bins 3 --position-resolution 0.75`で1.667 mmにできるが、これは分解能仮定を明示した診断である。対称な奇数binを使用し、中央の独立bin（u±7.5 mm、既定のv±2.5 mm; v=3binでは±0.8333 mm）が基準。外側binの一部はmarginで取得範囲から外れるので予測数も参照する。`--position-resolution 1.5` mm（既定）は細かすぎるbinを拒否するための**仮定の下限**であり、測定済み分解能ではない。bin幅は下限の2倍以上。予測track RMSをCrossingsに保持し、分解能・bin migrationの精密評価は別途必要。小規模検証で分解能を実証したとは解釈しない。

全sampleのADC binは共通（−200〜3000、8 ADC/bin）。histogramは絶対hit数、平均ADC・分母・flowを保持し、shape規格化はしない。Landau密度とGaussianの数値畳み込みは既存MIP fitの実装を再利用し、ROOT Landau mode displacementを補正したLandau MPVを求める。幅・sigmaも自由にfitして誤差を保存する。畳み込み後の最大位置とは違う。

threshold定数`expErfThre`の座標はHG−pedestal、補正前。各binで補正係数fの実測最大/最小を保存する。data/MC共通fit範囲の下端を`max[(threshold+3sigma)*fmax]`、上端を`min[(2600−pedestal)*fmin,2400]`としてwhole ADC binへ揃える。低ADC turn-onを一つの平均thresholdに置き換えず、保守的なplateauだけをfitする。4sigma下端のfitも行い、5%以内のMPV安定性を要求する。thresholdがchip補完・欠損なら未判定。

既定は区間内100hit以上、Minuit成功、共分散品質3、正NDF、MPV相対誤差15%以内、pooled Pearson χ²/NDF<=3、境界に張り付かないこと、MPV−2errorが下端+1binより上、4sigma fitも同じ品質を満たすことを要求する。低切断点に隠れたピークを採用しない。低統計・fit失敗・共通support不足は状態とNaNを保存し、補間しない。外部threshold自体のモデル誤差はこれだけで除去できない。

絶対MPVと`MPV(bin)/MPV(central)`をともに保存。異なるbinはdisjointなのでfit統計誤差を独立として伝播し、central/centralは1±0（定義）とする。共通threshold・温度・trackのsystematicと、異なる相対点が共有する基準MPVによる相関はこの1点ごとのerrorに含めない。比をさらにfitする際は共通基準の相関に注意する。

横error barは位置binの幅であって測定位置の誤差ではない。fit曲線はhistogramと異なる青/橙の破線で、MPV・幅・sigmaと状態をlegendに記載する。

## 交絡調整と独立検証

学習/検証は保存event identityをhashして二分し、同一eventの全channelを同じ側に置く。`--split-mode run --validation-runs 44`なら指定runを丸ごと検証へ回す。入力順・ADCは分割に使わない。

各axis別、channel別、split別に、run・温度bin（20℃を中心とする2℃）・sx/sy bin（中心0、0.02）・もう一方の位置binで層別する。

- `unadjusted`: 全予測通過。
- `position_matched`: sample内部で全位置binに共通するstratumを取り、各位置から最小予測件数だけをidentity hash順に選ぶ。
- `data_mc_matched`: 上記をdata/MC両方・全位置binに要求して最小件数を選ぶ。

ADC、target hitの有無を選択に使わず、分母を含めて条件を揃える。反面、どこかのbinが空ならそのstratumは使えない。data温度とMC20℃が重ならない場合は共通supportなしとなる。取得割合のSEは単純binomial近似で、少数例や0/1端点での厳密な信頼区間ではない。bin内の残存交絡はあり、共通supportを温度の外挿で作らない。`matching.tsv`に元件数とquota、`bins.tsv`に調整後の予測・取得件数を保存する。数は整数で、hit重みをfitへ持ち込まない。

SiPM10µm（layer4–27）・15µm（その他）の集計では、等channel加重とfit hit数加重を分ける。群平均はchannel MPVの平均で、混合分布をfitしたMPVではない。channel数、hit数、fit誤差、channel間SDと平均のSEを記録する。data/MC群比は共通の有効channel集合で別途出力する。

候補mapは学習側の`A_current`の相対MPVの逆数だけから作成する。既定は`position_matched`のu方向、`--map-axis`と`--map-condition`で明示変更可能。u/vの積や任意外挿は作らない。区間内のpiecewise constant候補として、測定できたbinかつ学習側に存在したrun・温度・角度・他方位置stratum内に限り独立検証eventのADCへ適用し、同じcovered hit集合で補正前後を比較する。未測定binはfactor=1で埋めず、未適用件数を出す。基準が測定できなければmap全体が未判定。区間内の連続変動、bin migration、温度外挿の妥当性は保証しない。map統計誤差はTSVに保存し、検証平均のSEは固定したmapに条件付けた統計誤差である。

**muonで得た応答mapをelectronに適用するには、digitizationのthreshold・SiPM応答より前の取り扱いを検討する必要がある。** 保存されたelectron ADCへ単純に掛けると未取得hitを回復できない。本解析は候補mapを本番digitizerへ適用しない。

## 出力と速度

基点は`analysis/result/channel_response/strip_position_response/<label>/`。Slurmログは`analysis/log/channel_response/strip_position_response/`。CoGを使わないので形式的なCoG階層を置かない。

| ファイル | 内容 |
|---|---|
| crossings.root | 保存identity/full CellID・対象layer除外track・予測/未取得・温度・3段階ADC、Events |
| histograms.root / bins.tsv | 位置別絶対ADC、flow、分母・平均・取得割合、温度/run/角度分布 |
| fits.root / fits.tsv / fit_status.tsv | 主/4sigma fit曲線、パラメータ・全誤差・状態、絶対/相対MPV |
| comparison.tsv | data−MC、data/MC、誤差、絶対平均差、両sampleの分母 |
| groups.tsv / common_channel_ratios.tsv | 10/15µm、等channel/ hit加重、channel間SD、共通channel比 |
| matching.tsv / channel_coverage.tsv | stratum別取得前件数とquota、未測定channel |
| candidate_map.tsv | 学習側のみの候補係数、誤差、未測定・非外挿状態 |
| validation.root / validation.tsv | 独立検証側の同一covered hitの補正前後分布・平均差・誤差・coverage。validation_position.tsvは位置別平均 |
| parameters.tsv / production_audit.tsv / event_audit.tsv / event_summary.tsv | 定数・geometry・由来の確認範囲、IDjoin、母集団の絶対hit数/energy |
| inputs/ / run.log / fit_worker*.log | manifest・snapshot・hash・保存source・入力stat・実行記録 |

PNGは`figures/<plot-label>/<data|mc|comparison>/<condition>/<energy>GeV/layerN/chipN/`。群図・layer mapはそれぞれcomparisonの下。図はPNGのみで、PDFなし。

- `*_adc.png`: 位置bin×ADC段階の絶対分布とfit。train/validationを別図にする。
- `*_response.png`: MPV、相対MPV、Landau幅、Gaussian sigma、取得割合、予測通過数。pedestal/path/temperatureの3段階を別図で表示し、全段階をfits.tsvにも保存。
- `*_ratios.png`: 絶対/相対MPVのdata/MC、平均ADC差、絶対hit数差。
- `*_confounding.png`: u/vと温度/run/sx/syの予測通過数。対応data/MCは同じ色軸。
- `*_before_after.png` / `*_validation_position.png`: 独立検証側の補正前後分布と位置別平均ADC。未測定mapならその状態を表示。
- `group*_response.png`: 等channel加重とhit加重の絶対・相対群平均、channel間SD。SDの図の誤差は正規近似SD/sqrt(2(N−1))で、fit誤差とは異なる。
- `candidate_maps/.../layerN.png`: 実strip座標上の候補係数。同じlayerでx/y同縮尺、data/MC共通色軸。白抜きは未測定。

PNGは既定で予測件数上位12channel、数値とhistogramは全収集channel。`--plot-limit 0`で全channel、`--plot-channels ID,...`で明示指定。図選択は数値集計や群平均を変えない。

raw/calibrationはscalar IDの事前走査後、必要vector branchだけをevent単位で読み、channelごとの再走査はしない。AggregateはCrossingsを一度だけ読み、予測レコードをメモリ内でchannel別に集計する。全量ではCrossings数に比例したメモリが必要なので、割当てに応じてrun/対象channelを制限する。fitは割当CPU以内の独立processに分け、別ROOTへ書いて成功後haddする。描画は1 processで数値ROOTをread-onlyで読む。plots-onlyは元入力・Crossings・定数を開かず、前後の数値checksum一致を検査する。

`--reuse-collection`はこの解析のCOLLECTED、schema、保存した収集設定と現行収集source hashを要求する。`--refit-from`は完了したhistogramと保存設定を引き継ぎ、位置bin・split・matching条件の変更を拒否する。再利用ファイルのhashを前後で照合する。新規完了結果にはhistogram/設定のanalysis_products.sha256も保存する。入力選択・mask・track条件を変えた再利用を拒否し、定数snapshotを再利用する。旧muon/mip Crossingsはschemaと選別が違うため再収集を選ぶ。共通ヘッダー・他解析コードへの書込みは不要。

## 小規模検証

2026-10-13、ROOT 6.36.08 / g++ -O2でビルドし、以下を実行した。**修正済みの数値は `verified_100GeV/`、現在の図はその `figures/geometry_checked/` を参照する。**

```bash
bash run.sh --label validation_100GeV --runs 40,44 --max-events 1000 --plot-limit 4
# 上記の実行時の位置設定はu=3, v=3, resolution-floor=0.75 mm。
# 現在の既定v=1でなく、この診断条件を新labelで再収集する場合:
bash run.sh --label position_diagnostic --runs 40,44 --max-events 1000 \
  --u-bins 3 --v-bins 3 --position-resolution 0.75 --plot-limit 2
bash run.sh --label verified_100GeV --refit-from \
  ../../result/channel_response/strip_position_response/validation_100GeV --plot-limit 2
bash run.sh --label verified_100GeV --plots-only --plot-label geometry_checked --plot-limit 2
bash submit.sh --label final_dry_run --runs 40,44 --max-events 1000 \
  --u-bins 3 --v-bins 3 --position-resolution 0.75 --dry-run
```

実行済みlabelは上書き不可。再実行では新しいlabel/plot-labelを指定する。

- data/MC × run40/44、それぞれ1000 event、合計4000 event。全4000件を保存IDで一意に照合し、欠落・並べ替え・重複はこの標本では0件。
- 6300channelを対象に収集し、どちらかのsampleで予測通過があったchannelは2694。予測通過data=20,942、MC=50,891、共通3段階ADCのusable hitはdata=14,966、MC=43,677。全71,833通過と58,643 hitについて集計ROOTとの一致、flowを含む件数、3段階の同一hit集合を検査した。
- 実入力fitの全581,904 status行は統計不足112,977 / 予測なし468,927。有効な位置別MPVは0件。候補map・実dataでの補正効果は未測定であり、0応答の測定とはしない。全量・全energyの解析結果ではない。
- 30layerについて対象layerのhitを大きく変更しても除外trackが不変であること、u/vの向き、通過長・温度補正を一時C++検証で確認した。
- 人工Landau–Gaussian（MPV350, Landau幅35, Gaussian sigma45 ADC、seed131026、30,000 hit）はMPV350.0099±0.6242 ADCを回収。Gaussian sigma44.6343 ADCの保存/読込とthreshold sigma8 ADCの独立性を検査した。隠れたMPVはMPV_hidden_or_boundary、3hitはinsufficient_statisticsで非採用。
- `synthetic_verified/`は人工位置応答data=(0.9,1.0,1.1)、MC=1の工程検証専用。学習側でdata相対MPV=(0.90484,1,1.09881)を回収。独立検証data 5701 covered hitの位置別平均ADCの最大差は147.349→9.990 ADC。**これは実験dataでの補正効果ではない。** map_verification.tsv、map_validation_check.tsvに数値と判定を保存。
- plots-onlyの数値ROOT/TSV checksum不変、refitでの元histogram checksum不変、u-bin変更を伴う不適合なrefitの拒否、存在しない5 GeV muon入力の拒否、Slurm完全コマンドのdry-runを確認。実ジョブは投入していない。

初期工程の`smoke_100GeV/`、`validation_100GeV/`、`coarse_validation/`、`synthetic_validation/`は由来・分布を保持している。初期Fits ROOTではthreshold sigmaとGaussian sigmaのbranch名が衝突していたため、そのfit表・幅の図は現在のbaselineにしない。`verified_100GeV`と`synthetic_verified`では `threshold_sigma` / `gaussian_sigma`を分離し、実行時に必要branchを検査する。元のCrossings・histogramとMPV計算の定義はこの修正で変更していない。

未検証: 実入力での十分な位置別MPV統計、独立なtrack位置分解能とbin migration、実runを丸ごとhold-outした補正効果、生成時CT/MIP設定の確証、electronへの効果、全muon run・全energy、実Slurm投入。現在のproducer source控えは過去productionの同一性を証明しない。

一時C++/ROOT検証ソースと実行コマンドは結果のinputs/verification_sources.tar.gzに再現資料として保存し、作業用.scratch内の検証コード・実行ファイルは完了後に削除する。恒久的なテストや共通処理の変更は追加していない。
