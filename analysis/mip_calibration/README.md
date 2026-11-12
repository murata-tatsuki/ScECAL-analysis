# Data muon MIP calibration

DecodeのHG ADCから、dataだけでchannelごとのLandau MPV候補と採否を求めます。C++17 / ROOT / shell。旧 `analysis/mip` と `mip_new` は変更しません。本番定数の置換・MC生成・electron再解析は行いません。

## 実行

```bash
cd /home/murata_t/ScECAL_BeamTest/analysis/mip_calibration
bash run.sh --dry-run
bash run.sh --label all_muons
bash submit.sh --label all_muons_slurm
# 小規模収集。Nは各入力ファイルから等間隔に選ぶevent数。
bash run.sh --label pilot --runs 40,44 --max-events 3000 --collect-only
# 完了したhistogramから別labelへ再fit。decodeを読み直さない。
bash run.sh --label pilot_fit --refit-from ../result/mip_calibration/pilot
# channel制限。track自体は他の全channelから構成する。
bash run.sh --label selected --runs 40 --max-events 3000 --channels 20013,420013
# 同じlabelの保存済みfitからPNGだけ再生成。
bash run.sh --label pilot_fit --plots-only --channels 20013,420013
# 既存channel図を維持し、6x6のchip一覧PNGだけを追加・再描画。
bash run.sh --label all_muons --chip-plots-only
# 温度補正なしの診断を別labelに保存。
bash run.sh --label pilot_none --refit-from ../result/mip_calibration/pilot --temperature-mode none
```

既定は全run・全event・全6300接続channel、4 fit process。Slurmは4 CPU/16 GB/24時間、`--partition / --mem / --time / --workers`で変更します。投入時のdry-runは設定とパスを検査し、ROOTの内容監査は実行時に行います。ジョブを実行中のlabelへの別投入は出力directory/lockで拒否します。既存labelの解析数値は上書きしません。失敗結果を消さず、別labelを使用してください。

`--input-manifest FILE`は`run<TAB>絶対decodeパス`、任意の`run/path`ヘッダーまたは`#`コメント。`--runs`も適用します。同一canonical pathの重複、run番号不一致、壊れた/recovered ROOT、必要branchの欠損を拒否します。入力が更新されたら収集完了にしません。`--channels`は物理CellID（layer*100000+chip*10000+channel）で指定します。制限した収集から未収集channelを再fitしても測定値は生成されません。

`--refit-from`では元の収集設定と定数snapshotを引き継ぎます。track、入力、event、factor-step等を変えるには再収集が必要です。fit条件・温度mode・channel・worker数は変更できます。MPV相対誤差上限は再fitでも新しい既定値10%を使います。旧5%を再指定する場合は`--max-relative-error 0.05`を明示してください。`--plots-only`は数値ROOTをread-onlyで読み、前後のchecksumを照合します。`--no-plots`は数値のみ、`--collect-only`は集計まで。全optionは`--help`を参照。

## 入力監査と定義

- data decode: `/home/murata_t/data_beamtest/ECAL_data/analysed/2023/sps/decode/mu-/100GeV/`
- pedestal、旧MIP、HLは現行 `Analysis_edit/share/`。thresholdは `analysis/result/threshold/threshold.root` の `expErfThre`。旧MIP/HLは比較および既存定数loaderのため読みますが、track、ADC、fit初期値・制約には使いません。
- 有効channelはlayer 0–29、chip 0–4は36、chip 5は30、計6300。旧ROOTの6720行は32*210。layer30/31は測定対象に加えず、6300行の中でも未測定channelはゼロの不採用印だけを出します。
- `Run_Num,Event_Time,CycleID,TriggerID`を全入力のscalar passで確認し、重複identityは全出現を除外。vector branchは一度だけ集計。decodeだけを使うのでcalib/SSAとのentry順joinはありません。full CellIDのmemory bitsはCrossingsに保持し、同一物理channelの複数HitTag=1は応答fitから除外します。
- HitTag=1のstrip中心による非加重のtrack。奇数layerでX、偶数でY。**対象layer全体を除去してから**推定。既定は各投影7layer以上、各投影角0.05 rad、RMS 3 mm以下。既存muon解析から `TrackGeometry.hh`へ抽出した同一処理を共有します。energy/MIP cutはありません。他layerの取得thresholdに条件付けられており、完全に無条件なmuon母集団ではありません。
- geometryは実際のdata readerの `EBUdecode.cxx`をリンク。pitch 45.4/5.3 mm、active strip 45*5*2 mm。DetectorConstructionのcrystalH=(2+2*ESR)とも照合。長手1 mm、幅0.25 mmの余白を含め全厚通過する直線を要求し、`length=2*sqrt(1+sx^2+sy^2)` mm。追跡誤差や散乱によるsystematicはfit誤差に含みません。
- `A=HG_Charge-float(pedestal)`、`A_path=A*2/length`。HG<2600の実測HGを使用。HG/LG切替、HG=4000の専用フラグ、12bit飽和HG/LG>=4095を区別し、LGでHGを置換しません。`forced_HG4000`という診断列は値が4000の件数です。実dataのDecode.cxxは12bit値を保存するため、その値だけでMCのような強制置換だったとは断定しません。GainTagはHL-modeの保存bitとして数え、HG/LGの両branchがあるためこれでgainを推測しません。crosstalkを取り除きません。
- 予測stripにhitがなくても分母に入り、`observed/predicted`および`unique/predicted`を保存。これは他layerのtrack条件下の取得割合です。trackの割当て誤差・隣接strip・重複を含み、hardware単独のefficiencyではありません。SEは単純binomial近似です。

温度はdecodeの16センサーから、現行data readerの関数を直接呼びstrip中心に逆距離補間します。現行 `C=1-(T-20)*alpha`、alphaは10 µmで1.6/135、15 µmで3.5/230。空センサーは現行reader同様20℃と記録しますが、温度欠損・異常のある通過はcurrent modeのfit入力から除外します。補正式再現のcurrent_adc診断histogramには空センサー時の20℃代入も残します。none modeは温度を必要としないため、その統計数は異なる場合があります。

校正既定は`A_current=A_path*C`のLandau MPVです。**現行readerの座標・補正式を再現する値であり、物理的な20℃への換算が独立に確立した値ではありません。** 補正なし`A_path`も同じhit集合で保存し、反対modeのfitを診断として出します。逆補正式へは切り替えません。`none` modeのMPVは測定温度混合の値で、物理的な20℃値ではありません。

## thresholdとfit

`threshold_finding.cc`は `HG_Charge-pedMeanHigh` をfitするのでthreshold/sigmaは**pedestal差引き後、温度・通過長補正前のADC**です。主fitの効率は `0.5*erfc(-(x/f-threshold)/(sqrt(2)*sigma))`、`f=(2/length)*C`（noneではC=1）。各runとf区間ごとに分布と独立規格化を持ち、共通のLandau MPV・Landau幅・Gaussian幅を同時fitします。run間で同じMPVという仮定はdependencies表で確認してください。

fは既定0.005間隔に丸め、ADCとthreshold/sigmaに**同一の丸めた係数**を適用します。各区間の実際のf最小最大もROOTに保存します。丸めによる相対誤差は`0.0025/f`以内。raw/path/currentの診断histogramとCrossingsには丸め前の値を保持。必要なら`--factor-step`を細かくして再収集します。別run・別係数のthresholdを一本の平均曲線に置き換えません。

主fitはpedestal差引き後・温度/通過長補正前のADCで40以上を対象とします（`--fit-min-adc`）。全条件でこの下限を満たす共通の補正後bin境界へ切り上げ、その区間内で条件ごとに規格化します。低ADCの別成分をLandau–Gaussianの裾へ無理に当てはめないためのfit範囲指定であり、元のhistogram・取得効率は変更しません。20/60 ADCへの変更も採用前に比較します。

fitはbin積分を3点Gauss求積したPoisson devianceの最小化。条件別normalizationは解析的にprofileし、幅・sigma・MPVは自由、MC/旧MIP/electron energyからの拘束はありません。Landau MPVはROOT Landau mode shift -0.22278298を補正したパラメータ。畳み込み後の最大位置は別の`convolution_peak`です。主fitとMPV profileは初期値をdata分布の最大binから3通り作り、収束した中で最小devianceの解を使います。条件変更fitは主fitの解から始め、採用候補の品質や安定性に問題が出た場合に3初期値で再試行します。Gaussian幅/Landau幅>8では積分変数をLandau側へ切り替え、極端に狭いLandau分布も分解できる区分Gauss–Legendre積分を使います。旧来の1600点打切りによる積分の凹凸を避けます。

旧FitMIP.Cの`fixed_width`はSetParameter/SetParLimitsのみでFixParameterではないため、そのAPIは流用しません。

χ²は条件ごとに隣接binを期待数5以上へpoolしたPearson **χ²そのもの**。NDFはpool数−条件別規格化数−3。fitのPoisson devianceとは別値です。低統計条件もfitに使いますが、NDFが不足すれば採用しません。

既定の直接採用条件:

- fit区間200hit以上、Minuit status=0、完全な共分散、正のNDF。
- MPV相対誤差≤10%、MPVと各幅の相関≤0.98、パラメータ境界に張り付かない。
- χ²/NDF≤2（10 µm）、≤1.7（15 µm）で現行readerと一致。
- `MPV-2*error`が全有効条件の`(threshold+4*sigma)*f+1bin`と実際のfit下限+1binの両方より上。ピークがturn-onに隠れたfitは収束しても未判定。
- threshold+2/3/4σ以上だけのfit、bin幅2倍、上限0.85倍、初期MPV0.7/1.3倍、threshold±sigma（負shift側はwidth1.5倍）、畳み込み求積2倍、低ADC下限±20 ADCが品質条件を満たし、MPV変動≤5%。温度mode比較はこの安定性判定に混ぜません。
- MPVを0.1 ADC、主fitの0.25/0.5/0.8/1.2/1.5/2倍に固定し、両幅を複数初期値から再fit。全点でΔdeviance≥3.84を要求。平坦なら`profile_not_identifiable`、主fitより低い別解があれば`profile_lower_minimum`で未判定。
- chip平均で補完したthresholdは未判定。current modeでは温度が正常な通過だけをfitする。
- 期待数0.001未満のbinにあるhitがmax(3,0.5%)以上ならthresholdモデル不適合として不採用。poolしたχ²だけでは見逃す極端な低確率領域を監査する。

`decision=1`はこの範囲での直接測定、0は未判定、2はfit失敗/χ²不適合。固定幅による救済はありません。主fitの下限到達は`MPV_at_lower_bound`、`landau_width_unresolved`、`gaussian_width_unresolved`に分けます。Gaussianを完全に外したLandau単独fit（sigma=0、自由パラメータ2個）はvariant15として比較保存しますが、その値を自動採用しません。独立pedestal幅は電子ノイズの一部を表すので、MIP分布のGaussian幅と同一視して固定しません。外部thresholdはelectron ADC spectrumのモデルfitであり独立なhardware turn-on測定ではなく、既存調査で高χ²が多数あることにも注意してください。本解析はthresholdモデルを真と証明するものではありません。採用をturn-onから離れた領域と診断fitの整合に制限します。run、温度、位置の応答変化が大きければ共通MPVの物理的解釈を再検討する必要があります。

## MPV・thresholdの詳しい診断

```bash
# 完了した全run集計から、MPV境界・Gaussian境界・良好なchannelを比較。
bash diagnose.sh ../result/mip_calibration/all_muons fit_diagnosis_custom 2440010,400020,1520021
# 修正後の手順で全channelを別labelへ再fitする例。
bash submit.sh --label all_muons_refit --refit-from ../result/mip_calibration/all_muons --max-relative-error 0.10 --fit-min-adc 40
```

`Diagnose.cc`は元の全ADC、20/40/60 ADC以上、threshold+2/3/4σ以上、threshold±σ・幅1.5倍、Landau単独、診断用に負MPVまで許したfitを比較します。負MPVは模型不適合・縮退の診断専用で、校正ROOTへ出しません。元の全ADC/40 ADC以上/plateau3についてMPV固定走査を行い、幅を毎回再fitします。比較区間が異なるfit間のdevianceの絶対値は直接比較しません。Δdevianceは同じデータ範囲・モデル内の比較です。

出力は`diagnostics.root`、`models.tsv`（`hypothesis`が比較条件名、`model`が数値model識別子）、`profiles.tsv`、`runs.tsv`、`noise.tsv`と`figures/cell<ID>_diagnosis.png`です。`runs.tsv`は40 ADC以上でthresholdを−2/−1/0/+1/+2σ動かしたrun別fit。run別だけは50 hit以上で診断し、少数runの値を校正に採用しません。`noise.tsv`は独立pedestal ROOTの`ChnLevel/PedHighSig`です。thresholdの変化はmuon形状変化と縮退し得るため、run別の最良shiftだけを新thresholdとして採用しません。

通常fitでもvariant12に低ADCを含む元モデル、13/14に下限−20/+20 ADC、15にGaussianなしを保存します。12/15はモデル診断で、13/14は採用の安定性検査です。PNGの赤線、6×6 chip一覧、互換ROOTの10 branchは維持しています。

## 出力

結果は `analysis/result/mip_calibration/<label>/`、ジョブログ・投入記録は `analysis/log/mip_calibration/`。channel図は `figures/layer<l>/chip<c>/cell<ID>_fit.png`、2250*1500、PNGのみ。

| 出力 | 内容 |
|---|---|
| mip.root / MIP_Fit | 現行互換10branch。Int_t: CellID,NDF,chnEntries。Double_t: LandauMPV,LandauWidth,GauSigma,ChiSquareと3誤差。ADC単位。 |
| histograms.root | 生HG、pedestal/path/current ADC、run/f区間のfit入力、Channels、Conditions、補正式 |
| crossings.root | 対象hitのない通過を含むidentity・track・長さ・温度・ADC。新規収集時のみ |
| fits.root | Diagnostics、FitVariants、fitに使ったhistogram、保存TF1、混合効率を含む期待bin数のTGraph、profile |
| fits.tsv / fit_variants.tsv | 候補MPV・幅・誤差、fit範囲、bin幅、共分散品質、相関、χ²、安定性、採否理由。raw_fit_lowは補正前ADC下限。boundary_maskのbit0/1/2はLandau幅/MPV/Gaussian幅の下限、bit3/4/5は上限。model=0はLandau–Gaussian、1はLandau単独。 |
| channels.tsv | 条件付き取得効率、除外内訳、旧個別/平均補完と新採否・実効MIPの比較 |
| summary.tsv | 10/15 µm、layer/chip、旧個別/補完別の集計。-1は全体 |
| dependencies.tsv | run・温度1℃・長手5 mm・幅1 mm・通過長0.001 mm別の件数と補正前後の平均ADC（MPVではない） |
| reader_status.tsv / schema.tsv | 現行readerの採用/補完数、補完平均、読込み可能性、branch型 |
| event_audit.tsv | manifest/run照合、処理event、identity重複、通過予測件数 |
| inputs/ / run.log | manifest、設定、コマンド、定数snapshot、実行binary・source archive、現行reader source、checksum、入力stat |

不採用の互換行は`LandauMPV=0,NDF=0`。候補値はDiagnostics/TSVに残し、平均補完値を実測値として書きません。現在のCalibration.cxx/Extract.cxxは両SiPM群に1個以上の有効donorがなければthrowします。小規模結果でdonorが足りない場合、ROOT schemaが正しくてもそのまま校正には入力できません。`reader_loadable=0`とNaN補完値を出します。両群にdonorがある場合は生成したROOTを共通reader helperで再読込みし、各channel採否・実効値を照合します。読込み不可でもdonorのある群には群内補完した場合の参考実効値を表記しますが、ファイル全体を読めることは意味しません。採用集合が片寄れば平均補完も片寄るため、`reader_loadable=1`は物理的な校正妥当性の保証ではありません。

fitパネルはdata分布を濃い青、fit曲線を赤で表示します。個別channel図の右上・左下とchip一覧に共通です。個別図の左上は温度補正前後のdata分布の比較で、青系2色を使います。

各layer直下の `figures/layer<l>/layer<l>_chip<c>.png` は、個別channel図の右上（主fit）を `Divide(6,6)` で並べた一覧です。channel 0–35を左上から横方向に並べ、採否・統計数に関係なく全接続channelを表示します。chip 5のchannel 30–35は未接続と表示します。fit曲線が保存されていないchannelも、同じ温度modeの保存済みrun/factor分布をfit範囲・bin幅に合わせて集計し、ヒストグラムと理由を表示します。各padの描画面積は既存2x2図の1padと同じで、canvasは縦横3倍（現在のROOTでは6738×4416 pixel）です。

通常の描画・`--plots-only`ではchip一覧も生成します。`--chip-plots-only`は一覧だけを生成し、既存の個別channel図・全体図を維持します。`--channels`を指定した場合、指定channelを含むchipの全接続channelを一覧に表示します。いずれも保存済み数値とfit曲線を読み、再fitしません。

個別の空channelのPNGは既定で省略しますが、ROOT/TSVには全6300channelを保持し、明示channel指定時には未判定図も出ます。再描画はTF1の再fitをせず保存期待値graphを使用します。workerは独立ROOTに書き、成功後にhadd。`.scratch/mip_calibration_runs/`の今回のshardだけを統合成功後に削除し、失敗時は再開調査のため場所をrun.logに残します。

## 参照と限界

監査したreader: `/home/murata_t/data_beamtest/analysis/ECAL_Analysis_LCIO/src/Calibration.cxx`の`MIPFactorExtract_simulation`と温度/HG処理、`Analysis_edit/src/Extract.cxx`の`MIPPeakValueExtract`、`MCDigi.cxx`のCT前一次光量。digitizerでのCT補正と実測ADCのcrosstalk込みMPVは別です。

比較・再利用: `analysis/mip/{MIP_calibration.cc,FitMIP.C}`、`analysis/mip_new/MIP_calibration.cc`（差は図の保存先）、`analysis/channel_response/muon_mip_response/`、`investigations/e_ssa_threshold_review_20260930/CurrentParameters.hh`、`docs/research/e-ssa-threshold.md`。旧MIPはSSA Event_Numとdecode TriggerID、SSA-tag付きCellIDを照合し、温度補正はコメントアウトされていました。今回はその選別を引き継ぎません。保存した現行producer sourceは過去入力の生成時sourceを証明するものではありません。

## 実装時の小規模検証（2026-10-13）

`analysis/result/mip_calibration/validation/`に実dataの検証結果を保存しました。run40/44の合計68,362 eventから10代表channelだけを集計し、通過予測6,602件。15 µm/10 µm、現行個別/平均補完、低統計/未取得、threshold近傍を含みます。CellID 10011は357.891±10.757 ADCで直接採用、他9channelは未判定または不採用でした。10 µm donorがないので、この小規模mip.rootの`reader_loadable`は0です。全channel校正の結果ではありません。

- `validation_checks.tsv`: 既知MPV人工分布の回収、threshold下に隠れたMPVと低統計の非採用、保存TF1の評価、30layerの対象除外、通過長、37条件のthreshold座標系、実decode 392件のidentity/full CellID/ADC/温度/長さ照合。
- `synthetic_validation.tsv`: ROOT TRandom3 seed=131026。既知MPV 500/140 ADC、Landau幅55/16、Gaussian幅65/22、2条件f=0.90/1.05の検証。十分な統計のMPVは500.480、500.302、139.853 ADC。hidden-peakはthreshold180、MPV140のケース。人工入力・専用検証コードは一時領域で処理し、完了後に削除。
- `synthetic_reader_check.tsv`: **実際のCalibration.cxxとExtract.cxx本体**をリンクして読込み、6300接続channelの実効定数差0。不採用行が平均補完されることを確認。data readerはlayer30/31も内部で補完するため内部補完数は420多くなりますが、出力にその測定行は追加しません。
- `production_reader_check.tsv`: 実data小規模ROOTでは、両readerとも10 µm donor不足で意図どおりthrow。架空donorや平均値の実測扱いで読込みを通していません。
- `inputs/plot_checksums_{before,after}.txt`: PNGのみ再描画でhistograms.root/fits.root/mip.rootが不変。

`validation_collect/`は代表選定に使った各run3000 event・全channelの保存分布。`validation_data/`は温度欠損時の20℃代入を含めた初期診断で、最終的な「温度不明の通過をcurrent fitから除外する」検証は`validation/`です。`validation_no_temperature/`は最終分布の2channelを再fitした補正なし比較で、decode再読込みを行いません。これらの科学的診断分布・条件・由来は保持します。本番校正定数の置換、全electron energyの再解析、全muon runの本実行はしていません。

## fit方法の再検討と代表channel検証（2026-10-13）

元の全run結果`all_muons/`を読み取り、9代表channelで低ADC除外、MPV固定走査、threshold±σ/幅変更、run別threshold走査を比較しました。初回の`fit_diagnosis/`と積分方式を変更した`fit_diagnosis_stable_integral/`は比較記録として保持しています。この初期2診断のLandau単独指定はROOTが以前のsigma下限0.3を保持したため、厳密にはsigma=0.3の近似です。最終コードでは固定変数をlimit設定前に定義し、`fit_method_validation/`のvariant15と詳しい診断では厳密なsigma=0を確認しています。

- 低ADCにある小成分は、主MIP成分のLandau–Gaussian×効率では説明できない例がありました。代表のpedestal sigmaは約3–5 ADCですが、20–40 ADCにも少数hitがあります。成分の発生原因は確定していません。40 ADC以上を主fitにし、20/60 ADCとplateau fitの整合を求める方法を採用しました。
- CellID **350001**: 低ADC除外後はMPV **545.30±8.23 ADC**。20/40/60 ADC・plateauで安定し、新たに採用。Gaussianを外すと同一区間のdevianceが約52.4増えるので、幅を0にする方法は適切ではありません。
- CellID **2440010**: 40 ADC以上でも正のMPVに安定した最小値がなく、負のMPVまで許すと解が負側へ動きます。threshold変更でも大きく変動し、未採用を維持します。
- CellID **400020**: GaussianなしとGaussian幅下限のfitのdeviance差は約0。Gaussian幅を決められていません。threshold変更でMPVが約59–93 ADCへ動くため、幅固定による救済は行いません。
- CellID **1520021**: 良好な10 µmの対照。MPV **135.27±2.33 ADC**で採用を維持します。
- run別では最良threshold shiftの向きが揃わず、低統計やMIP形状との縮退もあります。runごとのhardware threshold変化を確定できる証拠とは扱わず、既存thresholdを任意に置換していません。

最終コードの`fit_method_validation/`は、診断例に加えて固定seed=131026で元の誤差超過・threshold不適合・境界・採用済みから選んだ**49 channel**を検証しています。採用は**8→12**（15 µm **5→9**、10 µm **3→3**）。元の採用8 channelは全て維持。新規採用は140029、210024、250026、350001で、今回の増加は低ADCモデル不適合の解消によります。この標本では10%への緩和だけによる増加はありません。全6300 channelでの増加数はまだ再実行していないため未確定です。

`method_comparison.tsv`が元結果との比較、`implementation_comparison.tsv`が初回検証と最終コードの数値比較です。最終化で主MPVと採用集合は不変、Gaussianなしの48 fitはsigma=0になりました。`inputs/fit_trial_sources.tar.gz`に初回のコード、`inputs/analysis_sources.tar.gz`に最終コードを保存しています。最終検証用ROOTは6300行・互換10 branch、採用12/非採用6288のreader伝達を確認。49 channelだけをfitした検証ROOTであり、全channel校正結果ではありません。

人工分布の`synthetic_validation.tsv`はTRandom3 seed=131026、2条件f=0.9/1.05、各30000発生の5ケースです。既知MPV140/450 ADC、低ADC混入、threshold下のMPV、Gaussianなしを検証。`quadrature_validation.tsv`の175条件では独立な適応積分との最大絶対差5e-12以下でした。確認用一時コードは完了後に削除し、数値結果を保持します。元の`all_muons/`の数値ROOT、収集条件、本番校正定数は変更していません。
