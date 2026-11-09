# ScECAL 2023 PS・SPS 博士論文用引き継ぎ

## 0. 対象・調査方法・読み方

作成開始は **2026-10-24 09:16:34 Asia/Tokyo**。NFS待機中のagent再起動により中断し、同日18時台に未完成を確認して採取・組立・検証を再開した。採取時刻は各ファイルの記録を優先する。2026-10-26に再点検し、保存済みの採取・ROOT読出し記録から組立と最終検証を完了した。本書の「current／現行」は各sourceの採取時点を指し、10月26日の最新sourceを再取得した意味ではない。対象は2023年PS/SPSのScECALである。2022年SPS・2025年KEKの測定結果、SiECALを含むGNN/PFA研究は対象外。後年の装置紹介や旧年の定数が関係する場合は「依存・参考」として区別した。著者の原稿そのものは確認していない。

元のコード・定数・データ・結果は変更していない。調査は限定したdirectoryの一覧、コード・config・保存文書・表・ログの読出し、保存ROOTのobject/branch/fitの抽出、既存図の無改変コピーである。simulation、digitization、全event処理、再fit、性能の再集計を実行していない。スクリプトに含まれるjob投入・削除命令も実行していない。抽出スクリプトは `scripts/`、採取記録は `audit/` にある。保存ROOT読出しは新しい測定ではない。

- `evidence_status=確認済み` は記載・実装・保存値を直接読んだ意味であり、物理的正しさや再現性の認証ではない。`推定` は状況証拠による対応、`不明・未検証` は根拠不足。
- `status=暫定` は論文への採用判断待ち。新しいmtimeや良好なfitを理由に最終採用とはしない。
- `results.csv` の空欄は未保存・未確認でありゼロではない。`valid=unknown` は品質未確定、`false` は明示された除外・不正fit、`true` は限定された保存判定を満たす意味である。
- current sourceと生成時sourceは別物である。生成時のsourceが不明な結果にcurrent hashを割り当てていない。採取は逐次でatomic snapshotではない。ファイルごとのmtime・採取時刻・読出し中の変化をmanifestに保存する。大容量ROOT/rawの全体hashは取得していない。

先に本書を読み、値は `results.csv`、試料と用途は `runs.csv`、図は `figures.csv`、条件の参照は `manifest.yaml` を使う。以下の `S_CAL` 等は [SOURCES.md](SOURCES.md) に絶対path、同梱実体への相対リンク、SHA-256を示すsource aliasである。本文の行番号は今回採取したcurrent sourceに対するもの。過去の生成版を示すものではない。

主要な留保は、実機32層の説明と解析・MCの30層の違い、実際の定数tableの生成来歴不足、温度係数と基準の未確定、旧版と修正版MCの混在、低エネルギーfitの不安定性、独立な検証sample・系統誤差の不足である。

## A. 所在と処理経路

|置き場|実体と役割|
|---|---|
|`/home/murata_t/ScECAL_BeamTest`|`/net/megsrv01/export/megraid01/users/murata_t/scecal/ScECAL_BeamTest`。個別較正、診断、resolution解析、研究記録、図表。大文字表記は実際のdirectory名に従う。|
|`/home/murata_t/data_beamtest`|`/net/megsrv01/export/megraid01/users/data_beamtest`。raw/decode/calib、基本解析、Geant4/digitization。root全体は単一Git repositoryではなく、analysis・simulationに独立したGitがある。|

`/megraid01/users/...` は同じNFS実体への別名である。Git branch/commit/状態取得の成否はmanifest.repositories参照。statusのtimeoutをcleanと解釈してはいけない。主要ファイルのhashをGit情報に併記する。

実データの確認できた経路は次のとおり。

```text
ECAL_data/data/2023/{ps,sps}/{particle}/{beam label}/*.dat
  → analysis/ECAL_Analysis_LCIO/run/decode.sh → decode.cxx → src/Decode.cxx
  → ECAL_data/analysed/2023/{ps,sps}/decode/... : Raw_Hit
  → run/calibration*.sh → src/Calibration.cxx + run/global_config の定数
  → .../calib/... : Calib_Hit
  → ScECAL_BeamTest/analysis/resolution/SingleEnergyAnalysis.cc
  → analysis/result/resolution/{data,simulation/threshold}/{CoG}mm/*.root
  → MultiEnergyAnalysis.cc → comparsion/comparison_5mm.root、figures/comparsion/dataThre
```

simulationは `simulation/CEPCScECAL_SML_Portable_update_new/CEPCScECAL` のgeometry・generator・stepping/event集約から `generate/` にtruthを保存し、`Analysis_edit/main.cxx` と `src/MCDigi.cxx` がRaw_Hitを作る。以後Calibrationと評価コードを共用する。共用コードでも、入力定数、温度、threshold、tail除外、beam mapが同一だった証拠にはならない。

2026-10-12のresolution比較は `S_RESLOG` の開始03:58:52、終了06:25:12、exit 0と入力pathから裏付けられる。一方 `run/global_config` の現在値や `run_simulation/global_config` の最後の代入は、過去の全calibを生成した証拠ではない。旧 `calib_simulation`、`default`、`custom/hl_intercept`、`e-_ssa`、`threshold`、`threshold_MPVwoCT` を同じ版として扱わない。

rawから全段階を一対一で結ぶ完全なproduction manifestは見つかっていない。`ECAL_data/analysed/2023/sps/summary.txt` は「latest」や重複するcalib2ラベルを含み、現在のconfigと一致しない。これは古いメモとして併記する。既存研究索引 `docs/research/{README,reassessment,raw-calibration,muon-summary,hg-lg-tail,e-ssa-threshold}.md` は探索の入口であり、個々の保存物より優先する採用記録ではない。

探索は上記2つのtree、関連する `investigations/result`、run/config/ソース、較正・resolution・channel_responseの保存summary、2023 rawの粒子directoryを中心に行った。rawの内容を全件読まず、directory探索は深さを限定した。別diskのbackup・DAQ記録は限定探索であり完全な電子logbookを確認できていない。「見つからない」は不存在を意味しない。

## B. 実機・campaign・run

### B.1 装置の版を分ける

2020年設計論文（文献L1、PDF pp.3–4、印刷pp.2–3）は30 EBU、45×5×2 mm³ strip、210 strip/EBU、6個の36-channel SPIROC2E、16 analogue memoryを説明する。これは最初の2層の開発段階の資料であり、2023年as-builtの証明ではない。

2024/2025年の対象prototype紹介は32層・6720 channel、30 single-SiPM層と2 double-SiPM層、有感面約210×225 mm²、Cu/W absorber 3.2 mm、約23.4 X0を記す（L2/L4。資料の適用範囲は文献表参照）。一方、今回の物理解析は層0–29、6300 physical channel、MCも30層。32と30は一方を誤記として消さず、追加2層の位置・readout・除外理由・campaign間変更を著者確認事項にする。

current `S_MAP` のEBUdecode L17–50は6×36のmappingでchip 5のchannel 30–35を非接続とし、210 channel/層を定義する。pitchは5.3/45.4 mm。偶数層は(-y0,-x0)、奇数層は(-x0,-y0)、zは偶数 `1+(layer/2)*19.9` mm、奇数 `12.2+((layer-1)/2)*19.9` mm。memory番号は物理座標に反映しない。これを機械実測寸法と同一視しない。

SiPMを10 µm型と15 µm型に分け、current較正・digitizationでは層0–3,28–29を15 µm、4–27を10 µmとする。pixel数はMCで10000/4489。2020設計資料は10/15 µm、約10000/4500 pixelを説明するが、2023実装の全型番・serial・交換履歴までは確認できない。メーカー型番、bias/overvoltage、shaping time、ADCの実効校正、DAQ clock周期は一般値で補っていない。

### B.2 campaignと運転条件

Murata CALOR2024資料（L6、slide 6）はSPS H2（site 887）2023-04-26〜05-10、PS T9（site 157）05-17〜31を記す。2023 CALICE資料（L8、slide 7）は準備を含む可能性のある04-24、05-16開始を記す。PS raw pedestalは05-16から存在する。期間を一つに統一せず、setup/取得期間の定義を確認する。

L8 slides10/18はmomentumとしてSPSのμ100、e10–250、π10–120/350、p350 GeV/c、PSのμ10、e0.5–5、π1–15 GeV/cを記す。L6等の概略範囲と区別する。温度はL8 slides11/19にSPS17–23°C、PS21–26°Cと報告されるが、run別sensor温度や較正基準温度ではない。slideの範囲だけで各runのbeam設定を決めない。粒子directoryはe−、μ−、π−、SPSにはp+も存在する。資料にあるp−表記との対応は未解決で、pion/proton性能を測定済みとはしない。

L8（Yunlong Zhang, 2023 Beam Test of Sci-W ECAL and AHCAL Prototypes）のslide9はECAL 32 DIF/AHCAL 40 DIFとTLU同期、slides27–28はCherenkovによる粒子選別の使用を報告する。slide29には2023年前の2層再製作、温度dataの順序・autogain threshold・firmwareの更新がある。ただし今回のROOT sampleごとのtrigger/PID適用、効率、clock周期、設置角度・位置、beam spread、bias、sensor絶対校正は未確認。TOF使用は確認できない。装置/beamline写真が発表資料に存在しても、個別runの設定証明とは分ける。current MCのupstream scintillatorは外部triggerの実測geometryを自動的に再現したものとはしない。

### B.3 runs.csvの行の単位

rawファイルまたは解析input manifestの1行を「sample × 用途/condition × input path」の索引とした。同じrunのraw/decode、異なるtemperature/DAC group、dataと対応MCは別行になる。事象数を行間で足し合わせてはいけない。`run_id` は実データファイル等に明記された番号だけを記載し、MCは空欄と `associated_data_run` で対応させる。directory由来の「100GeV」は `beam_label` に保持し、energyかmomentumか未確定なら両物理量欄は空欄にする。chargeも粒子記録からの符号だけで実測PIDを意味しない。

raw取得の列挙、LED/gain入力、muon MIP入力、electron評価を区別する。pionはraw取得に加えL8 slides16/31に予備的response・contamination/absorber比較図がある。今回のサーバー調査ではそれに対応する採用可能なpion response/resolution一式を特定していない。未使用と断定せず「用途・性能評価未確認」とする。gainのelectron/muonはgain cross-check、muonはMIP・closure・beam mapに使われており、すべてを独立な性能sampleとして数えない。

## C. 較正手法と定数の来歴

### C.1 実際に参照されるtable

current `S_CONFIG` の最後のfactor_dirは `data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Analysis_edit/share`。pedestal=`pedestal2023_SPS.root`、HG/LG=`all_hl_electron2023_hlratio_v3.root`、MIP=`all_auto_muon_v4_trackfit.root`。digitizationのgainは `20230504_1425_Calib_LED_all.root`、thresholdは `ScECAL_BeamTest/analysis/result/threshold/threshold.root`。ROOT tree/schemaと小容量の全channel table抽出は `evidence/root_exports/` に収録する。全ROOTのコピー・全体hashは不要な大容量I/Oを避けて省略したものがある。

file名が2023/SPSであること、readerが読むこと、実際の生成run・選択・fitの来歴が確定していることは別である。特にPSへのSPS名定数の適用期間、較正温度、生成版とrun対応は未確定。コメント中のcosmic2020定数はcurrentの有効代入ではないが、旧 `calib_simulation` メモでは他年由来の依存として現れる。

### C.2 Pedestal

`S_PED` Pedestal.cxx L147–175はHitTag=0のnon-hit entry、Event_Time=0除外等を扱い、channel別HG/LG histogramを蓄積する。memory別fillはコメント化された箇所がある。L286–308にはmean±2 RMSからのGaussian fitがある。専用pedestal rawもPS/SPS双方に存在するが、現在のtableをそのrunから作ったと断定する実行記録は不足している。

current reader `S_CAL` L369–399はChnLevelのCellID、PedLowMean/PedLowSig/PedHighMean/PedHighSigを読み、memory別定数を使わない。保存treeの値branchはvector<double>、captured readerはvector<float>で読み、Extract.cxx L57–58はprocessCaliに丸めを合わせる意図を明記する。保存schemaとreader型を混同しない。digitizationはmeanを加えsigmaをGaussian幅に使い、reconstructionはmeanを引く。sigmaはvarianceやfit誤差ではない。table生成時のfit成功率・run温度・更新間隔は不明として残す。

### C.3 LED SiPM gain

`S_GAIN_README` と `S_COMB` は現在のgain候補抽出を説明する。memory 0、HitTag/GainTag、run/channel/temperature 0.5°C groupを分け、温度補正そのものは行わない。ADC 250–1400、quantile 0.5–99.5%、最大400 ADC幅、最低120 entry、gain探索6–40 ADC/p.e.という設定を持つ。LED DACの動作設定はinput manifestにunconfirmedが残る。

CombFit.ccのcommon spacingモデルはpeak中心 `phase+(j−1)g`、Gaussian幅 `g*width_fraction`。multi_gaussianではpeak面積を別々にprofileし、envelope_combではenvelopeで拘束する。bin積分Poisson deviance、背景、seedの半分/倍、Minuit収束/BICを扱う。spacingのfit error、個々のpeak幅、physical gain fluctuationは別概念である。bootstrap=0時の未評価stat_errorはNaNでありゼロではない。FFTのnull surrogate 199回はbootstrap誤差評価ではない。

full_grouped保存summaryは6,173,272 method rows、274,848 candidates、9,259,908 method pairs、317,239 cross-sample pairs、missing=[]、partial=false。これは解析tableの行数でありraw event数や有効channel数ではない。NUMERIC_COMPLETEを根拠に数値処理完了とするが、最終較正採用ではない。

selectionはDAC間・method間の整合性を使い、原則multi-GaussianとFFTの対称相対差5%以内、温度幅0.5°C、複数DACの中央値を用いる。strict版の3 DACと `full_grouped_relaxed_2dac` は別条件。後者の109条件×6300 channelからselected 7304 channel-conditionsという集計を、7304 unique channelと説明しない。gainErrorは条件付きfit error、ばらつき、stability等の最大値を含む品質幅で、統計1σではない。旧T_Event schemaのCellIDはlayer×10^6であり、通常のlayer×10^5と取り違えない。新しいselection結果がproduction gain fileへ採用された証拠はない。

### C.4 HG/LG intercalibration

採用table readerは `S_CAL` L465–530 / `S_EXTRACT` L103–181。実効変換は `LG=s HG+b` を前提に `A=1/s`, `B=−b/s`, `HG_equivalent=A LG+B`。finiteなsを閉区間[0.02,0.05]、finiteなbで受理する。現行readerはχ²による品質cutをしていない。

fallbackは30層の物理channel内のunique donorを10/15 µm別に集め、AとBを同じ母集団で算術平均する。`1/mean(s)`やoffset=0ではない。digitization側は `s=1/mean(A)`, `b=−mean(B)/mean(A)` に戻して整合させる。switch用XMaxはfiniteかつ>600ならXMax−600、欠測はvalid値の最小値を使う。

生成候補 `InterCalib.cxx` はraw HGを横軸・raw LGを縦軸にし、lower-case branch、slope 0.023–0.037、intercept 300–500等の制約を持つ。`ECAL/hl_ratio`にも別fitがある。しかし採用v3のuppercase schemaとの一致・pedestal subtraction・fit範囲・生成runを証明するproducerを特定できていない。これらのコードが存在するだけで、v3の生成方法として引用しない。

### C.5 MIP

production tableはMIP_FitのMPV等を読む。MPVはLandau成分のmode parameterであり、Landau-Gaussian convolutionの最大位置そのものとは限らない。旧 `analysis/mip` はSSA等を使うthreshold-aware fitを持つが、all_auto_muon_v4_trackfit.rootの完全な生成来歴は不明。旧FitMIP.Cのfixed_width名だけでwidth固定と判断しない（初期値/範囲設定とFixParameterは異なる）。

current reader `S_CAL` L401–463はMPV>0、finite、NDF>0、χ²>=0、χ²/NDF<=2（10 µm）または1.7（15 µm）。30層の物理channelからSiPM群別にdonor meanを作り欠測を補完する。過去の−10 sentinelをNDF=−10,χ²=−10の比で受理した旧問題はcurrent sourceでは修正されている。古いauditの平均81.576853/397.916766 ADCとclean候補97.098354/428.913481 ADCは旧tableに対する診断値であり、現在の最終性能を表さない。

新しい `analysis/mip_calibration` はtarget層を除くmuon track、path length補正（厚さ2 mm/track length）、threshold turn-onを含むLandau-Gaussianのbin積分Poisson deviance、複数seed、profile/安定性診断を用いる。外部thresholdは独立hardware測定ではなくelectron spectrum fitである。raw ADC下限40、最低200 fit hit、relative MPV error<=0.10、full covariance、正のNDF、correlation<=0.98、χ²/NDF<=2/1.7、turn-onからの距離、variant間5%等の基準は同梱README/config/sourceを参照する。threshold欠測補完や同定不能は直接測定として採用しない。

`all_muons_refit/COMPLETE` は2026-10-14 03:40:48+09:00。保存reader_statusは10 µmでaccepted 51 / fallback 4989、donor mean 129.324361056 ADC、15 µmでaccepted 588 / fallback 672、mean 463.311470672 ADC。ready_for_reader_review=1は論文採用の許可ではない。READMEの「全run未実行」は保存COMPLETEより古い記述。厳しい基準で多くがfallbackになることが主要な制約であり、全channelを測定済みと書けない。

## D. 再構成、温度、selection

### D.1 ADCからhit energy

current `S_CAL` Calibration() L560–582,659–709の順序は以下である。dataでも名前にsimulationを含むreaderを呼ぶ。

```text
HitTag == 0 → skip
layer = CellID/100000; chip = (CellID%100000)/10000; channel=CellID%100
chip==5 && channel>29 → skip
qH = rawHG − pedestalHG; qL = rawLG − pedestalLG
T = temperature vectorが空なら20°C、それ以外はtempReconstruction(layer,x,y,sensors)
alpha = 3.5/230 [15 µm] または1.6/135 [10 µm] /°C
C(T) = 1 − (T−20)*alpha
C*qH < 10 ADC → skip               # 等号は残す。LGを選ぶ場合にも先に判定
h = C*qH; l = C*qL*A + B           # offset Bは温度倍率の外
q = (rawHG < 2600) ? h : l         # 2600ちょうどはLG
E_hit [MeV] = 0.305*q/MIP_MPV
TotalEnergyDep = layer<30のE_hitの和
```

HG/LG switchはADC dynamic rangeの選択であり、SiPM pixel saturationの逆補正ではない。current pipelineに完全なSiPM saturation unfoldingを確認していない。Calib_HitはHit_Energy/HG_Energy/LG_Energy、XYZ、NewTemperature等を保存し、Event_Num=TriggerID。TotalEnergyDepはactive scintillatorのMIP換算信号和であり、入射energy、吸収体も含む全calorimeter energyとは異なる。

### D.2 温度とevent building

`S_MAP` tempReconstruction L74–105は16 sensorの座標を層向きに回転し逆距離加重する。current accumulatorは0初期化済み。sensorとの距離0はその値を返す。layer0–31、finite座標、sensor数16、finite温度を検証し異常はthrow。空vectorの20°Cはcallerのfallbackであり、異常vector全部を黙って20°Cに置換するものではない。係数の符号・20°C基準と較正tableの温度、sensor時刻対応の物理妥当性は別途未検証。MCは20°Cの固定16 sensorを出すのでC=1となる。

normal decodeの `S_DECODE` CatchSPIROCBag/FillChipBuffer/Decode（L263–352,367以降,543–624）はFA5A packet、36+36 ADC、BCID/chip、trigger、16温度を扱い、trigger差>30000で65536周期を補う。raw2Root normal entrypointとAutoGainDecodeを取り違えない。readout memory・重複hit・DAQ欠損を全runで監査したcut flowは未取得。event grouping実装があることはexternal detector同期・pile-up除去を実証しない。

### D.3 評価selection

current `S_SINGLE` は層<30、負energy除外、static bad_channels mask、層9/10のenergy-weighted CoG、指定window（5/10/20/200 mm）、HG/LG tail flag等を使う。beforeEventCutでも既にmask/CoG/tail処理を経ており「全取得event」の意味ではない。afterEventCutはさらにenergy–Nhitの事象選択後。current eventCutConditions L115–118は ((E−μE)/(3σE))²+((N−μN)/(3σN))²>1 を除外し、境界=1を残す。正確な境界と保存時のcommandはsource/ログを優先し、本文だけで再実装しない。

保存data100 GeV/5 mmのmetadataにはhl_tail_excluded=1、status=complete、after-CoG event数8348がある。tail源は `analysis/result/channel_response/hg_lg_tail/threshold/e_ssa`。common channelのsample別medianから残差を標準化し、低側tailを複数channel（>=3）持つeventをflagする。これはparticle PIDの代用ではない。

muon較正・closureはtrack層数、角度、残差、境界距離等の専用selectionを持つ。closure保存configにはmin_layers=7、max_angle=0.05、max_rms=3 mm等がある。electron、muon、pionで同一selectionと説明しない。全raw→decode→calib→最終selectionの一貫した事象数、PS/SPS別PID・contamination・複数粒子cutは未確認。保存histogram entriesはそのhistogramにfillされた数だけである。

## E. Simulationとdigitization

### E.1 Geant4 production

`S_G4MAIN` CEPCScECAL.cc L66–86はtime(0)+pid seed、single G4RunManager、**QGSP_BERT**を使う。custom PhysicsListの呼出しはコメントであり、そこにあるcut値を有効設定として採らない。現在のsetup環境はGeant4 11.1.3/ROOT 6.36.08等を参照するが、全保存sampleのbuild versionを保証しない。production時Geant4 version・cut値・seed・source hashの不足はmanifestに残す。

`S_GEOM` DetectorConstruction.cc L138–221,346–358: 30層、42×5 strip、scintillator 45×5×2 mm³、ESR 0.08 mm、kapton 0.05 mm、absorber 230×230×3.2 mm³、CuW密度16.17 g/cm³・Cu15%/W85%、PCB copper 0.24 mm+FR4 1.56 mm、air gap 0.64 mm等。別CuW2 materialが定義されてもplacementがコメントなら実装geometryに数えない。upstream Al 1 mm、scintillator厚41 mm中心z=−30 mmがある。実機trigger detectorの寸法証明ではない。reco centerとgeometry centerにも簡略化の可能性がある。

`S_STEP` SteppingAction/`S_EVENT` EventActionはRealCrystalのenergy depositをevent/channelごとに加算してtruthへ保存する。digitizationのStepDigiという名前でも、呼ばれる対象は各Geant stepではなく集約済みchannel deposit。optical photon追跡を用いた光収集/PDE simulationではない。

`generator_sps_e.sh` はdataのSSA trackFit由来beam mapからGetRandom2、z=−60 mm、方向(0,0,1)、GPS energyを使う。GPS energyはkinetic energy。予定100 file×1000 eventのloopは保存sampleの実数ではない。beam mapをdataで調整したMCは、そのbeam分布に関する独立検証ではない。別 `generator_e.sh` current point_noTriggerSciは(-6.7,4.9,-5) mm、明示seed系列、100 GeV等であり、e_ssa旧productionへ遡及適用しない。z=−5開始は上流scintillatorを通過しない設定で、geometryそのものを削除した意味ではない。

### E.2 current digitizationの式と順序

`S_DIGI` MCDigi.cxx StepDigi L199–295、PixelCorrection L175–195、PixelAva L141–170、Digi L304以降の値を以下に示す。生成時sourceが保存された小規模比較以外には、current値の使用を一律には認定しない。

|parameter|10 µm (layer4–27)|15 µm (layer0–3,28–29)|由来・留保|
|---|---:|---:|---|
|MIP energy [MeV]|0.305|0.305|コード定数。独立な測定誤差は未保存|
|gain [ADC/p.e.]|table、<=1なら8|table、<=1なら22.5|legacy gain schema、fallback便宜値|
|PDE|0.1|0.25|実装値。run別独立測定との対応未確認|
|crosstalk probability|0.1|0.15|隣接pixel cascadeの確率。メーカー測定そのものではない|
|pixel数|10000 (100²)|4489 (67²)|MCの格子モデル|
|光量uniform半幅|0.042|0.042|共通調整値。RMSは半幅/√3であり0.042そのものではない|
|ADC signal sigma|3√Nfired ADC|3√Nfired ADC|sigmaでありvarianceではない。gain fit errorとも別|

current light yieldは `L=(1−pCT)*M/G`。`Nγmean=floor((edep/0.305)*L/PDE)` → Poisson → uniform[(1−0.042)n,(1+0.042)n] → round → 各photonにBernoulli(PDE) → random pixel occupancy → crosstalk cascade → Nfired。occupancyは1回のchannel depositごとにresetし、同じpixelの再発火・recovery timeはない。隣接8方向から1方向を選びcascadeし、端で失う。位置依存は一様なpixel選択とglobal幅による近似で、strip内の実測光収集mapではない。

signal ADCをGaussian(mean=Nfired*G, sigma=3√Nfired)で揺らし、HGはさらにpedestal mean/sigmaを加えてround、LGはsignal/A+b+pedestalLGをGaussianでround。HGがchannel highSwitchを超えると4000 sentinel、LG>3000は3000へclipする。ApplyThresholdはGaussian(threshold+pedestal,sigma)を引き、ADC>=drawを通す。その後HG<pedestalを除外する。truth-hit channelだけをloopするのでnoise-only channelは生成しない。

DigiのgRandom->SetSeed(0)は固定seed再現ではない。出力Run_Num=0、BCID=−1、GainTag/HitTag=1、Hit_Time=−1、温度20°C等はsynthetic metadata。data runとして解釈しない。

### E.3 CT・threshold・調整の解釈

旧実装はL=M/G、currentは(1−pCT)M/Gである。M/Gは既に相関SiPM応答を含む可能性があり、追加CTとの二重計上は定数の定義次第である。currentの係数は単純な分岐の平均を意識した規格化だが、有限格子・端・occupancy・thresholdを含む最終responseが厳密に等しい証明ではない。primary photons、effective p.e.、fired pixelsを同じ変数として記述しない。

2026-10-13の小規模digitization比較は各energy 5/20/100 GeVのtruth file index 0,5等を使い、sources/binary/input/commandを保存している。旧MC、CT規格化変更、pixel数変更14400/19600、threshold+5/+10等を別conditionとして保存する。これらの感度比較は新しいpixel実測値や最適値の確定ではない。今回この比較を再実行していない。

threshold.rootのexpErfThreはCellID/threshold/sigma等。readerはchip meanで欠測を補完する箇所を持ち、χ²品質cutを読まない。独立hardware turn-on測定ではない。PDE/CT/光量幅/ADC fluctuationの決定sample、手動調整履歴、parameter covariance、独立validationの分離は不足している。合わせたdataへのagreementだけを予測能力の証拠にしない。

## F. 保存済み結果と限界

### F.1 評価量

energy sumの基本単位はMeV（表によってGeVへ表示済み）。入射energy GeVとのlinear scaleやsampling fractionを暗黙に導入しない。保存meanがGaussian μかhistogram meanか、resolutionがGaussian σ/μかRMS/meanかをmetric名・source_locatorで区別する。before/after selectionも別conditionである。

`S_SINGLE` のGaussian fitは保存histogramから初期mean/RMSを取り、概ねmax(15,mean−3sigma)〜mean+3sigma、次に±1.5sigmaという二段階を持つ。TF1自身の保存rangeだけでは実際の最後のfit区間を表さない場合がある。fit statusやcovarianceが保存されていないものはunknown。`S_MULTI` はσ/μの誤差をμとσの誤差から伝播するが、完全なfit covarianceの寄与・beam spread補正・全系統誤差を評価した証拠は不足する。

parameterizationは `sqrt(a²/E+c²)`、追加noise版 `sqrt(a²/E+b²/E²+c²)`。current MultiEnergyでは25<E<35、75<E<85をfitから除く処理があり、点の表示とfit使用を区別する。保存200 GeV dataは値/fit除外の明示記録がある。保存curveのparameterを最終採用値としない。

### F.2 結果群を混ぜない

|condition群|保存物・支持できること|限界|
|---|---|---|
|resolution Oct12|data/MC thresholdのCoG別保存histogram・fit・比較図・実行log。各5 mm energyの保存fitを数値抽出|生成calib/MCの完全な版とrun対応は未確定。data200除外、低energy不安定fitあり|
|Oct13 integrated response|5/20/100 GeVの条件付きactive signal比較、central/peripheral・hit数・RMS/Gaussian診断|CoG/energy window等を伴う診断。full energy scanの最終linearityではない|
|muon_mip_response/threshold_MPVwoCT|common/matched group別MIP closureとdata/MC比較|条件付きsample。旧Oct1 closureと混ぜない。fallback/individualを分ける|
|all_muons_refit|新MIP候補のfit判定・reader fallback集計・品質図|直接測定channelが少なく、全channel較正の完成を意味しない|
|gain full_grouped/selection|多method・DAC・温度の候補と品質比較|production table採用の証拠なし、統計誤差と品質幅が異なる|
|20260929 saved_fit_values|旧electron保存fitの転記、old condition scan|table内pathと現在のscript参照pathが異なり、元ROOTの移動/欠落・campaign・生成版を要確認。この旧tableは同梱するがresults.csvの2023性能点には統合していない|

Oct13 response_comparison.tsvの保存値の例：5 GeVでdata mean=0.142573793771 GeV、MC=0.179246892919 GeV、比=0.795405…、20 GeVでは0.72766446475/0.712937813645 GeV、100 GeVでは3.71219478559/3.39940923478 GeV。100 GeVのGaussian σ/μはdata0.04187790、MC0.02402077、RMS/meanはdata0.05018997、MC0.02633910。これらは保存されたselection条件付き信号であり、入射100 GeVの再構成energyとして3.71 GeVを無補足に書かない。未保存のuncertaintyを補わない。

data/MC差にはbeam profile・threshold・MIP補完・温度・HG/LG residual・sample selectionが関係し得る。原因寄与を分解した系統誤差budgetや独立validationは未完成。低energyでGaussian fitとRMSが大きく異なる保存例は、良好な分解能測定として採用しない。afterEventCutで細くなった分布だけからdetector intrinsic resolutionを結論しない。

results.csvは新規計算せず保存数値をmetric単位で転記したもの。大量channel定数は別JSON/tableとしてschemaとともに保存する。graphから目測した値は使わない。保存比較graphはGaussian σ/μを直接転記し、linear/logおよび同じ点の複数panelを独立測定として重複登録しない。追加列x_value/x_unit/x_definitionは保存graphの横軸とその定義を保持する。値が見つからないresolution parameter/linear residual等はunknownのままにし、異なる旧版をつないでenergy scanを埋めていない。

## G. 図・引用・貢献

`figures.csv` のfile行は既存図/資料、request行は未取得・必要内容である。filenameから推定したinput対応は推定と記す。保存PNG/PDFを無改変で同梱し、新しい性能図は作成していない。較正4例（pedestal/gain/HG-LG/MIP）、stack/mapping/DAQ/setup、linearity等で確実な対応がないものはrequestに残す。public PDF内の旧較正図を2023 runの測定例に流用しない。

|ID|一次資料・識別子|対象と参照箇所|
|---|---|---|
|L1|Yazhou Niu, *Design of Sc-ECAL Prototype for CEPC and Performance of First Two Layers*, JINST 15 (2020) C05036. DOI [10.1088/1748-0221/15/05/C05036](https://doi.org/10.1088/1748-0221/15/05/C05036), [arXiv:2002.01809v2](https://arxiv.org/abs/2002.01809v2)|PDF pp.3–4/fig.1は設計・readout、pp.5–6は初期2層の較正例。2023実機全体・2023性能の代替ではない。|
|L2|Jianbei Liu et al., *CALICE scintillator-SiPM calorimeter prototypes: R&D highlights and beamtests*, NIM A 1072 (2025) 170191. DOI [10.1016/j.nima.2024.170191](https://doi.org/10.1016/j.nima.2024.170191)|対象prototypeと2022/2023 campaign overview。両年を含むため2023結果へ無条件に帰属させない。|
|L3|Yong Liu, TIPP2023, 2023-09-06, *CERN Beamtests of CALICE Scintillator Calorimeter Prototypes*（配布PDF）|slide12 LED図のcreditにTatsuki Murata、slide24 campaign概要。正式な頁題は同梱PDFを優先。|
|L4|Wataru Ootani, *Status of ScW-ECAL*, DRD6, 2025-04-01（配布PDF）|slide4の32層/30 single+2 double等は後年の装置紹介。2025測定結果を今回のresultには入れない。|
|L6|Tatsuki Murata, *Status of high granular Sc-ECAL for future electron positron colliders*, CALOR2024（配布PDF）|slide6 campaign期間/setup写真。発表内容の担当範囲と、図の原データ所有者は別途確認。|
|L7|Tatsuki Murata, ILC-Japan会議配布資料（原PDFの表紙を参照）|campaign写真・説明の候補。日時/正式題名が未照合の項目は引用前に確認。|
|L8|Yunlong Zhang, *2023 Beam Test of Sci-W ECAL and AHCAL Prototypes*, CALICE 2023、BeamTest-20230928.pdf|slides7–11,18–19,27–29にschedule・DAQ・beam・温度・Cherenkov・更新。footerの日付表記に混在があり、正式会議日付は要照合。|

URL・取得可否・hashはmanifest.references / sourcesに収録する。学術引用は方法・装置の出典、source hashやjob logは解析来歴の証拠であり互いに代用しない。公開可能性・図の再利用許諾は採用時に確認する。

著者の貢献は、L3のLED解析図credit、本人名の発表資料から当該発表/図への関与を確認できる範囲に限定する。git authorやhome directory名だけで全較正・DAQ・simulationを本人の独創的成果としない。継承コード・共同開発・本人の修正/検証・共同者の測定を分けた貢献文の確定は著者に委ねる。

## H. 原稿対応・疑義・著者への質問

|項目|観察した実装と使用証拠|影響・未解決・次の確認|
|---|---|---|
|温度補間|current accumulatorは初期化済み。旧研究記録は旧不具合を扱う|どのcalib結果が修正前/後か、係数符号・基準温度・定数生成温度を確認。旧結果をcurrent修正済みとして説明しない|
|CT二重計上|current L=(1−pCT)M/G、旧L=M/G。Oct13小規模比較はsourceを保存|M/Gの定義、実測CT、response scaleを合わせたsampleを確認。current変更が全MCに反映済みとはしない|
|HG/LG|readerのA=1/s,B=−b/s、group fallbackは確認|v3 producerのraw/ped-sub axes、fit範囲、run、温度を特定。offsetと温度の非可換性がresponseへ影響|
|MIP|旧sentinel問題はcurrent readerで修正、新refitも保存|新候補採用の可否、fallback多数の解釈、thresholdモデルとrun安定性を確認|
|30層と32層|設計/実機紹介32と解析/MC30を併記|2023 as-built inventory、追加2層の位置/役割、除外の理由を確認|
|数値0.305/2600/10/0.042/3等|current箇所はC–E節に特定|source上の存在と過去production使用を分ける。由来・不確かさを独立に確認|
|resolution|保存値/figure/log、data200除外、不正低energy fitあり|どのconditionを論文候補にするか。fit安定性・beam spread・systematics未評価を結論の制約にする|

第3章 `scecal.tex` は設計、strip/mapping/readout、MCと実機の相違まで根拠付きで書ける。as-built32層・型番・DAQ/trigger・bias詳細は追加確認が必要。第4章 `beamtest.tex` はpipeline、current較正/再構成の式、保存run索引を説明できるが、実行版の対応・実験条件・PID/cut flow・最終定数は採用判断待ち。第5章 `beamtest_evaluation.tex` はdigitization順序、保存条件別の診断とdata/MC差を記述できるが、完成した精密性能測定・独立validationとしてまとめない。第9章の考察と第10章の結論では、調整/検証sample共有、fallback、低energy fit、系統誤差不足を落とさない。

著者への優先質問：

1. 論文の採用候補をどのproduction/condition（旧e_ssa、threshold、CT規格化変更、新MIP等）に置くか。明示的な採用記録はあるか。
2. 2023 PS/SPSの電子logbook、run list、as-built層構成、bias、trigger/PID、変更履歴を保管した場所と、装置/DAQ担当者は誰か。
3. pedestal2023_SPS、HL ratio v3、MIP v4 trackfit、20230504 gainの生成command/source/run/温度とPSへの適用根拠はどれか。
4. 本人担当・継承・共同開発の境界、転載する図のcredit、最終的に主張したい到達点をどう定義するか。

これらは資料だけでは決められない事項であり、今回の資料作成を停止するための事前確認にはしていない。

## I. 移送・検証・更新

snapshotと同じ親にある `.tar.gz` と隣の `.sha256` を対で移送する。受領先hostを仮定しない。最初にarchive checksumを確認し、thesis repositoryの `docs/research/snapshots/<snapshot_id>/` に展開して `sha256sum -c SHA256SUMS.txt` を実行する。snapshot原本は変更せず、`docs/research/scecal_context.md` を別に作り、本書とIDを参照する。

```bash
sha256sum -c scecal_2023_20261024_091634.tar.gz.sha256
mkdir -p docs/research/snapshots
tar -xzf scecal_2023_20261024_091634.tar.gz -C docs/research/snapshots
cd docs/research/snapshots/scecal_2023_20261024_091634
sha256sum -c SHA256SUMS.txt
```

CSVはUTF-8、空欄はmissing、validはtrue/false/unknown。各result_id→condition_id→sample/calibration/reconstruction/selection/evaluationとsource locatorをmanifestで照合する。元サーバーの絶対pathは来歴であり、閲覧には同梱relative pathを使う。抜粋/抽出は元ファイル全体と区別される。大容量ROOT/raw、実行binaryや外部softwareを全部は同梱しておらず、執筆環境で元解析を再現実行した証明ではない。

将来のsnapshotでは、source hashだけでなく実行source bundle、定数fileとschema/quality/fallback、input run/seed、温度、selection、result validity、source移動/欠落、図と数値の対応、採用判断の変更を比較する。同じcondition_idを異なる条件に再利用しない。追加解析が必要なら別作業として目的と適用範囲を定義し、本snapshotの保存値を上書きしない。

最終組立の補足：全channel定数はexport JSONとschemaへ、保存比較graph等はresults.csvへ収録する。linear/logおよび同じ点の複数panelは独立測定として重複させない。L3のPDFは取得失敗のためURLと参照箇所を保持する。
