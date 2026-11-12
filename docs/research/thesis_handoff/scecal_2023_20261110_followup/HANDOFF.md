# ScECAL 2023年ビーム試験の追加調査と執筆引継ぎ

2026年11月10日、移設先サーバーで読める解析資料を調べた追加調査。対象は2023年CERN PS・SPSのScECAL electron・muon・pion。10月12日の分解能比較、10月13日の診断、新MIP候補を分け、保存結果に結び付く条件を確認した。博士論文の第3章（装置）、第4章（実験・較正・再構成・選択）、第5章（simulation・digitization・応答・比較）への引継ぎを目的とする。

最大の訂正点は、10月12日の147分解能点が **tail除去・CoG cut後、ellipse cut前のbeforeEventCut** を用いること。10月13日の応答は別の抽出・energy窓を持つ診断である。新MIPは28 run・855,722 decode entriesを用いる候補だが、個別受理channelは639/6300で、正式採用は未確認。

元コード・定数・データ・保存結果は変更していない。新規simulation、digitization、全event再処理、再fitは行っていない。保存ROOTのgraph・TF1・histogram・primary branch、既存表・log・保存sourceの読出しによる調査である。この引継ぎ作成で追加したのは本directoryと同名ZIPだけ。

「確認済み」は記録・実装・保存値を直接読んだ意味で、再現性・物理的妥当性の保証ではない。「推定」「不明・未検証」と区別する。「採用済み／暫定／著者判断待ち」は別の軸であり、今回の主要系列は暫定・著者判断待ち。valid=trueは限定された判定を通る意味、falseは明示除外等、unknownは品質未確定。空欄・欠点を0としない。

旧snapshotは2026-10-24開始、10-26に組立・検証されたもので、source採取は逐次・非atomic。この追加調査は旧snapshot全文の置換ではない。執筆側のsnapshot ID scecal_2023_20261026_164835_chat とも別である。

## 受渡しファイル

- PROMPT.txt：元の執筆Codexに貼る指示。
- original_request.txt：今回のQ1〜Q8の依頼全文。
- tables/：保存値を既存表・保存抽出から転記したCSV。値を再fit・補間しない。
- manifest.json：同梱fileのSHA-256、CSV行数、収録範囲。

CSVの列名・source・object locatorを残して利用する。数値表の詳細はCSVを優先し、本文の表示桁を新しい測定精度と解釈しない。元のROOT、raw、旧snapshot全体、公開PDF本体はこのZIPに同梱していない。調査元の絶対pathが受取環境で存在するとは限らない。

|CSV|データ行数|収録範囲|
|---|---:|---|
|resolution_points.csv|147|8系列の保存Gaussian sigma/muとgraph誤差|
|response_5mm.csv|38|19 label×data/MCのGaussian μ/σ・誤差・before entries|
|resolution_counts.csv|152|19 label×data/MC×4 CoGのbefore/after事象数。先行チャットの直接ROOT読出し表から転記し、5 mm全38件を応答CSVと照合|
|resolution_curve_fits.csv|16|8系列×noise有無の保存TF1 parameter/error/χ²/NDF|
|integrated_response.csv|6|10/13の3 label×data/旧MC。後続診断にだけあるGaussian誤差は別source列|
|mip_reader_status.csv|4|旧新版×10/15 µmの受理・fallback。旧平均は保存集計値未発見のためCSVでは空欄|

旧MIP平均について本文に記した値は、先行調査で保存channel tableをreader条件に従って読出し集計した値である。CSVは保存済み集計表からの直接転記を優先しているため、旧平均を空欄とする。新しい再集計や補間で埋めていない。

## 移設先と参照略号

|略号|旧path|新サーバーの読める実体|
|---|---|---|
|R|/home/murata_t/ScECAL_BeamTest|/gpfs/group/ilc/users/murata/ScECAL_BeamTest|
|D|/home/murata_t/data_beamtest|/gpfs/group/ilc/users/murata/data_beamtest|
|A|旧D/analysis/ECAL_Analysis_LCIO|/gpfs/group/ilc/users/murata/data_beamtest/analysis/ECAL_Analysis_LCIO|
|S|旧D/simulation/CEPCScECAL_SML_Portable_update_new|/gpfs/group/ilc/users/murata/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new|
|B|旧R/docs/research/thesis_handoff/scecal_2023_20261024_091634|/gpfs/group/ilc/users/murata/ScECAL_BeamTest/docs/research/thesis_handoff/scecal_2023_20261024_091634|

以下のR/A/S/B等を含む参照は、この表の絶対pathへ展開する。

- C = R/analysis/result/resolution/comparsion/comparison_5mm.root。旧R/result/resolution/comparsion/comparison_5mm.rootの現在の読める結果。R/result自体は今回の確認時には存在せず、analysis/resultに実体がある。
- SC = B/evidence/captured/src_f4e156057998_SingleEnergyAnalysis.cc。
- ME = R/analysis/resolution/MultiEnergyAnalysis.cc。
- I = R/investigations/result/e_ssa_threshold_review_20260930/20261013/integrated。
- DI = R/investigations/result/e_ssa_threshold_review_20260930/20261013/digitization。
- M = R/analysis/result/mip_calibration/all_muons_refit。

対応根拠は研究索引・HANDOFF・保存source・log・結果の配置、および今回の現物照合。Cの147点すべてのx/y/ex/eyがB/evidence/comparison_canvas.jsonと完全一致した。ME、A/src/Calibration.cxx、S/Analysis_edit/src/MCDigi.cxxはB保存写しとSHA-256一致。現SingleEnergyAnalysis.ccはB保存写しと相違する。これらは生成binaryとの完全対応を保証しない。

Bのresults.csvは5871行、runs.csvは3400行。metric・用途・入力path別の索引であり、独立測定数やevent数ではない。執筆側のresults.csv 40行はfit係数・較正診断も含む抜粋で、「147点中40点」ではない。runs.csvの4例を全muon入力と扱わない。

M/event_audit.tsv、histograms.root、dependencies.tsvは旧NFS絶対pathへのsymlinkで、移設先のR/analysis/result/mip_calibration/all_muons/に読める引継ぎ元がある。リンクは修復していない。旧サーバーへの接続は前提としていない。

## Q1 各解析系列と使用データ

10月12日はelectron性能比較候補、10月13日は応答・ADC・幅の原因診断、新MIPは較正候補。正式採用の記録は確認できなかった。beamのGeV表記は以下ではfile/directory/graphのnominal label。実験資料のbeam設定は主にmomentum [GeV/c]。dataの実測momentum・実測energyと混同しない。MCのPrimaryEnergyは保存generatorでGPS kinetic energyに対応する。

### 10月12日のelectron入力

根拠：R/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/<E>GeV/inputs/data.tsv全行、同mc.tsv、resolutionのtail_source metadata、単energy log。tail manifest・file数・checked数が整合するが、10/12 CLIの全入力file一覧そのものの完全な恒久production manifestは未確認。analysis/tmp/analyze_data_*.txtは10/29後続実行で上書きされているため10/12の直接証拠にしない。

|label [GeV]|PS/SPS|e− data run|
|---:|---|---|
|0.5|PS|86,91|
|1|PS|84,85|
|2|PS|80,81,82,83|
|3|PS|78,79|
|4|PS|77|
|5|PS|63,64,75,76|
|10|SPS|110,111|
|20|SPS|96,97|
|30|SPS|137,138|
|40|SPS|98,99|
|50|SPS|135,136|
|60|SPS|100,101,102|
|70|SPS|133,134|
|80|SPS|103,104|
|100|SPS|105,106,107|
|120|SPS|108,109|
|150|SPS|184,185|
|200|SPS|186,187|
|250|SPS|188,189|

MCはS/Result_MC/calib/e-/e_ssa/threshold/<E>GeV/<E>GeV_*.root。各energyでchecked 100000 calibrated entries、tail除外0の保存記録を確認。生成取得数と区別する。100 GeVのdata生成logはR/analysis/log/resolution/singleEnergy/test-1326283.out:1–7、MCはtest-1326303.out:1–5。比較logはR/analysis/log/resolution/compare/test-1326323.out:1–9,4694、2026-10-12 03:58:52–06:25:12、exit 0。

raw取得にはPS5 GeVの87–90やSPS100 GeVの128・129もあるが、この比較へ使用したと補完しない。

### 10月13日の診断

dataは各file最大2000、MCは各file最大1000事象を、fileを10区間に分け各区間中央から抽出する。全取得eventの再評価ではない。保存sourceはI/inputs/analysis_source/IntegratedReview.cc:58–60。

|label|data campaign/run|MC file index|抽出data/MC|CoG20・tail後data/MC|energy窓後data/MC|
|---:|---|---|---:|---:|---:|
|5|PS 63,64,75,76|0,26,44,62,80,99|8000/6000|3073/2825|2733/2825|
|20|SPS 96,97|同上|4000/6000|1784/2740|1756/2740|
|100|SPS 105,106,107|同上|6000/6000|5048/4541|4999/4541|

全てelectron e_ssa/threshold。根拠：20261013/<E>GeV/baseline/inputs/{data,mc}.tsv、I/validation.tsv:2–7、I/summary.tsv:3,10,17,24,31,38,45,52,59,66,73,80。

### muon と pion

新MIP候補はSPS μ−100 GeV label、run29,30,31,33–57の28 run（M/inputs/manifest.tsv:1–28）。引継ぎ元all_muons/event_audit.tsv:2–29を単純合計するとdecode entries・sampled・unique_eventsは各855722、duplicate_identity_eventsとinvalid_temperature_crossingsは0。crossingsは6542834件だがtrack/channel交差件数でありevent数ではない。PS μ−10 GeV rawは存在するが、この28 runには含まれない。

muonはMIP較正、closure、channel応答・gain比較、beam map等に使う。これらを全て独立な性能sampleと数えない。

pionはraw取得と発表資料の混入診断を確認。例：PS π−1 GeV run31,33–36,43–45,53、SPS π−350 GeV run167–172,190,191。B/runs.csv:44–87,238–288および移設raw実体。論文採用可能なpion応答・分解能の完成系列は今回特定していない。未実施とは断定しない。2023実験資料pp.16,27–28にelectron/muon混入とCherenkov診断がある。

## Q2 結果ごとの較正と処理版

P=pedestal2023_SPS.root、H=all_hl_electron2023_hlratio_v3.root、M旧=all_auto_muon_v4_trackfit.root。

|系列|使用が確認できた条件|限界|
|---|---|---|
|10/12比較|既存単energy ROOTを読む。比較段階ではADC較正しない|比較logから元calib生成時の三定数を確定できない|
|10/12 tail判定側|P/H/M旧を使用|tail側の定数記録。元calib生成時の直接証拠とは区別|
|10/13 integrated|P/H/M旧と保存温度から保存HG/LG energyが再現される|table生成run・温度の物理的妥当性とは別|
|新MIP|all_muons収集histogramを再fit、定数snapshotを引継ぐ|新tableのproduction採用は未確認|
|旧SiPM gain|20230504_1425_Calib_LED_all.rootがCT診断commandに明示|新しい多手法gain選択への置換証拠なし|

根拠：R/analysis/result/channel_response/hg_lg_tail/threshold/e_ssa/100GeV/inputs/configuration.txt:9,18,20、I/inputs/analysis_source/IntegratedReview.cc:29–30,54–64、I/inputs/checksums.txt:8–10、M/inputs/command.sh:1、configuration.tsv:4,10,23,25,27–28、DI/inputs/ct_digitization_study.sh:54–55。

I/validation.tsv:2–7の保存closureは、data抽出8000/4000/6000、MC各6000事象。HG/LG最大差は約3.6e−15〜2.8e−14 MeV、選択信号和差は0。今回の再実行ではない。新MIPではなく旧trackfit＋現行補正式との一致を支持する。

生成来歴について、Pの生成run/温度、H v3のproducer/run/温度/fit範囲、M旧の生成run/選択/温度、旧LED gainの使用run/DAC/producerは未確定。候補コードの存在やschema一致だけで当該table生成法と断定しない。根拠はB/HANDOFF.md:82–113と保存source・table。

readerのMIP品質はfinite MPV>0、NDF>0、χ²>=0、χ²/NDF<=2（10 µm）/1.7（15 µm）。HLはfinite slope∈[0.02,0.05]、finite intercept、χ² cutなし。HL変換はA=1/s、B=−b/s。fallbackは30層接続channelの同SiPM群unique donorの算術平均。

|MIP版|SiPM µm|対象ch|個別受理|fallback|群平均ADC|
|---|---:|---:|---:|---:|---:|
|旧table＋10/10改訂reader|10|5040|2891|2149|97.0983536400|
|同|15|1260|987|273|428.913481257|
|新refit|10|5040|51|4989|129.324361056|
|同|15|1260|588|672|463.311470672|

新値はM/reader_status.tsv:2–3、旧受理数はM/summary.tsv:2,173。旧平均はB/evidence/root_exports/694fac608c4c_MIP_Fit.json全行をI保存CurrentParameters.hh:35–45,73–102の条件で読出し集計。10/12 runtimeの実使用fallback logとはしない。旧auditの81.576853/397.916766 ADCはsentinel等を扱う別reader条件で、混同しない。

HL保存tableの同条件読出しでは10 µmの受理4988/fallback52、平均A=36.0741284954、B=166.851219372 HG ADC。15 µmは1239/21、A=33.8609340907、B=165.864435743。根拠B/evidence/root_exports/c8bb349855ae_InterCalib.json、I保存CurrentParameters.hh:49–102。1/mean(slope)ではない。

現在A/run/global_config:29–30,96–100はPS/SPS切替と無関係に同三定数を指定。10/13 PS5 GeVも同三定数でclosureする。しかしPS全期間への共有の物理根拠・著者採用・全過去calib生成履歴は未確定。

新MIPは2026-10-14 00:15:25開始、03:40:48完了（M/run.log:1–3,184–185）。reader_loadableやready_for_reader_reviewは正式採用ではない。新gain多手法selectionも別候補（R/analysis/sipm_gain_selection/README.md:3,11–25）。

D/ECAL_data/analysed/2023/sps/summary.txt:39–43の旧calib_simulationは2020年由来三定数を記す。別較正版の依存として区別し、全結果へ割り当てない。同メモのlatest・重複calib2表記も確定production manifestには使えない。

## Q3 実際の選択

10/12はCalib_Hit入力。上流calibrationにHitTag=0、非接続channel、補正HG<10 ADCの除外、rawHG<2600でHG／2600以上でLGの処理がある（A/src/Calibration.cxx:659–709）。10/13抽出ではclosureを確認。10/12全calib生成版への対応は未完。

単energy評価の順序は以下。

1. 既存tail event maskを先に適用。
2. layer0–29を集計、Hit_Energy<0除外。0はこの条件では残る。
3. only_best=trueならstatic bad_channelsを除外。
4. layer9と10それぞれでenergy-weighted X/Y CoGを作り、両層の全4座標にabs(coord)<rを要求。
5. ここまでがbeforeEventCut。比較graphはこのGaussian fitを使う。
6. afterEventCutだけに、((E−μE)/(3σE))²+((Nhit−μN)/(3σN))²<=1を要求。μ/σはbeforeのfit。境界=1は残る。

根拠SC:42,45–51,115–118,178–203,423–481、ME:1028,1125–1174。r=5/10/20/200 mm、原点中心の正方形、半幅r・全幅2r、境界除外。200 mmは図ではnocutだがコード上±200 mm条件。tailのshower-cog20設定をresolution全CoGへ重ねない（HLTailSelection.hh:76はselected branchを読まない）。

tail定義はcommon channelでsample別medianを引いたHG/LG residualを標準化し、z<−5のdistinct channelが3以上でevent flag。R/analysis/resolution/HLTailSelection.hh:47–81と各energy tail inputs/configuration.txt。event identityはinput manifest coverageとRun_Num/Event_Time/Event_Numの一致を要求する（同:88–113）。PID cutではない。

static maskは2517 unique IDs、そのうちlayer0–29は2097。保存source配列の数であり、2023実機dead channel数ではない。B保存生成script src_ef2787785ec2_analysis_paralell_allCondition.sh:38はonly_best=1。SCは10/24採取で10/12 binary hashとの完全対応は未証明。metadata/logと一致する範囲の使用証拠として扱う。

trigger/PIDのrun別条件・効率・混入率、独立readout品質cutの有効適用は未確認。2023実験資料pp.27–28のCherenkov使用と個別ROOTの使用証拠を区別する。first-layer multiplicity除外はSC:463でコメント。

100 GeV/5 mmの保存数はdata checked105991、tail除外1972、before8348、after7705。MC checked100000、tail除外0、before10177、after9906。全energy・CoGの数は保存ROOTのbeforeEventCut/edep_<int(E)>GeVとafterEventCut同objectで確認した。data200 GeVは200/10/5 mmがincompleteでafter histogram無し（0としない）、20 mmはcompleteで49735→34515だが比較policyは全data200を除外。

10/13は各fileの部分抽出。CoG20は両層の全4座標を要求し、追加energy窓はellipseではない。bad mask前の全接続channel正energy和を、基準MC中央値の0.5倍超〜1.5倍未満に制限する。基準中央値は同6 MC files6000事象の接続channel符号付きenergy和から作る。5/20/100 GeV labelの保存中央値は0.190435085965／0.759516911546／3.65231923912 GeV。各energy baseline/reference.tsv:2、I保存PrepareReference.cc:13–18、IntegratedReview.cc:63–68。

## Q4 保存応答と分解能

tables/resolution_points.csvは全147点。C/resolution canvas pad1のTGraphErrors primitive indices1,3,5,7,9,11,13,15。系列順はdata200/MC200/data20/MC20/data10/MC10/data5/MC5 mm、点数18/19/18/19/17/19/18/19。重複panel/log図を別測定として数えない。metricはGaussian sigma/mu、無次元、全保存点valid=unknown、暫定。

除外：data200 GeV全CoGはC/resolution_value_exclusions、比較out:4–7、err:24586–24592により明示除外。data0.5 GeV/10 mmはμ=−0.915225 MeV、σ=36.283993 MeVで正mean条件を満たさず欠点（err:4、ME:139–148）。30/80 GeVは点を残すがcurve fitから除外。150/200/250はcurve fit範囲外。

tables/response_5mm.csvは全19 label×data/MC=38行、beforeEventCut/fit_gaus parameter1/2とparameter error。histogram mean/RMSではない。単位MeV。他CoGの平均応答を網羅した表ではない。data200の保存μ=36629.277、σ=10341.85 MeVは失敗・未完了fileの値でvalid=false。その他はfit status/covariance未保存でunknown。data0.5/1 GeVはN=8/29でμ誤差がμを上回るため、良好な性能測定として採用しない。

信号はactive scintillatorのMIP換算和で、入射energy推定値とは別。100 GeV/5 mmのdata μ=3913.0788±1.35059 MeV、σ=82.477719±1.70316 MeV、MC μ=3439.2482±0.93273、σ=63.883788±1.17113 MeV。10/13診断の3.712 GeVとはselectionが異なる。

Gaussian fitはまずmax(15,μseed−3σseed)〜μseed+3σseed、次に第一fitのμ±1.5σ（下端15 MeV）。200 GeVはseed μ=8000、σ=300 MeV。SC:581–605。TF1保存rangeだけでは最終fit区間を示さない。

10/12分解能誤差はsqrt[(δσ/μ)²+(σδμ/μ²)²]、μ–σ covariance無し、x error0（ME:1168–1174）。beam spread差引き・総合systematicsを含む証拠はない。exit0はfit成功・品質の保証ではない。

tables/resolution_curve_fits.csvはnoise無しsqrt(a²/E+c²)とnoise有りsqrt(a²/E+b²/E²+c²)の各8系列。Cのresolution/resolution_with_noise canvas pad1、res_func_low_0〜7／res_func_low_noise_0〜7。a[√GeV]、b[GeV]、c無次元、fit status unknown。fit範囲0.1–130 GeV、25<E<35と75<E<85をRejectPointで除外（ME:1337–1382,1457–1503）。保存rangeは後に描画用変更がありfit範囲とは別。noise無しcは符号無拘束、凡例はabs(c)。noise有りa∈[0,2],b∈[0,5],c∈[0,1]。大きなχ²・極端なerrorがあり無留保に最終係数として引用しない。

### 10月13日 integrated

tables/integrated_response.csvは全3 label×data/旧MC=6行。平均は選択eventの算術平均信号和[GeV]、RMSは母集団分散からのRMS/算術平均、Gaussian列は別fitのsigma/mu。

|label|sample|N|算術平均GeV|RMS/mean|Gaussian sigma/mu|
|---:|---|---:|---:|---:|---:|
|5|data|2733|0.142573793771|0.124471307782|0.0863037613705|
|5|old MC|2825|0.179246892919|0.0762337623799|0.0743273870002|
|20|data|1756|0.727664464750|0.0839483362257|0.0716787115578|
|20|old MC|2740|0.712937813645|0.0430572440382|0.0412421404612|
|100|data|4999|3.712194785590|0.0501899703144|0.0418778962159|
|100|old MC|4541|3.399409234780|0.0263390983202|0.0240207732510|

I/summary.tsv:17,24,45,52,73,80、response_comparison.tsv:2–4。data/MC平均比は0.795405／1.02066／1.09201。平均/RMSの誤差は未保存。

元integrated TSVにGaussian誤差も無いが、後続DI/summary.tsv:10–11,22–23,34–35のfrozen_old_cog20_energy_window/full_contextが同じevent表・選択でN/平均/RMS/比一致。保存Gaussian誤差はdata/MCで5 GeV:0.0030064643698428/0.0026183376709879、20:0.0035864099705673/0.0014276409528524、100:0.0011869526962809/0.00074233035558179。全fit status0だが物理的有効性はunknown、採用待ち。この誤差はμ–σ covariance込み（DI/inputs/reviewer_source/investigations/e_ssa_threshold_review_20260930/DigitizationImpactReview.cc:28,47,55）。今回の再fitではなく既存保存値の読出し。

診断Gaussian fitは240 bins、平均の0.3〜1.7倍histogram、mean±3 RMS、次にGaussian μ±1.5|σ|、下端0.015 GeV（I保存ResolutionImpact.cc:34–36、DI保存DigitizationImpactReview.cc:28）。μ/σ各々・最終数値fit区間は当該表に未保存。

### 後続の別系列

10/29のstrip選択はhitのabs(X),abs(Y)<22.5 mm、CoG event cut無効。R/analysis/log/resolution/singleEnergy/test-1350144.out:1–5とROOT metadataで確認。100 GeVのdata before104019→after82189、MC100000→77763。data/channel22.5mm/200mmとmc/threshold/channel22.5mm/200mmにdata19 ROOT・MC9 ROOTを確認（全completeは未検証）。10/12 CoG147点と混合しない。現sourceのstrip対応・分母保護・fit下限変更等を10/12へ遡及しない。

## Q5 実機32層と解析30層

B/assets/literature/calice_beamtest2023.pdf p9は2023 ECAL32 DIF、AHCAL40 DIF、TLU同期を記す。全層のrun別稼働率までは確定しない。

B/assets/literature/murata_calor2024.pdf p5には「最後の2層にdouble-SiPM readout部」「90 mm strip、標準45 mm」と明記。p3はghost hit問題への対策・位置情報の目的。型番S12571-010P/S12571-015P、pitch10/15 µm、pixel10000/4489もp5。これは2024年の当該prototype紹介で、2023全runのas-built・稼働証明とは分ける。

解析の0-based layer0–3,28–29は15 µm、4–27は10 µm。A/src/Calibration.cxx:25,685、MCDigi:226–243。calibrationは32層hitを保存・monitorしつつTotalEnergyDepだけlayer<30を加算（Calibration.cxx:696–715）。resolutionも0–29、MC採取sourceも30層。除外範囲は確定するが、追加2層の不動作・別目的・single-SiPMだけの評価等、物理的除外理由は不明。

後年ootani_2025.pdf p4の32層、30single＋2double、6720chは参考。2020設計論文は30層開発・初期2層の別時点。2025測定値や2020装置性能を2023結果へ流用しない。

CALOR2024 p9はcosmic→beamで全channel biasを+0.5 Vと記す。絶対bias・channel別Uov・PS/SPS間変更は未確定。R/analysis/sipm_gain_pool/README.md:3はユーザー確認に基づきPS内・SPS内それぞれbias/ASIC共通と宣言するが、PS=SPSではない。

2023資料p29は2層再製作、温度data順序修正、autogain threshold変更、firmware/hardware更新。再製作2層をdouble-SiPM2層と同一視できる証拠は未確認。

## Q6 温度補正

C(T)=1−α(T−20)、α10=1.6/135 ℃⁻¹、α15=3.5/230 ℃⁻¹。pedestal差引きHG/LGにCを掛け、LG offsetはCの外、次にMIP除算（A/src/Calibration.cxx:677–705）。T>20で信号を小さくする負符号は実装確認済み。

20℃はC=1の参照値。係数の一次測定・fit・出典、符号の物理根拠、旧MIP生成温度との整合は未検証。分子をADC/℃、分母を20℃MIP等と断定しない。新MIPは28run温度混合を現行Cで補正したfit。保存Collect sourceには物理的20℃外挿の検証ではない旨がある（B/evidence/captured/src_73556847c608_Collect.cc:63–86）。

normal decodeは科学data packetの16温度を同eventへ格納し128で割って℃化する。外部温度logとのtimestamp joinではない（A/src/Decode.cxx:204,263,415,467–474,604–608）。空間は16sensorから層向きを考慮した逆距離補間（A/src/EBUdecode.cxx:74–105）。run別sensor順序・絶対時刻・sensor較正までは未監査。空vectorはcallerで20℃、MCは固定20℃でC=1。

10/13 CoG20・tail後・energy窓前のevent平均温度は5/20/100 GeV dataで28.9740384955／23.5131122134／22.9929751002℃、MCは20℃。Nは3073/1784/5048と2825/2740/4541。event内energy加重温度をevent平均した量で、全run温度ではない（I/summary.tsv:3,10,31,38,59,66）。

10/13抽出では式使用をclosureが支持。10/12全calib生成時への適用対応は未完。別の温度感度試験はUov20=4.5/4.0 V、降伏電圧温度係数60mV/℃等を仮定する（R/docs/research/e-ssa-threshold.md:108–110）。これを測定済み絶対Uovや係数の一次出典としない。

2023資料SPS17–23℃/PS21–26℃、CALOR2024のbeam25–29℃、診断PS5 GeV約28.97℃は対象・期間・sensor量が揃わないため一つに統合しない。

## Q7 比較MCのproduction

旧e−/e_ssa/thresholdが10/12比較および10/13 integratedの旧MC。旧保存source・記録は光量M/Gを支持するが全fileと生成binary/sourceの完全対応は未確定。10/13 CT修正は(1−p)M/GとCT cascadeを使う別診断。

DI/inputs/manifest.tsv:2–7とct_digitization_study.sh:27–31により、CT修正truthは各energy file0/99・各1000事象。「indices0 and5」は6file manifest位置で実file名0/5ではない。

S/Result_MC/generate/e-/e_ssa/{5,20,100}GeV/<E>GeV_0.root のMC_Truth各1000 entryについてprimary branchを読んだ。全PDG11、PrimaryEnergy各5000/20000/100000MeV、z=−60mm、方向(0,0,1)。GPS kinetic energyに対応（B保存PrimaryGeneratorAction.cc:94–109）。少なくともこれらは単色で、11/10追加のmomentum spreadを持たない。

|項目|確認内容|証拠と限界|
|---|---|---|
|Geant4/physics|11.1 patch3、QGSP_BERT|S/CEPCScECAL/generate/out/test-1303375.out:1,10,18。旧10/10の100GeV campaign log。個々fileのjob直接対応未確定|
|beam|data SSA trackFitのh_hitmap_smoothをGetRandom2、z−60、方向z|B保存PrimaryGeneratorAction.cc:85–91,115–128、旧log:12804でmap読込成功。全energyのmap来歴未完|
|geometry|30層、42×5strip、45×5×2mm³ scint、230×230×3.2mm³ CuW、密度16.17、Cu15%W85%|B保存DetectorConstruction.cc:138–141,173–221。全旧binary対応ではない|
|受動材|ESR0.08、Kapton0.05、PCB Cu0.24＋FR4 1.56、airgap0.64mm|同source|
|上流|Al1mm中心z−2、scint41mm中心z−30|同:346–358。実機trigger実測geometryの証明ではない|
|hit|step edepをevent/channelへ集約後digitize|B保存SteppingAction.cc:99–105、EventAction.cc:158–172。StepDigiは各Geant step独立ではない|
|乱数|Geant4 time+pid、digitizer SetSeed(0)|固定seed完全再現条件ではない|

旧10/11 digitization logはS/Analysis_edit/job/scan/e-/e_ssa/5GeV/out/5GeV_0.root-1320998.out、20GeV/out/20GeV_0.root-1320198.out、100GeV/out/100GeV_0.root-1319598.out。各:3–20でthreshold.root読込、pedestal6300、MIP受理15/10µm=987/2891・補完2422、gain records1260/5040、threshold records1109/4853を確認。threshold群平均15/10µm=195.254/69.7125ADC、保存sigma平均15.7308/8.83222ADC。引数表示threshold0をmap無効化と解釈しない。

point/point_noTriggerSciは別beam系列。z−5開始では上流scintを通らないがgeometry削除を意味しない。11/10 commit7fb610dのmomentum spread追加や現在macro/software環境を旧productionへ遡及適用しない。

## Q8 Parameter決定と独立性

|parameter|10/15µm|定義と来歴|
|---|---|---|
|PDE|0.10/0.25|固定実装。独立測定/外部引用/調整sample対応不明|
|CT|0.10/0.15|隣接pixel cascade確率。run別実測値・誤差未確定|
|pixel|10000/4489|100²/67²。recovery・再発火無し。14400/19600試験は感度診断|
|光量幅|0.042/0.042|Poisson後photon数へのuniform相対半幅。sigmaではない。元測定・調整energy不明|
|ADC幅|3sqrt(Nfired)ADC|Gaussian sigma。係数3の測定・調整根拠不明|
|gain|旧LED table|CT診断command明示。旧log平均15/10µm=23.4165/13.849ADC/p.e.。新selection採用不明|
|MIPenergy|0.305MeV|コード定数。beam energy sampling補正ではない|
|threshold|channel map|electron spectrum fit、chip平均補完。独立hardware turn-onではない|
|光量norm|旧M/G、診断(1−p)M/G|CT平均増幅を意識した変更。有限pixel/threshold後MPV一致保証ではない|

根拠：R/analysis/result/channel_response/muon_mip_response/default/inputs/current_MCDigi.cxx:226–240、DI/inputs/digitizer_source/MCDigi.cxx:141–195,226–262、Extract.cxx:242–302,306–405、MCDigi.h:21、DI/inputs/ct_digitization_study.sh:54–55。

threshold producer候補は同じ2023 PS/SPS electron全energy群を使う（R/analysis/threshold/threshold_allGeVs.sh:6–7,49–61）。HitTag後ped差引HGを（指数成分＋定数）×erf turn-onでfit（threshold_finding.cc:145–160,1316–1354）。expErfThreのthresholdErr/sigmaErr/chi2ndfは保存するがdigitizerはthreshold/sigmaだけ読みfit誤差伝播・χ²品質cutしない。sigmaはturn-on幅。

data由来beam map・thresholdは比較dataと完全独立とは言えない。全run重複一覧・holdout記録未確認。PDE/CT/0.042/係数3の決定energy・観測量・目的関数・共分散を結ぶ記録は未確認。CT修正は同truthによるpaired感度試験で独立beam検証ではない。光量×1.10/1.05もmuon MPV不足を参考にした診断候補で最終決定ではない。

integrated表は信号和をそのまま比較。shape図は各sampleの平均でEを割り1/N規格化するため絶対応答差が消える（I保存IntegratedReview.cc:63–65,89–96）。10/12のunit-area normalization・rebinをenergy scale較正と呼ばない。前段の調整履歴は未確定。

PDE/CT/光量幅/ADC係数/threshold補完等の総合systematics budget、旧分解能のbeam spread差引きは未確認。限定的bootstrapはある：DI/comparison_methodology.tsv:5、DI/model_scan/model_methodology.tsv:6、joint_response/digitization/joint_methodology.tsv:6。400回paired100-entry-block、固定2truth files内で層別、dataも固定source内で別再標本化、16–84%区間。file間/model間の系統差を覆わない。

## 執筆への適用と残る判断

第3章はstrip/readout/mapping、装置紹介の32層と解析30層を区別して書ける。第4章は保存結果に結び付くtail/CoG/ellipse順序、10/13旧三定数closure、新MIP候補方法とfallbackを記述できる。第5章は旧MCとCT修正診断を分け、保存値と条件付きdata/MC差を示せる。

結果の説明に残る不足は、旧定数生成run/温度、10/12入力calib生成版、各MCとsource/binary/geometry/seedの一意対応、run別trigger/PID/readout/cutflow、温度係数の一次測定、parameter調整sampleと独立検証、総合systematics。見つからないことを不存在・未実施と断定しない。

著者判断が必要なのは、掲載主系列・CoG・before/after、data200と低energy不安定fit・大χ²曲線の扱い、新MIP/新gain/CT修正の採否とfallback許容、as-built追加層の稼働・除外理由と絶対bias/PID、温度の物理基準・符号と旧MIP温度整合。今回の引継ぎは採用の承認ではない。
