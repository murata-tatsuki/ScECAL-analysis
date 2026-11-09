ScECAL muon MIP規格化・校正の再調査：全28 run統合結果
============================================================

作成日：2026年10月1日。対象：2023 SPS、100 GeV mu−のdataと対応MC。解析結果：all_muon_runs_20260930。Slurm投入ジョブ：1299387。保存ログのCompleted、STATUS.jsonのcompleted=true、出力ROOTの読み取り、PDFのページ数を確認済み。本レポートの作成で再fit・再digitization・本番校正コードの変更は行っていない。

1. 当初確認したかったことと、今回の到達点
============================================================

electron showerでは光量・SiPM飽和・shower形状の影響が重なるため、まずmuonで「MC edep → digitization後HG ADC → 校正後MIP応答」を確認することが目的だった。入力MIP値とはdigitizerが採用するADC/MIPの定数であり、生成edepの基準0.305 MeV/MIPとは区別する。

当初の確認項目	今回の実施内容・到達点
data・MCのtrack／位置／通過長をそろえる	共通track cut、2 mm通過長補正、run・channel・位置・再構成傾きの共通bin重みを適用。MCにない真の入射角分布までは再現していない。
生成MPV／入力MIPをchannel別に測る	同じLandau×Gaussian fit手順で評価。MCで2,125 channelが有効。入力より高い傾向を確認。
個別／fallback、10／15 µmを分ける	生成側と校正側の採用区分を別々に保存。図の色分けは旧校正側の区分。
校正後にdata・MCが同じ基準へ戻るか	共通304 channelで比と各試料の絶対MPVを比較。15 µm個別に差が残り、fallbackでは双方が1からずれる例を確認。
不正なMIP fitを除いたfallback平均	正の有限MPV、NDF>0等を要求した平均を算出。旧平均への不正値混入を確認。
採用値変更と同じhit集合のenergy変化	data・MCとも2,422 channelで採用値が変化。同じ保存hit集合でdata −2.18%、MC +28.82%。
crosstalkを含む最終応答が入力に合う一次光電子数	規格化再検討の根拠は得られた。一次光電子数変更・crosstalk ONでの再digitization・再検証は未実施。

結論は二段階である。校正定数の採用条件に問題があることは具体的に確認できた。一方、その問題を取り除いてもMIP規格化の不一致が残り、特に15 µmの個別校正channelでdataとMCの差が見える。残る差の原因をcrosstalkだけに特定したわけではない。

2. 入力・統合範囲・処理の健全性
============================================================

対象run：29, 30, 31, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57。計28 run。Run32の入力は存在せず、発見された28 runからの除外は0。--runs all --max-tracks 0 --allow-recoveredで処理した。「全件」は保存trackの全entryであり、保存されていないtrackを復元したという意味ではない。

項目	data	MC
保存track数	204,215	1,196,321
選択track数	120,478	805,256
幾何選択後MIP hit数	2,158,464	18,799,042
共通binで重み>0のhit数	372,505	3,409,990
共通binに残る割合	17.26%	18.14%
共通binの重みの総和	363,196	363,196

raw・校正後・trackのイベント／hit対応について、欠落・曖昧な対応は0件。旧校正候補を使った保存HG値の再現不一致も0件。総energyの比較は選択されたdata 120,478イベント、MC 805,256イベントすべてで成立した。出力flat.root、matched_hits.root、spectra.rootは復旧なしで読み取り可能。

入力MC trackFit3Dは28 run中27 runでROOTの復旧が必要だった。Run44だけは復旧警告なし。読み取れたentryを診断対象に含め、diagnostic_recovered_inputs=trueとしている。解析処理の完了と、入力ファイルの完全性の保証は別である。

参照定数：all_auto_muon_v4_trackfit.root（data校正・MC校正・digitization候補）、pedestal2023_SPS.root、threshold.rootのexpErfThre。現在のソース・定数はinputs/にhash付きで記録されているが、過去のMC生成時の版が完全に同一だったことを証明したものではない。

3. 何をそろえて、どうfitしたか
============================================================

raw、校正後、trackは各runのイベントIDと時刻で対応付けた。trackFit3Dの最終パラメータが有限、NDF>0、各向き7層以上、投影角の絶対値<0.05 rad、sqrt(保存目的関数/NDF)<3 mmを要求した。後者は統計的なχ²そのものではなく、距離に基づくtrack品質の尺度である。

MIP用hitはtrack採用channel、raw HG<2600、重複のない物理channel等を要求し、45×5×2 mmのstripを厚さ方向に通過する条件を適用した。端から長手1 mm、幅0.25 mmの余白を設けた。通過長は最終傾きから L=2√(1+sx²+sy²) mm と再計算し、ADC・edepに2/Lを掛けて2 mmにそろえた。

run、channel、strip内長手位置、幅方向位置、sx、syのbinごとに、data・MCとも3 hit以上ある領域だけを採用した。bin幅は5 mm、1.5 mm、0.01、0.01。各試料の重みはmin(Ndata,Nmc)/Nsであり、ADCやenergyを使って応答を一致させてはいない。その後、全runの重み付きhitをchannelごとにまとめて一つの分布をfitした。run別MPVを平均した解析ではない。

fit関数はLandauとGaussianの畳み込み。ここでMPVと呼ぶのはLandau MPVパラメータであり、畳み込み曲線自体の最大位置はconvolution_modeとして別保存している。data・MCで関数・bin・fit範囲の決定手順を共通化した。

fit下端はthreshold+3σ等から決めた共通の効率plateau範囲を使う。fit範囲内の有効統計100以上、fit成功、良好な共分散、正で有限のMPV、NDF>0、χ²/NDF≤3、パラメータ境界でないこと、MPVがfit範囲内にあること等を要求した。MPVがthreshold付近にあるchannelは未判定となり、有効channelの分布には選択による偏りがあり得る。

ADC段階は A=(HG−pedestal)×2/L。校正後は現行温度係数を掛けたA20のfit MPVを採用MIP定数Cで割って、M=MPV(A20)/C とした。温度補正は現行の1−(T−20)×係数（10 µm：1.6/135、15 µm：3.5/230）を再現した。共通bin重み付け後の平均温度はdata 22.823℃、MC 20.000℃であり、温度そのものを同じ分布にしたわけではない。温度補正式と入力定数の基準温度の妥当性は、残る確認項目である。

4. beamの入射角・位置について確認できたこと
============================================================

保存MC truthの28 run、計2,800,000 primaryはすべて方向(0,0,1)、source z=−60 mmだった。MC生成にdataのbeam angle分布は入っていなかった。現在の生成コードはrun別h_hitmap_smoothからx,yをサンプリングする。この位置mapの作成コードはdata trackのz=0切片にpitch幅の一様な揺らぎを加えている。

したがって、「入射位置の分布を調整していた」という理解は概ね対応するが、各イベントのinitial hitを一対一に一致させた方法ではない。真の入射角分布や位置・角度相関も再現していない。共通の再構成傾きbinでの比較はこの差を抑える手段であり、MCに存在しない真の角度分布を作る手段ではない。

重み付け後にもbin内部の差は残る。例えばsyの平均はdata 0.001334、MC 0.000757、標準偏差はそれぞれ0.003482、0.002708だった。通過長の平均はdata 2.0000167 mm、MC 2.0000109 mmであり、この選択での平均通過長差そのものは7〜12%の応答差を説明する大きさではない。ただし角度・位置に伴う選択効果全体を否定したわけではない。

5. 生成されたMIP MPVは入力値を再現したか
============================================================

指標R_ADC=MPV(A_MC)/C_digi。1なら入力定数を再現し、1.10なら生成されたADC MPVが入力より10%高い。下表はMC ADC fitが有効な全channelの集計であり、dataでもfitできたchannelに限定した集計ではない。各channelを同じ重みで数えた中央値を示す。

種類	旧校正側の区分	有効channel数	R_ADC中央値	1からの差
10 µm	individual	764	1.074869	+7.49%
10 µm	fallback	649	1.066786	+6.68%
15 µm	individual	638	1.118789	+11.88%
15 µm	fallback	74	1.084395	+8.44%

重要な区分：closure_summaryの青／橙は旧校正側のindividual／fallbackである。「生成側でfallbackを使った」という区分とは別。今回、生成側individualは2,124 channel、生成側fallbackは1 channel。後者は10 µmかつ旧校正側individualに含まれ、R_ADC=1.093894。生成側individualだけに限った10 µm・旧校正individualは763 channel、中央値1.074581である。図の764 channel・1.074869との差はこの1 channelによる。生成側fallbackの応答一般を1 channelから評価することはできない。

両種類とも、入力より高いADC MPVを持つ傾向が確認できる。特に15 µm・個別校正の中央値は約+11.9%。ただしchannel間のばらつきであり、この中央値の差を、そのまま全channel共通の確定補正係数や統計的有意度として扱わない。

種類	ADC・truth共通channel数	edep MPV / 0.305中央値	R_ADC中央値	R_ADC / (edep MPV / 0.305)中央値
10 µm	975	0.990231	1.068258	1.080292
15 µm	588	0.983165	1.115314	1.134207

同じchannel集合でedep MPVは基準より約1〜2%低い。一方ADC応答は高いため、「muonのedep自体が高いだけ」という説明は、この選択結果からは支持されない。edepから光電子数・ADCへの変換やMPVの規格化を調べる根拠になる。最後の列はchannelごとの比の中央値であり、前の列の中央値同士を割った値とは限らない。

このtruth分布も検出・選択済みhitに条件付けられており、無条件の全生成muonのdE/dx分布ではない。また、crosstalk ON/OFFを切り替えた対照実験ではなく、現行のcrosstalkを含む想定のdigitization出力に対する診断である。光電子数の定義・gain・飽和・threshold・fit形状等のどれが支配的かは、これだけでは分離できない。

6. 校正後にdataとMCが1 MIPへ戻るか
============================================================

以下はdata・MCの両方で有効なfitが得られた同じ304 channelを使う。clean校正は両試料に同じ方針で適用した。各列はそれぞれの中央値であるため、dataの中央値／MCの中央値と、channelごとのdata／MC比の中央値は一般に一致しない。

種類	旧校正側の区分	共通channel数	data clean MPV [MIP]	MC clean MPV [MIP]	data / MC中央値
10 µm	individual	28	1.032682	1.009748	1.029582
10 µm	fallback	26	1.227777	1.226047	1.054522
15 µm	individual	242	0.970505	1.127849	0.857001
15 µm	fallback	8	0.712039	0.783775	0.938319

15 µm・個別校正242 channelでは、data約0.971 MIPに対してMC約1.128 MIP、data/MCの中央値約0.857である。fallbackだけが不一致の原因ではないことを示す。10 µm・個別校正28 channelではdata約1.033、MC約1.010で比較的近い。全MC集計の約1.075と異なるのは、dataでもfitできたchannel集合に限定されるためであり、10 µm全体が一様に+7.5%と断定できない。

10 µm・fallback26 channelではdata約1.228、MC約1.226 MIPとなる。dataとMCが近くても双方が1 MIPに戻っていない。逆に15 µm・fallback8 channelでは双方が1未満である。少数channelのため一般化は慎重に行う必要があるが、単一のsensor種平均を各channelに代入する方法では絶対MIP応答を保証できないことが見える。

共通校正係数Cで割った場合、(MPV_data/C)/(MPV_MC/C)=MPV_data/MPV_MCとなり、Cは相殺される。同じhit・fit結果を使う今回の比較では、fallback平均を変えてもchannelごとのdata/MC比は基本的に変わらない。校正修正の評価には比だけでなく、各試料のMPVが1へ戻るかを必ず併記する。

7. 負のMIP値、fallback平均、採用値の問題
============================================================

legacy校正はχ²/NDFの上限だけを見ており、MPV>0やNDF>0の検査がない。今回、負の個別校正値として採用された565 channelはすべてMPV=−10、NDF=−10、χ²/NDF=1だった。不正なfit情報でも比が1になるため、10 µmの上限2、15 µmの上限1.7を通る。その値はfallback平均にも、そのchannel自身の校正にも使われていた。

一方、MC生成側のExtractは負MPVを除外するが、同じχ²/NDF cutは適用しない。生成側にはMIP≤5のruntime代替値（10 µm：90 ADC、15 µm：400 ADC）もある。同じ定数ファイルを使っていても、生成側と校正側で採用値が一致するとは限らない。

種類	旧校正平均：donor数	うち非正値	旧校正平均 [ADC/MIP]	clean donor数	clean平均 [ADC/MIP]	生成側平均 [ADC/MIP]
10 µm	3381	490	81.576853	2891	97.098354	97.340297
15 µm	1062	75	397.916766	987	428.913481	406.655638

cleanはMPV>0、NDF>0、MPV・Landau幅・Gaussian幅・χ²が有限、χ²≥0、従来のχ²/NDF上限を要求し、接続channelの正常値だけでsensor種別平均を作る。この方針をdata・MCの両方へ適用した。今回の定数監査では、問題の565 channelは非正MPV・非正NDF・負χ²が同時に現れているため、それぞれの条件による影響を独立に分解した結果ではない。

採用値の変更はdata 2,422、MC 2,422 channel。種類別には10 µm 2,149、15 µm 273。負の個別値の代替に加え、旧fallback channelの平均変更も含む。従来の本番校正コードや定数ファイルを更新したのではなく、この解析内で比較用に再計算した値である。

種類	試料	共通fallback channel数	旧MPV [MIP]	clean MPV [MIP]
10 µm	data	26	1.461384	1.227777
10 µm	mc	26	1.459326	1.226047
15 µm	data	8	0.767505	0.712039
15 µm	mc	8	0.844829	0.783775

10 µmでは約1.46から約1.23へ改善するが、なお1より高い。15 µmは元々1未満であり、分母が大きくなることでさらに下がる。不正値を除いた平均は校正の健全性のために必要だが、その平均が各channelの正しいMIP値になるとは限らない。

8. 同じhit集合で総energyはどの程度変わるか
============================================================

ここではMIP fit用の共通bin重み・HG限定・位置選択を追加せず、選択muonイベントの保存Calib_Hitにあるlayer 0〜29の全接続hitを固定した。HG/LG切替、threshold、crosstalk、既存前段のhit選択は変えていない。追加のbad-channel maskは適用していない。したがって、通常解析で別のmaskを使う場合や、hit選択から再実行する場合の変化率とは同じとは限らない。

試料	イベント数	旧総energy [MeV]	clean総energy [MeV]	差 [MeV]	相対変化
data	120,478	1,423,160.765	1,392,093.990	-31,066.775	-2.1829%
MC	805,256	8,803,587.405	11,340,615.926	+2,537,028.521	+28.8181%

試料	既存fallback平均の変更寄与 [MeV]	負の個別校正値の置換寄与 [MeV]
data	-113,685.812	+82,619.037
MC	-781,486.327	+3,318,514.848

MCの+28.82%は主として負energyの解消による。旧校正では正のADCを負のMIP定数で割るため、負energyが総和を打ち消していた。MCの負の個別校正channel群は旧合計−3,117,289.746 MeVからclean合計+201,225.102 MeVへ変わった。一方、既存fallbackの分母が増える効果はenergyを下げる。この二つの反対方向の寄与の合計が総変化となる。光量が28.82%増えた、あるいはそのままelectron energy scaleも28.82%変わる、という結果ではない。dataとMCの総energy値同士も、イベント数・選択分布が異なるため直接比較しない。

9. 図の意図と、数値を追うためのファイル
============================================================

図・ファイル	読み取る内容
closure_summary.png / .pdf / .svg	10／15 µm別のMC ADC MPV／入力MIP分布。旧校正側individual／fallbackで区分。1を中心としたずれとchannel間の幅を見る。
mc_channel_map.png / .pdf / .svg	横軸36×chip+channel（0〜209）、縦軸layer。青は1未満、白は1付近、赤は1超、灰色は有効fitなし。XY実空間mapではない。色は0.5〜1.5で飽和し、範囲外件数は図に表示。
data_mc_pairs.pdf	共通304 channelのclean data/MC比。全16ページ。20 channel／ページでlayer・chip・channelを表示。誤差棒はdata・MC独立を仮定したfit統計誤差だけ。
data_mc_pairs_page_001.png等	比較PDFのページ別PNG／SVG。data_mc_pairs.png/svgは1ページ目のみ。
channel_axis.tsv	全体連番210×layer+36×chip+channelとCellIDの対応。chip 5はchannel 0〜29、他chipは0〜35。
channel_closure.tsv / fits.tsv	channel別応答・採用区分・fit状態と不可理由。
group_summary.tsv / calibrated_group_summary.tsv	区分別中央値。全MC有効集合と共通data/MC集合を混同しない。
*_fallback_means.tsv / *_changed_channels.tsv	平均のdonor、旧／新採用値と変更channel。
energy_summary.tsv / channel_energy.tsv / event_energy.tsv	同一hit集合の総energy、channel寄与、イベント別比較。
inputs/ / STATUS.json / visualization.json	入力、ソースhash、run採用、復旧状態、設定、描画元の記録。

10. 現時点の限界と、次に確認すべきこと
============================================================

比較可能channelは304で、全6,300接続channelを代表する保証はない。解析対象4,506 channelのdata ADC20 fitでは、有効321、統計不足3,491。残りにもplateau範囲外・パラメータ境界などがある。有効なdata channelのうちMCと共通になる304 channelを比較図に使っている。共通binに残るhitはdata約17%、MC約18%である。無効fitは応答ゼロや検出器故障を意味しない。

今回の各中央値には入力MIP定数の誤差、温度補正、位置・角度分布差、threshold、fit範囲・モデルの系統誤差を付けていない。図のfit統計誤差だけで検出器全体の不一致の有意度を確定しない。また、現行ソースの再現と保存HGの整合は確認したが、過去のMC生成条件を完全に復元したわけではない。

次の優先事項1：校正の採用条件を統一する。MPV=−10、NDF=−10等を個別校正にも平均にも入れず、data・MCで同じ基準を使う。各fallback channelには、平均の代入だけで1 MIPが保証されないことを踏まえ、十分なmuon統計での個別校正や局所的な代替定数を検討する。今回のレポート作成では本番コードを変更していない。

次の優先事項2：同じtruth・channel集合でdigitizationの規格化を検証する。入力MIP値、gain、一次光電子数、crosstalk増幅後の光電子数、飽和補正の定義を明示し、crosstalkを有効に保ったまま一次光電子数を変更して再digitizationする。逆数1/R_ADCは初期試行値にできるが、MPVは分布形状や選択で変わるので確定補正として一括適用しない。各channelのADC MPV／入力MIPと校正後1 MIPへの復帰を再確認する。

次の優先事項3：残る系統依存を調べる。data温度22.8℃とMC20℃に対する現行温度補正、定数の基準温度、thresholdとfit範囲、共通bin幅を変えた安定性を確認する。MCのbeam angle／位置相関をdataに近づけた検証と、復旧不要なtrack入力の再作成も必要。これは今回未実施である。

今回見えたものは、(a)負の不正fit値がχ²/NDFだけの判定を通る具体的な校正問題、(b)正常平均に置き換えても残る絶対MIP規格化のずれ、(c)同じchannelで比較しても残る15 µm個別校正のdata/MC差である。元の疑問を調べるための診断は成立したが、最終的な光量定義の修正とその再生成検証はこれからである。

付録：再現用の結果参照先
============================================================

/home/murata_t/ScECAL_BeamTest/analysis/result/channel_response/mip_closure/all_muon_runs_20260930

解析コード：/home/murata_t/ScECAL_BeamTest/analysis/channel_response/mip_closure/。再描画：python3.12 visual_report.py RESULT_DIRECTORY。本レポートの数値は保存済みTSVから集計し、読み込んだ入力のSHA-256は [sources.json](/home/murata_t/scecal_muon_summary_20261001/detailed_summary_ja_20261001.sources.json) に記録した。

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](../../analysis/result/channel_response/mip_closure/all_muon_runs_20260930/)。

[HTML版](/home/murata_t/scecal_muon_summary_20261001/detailed_summary_ja_20261001.html)。
