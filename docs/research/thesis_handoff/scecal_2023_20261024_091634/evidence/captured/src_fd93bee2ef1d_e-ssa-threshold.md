# ADC・hit数・分解能の追加統合調査（2026-10-13）

**CT補正後のmuon応答を参考に光量・threshold・大信号応答を組み合わせて変えても、100 GeVの相対分散差の縮小は約3%にとどまった。一方、位置・runに対応する変動を除く診断では約55%縮む。残る差にもlayer間の相関が大きく寄与する。個別MIPを検証できないchannelが電子信号の大部分を担っているため、MIP不均一の影響を小さいと結論することはできない。ADC・hit数・resolutionが同時に一致したモデルは、まだ得られていない。**

追加結果は既存結果内の [joint_response](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/) に統一した。5/20/100 GeVのみ。dataは既存の低LG tail event除外後、元のCoG20・energy窓を固定した2733/1756/4999 event。MCのCT補正・組合せ試験は同じ2 truth files×1000 eventから、元の選別を固定した937/910/1516 eventを使用した。元のproduction e_ssa/thresholdと、CT入力を補正した診断MCを区別して保存している。

以下の「分散差の縮小率」は、R=RMS/meanとして
1−(R_data,変更後²−R_MC,変更後²)/(R_data,基準²−R_MC,基準²)。
Gaussian中心幅の改善率ではない。表のenergyは校正後の信号和であり、入射beam energyではない。

## 1. muonから制約できるMIP校正の範囲

data muonのunmatched / temperature_adc fitについて、fit有効、2/3/4σ下限でMPVが3%以内、MPV相対誤差≤5%、有効hit数≥100を要求した。MC fitの成立はdata channel採用条件にしていない。採用できたchannelだけMIPをdata muon MPVに置き換え、同じelectron ADC・hit・event集合を再校正した。

| beam | data信号の平均補完MIP channel由来割合 | muonで検証できたdata信号割合 | data RMS/mean：元→再校正 | 分散差の縮小：dataのみ変更 |
|---|---:|---:|---:|---:|
| 5 GeV | 82.44% | 11.50% | 12.447→12.417% | +0.79% |
| 20 GeV | 85.42% | 7.30% | 8.395→8.360% | +1.10% |
| 100 GeV | 85.33% | 3.92% | 5.019→5.023% | −0.22% |

100 GeVの10 µm・平均補完区分で検証できた信号は、**その区分の0.331%**のみ。10 µm平均補完channelのdata fitは1340 channelがmpv_near_plateau_boundaryで、3σ fitがokの288 channelから安定性・誤差条件を通るのは15 channelだった。単なる処理失敗ではなく、threshold近傍・分布切断によってMPVを確定しにくい。

**82–85%が誤校正と判明したわけではない。大きな信号割合が独立に未検証である、という意味。** 少数の採用channelの倍率を、検証できない全channelへ適用する根拠はない。data/MC双方を同じ係数で再校正する診断も行ったが、MCの光量入力は変更していないため、自己整合した再digitizationではない。

3試料すべてで安定した個別MIP channelの比較：

| SiPM | channel数 | MC/data MPV：元→CT入力補正後 | Gaussian σ/MPV：data / 元MC / CT補正MC |
|---|---:|---:|---:|
| 10 µm | 30 | 1.028→0.913 | 0.274 / 0.335 / 0.373 |
| 15 µm | 609 | 1.126→0.950 | 0.227 / 0.191 / 0.213 |

10 µmではmuonの相対Gaussian幅が既にMCで大きい。electron total energyの幅がMCで小さいことを理由に、一律に独立noiseを増やす修正は両方の分布を整合させない。Gaussian convolution幅は電子回路noiseそのものではなく、光電子統計・edep・通過長・fit相関も含む。

根拠：[coverage.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/muon_calibration/coverage.tsv)、[impact.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/muon_calibration/impact.tsv)、[muon_group_constraints.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/muon_calibration/muon_group_constraints.tsv)、[100 GeV再校正図](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/muon_calibration/comparison/20mm/100GeV/muon_calibration_impact.png)。

## 2. 光量・threshold・大信号応答を同時に変更した実digitizer試験

CT入力補正済みMCを基準に、一次光量を10 µm×1.10・15 µm×1.05にした。これは上記MPV差を参考にした候補で、変更後muonのMPV・幅・取得率を完全に検証した倍率ではない。CT確率・1 p.e. gain・追加noise・校正MIPは保持し、20℃固定で比較した。

| beam | data平均信号和 [GeV] | MC平均：CT基準→光量変更 | data平均hit数 | MC平均hit数：CT基準→光量変更 | data R | MC R：CT基準→光量変更 |
|---|---:|---:|---:|---:|---:|---:|
| 5 | 0.14257 | 0.16048→0.17602 | 109.98 | 113.47→115.69 | 12.447% | 7.814→7.755% |
| 20 | 0.72766 | 0.63973→0.70232 | 262.25 | 276.10→282.41 | 8.395% | 4.211→4.184% |
| 100 | 3.71219 | 3.06943→3.35857 | 622.10 | 687.64→701.93 | 5.019% | 2.640→2.636% |

20/100 GeVの平均energyは改善するが、5 GeVは逆方向で、hit過剰も増える。一つの光量倍率で全energyの不一致を解消することはできない。

100 GeVではさらにthresholdをHG ADCで+10、有効10 µm pixel数を10000→19600にする診断を組み合わせた。有効pixel数は回復・大信号非線形性の感度を見る代用であり、物理的pixel数を変更する推奨ではない。

| 条件 | 平均信号和 [GeV] | 平均hit数 | RMS/mean | Gaussian σ/mean |
|---|---:|---:|---:|---:|
| data | 3.71219 | 622.10 | 5.019% | 4.188% |
| CT補正MC基準 | 3.06943 | 687.64 | 2.640% | 2.382% |
| 光量変更 | 3.35857 | 701.93 | 2.636% | 2.398% |
| 光量＋threshold | 3.35459 | 682.28 | 2.631% | 2.370% |
| 光量＋大信号応答 | 3.45619 | 702.05 | 2.766% | 2.618% |
| 光量＋大信号応答＋threshold | 3.45222 | 682.50 | 2.738% | 2.531% |

光量＋大信号応答の分散差縮小は3.74%［16–84%区間3.12–4.42%］、3項目同時では2.91%［2.40–3.46%］。400回のpaired block bootstrapによる区間で、MC 2 filesを超える系統差は含まない。

3項目同時でも平均energyはdataより7.00%低く、hit数は9.71%多い。中央hit数はdata191.13 / MC191.53で近いが、周辺は430.97 / 490.97でMCが13.92%多い。平均Nhitを合わせるには、周辺のchannel別取得効率・thresholdの立ち上がり・shower広がりを別に調べる必要がある。

共通10 µm個別MIP channel・同じgainで生存した同一hitでは、光量変更MC/CT基準ADC中央値は低HG 1.099–1.102。大信号応答変更を足しても低HGは1.101–1.103でほぼ同じだが、truth edep30–100 MeVのLGは1.093→1.138、≥100 MeVは1.083→1.177となる。大信号応答は低信号ADCとある程度独立に調整できる。ただしこれはMC内の同じhit比較であり、dataに対する飽和モデルの正しさを検証したものではない。

根拠：[joint_summary.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/digitization/joint_summary.tsv)、[joint_impact.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/digitization/joint_impact.tsv)、[ADC比較](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/digitization/adc_controls.tsv)、[100 GeV分布](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/digitization/100GeV/joint_response.png)。

## 3. 追加：位置・run依存を差し引いた後のlayer相関

各layerのenergyをshower CoGの2次式とdata run番号／MC file番号で回帰し、100-entry blockによる5foldの未使用eventへ適用した。layerごとの予測依存分を加法的に差し引き、各layer平均を元の値へ戻した。したがって平均は保存されるが、layer残差は負にもなり得る。この診断量を校正energyとして使うことはできない。従来のenergyを応答関数で割る位置補正とは異なる、相関の分解である。

| beam | data R：元→位置・run依存除去 | CT補正MC R：元→位置・file依存除去 | 分散差縮小 | 残った差のうちlayer間共分散項 |
|---|---:|---:|---:|---:|
| 5 | 12.447→12.021% | 7.814→7.657% | 8.53%［6.59–10.17%］ | 85.8% |
| 20 | 8.395→6.975% | 4.211→3.908% | 36.72%［33.23–40.67%］ | 74.0% |
| 100 | 5.019→3.607% | 2.640→2.185% | 54.81%［50.94–59.01%］ | 69.6% |

3項目同時変更MCでも100 GeVは2.738→2.304%、dataは3.607%。そのモデルに対する分散差縮小は56.47%、残差の78.3%はlayer間共分散項。ADC応答の変更と位置・run依存の診断を合わせても差は残る。

共分散項の割合は Var(ΣE_layer)=ΣVar(E_layer)+Σ(layer≠layer')Cov の恒等式による分解で、原因の寄与率ではない。正の相関が強い場合だけでなく、shower前後の負の相関による打ち消しが弱い場合もtotal幅は広がる。100 GeVで前方0–9層と後方20–29層の相関は、位置・run依存除去後data −0.747 / CT補正MC −0.810。独立noise、共通gain変動、縦方向shower形状、校正の空間不均一を区別する必要がある。

100 GeVの依存除去後の68%中心幅はdata2.250% / CT補正MC1.928%で、RMS差より近い。分布中心とtailの両方を監視することが重要。残差を全event共通の独立乗算ゆらぎだけで埋めると仮定した逆算値は5/20/100 GeVで9.24/5.77/2.87%。これは測定値でも追加noiseの推奨値でもない。

CoGはshower自身から求めるため、この54.81%を真の入射位置不均一やleakの寄与と同定できない。区間は300回のfile内paired block bootstrapで、回帰係数は固定。係数再fit・run間系統差は含まない。

根拠：[summary.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/layer_covariance/summary.tsv)、[impact.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/layer_covariance/impact.tsv)、[100 GeV layer共分散図](../../investigations/result/e_ssa_threshold_review_20260930/20261013/joint_response/layer_covariance/comparison/mc_joint/20mm/100GeV/layer_covariance.png)。図は上段が元の分布、下段が位置・run/file依存を差し引いた残差、左からdata・MC・data−MC。色はCov(E_i,E_j)/〈total E〉²。data/MCの色範囲を共通にし、差の図も上下で共通にした。

## 改善に向けた優先順位

1. **中央10 µmのMIP・光量基準を検証する。** muon trackから対象stripを除いて通過位置・通過長を予測し、取得効率とMPVを別に求める。現行thresholdでMPVが隠れるchannelでは低threshold muonまたは適切なLED測定が必要。採用MIP、測定温度、個別値/平均補完を照合する。未検証channelへ少数の倍率を拡張しない。
2. **実測gain・CT・PDEとmuon ADC分布を同時に制約する。** 1 p.e. charge/幅、CT多重度を独立に使い、固定biasでの温度基準を統一する。光量は実digitizer後のmuon MPVで決め、幅・取得率も確認する。現行3√PE追加ADC幅の根拠を検証する。温度モデルのみの分散差縮小は既調査で約2–3%。
3. **周辺hit取得率と大信号非線形性を分けて測定する。** thresholdの値だけでなく立ち上がり曲線とchannel/chipの同時欠落を確認する。SiPM大信号はLED強度・pulse幅とADC非線形性から回復・飽和を制約する。中央hit数が近いことは中央energy応答の一致を意味しない。
4. **event相関とshowerの縦・横形状を確認する。** 独立したmuon/LEDによるstrip内位置応答map、layer/chipで共有される時間変動、run別beam条件、前後layerのenergy交換、後方/側方信号を確認する。光学応答不均一、校正、beam spread、leakをこの診断だけで一意に分離できない。truthのgood-channel edep和でもMC幅が狭いという既結果は、SiPM飽和除去だけでは足りない参考であり、leak率そのものではない。

最終確認は、同じ選別・mask・温度条件でADC MPVと幅、threshold近傍効率、中央/周辺Nhit、平均energy、RMS/68%幅/Gaussian中心幅、layer共分散を同時に行う。上記impactには相互依存があり、百分率を足し合わせない。

再現コードは [調査README](../../investigations/e_ssa_threshold_review_20260930/README.md) を参照。追加の主なC++はJointResponseReview.cc、MuonCalibrationImpactReview.cc、LayerCovarianceReview.cc。新旧入力の対応、layer和・全energy分散の一致、入力checksum不変を確認した。科学的ROOT/TSVと使用ソースを保存し、一時実行ファイルは削除した。productionコード・定数は変更していない。

---

# 温度依存CT・gainモデルのelectron分解能試算（2026-10-13）

**温度モデルと20℃規格化を整合させる変更は、100 GeV MCのRMS/meanを旧CT修正MCの2.640%から2.647–2.665%へ少し動かすが、dataとの差をほとんど埋めない。MCだけの相対分散差の解消は0.22–0.75%、dataにも同じ逆応答係数を適用する診断を合わせても1.75–2.73%。温度モデルはADC・光電子数・hit取得率の整合性を直すために必要だが、現時点のelectron resolution差の主因を解消する変更とは言えない。**

結果は既存20261013配下の [一次近似](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/) と [100 GeV実digitizer試験](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_pilot/) に保存した。全energyの再生成は行わず、5/20/100 GeVの一次近似・data校正診断と、100 GeVの2,000truth eventだけを使用した。

## 仮定と試験条件

ユーザー確認によりbeam testのbiasは固定。基準温度20℃、Uov(20℃)=10 µmで4.5 V、15 µmで4.0 V、降伏電圧の温度係数60 mV/℃を仮定した。gain比は1−α(T−20)、α=1.6/135、3.5/230。CT確率はp0[Uov(T)/Uov(20)]²、p0=0.10/0.15。PDEは一定、またはUov比に比例する2条件。実機の測定曲線をfitした最終モデルではなく、感度を見る仮定であり、この2条件が全ての物理的な可能性を挟むとは限らない。

使用中LEDファイルのchannelごとの測定温度を読み、gain20=LEDgain/[1−α(T_LED−20)]に換算した。光子量は20℃で(1−p0)×既存MIP/gain20/PDE20から一度決め、その後は温度を変えても固定。現行の10,000/4,489 pixel、CT再帰処理、3√Nfiredの追加ADC幅、threshold、共通MIP定数を保持した。前のmuon全run fitで残ったMPVのずれを合わせる追加規格化は行っていない。

## 5/20/100 GeV：観測温度幅を使った一次近似

同じCT修正済みMC event・hitにdataのevent平均hit温度を独立に割り当てる、全hitの温度が同じという相関の強い条件で計算した。下表は平均応答の温度依存を乗算した段階で、整合する20℃補正をまだ戻さない場合。threshold・光電子統計・CT分散・pixel飽和の再生成は含まない。

| energy | data RMS/mean | CT修正MC | gain＋CT、PDE一定～線形 | 相対分散差の解消 |
|---|---:|---:|---:|---:|
| 5 GeV | 12.447% | 7.814% | 7.841–7.854% | 0.44–0.66% |
| 20 GeV | 8.395% | 4.211% | 4.239–4.315% | 0.46–1.68% |
| 100 GeV | 5.019% | 2.640% | 2.647–2.661% | 0.23–0.62% |

選別後のdata平均hit温度の平均/event間RMSは28.976/0.099℃、23.512/0.320℃、22.992/0.104℃。run内だけのRMSは0.047/0.049/0.033℃。PDE線形条件で純粋なevent間温度変動によるRMS/meanの増加は0.0063/0.1043/0.0174 percentage points。上表には平均温度が20℃と異なることで10/15 µmの相対寄与が変わる効果も含まれる。

全ての分散差をこのモデルの共通温度変動で作るには、温度RMS約3.05/2.46/1.46℃が必要で、観測値より大幅に大きい。ただし平均hit温度は空間的な温度場を全部表さず、厳密な上限ではない。モデル化した平均応答を全て20℃へ戻すと、この乗算近似の温度変動は数値的に相殺される。

[一次近似の全条件](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/summary.tsv)、[温度幅](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/temperature_distribution.tsv)、[必要温度幅の換算](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/required_temperature_spread.tsv)。

## 100 GeV：実digitizerによる確認

旧MCの同じtruth file0/99から各1,000event、元のCoG20・energy窓を固定した1,516eventを比較した。data選別4,999eventから固定seedでshuffleし、2,000eventの30layer×16sensor温度matrixを非復元抽出。Run105/106/107は636/660/704event、layer・sensor間同時相関を保持し、MC truthとは独立に割り当てた。空のlayer温度349/60,000は現行と同じ20℃fallback。旧・新raw/calib ROOTのclosure、event対応を確認した。

| MC条件 | 平均signal和 [GeV] | 平均hit数 | RMS/mean | Gaussian中心幅 |
|---|---:|---:|---:|---:|
| 旧CT修正MC：LED測定gainを20℃として使用 | 3.06943 | 687.64 | 2.640% | 2.382% |
| gainを20℃へ換算、温度20℃固定 | 3.07733 | 686.96 | 2.679% | 2.346% |
| 実温度、PDE一定、全応答逆補正 | 3.08010 | 679.71 | 2.647% | 2.210% |
| 実温度、PDE線形、全応答逆補正 | 3.08495 | 672.12 | 2.665% | 2.379% |

Gaussian fitの誤差は約0.10–0.12 percentage pointで、今回の小規模sampleから中心幅の有意な改善・悪化は判定できない。68%幅も旧MC2.445%、修正後2.450–2.475%でほぼ同じ。

全応答逆補正はHG相当信号に1/[gain比×PDE比×(1−p0)/(1−p(T))]を掛ける。LGではpedestal差引きLG×HLratio＋HGoffsetを先に再構成し、全体を温度補正する。低占有平均の近似なので飽和hitへの完全な補正ではない。

新20℃固定基準から実温度＋全応答逆補正へのRMS変化はPDE一定−0.0315 point [−0.0483,−0.0140]、PDE線形−0.0133 point [−0.0316,+0.0072]。区間は300回のpaired 100-entry-block bootstrapの16–84%で、file/モデル系統差を含まない。**温度依存を追加する部分だけでは、MC幅が必ず大きくなるわけではない。** 旧MCからの小さい全体変化にはgain参照温度の変更も含む。独立なdigitizer乱数系列なので、小さい差に乱数揺らぎも含まれる。

一方、実温度を入れても現行の補正をそのまま使うと、RMSは2.704–2.762%へ増えるが平均signalは2.841–2.731 GeVへ低下する。温度応答と整合しない校正を重ねた結果であり、resolutionを近づける変更として採用しない。

[実digitizerの基準値](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_pilot/summary.tsv)、[impactと区間](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_pilot/impact.tsv)、[比較PNG](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_pilot/temperature_resolution.png)。

## data校正を同時に変える診断

dataもrawから同じ逆応答係数で再換算した。旧event・非負good hit・HG/LG選択は固定。MIP定数の実際の取得温度・補正履歴と実測gain勾配が未確定なので、以下は校正感度であって確定修正ではない。

| data RMS/mean | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| 現行 | 12.447% | 8.395% | 5.019% |
| gain逆補正 | 12.388% | 8.262% | 4.994% |
| gain＋CT逆補正、PDE一定 | 12.383% | 8.249% | 4.991% |
| gain＋CT＋PDE線形の逆補正 | 12.373% | 8.214% | 4.983% |

100 GeVでdata・MCを同じ条件に変更し、G=R_data²−R_MC²、impact=1−G_after/G_beforeとすると旧CT修正MC基準の解消率はPDE一定1.75%、PDE線形2.73%。Gaussian中心幅の改善率とは異なる。平均signalのdata側変化は100 GeVで+8.4/+12.9%、5 GeVでは+29/+47%と大きいため、reference温度・係数の正しさを確認せず本番校正へ入れない。

[同じhitのdata再換算](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/data_calibration/summary.tsv)、[paired区間](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/data_calibration/impact.tsv)、[closure](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/temperature_estimate/data_calibration/validation.tsv)。逆補正はLG切片を係数内へ入れるため、温度係数の方向だけでなく切片位置の変更も含む。

C++/shellは [既存調査ディレクトリ](../../investigations/e_ssa_threshold_review_20260930/) の `TemperatureResolutionEstimate`、`DataTemperatureImpact`、`TemperatureDigitizationReview` と対応runnerに保存。入力一覧・温度bank・使用source/binary・ROOT・TSV・PNG・実行記録を保持し、production入力・source・校正定数は変更していない。

---

# ADC・hit数・resolutionの統合再調査（2026-10-13）

**CT入力規格化の修正はmuon ADCを改善するが、electronのresolution差への効果は小さい。残る差は、中央の大信号応答、周辺のhit分布、event単位のlayer間相関を同時に調べる必要がある。現段階でADC・hit数・resolutionを同時に再現したモデルは確立していない。**

結果は [20261013の一式](../../investigations/result/e_ssa_threshold_review_20260930/20261013/) に統一した。`integrated/` は最新electron入力の再読取り、`muon_ct/` はCT修正前後のmuon代表run、`raw_adc/` はchannel別ADCと占有率、`event_components/` はhit数・位置・layer相関、`digitization/` は実際のproduction digitizerによるelectron再生成の小規模比較である。

## 比較条件と基準値

5・20・100 GeVだけを調査した。dataは2023 PS/SPS、electron MCは `Result_MC/{decode,calib}/e-/e_ssa/threshold/`。今回確認したelectronファイルは10月11日生成のCT入力修正前であり、10月13日の現行digitizerソースと混同しない。decode・calibを再読取りし、event ID、full CellID、校正式、更新時刻・sizeの不変を検証した。基準値は10月12日の比較と一致する。

主比較は既存bad-channel maskの非負good-hit energy和、layer9/10のCoG±20 mm、全接続channelの正energy和が基準MC中央値の0.5–1.5倍というenergy窓を使用した。dataだけproduction低LG tail-rich eventを除外し、MCには古いtail判定を適用しない。energy窓なしの結果も保存し、特に低energy側の裾と中心幅を区別した。

| 入射energy | 採用event data / MC | Gaussian中心幅 data / MC | RMS/mean data / MC | 平均hit数 data / MC | event平均E/Nhit [MeV] data / MC |
|---|---:|---:|---:|---:|---:|
| 5 GeV | 2733 / 2825 | 8.63 / 7.43% | 12.45 / 7.62% | 110.0 / 116.3 | 1.308 / 1.550 |
| 20 GeV | 1756 / 2740 | 7.17 / 4.12% | 8.39 / 4.31% | 262.3 / 284.0 | 2.798 / 2.518 |
| 100 GeV | 4999 / 4541 | 4.19 / 2.40% | 5.02 / 2.63% | 622.1 / 705.7 | 5.991 / 4.821 |

Gaussianは共通の二段階fitによる診断値で、本番resolution fitそのものとは区別する。event平均E/Nhitは「各eventのE/Nhitを平均した値」であり、平均E÷平均Nhitではない。energyはcalibrated signal和で、入射energyではない。

[基準値](../../investigations/result/e_ssa_threshold_review_20260930/20261013/integrated/response_comparison.tsv)、[hit数・分位点](../../investigations/result/e_ssa_threshold_review_20260930/20261013/event_components/distribution_summary.tsv)、[ADC・event量](../../investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/comparison/20mm/100GeV/event_moments.tsv)。

## CT修正：ADC規格化とresolutionへの影響を分ける

現行修正はcrosstalkの実装を残し、一次光量に10 µmで0.90、15 µmで0.85を掛けるもの。校正MIP定数はdata/MC共通のままである。`threshold_MPVwoCT`という名前はCT自体を無効化した意味ではない。

muon Run44を各20,000event抽出し、旧MC・新MCで同じtruthを使用した。位置・角度を三sampleでそろえた15 µm個別定数群のADC/MIP中央値はdata1.113、旧MC1.343、新MC1.158で、dataに対する中央値の過大は約20.6%から4.0%へ縮小した。中央68%幅/中央値はdata0.432、旧MC0.412、新MC0.424だった。ただしmatched dataは約469 weighted hitであり、channel別の三sample共通MPV fitに必要な統計が不足している。全runのmuon fit結果としては扱わない。

671,702 truth rowのCellID、edep、primary chordを旧/新で照合し、完全一致した。同じtruthに対する検出率は10 µm個別定数群で4.773±0.096 percentage point、15 µmで1.627±0.070 point低下した。誤差は同一生成eventを対応させたevent-cluster SE。trackから予測したstripのhit取得率には位置推定誤差も含まれるため、hardware単独の効率には読み替えない。

[muon ADC比較](../../investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct/adc_distribution.tsv)、[同じtruthの検出率](../../investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct/paired_truth_detection.tsv)、[muon PNG](../../investigations/result/e_ssa_threshold_review_20260930/20261013/muon_ct/matched_three_samples_adc.png)。

さらに実際のproduction digitizerとcalibratorのbinary・ソース・定数を保存し、electronの同じtruth各energy2ファイル、計2,000eventを再digitizationした。元のevent選別を固定し、再選別の違いも別表に保存した。既存production MCを上書きせず、診断ROOTだけをこの結果配下へ保存した。

| CT入力修正：同じ旧選別event | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| MC平均energyの変化 | −10.58% | −10.19% | −9.70% |
| MC平均hit数の変化 | −2.32% | −2.65% | −2.48% |
| MC RMS/meanの変化 | +0.128 point | −0.020 point | +0.045 point |
| data/MC相対分散差の減少 | 2.07% | −0.32% | 1.27% |
| 新MC平均energy [GeV] | 0.16048 | 0.63973 | 3.06943 |
| 新MC平均hit数 | 113.47 | 276.10 | 687.64 |

pilotの採用eventは937/910/1516で、基準6ファイルのMC全sampleと同一ではない。RMS変化のpaired 100-entry-block bootstrapの16–84%区間は、5 GeV [+0.069,+0.209] point、20 GeV [−0.057,+0.010] point、100 GeV [+0.030,+0.058] point。100 GeVの新MCはdataより平均energyが約17.3%低く、hit数は約10.5%多い。CT規格化だけでは両方を合わせられない。

低truth edepの検出率も変わる。100 GeV、10 µmのtruth 0.1–0.305 MeVでは旧47.51%→新40.91%、0.305–0.6 MeVでは93.63%→90.43%。0.6–1.2 MeVでは99.90%→99.79%でほぼplateauである。低信号で「両方に残ったhit」だけ比べるとADC比が1に近づくことがあり、これはthreshold選別の効果であって入力補正が消えた意味ではない。

[CT修正impactとbootstrap](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/impact.tsv)、[再生成の基準値](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/summary.tsv)、[truth帯別の取得率](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/truth_hit_acceptance.tsv)。

## 100 GeVの実digitization追加試験：飽和とthresholdを別々に変更

CT入力修正後を基準に、10 µmの有効pixel数を10,000→14,400/19,600にしたケース、pixel数を据え置いてthresholdを+5/+10 HG ADCにしたケースを生成した。15 µmのpixel数4,489、CT確率、一次光量の(1−p)、共通MIP校正、threshold幅は固定。 pixel格子を広げる変更にはCTの境界損失の変化も含まれる。これは回復・非線形性への感度を見る診断であり、実際のmicrocell数が14,400/19,600だと推定したものではない。thresholdは保存hitの後段cutではなく、digitizationのhit取得判定で変更した。

100 GeVの同じtruth 2ファイル×1,000event、旧MC選別を固定した1,516eventで比較した。

| 条件 | 平均energy [GeV] | 平均hit数 | RMS/mean | CT基準からの相対分散差の減少 |
|---|---:|---:|---:|---:|
| data：全調査sample、tail除外 | 3.71219 | 622.10 | 5.019% | — |
| CT入力修正後のMC基準 | 3.06943 | 687.64 | 2.640% | 0% |
| 有効pixel数14,400 | 3.11941 | 687.49 | 2.680% | 1.17% |
| 有効pixel数19,600 | 3.15103 | 687.38 | 2.749% | 3.24% |
| threshold +5 ADC | 3.06821 | 676.91 | 2.632% | −0.21% |
| threshold +10 ADC | 3.06581 | 666.63 | 2.650% | 0.30% |

pixel19,600で平均energyは2.658%増加するが、dataより約15.1%低いまま。RMS変化は+0.109 point、paired block bootstrapの16–84%区間は[+0.090,+0.127] point、分散差の減少は3.24% [2.54,3.95]。CT基準からの追加効果と、旧MCからの効果4.47%を混同しない。

threshold+10でhit数は3.056%減るが、dataより約7.2%多い。中央hit数はCT基準191.67→190.74（data191.13）とほぼ一定で、周辺は495.97→475.89（data430.97）と減る。低信号の周辺hitを主に削るが、平均energy変化は−0.118%にとどまる。RMS変化+0.010 pointの区間は[−0.003,+0.025] pointで、全幅への影響は小さい。Gaussian幅は2.382→2.184%、68%幅も2.445→2.366%へ狭くなるので、RMSだけ合わせることも不適切である。

同じeventの**全good-channel Geant4 edepをADC threshold・飽和なしで線形に足したtruth参照**は平均3.27635 GeV、RMS/mean2.857%、68%幅2.660%。これはdetectorのgain・MIP規格化をdataで測り直した値ではないが、飽和を取り除いた線形truthのevent幅にもdataとの差が残ることを示す。 truthの幅は一般的な上限ではない。各ケースは独立なRNG系列で生成しており、2ファイル内のblock bootstrapはseedやfile間の系統差を評価していない。`edep≤1.2 MeV`だけのsubtotalは別に保存し、total resolutionへ読み替えない。以上の4例だけで飽和・threshold全般を否定はできないが、どちらか1係数を調整してADC・Nhit・resolutionを同時に再現できる根拠は得られなかった。

[条件別基準値](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/model_summary.tsv)、[CT基準・旧MCの両方に対するimpact](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/model_impact.tsv)、[truth参照](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/truth_reference.tsv)、[比較PNG](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/model_scan.png)。

全6条件に同じphysical hitが残り、低信号では全条件HG、高信号では全条件LGかつLG<3000となる共通hitだけの追加ADC比較も行った。10 µm個別定数channelでpixel19,600/CT基準の同一hit HG換算ADC中央値比はtruth0.6–1.2 MeVで1.0000（74,340 hit）、1.2–3 MeVで1.0017（46,942 hit）。LGの30–100 MeVでは1.0375（3,840 hit）、100 MeV以上では1.0775（624 hit）だった。つまり低信号規格化を約0–0.2%以内で保ちつつ、高信号の応答を約3.75–7.75%動かしている。ただし全条件に残るhitに限定した条件付き比較であり、dataに対する飽和不足を直接測った値ではない。ADC表のq16/q84はhit分布の分位点で、推定値の信頼区間ではない。[共通hitのADC比較](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/adc_controls.tsv)、[共通hitの割合](../../investigations/result/e_ssa_threshold_review_20260930/20261013/digitization/model_scan/adc_control_coverage.tsv)。

## hit数と大信号応答：空間的な差が残る

100 GeVの中央stripの平均hit数はdata191.13、MC192.59で近い。一方、周辺は430.97/513.11でdataが約16.0%少ない。共通10 µm channelを等重みで比べた占有率比は中央0.985、周辺0.782。これは共通channelに限定した比の平均であり、全channelのhit和の比0.840とは母集団・重みが異なる。

100 GeVではevent hit数のRMS/meanがdata7.06%、MC3.60%、event平均hit-energyのRMS/meanが7.52%/3.47%。meanやNhitだけでなく、それらのevent変動と相関にも差がある。中央共通10 µm channelのenergy中央値比はdata/MC1.139、非飽和HG20℃の中央値比は1.021で、主な中央差はHGのMIP-scaleだけでは説明できない。ただしelectronのADC分布はshower edep・channel混合も含み、MIP gainの直接測定ではない。

[100 GeV channel比較](../../investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/comparison/20mm/100GeV/common_channels.tsv)、[空間占有率とHG4000 map](../../investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/comparison/20mm/100GeV/occupancy_hg4000_maps.png)、[event分布](../../investigations/result/e_ssa_threshold_review_20260930/20261013/event_components/comparison/cog20_energy_window/100GeV/event_distributions.png)。

## 校正・低信号cut・HG上限だけでは差を埋めない

同じhit集合をSiPM群平均MIPで換算する診断の相対分散差の減少は5/20/100 GeVで3.28/4.69/3.39%。100 GeVのacquisition-block bootstrap 95%区間は2.27–4.84%。個別MIPを平均値に置き換えることを本番修正として勧める結果ではなく、channel換算値の混合だけでは幅の大部分を説明できないことを示す。

保存済みhitに0.3 MeV cutを追加しても分散差は4.69/1.79/0.17%悪化した。この結果は取得時のthresholdモデルの正しさを証明しない。取得前に失われたhit、trigger turn-on、event/chipの取得相関は別に検証する。

LG≥3000はdata/MCとも見つからず、HG4000を早期に採用した異常subsetのenergy寄与はMCで0.0208/0.0170/0.0344%と小さい。HG4000は通常LGへの切替表示なので、それらを全削除するとMCの正規の信号を5/20/100 GeVで約20/51/74%失う。図や幅が近づいても物理的な補正ではない。これはSiPMのpixel飽和の検証とは別の問題である。

低LG tail-rich除外後の共通channelの残る低側5σ residual率はdata約0.055/0.038/0.053%、MCほぼ0であり、以前の0.66%を現状の主因として再利用しない。raw/calib/選別event identityとenergy式のclosureは通過した。

[同じhitの診断impact](../../investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/comparison/20mm/100GeV/impact.tsv)、[clipping・占有率](../../investigations/result/e_ssa_threshold_review_20260930/20261013/raw_adc/comparison/20mm/100GeV/occupancy_clipping.tsv)。

温度補正を外すだけの診断でも、相対分散差の減少は5/20/100 GeVで1.15/2.42/0.72%と小さい。一方、dataの平均energyは11.92/4.24/3.62%増加するので、温度は平均scaleの検証から除外しない。前述の20 GeVの温度回帰は別runへの外挿に失敗しており、幅を合わせる目的で追加係数を確定しない。

後方5 layerの信号割合>5%のevent率はenergy窓なしでdata/MCが5 GeV9.66/1.42%、20 GeV2.58/3.25%、100 GeV9.23/19.91%。これは後方showerのproxyであって、escaped energyを直接測ったleak率ではない。低energy側のevent tail、run差、shower depthの条件も分けて確認する。[後方信号率](../../investigations/result/e_ssa_threshold_review_20260930/20261013/integrated/leak_proxy_rates.tsv)、[温度等の診断impact](../../investigations/result/e_ssa_threshold_review_20260930/20261013/integrated/impact.tsv)。

## resolutionの差：layer間相関と位置に関連した変動

R=RMS/mean、G=R_data²−R_MC²としてimpactを1−G_after/G_beforeで定義した。これはGaussian resolutionそのものの改善率ではない。層energyの分散・共分散の和はevent分散と数値的に一致し、全分散差のうち層間共分散差の寄与は5/20/100 GeVで約82.2/75.4/78.3%だった。100 GeVではMCのlayer間反相関による打ち消しがdataより強い。これは分散の会計であり、共通電子noise・beam spread・位置応答・shower/leakのどれかを確定した比率ではない。

| eventを分けて学習・評価した応答診断 | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| Nhitで平均応答を説明 | 47.6% | −1.3% | 5.4% |
| shower CoG位置で平均応答を説明 | 8.3% | 30.1% | 55.0% |
| 位置＋Nhit＋後方信号proxy | 76.6% | 55.6% | 73.5% |

上表は相対分散差の減少率。run/file単位での外部評価も概ね再現するが、CoG・Nhit・後方比は同じshowerの信号から得る量なので因果寄与とはしない。energy窓なしの5/20 GeVでは一部モデルが極端な補正を出し、別runの20 GeV温度モデルは外挿に失敗した。これらを本番補正に採用しない。

100 GeVでMCのNhit分布をdataへ直接重み付けすると、MC有効統計が約55 eventまで低下し、最大重みは約122となる。見かけの幅一致は稀なMC eventに依存するため採用できない。位置の分布をそろえるだけでは分散差は改善せず、位置に対する応答の違い・showerとの相関を区別する必要がある。

同じhit集合へ独立な一定noiseを足すだけでRMS差を埋める換算では、1 hit当たり約1.64/3.05/5.47 MeVもの追加σが必要になる。これはthreshold・hit集合固定の概算で、信号依存noise全般を否定しない。一律pedestal noiseの増加を選ぶ根拠にはしない。共通event gain変動なら必要σは約9.81/7.20/4.27%となるが、幅を合わせるために逆算したbenchmarkであり、物理的な測定値ではない。

[共分散・必要変動規模](../../investigations/result/e_ssa_threshold_review_20260930/20261013/event_components/variance_components.tsv)、[交差検証impact](../../investigations/result/e_ssa_threshold_review_20260930/20261013/event_components/crossvalidated_impact.tsv)、[重み付けの有効統計](../../investigations/result/e_ssa_threshold_review_20260930/20261013/event_components/reweighting_impact.tsv)。

## 改善の順序と採用条件

1. **muonのADC再現を基準に低信号規格化を固定する。** CT修正後の全run解析でMPV・分位点幅・tail・位置/温度依存を同時に確認する。個別定数と平均補完channelを分け、fit不安定・統計不足を確定係数へ使わない。
2. **大信号応答を独立に検証する。** ADCのHG/LG接続とSiPMの飽和を分け、LED等の大信号応答から有効pixel数・回復時間を制約する。muon規格化を再変更してelectron平均だけを合わせない。
3. **周辺の占有率・hit数をADC条件と同時に再現する。** threshold turn-on、無効channel、chip単位の取得、strip内光収集、beam/showerの横方向分布を調べる。全hitへの独立間引きや稀なeventへの大重みでNhitだけを一致させない。
4. **event・layer間の相関を再現する。** muonの複数層の応答、run内温度・pedestal/LEDの共通変動、beam条件を使い、独立に測った変動をMCへ入れる。shower由来CoG補正は診断として保持し、因果確認の代わりにしない。
5. **最後に共通のresolution解析で検証する。** 同じmask・event選別・tail処理を使い、ADC/hitの平均・幅・tail、layer/位置profile、layer共分散、Gaussian中心幅とRMS/68%幅を同じモデルで確認する。今回の試験は5/20/100 GeVのみで、全energyへの本番変更はしていない。

## 再現と保存

コードは既存の [調査ディレクトリ](../../investigations/e_ssa_threshold_review_20260930/) に置いた。`IntegratedReview.cc` はevent表へ温度・LG寄与・time/trigger・各layer energyを追加した。新しいC++/shellは `MuonsCTReview`、`RawADCReview`、`EventDistributionReview`、`DigitizationImpactReview`。CT electron再生成は `ct_digitization_study.sh`、100 GeVの4条件診断は `digitization_model_scan.sh` でジョブ投入した。科学的ROOT・TSV、PNG、実行記録、入力checksum、使用ソース・binaryは結果配下に保持した。入力変更・event対応・calibration closureのチェックに通過しており、検証用の一時実行ファイルだけを削除した。

---

# data/MC不一致の統合再調査（2026-10-12）

最新のdata（calib更新10月10日）とMC `e-/e_ssa/threshold`（decode/calib更新10月11日）を読み直した。対象は5・20・100 GeV。温度、MIP、HG/LG、位置依存、shower leak、SSA形状、beam/thresholdの別解析も監査した。

**残る中心幅の差で最も大きな診断効果は、位置に関連したevent応答変動。平均応答にはenergy依存の符号反転があり、大信号の非線形性と周辺stripの信号・hit数差を別々に調べる必要がある。後方信号の多いeventは主にenergy分布の裾に関係し、温度だけを外す効果は幅の差には小さい。**

[今回の結果一式](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/)、[平均応答・幅の比較](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/response_comparison.tsv)、[impact](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/impact.tsv)、[最新保存resolutionの読取り](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/production_resolution.tsv)。10月2日・4日・5日の数値は、以下の履歴に保持するが最新入力の数値として再利用しない。

## 条件と現在の差

MCは各energyの辞書順に広く6ファイル、各1000event。現在の利用可能MCは各100ファイル。dataは全runから各2000eventを10区間に分散して抽出し、5/20/100 GeVで8000/4000/6000event。runの全event数による重み付けは行わない。

主比較は既存bad-channel maskの非負good-hit energy和。layer9/10のgood-energy CoGがx/yとも±20mm、全接続channelの非負energy和がMCの符号付きenergy中央値の0.5–1.5倍を要求する。dataだけ現在のfull-sample低LG tail-rich判定を除外する。MCへのtail vetoは行わない。MC更新はtail解析より後なので、古いMC判定は再利用しない。dataの判定は入力更新時刻と全eventのrun/time/triggerを照合した。

| beam | 採用event data / MC | 平均energy data / MC [GeV] | data/MC−1 | RMS/mean data / MC | Gaussian中心幅 data / MC |
|---|---:|---:|---:|---:|---:|
| 5 | 2733 / 2825 | 0.142574 / 0.179247 | −20.46% | 12.45 / 7.62% | 8.63 / 7.43% |
| 20 | 1756 / 2740 | 0.727664 / 0.712938 | +2.07% | 8.39 / 4.31% | 7.17 / 4.12% |
| 100 | 4999 / 4541 | 3.712195 / 3.399409 | +9.20% | 5.02 / 2.63% | 4.19 / 2.40% |

energyは保存されたcalibrated signal和で、beam energyそのものではない。Gaussian幅は240bin、範囲[0.3,1.7]×各sample平均、最初に平均±3RMS、次にfit平均±1.5sigmaでfitした診断値。全84fitのstatusは0。選別・抽出・binningは本番resolutionとは異なる。

最新の全入力による保存済みresolution（CoG20、beforeEventCut）でもGaussian幅は5/20/100 GeVでdata8.84/6.66/4.14%、MC7.59/4.08/2.29%。同じ傾向を確認した。採用run統計とenergy裾の違いにより、sampleのRMSと本番RMSを混同しない。

## 1. 位置に関連した応答変動：中心幅への大きな効果

各sampleのenergyを、layer9 CoG xとlayer10 CoG yの二次曲面、run別定数で説明し、予測相対応答で割った。file・100entry block単位の5fold交差検証で、学習に使っていないeventの補正だけを採用。event・hit・gain選択は固定。

R=RMS/mean、G=R_data²−R_MC²とし、impactを1−G_after/G_beforeと定義する。

| 診断 | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| 位置のみ：相対分散差の減少 | 8.5% | 30.3% | 54.8% |
| 位置＋run：同 | 8.0% | 34.6% | 55.1% |
| run平均差のみ：同 | −0.7% | 4.6% | 0.5% |
| 温度補正係数を1へ：同 | 1.1% | 2.4% | 0.7% |
| 位置＋run後 Gaussian幅 data / MC | 7.75 / 7.49% | 4.82 / 3.78% | 2.07 / 1.96% |
| 位置＋run後 RMS/mean data / MC | 12.01 / 7.43% | 7.05 / 3.96% | 3.60 / 2.18% |

100 GeVの中心部分はかなり近づくがRMSの裾は残る。5 GeVのRMS差は位置補正で大きくは減らない。位置モデルの不確かさを含む信頼区間は今回は求めていない。

CoGはenergyから再構成した量なので、この改善にはstrip内光収集・channel間応答の不均一性に加え、shower形状・containment・再構成との相関が含まれる。strip固有の不均一性の因果寄与と同一視しない。外部track位置、または独立なlayer群の位置を用いた確認が次に必要。

[100 GeVの分布](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/comparison/20mm/100GeV/resolution_position_impact.png)。左右とも各sampleの平均でenergyを規格化しているので、平均応答差を示す図ではない。

layer共分散を分解すると、主比較の相対分散差のうち非対角項の差が約82/75/78%を占める。共分散の和はevent分散を再現することを検証した。これは分散の数学的な分解であり、独立した原因の寄与率ではない。位置/run/beam等と重複するため上表へ加算しない。[covariance.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/covariance.tsv)

## 2. 大信号応答と周辺strip：平均差の符号反転

中央は物理strip中心の|x|,|y|<22.5mm。実際のhit位置を長いstripの内部まで知っている定義ではない。

| beam | 中央energy data/MC | 中央hit数 data/MC | 周辺energy data/MC | 周辺hit数 data/MC |
|---|---:|---:|---:|---:|
| 5 | 0.817 | 0.974 | 0.713 | 0.904 |
| 20 | 1.055 | 0.979 | 0.895 | 0.879 |
| 100 | 1.161 | 0.992 | 0.817 | 0.840 |

100 GeVは中央hit数がほぼ同じなのにdataのenergyが16%高い。一方で周辺はenergy・hit数ともdataが少ない。20 GeVの平均が近いことは、中央と周辺の応答が一致したことを意味しない。

同じevent集合へ共通のcalibrated hit cutを0→0.3MeVにしても、data/MC平均比は20 GeV1.0207→1.0210、100 GeV1.0920→1.0927。低信号hitを後段で落とすだけでは主な平均差は解消しない。これは取得時に失ったhitの復元やhardware threshold変更の試験ではない。[hit_cuts.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/hit_cuts.tsv)

最新MCのtruthとdecodeも各energyの*_0.rootを1000eventずつ対応付けた。現在の採用MIPを用いた線形予測ADC `(truth edep/0.305)×MIP` に対する、pedestalとHL切片を処理したLG相当HG信号の比は以下。

| SiPM / truth edep | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| 10μm / 1–2MeV | 1.1037 | 1.1041 | 1.1047 |
| 15μm / 1–2MeV | 1.1482 | 1.1479 | 1.1455 |
| 10μm / 30–100MeV | 統計不足 | 1.0525（945hit） | 1.0168（20875hit） |
| 10μm / ≥100MeV | なし | なし | 0.9403（3989hit） |

比はratio of sums、取得されたHitTag=1、LG<3000のhitに限定。低信号で約10.4/14.5–14.8%の増幅、高信号で線形予測より低い応答が確認できる。現行MCDigiは実測MIP/gainを元に光量を設定しcrosstalk10/15%を追加、MIP÷1.1/1.15はコメントアウト。実測MIPがcrosstalk込みなら二重計上候補だが、electronのhit比からmuon MPVやtotal energyの補正量は確定できない。muonのraw ADC→MIP再現と有効pixel数・回復を含む大信号応答を検証する必要がある。[truth_scale.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/truth_scale.tsv)

## 3. 後方信号の多いevent：中心幅と裾を分ける

漏れの診断では上記energy窓を外し、低LG tail除外・good mask・CoG20を維持して再集計。track/粒子識別cutは追加していない。「最終5層の信号割合」はlayer25–29のpositive good energy / totalであり、真の漏出energy比ではない。

| beam | 最終5層割合>5% data / MC | RMS/mean：cut前→後 data | 同 MC | 相対分散差の減少 |
|---|---:|---:|---:|---:|
| 5 | 9.66 / 1.42% | 31.97→16.43% | 7.62→7.62% | 78.0% |
| 20 | 2.58 / 3.25% | 12.64→11.65% | 4.31→4.27% | 16.8% |
| 100 | 9.23 / 19.91% | 10.21→7.89% | 2.63→2.36% | 41.8% |

Gaussian中心幅はdataで5 GeV8.59→8.67%、20 GeV7.14→6.87%、100 GeV4.25→4.04%。MC100 GeVも2.40→2.28%。RMSの改善に比べ中心幅の差はあまり縮まらない。100 GeVで後方信号の多いeventはむしろMC側が多く、「dataで漏れが多いから中心resolutionが悪い」という単純な説明は支持されない。

元のenergy窓を維持した主比較では、この後方cutによる分散差減少は5/20/100 GeVで3.6/−1.3/−1.1%。母集団の定義でimpactが大きく変わることが確認できる。5 GeVの低energy側にはshower以外の成分やtrack様eventも混ざり得るため、上表の改善を全てshower leakへ帰属しない。別chatの旧shower_leak解析はtrack-v2除外を含むため今回と直接比較しない。[leak_proxy_rates.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/leak_proxy_rates.tsv)

今回のouter_fine_strip cutは、odd layerのX / even layerのYの細分方向で|座標|>103.35mmのstripにある信号割合を用いる。長軸方向の外縁を含まず、旧shower_leakのstripと10mm外縁帯の面積重なりとは定義が違う。この指標>5%は5 GeVdataで1/3073event、他は0eventで、代表sampleの主要因には見えない。横方向の真の漏れが無いことを証明した値ではない。

## 4. 温度：平均scaleへの感度は残る

補正係数Cだけ1に戻し、event/hit/gainとLG切片を固定した診断でdata平均は5/20/100 GeVで+11.92/+4.24/+3.62%。data/MC比は約0.890/1.064/1.132。5 GeVの不足の一部に対応する感度はあるが、20/100 GeVのdata過剰は大きくなる。幅の分散差減少は1.1/2.4/0.7%なので中心resolution差の主因としては小さい。

現在の温度再構成は16sensor、1/d重み、加算先0初期化・有限値・ゼロ距離等の検証へ修正済み。旧PS温度解析では保存温度と独立再構成は一致したが、MIP取得温度の由来は未確認、温度別MIP fitは統計不足、HG/LG residual傾きはrunと温度が交絡していた。現在も温度係数を独立に測定できた根拠にはしない。20℃基準を変える前に、MIP定数が何℃の応答を表すかを確認する。

旧解析参照：`analysis/result.bak/channel_response/temperature_calibration/ps_fast_png_v2/{summary.json,closure.tsv,reference_audit.tsv}`。保存後の不具合修正と入力更新を受けて、旧bugを現行原因として扱わない。

## 5. beam/SSA/threshold：別条件の証拠と現行条件を分ける

10月12日のbeam_quantitativeは100 GeVのpointとpoint_noTriggerSciのみ。通常e_ssaの同解析は未生成。lowLG除外・mask・CoG・energy窓は無く、data1355、MC1500eventなので主比較へ数値を合成しない。

それでもSSA前strip数はdata710.86に対しpoint MC865.12（+21.70%）、noTriggerSci862.72（+21.36%）。共通0.5/1/2MIP cut後にもpointの差は+21.38/+22.78/+24.40%。入口重心を引いた後のenergy横幅もlayer5/10/15/20/25でpoint MCがdataより約17/31/14/16/10%広い。低信号thresholdや入口中心だけで形状差は消えない。

truthに保存されたprimary開始位置はpoint z=−60mm、noTriggerSci z=−5mm。保存geometryでは41mm scintillatorがz=−50.5〜−9.5mmにあり、後者は下流開始によって通過を省いた比較（空気の距離も変わる）。strip数変化は−0.277%、横幅変化はlayerごとに約−1〜+2%で一貫した改善はない。通常e_ssaのresolutionへの物質効果を測った値ではない。

参照：`analysis/result/beam_quantitative/beam_scan/summary/scan/100GeV/{point,point_noTriggerSci}.{metrics.tsv,threshold.metrics.tsv,threshold.json}` と対応するevents ROOT・inputs。

現行e_ssaの5/20/100 GeVのdigitization代表ログではthreshold map読み込みに成功。338/6300channel（5.37%）をchip平均で補完。読込失敗fallbackを現在の主因にはしない。

一方、使用中の全energy合算threshold fitではχ²/ndf>2が5751/5962channel（96.46%）。同条件の各energy別fitでは4.11/24.71/44.93%。合算値とenergy別値の差の絶対中央値は10μmで7.76/2.78/1.25ADC、15μmで21.10/9.04/5.72ADC。合算signal spectrumのモデル検証と、独立した取得効率のturn-on確認が必要。高χ²はthreshold位置の誤りの証明ではなく、total resolutionへのimpactは未測定。5GeV別fitは2020channelのみで単純全置換は避ける。参照：`analysis/result/threshold/threshold{,_5GeV,_20GeV,_100GeV}.root` の `expErfThre`。

旧SSA shapeは10月5日・旧入力の結果。旧chip番号誤読による親strip欠落は現行StripSplitで10月6日に修正済み。旧SSAのbeam幅や親energy欠損を現在の原因・impactへそのまま流用しない。

## 修正済み事項と検証

現行Calibration/ExtractはMIPの正・有限MPV、NDF>0、有限非負χ²、接続channelのみの平均を採用。新しい定数抽出helperを現行Extract本体と比較し、6300channelのMIP/gain/pedestal差0、offset差は数値誤差1.14e−13ADC。現在のMIP fallbackは2422/6300channel（38.4%）で、これは旧「digi/calib不一致2422」と同じ意味ではない。平均置換が実際の各channelのMIP応答を表すかは別の課題。MCのdigitizationとcalibrationが同じ定数を使うため、数式closureだけではdataの残るchannel間応答差を保証できず、位置依存の原因候補に含める。

最新decode/calibのevent ID、full CellID、HG/LG式とHG=2600の採用を**13,557,919hit**で照合し、保存energyとの差は最大2.85e−14MeV。処理式の不一致は見つからない。[validation.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/validation.tsv)

共通channel（両sampleのHG800–2200で各100hit以上）の重みをそろえたEL−EHはdata5/20/100で+0.354/+0.252/−0.00375MeV、MC−0.00080/−0.00050/+0.000086MeV。data PS/低energy側のoverlapには残差があり、測定pedestal・HL比・温度条件をさらに確認できるが、MCの既対応intercept問題を再度原因としていない。good hitのLG3000到達は3energyとも0、MCのHG4000はLGを採用する表示・切替処理で、HGヒスト上限をenergy欠損と同一視しない。[gain_overlap.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261012/integrated/gain_overlap.tsv)

最新本番resolutionの5GeV/CoG5では、beforeのRMS23.37%に対し保存Gaussian幅46.79%、afterでもRMS8.40%に対しGaussian13.10%。統計は936/812event、sigma fit誤差14.47/2.16MeVで、CoG20/10と異なるfitの不安定さが見える。この点は物理的な差の根拠として使わず、binning・fit範囲・fit状態を確認する。今回のCoG20診断fitは全て成功。

再現コードは既存調査ディレクトリの [IntegratedReview.cc](../../investigations/e_ssa_threshold_review_20260930/IntegratedReview.cc)、[CurrentParameters.hh](../../investigations/e_ssa_threshold_review_20260930/CurrentParameters.hh)、[integrated_review.sh](../../investigations/e_ssa_threshold_review_20260930/integrated_review.sh)。入力一覧・前後の更新時刻/size・定数/コードchecksum・data tail checksum・実行記録は結果配下に保持。元ROOT・本番校正は変更していない。summaryのbaseline_*列は各selectionのnominal母集団の参考量であり、cut後の値ではない。

次の優先順は、(1)独立な位置からのlayer/channel応答、(2)muon MIP closureとSiPM大信号応答、(3)低energy側eventの深さ・track/粒子成分・後方割合、(4)threshold turn-onと補完channel、(5)温度係数の同一channel・run内検証。各診断impactは重複するため加算できない。

---

# 分解能差の追加調査（2026-10-04）

指定MC `/home/murata_t/data_beamtest/simulation/CEPCScECAL_SML_Portable_update_new/Result_MC/calib/e-/e_ssa/threshold` と、productionの低LG tail-rich判定を除外したdataを比較した。対象は5・20・100 GeVのみ。

**新しい所見は、位置に関連したevent応答の変動が分解能差のかなりの部分を担うこと。位置・run依存を別eventで学習して除く診断では、RMSで定義したdata/MC相対分散差が5/20/100 GeVで10/37/53%減った。Gaussian coreの差は特に20/100 GeVで大きく縮む一方、dataの分布の裾は残る。**

## 入力と再利用範囲

- 指定MCのrealpathは10月2日の調査入力と同一。今回、選択した全data/MCのdecode・calibのmtime/size、校正定数・digitizer/Extractソース・production tail判定のSHA256を照合し、変更がないことを確認した。
- 現在のMCファイル数は5/20/100 GeVで100/68/100。辞書順に広く6ファイルを選ぶ現行ルールでも、同じ6ファイルになる。全energy・全eventの再解析はしていない。
- 上記を確認して、10月2日の検証済み `events.tsv` を読み直し、新しい位置/run補正、Gaussian fit、impact計算を実施した。今回は元ROOTを再走査したという意味ではない。
- nominal非負good-hit和、既存shower選別、CoG±20 mm、低LG tail-rich除外を維持。data/MC採用数は5 GeV 2753/2870、20 GeV 1755/2719、100 GeV 5001/4566。MC側には該当tail-rich eventがない。
- 補正診断は同じevent・同じhit集合のまま実施し、元ROOTや本番校正を変更していない。

[追加結果一式](../../investigations/result/e_ssa_threshold_review_20260930/resolution_impact/)、[RMS impact表](../../investigations/result/e_ssa_threshold_review_20260930/resolution_impact/impact.tsv)、[Gaussian fit安定性](../../investigations/result/e_ssa_threshold_review_20260930/resolution_impact/gaussian_scan.tsv)。

## 1. RMSとGaussian coreを区別する

| beam | baseline RMS/mean data / MC | baseline Gaussian sigma/mean data / MC | 位置+run補正後 RMS/mean data / MC | 同 Gaussian sigma/mean data / MC |
|---|---:|---:|---:|---:|
| 5 GeV | 12.34 / 7.86% | 8.84 / 8.06% | 11.77 / 7.59% | 8.11 / 7.49% |
| 20 GeV | 8.70 / 4.62% | 7.63 / 4.17% | 7.18 / 4.14% | 4.96 / 4.04% |
| 100 GeV | 5.30 / 3.03% | 4.52 / 2.84% | 3.69 / 2.18% | 2.19 / 1.92% |

Gaussian fitは同一の選別sampleへ、平均±3RMSの初回fit、得られた平均±1.5sigmaの2回目fitを行う診断。既存resolution解析の二段階fitにならうが、binningとsample/選別はこの調査固有なので、本番resolution図の数値そのものではない。主表は240bin、範囲[0.3,1.7]×各sample平均。

補正後のGaussian比のfit統計誤差はdata/MCで、5 GeV ±0.28/0.27 percentage point、20 GeV ±0.23/0.15、100 GeV ±0.053/0.054。120/240/360bin、最終fit範囲±1.5/2sigmaの72fitすべて成功した。補正後の範囲はdata/MCで5 GeV 8.06–8.47/7.35–7.53%、20 GeV 4.91–5.02/3.98–4.09%、100 GeV 2.17–2.24/1.92–2.02%。fit設定による変化も残るため、単一のfit誤差だけで一致・不一致の有意度を決めない。

中央68%幅でも、補正後data/MCは5 GeV 8.92/7.53%、20 GeV 5.39/3.96%、100 GeV 2.36/2.02%。位置補正で中心は近づくが、完全一致とはいえない。

## 2. 原因別の簡単なimpact

RをRMS/meanとし、`G=R_data^2-R_MC^2`、改善率を`1-G_after/G_before`と定義する。これは平均energyの変化率やresolutionそのものの改善率ではない。

| 診断変更 | 5 GeV：分散差の減少率 | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| 温度係数を1に戻す | 1.4% | 2.3% | 0.8% |
| MIP採用条件を健全化（双方） | 3.6% | 6.0% | 7.1% |
| 校正MIP→digi MIP置換（双方） | 4.7% | 5.3% | 6.3% |
| run平均差のみ補正、交差検証 | −1.4% | 4.4% | 0.5% |
| 位置依存のみ補正、交差検証 | 11.6% | 32.7% | 52.5% |
| 位置+run補正、交差検証 | 10.4% | 36.7% | 52.8% |

各変更は独立にbaselineへ適用しており、寄与を加算できない。MIPや温度の絶対energyへの大きな影響を、そのまま分解能差への大きな影響と解釈しない。

位置補正はlayer9 CoG x、layer10 CoG yの二次曲面 `1,x,y,x^2,xy,y^2`、run補正はrun別の定数項。data/MCを別々に学習する。fileと100-entry blockから決めた5foldのうち4foldで係数を求め、残り1foldに適用し、全eventについてheld-out値を使用した。energyは予測相対応答で割る。MCはRun_Num=0のためrun定数は一つ。

100-entry blockを単位に400回bootstrapした参考範囲（16–84 percentile）は、位置+runの減少率で5 GeV 8.1–13.3%、20 GeV 32.5–41.1%、100 GeV 48.6–57.6%。学習済み補正値を固定した再標本化なので、補正モデル選択・再学習の不確かさをすべて含む区間ではない。

前回のCoG matchingは「位置の分布をdata/MCでそろえる」処理。今回は「各sampleの位置に対する平均応答変化を除く」処理であり、確認している内容が異なる。CoGはenergyから再構成した量なので、strip固有の不均一性だけでなくshower構造や再構成との相関も含む。外部trackがない現段階で、位置補正の改善をすべてstrip応答の因果寄与とは断定しない。

## 3. 残る原因候補

1. **位置に関連した応答変動：定量的に最も大きい。** strip内の光収集、channel間の応答差、layerごとの位置依存、shower/漏れとの相関を優先する。次は外部trackの位置、または独立なlayer群で推定した位置に対して同様の変動があるか確認する。
2. **低LG異常以外のevent分布の裾：特に5 GeVで重要。** 位置+run補正後でもRMSは11.77/7.59%だが、Gaussian coreは8.11/7.49%。各sampleのfit平均−3sigmaより低いeventはdata約1.49%、MC約0.035%、+3sigma以上は約1.93/0.38%。beam成分、上流損失、containment、通常eventの取得条件を確認する。これは処理済みの低LG tail判定とは別のenergy分布の裾であり、追加除外を推奨したものではない。3sigmaの絶対閾値は各sampleで異なる。
3. **残る相関変動・beam幅・shower再現。** 位置+run補正後にもRMS差があり、独立な相対揺らぎ1項へ換算すると5/20/100 GeVで9.00/5.87/2.98%。これはbeam spreadの測定値ではない。既存のlayer共分散の差と位置依存は重複し得るため、両方を別原因として加算しない。
4. **大信号応答・飽和：未分離。** 100 GeVのLG採用hitのenergy寄与はdata82.2%、MC77.0%。HG/LGに分けた分散ではLG側の違いが大きいが、LGは高信号hitを選んだ集合なので、これだけでLG electronicsの問題とはいえない。位置依存をそろえた後で大信号の非線形性を比較する必要がある。
5. **温度・MIPの誤適用：点検は必要だが、今回の幅差の主因としては小さい。** 固定hit集合での単独診断の改善は上表の通り。除外channelの復帰などを含む全面再校正の評価ではない。

## 実装・検証

既存調査ディレクトリへ `ResolutionImpact.cc` と `resolution_impact.sh` を追加し、Makefileと実行方法を更新した。再実行は `bash investigations/e_ssa_threshold_review_20260930/resolution_impact.sh`。入力が更新されていれば過去のevent表を再利用せず停止するため、その場合は先にraw/calib監査を更新する。

baselineのevent数・平均・RMSを既存独立集計と照合し、平均の差は2.5e-12 GeV未満、RMS比の差は2.1e-13未満。校正変更の値は元の保存event列を使用。補正eventのfile/entry同一性、Gaussian fit成功、shell構文、代表PNGを確認した。原データ・本番解析は変更していない。一時テスト用ファイルは新規作成していない。

- [位置補正前後の100 GeV](../../investigations/result/e_ssa_threshold_review_20260930/resolution_impact/comparison/20mm/100GeV/figures/resolution_position_impact.png)
- [位置依存の20 GeV](../../investigations/result/e_ssa_threshold_review_20260930/resolution_impact/comparison/20mm/20GeV/figures/energy_position_profile.png)
- [再現コード](../../investigations/e_ssa_threshold_review_20260930/ResolutionImpact.cc)

---

# 最新入力による再調査（2026-10-02、5 / 20 / 100 GeV）

対象は `Result_MC/{decode,calib}/e-/e_ssa/threshold/` と2023 PS/SPS data。10月2日に指定先を確認し、保存ROOTを読み直した。指定先のMC decode/calibは9月30日生成分であり、5/100 GeVの切片問題への対応は既に含まれる。低LG異常eventは、同じ入力を全event処理した既存 `hg_lg_tail/threshold/e_ssa/<energy>GeV` の判定で除外した。除外済みROOTへの置換ではなく、解析時のevent除外である。

**主な結論：平均応答の差にはenergy依存性があり、20 GeVで平均がほぼ一致しても分解能差は残る。MIP定数の健全化は絶対scaleを約14%動かすが、data/MC比は約0.2 percentage pointしか動かない。分解能差ではlayer間のenergy相関の違いが大きく、単純な一律校正では解消しない。**

- [今回の結果一式](../../investigations/result/e_ssa_threshold_review_20260930/20261002/)
- [影響度の数値表](../../investigations/result/e_ssa_threshold_review_20260930/20261002/impact_summary.tsv)
- [共分散の寄与](../../investigations/result/e_ssa_threshold_review_20260930/20261002/covariance_impact.tsv)
- [コード・実行方法](../../investigations/e_ssa_threshold_review_20260930/README.md)

## 条件と結果

各energyでMCを辞書順に広く6ファイル選択、各1000 event。利用可能ファイルは5/100 GeVが各100、20 GeVが68。dataは全runから各2000 eventを10区間に分散して抽出した。5 GeVは4run=8000 event、20 GeVはRun96/97の4000 event、100 GeVは3run=6000 event。抽出数はrunの全event数では重み付けしていない。

比較量は既存bad-channel maskを使った非負good hitのenergy和。全接続channelの非負energy和がMCの符号付きenergy中央値の0.5–1.5倍、layer9/10のgood-energy CoGがx/yとも±20 mmという従来のshower選別を維持した。その後、full-sample tail判定の異なる3 channel以上がz<−5のeventを除外。tail判定の `selected` は使わず、event IDを照合した `tail_rich` だけを使用した。

| beam | data/MC採用event数 | data平均 [GeV] | MC平均 [GeV] | data/MC−1 | data RMS/mean | MC RMS/mean |
|---|---:|---:|---:|---:|---:|---:|
| 5 GeV | 2753 / 2870 | 0.164539 | 0.213956 | −23.10% | 12.34% | 7.86% |
| 20 GeV | 1755 / 2719 | 0.844894 | 0.851639 | −0.79% | 8.70% | 4.62% |
| 100 GeV | 5001 / 4566 | 4.314263 | 4.031946 | +7.00% | 5.30% | 3.03% |

平均の標本SEMはdata/MCで、5 GeV 0.000387/0.000314、20 GeV 0.001755/0.000755、100 GeV 0.003232/0.001808 GeV。RMS/meanは選別後分布の診断量であり、Gaussian fitによる正式なresolutionではない。

非Gaussian tailへの依存も確認した。中央68%幅 `(q84−q16)/(2×median)` はdata/MCで、5 GeV **9.73/7.71%**、20 GeV **7.60/4.25%**、100 GeV **4.64/2.88%**。RMSの外側tailだけでなく、中心部分にも差がある。[robust_width.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261002/robust_width.tsv)

## 1. decode→calibrationの数式・対応は正常

3 energy、data+MCの計13,556,786 hitでraw/calibのevent IDとfull CellIDを照合した。保存HG/LG energyの再計算との差は最大2.28e−13 MeV未満、raw HG=2600による採用gainも一致。decode→calibrationのhit残存率は調査項目にしていない。

低LG異常除外後、共通channel（各sampleのHG 800–2200で100 hit以上）の重みをそろえたMCのHG/LG energy差は、5/20/100 GeVで −0.00075/−0.00109/+0.00043 MeV。MC切片の欠落を再び原因として扱う根拠はない。good hitのLG=3000 ADC上限到達は3 energyともdata/MCで0。HG=4000置換のhitはLGを採用しており、HGの表示上限だけをenergy欠損と解釈しない。

一方、dataの同じ重み付きHG/LG差は +0.416/+0.295/−0.0043 MeV（共通channel数170/358/1486）。5/20 GeVのdata内部の接続には残差がある。これはMCの既対応切片問題と区別し、dataのpedestal・gain比・温度条件で検証する。overlap領域の差をHG=2600の段差や全LG hitへの補正量へ外挿してはいけない。

参照：[gain_summary.tsv](../../investigations/result/e_ssa_threshold_review_20260930/20261002/gain_summary.tsv)、各energy `raw/validation.tsv` および `sensitivity/comparison/20mm/<energy>GeV/validation.tsv`。`raw` は除外前の監査、`sensitivity` のchannel表・図・主結果はtail除外後。

## 2. MIP・光量：絶対scaleへの感度とdata/MC差への影響を区別する

今回も6300接続channel中、校正MIPが負のもの565、digi/calib不一致2422、校正fallback1857。主比較の非負good energyのうちfallbackの寄与はdataで約85–87%、MCで約84–87%。

同じ採用event・同じ非負good hit集合について、以下の診断的な係数置換を行った。元のROOT・校正定数は変更していない。除外channelや負energy hitは復活させず、eventを再選別しない。

| 診断 | 5 GeV | 20 GeV | 100 GeV |
|---|---:|---:|---:|
| MCだけ校正MIP→digi MIP：MC平均の変化 | −15.98% | −16.34% | −15.82% |
| data/MC両方に同じ置換：data/MC比 | 0.7636 | 0.9863 | 1.0646 |
| data/MC両方のMIP品質を健全化：data平均の変化 | −13.54% | −13.97% | −13.96% |
| 同：MC平均の変化 | −13.37% | −13.83% | −13.78% |
| 同：data/MC比 | 0.7675 | 0.9904 | 1.0678 |

健全化はMPV>0、finite、NDF>0、chi2>=0および既存chi2/NDF cutを要求し、平均置換値を再計算したもの。この診断のdata/MC比の変化は−0.16/−0.17/−0.23 percentage pointであり、5/100 GeVの相対差を大幅に縮めない。分解能差も残る（健全化後data/MC RMS/meanは5 GeV 12.11/7.71%、20 GeV 8.40/4.42%、100 GeV 5.03/2.78%）。

MC truth→ADCでは、MIP一致10 μm channelのtruth 1–2 MeV/hitにおけるHG/線形予測が5/20/100 GeVで1.1025/1.1035/1.1031。約10%の低信号増幅が共通にある。実測MIPから求めた光量にcrosstalkを再度加えている可能性は残るが、electron hitの比はmuonのMIP MPV一致を検証したものではない。補正を確定するにはmuon closureが必要で、total energyを単純に10%下げることは推奨しない。

## 3. PSの温度補正は平均に大きく効くが、幅の差はほぼ残る

温度補正係数だけを1に戻し、採用hit・gain・LG切片を固定する診断では、data平均が5/20/100 GeVで **+11.89/+4.24/+3.62%**。data/MC比は **0.8604/1.0341/1.1087**。5 GeVの不足は約23%から約14%に縮むが、20/100 GeVではdata側を高くする。

同じ診断のdata RMS/meanは **12.34→12.28%、8.70→8.63%、5.30→5.28%**。今回のevent間の幅差は温度補正を外すだけではほとんど消えない。energy重み付き温度はdata約29.0/23.5/23.0℃、MC20℃。基準温度・MIP取得時の条件と、同一channelの温度依存で校正式を検証する必要がある。

## 4. 中央大信号と周辺hit数の差はenergy依存

中央はstrip中心の|x|,|y|<22.5 mm。以下はtail除外後。

| beam | 中央energy data/MC | 中央hit数 data/MC | 周辺energy data/MC | 周辺hit数 data/MC | 中央の差 / MC total | 周辺の差 / MC total |
|---|---:|---:|---:|---:|---:|---:|
| 5 | 0.7873 | 0.9723 | 0.6978 | 0.8978 | −16.93% | −6.16% |
| 20 | 1.0207 | 0.9749 | 0.8833 | 0.8756 | +1.64% | −2.43% |
| 100 | 1.1301 | 0.9913 | 0.8187 | 0.8402 | +10.50% | −3.50% |

100 GeVは中央のほぼ同じhit数に対しdataのenergyが13%大きい。20 GeVは中央と周辺の差が部分的に相殺してtotalが近い。したがって20 GeVの平均一致は、全channelの応答が一致したことを意味しない。100 GeVの大信号ではSiPM飽和、有効pixel数、非線形応答、strip内応答を優先して確認する。周辺ではthreshold/検出効率とshower形状・上流物質を分ける。周辺hit数比をそのまま検出効率とは呼ばない。

共通のcalibrated hit cutを0→0.3 MeVに上げても、20 GeVのdata/MC比は約0.992、100 GeVは約1.070のまま。0.3 MeV以下のhitだけが主な平均差を作っているわけではない。ただし、この後段cutはhardware thresholdの変更試験ではない。

## 5. 新しい所見：layer間共分散が相対分散の差の約8割を占める

採用eventごとに30 layerのenergyを保持し、次式でtotalの分散を分解した。

`Var(sum E_layer) = sum Var(E_layer) + sum_{i!=j} Cov(E_i,E_j)`

data/MCそれぞれのtotal平均の二乗で規格化した後に差を取ると、その差のうち非対角共分散項が占める割合は **5 GeV 79.9%、20 GeV 84.3%、100 GeV 83.8%**。MCではlayer間の増減が相殺する負の共分散の合計が強く、dataは相殺が弱い。個々のlayerの幅だけを見ると、totalの分解能差を見落とす。

これは数学的な分散分解であって、「特定の物理原因が8割を説明した」という因果寄与ではない。beam energy幅、位置依存応答とshower位置の揺らぎ、event共通のgain変動、shower/leakageの相関再現などが候補。独立Gaussian noiseを各hitに足す案より、相関を実測と照合する優先度が高い。

- CoGを5×5 mmの共通binで重み付けしても、非対角項の割合は78.3/83.6/79.0%。単なるCoG分布の混合では消えない。
- より中心（layer9 x、layer10 yが±5 mm）へ限定すると幅は特に20/100 GeVで狭くなるが、data/MCは20 GeV 5.91/4.19%、100 GeV 4.10/2.13%で差が残る。この狭い条件の統計はdata/MCで163/233、524/574 event。
- run別にも幅差がある。20 GeVのRun96/97は8.91/8.21%、MC4.62%。100 GeVのRun105/106/107は6.00/4.92/4.92%、MC3.03%。run平均の混合だけでも説明できない。
- MC truthのPrimaryEnergyは3 energy各6000 eventすべて指定値と一致。実dataのbeam energy spreadは今回未測定。

もし全幅差を独立な相対揺らぎ1項で表すなら必要な大きさは **9.51/7.37/4.35%**。これはRMSに基づく換算量であり、beam spreadの測定値ではない。5 GeVは中央68%幅とRMSの差も大きく、別のevent tail・beam成分を含めて調べる必要がある。

[共分散の数値](../../investigations/result/e_ssa_threshold_review_20260930/20261002/covariance_impact.tsv)。layerの平均・幅の図と、各run/条件の共分散分解は `sensitivity/comparison/20mm/<energy>GeV/{figures/layer_response.png,covariance_summary.tsv}`。

## 次に調べる優先順位

1. **分解能：event単位の相関の出所。** beam energy/入射位置・角度の実測、layer/chip間の共通応答変動、showerの前後layer間の相関を比較する。20 GeVは平均差が小さいため、この検証に向く。
2. **平均：PSの温度基準とMIP closure。** 5 GeVの約23%不足をdata校正とMC低信号規格化に分ける。同じMIP修正をdata/MC両方へ適用した影響を評価する。
3. **高energy：中央の大信号応答。** 低信号の規格化を固定した上で、飽和・strip内応答を調べる。周辺のhit不足は別に検証する。

## 再現・検証・保存

既存調査コードを20 GeVとlabelに対応させた。`run.sh <5|20|100> <label>`、完了済みの入力監査に対する追加処理は`extra.sh <energy> <label>`、集計は`summarize.sh <label>`。実行済みlabelは`20261002`。完成したbaselineがある場所への`run.sh`再実行は停止する。

今回の5 GeVでは旧HLInterceptStudyの不要なchannel histogram保存でI/O待ちが再発したため、自分の処理だけを終了し、`PrepareReference.cc`で必要なMC event基準値のみをメモリ上でROOT化して保存した。20/100 GeVは旧baseline生成が完了しており、共有するpaired_eventsの定義は同じ。入力ROOTを読み直す本体のraw/calib・truth・感度検証は3 energyとも完了した。今後のrun.shも軽量なPrepareReferenceを利用する。

- 原データ、MC、本番digitizer/calibration、校正定数は変更していない。
- raw/calib対応、gain選択、保存energyの再現を全サンプルで検証。最新Sensitivityのkeep集計は別実装InspectResponseの集計と一致。
- production tailのfull input path・全entry数・Run/Event_Time/Event_Numを照合。入力mtimeがtail判定より新しければ停止する。選択した入力のsize/mtimeは調査前後で不変。
- layer平均の和・共分散全要素の和が、event平均・分散を再現することを、keep/exclude/位置matching/中心限定/run別の全条件で検証。
- PNGを目視確認。図・表・ROOT・入力一覧・チェックサム・再現用ソースは今回の結果内に保存。PNGのみでPDFは出力していない。
- 一時I/O再試行の不完全ROOTと`.scratch`は削除済み。実行経緯を示すログは結果内に保持。

---

# 前回調査（2026-09-30、履歴）

以下は旧サンプル内でtailを再定義した除外前中心の報告。今回の主結果と混在させない。

# e_ssa/threshold: decodeからcalibrationまでのdata/MC比較

2026-09-30。MCは指定された `Result_MC/{decode,calib}/e-/e_ssa/threshold/`。
以前の `e-/sps/threshold/` や `e-/e_ssa_noHLintercept/threshold/` の結果と区別する。
全energyは実行せず、100 GeVと5 GeVだけを確認した。入力ROOT・本番のdigitization/calibration・解析の既定入力設定は変更していない。

**現時点で、差は一つの原因に集約できない。旧MCのHG/LG切片欠落は今回のMCでは解消している一方、MIP/光量の規格化、PSの温度・HG/LG校正、大信号応答、周辺hit効率、data固有のevent相関異常を分けて調べる必要がある。**

## 入力と比較条件

- 100 GeV data: Run105/106/107から各2000 event、計6000。5 GeV data: Run63/64/75/76から各2000、計8000。
- dataは各ファイルを10区間に分け、各区間中央から抽出。MCは各energy100ファイルの辞書順で広く6ファイルを選び、各1000 event、計6000。
- 実ファイル一覧は各energyの `baseline/inputs/{data,mc}.tsv`、抽出記録は `raw/inputs/sampled_files.tsv`。`baseline` はMC全6000 eventの整合性・参照値取得に使う。そこにあるdata先頭6000 eventの補助集計は、本文のdata比較値には使っていない。
- 主比較は非負good hitの和。比較用bad-channel maskは `analysis/cog_channel_scan/SingleEnergyAnalysis.cc`。中央はstrip中心の |x|,|y|<22.5 mm。
- 本文の主比較は `calibrated` / `sensitivity` のshower+CoG20: 全接続channelの非負energy和がMCの符号付きenergy中央値の0.5–1.5倍、layer9/10のCoGがx/yとも±20 mm。
- `tail` は既存hg_lg_tailの定義（good非負energy和と、そのMC中央値を使う窓）。100 GeVの選別eventは主比較と同一。5 GeVのdataは主比較2796、tail側2793と3 event異なる。数値を混ぜない。
- RMS/meanは選別後の標本分布の診断値。正式なGaussian fitによる分解能ではない。

| 入射energy | data/MC 選別event数 | data平均 [GeV] | MC平均 [GeV] | data/MC | data RMS/mean | MC RMS/mean |
|---|---:|---:|---:|---:|---:|---:|
| 100 GeV | 5118 / 4566 | 4.31008 | 4.03195 | 1.06898 | 5.364% | 3.030% |
| 5 GeV | 2796 / 2870 | 0.164522 | 0.213956 | 0.76895 | 12.334% | 7.863% |

参照: [100 GeV](../../investigations/result/e_ssa_threshold_review_20260930/100GeV/sensitivity/summary.tsv)、[5 GeV](../../investigations/result/e_ssa_threshold_review_20260930/5GeV/sensitivity/summary.tsv)。

## 1. 校正式は整合しており、今回のMCに切片欠落は見られない

raw/calibのRun・Event_Time・Trigger/Event_Numとfull CellIDを対応付けた。100 GeVでdata 3837758 hit、MC 5212106 hit、5 GeVでdata 623861、MC 789337 hitを検証。保存HG/LG energyと再計算値の不一致は0、最大誤差2.28e-13 MeV未満。raw HG=2600によるHG/LG選択も一致した。

raw HG 800–2200、両sample各100 hit以上のgood正MIP channelを選び、channel重みmin(Ndata,Nmc)で比較した結果:

| energy | 共通channel | data平均 E_LG−E_HG [MeV] | MC平均 E_LG−E_HG [MeV] |
|---|---:|---:|---:|
| 100 GeV | 1285 | −0.04688 | +0.000188 |
| 5 GeV | 310 | +0.35148 | −0.000552 |

MCのraw LG残差 `L−H/A−b` も100 GeV +0.00404 ADC、5 GeV −0.00597 ADCとほぼ0。今回の保存MCはraw切片を含むHG/LG関係と整合する。旧MCに対する「calibrationの切片を0にする」という対処を今回へそのまま適用してはいけない。実際、100 GeVの同じ選別eventでMCの切片だけを外すと平均4.03195→3.94304 GeVとなり、dataとの差は拡大する。

一方、5 GeV dataのHG/LGずれは残る。channel中央値を重み付けしたraw残差はdata +2.479 ADC、MC −0.0073 ADC。一定のraw offsetだけでも校正後の接続に効く。PSのpedestal・HL係数・温度依存をSPSと分けて再確認する価値が高い。HG 800–2200での差を、そのままHG/LG切替点2600 ADCでの段差量とは解釈しない。

参照: [gain_summary.tsv](../../investigations/result/e_ssa_threshold_review_20260930/gain_summary.tsv)、各energy `raw/validation.tsv`, `tail/channels.tsv`。

## 2. MCの光量規格化とMIP定数の整合性が残る

MC truthとdecodeを各energy6000 event対応付けた。PrimaryEnergyと指定入射energyは一致、decoded hitはtruthのCellIDに対応する。
10 um SiPM、truthの1–2 MeV/hit領域で、予測HG ADCを `(truth E/0.305) × digitization MIP` とすると:

| energy | digi/calib MIP一致channelのHG/予測 | 不一致channelのHG/予測 | 不一致channelの補正LG/calib予測 |
|---|---:|---:|---:|
| 100 GeV | 1.1031 | 1.1038 | 1.3193 |
| 5 GeV | 1.1025 | 1.1039 | 1.3341 |

HGはraw<2400に限定、LGはraw切片を除いてHG ADCへ換算。1–2 MeVの10 um領域ではHG利用可能率が高く、HG上限による大きな切断を避けている。15 umの同様なHG比は上限選択のbiasが大きいため、ここには混ぜない。

現行MCDigi.cxxは実測MIP/gainから光量を求め、その後にcrosstalkを入れる。MIPをcrosstalk分で割る処理はコメントアウトされている。保存MCに約10%の低信号増幅があることは確認できた。実測MIPがcrosstalk込みの値なら二重計上の可能性があり、muon/MIPのtruth→ADC→校正の閉じが最優先。ただし現在のソースが生成時と完全に同じだったとの履歴保証はない。

さらに接続6300 channel中2422でdigi/calibのMIPが異なり、校正fallbackは1857。100 GeVの選別後good energyのfallback寄与はdata87.4%、MC86.2%。信号を多く担うchannelで個別値と平均値の採否が異なる。これは設計上の「真の応答」と「推定校正」の差を表す場合もあるため、定数差だけでバグとは断定しない。

同じhit集合のMC energyを校正MIPからdigi MIPへ換算する感度では、100 GeV 4.03195→3.39420 GeV、5 GeV 0.213956→0.179766 GeV。5 GeVの差は縮まるが100 GeVでは拡大する。dataにも同じ係数置換をすればdataも動く。これは修正案の効果を確定したものではなく、一律スケール調整で両energyを合わせられない証拠。

参照: [MIP採用値](../../investigations/result/e_ssa_threshold_review_20260930/mip_audit.tsv)、[100 GeV truth response](../../investigations/result/e_ssa_threshold_review_20260930/100GeV/truth/response.tsv)、[5 GeV truth response](../../investigations/result/e_ssa_threshold_review_20260930/5GeV/truth/response.tsv)。category 0/1は15 umでMIP一致/不一致、2/3は10 umで一致/不一致。

## 3. 負のMIP fit値が校正に採用される問題も残る

6300 channel中565で校正MIPが負。失敗fitのMPV=chi2=NDF=-10でもchi2/NDF=1となり、採否判定を通る旧来の挙動が現在の保存energyにも整合する。負値はfallback平均にも影響する。

選別後の全接続channelの負energy和は、100 GeVでdata −0.01656、MC −0.18983 GeV/event、5 GeVでdata −0.001341、MC −0.009562。大部分が負MIP channel由来。符号付き全hit和をそのまま比較すると強く歪む。

本文の主比較は非負good hit和なので、負energyを直接足しているために主比較のMCが低いわけではない。それでも、MIP quality cutとfallbackの健全化は必要。負fit採用・平均への混入を正値/finite/NDF等で防ぐ案を、data/MCの両方で評価すべき。

## 4. 100 GeVは中央の大信号と周辺hit数で差の向きが異なる

| 領域 | Nhit data/MC | energy data/MC |
|---|---:|---:|
| 全good channel | 0.8811 | 1.0690 |
| 中央 | 0.9904 | 1.1289 |
| 周辺 | 0.8399 | 0.8183 |

中央のhit数はほぼ同じだが、energyはdataが12.9%高い。中央hit数190–191の同じbinでも、data3.6943、MC3.2662 GeVと差が残る。単なるhit数差ではない。選別後raw LG分布の大ADC側とLG採用hitの高energy側でもdataが厚い。

周辺はMCのhitが約19%多い（data/MC=0.840）。threshold/検出効率、shower横方向・縦方向の分布、beam/上流物質等が候補。ただしdataにはtruth分母がないため、これをそのままdataの検出効率と呼べない。

大信号側の候補はSiPM飽和モデル、有効pixel数、光量/gain/校正非線形性、strip内応答の再現。MC内部では、MIP一致10 umの補正LG/線形digi予測はtruth1–2 MeVで1.104、30–100 MeVで1.012、100 MeV以上で0.931に下がり、飽和は効いている。ただしdataとの差が飽和の強さだけによるかは未確定。

LG=3000 ADCの上限到達hitはgood正MIP channelの今回sampleでdata/MCとも0。HG=4000の置換はMCに多数あるが、校正ではLGを使うため、HG図の途切れだけでenergy低下を説明しない。低信号での早期HG置換設定は継続課題だが、`raw/channels.tsv` の `HG4000_impliedHG_lt2600` は旧切片無しMC用の近似で、新MCのraw切片を戻していないため、その件数を物理的な早期切替率として採用しない。

参照: [region_summary.tsv](../../investigations/result/e_ssa_threshold_review_20260930/region_summary.tsv)、[固定hit数比較](../../investigations/result/e_ssa_threshold_review_20260930/100GeV/calibrated/conditional_energy.tsv)、[raw ADC/LG energy図](../../investigations/result/e_ssa_threshold_review_20260930/100GeV/raw/figures/raw_adc_and_lg_energy.png)。

## 5. 位置分布をそろえても主な差は残る

今回のMC truthの入射位置は100 GeVで平均(x,y)=(-6.54,+4.85) mm、RMS=(15.00,13.61) mm、5 GeVで(-3.49,-7.79)、RMS=(22.06,23.32) mm。旧sps MCの中心入射を前提にした議論は流用しない。

100 GeVのshower選別後CoG20通過率はdata5118/5272=97.1%、MC4566/5994=76.2%。MCのCoG分布はdataより広い。
そこでCoG20後、layer9のxとlayer10のyを5×5 mm区画で分け、双方20 event以上ある区画に共通のmin(Ndata,Nmc)重みを与えた。

| energy | data平均 → matched [GeV] | MC平均 → matched [GeV] | matched RMS/mean data / MC |
|---|---:|---:|---:|
| 100 | 4.3101 → 4.3200 | 4.0319 → 4.0363 | 5.500% / 2.999% |
| 5 | 0.16452 → 0.16495 | 0.21396 → 0.21425 | 12.086% / 7.865% |

5 mm区画で揃えたCoGだけでは差を説明できない。ただしCoGは入射位置そのものではなく、strip応答・shower揺らぎを含む。入射角・全layerの分布まで一致させた検証ではない。
100 GeVのdataをrun別にしてもRMS/meanはRun105 6.05%、106 4.95%、107 5.03%で、MC3.03%より広い。run平均の混合だけでも説明できない。

参照: 各energyの `sensitivity/{summary,runs,cog_weights,events}.tsv`。

## 6. PSの温度補正とHG/LG接続は独立に再確認が必要

energy重み付き温度は100 GeV data約23℃、5 GeV data約29℃、MC20℃。同じ選別event・同じ非負hit集合のまま、保存式の温度係数Cだけを1に戻す診断では:

- 100 GeV data: 4.31008→4.46605 GeV（+3.62%）。MCとの差は拡大。
- 5 GeV data: 0.164522→0.184076 GeV（+11.89%）。MC0.213956との差は縮まるが、約14%残る。

温度補正はPSの差に無視できない。ただしMIP測定時の基準温度・符号が正しいか確認せず、補正を消すことを推奨するものではない。
PSにもSPS名のpedestal/HL/MIPファイルを用いた保存校正が再現された。ファイル名だけで誤適用とは断定できないが、PSでのHG/LG raw残差+2.48 ADCと校正後差+0.351 MeVは、温度と係数の適用条件を確認する具体的な動機になる。

HitTag=0のpedestal±5σ内baselineについて、good channelの残差中央値はHG/LGで100 GeV(+1.79,+4.49) ADC、5 GeV(+1.37,+3.92) ADC。trigger・小信号・切断biasがあるので、これをそのまま新pedestalとして使えない。MCにはHitTag=0がなく、同じ母集団の比較もできない。

## 7. dataの低LG tailは実在するが、分解能差の大部分は残る

| energy | 共通channel | data低tail率 | MC低tail率 | 選別後data tail-rich | data RMS/mean除去前→後 | MC RMS/mean |
|---|---:|---:|---:|---:|---:|---:|
| 100 | 1285 | 0.6551% | 0 / 1257399 hit | 116 / 5118 | 5.364→5.298% | 3.030% |
| 5 | 310 | 0.3877% | 0 / 151446 hit | 43 / 2793 | 12.299→12.301% | 7.863% |

率はchannel重みをそろえた選別前hit、tail-richは選別後eventで異なる3 channel以上にz<−5がある場合。100 GeVの選別後負LG hit705個は全て116 eventに集中。5 GeVでも14個全てが43 eventに含まれる。中央値による切片除去後にも残るため、一定の切片差ではない。

electronics/hold/readout/decodeのevent共通異常が候補だが、元binary照合は未実施。独立Gaussian noiseを全体に足すだけではこの集中を再現しない。tailを除いた後のdataにも大きな幅差が残り、通常eventの応答変動・beam energy幅等も別に調べる必要がある。MCのPrimaryEnergyは今回のsampleで全event指定値と一致しており、data beamのenergy幅を表す変動はここにはない。ただしdataのbeam energy spreadの実測値は今回確認しておらず、寄与量は未定。

## 次の優先順位

1. **muon/MIPの規格化を閉じる。** digitizationのMIP・crosstalk・SiPM gainと、calibrationのMIP採用条件を同じ基準で検証。負fitとfallback平均の問題を分離する。
2. **PSの校正条件を点検する。** 温度基準・符号、PS時のpedestal/HL直線をchannelごとに確認。SPSと共通の固定値が妥当か、HG採用域とLG採用域を分ける。
3. **100 GeVの中央大信号と周辺小信号を独立に検証する。** 校正が閉じた状態で飽和/光量とthresholdを個別に変える。global gainやHG上限だけで合わせない。
4. **通常eventの幅差を追う。** run分離・CoG matching・tail除去後にも残る。beam energy spread、channel/strip内応答の固定不均一性、event相関noise等を、実測情報と照合する。
5. **LG tailの取得条件を調べる。** 同一eventのlayer/chip/BCIDと元binaryを少数例で比較する。

## 再現と保存物

```bash
cd ~/ScECAL_BeamTest/investigations/e_ssa_threshold_review_20260930
# 初回実行用。完了済みbaselineがある場合は上書きせず停止する。
bash run.sh 100
bash run.sh 5
# truth・感度診断のみ更新
bash extra.sh 100
# 既存集計の表を更新
bash summarize.sh
```

既存の `InspectRawCalibration`, `InspectResponse`, `InspectTopology`, `HLInterceptStudy`, `HGLGTailStudy` を利用する。追加コードはC++の `TruthResponse.cc`, `Sensitivity.cc` とshell。`../result/e_ssa_threshold_review_20260930/<energy>GeV/` のbaseline/raw/calibrated/tail/truth/sensitivityに分類した。PNGのみでPDFは生成していない。一時tail cacheは終了時削除済み。

注意: `analysis/channel_response/adc_energy/run.sh` は現時点で説明文がe_ssaでも実際の既定入力がe_ssa_noHLinterceptになっていた。今回の入力はすべてmanifestでe_ssa/thresholdを明示した。この既定値は変更していない。

## 検証

- Raw/calib event ID・CellID照合、校正式の閉じ、HG/LG選択を両energyで確認。
- tailのhit/event別集計一致を両energyで確認。
- calibratedのMC符号付き和はbaselineと完全一致、切片0診断和の最大差1.78e-15 GeV。
- 追加Sensitivityの選別event数・平均は独立したInspectResponseと一致。
- truth/decodeのentry数、MC event identity、CellID対応、PrimaryEnergyを両energy各6000 eventで確認。
- raw ADC、event energy、CoGの代表PNGを目視確認。

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](../../investigations/result/e_ssa_threshold_review_20260930/)。
