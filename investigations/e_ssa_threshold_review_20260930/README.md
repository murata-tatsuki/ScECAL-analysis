> **追加統合調査（2026-10-13）：** [報告と定量結果](../../docs/research/e-ssa-threshold.md)、[追加結果](../result/e_ssa_threshold_review_20260930/20261013/joint_response/)。muonによる校正検証coverage、光量・threshold・大信号応答の組合せ試験、位置・run依存除去後のlayer共分散を追加した。既存20261013は完了済みで上書きしない。

追加解析は下記の順序・前提で実行する。LABELは既存のintegrated・digitization比較が完了した結果ラベル。

~~~bash
# CT入力補正の比較完了後。12 task：光量5/20/100、他の組合せ100 GeV
bash joint_digitization_study.sh LABEL submit
# 全task完了後に比較
bash run_joint_response_review.sh LABEL

# default・threshold_MPVwoCTのmuon解析とCT比較の完了後
bash run_muon_calibration_impact_review.sh LABEL

# joint比較完了後。5/20/100のlayer共分散と位置/runの交差検証
bash layer_covariance_review.sh LABEL

# PNGの再描画のみ。数値ROOT/TSVは変更しない
bash replot_joint_response.sh LABEL
bash layer_covariance_review.sh LABEL --plots-only
~~~

出力はjoint_response内のdigitization/、muon_calibration/、layer_covariance/。主な表は各々joint_summary.tsv・joint_impact.tsv、coverage.tsv・impact.tsv、summary.tsv・impact.tsv・covariance.tsv。layer_covarianceの対象別event表はdataまたはmc→試料→条件→20mm→energy、比較図はcomparison→MC条件→20mm→energyに保存する。図はPNGのみ、2250×1500。

LayerCovarianceReview.ccは、元のevent選別・CoGを固定し、MC calibの必要な4 branchのみ再読取りする。layerごとの位置・run/file予測を未使用blockへ適用する加法残差であり、校正energyや真のleak補正ではない。layer平均を保ち、共分散和がtotal相対分散と一致することを検証する。300回のblock bootstrapは回帰係数固定で、run間系統差を含まない。再描画ではdata/MCおよび上下の共通色範囲を使い、数値ファイルのSHA不変を検証する。

---

> **温度モデルの小規模試算（2026-10-13）：** [結果と解釈](../../docs/research/e-ssa-threshold.md)。`bash temperature_resolution_estimate.sh LABEL` は保存済みevent表による5/20/100 GeVの一次近似、`bash data_temperature_impact.sh LABEL` は同じ選別dataのrawを再読取りする逆応答係数の感度診断。結果は `digitization/temperature_estimate/` とその `data_calibration/`。100 GeVのみの実digitizer試験は `temperature_digitization_pilot.sh` のhelp/usageを参照し、3条件×2truth filesの完了後にreviewする。結果は `digitization/temperature_pilot/`。既存20261013は実行済み、production変更はしていない。CT/PDE曲線とgain勾配は感度を見る仮定で、実測fitによる最終モデルではない。

# ADC・hit数・resolutionの統合再調査（2026-10-13）

[調査報告](../../docs/research/e-ssa-threshold.md)、[結果一式](../result/e_ssa_threshold_review_20260930/20261013/)。
最新electron data/MCの5・20・100 GeVを再読取りし、raw ADC、channel占有率、event hit数、layer共分散を比較した。さらにCT入力規格化を修正したelectron MCを同じtruthから各energy2,000eventだけ再生成した。大信号応答とthresholdの追加試験は100 GeVだけである。

結果は同じラベル配下の `integrated/`、`muon_ct/`、`raw_adc/`、`event_components/`、`digitization/` にまとめる。図はPNGのみ。production入力・本番digitizer・校正定数を上書きしない。科学的ROOT/TSV、入力一覧、checksum、使用ソース・binary、実行記録は結果内に保持する。

新しい結果ラベルで再現する場合は、このディレクトリから次の順に実行する。以下の `LABEL` は既存結果のない任意の英数字ラベルに置き換える。`20261013` の結果は既に生成済みなので再投入は不要。

```bash
# 最新data/MCの入力監査とevent表（Slurm投入、完了を待つ）
bash integrated_review.sh LABEL

# integrated完了後のevent量・channel ADC比較
bash event_distribution_review.sh LABEL
bash raw_adc_review.sh ../../investigations/result/e_ssa_threshold_review_20260930/LABEL all

# muon Run44の旧/new CT入力比較（各20,000event）
bash run_muons_ct_review.sh ../../investigations/result/e_ssa_threshold_review_20260930/LABEL/muon_ct

# 同じelectron truthのCT入力修正試験（6タスク、完了を待つ）
bash ct_digitization_study.sh LABEL submit
bash run_digitization_impact_review.sh LABEL

# 100 GeVのみの有効pixel数・threshold診断（8タスク、完了を待つ）
bash digitization_model_scan.sh LABEL submit
bash run_digitization_impact_review.sh LABEL modelscan

# modelscan完了後、全条件に共通する同じhitの低/高信号ADC比較
bash run_digitization_impact_review.sh LABEL adccontrols
```

`integrated_review.sh` とdigitizationの `submit` はジョブを投入して戻る。後続処理は先行ジョブの完了後に実行する。CT比較は新しいraw/calibの全COMPLETEを確認し、modelscanはCT比較完了と8件の新しいraw/calibのCOMPLETEを確認する。完了済みのCT比較・modelscan・通常のraw ADC出力は上書きせず停止する。raw ADCのPNGだけ再描画する場合は第3引数に `--plots-only` を付ける。

主な数値表は `integrated/response_comparison.tsv`、`raw_adc/comparison/20mm/<energy>GeV/{event_moments,common_channels,impact}.tsv`、`event_components/{distribution_summary,variance_components,crossvalidated_impact,reweighting_impact}.tsv`、`muon_ct/{adc_distribution,paired_truth_detection}.tsv`、`digitization/{summary,impact,truth_hit_acceptance}.tsv`。modelscanの追加表・図は `digitization/model_scan/{model_summary,model_impact,truth_reference,adc_controls,adc_control_coverage}.tsv` と `model_scan.png/root` に保存する。

impactは特記しない限り `R=RMS/mean` の相対分散差 `R_data²−R_MC²` の減少率。Gaussian中心幅の改善率とは異なる。shower CoGや後方信号による診断は関連した変動を示すもので、独立な入射位置や真の漏れenergyの測定ではない。100 GeVのMC→data Nhit重み付けは有効統計が約55eventに落ちるため、見かけの分解能一致を採用しない。

---

> **2026-10-12統合再調査：** 最新data/MCを読み直す入口は `bash integrated_review.sh <結果ラベル>`。5/20/100 GeVのみ、dataのproduction低LG tail除外、MCは現行e_ssa/thresholdを使用する。同じ結果ラベルは上書きしない。実行済みラベルは `20261012`。
>
> `IntegratedReview.cc` は現在の定数採用を `CurrentParameters.hh` で再現し、raw/calib closure、energy窓あり/なし、位置・runの交差検証、温度係数、後方信号割合、hit cut、layer共分散、保存済み本番resolutionを比較する。最新truth/ADCは各energyの*_0.rootを全event読み、`truth_scale.tsv`へ保存する。truth診断単独は `IntegratedReview --truth-only <出力先>`。
>
> 結果：`investigations/result/e_ssa_threshold_review_20260930/<ラベル>/integrated/`。`summary.tsv` の `baseline_*` 列はnominal母集団の参考量、`models.tsv` のsample_selectionはenergy窓あり/なしを区別する。`response_comparison.tsv`、`gain_overlap.tsv`、`leak_proxy_rates.tsv`は同じ集計から作る比較表。図はPNGのみ。後方割合は真の漏れ比、outer_fine_stripは長軸も含めた外縁帯を表すものではない。
>
> [最新報告](../../docs/research/e-ssa-threshold.md)、[20261012結果](../result/e_ssa_threshold_review_20260930/20261012/integrated/)。旧10月2日/4日の入力再利用コードは当時の校正採用規則を保持しており、最新MCへそのまま適用しない。

> **2026-10-04 分解能impact追加：** `bash resolution_impact.sh` で検証済み5/20/100 GeV event表から位置・run・MIP・温度の幅への影響を比較します。指定MCの実体、元decode/calibの更新情報、校正定数、production tail判定を照合してから実行します。更新があれば停止します。
>
> 結果は [resolution_impact](../result/e_ssa_threshold_review_20260930/resolution_impact/)、結論・定義は [調査報告](../../docs/research/e-ssa-threshold.md)。`summary.tsv` はRMS/68%幅/Gaussian幅、`impact.tsv` は相対分散差の減少率、`gaussian_scan.tsv` はbin/fit範囲の感度、`components.tsv` はHG/LGおよび中央/周辺の分散分解です。位置補正は100-entry blockを単位に5foldで交差検証し、全energyへの適用や本番校正の変更はしていません。

> **2026-10-02更新：** 5 / 20 / 100 GeVを最新 `e_ssa/threshold` で再調査しました。最新の結論は [調査報告](../../docs/research/e-ssa-threshold.md)、結果は [20261002](../result/e_ssa_threshold_review_20260930/20261002/) を参照してください。下記の9月30日の数値は履歴です。
>
> 現行の実行方法：`bash run.sh 20 <label>`。既存baselineからの追加診断は `bash extra.sh 20 <label>`、3 energyの集約は `bash summarize.sh <label>`。実行済みlabelは `20261002`。未完了の新labelで5・20・100を各々実行してください。
>
> 最新の主比較は全sampleのproduction HG/LG tail判定を用いたevent除外後です。`Sensitivity.cc`は入力同一性・calibration closure・layer共分散の和を検証します。条件別出力は `<label>/sensitivity/{data,mc}/<条件>/<CoG>mm/<energy>GeV/`、比較表・図は `<label>/sensitivity/comparison/20mm/<energy>GeV/`。`PrepareReference.cc`は必要なMC event基準のみを作成し、個別channel画像は生成しません。

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
