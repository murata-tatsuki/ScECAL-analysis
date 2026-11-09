# impact studyを踏まえた追加調査（2026-09-30）

切片の問題だけではdata/MC差は説明できない。追加調査では、(1)負のMIP定数を採用する校正上の問題、(2)MIP定数の採用方針とMC光量規格化、(3)イベント構成と位置分布、(4)中央の大信号と周辺hitの差を区別できた。SiPM飽和が残差の唯一の原因だとは、まだ断定できない。

対象は `hl_intercept/threshold/production_v1` の旧MC。最新ソースにはHL切片のdigitizationへの追加があるが、その変更後のMCの評価ではない。元のcalibration・digitization・入力ROOTは変更していない。

**切片0でMC内部のHG/LG接続は改善するが、dataとの一致は保証されない。**

100 GeVの正のMIP定数・共通good channelのうち、raw HG 800–2200のdata/MCが各100 hit以上ある3,530 channelを集計した。重みはchannelごとの両sample hit数の小さい方。
平均 `E_LG-E_HG` は data −0.0470 MeV、元MC +0.7077 MeV、切片0 MC −0.000208 MeV。
MCの切片由来の段差を取り除く判断は支持される。一方、data側の残差形状は別途調べる必要がある。これはイベントのCoG選別をしないgain overlap領域の診断であり、総energyの比較表とは選択が異なる。
参照: [gain_consistency_100GeV.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/gain_consistency_100GeV.tsv)。

**1. 確認できた実装上の問題：失敗したMIP fitが校正に入る。**

`Calibration.cxx` の `MIPFactorExtract_simulation` は主にχ²/NDFで採否を決め、MPVやNDFの正値を確認していない。
入力には `LandauMPV=ChiSquare=NDF=-10` の失敗値がある。この比は1なのでquality cutを通る。
接続6,300 channel中565 channelで校正MIPが負のままになる（10 µm:490、15 µm:75）。
例はCellID 10030,10031,10032。負の定数はfallback用平均の計算にも入る。

実際の保存hitを調べると、負energyのほとんどがこの負MIPのchannelに対応する。
100 GeV、後述のshower+CoG20選択で、イベント当たりの負energy寄与はdata −0.0166 GeV、MC −0.1994 GeV。
250 GeVではdata −0.0378 GeV、MC −0.5136 GeV。単なるpedestal noiseの負側とは扱えない。

これは**impact studyで符号付きの全hit和を使った場合の平均・分布を歪める原因**。
既存comparisonは負hitを除くため、これだけで既存の正値hit比較の不一致を説明できるわけではない。
対処はMPV>0、NDF>0、finiteなどを採用条件にし、平均も正常なfitだけから計算すること。
その際、data/MCの同じchannel・同じ定数方針で影響を評価する。今回、本番コードの修正はしていない。
参照: [mip_constants_audit.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/mip_constants_audit.tsv)、各energyの `validation.tsv` と `event_metrics.tsv`。

**2. MCの低信号スケールに、切片とは別の問題候補がある。**

現在のdigitizationは負MPVだけを除外して個別MPVを使う一方、calibrationはχ²/NDFで個別値を棄却し、group平均に置き換える。
その結果、6,300 channel中2,422 channelで有効MIP定数が異なる。校正のfallbackは1,857 channel。

| CellID | calibration MIP [ADC] | digitization MIP [ADC] | digi/calib |
|---|---:|---:|---:|
| 240008 | 397.917 | 296.381 | 0.745 |
| 420008 | 81.577 | 87.598 | 1.074 |
| 820008 | 81.577 | 105.272 | 1.290 |
| 920012 | 81.577 | 99.203 | 1.216 |

100 GeVのgood hitのenergyのうち、data約87.4%、MC約87.1%がfallback MIPのchannelから来る。
channel数の割合だけでなく、信号を担うchannelの定数品質が重要である。
ただし、digiが真の応答、calibrationが推定値を表す設計なら定数差自体は必ずしもバグではない。
今回のようにquality cutで棄却した実測値をdigiでそのまま信頼してよいかを確認する必要がある。

さらにMC truthとraw ADCを照合した。100 GeVの先頭1,000イベント、truth hit E=1–2 MeV、10 µm channelについて：

| channel分類 | HG ADC / digi MIPによる予測 | HG ADC / calib MIPによる予測 |
|---|---:|---:|
| digi/calib MIP一致 | 1.102 | 1.102 |
| digi/calib MIP不一致 | 1.103 | 1.329 |

予測は `(truth E / 0.305) × MIP_ADC`。HG pedestalを差し引き、raw HG<2400で比較した。
1 GeVの先頭10,000イベントでも、MIP一致channelで約1.105、不一致channelではcalib予測比約1.363となる。
この増幅は保存済みMCに実在する。
`MCDigi.cxx` は `lightYieldperMIP=MIPPeakValue/SiPMGain` の後にpCrosstalk=0.1/0.15のPixelCorrectionを行い、MIP値を1.1/1.15で補正する行はコメントアウトされている。
実測MIP MPVがcrosstalk込みであるなら、同じ増幅を二度加えている可能性が高い。

次は、digi直後のmuon/MIP ADCのMPVが入力MIP値を再現するかを確認するべき。
crosstalkを単に消すのではなく、実測MIPに合う一次光電子数へ光量を規格化し、crosstalkと飽和は別に評価する。
この低信号MC過大応答はPS側の差と向きが合うが、高energyの「data中央energyが大きい」ことの単独説明にはならない。
参照: [100GeV/truth_response.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/100GeV/truth_response.tsv)、[1GeV/truth_response.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/1GeV/truth_response.tsv)、`truth_channels.tsv`。

**3. 全イベント平均は背景・低energy成分の混合に強く依存する。**

全イベントで `E_sum < 0.1 × median(MC nominal E_sum)` の割合：

| Beam | data | MC |
|---|---:|---:|
| 1 GeV | 30.35% | 0.31% |
| 5 GeV | 35.67% | ほぼ0% |
| 100 GeV | 10.96% | 0% |
| 200 GeV | 48.03% | 0% |
| 250 GeV | 13.24% | 0% |

200 GeVはRun186で48.21%、Run187で47.86%。一つのrunだけの異常ではない。
サンプルで調べた低energyイベントは1,871件中、hitが完全に空のものは15件。
682件は「20 layer以上に正のgood hitがあり、good hit数がactive layer数の2倍以下」というtracklike条件を満たす。
平均hit数31.7、active layer数15.5で、単なる空イベントだけではない。
muonなどの混入、showerを作らない通過粒子、trigger/beam条件の調査が必要だが、粒子同定はこの解析だけではできない。
1 GeVでは負MIPによるenergy相殺も低energy成分に含まれる。低energy割合をそのまま背景粒子率とは呼ばない。

100 GeVの符号付き平均はdata 3.911、MC 4.135 GeVだが、主ピークはおおよそdata 4.456、MC 4.206 GeVで差の向きが逆になる。
主ピークGaussianは非Gaussian性によるχ²/NDFが大きい場合があるため、正式な応答や分解能の最終fitではない。
参照: [event_populations.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/event_populations.tsv)、[data_runs.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/data_runs.tsv)、[低energy割合](../../investigations/result/response_diagnostics/threshold/production_v1/figures/low_energy_fraction_vs_beam.png)。

**4. PSの位置分布はMCと大きく異なる。**

1 GeVのshower選択後、dataのlayer9/10のCoGは広く分布し、MCは中心に集中する。
同じ診断sampleでCoG20を通る割合はdata 117/1599=7.3%、MC 9477/19834=47.8%。
5 GeVはdata 64.6%、MC 79.9%。5 GeV dataの平均CoGはlayer9 x≈−4.8 mm、layer10 y≈−9.3 mm、MCはほぼ0。
これは単純な一定倍率ではなく、beam位置・幅・角度、上流での散乱、shower分布等も合わせて調べる必要があることを示す。
CoG差は応答やthresholdにも依存するので、そのまま入射beam profileそのものと同一視はしない。

さらにPS dataの保存温度は約27–29℃、MCは20℃。実装されている10 µmの温度係数ではPSに約9–11%の減少が入る。
温度補正の基準と符号は元のMIP fit条件と合わせて確認すべきであり、温度を一律20℃に置換すれば正しいという意味ではない。
この量だけでは、CoG選別後にも残る1/5 GeVの25–33%程度の差は説明しきれない。
参照: [1GeVの位置分布](../../investigations/result/response_diagnostics/threshold/production_v1/1GeV/figures/shower_position.png)、[5GeV/topology.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/5GeV/topology.tsv)。

**5. 中央の大信号差と周辺のhit数差が残る。**

比較条件をそろえて、負hit・共通bad channelを除外、shower選択、layer9/10のCoG20を要求した。
中央はchannel中心の `|x|,|y|<22.5 mm`。下表は同じイベント選択での平均の比。

| Beam | 中央Nhit data/MC | 中央E data/MC | 中央E data/切片0 MC | 周辺E data/MC |
|---|---:|---:|---:|---:|
| 30 GeV | 0.974 | 1.033 | 1.071 | 0.897 |
| 100 GeV | 0.991 | 1.105 | 1.132 | 0.800 |
| 250 GeV | 1.012 | 1.152 | 1.170 | 0.714 |

100 GeVの全good hit数はdata/MC≈0.873だが、中央だけならほぼ一致する。
Nhitを同じbinにそろえても中央energy差は残る。100 GeV、Nhit=190–191のbinではdata 3.694、MC 3.321、切片0 MC 3.243 GeV。
大信号側の寄与も異なり、100 GeVのgood energyのうち1 hitが30 MeV以上の成分はdata61.5%、MC53.7%。
古い12 MeVまでの1 hit histogramでは、この重要な領域を十分比較できない。

中央と周辺の差が逆向きなので、global gainを一つ変えるだけでは合わない。
次の候補は、位置・layer分布をさらに合わせた上での大信号応答（SiPM飽和、LG応答、MIP定数）と、周辺のthreshold/検出効率。
同じNhit・CoG20でも位置分布とlayer分布は厳密には一致していないため、「SiPM飽和が原因」との断定はまだできない。
参照: [controlled_response_ratios.tsv](../../investigations/result/response_diagnostics/threshold/production_v1/controlled_response_ratios.tsv)、[100GeVの同じhit数での比較](../../investigations/result/response_diagnostics/threshold/production_v1/100GeV/figures/energy_at_fixed_nhit.png)、[layer別比較](../../investigations/result/response_diagnostics/threshold/production_v1/100GeV/figures/layer_response.png)。

**優先する次の検証**

1. MIP fitの失敗値を採用しない共通ルールを定義し、正常値だけでfallbackを計算する。まず同じhit集合上の校正差を比較し、その後threshold込みのMCへ進む。
2. muon/MIPのtruth→digitization→calibrationで規格化を閉じる。入力MIPと生成MIP ADCが一致するように光量・crosstalkの定義を統一する。
3. PSのbeam/shower位置分布と200 GeVの低energy成分を調べ、電子選別と同じ空間条件でdata/MCを比較する。
4. その状態で中央の大信号・周辺の小信号を分け、飽和とthresholdを独立に変える。現時点でSiPM pixel数だけを調整して一致させるのは早い。

**集計範囲と検証**

全19energyのイベント和・run別集計は全イベント。hit単位の詳細は0.5,1,5,30,100,200,250 GeV。
PS dataは全イベント、SPS dataは各ファイル2000イベントを10区間から抽出。MCは各energy計20000イベントを全ファイルから抽出した。
位置診断は1,5,100,200 GeV、data各ファイル2000イベント。厳密な抽出位置・選択はコードREADMEを参照。
PSのshower energy cutはdata主ピークを部分的に切るため、選択依存がある。ここでの条件付き平均をそのまま全電子の応答としない。
詳細MCの符号付き総energyは元のpaired_eventsと完全一致、切片0総和の差は最大3.6×10^-15 GeV。
HG/LG分岐が曖昧なhit、memory付きIDによるmask差は今回のsampleでは0。
signed=positive+negative、channel和とイベント平均の一致も検証した。
校正定数・コード・comparison maskの記録は `../../investigations/result/response_diagnostics/threshold/production_v1/inputs` に保存した。

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](../../investigations/result/response_diagnostics/threshold/production_v1/)。
