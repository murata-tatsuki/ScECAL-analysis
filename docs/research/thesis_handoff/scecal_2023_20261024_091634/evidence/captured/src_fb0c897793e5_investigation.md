# 100 GeV: raw ADCからcalibrationまでの追加調査

2026-09-30。今回は100 GeVのみ。data 3 run × 2000 events、MC 6 files × 1000 events、各6000 events。
全energyの再解析・再digitization・校正定数の変更はしていない。

**新しく確認できた最も明瞭な差は、dataで複数channelのLGが同じイベント内で低く出る現象。**
MCのHG/LG相関にはほぼ存在しない非Gaussianの低側tailがある。
ただし、それだけでは総energyのdata/MC差・分解能差は説明できない。

## 1. raw→calibrationの変換式は、data/MCとも保存値と一致

rawとcalibをRun_Num・Event_Time・TriggerID/Event_Numおよびfull CellIDで対応付け、
保存温度・pedestal・MIP・HL定数からEH/ELを再計算した。

| sample | 照合hit数 | 不一致hit | HGの最大差 [MeV] | LGの最大差 [MeV] |
|---|---:|---:|---:|---:|
| data | 3,837,758 | 0 | 1.42×10^-14 | 5.68×10^-14 |
| MC | 5,230,569 | 0 | 1.42×10^-14 | 2.27×10^-13 |

HG/LG選択も両者ともraw HG=2600で切り替わり、保存Hit_Energyと一致した。
したがって「dataだけ別の変換式を使う」「現在参照している定数では保存energyを再現できない」という差は、このsampleにはない。
これは定数や温度補正が物理的に正しいことまで保証するものではない。
[validation.tsv](../../investigations/result/raw_calibration/100GeV/validation.tsv)、[runs.tsv](../../investigations/result/raw_calibration/100GeV/runs.tsv)。

## 2. 切片差を除いても、dataのHG/LG相関に低LGのtailが残る

HGの非飽和側raw 800–2200 ADC、正のMIP、comparisonの共通good-channel maskを使った。
両sampleで各100 hit以上ある1,278 channelを比較。
残差は `r=(LG-pedLG)-(HG-pedHG)/A-b`、Aはcalibration側HL倍率。
dataは定数ファイルのb、旧MCはb=0を使い、既知の切片欠落を分離した。

| 指標 | data | MC |
|---|---:|---:|
| 共通channelのoverlap hit数 | 1,002,694 | 1,277,475 |
| 中心68%幅のchannel中央値 [LG ADC] | 2.729 | 2.532 |
| RMSのchannel中央値 [LG ADC] | 3.504 | 2.556 |
| channelごとの中央値から5σ以上のtail | 6,666 hit | 1 hit |
| うち負側 / 正側 | 6,623 / 43 | 0 / 1 |
| channelの重みをそろえたtail率 | 0.6563% | 0.000114% |

σは `sqrt(pedSigmaLG^2+(pedSigmaHG/A)^2+1/12)`。1/12はLG ADC丸めの分散近似。
重みはchannelごとに `min(Ndata,Nmc)` とし、両sampleに同じ重みを用いた。
同じchannelでもADC分布まで完全に同一にした比較ではない。

中心部分の幅の差は約8%で、RMSの大きな差は主として少数のtailによる。
例: CellID820008のRMSはdata8.79 / MC2.71 ADCだが、中心68%幅は3.04 / 2.67 ADC。
このchannelは両sampleのoverlapにGainTag=0がなく、gainタグの混合だけでは説明できない。
GainTag=1のみをrunごとに分けても、dataの広い残差は残る。

単にMCのGaussian noise幅を全体に大きくする変更は、中心成分を広げすぎる可能性がある。
まずtailの生成条件・同時発生を調べ、中心成分とtailを分ける必要がある。
[比較表](../../investigations/result/raw_calibration/100GeV/comparison_summary.tsv)、[channel別](../../investigations/result/raw_calibration/100GeV/robust_residuals.tsv)、[run別](../../investigations/result/raw_calibration/100GeV/residuals_by_run.tsv)。

![HG/LG残差の例](../../investigations/result/raw_calibration/100GeV/figures/hg_lg_residual_tails.png)

図は±120 ADCを表示。統計表は表示範囲外のtailも含む。

## 3. 低LGのtailはイベント内で強く集中し、LG採用時にも現れる

イベント診断ではraw HG=800–2200の各hitについて、上の直線から
`max(15 ADC,5*pedSigmaLG)` より下を「低LG外れ値」と定義した。
この定義は前節のchannel/sample中央値基準の5σと異なる。
同一イベントに3 hit以上を持つものを診断上の「同時発生イベント」とした。

- data全6000イベント中117イベント。MCは0イベント。
- shower+CoG20後のdata5118イベント中112イベント（2.19%）。MC5015イベントでは0。
- dataの正MIP・good channelで、LG採用時のenergyが負となるhitは全体719、選別後705。
  全719 hitが上記117イベント、選別後705 hitが上記112イベントに集中した。
- MCの同条件の負LG hitは全体21、選別後16。これをdataの同時発生tailと同じ機構とは断定しない。

具体例はRun105のraw/calib entry4973（BCID148）。overlap229 hitのうち123 hitが低LG外れ値、
LG採用100 hitのうち14 hitが負energy。このイベントの隣entry4972/4974にはどちらもoverlapの低LG外れ値がない。
単独channelの恒常的なHL ratioのずれや、独立なGaussian noiseだけでは説明しにくい構造である。
ただし、electronics・hold/readout条件・decodeのどれが原因かはまだ確定していない。

BCID<100だけに限定した問題ではない。dataの低LG外れ値6285 hit中6210 hitはBCID>=100。
次の調査では、該当イベントのchip/layer別分布と元binaryのHG/LG wordを確認し、
decodeの並び違いか、取得時のLG低下かを切り分けるのが有効。
MCのBCIDは-1固定なので、dataと同じタイミング条件の再現は現在の保存情報では行えない。

### total energyへの影響の大きさ

イベント集合を診断上分けた結果。energyは正のgood hit和。

| sample / 選択 | events | 平均 [GeV] | RMS [GeV] | RMS/mean |
|---|---:|---:|---:|---:|
| data: shower+CoG20 | 5118 | 4.31008 | 0.23120 | 5.364% |
| data: 同時発生112イベント | 112 | 4.12742 | 0.26787 | 6.490% |
| data: 上記を除く | 5006 | 4.31417 | 0.22865 | 5.300% |
| MC: shower+CoG20 | 5015 | 4.12170 | 0.11533 | 2.798% |

同時発生イベントは平均energyが低い。ただし、除外によるdata RMS/meanの変化は
5.364%→5.300%（相対約1.2%）で、残るdata/MCの分解能差は大きい。
従って、この現象は実在する差だが、それだけがtotal resolutionの原因ではない。
これは同じカット後の標本RMSによる診断であり、Gaussian core fitによる正式な分解能測定ではない。
イベントの物理条件も完全には一致していないため、この除外を本番カットとして採用したわけではない。
LGが低く出る現象はdata energyを下げる向きであり、data中央energyがMCより高いことの単独原因にもならない。
[event_readout.tsv](../../investigations/result/raw_calibration/100GeV/event_readout.tsv)、[集計](../../investigations/result/raw_calibration/100GeV/event_readout_summary.tsv)、[低LGの実例](../../investigations/result/raw_calibration/100GeV/lg_dropout_examples.tsv)。

![イベント内での同時発生](../../investigations/result/raw_calibration/100GeV/figures/event_low_lg_coincidence.png)

## 4. MCのHG上限処理には低信号での不自然な切替がある

現行実装は `HG > XMax-600` でHGを4000に置換し、校正側はHG>=2600でLGを選ぶ。
接続6300 channel中2531 channelで設定値XMax-600が2600未満、うち正MIP・good channelは1119。
HL定数欠損73 channelの既定値は、全channel最小の303 ADCになっている。

実際の保存MCでも、HG=4000なのにLGから推定したHGが2600未満のhitが、正MIP・good channelで2589 hitあった。
例はCellID500000: HG上限303 ADC、sampleの全60 hitがHG=4000。
ただし、この判定はLG noiseも含む推定で、境界に近いhitには曖昧さがある。
上限が303など極端に小さい例は、境界付近のnoiseだけでは説明できない。

HL ratioの欠損時fallbackもdigiは全channel平均、calibrationはSiPM group別平均で、73 channelに差がある。
ただし共通goodかつ正MIPの範囲では該当するのはCellID500000の1 channelだけで、このsampleのenergy寄与は約21.8 MeV/6000 events。

低HGでLG側へ移る影響には、既知の切片offsetとLG noiseが含まれる。
切片offsetだけの平均寄与は、shower+CoG20のMCで約0.289 MeV/event（totalの約0.007%）と小さい。
MCのHG4000の山自体をdataとの差の原因とみなすのではなく、実際にLGへ切り替わる位置と影響を見るべき。
[early_switch_channels.tsv](../../investigations/result/raw_calibration/100GeV/early_switch_channels.tsv)、[tail_event_effect.tsv](../../investigations/result/raw_calibration/100GeV/tail_event_effect.tsv)。

## 5. pedestalと温度補正順序

dataのHitTag=0で、HG/LGとも指定pedestalの±5σ内を使うと、4157 good正MIP channelのbaseline差中央値は
HG +1.79 ADC、LG +4.49 ADCだった。
現在の固定pedestalとbeam中のHitTag=0分布が一致しているとは限らない。
ただし、HitTag=0にはtrigger選択・小信号・±5σ切断のbiasがある。
この数字をそのまま新しいpedestalとして差し引くことはできない。HL切片との整合性も必要。
MCはHitTag=1しか保存しないので、同じHitTag=0分布の直接比較はできない。

また、現在のLG校正は `f*(C*A*qLG+B)`、温度補正をHL変換後に行う式は `f*C*(A*qLG+B)`。
両式の差は `f*(C-1)*B`。shower+CoG20のdataで後者に変えた場合の総energy差は
平均−3.22 MeV/event（約−0.075%）、MCは温度20℃なのでほぼ0。
これは順序を変えるべきと断定したものではなく、現在の切片定義の下での感度評価。
今回の大きなdata/MC差の主因となる量ではない。[event_summary.tsv](../../investigations/result/raw_calibration/100GeV/event_summary.tsv)。

## 6. 今回のsampleで主要因とならなかったもの

- LG ADCの3000への張り付き: 正MIP・good channelのdata/MCとも0 hit。LG ADC>3000も0。
  少なくとも今回の100 GeV sampleでは、この上限による差は観測していない。高energyは未調査。
- 保存Hit_EnergyのHG/LG分岐の実装差: 見つからない。
- 校正定数の読み違いによるraw→saved energyの不一致: 見つからない。
- 強い正側HG/LG不一致（HG採用中に補正後LG>2*HG、差>3 MeV）: data全体34 hit、選別後3 hit、MC0。
  選別後に仮にLGへ置換しても平均の変化は約0.0045 MeV/eventで小さい。
  今回重要なのは、この正側よりも負側のイベント内同時発生。

## 次の優先順位

1. Run105 entry4973などを使い、低LGがchip/layer・読み出し順・取得条件のどれにまとまるか調べる。
   元binaryが使えるならdecode前後を数イベントだけ照合する。単純にMC noiseを増やす修正は先にしない。
2. HG/LG定数の評価では非Gaussian tailを分離し、中心の傾き・切片・温度依存を確認する。
3. MCのHG上限の品質条件を見直し、低統計・欠損channelを極端な最小値で埋めない方法を検討する。
   ただし本sampleで定量化できた平均影響は小さい。
4. 以前確認したMIP定数の採用差・光量/crosstalk規格化・中央大信号応答・beam位置の差の調査を継続する。
   今回のraw由来のtailだけで大きな分解能差は埋まらない。

詳細な定義と再実行方法は [README](../../investigations/raw_calibration/README.md) を参照。
raw/calibデータ、校正定数、本番digitizationソースは変更していない。

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](../../investigations/result/raw_calibration/100GeV/)。
