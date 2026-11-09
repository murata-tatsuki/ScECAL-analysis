# HG–LG tail生成結果の確認（2026-09-30）

既存の `analysis/result/channel_response/hg_lg_tail/threshold/default/` を確認した。元ROOTの再解析・ジョブ投入・解析コード変更は行っていない。数値一覧は `../../investigations/result/hg_lg_tail_review_20260930/completed_summary.tsv`。tail率はchannel構成をそろえた選別前hitの率、分解能・集中event率は選別後の値であり、分母が異なる。

## 完了状況と整合性

- 完了：0.5, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50 GeV。各5枚のPNGと最終ROOT・TSVがあり、run.logにFinishedがある。完了分のログにError/Warningは見つからない。
- 未完了：60, 70, 80, 100, 120, 150, 200, 250 GeV。全8点で一時ROOT書込み時の `No space left on device` と、最後に `Cannot read cached event` を確認。`/tmp/scecal_hg_lg_tail.*/events.root` の容量不足による失敗。途中のevents.tsvや.root.tmpは最終結果として扱わない。
- 全完了点でvalidation.tsvのevent/channel別tail数・対象hit数は一致。ADCから再計算したHG/LG energyと保存値の最大差は4.55e-13 MeV以下。これは実装・入力の整合性であり、校正係数の物理的妥当性を保証するものではない。
- 50 GeVの最終ROOTを開き、events=219817、tail_hits=113331、選別後data=95787、tail-rich data=2586を確認。events.tsvの独立集計とも一致。
- 50 GeVのoverview・channel map・event相関・energy、および5/0.5 GeVのenergy PNGを目視確認。

## 判定条件

raw HG 800–2200 ADCの範囲で、data/MC双方に100 hit以上あるgood channelを使用。sample・channelごとの残差中央値を引き、pedestal noiseと量子化に基づく基準幅に対してz < -5を低LG tailとする。一定の切片差は除かれるが、傾き差・信号依存性は残る。同一event内の異なる物理channelに3個以上tailがあればtail-rich。

event選別はlayer 9/10でCoGのx/yが±20 mm、かつ正のgood energy合計がMC中央値の0.5–1.5倍。event energyは保存済み非負Hit_Energyの和であり、入射energy自体ではない。負LG energyを加算して総energyを下げているわけではない。今回のRMS/meanはこの診断選別での値で、正式解析の分解能と無条件に同一視しない。

## 確認できた現象

1. dataに低側tailがあり、MCでは極めて少ない。20/30/40/50 GeVのdata率は0.527/0.533/0.665/0.725%。MCは同範囲で0.000023–0.000061%程度。50 GeVのz分布は中心付近が近くても負側に長いtailが残り、一定の切片差では説明できない。
2. tailと負LG energyは同じevent群に集中する。50 GeVで選別後95787 event中2586（2.700%）がtail-rich。そこに選別後の低tail hitの89.79%、負LG energy hitの10161/10164（99.97%）が含まれる。負LGを持つeventは1718、そのうち1715がtail-rich。1 eventの低tail物理channel数は最大100。MCでは完了11点ともtail-rich eventは0。低HG領域で検出したtailと、別のLG採用hitの異常が同じeventに現れるという関連であり、tail hitそのものがLG採用されている意味ではない。
3. tail-rich除去だけではdata/MCの幅差は解消しない。RMS/meanは20 GeVで9.057→8.985%（MC4.419%）、30 GeVで7.216→7.158%（MC3.832%）、40 GeVで6.660→6.573%（MC3.516%）、50 GeVで6.078→6.005%（MC3.292%）。50 GeVの平均energyは除去により2.213709→2.216308 GeV（+0.117%）で、MC2.142146 GeVとの差も残る。tail-rich群の平均は2.120063 GeVでtail-poor群より低い。これは診断的なevent除去であり、異常hitの復元・補正の効果を測ったものではない。
4. PS 2–5 GeVでは選別後平均energyのdata/MC比が0.733–0.755。5 GeVのtail-rich除去で平均は0.164910→0.164980 GeV、MCは0.220908 GeVなので、約25%の平均応答差はこのtailでは説明できない。

## 解釈上の制約と次に切り分ける事項

- 0.5/1 GeVの選別通過率はdata 1.63/2.42%、MC 19.47/47.42%。0.5 GeVのdataは151 event、共通channelは74しかない。0.5/1 GeVでtail-richが0でも現象の不在とは言えない。PSではCoG条件・energy窓・粒子選別・入射位置による影響を分ける必要がある。2–5 GeVの平均差にも選別の違いが含まれ得る。
- PS入力にもSPS名のpedestalと共通のMIP/HL係数が使われている。closureは成立しているが、PS/SPSで係数・温度・MC条件が適切かは別の確認項目。ファイル名だけを根拠に校正ミスとは断定しない。
- 10 GeVのdata tail率は0.060%で、20 GeVの0.527%より大幅に低い。energyごとに共通channel・入力run・照射条件が異なるので、単純なenergy依存と解釈しない。共通channelとHG分布を固定したrun/BCID別比較が次の候補。
- 50 GeVのchannel mapはchip 2付近の帯が目立つ。全layerをまとめた選別前dataのchip 2の率は0.867%、他chipは0.137–0.335%。chip 2には比較範囲hitの76.6%が集中しているため、配置・信号量をそろえずに電子回路の故障箇所とは断定できない。
- event単位の相関は独立Gaussian noiseだけでは表しにくい。次は代表1点の50 GeVで、同一event内のLG残差の符号・大きさとlayer/chip/BCIDのまとまりを見れば、baseline変動、readout異常、ADCの取り違え等を切り分けられる。原因そのものは今回の図だけでは確定できない。
- 再実行前には一時領域を容量のある保存先へ変更し、同時実行数を制限する必要がある。既存コードはcache.Fill/Write失敗で直ちに停止せず、再読込み時に失敗を検出している。書込み時のエラー検出も改善候補。今回再投入はしていない。

50 GeVの主な参照：
- [HG–LG残差と多重度](../../analysis/result/channel_response/hg_lg_tail/threshold/default/50GeV/figures/hg_lg_tail_overview.png)
- [event内相関](../../analysis/result/channel_response/hg_lg_tail/threshold/default/50GeV/figures/event_tail_correlations.png)
- [event energy](../../analysis/result/channel_response/hg_lg_tail/threshold/default/50GeV/figures/event_energy.png)
- [数値集計](../../analysis/result/channel_response/hg_lg_tail/threshold/default/50GeV/summary.tsv)

資料の参照先: [索引](README.md)、[関連する図・表・実行記録](../../investigations/result/hg_lg_tail_review_20260930/)。

配置確認（2026-10-02）: 旧 `hg_lg_tail/threshold/default/50GeV/` の図3件と `summary.tsv` は、今回の整理前から見つかりませんでした。保存済みの集計表は上記の既存保存先にあります。別条件の結果への差し替えは行っていません。
