# Rawからcalibrationまでのdata/MC調査

これは代表energyで原因を切り分けるための調査。通常の解析コードは `analysis/` に残す。
今回の対象は既存の100 GeV data/MC。全energyの実行・ジョブ投入は行わない。

```bash
cd ~/ScECAL_BeamTest/investigations/raw_calibration
# 新しい出力先へ100 GeVの調査を実行（C++/ROOT + shell、PNGのみ）
bash run.sh 100 ../result/raw_calibration/100GeV_repeat
# 保存済みsampleの外れ値・BCID診断を更新
bash tails.sh 100
# 保存済み集計から概要図・表を更新し、2回の読み取り結果を照合
./SummarizeRaw ../result/raw_calibration/100GeV
```

`run.sh` の既定energyは100。`all`は受け付けない。既存のchannels.tsvがある出力先を上書きしない。
`tails.sh` と `SummarizeRaw` は対応する診断結果のみを更新する。
全energyの解析が別途依頼された場合は `analysis/channel_response/` のsubdirectoryとして作成する。

入力は `analysis/result/impact_studies/hl_intercept/threshold/production_v1/<energy>GeV/inputs/` のraw TAB calib対応表。
100 GeV dataは3 runから各2000イベント、MCは対応表の先頭6ファイルから各1000イベント。
各ファイルを10区間に分け、区間中央の連続イベントを取る。先頭だけを調べる方法ではない。
MCファイルは対応表の順序（辞書順）。入出力、実際のファイルとsample数はinputsに保存する。
イベントのRun_Num・Event_Time・TriggerID/Event_Numの一致とfull CellIDの対応を確認する。
hit残存率の解析はしない。

- `InspectRawCalibration`: pedestal・MIP・HL定数・保存温度からEH/ELを再計算し、保存値との閉じを検証。
  raw ADC、HG/LG切替、LGの上限、HG/LGの線形回帰、HitTag=0のbaselineを集計。
- `InspectTails`: 同じsampleを読み、HG/LGの残差を中央値・中心68%幅・RMS・5sigma外れ値に分解。
  gainタグ・BCID・低LGの外れ値、イベント内での同時発生を調査する。
- `SummarizeRaw`: 同じchannel集合で比較表・PNGを作成。2回のsample数・energy平均の一致と校正式の閉じを検証。

結果は `investigations/result/raw_calibration/<energy>GeV/`。
今回の解釈と限界は [調査本文](../../docs/research/raw-calibration.md) を参照。

## 選別と定義

生の校正式の検証は接続channelの全校正hit。原則layer0–29、chip5のchannel30–35は除外。
raw ADCからEH/ELを再現する定数は従来のStudyCommon.hhを共用し、既存の負MIP/fallbackの挙動も再現する。
実際のdata/MC差を比較する外れ値診断は正のMIPとcomparisonの共通good-channel maskを使用。
HG/LG比較はraw HG=800–2200 ADC。残差は `LG-pedLG-(HG-pedHG)/A-b`。
dataは定数ファイルのb、旧MCはb=0を使用し、既知の切片差を除いて形を比較する。
channel/sampleごとの中央値から `5*sqrt(sigmaLG^2+(sigmaHG/A)^2+1/12)` を超えるものを5sigma外れ値とする。
1/12はLG ADCの丸めの分散近似。HG丸めの寄与は小さく、この定義では省略。
共通比較は両sample各100 hit以上。同じchannelに同じ重み `min(Ndata,Nmc)` を与えた外れ値率も保存。

HitTag=0のbaselineは両ADCが指定pedestalの±5sigma内のもの。trigger条件・小信号・切断によるbiasがあるため、
この平均を新しいpedestalと決めつけない。MCはHitTag=1のみを保存するので同等のHitTag=0比較はできない。

イベントenergyは正のgood hitの和。showerは正の全hit和が既存MCの符号付き総energy中央値の0.5–1.5倍、
CoG20は全good hitから求めたlayer9/10それぞれの|x|,|y|<20mm。beam条件は完全には一致していない。
`tail_event_effect.tsv` のLG置換・疑わしいhit除外は強い正側不一致だけを使った診断上の仮定であり、
正しい再校正や新しい推奨カットを意味しない。負側の問題はlg_dropout_counts/event_readoutを確認する。

`lg_dropout_counts.tsv` の負残差は中央値でなく既存HL直線から max(15ADC,5*pedSigmaLG) 下側。
HGが飽和・4000へ置換された領域の残差は物理的に比較できないため、LG使用時の負残差件数だけで異常率を比較しない。
正のMIPかつLG採用時にEL<0となる件数と、未飽和overlap領域の外れ値を区別する。

閾値ファイルは現行digitizationスクリプトが指定するthreshold.rootを参考情報として記録する。
既存MC生成当時の閾値ファイルが同じだったと保証するものではなく、本調査では効率の再fit・再digitizationは行わない。
ROOT・calibration・digitizationの本番ファイルは変更しない。
